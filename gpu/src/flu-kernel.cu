#include <cstdio>

#include "flu-kernel.h"
#include "block-size.h"

__global__ void flu_kernel(
    uint32_t simulation_step,                                                   // RO
    uint32_t population_size, uint32_t connection_count,                        // RO
    float mean_duration_inc, float stddev_duration_inc,                         // RO
    IND_TYPE* network_m, CON_TYPE* connect_m, CON_TYPE* beta_adjust,            // RO
    IND_TYPE* inf_q_in, TIMER_TYPE* inf_timers_in, uint32_t inf_queue_size,     // RO
    INF_TYPE* status_v,                                                         // R/W
    IND_TYPE* exp_q_out, TIMER_TYPE* exp_timers_out, uint32_t* exp_queue_count, // WO
    IND_TYPE* inf_q_out, TIMER_TYPE* inf_timers_out, uint32_t* inf_queue_count, // WO
    uint32_t* transmission_count,                                               // RW
    TxEvent* transmission_events,                                               // WO
    curandState* prngStates)                                                    // R/W
{
    uint32_t infected_index = (blockDim.x * blockIdx.x) + threadIdx.x;
///    printf("Block dim (%02d, %02d), block index (%02d, %02d), thread index (%02d, %02d) => index %02d, neighbor %02d\n",
///        blockDim.x, blockDim.y, blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y, infected_index, neighbor_index);

    // Entire blocks of threads are scheduled, so the "last" block might have some threads
    // with no work to do.
    if (infected_index < inf_queue_size)
    {
        IND_TYPE source_individual = inf_q_in[infected_index];
        IND_TYPE* neighbors = network_m + (connection_count * source_individual);
        CON_TYPE* susceptibility = connect_m + (connection_count * source_individual);
        for (uint32_t neighbor_index = 0; neighbor_index < connection_count; ++neighbor_index)
        {
            IND_TYPE target_individual = neighbors[neighbor_index];
            // TODO - ignore status_v[target_individual] here to allow super-infection?
            if (STATUS_SUSCEPTIBLE==status_v[target_individual])
            {
                curandState* prngState = prngStates + source_individual;
                CON_TYPE draw = CON_TYPE(curand_uniform(prngState));
                CON_TYPE beta_modifier = beta_adjust[target_individual];
                if (draw < (susceptibility[neighbor_index]*beta_modifier))
                {
                    INF_TYPE infected = atomicExch(&status_v[target_individual], INF_TYPE(STATUS_EXPOSED));
                    if (infected == STATUS_SUSCEPTIBLE)
                    {
                        uint32_t queue_index = atomicAdd((int*)exp_queue_count, 1);

                        exp_q_out[queue_index] = target_individual;
                        int32_t duration = __float2int_rn(curand_log_normal(prngState, mean_duration_inc, stddev_duration_inc));
                        exp_timers_out[queue_index] = TIMER_TYPE(max(duration, 1)); // At least one day.

                        uint32_t event_index = atomicAdd((int*)transmission_count, 1);
                        transmission_events[target_individual].source       = source_individual;
                        transmission_events[target_individual].time_exposed = simulation_step;
                    }
                }
            }
        }

        TIMER_TYPE remaining = inf_timers_in[infected_index] - 1;
        if (remaining > 0)
        {
            uint32_t queue_index = atomicAdd((int*)inf_queue_count, 1);
            inf_q_out[queue_index] = source_individual;
            inf_timers_out[queue_index] = remaining;
        }
        else
        {
            // Increment number of recovered individuals
            status_v[source_individual] = STATUS_RECOVERED;
            transmission_events[source_individual].time_recovered = simulation_step;
        }
    }
}

__global__ void flu_exposed(
    uint32_t simulation_step,                                                   // RO
    float mean_duration_inf, float stddev_duration_inf,                         // RO
    IND_TYPE* exp_q_in, TIMER_TYPE* exp_timers_in, uint32_t exp_queue_size,     // RO
    INF_TYPE* status_v,                                                         // R/W
    IND_TYPE* exp_q_out, TIMER_TYPE* exp_timers_out, uint32_t* exp_queue_count, // WO
    IND_TYPE* inf_q_out, TIMER_TYPE* inf_timers_out, uint32_t* inf_queue_count, // WO
    TxEvent* transmission_events,                                               // WO
    curandState* prngStates)                                                    // R/W
{
    uint32_t exposed_index = (blockDim.x * blockIdx.x) + threadIdx.x;
    if(exposed_index < exp_queue_size)
    {
        IND_TYPE source_individual = exp_q_in[exposed_index];
        TIMER_TYPE remaining = exp_timers_in[exposed_index] - 1;
        if (remaining > 0)
        {
            uint32_t queue_index = atomicAdd((int*)exp_queue_count, 1);
            exp_q_out[queue_index] = source_individual;
            exp_timers_out[queue_index] = remaining;
        }
        else
        {
            // Increment number of infectious individuals
            uint32_t queue_index = atomicAdd((int*)inf_queue_count, 1);
            inf_q_out[queue_index] = source_individual;
            curandState* prngState = prngStates + source_individual;
            int32_t duration = __float2int_rn(curand_log_normal(prngState, mean_duration_inf, stddev_duration_inf));
            inf_timers_out[queue_index] = TIMER_TYPE(max(duration, 1)); // At least one day.
            status_v[source_individual] = STATUS_INFECTIOUS;
            transmission_events[source_individual].time_infectious = simulation_step;
        }
    }
}

// for now all importations are exposures and the time remaining has not been shorten
__global__ void flu_importation(
    uint32_t simulation_step, uint32_t population_size,                         // RO
    uint32_t* importations,                                                      // R/W
    float mean_duration_inc, float stddev_duration_inc,                         // RO
    INF_TYPE* status_v,                                                         // R/W
    IND_TYPE* exp_q_out, TIMER_TYPE* exp_timers_out, uint32_t* exp_queue_count, // WO
    uint32_t* transmission_count,                                               // RW
    TxEvent* transmission_events,                                               // WO
    curandState* prngStates)                                                    // R/W
{
    uint32_t individual_index = (blockDim.x * blockIdx.x) + threadIdx.x;

    if(individual_index < population_size)
    {
        if(status_v[individual_index] == STATUS_SUSCEPTIBLE)
        {
            int32_t importsLeft = atomicAdd((int*)importations, -1);
            if(importsLeft > 0)
            {
                uint32_t queue_index = atomicAdd((int*)exp_queue_count, 1);

                curandState* prngState = prngStates + individual_index;
                exp_q_out[queue_index] = individual_index;
                int32_t duration = __float2int_rn(curand_log_normal(prngState, mean_duration_inc, stddev_duration_inc));
                exp_timers_out[queue_index] = TIMER_TYPE(max(duration, 1)); // At least one day.

                uint32_t event_index = atomicAdd((int*)transmission_count, 1);
                status_v[individual_index] = STATUS_EXPOSED;
                transmission_events[individual_index].time_exposed = simulation_step;
            }
        }
    }
}

__global__ void flu_vaccinate(
    CON_TYPE* beta_adjust, float effectiveness, float coverage,
    uint32_t population, curandState* prngStates)
{
    uint32_t individual = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( individual < population )
    {
        curandState* prngState = prngStates + individual;
        float draw = curand_uniform(prngState);
        if (draw < coverage)
        {
            beta_adjust[individual] = 1.0 - effectiveness;
        }
    }
}

__global__ void initialize_prng(curandState* prngState, uint32_t population, uint32_t seed)
{
    uint32_t individual = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( individual < population )
    {
        memset(prngState + individual, 0, sizeof(curandState));
        //curand_init(seed, individual, 0, prngState + individual);    // VERY SLOW!
        curand_init(individual + seed, 0, 0, prngState + individual);
    }
}

__device__ bool contains(IND_TYPE* neighbors, uint32_t count, IND_TYPE neighbor)
{
    for (uint32_t i = 0; i < count; ++i)
    {
        if ( neighbors[i] == neighbor )
        {
            return true;
        }
    }

    return false;
}

// TODO - determine appropriate R0 and thus weights below.
// Consider "Modeling influenza epidemics and pandemics: insights into the future of swine flu (H1N1)"
//          "https://www.ncbi.nlm.nih.gov/pubmed/19545404"
// "The R0 for novel influenza A (H1N1) has recently been estimated to be between 1.4 and 1.6.
//  This value is below values of R0 estimated for the 1918-1919 pandemic strain (mean R0 approximately 2: range 1.4 to 2.8)
//  and is comparable to R0 values estimated for seasonal strains of influenza (mean R0 1.3: range 0.9 to 2.1)."

__global__ void setup_network(
    IND_TYPE* network, CON_TYPE* connections,
    uint32_t population_count, uint32_t connection_count,
    float infectious_duration, float r0, float r0_std,
    curandState* prngStates)
{
    uint32_t individual = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( individual < population_count )
    {
        IND_TYPE* neighbors = network +     (connection_count * individual);
        CON_TYPE* weights   = connections + (connection_count * individual);
        curandState* prngState = prngStates + individual;
        float Rcurr = curand_log_normal(prngState, r0, r0_std);
        for (uint32_t connection = 0; connection < connection_count; ++connection)
        {
            uint32_t neighbor;
            do
            {
                neighbor = curand(prngStates + individual) % population_count;
            }
            while (contains(neighbors, connection, neighbor));

            neighbors[connection] = neighbor;
            // TODO - consider distribution around this value? Gaussian/normal?
            weights[connection] = CON_TYPE(Rcurr / (infectious_duration * connection_count));
        }
    }
}

__global__ void initialize_vectors(
    INF_TYPE* status_v, CON_TYPE* beta_adjust, TxEvent* transmission_events, uint32_t population_count)
{
    uint32_t individual = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( individual < population_count )
    {
       status_v[individual]                            = STATUS_SUSCEPTIBLE;
       beta_adjust[individual]                         = 1.0f;
       transmission_events[individual].source          = -1;
       transmission_events[individual].time_exposed    = -1;
       transmission_events[individual].time_infectious = -1;
       transmission_events[individual].time_recovered  = -1;
    }
}
