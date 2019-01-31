#include <cstdio>

#include "flu-kernel.h"
#include "block-size.h"

__global__ void flu_kernel(
    uint32_t simulation_step,                                                   // RO
    uint32_t population_size, uint32_t connection_count,                        // RO
    float mean_duration,                                                        // RO
    IND_TYPE* network_m, CON_TYPE* connect_m,                                   // RO
    IND_TYPE* infected_q_in, TIMER_TYPE* timers_in, uint32_t queue_size,        // RO
    INF_TYPE* infected_v, MASK_TYPE* recovered_v,                               // R/W
    IND_TYPE* infected_q_out, TIMER_TYPE* timers_out, uint32_t* queue_count,    // WO
    uint32_t* transmission_count,                                               // RW
    TxEvent* transmission_events,                                               // WO
    curandState* prngStates)                                                    // R/W
{
    uint32_t infected_index = (blockDim.x * blockIdx.x) + threadIdx.x;
///    printf("Block dim (%02d, %02d), block index (%02d, %02d), thread index (%02d, %02d) => index %02d, neighbor %02d\n",
///        blockDim.x, blockDim.y, blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y, infected_index, neighbor_index);

    // Entire blocks of threads are scheduled, so the "last" block might have some threads
    // with no work to do.
    if (infected_index < queue_size)
    {
        IND_TYPE source_individual = infected_q_in[infected_index];
        IND_TYPE* neighbors = network_m + (connection_count * source_individual);
        CON_TYPE* susceptibility = connect_m + (connection_count * source_individual);
        for (uint32_t neighbor_index = 0; neighbor_index < connection_count; ++neighbor_index)
        {
            IND_TYPE target_individual = neighbors[neighbor_index];
            // TODO - ignore infected_v[target_individual] here to allow super-infection?
            if (!infected_v[target_individual] && !recovered_v[target_individual])
            {
                curandState* prngState = prngStates + source_individual;
                CON_TYPE draw = CON_TYPE(curand_uniform(prngState));
                if (draw < susceptibility[neighbor_index])
                {
                    INF_TYPE infected = atomicExch(&infected_v[target_individual], INF_TYPE(1));
                    if ( !infected )
                    {
                        uint32_t queue_index = atomicAdd((int*)queue_count, 1);

                        infected_q_out[queue_index] = target_individual;
                        int32_t duration = curand_poisson(prngState, mean_duration);
                        timers_out[queue_index] = TIMER_TYPE(max(duration, 1)); // At least one day.

                        uint32_t event_index = atomicAdd((int*)transmission_count, 1);
                        transmission_events[event_index].source    = source_individual;
                        transmission_events[event_index].target    = target_individual;
                        transmission_events[event_index].time_step = simulation_step;
                    }
                }
            }
        }

        TIMER_TYPE remaining = timers_in[infected_index] - 1;
        if (remaining > 0)
        {
            uint32_t queue_index = atomicAdd((int*)queue_count, 1);
            infected_q_out[queue_index] = source_individual;
            timers_out[queue_index] = remaining;
        }
        else
        {
            recovered_v[source_individual] = 1;
            // Increment number of recovered individuals
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
        curand_init(individual, 0, 0, prngState + individual);
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
    float infectious_duration, float r0,
    curandState* prngStates)
{
    uint32_t individual = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( individual < population_count )
    {
        IND_TYPE* neighbors = network +     (connection_count * individual);
        CON_TYPE* weights   = connections + (connection_count * individual);
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
            weights[connection] = CON_TYPE(r0 / (infectious_duration * connection_count));
        }
    }
}