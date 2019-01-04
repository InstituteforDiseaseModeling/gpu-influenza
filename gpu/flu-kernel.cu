#include <cstdio>

#include "flu-kernel.h"
#include "block-size.h"

__global__ void flu_kernel(
    uint32_t population_size, uint32_t connection_count, IND_TYPE* network_m, CON_TYPE* connect_m,  // RO
    IND_TYPE* infected_q_in, TIMER_TYPE* timers_in, uint32_t queue_size,                            // RO
    INF_TYPE* infected_v, MASK_TYPE* recovered_v,                                                   // R/W
    IND_TYPE* infected_q_out, TIMER_TYPE* timers_out, uint32_t* queue_count,                        // WO
    curandState* pPrngStates)                                                                       // R/W
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
            // TODO - ignore infected_v[target_individual] here to track super-infection?
            if (!infected_v[target_individual] && !recovered_v[target_individual])
            {
                curandState* pPrngState = pPrngStates + source_individual;
                float draw = curand_uniform(pPrngState);
                if (draw < susceptibility[neighbor_index])
                {
                    INF_TYPE infected = atomicExch(&infected_v[target_individual], INF_TYPE(1));
                    if ( !infected )
                    {
                        uint32_t queue_index = atomicAdd((int*)queue_count, 1);

                        infected_q_out[queue_index] = target_individual;
                        int32_t duration = curand_poisson(pPrngState, 3);
                        timers_out[queue_index] = TIMER_TYPE(duration+1);   // At least one day.
                    }
                    // TODO - should still note that individual (super) infected target
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
        }
    }
}

__global__ void initialize_prng(curandState* pPrngState, uint32_t population, uint32_t seed)
{
    uint32_t individual = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( individual < population )
    {
        memset(pPrngState + individual, 0, sizeof(curandState));
        //curand_init(seed, individual, 0, pPrngState + individual);    // VERY SLOW!
        curand_init(individual, 0, 0, pPrngState + individual);
    }
}
