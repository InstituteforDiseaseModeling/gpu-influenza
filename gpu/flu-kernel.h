#pragma once

#include <cstdint>
#include <cuda_runtime.h>
#include <curand_kernel.h>

__global__ void initialize_prng(curandState* pPrngState, uint32_t population, uint32_t seed);

#define IND_TYPE    uint32_t
#define CON_TYPE    float
#define TIMER_TYPE  uint8_t
#define INF_TYPE    uint32_t
#define MASK_TYPE   uint8_t

__global__ void flu_kernel(
    uint32_t population_size,       // RO
    uint32_t max_connections,       // RO
    IND_TYPE* network_matrix,       // RO
    CON_TYPE* connections_matrix,   // RO
    IND_TYPE* infected_queue_in,    // RO
    TIMER_TYPE* timers_in,          // RO
    uint32_t infected_size,         // RO
    INF_TYPE* infected_vector,      // RW
    MASK_TYPE* recovered_vector,    // RW
    IND_TYPE* infected_queue_out,   // WO
    TIMER_TYPE* timers_out,         // WO
    uint32_t* infected_count,       // WO
    curandState* pPrngStates        // RW
);
