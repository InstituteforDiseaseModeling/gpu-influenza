#pragma once

#include <cstdint>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#include "sim-types.h"

__global__ void initialize_prng(curandState* prngState, uint32_t population, uint32_t seed);

__global__ void setup_network(
    IND_TYPE* network, CON_TYPE* connections,
    uint32_t population_count, uint32_t connection_count,
    float infectious_duration, float r0,
    curandState* prngStates);

__global__ void flu_kernel(
    uint32_t simulation_step,       // RO
    uint32_t population_size,       // RO
    uint32_t max_connections,       // RO
    float mean_duration,            // RO
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
    uint32_t* transmission_count,   // RW
    TxEvent* transmission_events,   // WO
    curandState* prngStates         // RW
);
