#pragma once

#include <cstdint>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#include "sim-types.h"

__global__ void initialize_prng(curandState* prngState, uint32_t population, uint32_t seed);

__global__ void flu_vaccinate(CON_TYPE* beta_adjust, float effectiveness, float coverage,
		uint32_t popualtion, curandState* prngStates);

__global__ void initialize_vectors(INF_TYPE* status_v, CON_TYPE* beta_adjust, TxEvent* transmission_events,
	       	uint32_t population_count);

__global__ void setup_network(
    IND_TYPE* network, CON_TYPE* connections,
    uint32_t population_count, uint32_t connection_count,
    float infectious_duration, float r0, float r0_std,
    curandState* prngStates);

__global__ void flu_kernel(
    uint32_t simulation_step,       // RO
    uint32_t population_size,       // RO
    uint32_t max_connections,       // RO
    float mean_duration_inc,        // RO
    float std_duration_inc,         // RO
    float mean_duration_inf,        // RO
    float std_duration_inf,         // RO
    IND_TYPE* network_matrix,       // RO
    CON_TYPE* connections_matrix,   // RO
    CON_TYPE* beta_adjust,          // RO
    IND_TYPE* exposped_queue_in,    // RO
    TIMER_TYPE* exp_timers_in,      // RO
    uint32_t exposed_size,          // RO
    IND_TYPE* infected_queue_in,    // RO
    TIMER_TYPE* inf_timers_in,      // RO
    uint32_t infected_size,         // RO
    INF_TYPE* status_vector,        // RW
    IND_TYPE* exposed_queue_out,    // WO
    TIMER_TYPE* exp_timers_out,     // WO
    uint32_t* exposed_count,        // WO
    IND_TYPE* infected_queue_out,   // WO
    TIMER_TYPE* inf_timers_out,     // WO
    uint32_t* infected_count,       // WO
    uint32_t* transmission_count,   // RW
    TxEvent* transmission_events,   // WO
    curandState* prngStates         // RW
);
