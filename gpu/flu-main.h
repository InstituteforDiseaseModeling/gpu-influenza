#pragma once

#include <cstdint>
#include <string>

#include <curand.h>

#include "sim-types.h"

#define DEFAULT_POP     (32)
#define DEFAULT_CONS     (4)
#define DEFAULT_INFS     (1)
#define DEFAULT_STEPS   (10)
#define DEFAULT_SEED    (20190101)

#define DEFAULT_INFECTIOUS_DURATION_MEAN    (3.0f)
#define DEFAULT_R0                          (1.5f)

#define DEFAULT_OUTPUT  ("transmissions.bin")

struct Configuration {
    uint32_t    population_size;
    uint32_t    max_connections;
    uint32_t    initial_infections;
    uint32_t    simulation_steps;
    uint32_t    prng_seed;
    float       infectious_duration;
    float       r0;
    bool        quiet;
    bool        debug;
    std::string input_filename;
    std::string output_filename;
    std::string network_filename;

    Configuration()
        : population_size(DEFAULT_POP)
        , max_connections(DEFAULT_CONS)
        , initial_infections(DEFAULT_INFS)
        , simulation_steps(DEFAULT_STEPS)
        , prng_seed(DEFAULT_SEED)
        , infectious_duration(DEFAULT_INFECTIOUS_DURATION_MEAN)
        , r0(DEFAULT_R0)
        , quiet(false)
        , debug(false)
        , input_filename("")
        , output_filename(DEFAULT_OUTPUT)
        , network_filename("")
    {}
};

struct Counters {
    uint32_t    infected_previously;
    uint32_t    infected_currently;
    uint32_t    recovered;
    uint32_t    transmissions;
};

struct SimulationState {
    uint32_t        population_size;
    uint32_t        max_connections;
    IND_TYPE*       network_m;
    CON_TYPE*       connect_m;
    curandState*    prngStates;
    IND_TYPE*       queue_a;
    IND_TYPE*       queue_b;
    TIMER_TYPE*     timers_a;
    TIMER_TYPE*     timers_b;
    INF_TYPE*       infected_v;
    MASK_TYPE*      recovered_v;
    TxEvent*        tx_events;
    Counters*       counters;

    SimulationState()
        : population_size(0)
        , max_connections(0)
        , network_m(nullptr)
        , connect_m(nullptr)
        , prngStates(nullptr)
        , queue_a(nullptr)
        , queue_b(nullptr)
        , timers_a(nullptr)
        , timers_b(nullptr)
        , infected_v(nullptr)
        , recovered_v(nullptr)
        , counters(nullptr)
    {}
};
