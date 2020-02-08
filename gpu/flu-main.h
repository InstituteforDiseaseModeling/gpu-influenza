#pragma once

#include <cstdint>
#include <string>

#include <curand.h>

#include "sim-types.h"

#define DEFAULT_POP     (32)
#define DEFAULT_CONS     (4)
#define DEFAULT_EXPS     (0)
#define DEFAULT_INFS     (1)
#define DEFAULT_STEPS   (10)
#define DEFAULT_SEED    (20190101)

#define DEFAULT_EXPOSED_DURATION_MEAN       (1.609438f)  // log(5)
#define DEFAULT_EXPOSED_DURATION_STD        (0.4054651f) // log(1.5)
#define DEFAULT_INFECTIOUS_DURATION_MEAN    (1.609438f)  // log(5)
#define DEFAULT_INFECTIOUS_DURATION_STD     (0.4054651f) // log(1.5)
#define DEFAULT_R0                          (0.9162907f) // log(2.5)
#define DEFAULT_R0_STD                      (0.4054651f) // log(1.5)

#define DEFAULT_VACCINATIONS                (0) 
#define DEFAULT_VAC_TIMES                   (NULL) 
#define DEFAULT_VAC_EFFECTIVENESS           (NULL) 
#define DEFAULT_VAC_COVERAGE                (NULL) 

#define DEFAULT_OUTPUT  ("transmissions.bin")

struct Configuration {
    uint32_t    population_size;
    uint32_t    max_connections;
    uint32_t    initial_exposures;
    uint32_t    initial_infections;
    uint32_t    simulation_steps;
    uint32_t    prng_seed;
    float       exposed_duration;
    float       exposed_duration_std;
    float       infectious_duration;
    float       infectious_duration_std;
    float       r0;
    float       r0_std;
    bool        quiet;
    bool        debug;
    std::string input_filename;
    std::string output_filename;
    std::string network_filename;
    uint32_t    vac_num;
    uint32_t*   vac_times;
    float*      vac_effectiveness;
    float*      vac_coverage;

    Configuration()
        : population_size(DEFAULT_POP)
        , max_connections(DEFAULT_CONS)
        , initial_exposures(DEFAULT_EXPS)
        , initial_infections(DEFAULT_INFS)
        , simulation_steps(DEFAULT_STEPS)
        , prng_seed(DEFAULT_SEED)
        , exposed_duration(DEFAULT_EXPOSED_DURATION_MEAN)
        , exposed_duration_std(DEFAULT_EXPOSED_DURATION_STD)
        , infectious_duration(DEFAULT_INFECTIOUS_DURATION_MEAN)
        , infectious_duration_std(DEFAULT_INFECTIOUS_DURATION_STD)
        , r0(DEFAULT_R0)
        , r0_std(DEFAULT_R0_STD)
        , quiet(false)
        , debug(false)
        , input_filename("")
        , output_filename(DEFAULT_OUTPUT)
        , network_filename("")
	, vac_num(DEFAULT_VACCINATIONS)
	, vac_times(DEFAULT_VAC_TIMES)
	, vac_effectiveness(DEFAULT_VAC_EFFECTIVENESS)
	, vac_coverage(DEFAULT_VAC_COVERAGE)
    {}
};

struct Counters {
    uint32_t    exposed_previously;
    uint32_t    exposed_currently;
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
    CON_TYPE*       beta_factor;
    curandState*    prngStates;
    IND_TYPE*       queue_a_exp;
    IND_TYPE*       queue_b_exp;
    IND_TYPE*       queue_a_inf;
    IND_TYPE*       queue_b_inf;
    TIMER_TYPE*     timers_a_exp;
    TIMER_TYPE*     timers_b_exp;
    TIMER_TYPE*     timers_a_inf;
    TIMER_TYPE*     timers_b_inf;
    INF_TYPE*       status_v;
    TxEvent*        tx_events;
    Counters*       counters;

    SimulationState()
        : population_size(0)
        , max_connections(0)
        , network_m(nullptr)
        , connect_m(nullptr)
        , beta_factor(nullptr)
        , prngStates(nullptr)
        , queue_a_exp(nullptr)
        , queue_b_exp(nullptr)
        , queue_a_inf(nullptr)
        , queue_b_inf(nullptr)
        , timers_a_exp(nullptr)
        , timers_b_exp(nullptr)
        , timers_a_inf(nullptr)
        , timers_b_inf(nullptr)
        , status_v(nullptr)
        , counters(nullptr)
    {}
};
