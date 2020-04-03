#pragma once

#include <cstdint>
#include <map>
#include <string>

#include <curand_kernel.h>  // includes curand.h which includes cuda_runtime.h

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
#define DEFAULT_VAC_TIMES                   (nullptr)
#define DEFAULT_VAC_EFFECTIVENESS           (nullptr)
#define DEFAULT_VAC_COVERAGE                (nullptr)

#define DEFAULT_IMPORTATION_MAG             (0.0f)
#define DEFAULT_IMPORTATION_DOUBLING        (0.0990210f) // log(2)/7
#define DEFAULT_IMPORTATION_DURATION_MEAN   (0.9162901f) // log(2.5)
#define DEFAULT_IMPORTATION_DURATION_STD    (0.4054651f) // log(1.5)
#define DEFAULT_RESEED_TIME                 (0)

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
    float       import_mag;
    float       import_rate;
    float       import_duration;
    float       import_duration_std;
    uint32_t    reseed_time;

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
        , import_mag(DEFAULT_IMPORTATION_MAG)
        , import_rate(DEFAULT_IMPORTATION_DOUBLING)
        , import_duration(DEFAULT_IMPORTATION_DURATION_MEAN)
        , import_duration_std(DEFAULT_IMPORTATION_DURATION_STD)
        , reseed_time(DEFAULT_RESEED_TIME)
    {}
};

struct Counters {
    uint32_t    exposed_previously;
    uint32_t    exposed_currently;
    uint32_t    infected_previously;
    uint32_t    infected_currently;
    uint32_t    recovered;
    uint32_t    transmissions;
    uint32_t    importations;
};

struct SimulationState {
    uint32_t        population_size;
    uint32_t        max_connections;
    IND_TYPE*       network_m;
    CON_TYPE*       connect_m;
    CON_TYPE*       rx_modifiers;
    CON_TYPE*       tx_modifiers;
    curandState*    prngStates;
    IND_TYPE*       exposed_individuals_odd_time_steps;
    IND_TYPE*       exposed_individuals_even_time_steps;
    IND_TYPE*       infected_individuals_odd_time_steps;
    IND_TYPE*       infected_individuals_even_time_steps;
    TIMER_TYPE*     exposure_timers_odd_time_steps;
    TIMER_TYPE*     exposure_timers_even_time_steps;
    TIMER_TYPE*     infected_timers_odd_time_steps;
    TIMER_TYPE*     infected_timers_even_time_steps;
    INF_TYPE*       status_v;
    TxEvent*        tx_events;
    Counters*       counters;
    std::map<std::string, ATTR_TYPE*> attributes;

    SimulationState()
        : population_size(0)
        , max_connections(0)
        , network_m(nullptr)
        , connect_m(nullptr)
        , rx_modifiers(nullptr)
        , tx_modifiers(nullptr)
        , prngStates(nullptr)
        , exposed_individuals_odd_time_steps(nullptr)
        , exposed_individuals_even_time_steps(nullptr)
        , infected_individuals_odd_time_steps(nullptr)
        , infected_individuals_even_time_steps(nullptr)
        , exposure_timers_odd_time_steps(nullptr)
        , exposure_timers_even_time_steps(nullptr)
        , infected_timers_odd_time_steps(nullptr)
        , infected_timers_even_time_steps(nullptr)
        , status_v(nullptr)
        , counters(nullptr)
        , attributes()
    {}
};
