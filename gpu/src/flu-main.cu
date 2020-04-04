#include <cuda_runtime.h>

#include <chrono>
#include <cstdio>

#include "block-size.h"
#include "check-return.h"

#include <iostream>
#include <iomanip>
#include <random>

#include "flu-kernel.h"
#include "flu-main.h"

/*

OVERVIEW:

1) Create a network (matrix) representing connectivity from infecteds to susceptibles.
   Parameterize by total population size (N) and maximum number of connections (C).

2) Create a vector (size N) representing infectiousness state of individuals.

For each time step:

3) For each infected individual, test for transmission to connected susceptibles.

4) TODO - consider multiple "simultaneous" transmissions.

5) Record source (and destination) for each transmission.

6) Update infectious state, mark recovered individuals as recovered (excluded).

Lather, rinse, repeat.

7) Save transmission information.

*/

void ProcessCommandLine(int argc, const char* argv[], Configuration& config);

void AllocateState(SimulationState& state, Configuration& config);
void InitializeState(SimulationState& state, Configuration& config, std::mt19937_64& generator);
void SaveStateToFile(SimulationState& state, Configuration& config);
void WriteTransmissionEventsToFile( SimulationState& state, Configuration& config );

void DisplayNetwork( SimulationState& state );
void DisplayInfected( IND_TYPE* queue_exp, TIMER_TYPE* timers_exp, uint32_t count_exp,
                      IND_TYPE* queue_inf, TIMER_TYPE* timers_inf, uint32_t count_inf );

int main(int argc, const char* argv[])
{
    Configuration   config;
    SimulationState state;

    ProcessCommandLine(argc, argv, config);

    // prng_seed is also used on the GPU, lets flip the bits for the CPU seed.
    std::mt19937_64 generator(config.prng_seed ^ 0xFFFFFFFF);
    AllocateState(state, config);
    InitializeState(state, config, generator);
    uint32_t vaccination = 0;
    uint32_t vac_start = (config.vac_num > 0) ? config.vac_times[0] : config.simulation_steps;

    DisplayNetwork( state );
    DisplayInfected( state.queue_a_exp, state.timers_a_exp, state.counters->exposed_currently,
                     state.queue_a_inf, state.timers_a_inf, state.counters->infected_currently );

    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = (config.population_size + block.x - 1) / block.x;
    dim3 grid(gridX, 1, 1);

    IND_TYPE*   previously_infected;
    IND_TYPE*   currently_infected;
    IND_TYPE*   previously_exposed;
    IND_TYPE*   currently_exposed;
    TIMER_TYPE* previous_timers_inf;
    TIMER_TYPE* current_timers_inf;
    TIMER_TYPE* previous_timers_exp;
    TIMER_TYPE* current_timers_exp;

    auto start = clock();
    uint32_t extinction_duration = 0;
    using param_t = std::poisson_distribution<int32_t>::param_type;
    std::poisson_distribution<int32_t> import_distribution(config.import_mag);
    for (uint32_t step = 0; step < config.simulation_steps; ++step)
    {
        state.counters->infected_previously = state.counters->infected_currently;
        state.counters->infected_currently = 0;
        state.counters->exposed_previously = state.counters->exposed_currently;
        state.counters->exposed_currently = 0;

        previously_exposed = ((step & 1) ? state.queue_b_exp : state.queue_a_exp);
        currently_exposed  = ((step & 1) ? state.queue_a_exp : state.queue_b_exp);
        previous_timers_exp = ((step & 1) ? state.timers_b_exp : state.timers_a_exp);
        current_timers_exp  = ((step & 1) ? state.timers_a_exp : state.timers_b_exp);
        previously_infected = ((step & 1) ? state.queue_b_inf : state.queue_a_inf);
        currently_infected  = ((step & 1) ? state.queue_a_inf : state.queue_b_inf);
        previous_timers_inf = ((step & 1) ? state.timers_b_inf : state.timers_a_inf);
        current_timers_inf  = ((step & 1) ? state.timers_a_inf : state.timers_b_inf);

        // with potential concurrent kernels I may need to change this?
        if ( (state.counters->infected_previously + state.counters->exposed_previously) > 0 )
        {
            // have those in the exposed compartment age a day
            grid.x = max((state.counters->exposed_previously + block.x - 1) / block.x, 1);
            flu_exposed<<<grid, block>>>(
                step,
                config.infectious_duration, config.infectious_duration_std,
                previously_exposed, previous_timers_exp, state.counters->exposed_previously,
                state.status_v,
                currently_exposed, current_timers_exp, &state.counters->exposed_currently,
                currently_infected, current_timers_inf, &state.counters->infected_currently,
                state.tx_events,
                state.prngStates
            );

            // have the infectious attempt to transmit
            grid.x = max((state.counters->infected_previously + block.x - 1) / block.x, 1);
            flu_kernel<<<grid, block>>>(
                step,
                state.population_size, state.max_connections,
                config.exposed_duration, config.exposed_duration_std,
                state.network_m, state.connect_m, state.beta_factor,
                previously_infected, previous_timers_inf, state.counters->infected_previously,
                state.status_v,
                currently_exposed, current_timers_exp, &state.counters->exposed_currently,
                currently_infected, current_timers_inf, &state.counters->infected_currently,
                &state.counters->transmissions, state.tx_events,
                state.prngStates
            );

            import_distribution.param(param_t{config.import_mag*exp(config.import_rate*step)});
            state.counters->importations = import_distribution(generator);
            grid.x = (state.population_size + block.x - 1) / block.x;
            flu_importation<<<grid, block>>>(
                step, state.population_size, &state.counters->importations,
                config.import_duration, config.import_duration_std,
                state.status_v,
                currently_exposed, current_timers_exp, &state.counters->exposed_currently,
                &state.counters->transmissions, state.tx_events,
                state.prngStates
                );

            CheckReturn(cudaDeviceSynchronize());
            CheckReturn(cudaGetLastError());

            extinction_duration = 0;

            if (vac_start == step)
            {
                flu_vaccinate<<<grid, block>>>(state.beta_factor, config.vac_effectiveness[vaccination],
                config.vac_coverage[vaccination],
                state.population_size, state.prngStates);
                ++vaccination;
                if (config.vac_num > vaccination)
                {
                    vac_start = config.vac_times[vaccination];
                }
            }

            if ( config.debug )
            {
                uint32_t w = uint32_t(log10(config.simulation_steps)) + 1;
                std::cout << std::setw(w) << step << ": " << state.counters->exposed_currently << " exposed, " << state.counters->infected_currently << " infected." << std::endl;
            }

            if ( !config.quiet )
            {
                DisplayInfected( currently_exposed, current_timers_exp, state.counters->exposed_currently,
                currently_infected, current_timers_inf, state.counters->infected_currently );
            }
        }
        else
        {
            ++extinction_duration;
            if(extinction_duration > config.reseed_time)
            {
                std::cout << "No remaining exposed or infected for " << extinction_duration << " days at time step " << step << '.' << std::endl;
                break;
            }
        }
    }
    auto end = clock();
    std::cout << (1000 * uint64_t(end - start) / double(CLOCKS_PER_SEC)) << " milliseconds for simulation." << std::endl;

    WriteTransmissionEventsToFile( state, config );

    free(config.vac_times);
    free(config.vac_effectiveness);
    free(config.vac_coverage);

    return 0;
}

void ShowHelp(const char* name)
{
    std::cout << std::endl;
    std::cout << name << ':' << std::endl;
    std::cout << "\t--population:#       - population size [" << DEFAULT_POP << ']' << std::endl;
    std::cout << "\t--connections:#      - per agent network degree [" << DEFAULT_CONS << ']' << std::endl;
    std::cout << "\t--exposures:#        - number of initial exposed [" << DEFAULT_EXPS << ']' << std::endl;
    std::cout << "\t--infections:#       - number of initial infections [" << DEFAULT_INFS << ']' << std::endl;
    std::cout << "\t--time:#             - number of time steps [" << DEFAULT_STEPS << ']' << std::endl;
    std::cout << "\t--durationE:#        - mean exposure duration (days) [" << exp(DEFAULT_EXPOSED_DURATION_MEAN) << ']' << std::endl;
    std::cout << "\t--stddevE:#          - exposure duration std dev (days) [" << exp(DEFAULT_EXPOSED_DURATION_STD) << ']' << std::endl;
    std::cout << "\t--durationI:#        - mean infectious duration (days) [" << exp(DEFAULT_INFECTIOUS_DURATION_MEAN) << ']' << std::endl;
    std::cout << "\t--stddevI:#          - infectious duration std dev (days) [" << exp(DEFAULT_INFECTIOUS_DURATION_STD) << ']' << std::endl;
    std::cout << "\t--r0:#               - mean r0 [" << exp(DEFAULT_R0) << ']' << std::endl;
    std::cout << "\t--stddevr0:#         - r0 std dev [" << exp(DEFAULT_R0_STD) << ']' << std::endl;
    std::cout << "\t--seed:#             - random number seed [" << DEFAULT_SEED << ']' << std::endl;
    std::cout << "\t--input:filename     - name of file with network data [ none ]" << std::endl;
    std::cout << "\t--output:filename    - name of file for transmission data [ '" << DEFAULT_OUTPUT << "' ]" << std::endl;
    std::cout << "\t--network:filename   - name of file for saving network data [ none ]" << std::endl;
    std::cout << "\t--quiet              - do not display state during execution" << std::endl;
    std::cout << "\t--debug              - display debugging information" << std::endl;
    std::cout << "\t--vaccinations:#     - number of vaccination events [" << DEFAULT_VACCINATIONS << "]" << std::endl;
    std::cout << "\t--vacTimes:#,#       - time points to apply vaccination [none]" << std::endl;
    std::cout << "\t--vacEffect:#,#      - effectiveness of each intervention [none]" << std::endl;
    std::cout << "\t--vacCoverage:#,#    - coverage of each intervention [none]" << std::endl;
    std::cout << "\t--importMag:#        - magitude of the importations [" << DEFAULT_IMPORTATION_MAG << "]" << std::endl;
    std::cout << "\t--importDoubleTime:# - doubling time of the importations (days) [" << DEFAULT_IMPORTATION_DOUBLING << "]" << std::endl;
    std::cout << "\t--durationImport:#   - mean exposure time (left) for imported cases (days) [" << DEFAULT_IMPORTATION_DURATION_MEAN << "]" << std::endl;
    std::cout << "\t--stddevImport:#     - exposed time std dev of the imported cases (days) [" << DEFAULT_IMPORTATION_DURATION_STD << "]" << std::endl;
    std::cout << "\t--reseedTime:#       - maximum number of days without infection before ending sim (days) [" << DEFAULT_RESEED_TIME << "]" << std::endl;
    std::cout << std::endl;
    exit(-1);
}

void ValidateArguments(Configuration& config)
{
    bool okay = true;

    if ( config.population_size == 0 )            { std::cerr << "population size must be > 0"          << std::endl; okay = false; }
    if ( config.max_connections == 0 )            { std::cerr << "connections must be > 0"              << std::endl; okay = false; }
    //if ( config.initial_exposures < 0 )           { std::cerr << "initial exposed must be >= 0"         << std::endl; okay = false; }
    //if ( config.initial_infections < 0 )          { std::cerr << "initial infections must be >= 0"      << std::endl; okay = false; }
    if ( config.simulation_steps == 0 )           { std::cerr << "number of time steps must be > 0"     << std::endl; okay = false; }
    if ( config.exposed_duration < 0.0f )         { std::cerr << "exposed duration must be >= 0"        << std::endl; okay = false; }
    if ( config.exposed_duration_std < 0.0f )     { std::cerr << "exposed duration std must be >= 0"    << std::endl; okay = false; }
    if ( config.infectious_duration <= 0.0f )     { std::cerr << "infectious duration must be > 0"      << std::endl; okay = false; }
    if ( config.infectious_duration_std < 0.0f )  { std::cerr << "infectious duration std must be >= 0" << std::endl; okay = false; }
    if ( config.r0 <= 0.0f )                      { std::cerr << "r0 must be > 0"                       << std::endl; okay = false; }
    if ( config.r0_std < 0.0f )                   { std::cerr << "r0 std must be >= 0"                  << std::endl; okay = false; }
    if ( config.vac_num < 0.0f )                  { std::cerr << "number of vaccinations must be >= 0"  << std::endl; okay = false; }
    //for(int i=0; i<config.vac_num; i++){
    //    if ( config.vac_times[i] < 0 ){
    //       std::cerr << "The time of a vaccination must be >= 0" << std::endl; okay = false;
    //    }
    //}
    for(int i=0; i<config.vac_num; i++){
        if ( config.vac_effectiveness[i] < 0.0f || config.vac_effectiveness[i] > 1.0f){
            std::cerr << "The effectiveness of a vaccine must be between 0 and 1 (0% and 100%) "<< std::endl; okay = false;
        }
    }
    for(int i=0; i<config.vac_num; i++){
        if ( config.vac_coverage[i] < 0.0f || config.vac_coverage[i] > 1.0f){
            std::cerr << "The coverage of a vaccine must be between 0 and 1 (0% and 100%) "<< std::endl; okay = false;
        }
    }
    if ( config.import_mag < 0.0f )               { std::cerr << "magitude of importations must be >= 0"      << std::endl; okay = false; }
    if ( config.import_rate < 0.0f )              { std::cerr << "doubling time of importations must be >= 0" << std::endl; okay = false; }
    if ( config.import_duration < 0.0f )          { std::cerr << "imported exposed duration must be >= 0"     << std::endl; okay = false; }
    if ( config.import_duration_std < 0.0f )      { std::cerr << "imported exposed duration std must be >= 0" << std::endl; okay = false; }
    //if ( config.reseed_time < 0 )                 { std::cerr << "reseeding time must be >= 0"                << std::endl; okay = false; }

    if ( !okay ) { exit(-1); }

    return;
}

void PrintArray(float* array, int n);
void PrintArray(uint32_t* array, int n);
void ReadArray(uint32_t expectedSize, float* outputArray, const char* arg);
void ReadArray(uint32_t expectedSize, uint32_t* outputArray, const char* arg);

void ProcessCommandLine(int argc, const char* argv[], Configuration& config)
{
    int vac_param_num = 0;
    for (size_t index = 1; index < argc; ++index)
    {
        if (strncmp(argv[index], "--", 2) == 0)
        {
            const char* arg = argv[index] + 2;
            size_t len = strlen(arg);
            len = strchr(arg, ':') ? (strchr(arg, ':') - arg) : len;
            len = strchr(arg, '=') ? (strchr(arg, '=') - arg) : len;
            const char* value = arg + len + 1;

            if ( strncmp(arg, "help", len) == 0 ) { ShowHelp(argv[0]); }
            else if ( strncmp(arg, "population",       len) == 0 ) { config.population_size         = atoi(value); }
            else if ( strncmp(arg, "connections",      len) == 0 ) { config.max_connections         = atoi(value); }
            else if ( strncmp(arg, "exposures",        len) == 0 ) { config.initial_exposures       = atoi(value); }
            else if ( strncmp(arg, "infections",       len) == 0 ) { config.initial_infections      = atoi(value); }
            else if ( strncmp(arg, "time",             len) == 0 ) { config.simulation_steps        = atoi(value); }
            else if ( strncmp(arg, "durationE",        len) == 0 ) { config.exposed_duration        = log(atof(value)); }  // float
            else if ( strncmp(arg, "stddevE",          len) == 0 ) { config.exposed_duration_std    = log(atof(value)); }  // float
            else if ( strncmp(arg, "durationI",        len) == 0 ) { config.infectious_duration     = log(atof(value)); }  // float
            else if ( strncmp(arg, "stddevI",          len) == 0 ) { config.infectious_duration_std = log(atof(value)); }  // float
            else if ( strncmp(arg, "r0",               len) == 0 ) { config.r0                      = log(atof(value)); }  // float
            else if ( strncmp(arg, "stddevr0",         len) == 0 ) { config.r0_std                  = log(atof(value)); }  // float
            else if ( strncmp(arg, "seed",             len) == 0 ) { config.prng_seed               = atoi(value); }
            else if ( strncmp(arg, "quiet",            len) == 0 ) { config.quiet                   = true; }
            else if ( strncmp(arg, "debug",            len) == 0 ) { config.debug                   = true; }
            else if ( strncmp(arg, "input",            len) == 0 ) { config.input_filename          = std::string(value); }
            else if ( strncmp(arg, "output",           len) == 0 ) { config.output_filename         = std::string(value); }
            else if ( strncmp(arg, "network",          len) == 0 ) { config.network_filename        = std::string(value); }
            else if ( strncmp(arg, "vaccinations",     len) == 0 )
            {
                config.vac_num           = atoi(value);
                config.vac_times         = (uint32_t*)malloc(config.vac_num*sizeof(uint32_t));
                config.vac_effectiveness = (float*)malloc(config.vac_num*sizeof(float));
                config.vac_coverage      = (float*)malloc(config.vac_num*sizeof(float));
                ++vac_param_num;
            }
            else if ( strncmp(arg, "vacTimes",         len) == 0 ) { ReadArray(config.vac_num, config.vac_times, value); vac_param_num++;}
            else if ( strncmp(arg, "vacEffect",        len) == 0 ) { ReadArray(config.vac_num, config.vac_effectiveness, value); vac_param_num++;}
            else if ( strncmp(arg, "vacCoverage",      len) == 0 ) { ReadArray(config.vac_num, config.vac_coverage, value); vac_param_num++;}
            else if ( strncmp(arg, "importMag",        len) == 0 ) { config.import_mag              = atof(value); }
            else if ( strncmp(arg, "importDoubleTime", len) == 0 ) { config.import_rate             = log(2)/atof(value); }
            else if ( strncmp(arg, "durationImport",   len) == 0 ) { config.import_duration         = log(atof(value)); }  // float
            else if ( strncmp(arg, "stddevImport",     len) == 0 ) { config.import_duration_std     = log(atof(value)); }  // float
            else if ( strncmp(arg, "reseedTime",       len) == 0 ) { config.reseed_time             = atoi(value); }
            else
            {
                std::cerr << "Unknown command line argument '" << argv[index] << "'." << std::endl;
            }
        }
        else
        {
            std::cerr << "Unknown command line argument '" << argv[index] << "'." << std::endl;
        }
    }

    // check that all or none of the vaccinaton parameters were passed in
    if ((vac_param_num!=0) && (vac_param_num!=4))
    {
        std::cerr << "Received some vaccine parameters but not all the required parameters." << std::endl;
        exit(-1);
    }

    ValidateArguments(config);

    std::cout << "Population:                        " << config.population_size              << std::endl;
    std::cout << "Max connections: ................  " << config.max_connections              << std::endl;
    std::cout << "Initial exposures:                 " << config.initial_exposures            << std::endl;
    std::cout << "Initial infections:                " << config.initial_infections           << std::endl;
    std::cout << "Simulation duration (time steps):  " << config.simulation_steps             << std::endl;
    std::cout << "Mean exposed duration:             " << exp(config.exposed_duration)        << std::endl;
    std::cout << "Exposed duration std dev:          " << exp(config.exposed_duration_std)    << std::endl;
    std::cout << "Mean infectious duration:          " << exp(config.infectious_duration)     << std::endl;
    std::cout << "Infectious duration std dev:       " << exp(config.infectious_duration_std) << std::endl;
    std::cout << "R0:                                " << exp(config.r0)                      << std::endl;
    std::cout << "R0 std dev:                        " << exp(config.r0_std)                  << std::endl;
    std::cout << "PRNG seed:                         " << config.prng_seed                    << std::endl;
    std::cout << "Input (network) filename:          " << config.input_filename.c_str()       << std::endl;
    std::cout << "Output (transmissions) filename:   " << config.output_filename.c_str()      << std::endl;
    std::cout << "Network (save) filename:           " << config.network_filename.c_str()     << std::endl;
    std::cout << "Quiet mode:                        " << (config.quiet ? "on" : "off")       << std::endl;
    std::cout << "Debugging mode:                    " << (config.debug ? "on" : "off")       << std::endl;
    std::cout << "Number of vaccination events:      " << config.vac_num                      << std::endl;
    std::cout << "Times of vaccination events:       "; PrintArray(config.vac_times, config.vac_num);
    std::cout << "Effectiveness of vaccine:          "; PrintArray(config.vac_effectiveness, config.vac_num);
    std::cout << "Vacination coverage:               "; PrintArray(config.vac_coverage, config.vac_num);
    std::cout << "Magitude of importations:          " << config.import_mag                   << std::endl;
    std::cout << "Doubling time for importatons:     " << log(2)/config.import_rate           << std::endl;
    std::cout << "Mean imported exposed duration:    " << exp(config.import_duration)         << std::endl;
    std::cout << "Imported exposed duration std dev: " << exp(config.import_duration_std)     << std::endl;
    std::cout << "Days without infection before end: " << config.reseed_time                  << std::endl;
}

void PrintArray(float* array, int n)
{
    for(int i=0; i<n; i++)
    {
        std::cout << array[i] << ' ';
    }
    std::cout << std::endl;
}

void PrintArray(uint32_t* array, int n)
{
    for(int i=0; i<n; i++)
    {
        std::cout << array[i] << ' ';
    }
    std::cout << std::endl;
}

void ReadArray(uint32_t expectedSize, float* outputArray, const char* arg)
{
    char* pch;
    char* copy = strdup(arg);
    int position = 0;
    pch = strtok(copy, ",");
    while (pch != NULL)
    {
        outputArray[position] = atof(pch);
        pch = strtok(NULL, ",");
        position++;
    }
    if (position == 1)
    {
        for(int i = 1; i<expectedSize; i++)
        {
            outputArray[i] = outputArray[0];
        }
    }
    else if (position != expectedSize)
    {
        std::cerr << "Mismatch between number of vaccinations asked for and number of parameters given" << std::endl;
        exit(-1);
    }
    free(copy);
}

void ReadArray(uint32_t expectedSize, uint32_t* outputArray, const char* arg)
{
    char* pch;
    char* copy = strdup(arg);
    int position = 0;
    pch = strtok(copy, ",");
    while (pch != NULL)
    {
        outputArray[position] = atoi(pch);
        pch = strtok(NULL, ",");
        position++;
    }
    if (position == 1)
    {
        for(int i = 1; i<expectedSize; i++)
        {
            outputArray[i] = outputArray[0];
        }
    }
    else if (position != expectedSize)
    {
        std::cerr << "Mismatch between number of vaccinations asked for and number of parameters given" << std::endl;
        exit(-1);
    }
    free(copy);
}

void DisplayNetwork( SimulationState& state )
{

}

void DisplayInfected( IND_TYPE* exposed_queue, TIMER_TYPE* exposed_timers, uint32_t exposed_count,
                      IND_TYPE* infected_queue, TIMER_TYPE* infected_timers, uint32_t infected_count)
{
    std::cout << "Exposed (" << exposed_count << "): ";
    for (uint32_t i = 0; i < exposed_count; ++i)
    {
        std::cout << exposed_queue[i] << '/' << uint32_t(exposed_timers[i]) << ' ';
    }
    std::cout << std::endl;
    std::cout << "Infected (" << infected_count << "): ";
    for (uint32_t i = 0; i < infected_count; ++i)
    {
        std::cout << infected_queue[i] << '/' << uint32_t(infected_timers[i]) << ' ';
    }
    std::cout << std::endl;
}

void AllocateState(SimulationState& state, Configuration& config)
{
    CheckReturn(cudaMallocManaged(&state.network_m,      config.population_size * config.max_connections * sizeof(IND_TYPE)));
    CheckReturn(cudaMallocManaged(&state.connect_m,      config.population_size * config.max_connections * sizeof(CON_TYPE)));
    CheckReturn(cudaMallocManaged(&state.beta_factor,    config.population_size * sizeof(CON_TYPE)));
    CheckReturn(cudaMallocManaged(&state.prngStates,     config.population_size * sizeof(curandState)));
    CheckReturn(cudaMallocManaged(&state.queue_a_exp,    config.population_size * sizeof(IND_TYPE)));
    CheckReturn(cudaMallocManaged(&state.queue_b_exp,    config.population_size * sizeof(IND_TYPE)));
    CheckReturn(cudaMallocManaged(&state.queue_a_inf,    config.population_size * sizeof(IND_TYPE)));
    CheckReturn(cudaMallocManaged(&state.queue_b_inf,    config.population_size * sizeof(IND_TYPE)));
    CheckReturn(cudaMallocManaged(&state.timers_a_exp,   config.population_size * sizeof(TIMER_TYPE)));
    CheckReturn(cudaMallocManaged(&state.timers_b_exp,   config.population_size * sizeof(TIMER_TYPE)));
    CheckReturn(cudaMallocManaged(&state.timers_a_inf,   config.population_size * sizeof(TIMER_TYPE)));
    CheckReturn(cudaMallocManaged(&state.timers_b_inf,   config.population_size * sizeof(TIMER_TYPE)));
    CheckReturn(cudaMallocManaged(&state.status_v,       config.population_size * sizeof(INF_TYPE)));
    CheckReturn(cudaMallocManaged(&state.tx_events,      config.population_size * sizeof(TxEvent)));
    CheckReturn(cudaMallocManaged(&state.counters,       sizeof(Counters)));
}

void InitializePrngStates(SimulationState& state, Configuration& config);
void InitializeNetwork(SimulationState& state, Configuration& config);
void SeedInfections(SimulationState& state, Configuration& config, std::mt19937_64& generator);
void ZeroVectors(SimulationState& state, uint32_t simulation_steps);

void InitializeState(SimulationState& state, Configuration& config, std::mt19937_64& generator)
{
    state.population_size = config.population_size;
    state.max_connections = config.max_connections;
    state.counters->exposed_previously  = 0;
    state.counters->exposed_currently   = 0;
    state.counters->infected_previously = 0;
    state.counters->infected_currently  = 0;
    state.counters->recovered           = 0;
    state.counters->transmissions       = 0;
    state.counters->importations        = 0;

    InitializePrngStates(state, config);
    InitializeNetwork(state, config);
    ZeroVectors(state, config.simulation_steps);
    SeedInfections(state, config, generator);
}

void InitializePrngStates(SimulationState& state, Configuration& config)
{
    auto start = clock();

    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = (state.population_size + block.x - 1) / block.x;
    dim3 grid(gridX, 1, 1);

    initialize_prng<<<grid, block>>>(state.prngStates, state.population_size, config.prng_seed);
    CheckReturn(cudaDeviceSynchronize());
    CheckReturn(cudaGetLastError());

    auto end = clock();
}

void LoadNetworkFromFile(SimulationState& state, Configuration& config);
void InitializeNetworkOnGpu(SimulationState& state, Configuration& config);

void InitializeNetwork(SimulationState& state, Configuration& config)
{
    if ( config.input_filename != "" )
    {
        LoadNetworkFromFile(state, config);
    }
    else
    {
        InitializeNetworkOnGpu(state, config);
        if ( config.network_filename != "" )
        {
            SaveStateToFile(state, config);
        }
    }
}

#include <sys/stat.h>

void LoadNetworkFromFile(SimulationState& state, Configuration& config)
{
    auto start = clock();

    const char* filename = config.input_filename.c_str();

    // See if the file is the correct size.
    struct stat buffer;
    if ( stat( filename, &buffer ) == 0 )
    {
        size_t element_count = config.population_size * config.max_connections;
        size_t expected_size = size_t((element_count * sizeof(IND_TYPE)) + (element_count * sizeof(CON_TYPE)));
        if ( buffer.st_size == expected_size )
        {
            auto file = fopen( filename, "rb" );
            if ( file != nullptr )
            {
                auto count = fread( state.network_m, sizeof(IND_TYPE), element_count, file );
                if ( count != element_count )
                {
                    std::cerr << "Error reading from '" << filename << "'. Wanted " << element_count << " entries but only read " << count << '.' << std::endl;
                    exit(-1);
                }
                count = fread( state.connect_m, sizeof(CON_TYPE), element_count, file );
                if ( count != element_count )
                {
                    std::cerr << "Error reading from '" << filename << "'. Wanted " << element_count << " entries but only read " << count << '.' << std::endl;
                    exit(-1);
                }
                fclose( file );
            }
            else
            {
                std::cerr << "Could not open '" << filename << "' for reading." << std::endl;
                exit(-1);
            }
        }
        else
        {
            std::cerr << "'" << filename << "' is not the correct size. Expected " << expected_size << " but found " << buffer.st_size << std::endl;
            exit(-1);
        }
    }
    else
    {
        std::cerr << "Could not call stat() on '" << filename << "'." << std::endl;
        exit(-1);
    }

    // Update network weights with transmissibility information
    dim3 blocks(BLOCK_SIZE, 1, 1);
    dim3 grid((state.population_size + BLOCK_SIZE - 1) / BLOCK_SIZE, 1, 1);
    init_network_from_file<<<grid, blocks>>>(state.connect_m, state.population_size, state.max_connections, config.infectious_duration, config.r0, config.r0_std, state.prngStates);
    CheckReturn(cudaDeviceSynchronize());
    CheckReturn(cudaGetLastError());

    auto end = clock();
    std::cout << (1000 * uint64_t(end - start) / double(CLOCKS_PER_SEC)) << " milliseconds to initialize network." << std::endl;
}

void InitializeNetworkOnGpu(SimulationState& state, Configuration& config)
{
    auto start = clock();

    dim3 blocks(BLOCK_SIZE, 1, 1);
    dim3 grid((state.population_size + BLOCK_SIZE - 1) / BLOCK_SIZE, 1, 1);
    setup_network<<<grid, blocks>>>(state.network_m, state.connect_m, state.population_size, state.max_connections, config.infectious_duration, config.r0, config.r0_std, state.prngStates);
    CheckReturn(cudaDeviceSynchronize());
    CheckReturn(cudaGetLastError());

    auto end = clock();
    std::cout << (1000 * uint64_t(end - start) / double(CLOCKS_PER_SEC)) << " milliseconds to initialize network." << std::endl;
}

// TODO - consider doing this on the GPU
void SeedInfections(SimulationState& state, Configuration& config, std::mt19937_64& generator)
{
    // If we are using saved state, infections may already be present.
    if ( config.initial_infections == 0 && config.initial_exposures == 0) return;

    std::uniform_int_distribution<IND_TYPE> uniform_int(0, config.population_size - 1);
    std::lognormal_distribution<double> duration(config.exposed_duration, config.exposed_duration_std);

    for (size_t j = 0; j < config.initial_exposures; ++j)
    {
        IND_TYPE index;
        do
        {
            index = uniform_int(generator);
        } while (state.status_v[index] != 0);
        state.queue_a_exp[j] = index;
        state.timers_a_exp[j] = TIMER_TYPE(std::max(uint32_t(duration(generator)+0.5), uint32_t(1))); // At least one day of infectiousness, also the +0.5 is so the cast to int rounds correctly
        state.status_v[index] = STATUS_EXPOSED;

        if (config.debug)
        {
            std::cout << "Individual " << index << " placed at location " << j << " with exposure period " << uint32_t(state.timers_a_exp[j]) << std::endl;
        }
    }

    state.counters->exposed_currently = config.initial_exposures;

    duration = std::lognormal_distribution<double>(config.infectious_duration, config.infectious_duration_std);
    for (size_t j = 0; j < config.initial_infections; ++j)
    {
        IND_TYPE index;
        do
        {
            index = uniform_int(generator);
        } while (state.status_v[index] != 0);
        state.queue_a_inf[j] = index;
        state.timers_a_inf[j] = TIMER_TYPE(std::max(uint32_t(duration(generator)+0.5), uint32_t(1))); // At least one day of infectiousness, also the +0.5 is so the cast to int rounds correctly
        state.status_v[index] = STATUS_INFECTIOUS;

        if (config.debug)
        {
            std::cout << "Individual " << index << " placed at location " << j << " with infectious period " << uint32_t(state.timers_a_inf[j]) << std::endl;
        }
    }

    state.counters->infected_currently = config.initial_infections;
}

// TODO - could do this on the GPU with a tiny kernel
void ZeroVectors(SimulationState& state, uint32_t simulation_steps)
{
    dim3 blocks(BLOCK_SIZE, 1, 1);
    dim3 grid((state.population_size + BLOCK_SIZE - 1) / BLOCK_SIZE, 1, 1);
    initialize_vectors<<<grid, blocks>>>(state.status_v, state.beta_factor, state.tx_events, state.population_size);
    CheckReturn(cudaDeviceSynchronize());
    CheckReturn(cudaGetLastError());
}

void WriteTransmissionEventsToFile( SimulationState& state, Configuration& config )
{
    const char* filename = config.output_filename.c_str();
    auto file = fopen( filename, "wb" );
    std::cout << "Writing " << state.counters->transmissions << " transmission events to '" << filename << "'." << std::endl;
    auto count = fwrite( state.tx_events, sizeof(TxEvent), state.population_size, file );
    std::cout << "Wrote   " << count << " transmission events to '" << filename << "'." << std::endl;
    fclose( file );
}

void SaveStateToFile(SimulationState& state, Configuration& config)
{
    const char* filename = config.network_filename.c_str();
    auto file = fopen( filename, "wb" );
    size_t element_count = config.population_size * config.max_connections;
    std::cout << "Writing network edges to '" << filename << "'... ";
    auto count = fwrite( state.network_m, sizeof(IND_TYPE), element_count, file );
    std::cout << "Wrote " << count << " edges." << std::endl;
    std::cout << "Writing edge weights to '" << filename << "'... ";
    count = fwrite( state.connect_m, sizeof(CON_TYPE), element_count, file );
    std::cout << "Wrote " << count << " weights." << std::endl;
    fclose( file );
}
