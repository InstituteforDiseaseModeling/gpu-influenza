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
void InitializeState(SimulationState& state, Configuration& config);
void SaveStateToFile(SimulationState& state, Configuration& config);
void WriteTransmissionEventsToFile( SimulationState& state, Configuration& config );

void DisplayNetwork( SimulationState& state );
void DisplayInfected( IND_TYPE* queue_exp, TIMER_TYPE* timers_exp, uint32_t count_exp,
	       IND_TYPE* queue_inf, TIMER_TYPE* teimers_inf, uint32_t count_inf	);

int main(int argc, const char* argv[])
{
    Configuration   config;
    SimulationState state;

    ProcessCommandLine(argc, argv, config);

    AllocateState(state, config);
    InitializeState(state, config);
    uint32_t intervention = 0;
    uint32_t intervention_start = config.simulation_steps;
    if(config.intervention_num>intervention)
    {
        intervention_start = config.intervention_times[intervention];
    }

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

        grid.x = (state.counters->infected_previously + state.counters->exposed_previously + block.x - 1) / block.x;
        if ( grid.x > 0 )
        {
            flu_kernel<<<grid, block>>>(
                step,
                state.population_size, state.max_connections,
                config.exposed_duration, config.exposed_duration_std,
                config.infectious_duration, config.infectious_duration_std,
                state.network_m, state.connect_m, state.beta_factor,
                previously_exposed, previous_timers_exp, state.counters->exposed_previously,
                previously_infected, previous_timers_inf, state.counters->infected_previously,
                state.status_v, 
                currently_exposed, current_timers_exp, &state.counters->exposed_currently,
                currently_infected, current_timers_inf, &state.counters->infected_currently,
                &state.counters->transmissions, state.tx_events,
                state.prngStates
            );
            CheckReturn(cudaDeviceSynchronize());
            CheckReturn(cudaGetLastError());

	    if (intervention_start == step)
	    {
                flu_vaccinate<<<grid, block>>>(state.beta_factor, config.intervention_effectiveness[intervention],
				config.intervention_coverage[intervention],
				state.population_size, state.prngStates);
		intervention++;
                if(config.intervention_num>intervention)
                {
                    intervention_start = config.intervention_times[intervention];
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
            std::cout << "No remaining exposed or infected at time step " << step << '.' << std::endl;
            break;
        }
    }
    auto end = clock();
    std::cout << (1000 * uint64_t(end - start) / double(CLOCKS_PER_SEC)) << " milliseconds for simulation." << std::endl;

    WriteTransmissionEventsToFile( state, config );

    free(config.intervention_times);
    free(config.intervention_effectiveness);
    free(config.intervention_coverage);

    return 0;
}

void ShowHelp(const char* name)
{
    std::cout << std::endl;
    std::cout << name << ':' << std::endl;
    std::cout << "\t--population:#       - population size [" << DEFAULT_POP << ']' << std::endl;
    std::cout << "\t--connections:#      - per agent network degree [" << DEFAULT_CONS << ']' << std::endl;
    std::cout << "\t--infections:#       - number of initial infections [" << DEFAULT_INFS << ']' << std::endl;
    std::cout << "\t--time:#             - number of time steps [" << DEFAULT_STEPS << ']' << std::endl;
    std::cout << "\t--durationE:#        - mean exposure duration (days) [" << DEFAULT_EXPOSED_DURATION_MEAN << ']' << std::endl;
    std::cout << "\t--stddevE:#          - exposure duration std dev (days) [" << DEFAULT_EXPOSED_DURATION_STD << ']' << std::endl;
    std::cout << "\t--durationI:#        - mean infectious duration (days) [" << DEFAULT_INFECTIOUS_DURATION_MEAN << ']' << std::endl;
    std::cout << "\t--stddevI:#          - infectious duration std dev (days) [" << DEFAULT_INFECTIOUS_DURATION_STD << ']' << std::endl;
    std::cout << "\t--r0:#               - mean r0 [" << DEFAULT_R0 << ']' << std::endl;
    std::cout << "\t--stddevr0:#         - r0 std dev [" << DEFAULT_R0 << ']' << std::endl;
    std::cout << "\t--seed:#             - random number seed [" << DEFAULT_SEED << ']' << std::endl;
    std::cout << "\t--input:filename     - name of file with network data [ none ]" << std::endl;
    std::cout << "\t--output:filename    - name of file for transmission data [ '" << DEFAULT_OUTPUT << "' ]" << std::endl;
    std::cout << "\t--network:filename   - name of file for saving network data [ none ]" << std::endl;
    std::cout << "\t--quiet              - do not display state during execution" << std::endl;
    std::cout << "\t--debug              - display debugging information" << std::endl;
    std::cout << "\t--interventions:#    - number of interventions [" << DEFAULT_INTERVENTIONS << "]" << std::endl;
    std::cout << "\t--interveneTimes:#,# - time points to apply interventions" << std::endl;
    std::cout << "\t--effectiveness:#,#  - effectiveness of each intervention" << std::endl;
    std::cout << "\t--coverage:#,#       - coverage of each intervention" << std::endl;
    std::cout << std::endl;
    exit(-1);
}

void ValidateArguments(Configuration& config)
{
    bool okay = true;

    if ( config.population_size == 0 )            { std::cerr << "population size must be > 0"         << std::endl; okay = false; }
    if ( config.max_connections == 0 )            { std::cerr << "connections must be > 0"             << std::endl; okay = false; }
    if ( config.initial_infections == 0 )         { std::cerr << "initial infections must be > 0"      << std::endl; okay = false; }
    if ( config.simulation_steps == 0 )           { std::cerr << "number of time steps must be > 0"    << std::endl; okay = false; }
    if ( config.exposed_duration <= 0.0f )        { std::cerr << "exposed duration must be > 0"        << std::endl; okay = false; }
    if ( config.exposed_duration_std <= 0.0f )    { std::cerr << "exposed duration std must be > 0"    << std::endl; okay = false; }
    if ( config.infectious_duration <= 0.0f )     { std::cerr << "infectious duration must be > 0"     << std::endl; okay = false; }
    if ( config.infectious_duration_std <= 0.0f ) { std::cerr << "infectious duration std must be > 0" << std::endl; okay = false; }
    if ( config.r0 <= 0.0f )                      { std::cerr << "r0 must be > 0"                      << std::endl; okay = false; }
    if ( config.r0_std < 0.0f )                   { std::cerr << "r0 std must be >= 0"                 << std::endl; okay = false; }
    if ( config.intervention_num < 0.0f )         { std::cerr << "number of interventions must be >= 0"<< std::endl; okay = false; }
    for(int i=0; i<config.intervention_num; i++){
        if ( config.intervention_times[i] < 0.0f ){
           	std::cerr << "The time of an interventions must be >= 0"<< std::endl; okay = false;
       	}
    }
    for(int i=0; i<config.intervention_num; i++){
        if ( config.intervention_effectiveness[i] < 0.0f || config.intervention_effectiveness[i] > 1.0f){
           	std::cerr << "The effectiveness of an intervention must be between 0 and 1 (0% and 100%) "<< std::endl; okay = false;
       	}
    }
    for(int i=0; i<config.intervention_num; i++){
        if ( config.intervention_coverage[i] < 0.0f || config.intervention_coverage[i] > 1.0f){
           	std::cerr << "The coverage of an intervention must be between 0 and 1 (0% and 100%) "<< std::endl; okay = false;
       	}
    }

    if ( !okay ) { exit(-1); }

    return;
}

void ReadArray(uint32_t expectedSize, float* outputArray, const char* arg);
void ReadArray(uint32_t expectedSize, uint32_t* outputArray, const char* arg);

void ProcessCommandLine(int argc, const char* argv[], Configuration& config)
{
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
            else if ( strncmp(arg, "population",     len) == 0 ) { config.population_size         = atoi(value); }
            else if ( strncmp(arg, "connections",    len) == 0 ) { config.max_connections         = atoi(value); }
            else if ( strncmp(arg, "infections",     len) == 0 ) { config.initial_infections      = atoi(value); }
            else if ( strncmp(arg, "time",           len) == 0 ) { config.simulation_steps        = atoi(value); }
            else if ( strncmp(arg, "durationE",      len) == 0 ) { config.exposed_duration        = log(atof(value)); }  // float
            else if ( strncmp(arg, "stddevE",        len) == 0 ) { config.exposed_duration_std    = log(atof(value)); }  // float
            else if ( strncmp(arg, "durationI",      len) == 0 ) { config.infectious_duration     = log(atof(value)); }  // float
            else if ( strncmp(arg, "stddevI",        len) == 0 ) { config.infectious_duration_std = log(atof(value)); }  // float
            else if ( strncmp(arg, "r0",             len) == 0 ) { config.r0                      = log(atof(value)); }  // float
            else if ( strncmp(arg, "stddevr0",       len) == 0 ) { config.r0_std                  = log(atof(value)); }  // float
            else if ( strncmp(arg, "seed",           len) == 0 ) { config.prng_seed               = atoi(value); }
            else if ( strncmp(arg, "quiet",          len) == 0 ) { config.quiet = true; }
            else if ( strncmp(arg, "debug",          len) == 0 ) { config.debug = true; }
            else if ( strncmp(arg, "input",          len) == 0 ) { config.input_filename   = std::string(value); }
            else if ( strncmp(arg, "output",         len) == 0 ) { config.output_filename  = std::string(value); }
            else if ( strncmp(arg, "network",        len) == 0 ) { config.network_filename = std::string(value); }
            else if ( strncmp(arg, "interventions",  len) == 0 )
	    {
	        config.intervention_num = atoi(value);
		config.intervention_times = (uint32_t*)malloc(config.intervention_num*sizeof(uint32_t));
		config.intervention_effectiveness = (float*)malloc(config.intervention_num*sizeof(float));
		config.intervention_coverage = (float*)malloc(config.intervention_num*sizeof(float));
	    }
            else if ( strncmp(arg, "interveneTimes", len) == 0 ) { ReadArray(config.intervention_num, config.intervention_times, value); }
            else if ( strncmp(arg, "effectiveness",  len) == 0 ) { ReadArray(config.intervention_num, config.intervention_effectiveness, value); }
            else if ( strncmp(arg, "coverage",       len) == 0 ) { ReadArray(config.intervention_num, config.intervention_coverage, value); }
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

    ValidateArguments(config);

    std::cout << "Population:                       " << config.population_size              << std::endl;
    std::cout << "Max connections: ................ " << config.max_connections              << std::endl;
    std::cout << "Initial infections:               " << config.initial_infections           << std::endl;
    std::cout << "Simulation duration (time steps): " << config.simulation_steps             << std::endl;
    std::cout << "Mean exposed duration:            " << exp(config.exposed_duration)        << std::endl;
    std::cout << "Exposed duration std dev:         " << exp(config.exposed_duration_std)    << std::endl;
    std::cout << "Mean infectious duration:         " << exp(config.infectious_duration)     << std::endl;
    std::cout << "Infectious duration std dev:      " << exp(config.infectious_duration_std) << std::endl;
    std::cout << "R0:                               " << exp(config.r0)                      << std::endl;
    std::cout << "R0 std dev:                       " << exp(config.r0_std)                  << std::endl;
    std::cout << "PRNG seed:                        " << config.prng_seed                    << std::endl;
    std::cout << "Input (network) filename:         " << config.input_filename.c_str()       << std::endl;
    std::cout << "Output (transmissions) filename:  " << config.output_filename.c_str()      << std::endl;
    std::cout << "Network (save) filename:          " << config.network_filename.c_str()     << std::endl;
    std::cout << "Quiet mode:                       " << (config.quiet ? "on" : "off")       << std::endl;
    std::cout << "Debugging mode:                   " << (config.debug ? "on" : "off")       << std::endl;
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
    if (position != expectedSize)
    {
        std::cerr << "Miss-match between number of interventions asked for and number of parameters given" << std::endl;
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
    if (position != expectedSize)
    {
        std::cerr << "Miss-match between number of interventions asked for and number of parameters given" << std::endl;
	exit(-1);
    }
    free(copy);
}

void DisplayNetwork( SimulationState& state )
{

}

void DisplayInfected( IND_TYPE* exposed_queue, TIMER_TYPE* exposed_timers, uint32_t exposed_count,
	       IND_TYPE* infected_queue, TIMER_TYPE* infected_timers, uint32_t infected_count	)
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
void SeedInfections(SimulationState& state, Configuration& config);
void ZeroVectors(SimulationState& state, uint32_t simulation_steps);

void InitializeState(SimulationState& state, Configuration& config)
{
    state.population_size = config.population_size;
    state.max_connections = config.max_connections;
    state.counters->exposed_previously  = 0;
    state.counters->exposed_currently   = 0;
    state.counters->infected_previously = 0;
    state.counters->infected_currently  = 0;
    state.counters->recovered           = 0;
    state.counters->transmissions       = 0;

    InitializePrngStates(state, config);
    InitializeNetwork(state, config);
    ZeroVectors(state, config.simulation_steps);
    SeedInfections(state, config);
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
void SeedInfections(SimulationState& state, Configuration& config)
{
    // If we are using saved state, infections may already be present.
    if ( config.initial_infections == 0 ) return;

    // prng_seed is also used on the GPU, lets flip the bits for the CPU seed.
    std::mt19937_64 generator(config.prng_seed ^ 0xFFFFFFFF);
    std::uniform_int_distribution<IND_TYPE> uniform_int(0, config.population_size - 1);
    std::lognormal_distribution<double> duration(config.infectious_duration, config.infectious_duration_std);

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
