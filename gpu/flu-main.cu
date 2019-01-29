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
void DisplayInfected( IND_TYPE* queue, TIMER_TYPE* timers, uint32_t count );

int main(int argc, const char* argv[])
{
    Configuration   config;
    SimulationState state;

    ProcessCommandLine(argc, argv, config);

    AllocateState(state, config);
    InitializeState(state, config);

    DisplayNetwork( state );
    DisplayInfected( state.queue_a, state.timers_a, state.counters->infected_currently );

    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = (config.population_size + block.x - 1) / block.x;
    dim3 grid(gridX, 1, 1);

    IND_TYPE*   previously_infected;
    IND_TYPE*   currently_infected;
    TIMER_TYPE* previous_timers;
    TIMER_TYPE* current_timers;

    auto start = clock();
    for (uint32_t step = 0; step < config.simulation_steps; ++step)
    {
        state.counters->infected_previously = state.counters->infected_currently;
        state.counters->infected_currently = 0;

        previously_infected = ((step & 1) ? state.queue_b : state.queue_a);
        currently_infected  = ((step & 1) ? state.queue_a : state.queue_b);
        previous_timers = ((step & 1) ? state.timers_b : state.timers_a);
        current_timers  = ((step & 1) ? state.timers_a : state.timers_b);

        grid.x = (state.counters->infected_previously + block.x - 1) / block.x;
        if ( grid.x > 0 )
        {
            flu_kernel<<<grid, block>>>(
                step,
                state.population_size, state.max_connections,
                config.infectious_duration,
                state.network_m, state.connect_m,
                previously_infected, previous_timers, state.counters->infected_previously,
                state.infected_v, state.recovered_v,
                currently_infected, current_timers, &state.counters->infected_currently,
                &state.counters->transmissions, state.tx_events,
                state.prngStates
            );
            CheckReturn(cudaDeviceSynchronize());
            CheckReturn(cudaGetLastError());

            if ( config.debug )
            {
                uint32_t w = uint32_t(log10(config.simulation_steps)) + 1;
                std::cout << std::setw(w) << step << ": " << state.counters->infected_currently << " infected." << std::endl;
            }

            DisplayInfected( currently_infected, current_timers, state.counters->infected_currently );
        }
        else
        {
            std::cout << "No remaining infected at time step " << step << '.' << std::endl;
            break;
        }
    }
    auto end = clock();
    std::cout << (1000 * uint64_t(end - start) / double(CLOCKS_PER_SEC)) << " milliseconds for simulation." << std::endl;

    WriteTransmissionEventsToFile( state, config );

    return 0;
}

void ShowHelp(const char* name)
{
    std::cout << std::endl;
    std::cout << name << ':' << std::endl;
    std::cout << "\t--population:#     - population size [" << DEFAULT_POP << ']' << std::endl;
    std::cout << "\t--connections:#    - per agent network degree [" << DEFAULT_CONS << ']' << std::endl;
    std::cout << "\t--infections:#     - number of initial infections [" << DEFAULT_INFS << ']' << std::endl;
    std::cout << "\t--time:#           - number of time steps [" << DEFAULT_STEPS << ']' << std::endl;
    std::cout << "\t--duration:#       - mean infectious duration (days) [" << DEFAULT_INFECTIOUS_DURATION_MEAN << ']' << std::endl;
    std::cout << "\t--r0:#             - mean r0 [" << DEFAULT_R0 << ']' << std::endl;
    std::cout << "\t--seed:#           - random number seed [" << DEFAULT_SEED << ']' << std::endl;
    std::cout << "\t--input:filename   - name of file with network data [ none ]" << std::endl;
    std::cout << "\t--output:filename  - name of file for transmission data [ '" << DEFAULT_OUTPUT << "' ]" << std::endl;
    std::cout << "\t--network:filename - name of file for saving network data [ none ]" << std::endl;
    std::cout << "\t--quiet            - do not display state during execution" << std::endl;
    std::cout << "\t--debug            - display debugging information" << std::endl;
    std::cout << std::endl;
    exit(-1);
}

void ValidateArguments(Configuration& config)
{
    bool okay = true;

    if ( config.population_size == 0 )         { std::cerr << "population size must be > 0"      << std::endl; okay = false; }
    if ( config.max_connections == 0 )         { std::cerr << "connections must be > 0"          << std::endl; okay = false; }
    if ( config.initial_infections == 0 )      { std::cerr << "initial infections must be > 0"   << std::endl; okay = false; }
    if ( config.simulation_steps == 0 )        { std::cerr << "number of time steps must be > 0" << std::endl; okay = false; }
    if ( config.infectious_duration <= 0.0f )  { std::cerr << "infectious duration must be > 0"  << std::endl; okay = false; }
    if ( config.r0 <= 0.0f )                   { std::cerr << "r0 must be > 0"                   << std::endl; okay = false; }

    if ( !okay ) { exit(-1); }

    return;
}

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
            else if ( strncmp(arg, "population",  len) == 0 ) { config.population_size     = atoi(value); }
            else if ( strncmp(arg, "connections", len) == 0 ) { config.max_connections     = atoi(value); }
            else if ( strncmp(arg, "infections",  len) == 0 ) { config.initial_infections  = atoi(value); }
            else if ( strncmp(arg, "time",        len) == 0 ) { config.simulation_steps    = atoi(value); }
            else if ( strncmp(arg, "duration",    len) == 0 ) { config.infectious_duration = atof(value); }  // float
            else if ( strncmp(arg, "r0",          len) == 0 ) { config.r0                  = atof(value); }  // float
            else if ( strncmp(arg, "seed",        len) == 0 ) { config.prng_seed           = atoi(value); }
            else if ( strncmp(arg, "quiet",       len) == 0 ) { config.quiet = true; }
            else if ( strncmp(arg, "debug",       len) == 0 ) { config.debug = true; }
            else if ( strncmp(arg, "input",       len) == 0 ) { config.input_filename   = std::string(value); }
            else if ( strncmp(arg, "output",      len) == 0 ) { config.output_filename  = std::string(value); }
            else if ( strncmp(arg, "network",     len) == 0 ) { config.network_filename = std::string(value); }
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

    std::cout << "Population:                       " << config.population_size          << std::endl;
    std::cout << "Max connections: ................ " << config.max_connections          << std::endl;
    std::cout << "Initial infections:               " << config.initial_infections       << std::endl;
    std::cout << "Simulation duration (time steps): " << config.simulation_steps         << std::endl;
    std::cout << "Mean infectious duration:         " << config.infectious_duration      << std::endl;
    std::cout << "R0:                               " << config.r0                       << std::endl;
    std::cout << "PRNG seed:                        " << config.prng_seed                << std::endl;
    std::cout << "Input (network) filename:         " << config.input_filename.c_str()   << std::endl;
    std::cout << "Output (transmissions) filename:  " << config.output_filename.c_str()  << std::endl;
    std::cout << "Network (save) filename:          " << config.network_filename.c_str() << std::endl;
    std::cout << "Quiet mode:                       " << (config.quiet ? "on" : "off")   << std::endl;
    std::cout << "Debugging mode:                   " << (config.debug ? "on" : "off")   << std::endl;
}

void DisplayNetwork( SimulationState& state )
{

}

void DisplayInfected( IND_TYPE* infected_queue, TIMER_TYPE* timers, uint32_t count )
{
    std::cout << "Infected (" << count << "): ";
    for (uint32_t i = 0; i < count; ++i)
    {
        std::cout << infected_queue[i] << '/' << uint32_t(timers[i]) << ' ';
    }
    std::cout << std::endl;
}

void AllocateState(SimulationState& state, Configuration& config)
{
    CheckReturn(cudaMallocManaged(&state.network_m,   config.population_size * config.max_connections * sizeof(IND_TYPE)));
    CheckReturn(cudaMallocManaged(&state.connect_m,   config.population_size * config.max_connections * sizeof(CON_TYPE)));
    CheckReturn(cudaMallocManaged(&state.prngStates,  config.population_size * sizeof(curandState)));
    CheckReturn(cudaMallocManaged(&state.queue_a,     config.population_size * sizeof(IND_TYPE)));
    CheckReturn(cudaMallocManaged(&state.queue_b,     config.population_size * sizeof(IND_TYPE)));
    CheckReturn(cudaMallocManaged(&state.timers_a,    config.population_size * sizeof(TIMER_TYPE)));
    CheckReturn(cudaMallocManaged(&state.timers_b,    config.population_size * sizeof(TIMER_TYPE)));
    CheckReturn(cudaMallocManaged(&state.infected_v,  config.population_size * sizeof(INF_TYPE)));
    CheckReturn(cudaMallocManaged(&state.recovered_v, config.population_size * sizeof(MASK_TYPE)));
    CheckReturn(cudaMallocManaged(&state.tx_events,   config.population_size * sizeof(TxEvent)));
    CheckReturn(cudaMallocManaged(&state.counters,    sizeof(Counters)));
}

void InitializePrngStates(SimulationState& state, Configuration& config);
void InitializeNetwork(SimulationState& state, Configuration& config);
void SeedInfections(SimulationState& state, Configuration& config);
void ZeroVectors(SimulationState& state);

void InitializeState(SimulationState& state, Configuration& config)
{
    state.population_size = config.population_size;
    state.max_connections = config.max_connections;
    state.counters->infected_previously = 0;
    state.counters->infected_currently  = 0;
    state.counters->recovered           = 0;
    state.counters->transmissions       = 0;

    InitializePrngStates(state, config);
    InitializeNetwork(state, config);
    SeedInfections(state, config);
    ZeroVectors(state);
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
    setup_network<<<grid, blocks>>>(state.network_m, state.connect_m, state.population_size, state.max_connections, config.infectious_duration, config.r0, state.prngStates);
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
    std::poisson_distribution<uint32_t> duration(config.infectious_duration);

    for (size_t j = 0; j < config.initial_infections; ++j)
    {
        IND_TYPE index;
        do
        {
            index = uniform_int(generator);
        } while (state.infected_v[index] != 0);
        state.queue_a[j] = index;
        state.timers_a[j] = TIMER_TYPE(std::max(duration(generator), uint32_t(1))); // At least one day of infectiousness.
        state.infected_v[index] = 1;

        if (config.debug)
        {
            std::cout << "Individual " << index << " placed at location " << j << " with infectious period " << uint32_t(state.timers_a[j]) << std::endl;
        }
    }

    state.counters->infected_currently = config.initial_infections;
}

// TODO - could do this on the GPU with a tiny kernel
void ZeroVectors(SimulationState& state)
{
    memset(state.infected_v,  0, state.population_size * sizeof(INF_TYPE));
    memset(state.recovered_v, 0, state.population_size * sizeof(MASK_TYPE));
}

void WriteTransmissionEventsToFile( SimulationState& state, Configuration& config )
{
    const char* filename = config.output_filename.c_str();
    auto file = fopen( filename, "wb" );
    std::cout << "Writing " << state.counters->transmissions << " transmission events to '" << filename << "'." << std::endl;
    auto count = fwrite( state.tx_events, sizeof(TxEvent), state.counters->transmissions, file );
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