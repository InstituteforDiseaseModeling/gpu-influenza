#include <cuda_runtime.h>
#include <curand.h>

#include <chrono>
#include <cstdio>

#include "block-size.h"
#include "../common/check-return.h"

#include <iostream>
#include <iomanip>
#include <random>

#include "flu-kernel.h"

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

void ProcessCommandLine(int argc, const char* argv[], uint32_t& n, uint32_t& c, uint32_t& i, uint32_t& d, float& duration, float& r, uint32_t& seed, bool& quiet, bool& debug);
void AllocateAndInitializeNetwork(IND_TYPE** network_m, CON_TYPE** connect_m, uint32_t population_size, uint32_t max_connections, float duration, float r0, uint32_t prng_seed);
void AllocateQueuesAndTimers(IND_TYPE** queue_a, IND_TYPE** queue_b, TIMER_TYPE** timers_a, TIMER_TYPE** timers_b, uint32_t population_size);

void InitializeNetwork(uint32_t population_size, uint32_t max_connections, IND_TYPE* network, CON_TYPE* connections, float duration, float r0, uint32_t seed);
void SeedInfected(uint32_t n, uint32_t i, float mean_duration, IND_TYPE* infected_q, TIMER_TYPE* timers, INF_TYPE*infected_v, uint32_t& count, uint32_t seed);
void AllocatePrngState(curandState** pPrngState, uint32_t population);
void DisplayInfectedQueue( IND_TYPE* infected_queue, TIMER_TYPE* timers, uint32_t count );

bool quiet = false;
bool debug = false;

#define DEFAULT_POP     (32)
#define DEFAULT_CONS     (4)
#define DEFAULT_INFS     (1)
#define DEFAULT_STEPS   (10)
#define DEFAULT_SEED    (20190101)

#define DEFAULT_INFECTIOUS_DURATION_MEAN    (3.0f)
#define DEFAULT_R0                          (1.5f)

int main(int argc, const char* argv[])
{
    uint32_t population_size    = DEFAULT_POP;
    uint32_t max_connections    = DEFAULT_CONS;
    uint32_t initial_infections = DEFAULT_INFS;
    uint32_t simulation_days    = DEFAULT_STEPS;
    uint32_t prng_seed          = DEFAULT_SEED;

    float infectious_duration = DEFAULT_INFECTIOUS_DURATION_MEAN;
    float r0                  = DEFAULT_R0;

    ProcessCommandLine(argc, argv,
        population_size, max_connections, initial_infections, simulation_days,
        infectious_duration, r0,
        prng_seed, quiet, debug);

    IND_TYPE* network_m = nullptr;
    CON_TYPE* connect_m = nullptr;
    AllocateAndInitializeNetwork(&network_m, &connect_m, population_size, max_connections, infectious_duration, r0, prng_seed);

    IND_TYPE* queue_a = nullptr;
    IND_TYPE* queue_b = nullptr;
    TIMER_TYPE* timers_a = nullptr;
    TIMER_TYPE* timers_b = nullptr;

    AllocateQueuesAndTimers(&queue_a, &queue_b, &timers_a, &timers_b, population_size);

    INF_TYPE* infected_v = nullptr;
    CheckReturn(cudaMallocManaged(&infected_v, population_size * sizeof(INF_TYPE)));
    memset(infected_v, 0, population_size * sizeof(INF_TYPE));

    uint32_t infected_size = 0;

    SeedInfected(population_size, initial_infections, infectious_duration, queue_a, timers_a, infected_v, infected_size, prng_seed);
    DisplayInfectedQueue( queue_a, timers_a, infected_size );

    MASK_TYPE* recovered_v = nullptr;
    CheckReturn(cudaMallocManaged(&recovered_v, population_size * sizeof(MASK_TYPE)));
    memset(recovered_v, 0, population_size * sizeof(MASK_TYPE));

    uint32_t* infected_count = nullptr;
    CheckReturn(cudaMallocManaged(&infected_count, sizeof(uint32_t)));

    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = (population_size + block.x - 1) / block.x;
    dim3 grid(gridX, 1, 1);

    curandState* pPrngStates;
    AllocatePrngState(&pPrngStates, population_size);
    initialize_prng<<<grid, block>>>(pPrngStates, population_size, prng_seed);

    IND_TYPE* infected_q_in = queue_a;
    TIMER_TYPE* timers_in = timers_a;
    IND_TYPE* infected_q_out = queue_b;
    TIMER_TYPE* timers_out = timers_b;

    auto start = clock();
    for (uint32_t t = 0; t < simulation_days; ++t)
    {
        *infected_count = 0;

        grid.x = (infected_size + block.x - 1) / block.x;
        if ( grid.x > 0 )
        {
            flu_kernel<<<grid, block>>>(
                population_size, max_connections, network_m, connect_m,
                infected_q_in, timers_in, infected_size,
                infected_v, recovered_v,
                infected_q_out, timers_out, infected_count,
                pPrngStates
            );
            CheckReturn(cudaDeviceSynchronize());
            CheckReturn(cudaGetLastError());

            if ( debug )
            {
                uint32_t w = uint32_t(log10(simulation_days)) + 1;
                std::cout << std::setw(w) << t << ": " << *infected_count << " infected." << std::endl;
            }

            DisplayInfectedQueue( infected_q_out, timers_out, *infected_count );

            // Swap queues
            auto temp1 = infected_q_in;
            infected_q_in = infected_q_out;
            infected_q_out = temp1;

            // Swap timers
            auto temp2 = timers_in;
            timers_in = timers_out;
            timers_out = temp2;

            infected_size = *infected_count;
        }
        else
        {
            std::cout << "No remaining infected at time step " << t << '.' << std::endl;
            break;
        }
    }
    auto end = clock();
    std::cout << (1000 * uint64_t(end - start) / double(CLOCKS_PER_SEC)) << " milliseconds for simulation." << std::endl;

    return 0;
}

void ShowHelp(const char* name)
{
    std::cout << std::endl;
    std::cout << name << ':' << std::endl;
    std::cout << "\t--population:#  - population size [" << DEFAULT_POP << ']' << std::endl;
    std::cout << "\t--connections:# - per agent network degree [" << DEFAULT_CONS << ']' << std::endl;
    std::cout << "\t--infections:#  - number of initial infections [" << DEFAULT_INFS << ']' << std::endl;
    std::cout << "\t--time:#        - number of time steps [" << DEFAULT_STEPS << ']' << std::endl;
    std::cout << "\t--duration:#    - mean infectious duration (days) [" << DEFAULT_INFECTIOUS_DURATION_MEAN << ']' << std::endl;
    std::cout << "\t--r0:#          - mean r0 [" << DEFAULT_R0 << ']' << std::endl;
    std::cout << "\t--seed:#        - random number seed [" << DEFAULT_SEED << ']' << std::endl;
    std::cout << "\t--quiet         - do not display state during execution" << std::endl;
    std::cout << "\t--debug         - display debugging information" << std::endl;
    std::cout << std::endl;
    exit(-1);
}

void ValidateArguments(uint32_t n, uint32_t c, uint32_t i, uint32_t d, float m, float r)
{
    bool okay = true;

    if ( n == 0 )    { std::cerr << "population size must be > 0"      << std::endl; okay = false; }
    if ( c == 0 )    { std::cerr << "connections must be > 0"          << std::endl; okay = false; }
    if ( i == 0 )    { std::cerr << "initial infections must be > 0"   << std::endl; okay = false; }
    if ( d == 0 )    { std::cerr << "number of time steps must be > 0" << std::endl; okay = false; }
    if ( m <= 0.0f ) { std::cerr << "infectious duration must be > 0"  << std::endl; okay = false; }
    if ( r <= 0.0f ) { std::cerr << "r0 must be > 0"                   << std::endl; okay = false; }

    if ( !okay ) { exit(-1); }

    return;
}

void ProcessCommandLine(int argc, const char* argv[], uint32_t& n, uint32_t& c, uint32_t& i, uint32_t& d, float& m, float& r, uint32_t& s, bool& q, bool& dbg)
{
    for (size_t index = 1; index < argc; ++index)
    {
        if (strncmp(argv[index], "--", 2) == 0)
        {
            const char* arg = argv[index] + 2;
            size_t len = strlen(arg);
            len = strchr(arg, ':') ? (strchr(arg, ':') - arg) : len;
            len = strchr(arg, '=') ? (strchr(arg, '=') - arg) : len;

            if ( strncmp(arg, "help", len) == 0 ) { ShowHelp(argv[0]); }
            else if ( strncmp(arg, "population",  len) == 0 ) { n = atoi(arg + len + 1); }
            else if ( strncmp(arg, "connections", len) == 0 ) { c = atoi(arg + len + 1); }
            else if ( strncmp(arg, "infections",  len) == 0 ) { i = atoi(arg + len + 1); }
            else if ( strncmp(arg, "time",        len) == 0 ) { d = atoi(arg + len + 1); }
            else if ( strncmp(arg, "duration",    len) == 0 ) { m = atof(arg + len + 1); }  // float
            else if ( strncmp(arg, "r0",          len) == 0 ) { r = atof(arg + len + 1); }  // float
            else if ( strncmp(arg, "seed",        len) == 0 ) { s = atoi(arg + len + 1); }
            else if ( strncmp(arg, "quiet",       len) == 0 ) { q = true; }
            else if ( strncmp(arg, "debug",       len) == 0 ) { dbg = true; }
            else
            {
                std::cerr << "Unknown command line argument '" << argv[i] << "'." << std::endl;
            }
        }
        else
        {
            std::cerr << "Unknown command line argument '" << argv[i] << "'." << std::endl;
        }
    }

    ValidateArguments(n, c, i, d, m, r);

    std::cout << "Population:                       " << n << std::endl;
    std::cout << "Max connections: ................ " << c << std::endl;
    std::cout << "Initial infections:               " << i << std::endl;
    std::cout << "Simulation duration (time steps): " << d << std::endl;
    std::cout << "Mean infectious duration:         " << m << std::endl;
    std::cout << "R0:                               " << r << std::endl;
    std::cout << "PRNG seed:                        " << s << std::endl;
    std::cout << "Quiet mode:                       " << (q ? "on" : "off") << std::endl;
    std::cout << "Debugging mode:                   " << (dbg ? "on" : "off") << std::endl;
}

void AllocateAndInitializeNetwork(IND_TYPE** network_m, CON_TYPE** connect_m, uint32_t population_size, uint32_t max_connections, float duration, float r0, uint32_t prng_seed)
{
    auto start = clock();
    CheckReturn(cudaMallocManaged(network_m, population_size * max_connections * sizeof(IND_TYPE)));
    CheckReturn(cudaMallocManaged(connect_m, population_size * max_connections * sizeof(CON_TYPE)));
    InitializeNetwork(population_size, max_connections, *network_m, *connect_m, duration, r0, prng_seed);
    auto end = clock();
    std::cout << (1000 * uint64_t(end - start) / double(CLOCKS_PER_SEC)) << " milliseconds to allocate and initialize network." << std::endl;
}

bool Contains(IND_TYPE* list, size_t count, IND_TYPE target);

// TODO - determine appropriate R0 and thus "strength" below.
// Consider "Modeling influenza epidemics and pandemics: insights into the future of swine flu (H1N1)"
//          "https://www.ncbi.nlm.nih.gov/pubmed/19545404"
// "The R0 for novel influenza A (H1N1) has recently been estimated to be between 1.4 and 1.6.
//  This value is below values of R0 estimated for the 1918-1919 pandemic strain (mean R0 approximately 2: range 1.4 to 2.8)
//  and is comparable to R0 values estimated for seasonal strains of influenza (mean R0 1.3: range 0.9 to 2.1)."

void InitializeNetwork(uint32_t n, uint32_t c, IND_TYPE* network, CON_TYPE* connect, float duration, float r0, uint32_t seed)
{
    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<IND_TYPE> uniform_int(0, n-1);
    std::uniform_real_distribution<CON_TYPE> uniform_float(0.25f, 0.75f);

    for (size_t i = 0; i < n; ++i)
    {
        IND_TYPE* entries = network + (c * i);
        CON_TYPE* strengths = connect + (c * i);
        for (size_t j = 0; j < c; ++j)
        {
            IND_TYPE susceptible = uniform_int(generator);
            while ((susceptible == IND_TYPE(i)) || Contains(entries, j, susceptible))
            {
                susceptible = uniform_int(generator);
            }
            entries[j] = susceptible;
            // TODO - consider distribution around this value? normal_float?
            CON_TYPE strength = r0 / (duration * c);    // uniform_float(generator);
            strengths[j] = strength;
        }

        if ( debug && (i < 16) )
        {
            std::cout << std::fixed<< std::setprecision(4) << std::showpoint;
            uint32_t w = uint32_t(log10(n)) + 1;
            std::cout << std::setw(w) << i << ": ";
            for (size_t j = 0; j < c; ++j)
            {
                std::cout << std::setw(w) << entries[j] << '/' << strengths[j] << ", ";
            }
            std::cout << std::endl;
        }
    }
}

bool Contains(IND_TYPE* list, size_t count, IND_TYPE target)
{
    for (size_t i = 0; i < count; ++i)
    {
        if (list[i] == target)
        {
            return true;
        }
    }

    return false;
}

void AllocateQueuesAndTimers(IND_TYPE** queue_a, IND_TYPE** queue_b, TIMER_TYPE** timers_a, TIMER_TYPE** timers_b, uint32_t population_size)
{
    CheckReturn(cudaMallocManaged(queue_a, population_size * sizeof(IND_TYPE)));
    CheckReturn(cudaMallocManaged(queue_b, population_size * sizeof(IND_TYPE)));

    CheckReturn(cudaMallocManaged(timers_a, population_size * sizeof(TIMER_TYPE)));
    CheckReturn(cudaMallocManaged(timers_b, population_size * sizeof(TIMER_TYPE)));
}


void SeedInfected(uint32_t population_size, uint32_t initial_infections, float mean_duration, IND_TYPE* infected_q, TIMER_TYPE* timers, INF_TYPE* infected_v, uint32_t& count, uint32_t seed)
{
    std::mt19937_64 generator(seed ^ 0xFFFFFFFF);
    std::uniform_int_distribution<IND_TYPE> uniform_int(0, population_size - 1);
    std::poisson_distribution<uint32_t> duration(mean_duration);

    for (size_t j = 0; j < initial_infections; ++j)
    {
        IND_TYPE index;
        do
        {
            index = uniform_int(generator);
        } while (infected_v[index] != 0);
        infected_q[j] = index;
        timers[j] = TIMER_TYPE(std::max(duration(generator), uint32_t(1))); // At least one day of infectiousness.
        infected_v[index] = 1;
        std::cout << "Individual " << index << " placed at location " << j << " with infectious period " << uint32_t(timers[j]) << std::endl;
    }

    count = initial_infections;
}

void AllocatePrngState(curandState** pPrngState, uint32_t population)
{
    // CheckReturn(cudaMallocManaged(pPrngState, population * sizeof(curandState)));
    CheckReturn(cudaMalloc(pPrngState, population * sizeof(curandState)));
    // CheckReturn(cudaMemset(*pPrngState, 0, population * sizeof(curandState)));
}

void DisplayInfectedQueue( IND_TYPE* infected_queue, TIMER_TYPE* timers, uint32_t count )
{
    if ( quiet ) return;

    std::cout << "Infected (" << count << "): ";
    for (uint32_t i = 0; i < count; ++i)
    {
        std::cout << infected_queue[i] << '/' << uint32_t(timers[i]) << ' ';
    }
    std::cout << std::endl;
}