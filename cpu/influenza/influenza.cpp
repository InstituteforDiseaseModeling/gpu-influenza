// influenza.cpp : Defines the entry point for the console application.
//

#include <chrono>
#include <iostream>
#include <random>
#include <cstring>

// TODO - include attenuation factor for each agent

void ProcessCommandLine(int argc, const char* argv[], uint32_t& n, uint32_t& c, uint32_t& i, uint32_t& d, bool& quiet);
void InitializeNetwork(uint32_t n, uint32_t c, uint32_t* network_m, float* connect_m);
void DisplayNetwork(uint32_t n, uint32_t c, uint32_t* network_m, float* connect_m);
void SeedInfected(uint32_t n, uint32_t i, uint32_t* infected_q, uint8_t* timers, uint8_t*infected_v, uint32_t& count);
void DisplayState(uint32_t n, uint32_t* infected_v, uint8_t* timers, uint32_t count);
void Step(uint32_t n, uint32_t c, uint32_t* network_m, float* connect_m,        // RO
    uint32_t* infected_q_in, uint8_t* timers_in, uint32_t infected_size,        // RO
    uint8_t* infected_v, uint8_t* recovered_v,                                  // R/W
    uint32_t& susceptible_c, uint32_t& infectious_c, uint32_t& recovered_c,     // R/W
    uint32_t* infected_q_out, uint8_t* timers_out, uint32_t& infected_count);   // WO
void DisplayRecovered(uint32_t n, uint8_t* recovered_v);
void RadixSort(uint32_t* data, uint32_t count);
void Unique(uint32_t* data, uint32_t& count);

bool quiet = false;

int main(int argc, const char* argv[])
{
    uint32_t N = 32;
    uint32_t C = 4;
    uint32_t I = 1;
    uint32_t D = 10;

    ProcessCommandLine(argc, argv, N, C, I, D, quiet);

    std::cout << "Population:                       " << N << std::endl;
    std::cout << "Max connections: ................ " << C << std::endl;
    std::cout << "Initial infections:               " << I << std::endl;
    std::cout << "Simulation duration (time steps): " << D << std::endl;

    size_t total_bytes = 0;

    size_t size = N * C * sizeof(uint32_t);
    std::cout << "Allocating " << size << " bytes for network 'matrix'." << std::endl;
    uint32_t* network_m = (uint32_t*)malloc(size);
    total_bytes += size;

    size = N * C * sizeof(float);
    std::cout << "Allocating " << size << " bytes for connection 'matrix'." << std::endl;
    float* connect_m = (float*)malloc(size);
    total_bytes += size;

    InitializeNetwork(N, C, network_m, connect_m);
    if (!quiet && (N <= 128)) DisplayNetwork(N, C, network_m, connect_m);

    size = N * sizeof(float);
    std::cout << "Allocating " << size << " bytes for interventions vector." << std::endl;
    float* interventions = (float*)malloc(size);
    total_bytes += size;

    for (uint32_t i = 0; i < N; ++i)
    {
        interventions[i] = 1.0f;
    }

    std::cout << "Allocating " << 2 * N * sizeof(uint32_t) << " bytes for infected queues." << std::endl;
    uint32_t* queue_a = (uint32_t*)calloc(N, sizeof(uint32_t));
    uint32_t* queue_b = (uint32_t*)calloc(N, sizeof(uint32_t));
    uint32_t infected_size = 0;
    total_bytes += 2 * N * sizeof(uint32_t);

    std::cout << "Allocating " << 2 * N * sizeof(uint8_t) << " bytes for infection timers." << std::endl;
    uint8_t* timers_a = (uint8_t*)calloc(N, sizeof(uint8_t));
    uint8_t* timers_b = (uint8_t*)calloc(N, sizeof(uint8_t));
    total_bytes += 2 * N * sizeof(uint8_t);

    std::cout << "Allocating " << N * sizeof(uint8_t) << " bytes for infected vector." << std::endl;
    uint8_t* infected_v = (uint8_t*)calloc(N, sizeof(uint8_t));
    total_bytes += N * sizeof(uint8_t);

    SeedInfected(N, I, queue_a, timers_a, infected_v, infected_size);
    if (!quiet) DisplayState(N, queue_a, timers_a, infected_size);

    std::cout << "Allocating " << N * sizeof(uint8_t) << " bytes for recovered vector." << std::endl;
    uint8_t* recovered_v = (uint8_t*)calloc(N, sizeof(uint8_t));
    total_bytes += N * sizeof(uint8_t);

    std::cout << "Total bytes allocated = " << total_bytes << std::endl;

    uint32_t susceptible_c = N - I;
    uint32_t infectious_c = I;
    uint32_t recovered_c = 0;

    uint32_t* infected_q_input = queue_a;
    uint32_t* infected_q_output = queue_b;

    uint8_t* timers_in = timers_a;
    uint8_t* timers_out = timers_b;

    clock_t start = clock();
    for (size_t d = 0; d < D; ++d)
    {
        std::cout << "Time step " << d << ": " << susceptible_c << ',' << infectious_c << ',' << recovered_c << std::endl;
        uint32_t infected_count = 0;
        Step(N, C, network_m, connect_m,
            infected_q_input, timers_in, infected_size,
            infected_v, recovered_v,
            susceptible_c, infectious_c, recovered_c,
            infected_q_output, timers_out, infected_count);

        void* temp = infected_q_input;
        infected_q_input = infected_q_output;
        infected_q_output = (uint32_t*)temp;

        infected_size = infected_count;

        temp = timers_in;
        timers_in = timers_out;
        timers_out = (uint8_t*)temp;

        if (!quiet) DisplayState(N, infected_q_input, timers_in, infected_size);
    }
    clock_t end = clock();
    double elapsed = 1000 * uint64_t(end - start) / double(CLOCKS_PER_SEC);
    std::cout << "Elapsed time " << uint64_t(elapsed) << " milliseconds." << std::endl;

    if (!quiet) DisplayRecovered(N, recovered_v);
}

void ProcessCommandLine(int argc, const char* argv[], uint32_t& n, uint32_t& c, uint32_t& i, uint32_t& d, bool& q)
{
    for (size_t index = 1; index < argc; ++index)
    {
        if (strncmp(argv[index], "--", 2) == 0)
        {
            const char* arg = argv[index] + 2;
            size_t len = strlen(arg);
            len = strchr(arg, ':') ? (strchr(arg, ':') - arg) : len;
            len = strchr(arg, '=') ? (strchr(arg, '=') - arg) : len;
            if (strncmp(arg, "help", len) == 0)
            {
                std::cout << argv[0] << ':' << std::endl;
                std::cout << "\t--population:#  - population size" << std::endl;
                std::cout << "\t--connections:# - per agent network degree" << std::endl;
                std::cout << "\t--infections:#  - number of initial infections" << std::endl;
                std::cout << "\t--duration:#    - number of time steps" << std::endl;
                std::cout << "\t--quiet         - do not display state during execution" << std::endl;
                exit(-1);
            }
            else if (strncmp(arg, "population", len) == 0)
            {
                n = atoi(arg + len + 1);
            }
            else if ( strncmp(arg, "connections", len) == 0 )
            {
                c = atoi(arg + len + 1);
            }
            else if ( strncmp(arg, "infections", len) == 0 )
            {
                i = atoi(arg + len + 1);
            }
            else if (strncmp(arg, "duration", len) == 0)
            {
                d = atoi(arg + len + 1);
            }
            else if (strncmp(arg, "quiet", len) == 0)
            {
                q = true;
            }
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
}

bool Contains(uint32_t* list, size_t count, uint32_t target);

void InitializeNetwork(uint32_t n, uint32_t c, uint32_t* network, float* connect)
{
    auto start = clock();
    std::mt19937_64 generator(20181204);
    std::uniform_int_distribution<uint32_t> uniform_int(0, n-1);
    std::uniform_real_distribution<float> uniform_float(0.25f, 0.75f);

    for (size_t i = 0; i < n; ++i)
    {
        uint32_t* entries = network + (c * i);
        float* strengths = connect + (c * i);
        for (size_t j = 0; j < c; ++j)
        {
            uint32_t susceptible = uniform_int(generator);
            while ((susceptible == uint32_t(i)) || Contains(entries, j, susceptible))
            {
                susceptible = uniform_int(generator);
            }
            entries[j] = susceptible;
            float strength = 2.0f / (3 * c); // uniform_float(generator);
            strengths[j] = strength;
        }
    }
    auto end = clock();
    std::cout << "Initializing " << n << " individuals with " << c << " connections each ... ";
    std::cout << (uint64_t(end - start) * 1000 / double(CLOCKS_PER_SEC)) << " milliseconds." << std::endl;
}

bool Contains(uint32_t* list, size_t count, uint32_t target)
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

void DisplayNetwork(uint32_t n, uint32_t c, uint32_t* network, float* connect)
{
    for (size_t source = 0; source < n; ++source)
    {
        uint32_t* entries = network + (c * source);
        float* strengths = connect + (c * source);
        std::cout << source << ": ";
        for (size_t i = 0; i < c; ++i)
        {
            std::cout << entries[i] << '/' << strengths[i] << ' ';
        }
        std::cout << std::endl;
    }
}

void SeedInfected(uint32_t population_size, uint32_t initial_infections, uint32_t* infected_q, uint8_t* timers, uint8_t* infected_v, uint32_t& count)
{
    std::mt19937_64 generator(101110);
    std::uniform_int_distribution<uint32_t> uniform_int(0, population_size - 1);
    std::poisson_distribution<uint32_t> duration(3);

    for (size_t j = 0; j < initial_infections; ++j)
    {
        uint32_t index;
        do
        {
            index = uniform_int(generator);
        } while (infected_v[index] != 0);
        infected_q[j] = index;
        timers[j] = uint8_t(std::max(duration(generator), uint32_t(1)));    // At least one day of infectiousness.
        infected_v[index] = 1;
    }

    count = initial_infections;
}

void DisplayState(uint32_t n, uint32_t* infected, uint8_t* timers, uint32_t count)
{
    std::cout << "State: ";
    for (size_t i = 0; i < count; ++i)
    {
        std::cout << infected[i] << '/' << uint32_t(timers[i]) << ' ';
    }
    std::cout << std::endl;
}

void Step(uint32_t population_size, uint32_t connection_count, uint32_t* network_m, float* connect_m,        // RO
    uint32_t* infected_q_in, uint8_t* timers_in, uint32_t infected_size,        // RO
    uint8_t* infected_v, uint8_t* recovered_v,                                  // R/W
    uint32_t& susceptible_c, uint32_t& infectious_c, uint32_t& recovered_c,     // R/W
    uint32_t* infected_q_out, uint8_t* timers_out, uint32_t& infected_count)    // WO
{
    std::mt19937_64 generator(24240101);
    std::uniform_real_distribution<float> draw(0.0f, 1.0f);
    std::poisson_distribution<uint32_t> duration(3);

    infected_count = 0;

    for (uint32_t queue_index = 0; queue_index < infected_size; ++queue_index)
    {
        uint32_t source_individual = infected_q_in[queue_index];
        uint32_t* neighbors = network_m + (connection_count * source_individual);
        float* susceptibility = connect_m + (connection_count * source_individual);
        for (size_t j = 0; j < connection_count; ++j)
        {
            uint32_t target_individual = neighbors[j];
            if (!infected_v[target_individual] && !recovered_v[target_individual])
            {
                if (draw(generator) < susceptibility[j])
                {
                    if (!quiet) { std::cout << "Individual " << source_individual << " infected individual " << target_individual << '.' << std::endl; }

                    --susceptible_c;
                    ++infectious_c;
                    infected_v[target_individual] = 1;

                    infected_q_out[infected_count] = target_individual;
                    timers_out[infected_count] = uint8_t(std::max(duration(generator), uint32_t(1)));
                    ++infected_count;
                }
            }
        }

        uint8_t remaining = timers_in[queue_index] - 1;
        if (remaining > 0)
        {
            infected_q_out[infected_count] = source_individual;
            timers_out[infected_count] = remaining;
            ++infected_count;
        }
        else
        {
            if (!quiet) { std::cout << "Individual " << source_individual << " recovered." << std::endl; }
            infected_v[source_individual] = 0;
            recovered_v[source_individual] = 1;
            --infectious_c;
            ++recovered_c;
        }
    }
}

void DisplayRecovered(uint32_t n, uint8_t* recovered)
{
    std::cout << "Recovered: ";
    for (size_t i = 0; i < n; ++i)
    {
        std::cout << uint32_t(recovered[i]) << ' ';
    }
    std::cout << std::endl;
}

void RadixSort(uint32_t* data, uint32_t count)
{
    uint32_t a[256];
    uint32_t b[256];
    uint32_t c[256];
    uint32_t d[256];

    uint32_t* buffer = new uint32_t[count];

    memset(a, 0, 256 * sizeof(uint32_t));
    memset(b, 0, 256 * sizeof(uint32_t));
    memset(c, 0, 256 * sizeof(uint32_t));
    memset(d, 0, 256 * sizeof(uint32_t));

    for (uint32_t i = 0; i < count; ++i)
    {
        uint32_t value = data[i];
        uint32_t ia = value & 0xFF;
        a[ia] += 1;
        value >>= 8;
        uint32_t ib = value & 0xFF;
        b[ib] += 1;
        value >>= 8;
        uint32_t ic = value & 0xFF;
        c[ic] += 1;
        value >>= 8;
        uint32_t id = value & 0xFF;
        d[id] += 1;
    }

    for (uint32_t i = 1; i < 256; ++i)
    {
        a[i] += a[i - 1];
        b[i] += b[i - 1];
        c[i] += c[i - 1];
        d[i] += d[i - 1];
    }

    for (int32_t i = count - 1; i >= 0; --i)
    {
        uint32_t value = data[i];
        uint32_t ia = value & 0xFF;
        buffer[a[ia] - 1] = value;
        a[ia] -= 1;
    }

    for (int32_t i = count - 1; i >= 0; --i)
    {
        uint32_t value = buffer[i];
        uint32_t ib = (value >> 8) & 0xFF;
        data[b[ib] - 1] = value;
        b[ib] -= 1;
    }

    for (int32_t i = count - 1; i >= 0; --i)
    {
        uint32_t value = data[i];
        uint32_t ic = (value >> 16) & 0xFF;
        buffer[c[ic] - 1] = value;
        c[ic] -= 1;
    }

    for (int32_t i = count - 1; i >= 0; --i)
    {
        uint32_t value = buffer[i];
        uint32_t id = (value >> 24) & 0xFF;
        data[d[id] - 1] = value;
        d[id] -= 1;
    }

    delete[] buffer;

    return;
}

void Unique(uint32_t* data, uint32_t& count)
{
    uint32_t* a = data;
    uint32_t* b = data + 1;
    uint32_t* limit = data + count;

    while ((b < limit) && (*b != *a))
    {
        ++a;
        ++b;
    }

    while (b < limit)
    {
        while ((b < limit) && (*b == *a))
        {
            ++b;
        }
        while ((b < limit) && (*b != *a))
        {
            ++a;
            *a = *b;
            ++b;
        }
    }

    ++a;

    count = uint32_t(a - data);

    return;
}
