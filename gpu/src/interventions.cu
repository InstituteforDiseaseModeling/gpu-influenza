/*
** copyright banner
*/

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <map>
#include <vector>

#include "block-size.h"
#include "check-return.h"
#include "flu-main.h"
#include "interventions.h"

#define CSV_IO_NO_THREAD
#include "csv.h"    // https://github.com/ben-strasser/fast-cpp-csv-parser

BaseIntervention::BaseIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute, uint32_t attribute_value)
    : tag(tag)
    , duration(duration)
    , factor(factor)
    , attribute(attribute)
    , attribute_value(attribute_value)
{
    // This space intentionally left blank.
}

static std::map<uint32_t, std::vector<IIntervention*>> scheduled_interventions;
static std::vector<std::pair<uint32_t, IIntervention*>> triggered_interventions;

__global__ void age_interventions(
    uint32_t population_size,
    uint8_t* active_interventions,
    TIMER_TYPE* intervention_timers,
    uint32_t max_interventions,
    CON_TYPE* intervention_modifiers,
    CON_TYPE* tx_modifiers
)
{
    uint32_t individual_id = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( individual_id < population_size )
    {
        uint32_t active_ivs = uint32_t(active_interventions[individual_id]);                // active interventions for this node
        TIMER_TYPE* timers = intervention_timers + (individual_id * max_interventions);     // pointer to timers for this node's interventions
        CON_TYPE* modifiers = intervention_modifiers + (individual_id * max_interventions); // pointer to modifiers for this node's interventions
        CON_TYPE tx_modifier = CON_TYPE(1.0f);                                              // current tx_modifier, starts at 1.0
        uint8_t next_iv = 0;                                                                // next intervention slot (index)

        for ( uint32_t iv_index = 0; iv_index < active_ivs; ++iv_index )                    // iterate over active interventions
        {
            uint8_t timer = timers[iv_index] - 1;                                           // get this intervention's timer and decrement
            if ( timer > 0 )                                                                // if still active...
            {
                CON_TYPE modifier = modifiers[iv_index];                                    // get its modifier
                tx_modifier *= modifier;                                                    // update current tx_modifier
                timers[next_iv] = timer;                                                    // put remaining timer into next intervention slot (index)
                modifiers[next_iv] = modifier;                                              // copy this interventions modifier to next intervention slot (index)
                ++next_iv;                                                                  // increment index to next intervention slot (index)
            }
        }

        if ( next_iv != active_ivs )                                                        // if we no longer have as many active interventions...
        {
            active_interventions[individual_id] = uint8_t(next_iv);                         // update the active interventions counter for this node
            tx_modifiers[individual_id] = tx_modifier;                                      // update the node's tx_modifier
        }
    }
}

void process_interventions(uint32_t step, Configuration& config, SimulationState& state)
{
    uint32_t infectious_count = state.counters->infected_currently;
    printf("Processing interventions for time step %d (%d currently infected/infectious).\n", step, infectious_count);

    /*
    [age interventions] for each intervention j: 
        Read a[j], the attenuation for this intervention
        Decrement remaining_time for each active intervention and node
        For each active intervention on each node k that reaches 0:
            rx_mod[k] = rx_mod[k] / (1 – a[j])
    */

    dim3 block(BLOCK_SIZE, 1, 1);
    dim3 grid((state.population_size + block.x - 1) / block.x, 1, 1);

    age_interventions<<<grid, block>>>(
        state.population_size,
        state.active_interventions,
        state.intervention_timers,
        state.max_interventions,
        state.intervention_modifiers,
        state.tx_modifiers
    );

    /*
    [trigger new interventions] for each intervention:
        if ( ( trigger_type == “time” AND t == trigger_value ) ||
             ( trigger_type == “infectious_count” AND size(infectious) >= trigger_value ) )
            switch target_type:
                case “neighbor”:  neighbor_intervention(); break
                case “attribute”:  attribute_intervention(); break
                case “node_matching_attribute”: node_matching_intervention(); break
                case “neighbor_matching_attribute”: neighbor_matching_intervention(); break
    */

    auto it = scheduled_interventions.find(step);
    if (it != scheduled_interventions.end())
    {
        for (auto intervention : (*it).second)
        {
            intervention->Apply(step, config, state);
        }
    }

    for (size_t i = 0; i < triggered_interventions.size(); /**/)
    {
        auto& entry = triggered_interventions[i];
        if ( infectious_count >= entry.first )
        {
            IIntervention* intervention = entry.second;
            intervention->Apply(step, config, state);

            // Replace the current intervention with the last one in the vector...
            triggered_interventions[i] = triggered_interventions[triggered_interventions.size()-1];
            triggered_interventions.pop_back(); // Pop/delete the last intervention off the vector

            // Do not increment i, we (may have) just moved the last intervention to this slot.
        }
        else
        {
            ++i;
        }
    }
}

/*
attribute_intervention(): the intervention event is applied to all nodes that match an attribute value.
    read target_attribute, target_value, duration, and attenuation_factor
    for each node k in network:
        if attribute[k][target_attribute] == target_value:
            rx_mod[k] = rx_mod[k] * attenuation_factor
            add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

__global__ void attribute_intervention_kernel(
    uint32_t population_size,
    ATTR_TYPE* attributes,
    ATTR_TYPE target_value,
    CON_TYPE* tx_modifiers,
    float factor,
    uint8_t* active_interventions,
    TIMER_TYPE* intervention_timers,
    uint32_t max_interventions,
    uint8_t duration,
    CON_TYPE* intervention_effects
)
{
    uint32_t individual_id = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( individual_id < population_size )
    {
        if ( attributes[individual_id] == target_value )
        {
            tx_modifiers[individual_id] *= factor;

            uint32_t next_intervention = active_interventions[individual_id]++;
            TIMER_TYPE* timers = intervention_timers + (individual_id * max_interventions);
            timers[next_intervention] = duration;
            CON_TYPE* effects = intervention_effects + (individual_id * max_interventions);
            effects[next_intervention] = factor;
        }
    }
}

/*
neighbor_intervention(): The Intervention event is applied to all neighbors of each infectious node.
    read duration, and attenuation_factor
    initialize “implement_intervention” vector of size N with all zeros.
    for each infectious node k in network:
        get neighbors of node k 
        for each neighbor j:
            implement_intervention[j] = 1
    for each node k in implement_intervention:
        if implement_intervention[k] == 1 then:
            tx_mod[k] = tx_mod[k] * attenuation_factor
            add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

// Individuals might be neighbors of more than one infected individual, we identify them here so the intervention is only applied once.
// Actual tx_modifier change is applied separately.

__global__ void identify_neighbors_kernel(
    uint32_t infected_count,
    IND_TYPE* infected,
    IND_TYPE* network_m,
    uint32_t connection_count,
    uint8_t* intervene
)
{
    // Assumes that state.intervene has been zeroed.

    uint32_t infected_index = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( infected_index < infected_count )
    {
        uint32_t individual_id = infected[infected_index];
        IND_TYPE* neighbors = network_m + (connection_count * individual_id);
        for (uint32_t neighbor_index = 0; neighbor_index < connection_count; ++neighbor_index)
        {
            uint32_t neighbor_id = neighbors[neighbor_index];
            intervene[neighbor_id] = 1;
        }
    }
}

__global__ void intervention_kernel(
    uint32_t population_size,
    uint8_t* intervene,
    CON_TYPE* tx_modifiers,
    float factor,
    uint8_t* active_interventions,
    TIMER_TYPE* intervention_timers,
    uint32_t max_interventions,
    uint8_t duration,
    CON_TYPE* intervention_effects
)
{
    // Assumes that state.intervene has been set with identify_neighbors_kernel (above).

    uint32_t individual_id = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( individual_id < population_size )
    {
        if ( intervene[individual_id] != 0 )
        {
            tx_modifiers[individual_id] *= factor;

            uint32_t next_intervention = active_interventions[individual_id]++;
            TIMER_TYPE* timers = intervention_timers + (individual_id * max_interventions);
            timers[next_intervention] = duration;
            CON_TYPE* effects = intervention_effects + (individual_id * max_interventions);
            effects[next_intervention] = factor;
        }
    }
}

/*
node_matching_intervention(): In this case, the intervention event is applied to all nodes whose
attribute value matches the value of the same attribute in the infectious node.

read target_attribute, duration, and attenuation_factor
initialize an empty “attribute_values” vector.
for each infectious node k in network:
    if attribute[k][target_attribute] is not in attribute_values, then:
        add attribute[target_attribute] to attribute_values
for each node k in network:
    if attribute[k][target_attribute] is in attribute_values:
        // tx_mod[k] = tx_mod[k] * attenuation_factor
        connection_weights *= attenuation_factor
        add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

__global__ void attribute_collection_kernel(uint32_t infected_count, IND_TYPE* infected, ATTR_TYPE* attributes, uint8_t* scratch)
{
    uint32_t infected_index = (blockDim.x * blockIdx.x) + threadIdx.x;
    if ( infected_index < infected_count )
    {
        IND_TYPE individual_id = infected[infected_index];
        ATTR_TYPE attribute = attributes[individual_id];
        scratch[attribute] = 1;
    }
}

__global__ void node_matching_intervention_kernel(
    uint32_t population_size,
    ATTR_TYPE* attributes,
    uint8_t* scratch,
    CON_TYPE* tx_modifiers,
    CON_TYPE factor,
    uint8_t* active_interventions,
    TIMER_TYPE* intervention_timers,
    uint32_t max_interventions,
    uint8_t duration,
    CON_TYPE* intervention_effects
)
{
    uint32_t individual_id = (blockDim.x * blockIdx.x) + threadIdx.x;
    if ( individual_id < population_size )
    {
        ATTR_TYPE attribute = attributes[individual_id];
        if ( scratch[attribute] != 0 )
        {
            tx_modifiers[individual_id] *= factor;

            uint32_t next_intervention = active_interventions[individual_id]++;
            TIMER_TYPE* timers = intervention_timers + (individual_id * max_interventions);
            timers[next_intervention] = duration;
            CON_TYPE* effects = intervention_effects + (individual_id * max_interventions);
            effects[next_intervention] = factor;
        }
    }
}

/*
neighbor_matching_intervention(): In this case, the intervention event is applied to all neighboring
nodes whose attribute value matches the value of the same attribute in the infectious node.

read target_attribute, target_value, duration, and attenuation_factor
initialize “implement_intervention” vector of size N with all zeros.
for each infectious node k in network:
    get neighbors of node k 
    for each neighbor j:
        if attribute[j][target_attribute] == attribute[k][target_attribute]:
            implement_intervention[j] = 1
for each node k in implement_intervention:
    if implement_intervention[k] == 1 then:
        tx_mod[k] = tx_mod[k] * attenuation_factor
        add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

// Individuals might be neighbors of more than one infected individual, we identify them here so the intervention is only applied once.
// Actual tx_modifier change is applied separately.

__global__ void identify_matching_neighbors_kernel(
    uint32_t infected_queue_size,
    uint32_t* infected_queue,
    ATTR_TYPE* attributes,
    IND_TYPE* network_m,
    uint32_t connection_count,
    uint8_t* intervene
)
{
    // Assumes that state.intervene has been zeroed.

    uint32_t infected_index = (blockDim.x * blockIdx.x) + threadIdx.x;

    if ( infected_index < infected_queue_size )
    {
        IND_TYPE source_individual_id = infected_queue[infected_index];
        ATTR_TYPE attribute = attributes[source_individual_id];
        IND_TYPE* neighbors = network_m + (connection_count * source_individual_id);

        for (uint32_t neighbor_index = 0; neighbor_index < connection_count; ++neighbor_index)
        {
            IND_TYPE neighbor_id = neighbors[neighbor_index];
            if ( attributes[neighbor_id] == attribute )
            {
                intervene[neighbor_id] = 1;
            }
        }
    }
}

class AttributeIntervention : public BaseIntervention
{
public:
    AttributeIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute = "", uint32_t attribute_value = 0);
    virtual ~AttributeIntervention() {}

    virtual void Apply(uint32_t step, Configuration& config, SimulationState& state) override;
};

AttributeIntervention::AttributeIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute, uint32_t attribute_value)
    : BaseIntervention(tag, duration, factor, attribute, attribute_value)
{
    // This space intentionally left blank.
}

void AttributeIntervention::Apply(uint32_t step, Configuration& config, SimulationState& state)
{
    printf("Applying AttributeIntervention '%s' at time %d\n", tag.c_str(), step);

    dim3 block(BLOCK_SIZE, 1, 1);
    dim3 grid((state.population_size + block.x - 1) / block.x, 1, 1);

    ATTR_TYPE* attributes(state.attributes.at(attribute));

    attribute_intervention_kernel<<<grid, block>>>(
        state.population_size,
        attributes,
        attribute_value,
        state.tx_modifiers,
        factor,
        state.active_interventions,
        state.intervention_timers,
        state.max_interventions,
        duration,
        state.intervention_modifiers
    );

    CheckReturn(cudaDeviceSynchronize());
    CheckReturn(cudaGetLastError());
}

class NeighborIntervention : public BaseIntervention
{
public:
    NeighborIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute = "", uint32_t attribute_value = 0);
    virtual ~NeighborIntervention() {}

    virtual void Apply(uint32_t step, Configuration& config, SimulationState& state) override;
};

NeighborIntervention::NeighborIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute, uint32_t attribute_value)
    : BaseIntervention(tag, duration, factor, attribute, attribute_value)
{
    // This space intentionally left blank.
}

void NeighborIntervention::Apply(uint32_t step, Configuration& config, SimulationState& state)
{
    printf("Applying NeighborIntervention '%s' at time %d\n", tag.c_str(), step);

    CheckReturn(cudaMemset( state.intervene, 0, state.population_size * sizeof(uint8_t) ));

    IND_TYPE* infected = ((step & 1) ? state.infected_individuals_odd_time_steps : state.infected_individuals_even_time_steps);
    dim3 block(BLOCK_SIZE, 1, 1);
    dim3 grid((state.counters->infected_currently + block.x - 1) / block.x, 1, 1);  // iterate over currently infectious nodes
    identify_neighbors_kernel<<<grid, block>>>( state.counters->infected_currently, infected, state.network_m, state.max_connections, state.intervene );

    grid = dim3((state.population_size + block.x - 1) / block.x, 1, 1); // iterate over all nodes
    intervention_kernel<<<grid, block>>>(
        state.population_size,
        state.intervene,
        state.tx_modifiers,
        factor,
        state.active_interventions,
        state.intervention_timers,
        state.max_interventions,
        duration,
        state.intervention_modifiers
    );

    CheckReturn(cudaDeviceSynchronize());
    CheckReturn(cudaGetLastError());
}

class NodeMatchingIntervention : public BaseIntervention
{
public:
    NodeMatchingIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute = "", uint32_t attribute_value = 0);
    virtual ~NodeMatchingIntervention() {}

    virtual void Apply(uint32_t step, Configuration& config, SimulationState& state) override;
};

NodeMatchingIntervention::NodeMatchingIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute, uint32_t attribute_value)
    : BaseIntervention(tag, duration, factor, attribute, attribute_value)
{
    // This space intentionally left blank.
}

void NodeMatchingIntervention::Apply(uint32_t step, Configuration& config, SimulationState& state)
{
    printf("Applying NodeMatchingIntervention '%s' (match '%s', attenuate %f) at time %d\n", tag.c_str(), attribute.c_str(), float(factor), step);

    // 1. zero out scratch space
    // 2. mark scratch[attribute_value] = 1 for given attribute from infected individuals
    // 3. apply intervention to all individuals with same attribute_value (i.e., scratch[attribute_value] != 0)

    CheckReturn(cudaMemset(state.scratch, uint8_t(0), (1 << (sizeof(ATTR_TYPE)*8)) * sizeof(uint8_t)));

    dim3 block(BLOCK_SIZE, 1, 1);
    dim3 grid(max((state.counters->infected_currently + block.x - 1) / block.x, 1), 1, 1);  // the first kernel just inspects infectious individuals

    IND_TYPE* infected = ((step & 1) ? state.infected_individuals_odd_time_steps : state.infected_individuals_even_time_steps);
    ATTR_TYPE* attributes(state.attributes.at( attribute ));
    attribute_collection_kernel<<<grid, block>>>(state.counters->infected_currently, infected, attributes, state.scratch);

    grid = dim3(max((state.population_size + block.x - 1) / block.x, 1), 1, 1); // now consider everyone for matching attributes
    node_matching_intervention_kernel<<<grid, block>>>(
        state.population_size,
        attributes,
        state.scratch,
        state.tx_modifiers,
        factor,
        state.active_interventions,
        state.intervention_timers,
        state.max_interventions,
        duration,
        state.intervention_modifiers
    );

    CheckReturn(cudaDeviceSynchronize());
    CheckReturn(cudaGetLastError());
}

class NeighborMatchingIntervention : public BaseIntervention
{
public:
    NeighborMatchingIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute = "", uint32_t attribute_value = 0);
    virtual ~NeighborMatchingIntervention() {}

    virtual void Apply(uint32_t step, Configuration& config, SimulationState& state) override;
};

NeighborMatchingIntervention::NeighborMatchingIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute, uint32_t attribute_value)
    : BaseIntervention(tag, duration, factor, attribute, attribute_value)
{
    // This space intentionally left blank.
}

void NeighborMatchingIntervention::Apply(uint32_t step, Configuration& config, SimulationState& state)
{
    printf("Applying NeighborMatchingIntervention '%s' at time %d\n", tag.c_str(), step);

    // 1. zero out intervene vector (one entry per node)
    // 2. identify nodes who are neighbors of an infected node and have a matching attribute value, set intervene[node_id] = 1
    // 3. apply intervention to all identified nodes (i.e., intervene[node_id] != 0)

    CheckReturn(cudaMemset( state.intervene, 0, state.population_size * sizeof(uint8_t) ));

    uint32_t* infected_queue = ((step & 1) ? state.infected_individuals_odd_time_steps : state.infected_individuals_even_time_steps);
    ATTR_TYPE* attributes(state.attributes.at(attribute));

    dim3 block(BLOCK_SIZE, 1, 1);
    dim3 grid(max((state.counters->infected_currently + block.x - 1) / block.x, 1), 1, 1);  // the first kernel just inspects infectious individuals

    identify_matching_neighbors_kernel<<<grid, block>>>(
        state.counters->infected_currently,
        infected_queue,
        attributes,
        state.network_m,
        state.max_connections,
        state.intervene
    );

    grid = dim3(max((state.population_size + block.x - 1) / block.x, 1), 1, 1); // now consider everyone for matching attributes
    intervention_kernel<<<grid, block>>>(
        state.population_size,
        state.intervene,
        state.tx_modifiers,
        factor,
        state.active_interventions,
        state.intervention_timers,
        state.max_interventions,
        duration,
        state.intervention_modifiers
    );

    CheckReturn(cudaDeviceSynchronize());
    CheckReturn(cudaGetLastError());
}

enum
{
    ATTRIBUTE        = 1,
    NEIGHBOR         = 2,
    NODEMATCHING     = 3,
    NEIGHBORMATCHING = 4
};

#define TIME    ("TIME")
#define COUNT   ("COUNT")

bool ValidateIntervention( const std::string& type, const std::string& trigger, float factor, const std::map<std::string, uint32_t>& map, const std::string& filename );
void AddIntervention( IIntervention* intervention, std::string& trigger, uint32_t trigger_value );

void load_interventions(const std::string& filename)
{
    if ( filename.size() == 0 )
    {
        printf( "No campaign specified.\n" );
        return;
    }

    static std::map<std::string, uint32_t> enum_map =
        {
            { "1", ATTRIBUTE },        { "ATTRIBUTE", ATTRIBUTE },               { "ATTRIBUTES", ATTRIBUTE },          { "ATTR", ATTRIBUTE },
            { "2", NEIGHBOR },         { "NEIGHBOR", NEIGHBOR },                 { "NEIGHBORS", NEIGHBOR },
            { "3", NODEMATCHING },     { "NODEMATCHING", NODEMATCHING },         { "NODEMATCH", NODEMATCHING },
            { "4", NEIGHBORMATCHING }, { "NEIGHBORMATCHING", NEIGHBORMATCHING }, { "NEIGHBORMATCH", NEIGHBORMATCHING }
        };

    io::CSVReader<8, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>, io::throw_on_overflow, io::single_line_comment<'#'>> reader(filename);    // allow comment lines prefixed with '#'

    reader.read_header(io::ignore_extra_column, "tag", "type", "trigger", "trigger_value", "duration", "factor", "attribute", "attribute_value");

    std::string tag;
    std::string type_string;
    std::string trigger;
    uint32_t trigger_value;
    uint32_t duration;
    float factor;
    std::string attribute;
    uint16_t attribute_value;
    bool valid = true;

    while( reader.read_row( tag, type_string, trigger, trigger_value, duration, factor, attribute, attribute_value ) )
    {
        transform(type_string.begin(), type_string.end(), type_string.begin(), ::toupper);  // convert to uppercase to be case-insensitive
        transform(trigger.begin(), trigger.end(), trigger.begin(), ::toupper);              // convert to uppercase to be case-insensitive
        valid &= ValidateIntervention( type_string, trigger, factor, enum_map, filename );
        if ( !valid ) continue;

        IIntervention* intervention(nullptr);

        uint32_t type = enum_map.at( type_string );

        switch (type)
        {
            case ATTRIBUTE:
                printf("Parsed attribute intervention (%s, %s, %d, %d, %f)\n", tag.c_str(), trigger.c_str(), trigger_value, duration, factor);
                intervention = new AttributeIntervention(tag, duration, factor, attribute, attribute_value);
                break;

            case NEIGHBOR:
                printf("Parsed neighbor intervention (%s, %s, %d, %d, %f)\n", tag.c_str(), trigger.c_str(), trigger_value, duration, factor);
                intervention = new NeighborIntervention(tag, duration, factor, attribute, attribute_value);
                break;

            case NODEMATCHING:
                printf("Parsed node matching intervention (%s, %s, %d, %d, %f)\n", tag.c_str(), trigger.c_str(), trigger_value, duration, factor);
                intervention = new NodeMatchingIntervention(tag, duration, factor, attribute, attribute_value);
                break;

            case NEIGHBORMATCHING:
                printf("Parsed neighbor matching intervention (%s, %s, %d, %d, %f)\n", tag.c_str(), trigger.c_str(), trigger_value, duration, factor);
                intervention = new NeighborMatchingIntervention(tag, duration, factor, attribute, attribute_value);
                break;

            default:
                fprintf(stderr, "Unknown intervention type '%d' in campaign file '%s'.\n", type, filename.c_str());
                valid = false;
                break;
        }

        if ( valid )
        {
            AddIntervention( intervention, trigger, trigger_value );
        }
    }

    if ( !valid ) exit( 1 );
}

bool ValidateIntervention( const std::string& type, const std::string& trigger, float factor, const std::map<std::string, uint32_t>& map, const std::string& filename )
{
    bool valid = true;

    if ( (trigger != TIME) and (trigger != COUNT) )
    {
        fprintf( stderr, "Invalid trigger '%s' in intervention file '%s'.\n", trigger.c_str(), filename.c_str() );
        valid = false;
    }

    if ( map.find( type ) == map.end() )
    {
        fprintf( stderr, "Invalid intervention type '%s' in intervention file '%s'.\n", type.c_str(), filename.c_str() );
        valid = false;
    }

    if ( factor < 0.0f )
    {
        fprintf( stderr, "Invalid intervention factor %f in intervention file '%s'.\n", factor, filename.c_str() );
        valid = false;
    }

    return valid;
}

void AddIntervention( IIntervention* intervention, std::string& trigger, uint32_t trigger_value )
{
    if ( trigger == TIME )
    {
        scheduled_interventions[trigger_value].push_back( intervention );
    }
    else if ( trigger == COUNT )
    {
        triggered_interventions.push_back( std::pair<uint32_t, IIntervention*>( trigger_value, intervention ) );
    }
    else
    {
        fprintf( stderr, "Invalid trigger type, '%s', passed to %s.\n", trigger.c_str(), __FUNCTION__ );
        exit( 1 );
    }
}
