/*
** copyright banner
*/

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

struct IIntervention
{
// public:
    IIntervention() {}
    virtual ~IIntervention() {}
    virtual void Apply(uint32_t time_step, Configuration& config, SimulationState& state) = 0;
    virtual uint32_t GetDuration( void ) = 0;
};

struct BaseIntervention : IIntervention
{
//public:
    BaseIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute = "", uint32_t attribute_value = 0);
    uint32_t GetDuration( void ) { return duration; }

protected:
    std::string tag;
    uint32_t duration;
    float factor;
    std::string attribute;
    uint32_t attribute_value;
};

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
    uint16_t* attributes,
    uint16_t target_value,
    CON_TYPE* rx_modifiers,
    float factor)
{
    uint32_t individual_index = (blockDim.x * blockIdx.x) + threadIdx.x;
    if (individual_index < population_size)
    {
        if ( attributes[individual_index] == target_value )
        {
            rx_modifiers[individual_index] *= factor;
        }
    }
}

void attribute_intervention(uint32_t step, SimulationState& state, const std::string& attribute_name, uint16_t target_value, float attenuation)
{
    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = (state.population_size + block.x - 1) / block.x;
    dim3 grid(gridX, 1, 1);

    ATTR_TYPE* attributes(state.attributes.at(attribute_name));

    attribute_intervention_kernel<<<grid, block>>>( state.population_size, attributes, target_value, state.rx_modifiers, attenuation);

    CheckReturn(cudaDeviceSynchronize());
    CheckReturn(cudaGetLastError());
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
            rx_mod[k] = rx_mod[k] * attenuation_factor
            add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

__global__ void neighbor_intervention_kernel(
    uint32_t queue_size,
    IND_TYPE* queue,
    CON_TYPE* connect_m,
    uint32_t connection_count,
    float factor
)
{
    uint32_t queue_index = (blockDim.x * blockIdx.x) + threadIdx.x;
    if ( queue_index < queue_size )
    {
        uint32_t individual_id = queue[queue_index];
        CON_TYPE* connections = connect_m + (connection_count * individual_id);
        for (uint32_t neighbor_index = 0; neighbor_index < connection_count; ++neighbor_index)
        {
            connections[neighbor_index] *= factor;
        }
    }
}

void neighbor_intervention(uint32_t step, Configuration& config, SimulationState& state, float factor)
{
    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = (state.population_size + block.x - 1) / block.x;
    dim3 grid(gridX, 1, 1);

    IND_TYPE* queue = ((step & 1) ? state.infected_individuals_odd_time_steps : state.infected_individuals_even_time_steps);
    neighbor_intervention_kernel<<<grid, block>>>( state.counters->infected_currently, queue, state.connect_m, config.max_connections, factor );

    CheckReturn(cudaDeviceSynchronize());
    CheckReturn(cudaGetLastError());
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
        // rx_mod[k] = rx_mod[k] * attenuation_factor
        connection_weights *= attenuation_factor
        add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

__global__ void node_matching_intervention_kernel(

)
{

}

void node_matching_intervention(uint32_t step, SimulationState& state, uint32_t attribute, uint16_t duration, CON_TYPE factor)
{
    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = max((state.counters->infected_currently + block.x - 1) / block.x, 1);
    dim3 grid(gridX, 1, 1);

    // uint32_t* infected_queue = ((step & 1) ? state.infected_individuals_odd_time_steps : state.infected_individuals_even_time_steps);

    node_matching_intervention_kernel<<<grid, block>>>();
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
        rx_mod[k] = rx_mod[k] * attenuation_factor
        add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

__global__ void neighbor_matching_intervention_kernel(
    uint32_t infected_queue_size,
    uint32_t* infected_queue,
    ATTR_TYPE* attributes,
    IND_TYPE* network_m,
    uint32_t connection_count,
    CON_TYPE* rx_modifiers,
    CON_TYPE factor
)
{
    uint32_t infected_index = (blockDim.x * blockIdx.x) + threadIdx.x;
    if (infected_index < infected_queue_size)
    {
        IND_TYPE source_individual = infected_queue[infected_index];
        ATTR_TYPE attribute = attributes[source_individual];
        IND_TYPE* neighbors = network_m + (connection_count * source_individual);
        for (uint32_t neighbor_index = 0; neighbor_index < connection_count; ++neighbor_index)
        {
            IND_TYPE target_individual = neighbors[neighbor_index];
            if ( attributes[target_individual] == attribute )
            {
                rx_modifiers[target_individual] *= factor;
            }
        }
    }
}

void neighbor_matching_intervention(uint32_t step, SimulationState& state, const std::string& attribute_name, CON_TYPE factor)
{
    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = max((state.counters->infected_currently + block.x - 1) / block.x, 1);
    dim3 grid(gridX, 1, 1);

    uint32_t* infected_queue = ((step & 1) ? state.infected_individuals_odd_time_steps : state.infected_individuals_even_time_steps);
    ATTR_TYPE* attributes(state.attributes.at(attribute_name));

    neighbor_matching_intervention_kernel<<<grid, block>>>(
        state.counters->infected_currently,
        infected_queue,
        attributes,
        state.network_m,
        state.max_connections,
        state.rx_modifiers,
        factor
    );
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
    // TBD
    printf("Applying AttributeIntervention '%s' at time %d\n", tag.c_str(), step);
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
    // TBD
    printf("Applying NeighborIntervention '%s' at time %d\n", tag.c_str(), step);
    neighbor_intervention(step, config, state, factor);
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
    // TBD
    printf("Applying NodeMatchingIntervention '%s' at time %d\n", tag.c_str(), step);
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
    neighbor_matching_intervention(step, state, attribute, factor);
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

void load_interventions(const std::string filename)
{
    io::CSVReader<8, io::trim_chars<' ', '\t'>, io::no_quote_escape<','>, io::throw_on_overflow, io::single_line_comment<'#'>> reader(filename);    // allow comment lines prefixed with '#'
    reader.read_header(io::ignore_extra_column, "tag", "type", "trigger", "trigger_value", "duration", "factor", "attribute", "attribute_value");
    std::string tag;
    uint32_t type;
    std::string trigger;
    uint32_t trigger_value;
    uint32_t duration;
    float factor;
    std::string attribute;
    uint16_t attribute_value;
    bool error = false;
    while(reader.read_row(tag, type, trigger, trigger_value, duration, factor, attribute, attribute_value))
    {
        if ( (trigger != TIME) and (trigger != COUNT) )
        {
            fprintf( stderr, "Invalid trigger '%s' in intervention file '%s'.\n", trigger.c_str(), filename.c_str() );
            error = true;
            continue;
        }

        IIntervention* intervention(nullptr);

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
                error = true;
                break;
        }

        if ( !error )
        {
            AddIntervention( intervention, trigger, trigger_value );
        }
    }

    if ( error )
    {
        exit(1);
    }
}
