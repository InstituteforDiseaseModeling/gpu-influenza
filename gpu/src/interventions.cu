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

#include "csv.h"

struct IIntervention
{
// public:
    virtual ~IIntervention() {}
    virtual void Apply(uint32_t time_step) = 0;
};

struct IScheduledIntervention : public IIntervention
{
    IScheduledIntervention(uint32_t time) : scheduled_time( time ) {}
    virtual uint32_t GetScheduledTime() { return scheduled_time; }

protected:
    uint32_t scheduled_time;
};

struct ITriggeredIntervention : public IIntervention
{
    ITriggeredIntervention(uint32_t threshold) : trigger_value( threshold ) {}
    virtual uint32_t GetTriggerValue() { return trigger_value; }

protected:
    uint32_t trigger_value;
};

class DummyScheduled : public IScheduledIntervention
{
public:
    DummyScheduled(uint32_t time) : IScheduledIntervention(time) {};
    virtual ~DummyScheduled() {}

    virtual void Apply(uint32_t time_step) override;
};

void DummyScheduled::Apply(uint32_t time_step)
{
    printf("Applying scheduled intervention at time %d (scheduled for %d).\n", time_step, scheduled_time);
}

class DummyTriggered : public ITriggeredIntervention
{
public:
    DummyTriggered(uint32_t threshold) : ITriggeredIntervention(threshold) {};
    virtual ~DummyTriggered() {}

    virtual void Apply(uint32_t time_step) override;
};

void DummyTriggered::Apply(uint32_t time_step)
{
    printf("Applying triggered intervention at time %d (threshold %d).\n", time_step, trigger_value);
}

std::map<uint32_t, std::vector<IScheduledIntervention*>> scheduled_interventions;
std::vector<ITriggeredIntervention*> triggered_interventions;

// Testing/debugging
static bool initialized = false;

static void initialize(void)
{
    scheduled_interventions[13].push_back(new DummyScheduled(13));
    scheduled_interventions[42].push_back(new DummyScheduled(42));
    triggered_interventions.push_back(new DummyTriggered(31));
    triggered_interventions.push_back(new DummyTriggered(314));
    initialized = true;

    return;
}

void load_interventsion(std::string& filename)
{

}

void process_interventions(uint32_t step, uint32_t infectious_count)
{
    // Testing/debugging
    if (!initialized) initialize();

    printf("Processing interventions for time step %d (%d currently infected/infectious).\n", step, infectious_count);

    /*
    [age interventions] for each intervention j: 
        Read a[j], the attenuation for this intervention
        Decrement remaining_time for each active intervention and node
        For each active intervention on each node k that reaches 0:
            beta_modifier[k] = beta_modifier[k] / (1 – a[j])
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
            intervention->Apply(step);
        }
    }

    for (size_t i = 0; i < triggered_interventions.size(); /**/)
    {
        auto intervention = triggered_interventions[i];
        if ( infectious_count >= intervention->GetTriggerValue() )
        {
            intervention->Apply(step);
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
            beta_modifier[k] = beta_modifier[k] * attenuation_factor
            add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

__global__ void attribute_intervention_kernel(
    uint32_t population_size,
    uint16_t* attributes,
    uint16_t target_value,
    CON_TYPE* beta_modifier,
    float factor)
{
    uint32_t individual_index = (blockDim.x * blockIdx.x) + threadIdx.x;
    if (individual_index < population_size)
    {
        if ( attributes[individual_index] == target_value )
        {
            beta_modifier[individual_index] *= factor;
        }
    }
}

void attribute_intervention(uint32_t step, SimulationState& state, uint32_t attribute, uint16_t target_value, float attenuation)
{
    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = (state.population_size + block.x - 1) / block.x;
    dim3 grid(gridX, 1, 1);

    uint16_t* attributes(nullptr);
    switch (attribute)
    {
        default:
            break;
    }

    attribute_intervention_kernel<<<grid, block>>>( state.population_size, attributes, target_value, state.beta_factor, attenuation);

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
            beta_modifier[k] = beta_modifier[k] * attenuation_factor
            add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

__global__ void neighbor_intervention_kernel(
    uint32_t population_size,
    IND_TYPE* network_m,
    uint32_t connection_count,
    CON_TYPE* beta_modifier,
    float factor
)
{
    uint32_t source_individual = (blockDim.x * blockIdx.x) + threadIdx.x;
    if ( source_individual < population_size )
    {
        IND_TYPE* neighbors = network_m + (connection_count * source_individual);
        for (uint32_t neighbor_index = 0; neighbor_index < connection_count; ++neighbor_index)
        {
            IND_TYPE target_individual = neighbors[neighbor_index];
            beta_modifier[target_individual] *= factor;
        }
    }
}

void neighbor_intervention(uint32_t step, SimulationState& state, float factor)
{
    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = (state.population_size + block.x - 1) / block.x;
    dim3 grid(gridX, 1, 1);

    neighbor_intervention_kernel<<<grid, block>>>( state.population_size, state.network_m, state.max_connections, state.beta_factor, factor);

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
        beta_modifier[k] = beta_modifier[k] * attenuation_factor
        add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

__global__ void node_matching_intervention_kernel(

)
{

}

void node_matching_intervention(uint32_t step, SimulationState& state, uint32_t attribute, uint16_t duration, CON_TYPE factor)
{
    node_matching_intervention<<<grid, block>>>();
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
        beta_modifier[k] = beta_modifier[k] * attenuation_factor
        add (duration, attenuation_factor) to beta_modifier_tracker[k] vector.
*/

__global__ void neighbor_matching_intervention_kernel(
    uint32_t infected_queue_size,
    uint32_t* infected_queue,
    ATTR_TYPE* attributes,
    IND_TYPE* network_m,
    uint32_t connection_count,
    CON_TYPE* beta_modifiers,
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
                beta_modifiers[target_individual] *= factor;
            }
        }
    }
}

void neighbor_matching_intervention(uint32_t step, SimulationState& state, CON_TYPE factor)
{
    uint32_t blockX = BLOCK_SIZE;
    dim3 block(blockX, 1, 1);
    uint32_t gridX = max((state.counters->infected_currently + block.x - 1) / block.x, 1);
    dim3 grid(gridX, 1, 1);

    uint32_t* infected_queue = ((step & 1) ? state.queue_a_inf : state.queue_b_inf);
    ATTR_TYPE* attributes(nullptr);

    neighbor_matching_intervention_kernel<<<grid, block>>>(
        state.counters->infected_currently,
        infected_queue,
        attributes,
        state.network_m,
        state.max_connections,
        state.beta_factor,
        factor
    );
}