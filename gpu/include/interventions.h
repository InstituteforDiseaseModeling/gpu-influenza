/*
** copyright banner
*/

#pragma once

#include "flu-main.h"

void load_interventions(const std::string& filename);
void process_interventions(uint32_t step, Configuration& config, SimulationState& state);

struct IIntervention
{
// public:
    IIntervention() {}
    virtual ~IIntervention() {}
    virtual void Apply(uint32_t time_step, Configuration& config, SimulationState& state) = 0;
};

struct BaseIntervention : IIntervention
{
// public:
    BaseIntervention(const std::string& tag, uint32_t duration, float factor, const std::string& attribute = "", uint32_t attribute_value = 0);

protected:
    std::string tag;
    uint32_t duration;
    float factor;
    std::string attribute;
    uint32_t attribute_value;
};
