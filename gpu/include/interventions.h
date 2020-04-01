/*
** copyright banner
*/

#pragma once

#include "flu-main.h"

void load_interventions(const std::string filename);
void process_interventions(uint32_t step, Configuration& config, SimulationState& state);
