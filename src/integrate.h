// integrate.h
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "simulation.h"
#include "parameters.h"

void verlet_integration(SimulationData& sim_data,std::fstream &logfile);


void initialize_simulation(SimulationData& sim_data, Parameter& parameter,std::fstream& logfile);

#endif
// integrate.h

