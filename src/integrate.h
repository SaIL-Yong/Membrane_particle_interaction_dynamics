// integrate.h
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "simulation.h"
#include "parameters.h"

void verlet_integration(SimulationData& sim_data,std::fstream &logfile);

void calculate_forces(SimulationData& sim_data, Mesh& M1, Energy& E1,int current_iteration);

void initialize_simulation(SimulationData& sim_data, Parameter& parameter,std::fstream& logfile);

#endif
// integrate.h

