#include <iostream>
#include <fstream>
#include "meshops.h"
#include "energy.h"
#include "parameters.h"
#include "simulation.h"
#include "integrate.h"

int main() {
    // 1) Read parameters from run_parameters.txt
    Parameter parameter;
    parameter.readParameter();  // no argument version

    // 2) Prepare SimulationData
    SimulationData sim_data;

    // 3) Open logfile
    std::fstream logfile("logfile.txt", std::ios::out);
    if (!logfile.is_open()) {
        std::cerr << "ERROR: cannot open logfile.txt\n";
        return 1;
    }

    // 4) Initialize simulation (load meshes, set constants, etc.)
    initialize_simulation(sim_data, parameter, logfile);

    // 5) Run Verlet integration
    verlet_integration(sim_data, logfile);

    logfile.close();
    return 0;
}
