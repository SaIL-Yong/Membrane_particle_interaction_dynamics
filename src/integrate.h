#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <fstream>
#include "simulation.h"
#include "meshops.h"
#include "energy.h"

// ------------------------------------------------------
// Initialize everything: read OFF meshes, set SimulationData fields, etc.
// ------------------------------------------------------
void initialize_simulation(
    SimulationData& sim_data,
    const Parameter& parameter,
    std::fstream& logfile
);

// ------------------------------------------------------
// Verlet‐integration loop (advances positions+velocities over time).
// ------------------------------------------------------
void verlet_integration(
    SimulationData& sim_data,
    std::fstream& logfile
);

// ------------------------------------------------------
// Calculate all forces at the current iteration.
//   M1, E1 → fields for mesh 1;  M2, E2 → fields for mesh 2.
//   “current_iteration” is the loop index.
// ------------------------------------------------------
void calculate_forces(
    SimulationData& sim_data,
    Mesh&           M1,
    Energy&         E1,
    Mesh&           M2,
    Energy&         E2,
    int             current_iteration
);

#endif
