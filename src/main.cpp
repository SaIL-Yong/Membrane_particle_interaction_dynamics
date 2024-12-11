#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include "meshops.h"
#include "energy.h"
#include "parameters.h"
#include "rigidbody.h"
#include "simulation.h"
#include "integrate.h"
using namespace std::chrono;


Parameter parameter;
int main() {
  // initialization of simulaiton parameters
  readParameter();
  //SimulationData sim_data;
  SimulationData sim_data;
  // Log file setup
  std::fstream logfile("logfile.txt", std::ios::out);
    if (!logfile.is_open()) {
        std::cerr << "ERROR: cannot access logfile." << std::endl;
        return 1;
    }
  // Start timing simulation
  auto start = std::chrono::system_clock::now();
  std::cout << "Simulation Start:\n";

  // Initialize the simulation
  initialize_simulation(sim_data, parameter, logfile);
  //call the integration function
  verlet_integration(sim_data,logfile);
  
  // End timing simulation
  auto end = system_clock::now();
  auto duration = duration_cast<minutes>(end - start);
  logfile<<"Total run time: "<<duration.count()<<" mins"<<std::endl;
  logfile.close();
  return 0;
}
void readParameter()
{
  std::string line;
  std::ifstream runfile;
  runfile.open("run_parameters.txt");
  if (!runfile.is_open()) std::cout<<"ERROR: cannot access simlation parameter file."<<std::endl;
  getline(runfile, line);
  runfile >> parameter.iterations;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.dt;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.Kb;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.Ka;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.Kv;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.reduced_volume;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.tolerance;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.tolerance_flag;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.tolfrequency;  // in units of sim. time (iteration * timestep)
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.gamma;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.mass;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.kbT;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.particle_flag;
  getline(runfile, line);
  if (parameter.particle_flag) {
    getline(runfile, line);
    runfile >> parameter.particle_position;
    getline(runfile, line);
    getline(runfile, line);
    // if (line.compare("particle_coordinate") == 0) {
    //   runfile >> parameter.X0 >> parameter.Y0 >> parameter.Z0;
    //   parameter.particle_coord_flag = 1;
    //   getline(runfile, line);
    //   getline(runfile, line);
    // }
    // else parameter.particle_coord_flag = 0;
    runfile >> parameter.particle_coord_flag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.particle_radius;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.adhesion_strength;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.mass_particle_ratio;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.potential_range;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.r_equilibrium;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.angle_condition_flag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.forced_wrapping_flag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.wrapping_fraction;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.wrapping_bias_strength;
    getline(runfile, line);
  }
  getline(runfile, line);
  getline(runfile, parameter.meshFile);
  getline(runfile, line);
  getline(runfile, parameter.particleFile);
  getline(runfile, line);
  getline(runfile, parameter.outFile);
  getline(runfile, line);
  getline(runfile, parameter.resFile);
  getline(runfile, line);
  runfile >> parameter.logfrequency;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.dumpfrequency;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.resfrequency;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.mesh_reg_frequency;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.bondfrequency;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.vertex_smoothing_flag;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.delaunay_triangulation_flag;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.random_force_flag;
  runfile.close();
}

