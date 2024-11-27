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

  std::fstream logfile("logfile.txt", std::ios::out);
    if (!logfile.is_open()) {
        std::cerr << "ERROR: cannot access logfile." << std::endl;
        return 1;
    }
  // Initialize the simulation based on the number of vertices

    // Start timing simulation
    auto start = std::chrono::system_clock::now();
    std::cout << "Simulation Start:\n";

    // Log file setup

    // Initialize the simulation based on the number of vertices
  initialize_simulation(sim_data, parameter, logfile);

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


    /*
    logfile << "Simulation Start:\n";
    // Dynamic logging depending on flags
    if (sim_data.particle_flag) {
          logfile << "Iteration  Time  Area  Volume  ReducedVolume  BendingEnergy  AreaEnergy  VolumeEnergy  AdhesionEnergy  BiasedWrappingEnergy  PotentialEnergy  TotalEnergy  KineticEnergy  ParticleKineticEnergy  EnergyChangeRate  ForceResidual" << std::endl;
      } else {
          logfile << "Iteration  Time  Area  Volume  ReducedVolume  BendingEnergy  AreaEnergy  VolumeEnergy  PotentialEnergy  TotalEnergy  KineticEnergy  EnergyChangeRate  ForceResidual" << std::endl;
    }
 
    // Initialize the rigid body dynamics
    RigidBody body;
    body.calculate_properties(sim_data.V2, sim_data.mass_particle, sim_data.rotation_matrix, sim_data.idiag, sim_data.displace);
    std::cout << "Principal moments of inertia (Diagonal):\n" << sim_data.idiag << std::endl;
    sim_data.center_of_mass = body.getCenterOfMass();
    std::cout << "Center of Mass: " << sim_data.center_of_mass.transpose() << std::endl;

    //Eigen::Quaterniond current_quaternion;
    body.exyz_to_q(sim_data.rotation_matrix, sim_data.current_quaternion);
    std::cout << "Initial Quaternion (Identity): "
              << "w = " << sim_data.current_quaternion.w() << ", "
              << "x = " << sim_data.current_quaternion.x() << ", "
              << "y = " << sim_data.current_quaternion.y() << ", "
              << "z = " << sim_data.current_quaternion.z() << std::endl;

    // Display initial simulation output
    if (sim_data.particle_flag) {
        std::cout << "Iteration  ReducedVolume  BendingEnergy  AdhesionEnergy  TotalEnergy  EnergyChangeRate ForceResidual" << std::endl;
    } else {
        std::cout << "Iteration  ReducedVolume  BendingEnergy  TotalEnergy  EnergyChangeRate  ForceResidual" << std::endl;
    }

  */
  // Main simulation loop

  /*
  
    igl::readOFF(parameter.meshFile,sim_data.V1,sim_data.F1);
  igl::readOFF(parameter.particleFile,sim_data.V2, sim_data.F2);
  
  sim_data.numF = sim_data.F1.rows();
  sim_data.numV = sim_data.V1.rows();
  sim_data.distance_threshold = 0.1;
  sim_data.initialize(sim_data.numV);
  // screen and log output of simulation settings

  std::fstream logfile("logfile.txt", std::ios::out); // Open for writing and to append | std::ios::app
  if (!logfile.is_open()) {
        std::cerr << "ERROR: cannot access logfile." << std::endl;
        return 1;
    }

    // Write initial log
  logfile << "This is logfile for simulation" << std::endl;

  sim_data.iterations = parameter.iterations;
  sim_data.logfrequency = parameter.logfrequency;
  sim_data.dumpfrequency = parameter.dumpfrequency;
  sim_data.bondfrequency=parameter.bondfrequency;
  sim_data.resfrequency = parameter.resfrequency;
  sim_data.dt = parameter.dt;
  sim_data.time =0.0;
  sim_data.dtf=sim_data.dt/2.0;
  sim_data.tolerance = parameter.tolerance;
  //sim_data.force_residual,force_ratio;
  sim_data.tolerance_flag = parameter.tolerance_flag;
  sim_data.tolfrequency = parameter.tolfrequency;
  sim_data.tolsteps = floor(sim_data.tolfrequency /sim_data.dt);
  sim_data.tolmean_steps = floor(sim_data.tolsteps/sim_data.logfrequency);
  sim_data.etol;
  sim_data.etol.resize(floor(sim_data.iterations/sim_data.logfrequency));
  sim_data.etol.setZero();

  std::cout<<"Mesh info:"<<std::endl;
  std::cout<<"Number of vertices: "<<sim_data.numV<<" Number of faces: "<<sim_data.numF<<"\n"<<std::endl;
  std::cout<<"Max number of iterations: "<<sim_data.iterations<<std::endl;
  std::cout<<"Log output frequency: "<<sim_data.logfrequency<<std::endl;
  std::cout<<"Mesh dump frequency: "<<sim_data.dumpfrequency<<std::endl;
  std::cout<<"Restart save frequency: "<<sim_data.resfrequency<<std::endl;
  std::cout<<"Time step: "<<sim_data.dt<<std::endl;
  logfile<<"Mesh info:"<<std::endl;
  logfile<<"Number of vertices: "<<sim_data.numV<<" Number of faces: "<<sim_data.numF<<"\n"<<std::endl;
  logfile<<"Max number of iterations: "<<sim_data.iterations<<std::endl;
  logfile<<"Log output frequency: "<<sim_data.logfrequency<<std::endl;
  logfile<<"Mesh dump frequency: "<<sim_data.dumpfrequency<<std::endl;
  logfile<<"Restart save frequency: "<<sim_data.resfrequency<<std::endl;
  logfile<<"Time step: "<<sim_data.dt<<std::endl;

  if (sim_data.tolerance_flag) {
    std::cout<<"Convergence: ON, Tolerance: "<<sim_data.tolerance<<std::endl;
    std::cout<<"Tolerance check frequency: "<<sim_data.tolfrequency<<" time units\n"<<std::endl;
    logfile<<"Convergence: ON, Tolerance: "<<sim_data.tolerance<<std::endl;
    logfile<<"Tolerance check frequency: "<<sim_data.tolfrequency<<" time units\n"<<std::endl;
  }
  else {
    std::cout<<"Convergence: OFF\n"<<std::endl;
    logfile<<"Convergence: OFF\n"<<std::endl;
  }
  
  // paraemters for membrane properties
  sim_data.gamma = parameter.gamma;
  sim_data.mass = parameter.mass;
  sim_data.kbT = parameter.kbT;
  sim_data.Kb = parameter.Kb;
  sim_data.Kv = 0.0;
  sim_data.Ka = parameter.Ka;
  sim_data.Rv = 1.0;
  sim_data.area_target = 4*PI*sim_data.Rv*sim_data.Rv;
  sim_data.volume_target = 0.0;
  sim_data.rVol; // true reduced volume
  

  std::cout<<"Vesicle radius: "<<sim_data.Rv<<std::endl;
  std::cout<<"Membrane drag coefficient: "<<sim_data.gamma<<std::endl;
  std::cout<<"Membrane mass coefficient: "<<sim_data.mass<<std::endl;
  std::cout<<"Membrane bending modulus: "<<sim_data.Kb<<std::endl;
  std::cout<<"Membrane stretching modulus: "<<sim_data.Ka<<std::endl;
  logfile<<"Vesicle radius: "<<sim_data.Rv<<std::endl;
  logfile<<"Membrane drag coefficient: "<<sim_data.gamma<<std::endl;
  logfile<<"Membrane mass coefficient: "<<sim_data.mass<<std::endl;
  logfile<<"Membrane bending modulus: "<<sim_data.Kb<<std::endl;
  logfile<<"Membrane stretching modulus: "<<sim_data.Ka<<std::endl;

  if (std::abs(parameter.Kv) > EPS) {
    sim_data.rVol_t = parameter.reduced_volume;
    sim_data.Kv = parameter.Kv;
    sim_data.volume_target = sim_data.rVol_t*(4.0/3.0)*PI*pow(sim_data.Rv,3);
    std::cout<<"Target vesicle reduced volume: "<<sim_data.rVol_t<<std::endl;
    std::cout<<"Vesicle osmotic strength constant: "<<sim_data.Kv<<"\n"<<std::endl;
    logfile<<"Target vesicle reduced volume: "<<sim_data.rVol_t<<std::endl;
    logfile<<"Vesicle osmotic strength constant: "<<sim_data.Kv<<"\n"<<std::endl;
  }
  ///Flags
  sim_data.angle_flag = parameter.angle_condition_flag;
  sim_data.random_force_flag = parameter.random_force_flag;
  sim_data.particle_flag = parameter.particle_flag;
  sim_data.forced_wrapping_flag = parameter.forced_wrapping_flag;

  //file input/oupout
  sim_data.meshFile = parameter.meshFile;
  sim_data.outFile = parameter.outFile;

  if (sim_data.particle_flag) {
        igl::readOFF(parameter.particleFile, sim_data.V2, sim_data.F2);
        sim_data.numVp = sim_data.V2.rows();
        sim_data.numFp = sim_data.F2.rows();
        std::cout << "Number of particle vertices: " << sim_data.numVp << " Number of particle faces: " << sim_data.numFp << "\n";
        
        Mesh M2;
        M2.mesh_cal(sim_data.V2, sim_data.F2);
        sim_data.Rp = sqrt(M2.area_total / (4 * PI));
        sim_data.u = parameter.adhesion_strength;
        sim_data.U = (sim_data.Kb * sim_data.u) / (sim_data.Rp * sim_data.Rp);
        sim_data.rho = parameter.potential_range;
        sim_data.r_equilibrium = parameter.r_equilibrium;
        sim_data.rc = 10.0 * sim_data.rho;
        sim_data.mass_particle = sim_data.mass*parameter.mass_particle_ratio;
        sim_data.total_mass_particle = sim_data.mass_particle * sim_data.V2.rows();
        sim_data.particle_position = parameter.particle_position;

        igl::centroid(sim_data.V2, sim_data.F2, sim_data.COM);
        std::cout << "Particle position: " << sim_data.COM.transpose() << std::endl;
        std::cout << "Particle adhesion strength: " << sim_data.U << std::endl;
        std::cout << "Particle adhesion range: " << sim_data.rho << std::endl;
        std::cout << "Particle adhesion cutoff: " << sim_data.rc << std::endl;
        std::cout << "Distance threshold: " << sim_data.distance_threshold << std::endl;
        std::cout << "Particle mass: " << sim_data.total_mass_particle << std::endl;
        //std::fstream logfile("logfile.txt", std::ios::out);
        
        logfile << "Particle properties and interactions:" << std::endl;
        logfile << "Particle position: " << sim_data.COM.transpose() << std::endl;
        logfile << "Particle radius: " << sim_data.Rp << std::endl;
        logfile << "Particle adhesion strength: " << sim_data.U << std::endl;
        logfile << "Particle adhesion range: " << sim_data.rho << std::endl;
        logfile << "Particle adhesion cutoff: " << sim_data.rc << std::endl;
        logfile << "Distance threshold: " << sim_data.distance_threshold << std::endl;
        if (sim_data.angle_flag) {
                logfile << "Angle criterion: ON" << std::endl;
          } else {
                logfile << "Angle criterion: OFF" << std::endl;
        }    
    


    // parameters for forced wrapping
    if (parameter.forced_wrapping_flag) {
        sim_data.forced_wrapping_fraction = parameter.wrapping_fraction;
        sim_data.Kw = parameter.wrapping_bias_strength;
        double Area_w_t = sim_data.forced_wrapping_fraction * M2.area_total;  // Ensure M2 is defined and area_total is calculated earlier
        sim_data.Ew_t = -sim_data.U * Area_w_t;

        std::cout << "Forced wrapping fraction: " << sim_data.forced_wrapping_fraction << std::endl;
        std::cout << "Forced wrapping strength constant: " << sim_data.Kw << "\n";
        std::cout << "Particle Surface Area: " << M2.area_total << std::endl;
        std::cout << "Target adhesion energy: " << sim_data.Ew_t << std::endl;

        logfile << "Forced wrapping fraction: " << sim_data.forced_wrapping_fraction << std::endl;
        logfile << "Forced wrapping strength constant: " << sim_data.Kw << "\n";
        logfile << "Particle Surface Area: " << M2.area_total << std::endl;
        logfile << "Target adhesion energy: " << sim_data.Ew_t << std::endl;
      }
    }
    // Mesh regularization settings
    sim_data.v_smooth_flag = parameter.vertex_smoothing_flag;
    sim_data.delaunay_tri_flag = parameter.delaunay_triangulation_flag;
    sim_data.mesh_reg_frequency = parameter.mesh_reg_frequency;
    sim_data.bondfrequency = parameter.bondfrequency;
    std::cout << "Vertex smoothing: " << (sim_data.v_smooth_flag ? "ON" : "OFF") << std::endl;
    std::cout << "Delaunay triangulation: " << (sim_data.delaunay_tri_flag ? "ON" : "OFF") << std::endl;
    std::cout << "Mesh regularization frequency: " << sim_data.mesh_reg_frequency << "\n";

    logfile << "Vertex smoothing: " << (sim_data.v_smooth_flag ? "ON" : "OFF") << std::endl;
    logfile << "Delaunay triangulation: " << (sim_data.delaunay_tri_flag ? "ON" : "OFF") << std::endl;
    logfile << "Mesh regularization frequency: " << sim_data.mesh_reg_frequency << "\n";

  */