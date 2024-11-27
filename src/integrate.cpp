// integrate.cpp
#include "simulation.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <Eigen/Dense> // For matrix operations

#include "energy.h" 
#include "rigidbody.h"
#include "parameters.h"
#include "integrate.h"
#include "meshops.h"

void verlet_integration(SimulationData& sim_data,std::fstream &logfile) {
// main loop
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

  Mesh M1;
  Energy E1;
  //RigidBody body;
  int i=0;
  int toln = 0;
  for (i = 0; i < sim_data.iterations; i++){

    //Initial Integration Step
    M1.mesh_cal(sim_data.V1, sim_data.F1);
    E1.compute_bendingenergy_force(sim_data.V1,sim_data. F1,sim_data. Kb,sim_data. Force_Bending, sim_data.EnergyBending, M1);
    E1.compute_areaenergy_force(sim_data.V1, sim_data.F1, sim_data.Ka, sim_data.area_target, sim_data.Force_Area, sim_data.EnergyArea, M1);
    E1.compute_volumeenergy_force(sim_data.V1, sim_data.F1, sim_data.Kv, sim_data.volume_target, sim_data.Force_Volume, sim_data.EnergyVolume, M1);
    if(sim_data.particle_flag && i%sim_data.bondfrequency==0)igl::signed_distance(sim_data.V1, sim_data.V2, sim_data.F2, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, sim_data.signed_distance, sim_data.facet_index, sim_data.closest_points, sim_data.normals_closest_points);
    if(sim_data.particle_flag && i%sim_data.bondfrequency==0){E1.compute_adhesion_energy_force(sim_data.V1, sim_data.F1, sim_data.closest_points, sim_data.rho, sim_data.U,sim_data.r_equilibrium,sim_data.rc,sim_data.angle_flag,
                                    sim_data.particle_position, sim_data.Ew_t, sim_data.Kw,sim_data.Force_Adhesion,sim_data.Force_Repulsion,sim_data.signed_distance, sim_data.EnergyAdhesion,sim_data.EnergyBias, M1);}
    
    if (sim_data.random_force_flag) {E1.compute_random_force( sim_data.gamma,sim_data. kbT, sim_data.mass, sim_data.dt, sim_data.Force_Random);
                            E1.compute_drag_force(sim_data.velocity,sim_data.gamma, sim_data.mass, sim_data.Force_Drag);}
    //std::cout << "Force Adhesion" << Force_Adhesion << std::endl;
    sim_data.EnergyPotential= sim_data.EnergyBending + sim_data.EnergyArea + sim_data.EnergyVolume + sim_data.EnergyAdhesion + sim_data.EnergyBias;
    sim_data.EnergyKinetic= 0.5 * sim_data.mass * (sim_data.velocity.rowwise().squaredNorm().sum());
    sim_data.EnergyTotal =sim_data. EnergyPotential + sim_data.EnergyKinetic+sim_data.EnergyParticleKinetic;
    sim_data.Force_Total =sim_data. Force_Bending + sim_data.Force_Area + sim_data.Force_Volume + sim_data.Force_Adhesion + sim_data.Force_Random + sim_data.Force_Drag;
    
    //LAMMPS Integration
    sim_data.velocity += (sim_data.Force_Total / sim_data.mass) * 0.5 * sim_data.dt;
    sim_data.V1 += sim_data.velocity * sim_data.dt;

    //old_updater 
    //acceleration_half_step = Force_Total / mass;
    //velocity_half_step = velocity_half_step + 0.5 *dt* (acceleration_half_step - (gamma*velocity));// + Force_Random ;
    //V1 += velocity * dt + 0.5 * acceleration_half_step * (dt * dt);
    

    //  Rigid Body Calculations 
    //ForcesonParticleVertices
    if (sim_data.particle_flag) {
      sim_data.ForcesOnParticle = -sim_data.Force_Repulsion;
      // Calculate the acceleration of the center of mass based on the net force
      sim_data.particle_acceleration_com = (sim_data.ForcesOnParticle.colwise().sum()) / sim_data.total_mass_particle;
      // Update all vertex positions by translating with the velocity
      sim_data.particle_velocity_com += sim_data.particle_acceleration_com * sim_data.dtf;
      if (i % sim_data.logfrequency == 0) std::cout << "Particle Velocity: " << sim_data.particle_velocity_com.transpose() << std::endl;
      sim_data.V2.rowwise() += (sim_data.particle_velocity_com * sim_data.dt).transpose();
  
    
    //  Rigid Body Calculations (Initial Integration Step)
    
    body.calculate_torque(sim_data.ForcesOnParticle, sim_data.closest_points, sim_data.center_of_mass, sim_data.torque); // torque calculation, Tau = r x F
    if (i%sim_data.logfrequency==0)std::cout << "Torque: " << sim_data.torque.transpose() << std::endl;
    //calculate angular momentum
    body.angular_momentum(sim_data.torque, sim_data.dtf ,sim_data.ang_momentum);
    // //calculate angular velocity
    body.calculate_omega(sim_data.ang_momentum, sim_data.rotation_matrix, sim_data.idiag, sim_data.ang_velocity);

    if (i % sim_data.logfrequency == 0) std::cout << "Angular Velocity: " << sim_data.ang_velocity.transpose() << std::endl;

    // Update the quaternion
    body.update_quaternion(sim_data.current_quaternion, sim_data.ang_velocity, sim_data.dt, sim_data.new_quaternion);
    sim_data.current_quaternion=sim_data.new_quaternion;
    body.q_to_exyz(sim_data.new_quaternion, sim_data.rotation_matrix); 
    if (i % sim_data.logfrequency == 0) std::cout << "Rotation Matrix: \n" << sim_data.rotation_matrix << std::endl;
    // // Update the particle vertices based on the rotation matrix
    //v= vcm + omega x r
    //rotate the particle based on rotation_matrix
    body.calculate_center_of_mass(sim_data.V2,sim_data.F2,sim_data.center_of_mass);
    body.rotate_vertices(sim_data.V2,sim_data.center_of_mass,sim_data.displace,sim_data.rotation_matrix); // r_vector_new= Rotation_matrix * r_vector 
    
    //Particle Kinetic Energy
    sim_data.EnergyParticleKineticTranslation = 0.5 * sim_data.total_mass_particle * (sim_data.particle_velocity_com.rowwise().squaredNorm().sum());
    sim_data.EnergyParticleKineticRotation = 0.5 * sim_data.ang_velocity.transpose() * sim_data.idiag * sim_data.ang_velocity;
    sim_data.EnergyParticleKinetic = sim_data.EnergyParticleKineticTranslation + sim_data.EnergyParticleKineticRotation;

    //Rigid Body Calculations End
    }
   

    

   //Final Integration Step

    //Repeat the force calucaltion here//Final Integration Step
    M1.mesh_cal(sim_data.V1, sim_data.F1);
    E1.compute_bendingenergy_force(sim_data.V1, sim_data.F1, sim_data.Kb, sim_data.Force_Bending, sim_data.EnergyBending, M1);
    E1.compute_areaenergy_force(sim_data.V1, sim_data.F1, sim_data.Ka, sim_data.area_target, sim_data.Force_Area, sim_data.EnergyArea, M1);
    E1.compute_volumeenergy_force(sim_data.V1, sim_data.F1, sim_data.Kv, sim_data.volume_target, sim_data.Force_Volume, sim_data.EnergyVolume, M1);
    if(sim_data.particle_flag && i % sim_data.bondfrequency == 0) {
      igl::signed_distance(sim_data.V1, sim_data.V2, sim_data.F2, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, sim_data.signed_distance, sim_data.facet_index, sim_data.closest_points, sim_data.normals_closest_points);
    }


  
    if (sim_data.particle_flag && i % sim_data.bondfrequency == 0) {
          E1.compute_adhesion_energy_force(sim_data.V1, sim_data.F1, sim_data.closest_points, sim_data.rho, sim_data.U,sim_data.r_equilibrium,sim_data.rc,sim_data.angle_flag,
                                    sim_data.particle_position, sim_data.Ew_t, sim_data.Kw,sim_data.Force_Adhesion,sim_data.Force_Repulsion,sim_data.signed_distance, sim_data.EnergyAdhesion,sim_data.EnergyBias, M1);
        }
    if (sim_data.random_force_flag) {
      E1.compute_random_force(sim_data.gamma, sim_data.kbT, sim_data.mass, sim_data.dt, sim_data.Force_Random);
      E1.compute_drag_force(sim_data.velocity, sim_data.gamma, sim_data.mass, sim_data.Force_Drag);
        }
    sim_data.EnergyPotential = sim_data.EnergyBending + sim_data.EnergyArea + sim_data.EnergyVolume + sim_data.EnergyAdhesion + sim_data.EnergyBias;
    sim_data.EnergyKinetic = 0.5 * sim_data.mass * (sim_data.velocity.rowwise().squaredNorm().sum());
    sim_data.EnergyTotal = sim_data.EnergyPotential + sim_data.EnergyKinetic + sim_data.EnergyParticleKinetic;
    sim_data.Force_Total = sim_data.Force_Bending + sim_data.Force_Area + sim_data.Force_Volume + sim_data.Force_Adhesion + sim_data.Force_Random + sim_data.Force_Drag;
    sim_data.force_residual = sim_data.Force_Total.norm();

    //LAMMPS Integration
    sim_data.velocity += (sim_data.Force_Total / sim_data.mass) * 0.5 * sim_data.dt;

    //old_updater 
    //acceleration = Force_Total / mass;
    // // Update velocities with average acceleration
    // //velocity = velocity_half_step + 0.5 *dt*(acceleration - (gamma * velocity_half_step));// + Force_Random ;
    //velocity = 0.5 * (acceleration + acceleration_half_step) * dt;

    //ForcesonParticleVertices
    //if(particle_flag){E1.redistributeAdhesionForce(V2,F2,closest_points, Force_Repulsion, facet_index,ForcesOnVertices); 
    if(sim_data.particle_flag){
    sim_data.ForcesOnParticle = -sim_data.Force_Repulsion;

    sim_data.particle_acceleration_com = (sim_data.ForcesOnParticle.colwise().sum()) / sim_data.total_mass_particle;

    sim_data.particle_velocity_com += sim_data.particle_acceleration_com * sim_data.dtf; // update the particle velocity

    
    //  Rigid Body Calculations (Final Integration Step)
    body.calculate_torque(sim_data.ForcesOnParticle, sim_data.closest_points, sim_data.center_of_mass, sim_data.torque); // torque calculation, Tau = r x F

    body.angular_momentum(sim_data.torque, sim_data.dtf, sim_data.ang_momentum);

    body.calculate_omega(sim_data.ang_momentum, sim_data.rotation_matrix, sim_data.idiag, sim_data.ang_velocity);
    //Rigid Body Calculations End

    }
  
    sim_data.rVol = 6 * sqrt(PI) * M1.volume_total * pow(M1.area_total, -1.5);   

    if (i % sim_data.logfrequency == 0) {
      sim_data.EnergyChangeRate_log = (sim_data.EnergyTotal - sim_data.EnergyTotalold_log) / (sim_data.logfrequency * sim_data.dt);
      sim_data.EnergyTotalold_log = sim_data.EnergyTotal;
      sim_data.etol(toln++) = sim_data.EnergyChangeRate_log;

      // screen output
      if (sim_data.particle_flag)
      std::cout<<i<<"  "<<sim_data.rVol<<"  "<<sim_data.EnergyBending<<"  "<<sim_data.EnergyAdhesion<<"  "<<sim_data.EnergyTotal<<"  "<<sim_data.EnergyChangeRate_log<<"  "<<sim_data.force_residual<<"  "<<std::endl;
      else
      std::cout<<i<<"  "<<sim_data.rVol<<"  "<<sim_data.EnergyBending<<"  "<<sim_data.EnergyTotal<<"  "<<sim_data.EnergyChangeRate_log<<"  "<<sim_data.force_residual<<"  "<<std::endl;
      // logfile output
      logfile<<i<<"  ";
      logfile<<sim_data.time<<"  ";
      logfile<<M1.area_total<<"  ";
      logfile<<M1.volume_total<<"  ";
      logfile<<sim_data.rVol<<"  ";
      logfile<<sim_data.EnergyBending<<"  ";
      logfile<<sim_data.EnergyArea<<"  ";   
      logfile<<sim_data.EnergyVolume<<"  ";
      if (sim_data.particle_flag) {
      logfile<<sim_data.EnergyAdhesion<<"  "; 
      if (sim_data.forced_wrapping_flag) logfile<<sim_data.EnergyBias<<"  ";
      }
      logfile<<sim_data.EnergyPotential<<"  ";
      logfile<<sim_data.EnergyTotal<<"  ";
      logfile<<sim_data.EnergyKinetic<<"  ";
      logfile<<sim_data.EnergyParticleKinetic<<"  ";
      logfile<<sim_data.EnergyChangeRate_log<<"  ";
      logfile<<sim_data.force_residual<<std::endl;
    }

    if (i % sim_data.tolsteps == 0) {
      if (i != 0) {
        sim_data.EnergyChangeRate_avg = sim_data.etol(Eigen::seq(toln-1-sim_data.tolmean_steps, toln-1)).mean();

        if (std::abs(sim_data.EnergyChangeRate_avg) < sim_data.tolerance && sim_data.tolerance_flag) {
          std::cout<<"Energy change rate reaches the threshold."<<std::endl;
          std::cout<<"Simulation reaches equilibrium state."<<std::endl;
          break;
        }
      }
    }

    if (i % sim_data.dumpfrequency == 0) {
      char dumpfilename[128];
      sprintf(dumpfilename, "dump%08d.off", i);
	    igl::writeOFF(dumpfilename, sim_data.V1, sim_data.F1);
      if(sim_data.particle_flag){
      char dumpfilename_p[128];
      sprintf(dumpfilename_p, "particle%08d.off", i);
      igl::writeOFF(dumpfilename_p, sim_data.V2, sim_data.F2);}
	  }

    if (sim_data.v_smooth_flag || sim_data.delaunay_tri_flag) {
      if ((i+1) %sim_data.mesh_reg_frequency == 0) {
        if (sim_data.v_smooth_flag) sim_data.V1 = M1.vertex_smoothing(sim_data.V1, sim_data.F1);
        if (sim_data.delaunay_tri_flag) {
          igl::edge_lengths(sim_data.V1, sim_data.F1, sim_data.l);
          igl::intrinsic_delaunay_triangulation(sim_data.l, sim_data.F1, sim_data.l, sim_data.F1);
        }
      }
    }

    sim_data.time += sim_data.dt;

  

    if ((i+1) == sim_data.iterations) std::cout<<"Simulation reaches max iterations."<<std::endl;


    //main loop ends here

    // // logfile output
    // logfile<<i<<"  ";
    // logfile<<sim_data.time<<"  ";
    // logfile<<M1.area_total<<"  ";
    // logfile<<M1.volume_total<<"  ";
    // logfile<<sim_data.rVol<<"  ";
    // logfile<<sim_data.EnergyBending<<"  ";
    // logfile<<sim_data.EnergyArea<<"  ";   
    // logfile<<sim_data.EnergyVolume<<"  ";
    // if (sim_data.particle_flag) {
    //   logfile<<sim_data.EnergyAdhesion<<"  ";
    //   if (sim_data.forced_wrapping_flag) logfile<<sim_data.EnergyBias<<"  ";
    // }
    // logfile<<sim_data.EnergyTotal<<"  ";
    // logfile<<sim_data.force_residual<<std::endl;
    //logfile<<"Total run time: "<<duration.count()<<" mins"<<std::endl;
    // logfile.close();   
    igl::writeOFF(sim_data.outFile, sim_data.V1, sim_data.F1);
  }

}

void initialize_simulation(SimulationData& sim_data, Parameter& parameter,std::fstream& logfile){

//     Parameter parameter;

//   // initialization of simulaiton parameters
//   readParameter();
  //SimulationData sim_data;
  //SimulationData sim_data;
  igl::readOFF(parameter.meshFile,sim_data.V1,sim_data.F1);
  igl::readOFF(parameter.particleFile,sim_data.V2, sim_data.F2);
  
  sim_data.numF = sim_data.F1.rows();
  sim_data.numV = sim_data.V1.rows();
  sim_data.distance_threshold = 0.1;
  sim_data.initialize(sim_data.numV);
  // screen and log output of simulation settings

//   std::fstream logfile("logfile.txt", std::ios::out); // Open for writing and to append | std::ios::app
//   if (!logfile.is_open()) {
//         std::cerr << "ERROR: cannot access logfile." << std::endl;
//         return 1;
//     }

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

 
  // Initialize the simulation based on the number of vertices

}