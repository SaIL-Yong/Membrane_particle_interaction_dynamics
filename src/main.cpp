#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include "meshops.h"
#include "energy.h"
#include "parameters.h"
#include "rigidbody.h"
using namespace std::chrono;
int numV;                                               // number of vertices
int numF;                                               // number of faces
Eigen::MatrixXd V1,V2;                                      // matrix storing vertice coordinates
Eigen::MatrixXi F1,F2;
Parameter parameter;

int main() {
  // initialization of simulaiton parameters
  readParameter();
  igl::readOFF(parameter.meshFile, V1, F1);
  //igl::readOFF(parameter.particleFile, V2, F2);
  
  numF = F1.rows();
  numV = V1.rows();



  Eigen::VectorXi nearest;              // Nearest neighbor of each vertex in V1 in V2
  std::vector<std::pair<int, int>> bonds;

  /**
   * @brief The threshold distance used for determining if two particles are close enough to interact with each other.
   */
  double distance_threshold = 0.1;

  // screen and log output of simulation settings
  std::fstream logfile;
  logfile.open("logfile.txt",std::ios::out);
  if (logfile.is_open())
  {
    logfile<<"This is logfile for simulation"<<std::endl;
  } else {
    std::cout<<"ERROR: cannot access logfile."<<std::endl;
  }

  int iterations = parameter.iterations;
  int logfrequency = parameter.logfrequency;
  int dumpfrequency = parameter.dumpfrequency;
  int bondfrequency=parameter.bondfrequency;
  int resfrequency = parameter.resfrequency;
  double dt = parameter.dt, time = 0.0;
  double tolerance = parameter.tolerance;
  double force_residual,force_ratio;
  int tolerance_flag = parameter.tolerance_flag;
  double tolfrequency = parameter.tolfrequency;
  int tolsteps = floor(tolfrequency / dt);
  int tolmean_steps = floor(tolsteps/logfrequency);
  Eigen::VectorXd etol;
  etol.resize(floor(iterations/logfrequency));
  etol.setZero();

  std::cout<<"Mesh info:"<<std::endl;
  std::cout<<"Number of vertices: "<<numV<<" Number of faces: "<<numF<<"\n"<<std::endl;
  //if(particle_flag)std::cout<<"Number of particle vertices: "<<numVp<<" Number of paticle faces: "<<numFp<<"\n"<<std::endl;
  std::cout<<"Max number of iterations: "<<iterations<<std::endl;
  std::cout<<"Log output frequency: "<<logfrequency<<std::endl;
  std::cout<<"Mesh dump frequency: "<<dumpfrequency<<std::endl;
  std::cout<<"Restart save frequency: "<<resfrequency<<std::endl;
  std::cout<<"Time step: "<<dt<<std::endl;
  logfile<<"Mesh info:"<<std::endl;
  logfile<<"Number of vertices: "<<numV<<" Number of faces: "<<numF<<"\n"<<std::endl;
  logfile<<"Max number of iterations: "<<iterations<<std::endl;
  logfile<<"Log output frequency: "<<logfrequency<<std::endl;
  logfile<<"Mesh dump frequency: "<<dumpfrequency<<std::endl;
  logfile<<"Restart save frequency: "<<resfrequency<<std::endl;
  logfile<<"Time step: "<<dt<<std::endl;

  if (tolerance_flag) {
    std::cout<<"Convergence: ON, Tolerance: "<<tolerance<<std::endl;
    std::cout<<"Tolerance check frequency: "<<tolfrequency<<" time units\n"<<std::endl;
    logfile<<"Convergence: ON, Tolerance: "<<tolerance<<std::endl;
    logfile<<"Tolerance check frequency: "<<tolfrequency<<" time units\n"<<std::endl;
  }
  else {
    std::cout<<"Convergence: OFF\n"<<std::endl;
    logfile<<"Convergence: OFF\n"<<std::endl;
  }
  
  // paraemters for membrane properties]
  double gamma = parameter.gamma;
  double mass = parameter.mass;
  double kbT = parameter.kbT;
  double Kb = parameter.Kb;
  double Kv = 0.0;
  double Ka = parameter.Ka;
  double Rv = 1.0;
  double area_target = 4*PI*Rv*Rv;
  double volume_target = 0.0;
  double rVol; // true reduced volume

  std::cout<<"Vesicle radius: "<<Rv<<std::endl;
  std::cout<<"Membrane drag coefficient: "<<gamma<<std::endl;
  std::cout<<"Membrane mass coefficient: "<<mass<<std::endl;
  std::cout<<"Membrane bending modulus: "<<Kb<<std::endl;
  std::cout<<"Membrane stretching modulus: "<<Ka<<std::endl;
  logfile<<"Vesicle radius: "<<Rv<<std::endl;
  logfile<<"Membrane drag coefficient: "<<gamma<<std::endl;
  logfile<<"Membrane mass coefficient: "<<mass<<std::endl;
  logfile<<"Membrane bending modulus: "<<Kb<<std::endl;
  logfile<<"Membrane stretching modulus: "<<Ka<<std::endl;

  if (std::abs(parameter.Kv) > EPS) {
    double rVol_t = parameter.reduced_volume;
    Kv = parameter.Kv;
    volume_target = rVol_t*(4.0/3.0)*PI*pow(Rv,3);
    std::cout<<"Target vesicle reduced volume: "<<rVol_t<<std::endl;
    std::cout<<"Vesicle osmotic strength constant: "<<Kv<<"\n"<<std::endl;
    logfile<<"Target vesicle reduced volume: "<<rVol_t<<std::endl;
    logfile<<"Vesicle osmotic strength constant: "<<Kv<<"\n"<<std::endl;
  }

  // parameters for particle adhesion
  int particle_flag = parameter.particle_flag;
  int particle_position = parameter.particle_position;
  double Rp, u, U, rho, rc, X0, Y0, Z0, Ew_t, Kw,r_equilibrium,epsilon,sigma;
  int angle_flag;

  // Declaration and initialization of COM
  Eigen::Vector3d COM(0.0, 0.0, 0.0),center_of_mass;
  Eigen::MatrixXd signed_distance;  // Matrix to store signed distances
  Eigen::MatrixXi facet_index;  // Matrix to store facet indices
  Eigen::MatrixXd closest_points;  // Matrix to store closest points
  Eigen::MatrixXd normals_closest_points;  // Matrix to store closest normals (if needed)

  int random_force_flag = parameter.random_force_flag;

  if (particle_flag) {
    igl::readOFF(parameter.particleFile, V2, F2);
    int numFp = F2.rows();
    int numVp = V2.rows();
    std::cout<<"Number of particle vertices: "<<numVp<<" Number of paticle faces: "<<numFp<<"\n"<<std::endl;
    logfile<<"Number of particle vertices: "<<numVp<<" Number of paticle faces: "<<numFp<<"\n"<<std::endl;
    Rp = parameter.particle_radius;
    u = parameter.adhesion_strength;
    U = (Kb * u) / (Rp * Rp);
    rho =  parameter.potential_range;
    r_equilibrium=parameter.r_equilibrium;
    rc = 5.0*rho;
    angle_flag = parameter.angle_condition_flag;
    

    if (parameter.particle_position > 0) {
      std::cout<<"Particle position: outside"<<std::endl;
      logfile<<"Particle position: outside"<<std::endl;
    } 
    else if (parameter.particle_position < 0){
      std::cout<<"Particle position: inside"<<std::endl;
      logfile<<"Particle position: inside"<<std::endl;
    }
    
    // position of the particle
    // if (!parameter.particle_coord_flag) {
    //   X0 = 0.0, Y0 = 0.0, Z0 = V1.col(2).maxCoeff() + parameter.particle_position * (Rp + 1.0*rho);
    // }
    // else {
    //   X0 = parameter.X0, Y0 = parameter.Y0, Z0 = parameter.Z0;
    // }

    if (parameter.particle_coord_flag==0){//to check if the given mesh would be taken or added
    double Z0 = V1.col(2).maxCoeff() + 1.0*(parameter.particle_position* (V2.col(2).maxCoeff()+1.0*rho));
    Eigen::ArrayXd v2_col = V2.col(2).array();
    v2_col += Z0;
    V2.col(2) = v2_col.matrix();
    }
    std::string initialparticle="initialparticle.off";
    igl::writeOFF(initialparticle, V2, F2);   //storing initial particle mesh file
    igl::centroid(V2, F2, COM);



    std::cout<<"Particle position: "<<COM(0)<<", "<<COM(1)<<", "<<COM(2)<<std::endl;
    //std::cout<<"Particle radius: "<<Rp<<std::endl;
    std::cout<<"Particle adhesion strength: "<<U<<std::endl;
    std::cout<<"Particle adhesion range: "<<rho<<std::endl;   
    std::cout<<"Particle adhesion cutoff: "<<rc<<std::endl;
    std::cout<<"distance threshold: "<<distance_threshold<<std::endl;
    logfile<<"Particle position: "<<COM(0)<<", "<<COM(1)<<", "<<COM(2)<<std::endl;
    //logfile<<"Particle radius: "<<Rp<<std::endl;
    logfile<<"Particle adhesion strength: "<<U<<std::endl;
    logfile<<"Particle adhesion range: "<<rho<<std::endl;   
    logfile<<"Particle adhesion cutoff: "<<rc<<std::endl;
    logfile<<"distance threshold: "<<distance_threshold<<std::endl;
    if (angle_flag) {
      std::cout<<"Angle criterion: ON\n"<<std::endl;
      logfile<<"Angle criterion: ON\n"<<std::endl;
    }
    else {
      std::cout<<"Angle criterion: OFF\n"<<std::endl;
      logfile<<"Angle criterion: OFF\n"<<std::endl;
    }

    // parameters for forced wrapping
    Mesh M2;
    M2.mesh_cal(V2, F2);

    Ew_t = 0.0;
    Kw = 0.0;
    if (parameter.forced_wrapping_flag) {
        double chi = parameter.wrapping_fraction;
        Kw = parameter.wrapping_bias_strength;
        double Area_w_t = chi*M2.area_total;
        Ew_t = -U*Area_w_t;
        std::cout<<"Forced wrapping fraction: "<<chi<<std::endl;
        std::cout<<"Forced wrapping strength constant: "<<Kw<<"\n"<<std::endl;
        std::cout<<"Particle Surface Area: "<<M2.area_total<<std::endl;
        std::cout<<"target adhesion energy: "<<Ew_t<<std::endl;
        logfile<<"Forced wrapping fraction: "<<chi<<std::endl;
        logfile<<"Forced wrapping strength constant: "<<Kw<<"\n"<<std::endl;
        logfile<<"Particle Surface Area: "<<M2.area_total<<std::endl;
        logfile<<"target adhesion energy: "<<Ew_t<<std::endl;
    }
  }

  // mesh regularization flag
  int v_smooth_flag = parameter.vertex_smoothing_flag;
  int delaunay_tri_flag = parameter.delaunay_triangulation_flag;
  if (v_smooth_flag) {
    std::cout<<"Vertex smoothing: ON"<<std::endl;
    logfile<<"Vertex smoothing: ON"<<std::endl;
  }
  else {
    std::cout<<"Vertex smoothing: OFF"<<std::endl;
    logfile<<"Vertex smoothing: OFF"<<std::endl;
  }
  if (delaunay_tri_flag) {
    std::cout<<"Delaunay triangulation: ON"<<std::endl;
    logfile<<"Delaunay triangulation: ON"<<std::endl;
  }
  else {
    std::cout<<"Delaunay triangulation: OFF"<<std::endl;
    logfile<<"Delaunay triangulation: OFF"<<std::endl;
  }
  int mesh_reg_frequency = parameter.mesh_reg_frequency;
  if (v_smooth_flag || delaunay_tri_flag) {
    std::cout<<"Mesh regularization frequency: "<<mesh_reg_frequency<<"\n"<<std::endl;
    logfile<<"Mesh regularization frequency: "<<mesh_reg_frequency<<"\n"<<std::endl;
  }
  
  // Calculate the distances between each pair of vertices
  ParticleAdhesion P1;
  //P1.find_pairs(V1, F1, V2, F2, distance_threshold, bonds);
  Mesh M1;

  Energy E1;

  Eigen::MatrixXd Force_Area(numV, 3), Force_Volume(numV, 3), Force_Bending(numV, 3), Force_Adhesion(numV, 3),Force_Random(numV,3),Force_Repulsion(numV,3), velocity(numV, 3), Force_Total(numV, 3),
                  acceleration(numV,3),acceleration_half_step(numV,3),velocity_half_step(numV,3),ForcesOnVertices; //force components

  //Eigen::MatrixXd ForcesOnVertices;
  Force_Adhesion.setZero();
  Force_Random.setZero();
  Force_Repulsion.setZero();
  //ForcesOnVertices.setZero();
  velocity.setZero();
  velocity_half_step.setZero();
  
  Force_Total.setZero();
  //Force_Total_old.setZero();
  acceleration.setZero();
  acceleration_half_step.setZero();
 

  double EnergyVolume = 0.0, EnergyArea = 0.0, EnergyBending = 0.0, EnergyAdhesion = 0.0,  EnergyBias = 0.0,
         EnergyTotal = 0.0, EnergyTotalold_log = 0.0, EnergyChangeRate_log = 0.0, EnergyChangeRate_avg = 0.0;  //energy components
  Eigen::MatrixXd l;
  std::cout<<"Simulation Start:\n"<<std::endl;
  logfile<<"Simulation Start:\n"<<std::endl;
  auto start = system_clock::now();

  // initiate logfile output
  if (particle_flag) {
    if (parameter.forced_wrapping_flag)
      logfile<<"Iteration  Time  Area  Volume  ReducedVolume  BendingEnergy  AreaEnergy  VolumeEnergy  AdhesionEnergy  BiasedWrappingEnergy  TotalEnergy  EnergyChangeRate  ForceResidual"<<std::endl;
    else
      logfile<<"Iteration  Time  Area  Volume  ReducedVolume  BendingEnergy  AreaEnergy  VolumeEnergy  AdhesionEnergy  TotalEnergy  EnergyChangeRate  ForceResidual"<<std::endl;
  } else {
    logfile<<"Iteration  Time  Area  Volume  ReducedVolume  BendingEnergy  AreaEnergy  VolumeEnergy  TotalEnergy  EnergyChangeRate  ForceResidual"<<std::endl;
  }
  

  
  //RigidBodyCalculations and particl motions
  //if (particle_flag){
  RigidBody body;
  Eigen::Matrix3d rotation_matrix,idiag;
  Eigen::Matrix3d moment_of_inertia;
  Eigen::Matrix3d inverse_moment_of_inertia;
  Eigen::Vector3d torque,ang_momentum,ang_velocity,particle_velocity_com,particle_acceleration_com;
  Eigen::MatrixXd particle_velocities(V2.rows(),3);
  body.calculate_properties(V2,mass);///space_frame
  // Access and use the calculated properties
  std::cout << "Center of Mass: " << body.getCenterOfMass().transpose() << std::endl;
  std::cout << "Moment of Inertia: \n" << body.getMomentOfInertia() << std::endl;
  std::cout << "Inverse Moment of Inertia (if calculated): \n" << body.getInverseMomentOfInertia() << std::endl;
  body.diagonalize_inertia_tensor(body.getMomentOfInertia(), rotation_matrix,idiag);
  //std::cout << "Rotation Matrix: \n" << rotation_matrix << std::endl;
  Eigen::Quaterniond current_quaternion; // Initialize the quaternion
  body.exyz_to_q(rotation_matrix,current_quaternion);
  Eigen::Quaterniond new_quaternion;//= Eigen::Quaterniond::Identity();
  std::cout << "Initial Quaternion (Identity): "
               << "w = " <<current_quaternion.w() << ", "
               << "x = " <<current_quaternion.x() << ", "
               << "y = " <<current_quaternion.y() << ", "
               << "z = " <<current_quaternion.z() << std::endl;

  // Calculate the angular velocity
 // }
  
///

  // initiate screen output
  if (particle_flag) std::cout<<"Iteration  ReducedVolume  BendingEnergy  AdhesionEnergy  TotalEnergy  EnergyChangeRate ForceResidual"<<std::endl;
  else std::cout<<"Iteration  ReducedVolume  BendingEnergy  TotalEnergy  EnergyChangeRate  ForceResidual"<<std::endl;
  //P1.find_pairs(V1, F1, V2, F2, distance_threshold, bonds);
  //std::cout << "bond is updated." << std::endl;
  // main loop
  int i;
  int toln = 0;
  for (i = 0; i < iterations; i++)
  {
    M1.mesh_cal(V1, F1);
    E1.compute_bendingenergy_force(V1, F1, Kb, Force_Bending, EnergyBending, M1);
    E1.compute_areaenergy_force(V1, F1, Ka, area_target, Force_Area, EnergyArea, M1);
    E1.compute_volumeenergy_force(V1, F1, Kv, volume_target, Force_Volume, EnergyVolume, M1);
    if(particle_flag && i%bondfrequency==0)igl::signed_distance(V1, V2, F2, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, signed_distance, facet_index, closest_points, normals_closest_points);
    if(particle_flag){E1.compute_adhesion_energy_force(V1, F1, closest_points, rho, U,r_equilibrium,rc,angle_flag,
                                    particle_position, Ew_t, Kw,Force_Adhesion,Force_Repulsion,signed_distance, EnergyAdhesion,EnergyBias, M1);}
    if (random_force_flag)E1.compute_random_force(V1, gamma, kbT, mass, dt, Force_Random);
    //std::cout << "Force Adhesion" << Force_Adhesion << std::endl;
    EnergyTotal = EnergyBending + EnergyArea + EnergyVolume + EnergyAdhesion + EnergyBias;
    Force_Total = Force_Bending + Force_Area + Force_Volume + Force_Adhesion+ Force_Random;

    //acceleration = Force_Total/mass;
    acceleration_half_step = Force_Total / mass;
    //velocity_half_step = velocity_half_step + 0.5 *dt* (acceleration_half_step - (gamma*velocity));// + Force_Random ;

    V1 += velocity * dt + 0.5 * acceleration_half_step * (dt * dt);

    //V1 += velocity * dt + 0.5 * acceleration_half_step * (dt * dt);
        //ForcesonParticleVertices
    if(particle_flag){E1.redistributeAdhesionForce(V2,F2,closest_points, Force_Repulsion, facet_index,ForcesOnVertices); 
    /*std::ofstream file_force("particle_force.txt");
    if (file_force.is_open()) {
    file_force<< ForcesOnVertices<< std::endl;
    file_force.close();
    std::cout << "particle force successfully saved to file." << std::endl;
    }
    else {
    std::cout << "Error: cannot open particle force file." << std::endl;
    }*/
    } 

    //  Rigid Body Calculations 
    // body.calculate_center_of_mass(V2,F2,center_of_mass);
    // body.calculate_torque(ForcesOnVertices, V2, center_of_mass, torque); //torque calculation
    //     // Calculate the acceleration of the center of mass based on the net force
    // particle_acceleration_com =  ForcesOnVertices.colwise().sum() / V2.rows();

    // // Update the velocity of the center of mass based on the acceleration
    // particle_velocity_com = particle_acceleration_com * dt;
    //     //std::cout << "Particle Velocity: " << particle_velocity.transpose() << std::endl;

    // // Update all vertex positions by translating with the velocity
    // //V2.rowwise() += (particle_velocity_com * dt).transpose();
    // //calculate angular momentum
    // body.angular_momentum(torque, dt ,ang_momentum);
    // //calculate angular velocity
    // body.calculate_omega(ang_momentum, rotation_matrix, idiag, ang_velocity);
    // std::cout << "Angular Momentum: " << ang_momentum.transpose() << std::endl;
    // std::cout << "Angular Velocity: " << ang_velocity.transpose() << std::endl;

    // // Update the quaternion
    // body.update_quaternion(current_quaternion, ang_velocity, dt,new_quaternion); 
    // body.q_to_exyz(new_quaternion, rotation_matrix);
    // std::cout << "Rotation Matrix: \n" << rotation_matrix << std::endl;

    // //V= vcm + omega x r
    // body.update_vertex_velocities_positions(V2, center_of_mass,particle_velocity_com ,ang_velocity, dt,particle_velocities);
    // //std::cout << "Particle Velocity: " << particle_velocity.transpose() << std::endl;
    
    
    //Rigid Body Calculations End



    //Repeat the force calucaltion here
    M1.mesh_cal(V1, F1);
    E1.compute_bendingenergy_force(V1, F1, Kb, Force_Bending, EnergyBending, M1);
    E1.compute_areaenergy_force(V1, F1, Ka, area_target, Force_Area, EnergyArea, M1);
    E1.compute_volumeenergy_force(V1, F1, Kv, volume_target, Force_Volume, EnergyVolume, M1);
    if(particle_flag && i%bondfrequency==0){igl::signed_distance(V1, V2, F2, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, signed_distance, facet_index, closest_points, normals_closest_points);

   //std::cout << "Facet Index" << facet_index<<std::endl;
    std::ofstream outfile("signed_distance.txt");
   // Check if the file was successfully opened
    if (outfile.is_open()) {
      outfile << signed_distance << std::endl;
      outfile << facet_index << std::endl;
      outfile.close();
      } 
      else {
      std::cout << "Error: cannot open adhesion force file." <<std::endl;
      }
    }
    if(particle_flag){E1.compute_adhesion_energy_force(V1, F1, closest_points, rho, U,r_equilibrium,rc,angle_flag,
                                    particle_position, Ew_t, Kw,Force_Adhesion,Force_Repulsion,signed_distance, EnergyAdhesion,EnergyBias, M1);}
    if (random_force_flag)E1.compute_random_force(V1, gamma, kbT, mass, dt, Force_Random);
    EnergyTotal = EnergyBending + EnergyArea + EnergyVolume + EnergyAdhesion + EnergyBias;
    Force_Total = Force_Bending + Force_Area + Force_Volume + Force_Adhesion+Force_Random;
    force_residual = Force_Total.norm();
    if (random_force_flag){force_ratio=Force_Random.norm()/(Force_Bending + Force_Area + Force_Volume + Force_Adhesion).norm();
    std::cout<<"Force Ratio: "<<force_ratio<<std::endl;}

    acceleration = Force_Total / mass;
    // Update velocities with average acceleration
    //velocity = velocity_half_step + 0.5 *dt*(acceleration - (gamma * velocity_half_step));// + Force_Random ;
    velocity = 0.5 * (acceleration + acceleration_half_step) * dt;

    rVol = 6 * sqrt(PI) * M1.volume_total * pow(M1.area_total, -1.5);

    //ForcesonParticleVertices
    if(particle_flag){E1.redistributeAdhesionForce(V2,F2,closest_points, Force_Repulsion, facet_index,ForcesOnVertices); 
    /*std::ofstream file_force("particle_force.txt");
    if (file_force.is_open()) {
    file_force<< ForcesOnVertices<< std::endl;
    file_force.close();
    std::cout << "particle force successfully saved to file." << std::endl;
    }
    else {
    std::cout << "Error: cannot open particle force file." << std::endl;
    }*/
    } 

    //  Rigid Body Calculations 
    body.calculate_center_of_mass(V2,F2,center_of_mass);
    body.calculate_torque(ForcesOnVertices, V2, center_of_mass, torque); //torque calculation
        // Calculate the acceleration of the center of mass based on the net force
    particle_acceleration_com =  ForcesOnVertices.colwise().sum() / V2.rows();

    // Update the velocity of the center of mass based on the acceleration
    particle_velocity_com = particle_acceleration_com * dt;
        //std::cout << "Particle Velocity: " << particle_velocity.transpose() << std::endl;

    // Update all vertex positions by translating with the velocity
    //V2.rowwise() += (particle_velocity_com * dt).transpose();
    //calculate angular momentum
    body.angular_momentum(torque, dt ,ang_momentum);
    //calculate angular velocity
    body.calculate_omega(ang_momentum, rotation_matrix, idiag, ang_velocity);
    std::cout << "Angular Momentum: " << ang_momentum.transpose() << std::endl;
    std::cout << "Angular Velocity: " << ang_velocity.transpose() << std::endl;

    // Update the quaternion
    body.update_quaternion(current_quaternion, ang_velocity, dt,new_quaternion); ///simple euler update
    current_quaternion=new_quaternion;
    body.q_to_exyz(new_quaternion, rotation_matrix);
    std::cout << "Rotation Matrix: \n" << rotation_matrix << std::endl;

    //v= vcm + omega x r
    body.update_vertex_velocities_positions(V2,rotation_matrix ,center_of_mass,particle_velocity_com ,ang_velocity, dt,particle_velocities);
    //std::cout << "Particle Velocity: " << particle_velocity.transpose() << std::endl;
    
    
    //Rigid Body Calculations End
    
    
    if (i % logfrequency == 0) {
      EnergyChangeRate_log = (EnergyTotal - EnergyTotalold_log) / (logfrequency * dt);
      EnergyTotalold_log = EnergyTotal;
      etol(toln++) = EnergyChangeRate_log;

      // screen output
      if (particle_flag)
        std::cout<<i<<"  "<<rVol<<"  "<<EnergyBending<<"  "<<EnergyAdhesion<<"  "<<EnergyTotal<<"  "<<EnergyChangeRate_log<<"  "<<force_residual<<"  "<<std::endl;
      else
        std::cout<<i<<"  "<<rVol<<"  "<<EnergyBending<<"  "<<EnergyTotal<<"  "<<EnergyChangeRate_log<<"  "<<force_residual<<"  "<<std::endl;
      // logfile output
      logfile<<i<<"  ";
      logfile<<time<<"  ";
      logfile<<M1.area_total<<"  ";
      logfile<<M1.volume_total<<"  ";
      logfile<<rVol<<"  ";
      logfile<<EnergyBending<<"  ";
      logfile<<EnergyArea<<"  ";   
      logfile<<EnergyVolume<<"  ";
      if (particle_flag) {
        logfile<<EnergyAdhesion<<"  "; 
        if (parameter.forced_wrapping_flag) logfile<<EnergyBias<<"  ";
      }
      logfile<<EnergyTotal<<"  ";
      logfile<<EnergyChangeRate_log<<"  ";
      logfile<<force_ratio<<"  ";
      logfile<<force_residual<<std::endl;
    }

    if (i % tolsteps == 0) {
      if (i != 0) {
        EnergyChangeRate_avg = etol(Eigen::seq(toln-1-tolmean_steps,toln-1)).mean();

        if (std::abs(EnergyChangeRate_avg) < tolerance && tolerance_flag) {
          std::cout<<"Energy change rate reaches the threshold."<<std::endl;
          std::cout<<"Simulation reaches equilibrium state."<<std::endl;
          break;
        }
      }
    }

    if (i % dumpfrequency == 0) {
      char dumpfilename[128];
      sprintf(dumpfilename, "dump%08d.off", i);
	    igl::writeOFF(dumpfilename, V1, F1);
      char dumpfilename_p[128];
      sprintf(dumpfilename_p, "particle%08d.off", i);
      igl::writeOFF(dumpfilename_p, V2, F2);
	  }

    if (i % resfrequency == 0) igl::writeOFF(parameter.resFile, V1, F1);



    //velocity = Force_Total / gamma;
    //V1 += velocity * dt;

    if (v_smooth_flag || delaunay_tri_flag) {
      if ((i+1) % mesh_reg_frequency == 0) {
        if (v_smooth_flag) V1 = M1.vertex_smoothing(V1, F1);
        if (delaunay_tri_flag) {
          igl::edge_lengths(V1, F1, l);
          igl::intrinsic_delaunay_triangulation(l, F1, l, F1);
        }
      }
    }

    time += dt;

  }

  if ((i+1) == iterations) std::cout<<"Simulation reaches max iterations."<<std::endl;

  auto end = system_clock::now();
  auto duration = duration_cast<minutes>(end - start);
  //main loop ends here

  // logfile output
  logfile<<i<<"  ";
  logfile<<time<<"  ";
  logfile<<M1.area_total<<"  ";
  logfile<<M1.volume_total<<"  ";
  logfile<<rVol<<"  ";
  logfile<<EnergyBending<<"  ";
  logfile<<EnergyArea<<"  ";   
  logfile<<EnergyVolume<<"  ";
  if (particle_flag) {
    logfile<<EnergyAdhesion<<"  ";
    if (parameter.forced_wrapping_flag) logfile<<EnergyBias<<"  ";
  }
  logfile<<EnergyTotal<<"  ";
  logfile<<force_residual<<std::endl;
  logfile<<"Total run time: "<<duration.count()<<" mins"<<std::endl;
  logfile.close();
      
  igl::writeOFF(parameter.outFile, V1, F1);

  //Storing force components to text file after equilibrium
  std::ofstream file1("Adhesion_Force.txt");
  // Check if the file was successfully opened
  if (file1.is_open()) {
    file1 << Force_Adhesion << std::endl;
    file1.close();
    std::cout << "Adhesion force successfully saved to file." << std::endl;
  }
  else {
    std::cout << "Error: cannot open adhesion force file." <<std::endl;
  }
  std::ofstream file2("Bending_Force.txt");
  if (file2.is_open()) {
    file2 << Force_Bending << std::endl;
    file2.close();
    std::cout << "Bending force successfully saved to file." << std::endl;
  }
  else {
    std::cout << "Error: cannot open bending force file." << std::endl;
  }
  std::ofstream file3("Area_Force.txt");
  if (file3.is_open()) {
    file3<< Force_Area << std::endl;
    file3.close();
    std::cout << "Area force successfully saved to file." << std::endl;
  }
  else {
    std::cout << "Error: cannot open area force file." << std::endl;
  }
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