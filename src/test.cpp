#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include "meshops.h"
#include "energy.h"
#include "parameters.h"
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
  igl::readOFF(parameter.particleFile, V2, F2);
  
  numF = F1.rows();
  numV = V1.rows();
  int numFp = F2.rows();
  int numVp = V2.rows();

  Eigen::VectorXi nearest;              // Nearest neighbor of each vertex in V1 in V2
  std::vector<std::pair<int, int>> bonds;
  double distance_threshold = 0.20;

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
  double force_residual;
  int tolerance_flag = parameter.tolerance_flag;
  double tolfrequency = parameter.tolfrequency;
  int tolsteps = floor(tolfrequency / dt);
  int tolmean_steps = floor(tolsteps/logfrequency);
  Eigen::VectorXd etol;
  etol.resize(floor(iterations/logfrequency));
  etol.setZero();

  std::cout<<"Mesh info:"<<std::endl;
  std::cout<<"Number of vertices: "<<numV<<" Number of faces: "<<numF<<"\n"<<std::endl;
  std::cout<<"Number of particle vertices: "<<numVp<<" Number of paticle faces: "<<numFp<<"\n"<<std::endl;
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
  
  // paraemters for membrane properties
  double gamma = parameter.gamma;
  double Kb = parameter.Kb;
  double Kv = 0.0;
  double Ka = parameter.Ka;
  double Rv = 1.0;
  double area_target = 4*PI*Rv*Rv;
  double volume_target = 0.0;
  double rVol; // true reduced volume

  std::cout<<"Vesicle radius: "<<Rv<<std::endl;
  std::cout<<"Membrane drag coefficient: "<<gamma<<std::endl;
  std::cout<<"Membrane bending modulus: "<<Kb<<std::endl;
  std::cout<<"Membrane stretching modulus: "<<Ka<<std::endl;
  logfile<<"Vesicle radius: "<<Rv<<std::endl;
  logfile<<"Membrane drag coefficient: "<<gamma<<std::endl;
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
  if (particle_flag) {
    Rp = parameter.particle_radius;
    u = parameter.adhesion_strength;
    U = (Kb * u) / (Rp * Rp);
    rho =  parameter.potential_range;
    r_equilibrium=parameter.r_equilibrium;
    epsilon=parameter.epsilon;
    sigma=parameter.sigma;
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
    double Z0 = V1.col(2).maxCoeff() + (parameter.particle_position* V2.col(2).maxCoeff());
    Eigen::ArrayXd v2_col = V2.col(2).array();
    v2_col += Z0;
    V2.col(2) = v2_col.matrix();
    std::string initialparticle="initialparticle.off";
    igl::writeOFF(initialparticle, V2, F2);   //storing initial particle mesh file


    std::cout<<"Particle adhesion strength: "<<U<<std::endl;
    std::cout<<"Particle adhesion range: "<<rho<<std::endl;   
    std::cout<<"Particle adhesion cutoff: "<<rc<<std::endl;
    logfile<<"Particle position: "<<X0<<", "<<Y0<<", "<<Z0<<std::endl;
    logfile<<"Particle radius: "<<Rp<<std::endl;
    logfile<<"Particle adhesion strength: "<<U<<std::endl;
    logfile<<"Particle adhesion range: "<<rho<<std::endl;   
    logfile<<"Particle adhesion cutoff: "<<rc<<std::endl;
    if (angle_flag) {
      std::cout<<"Angle criterion: ON\n"<<std::endl;
      logfile<<"Angle criterion: ON\n"<<std::endl;
    }
    else {
      std::cout<<"Angle criterion: OFF\n"<<std::endl;
      logfile<<"Angle criterion: OFF\n"<<std::endl;
    }

    // parameters for forced wrapping
    Ew_t = 0.0;
    Kw = 0.0;
    if (parameter.forced_wrapping_flag) {
        double chi = parameter.wrapping_fraction;
        Kw = parameter.wrapping_bias_strength;
        double Area_w_t = chi*4.0*PI*Rp*Rp;
        Ew_t = -U*Area_w_t;
        std::cout<<"Forced wrapping fraction: "<<chi<<std::endl;
        std::cout<<"Forced wrapping strength constant: "<<Kw<<"\n"<<std::endl;
        logfile<<"Forced wrapping fraction: "<<chi<<std::endl;
        logfile<<"Forced wrapping strength constant: "<<Kw<<"\n"<<std::endl;
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

  Eigen::MatrixXd Force_Area(numV, 3), Force_Volume(numV, 3), Force_Bending(numV, 3), Force_Adhesion(numV, 3), velocity(numV, 3), Force_Total(numV, 3); //force components
  velocity.setZero();
  double EnergyVolume = 0.0, EnergyArea = 0.0, EnergyBending = 0.0, EnergyAdhesion = 0.0,  EnergyBias = 0.0,
         EnergyTotal = 0.0, EnergyTotalold_log = 0.0, EnergyChangeRate_log = 0.0, EnergyChangeRate_avg = 0.0;  //energy components
  Eigen::MatrixXd l;
  P1.find_pairs(V1, F1, V2, F2, distance_threshold, bonds);

  M1.mesh_cal(V1, F1);
  E1.compute_bendingenergy_force(V1, F1, Kb, Force_Bending, EnergyBending, M1);
  E1.compute_areaenergy_force(V1, F1, Ka, area_target, Force_Area, EnergyArea, M1);
  E1.compute_volumeenergy_force(V1, F1, Kv, volume_target, Force_Volume, EnergyVolume, M1);
  E1.compute_adhesion_energy_force(V1, F1, V2, F2, rho, U,r_equilibrium,epsilon,sigma, Force_Adhesion, bonds, EnergyAdhesion, M1);
  std::cout<<"w/o lennard jone potential force: "<<Force_Adhesion<<std::endl;

        

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
  runfile >> parameter.particle_flag;
  getline(runfile, line);
  if (parameter.particle_flag) {
    getline(runfile, line);
    runfile >> parameter.particle_position;
    getline(runfile, line);
    getline(runfile, line);
    if (line.compare("particle_coordinate") == 0) {
      runfile >> parameter.X0 >> parameter.Y0 >> parameter.Z0;
      parameter.particle_coord_flag = 1;
      getline(runfile, line);
      getline(runfile, line);
    }
    else parameter.particle_coord_flag = 0;
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
    runfile >> parameter.epsilon;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.sigma;
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
  runfile.close();
}
