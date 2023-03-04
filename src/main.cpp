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
Eigen::MatrixXd V;                                      // matrix storing vertice coordinates
Eigen::MatrixXi F;
Parameter parameter;

int main() {
  // initialization of simulaiton parameters
  readParameter();
  igl::readOFF(parameter.meshFile, V, F);
  numF = F.rows();
  numV = V.rows();

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
  double Rp, u, U, rho, rc, X0, Y0, Z0, Ew_t, Kw;
  int angle_flag;
  if (particle_flag) {
    Rp = parameter.particle_radius;
    u = parameter.adhesion_strength;
    U = (Kb * u) / (Rp * Rp);
    rho =  parameter.potential_range;
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
    if (!parameter.particle_coord_flag) {
      X0 = 0.0, Y0 = 0.0, Z0 = V.col(2).maxCoeff() + parameter.particle_position * (Rp + 1.0*rho);
    }
    else {
      X0 = parameter.X0, Y0 = parameter.Y0, Z0 = parameter.Z0;
    }

    std::cout<<"Particle position: "<<X0<<", "<<Y0<<", "<<Z0<<std::endl;
    std::cout<<"Particle radius: "<<Rp<<std::endl;
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

  Mesh M1;
  Energy E1;

  Eigen::MatrixXd Force_Area(numV, 3), Force_Volume(numV, 3), Force_Bending(numV, 3), Force_Adhesion(numV, 3), velocity(numV, 3), Force_Total(numV, 3); //force components
  velocity.setZero();
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

  // initiate screen output
  if (particle_flag) std::cout<<"Iteration  ReducedVolume  BendingEnergy  AdhesionEnergy  TotalEnergy  EnergyChangeRate  ForceResidual"<<std::endl;
  else std::cout<<"Iteration  ReducedVolume  BendingEnergy  TotalEnergy  EnergyChangeRate  ForceResidual"<<std::endl;

  // main loop
  int i;
  int toln = 0;
  for (i = 0; i < iterations; i++)
  {
    M1.mesh_cal(V, F);
    E1.compute_bendingenergy_force(V, F, Kb, Force_Bending, EnergyBending, M1);
    E1.compute_areaenergy_force(V, F, Ka, area_target, Force_Area, EnergyArea, M1);
    E1.compute_volumeenergy_force(V, F, Kv, volume_target, Force_Volume, EnergyVolume, M1);
    if (particle_flag) E1.compute_adhesion_energy_force(V, F, X0, Y0, Z0, Rp, rho, U, rc, angle_flag, particle_position, Ew_t, Kw, Force_Adhesion, EnergyAdhesion, EnergyBias, M1);

    EnergyTotal = EnergyBending + EnergyArea + EnergyVolume + EnergyAdhesion + EnergyBias;
    Force_Total = Force_Bending + Force_Area + Force_Volume + Force_Adhesion;
    force_residual = Force_Total.norm();

    rVol = 6 * sqrt(PI) * M1.volume_total * pow(M1.area_total, -1.5);

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
	    igl::writeOFF(dumpfilename, V, F);
	  }

    if (i % resfrequency == 0) igl::writeOFF(parameter.resFile, V, F);

    velocity = Force_Total / gamma;
    V += velocity * dt;

    if (v_smooth_flag || delaunay_tri_flag) {
      if ((i+1) % mesh_reg_frequency == 0) {
        if (v_smooth_flag) V = M1.vertex_smoothing(V, F);
        if (delaunay_tri_flag) {
          igl::edge_lengths(V, F, l);
          igl::intrinsic_delaunay_triangulation(l, F, l, F);
        }
      }
    }

    time += dt;

  }

  if ((i+1) == iterations) std::cout<<"Simulation reaches max iterations."<<std::endl;

  auto end = system_clock::now();
  auto duration = duration_cast<minutes>(end - start);

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
      
  igl::writeOFF(parameter.outFile, V, F);

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
  runfile >> parameter.vertex_smoothing_flag;
  getline(runfile, line);
  getline(runfile, line);
  runfile >> parameter.delaunay_triangulation_flag;
  runfile.close();
}
