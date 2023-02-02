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
double area_avg;                                        // average area of each triangle mesh
Eigen::MatrixXd V,V_new;                                      // matrix storing vertice coordinates
Eigen::MatrixXi F;
Parameter parameter;

int main(){
    readParameter();
    std::string filename = parameter.meshFile;
    std::string resfilename = parameter.resFile;
    std::string outfilename = parameter.outFile;
    igl::readOFF(filename, V, F);
    numF = F.rows();
    numV = V.rows();
    int iterations=parameter.iterations;
    int logfrequency=100;
    Eigen::MatrixXd Force_Area,Force_Volume,Force_Bending,Force_Adhesion,velocity,ForceTotal; //forces

    double dt=parameter.dt,time=0.0; //time step
    double gamma=parameter.gamma;
    double tolerance=parameter.tolerance;
    int tolerance_flag=parameter.tolerance_flag;
    float reduced_volume=parameter.reduced_volume;
    double EnergyVolume,EnergyArea,EnergyBending,EnergyTotal,EnergyAdhesion,EnergyTotal_new, EnergyChange,E_bias;  //energies
//parameters for particle adhesion
    float bending_modulus=parameter.Kb;
    float rp=parameter.particle_radious;
    float u=parameter.adhesion_strength;
    float rho=(parameter.potential_range)*rp;
    float U=(bending_modulus*u)/(pow(rp,2)) ;
    float rc=5*rho;
    float X=V.maxCoeff()-rp-rho*1,Y=0.0,Z=0.0;
    float chi=parameter.wrapping_fraction;
    float Area_w_t=chi*4.0*PI*(pow(rp,2));
    float Ew_t=-U*Area_w_t;
    float K_bias=10.0;

//parameters for particle adhesion

    Mesh M1;
    Energy E1;
    E1.compute_bendingenergy_force(V,F,Force_Bending,EnergyBending);
    E1.compute_areaenergy_force(V,F,Force_Area,EnergyArea);
    E1.compute_volumeenergy_force(V,F,reduced_volume,Force_Volume,EnergyVolume);
    E1.compute_adhesion_energy_force(V,F,X,Y,Z,rp,rho,u,U,rc,Ew_t,K_bias,Force_Adhesion,EnergyAdhesion,E_bias);
    EnergyTotal=EnergyBending + EnergyArea + EnergyVolume + EnergyAdhesion+E_bias;

    ForceTotal=Force_Bending+Force_Area+Force_Volume+Force_Adhesion;

    std::cout<<"Bending Energy: "<< EnergyBending<<std::endl;
    std::cout<<"\n Total Energy: "<< EnergyTotal<<std::endl;
    std::cout<<"\n reduced_volume: "<< reduced_volume<<std::endl;
    std::cout<<"\n U: "<< U<<std::endl;
    std::cout<<"\n rho: "<< parameter.potential_range<<std::endl;
    std::cout<<"\n rp: "<< rp<<std::endl;
    std::cout<<"\n u: "<< u <<std::endl;
    std::cout<<"\n cut-off range: "<< rc <<std::endl;
    V_new=V;
    std::fstream logfile;
    logfile.open("logfile.txt",std::ios::out);
    if(logfile.is_open()){

      logfile<<"This is logfile for simulation"<<std::endl;
      logfile<<"timestep: "<< dt<<std::endl;
      logfile<<"\n reduced_volume:  "<<reduced_volume<<std::endl;
      logfile<<"\n Bending Energy: "<< EnergyBending<<std::endl;
      logfile<<"\n Area Energy: "<< EnergyArea<<std::endl;
      logfile<<"\n Volume Energy: "<< EnergyVolume<<std::endl;
      logfile<<"\n Adhesion Energy: "<< EnergyAdhesion<<std::endl;
      logfile<<"\n Biased Energy: "<< E_bias<<std::endl;
      logfile<<"\n Total Energy: "<< EnergyTotal_new<<std::endl;
      logfile<<"\n number of vertices: "<< numV<<std::endl;
      logfile.close();
    }
  auto start = high_resolution_clock::now();
  for (int i=0; i< iterations; i++){
         velocity = ForceTotal/gamma;
         V_new += velocity*dt;
         Eigen::MatrixXd l;
         igl::edge_lengths(V_new,F,l);
         V_new=M1.vertex_smoothing(V_new,F);
         igl::intrinsic_delaunay_triangulation(l,F,l,F);
         E1.compute_bendingenergy_force(V_new,F,Force_Bending,EnergyBending);
         E1.compute_areaenergy_force(V_new,F,Force_Area,EnergyArea);
         E1.compute_volumeenergy_force(V_new,F,reduced_volume,Force_Volume,EnergyVolume);
         E1.compute_adhesion_energy_force(V_new,F,X,Y,Z,rp,rho,u,U,rc,Ew_t,K_bias,Force_Adhesion,EnergyAdhesion,E_bias);

         EnergyTotal_new=EnergyBending+EnergyArea+EnergyVolume+EnergyAdhesion+E_bias;
         ForceTotal=Force_Bending+Force_Area+Force_Volume+Force_Adhesion;

         EnergyChange=abs(EnergyTotal_new-EnergyTotal);
         EnergyTotal=EnergyTotal_new;

         time += dt;
         //std::cout<<"time \n "<< time<<"\n iteration \n"<<i<<std::endl;
         //std::cout<<"EnergyChange \n "<<EnergyChange<<std::endl;

        if(i%logfrequency ==0){
          std::cout<<"time "<< time<<"\n iteration "<<i<<std::endl;
          std::cout<<"\n EnergyChange:  "<<EnergyChange<<std::endl;
          std::cout<<"\n Bending Energy: "<< EnergyBending<<std::endl;
          std::cout<<"\n Area Energy: "<< EnergyArea<<std::endl;
          std::cout<<"\n Volume Energy: "<< EnergyVolume<<std::endl;
          std::cout<<"\n Adhesion Energy: "<< EnergyAdhesion<<std::endl;
          std::cout<<"\n Biased Energy: "<< E_bias<<std::endl;
          std::cout<<"\n Total Energy: "<< EnergyTotal_new<<std::endl;
          std::cout<<"\n number of vertices: "<< numV<<std::endl;
          //std::fstream logfile;
          logfile.open("logfile.txt",std::ios::app);
          if(logfile.is_open()){
            logfile<<"######"<<std::endl;
            logfile<<"time: "<< time<<std::endl;
            logfile<<"iteration: "<<i<<std::endl;
            logfile<<"\n EnergyChange:  "<<EnergyChange<<std::endl;
            logfile<<"\n Bending Energy: "<< EnergyBending<<std::endl;
            logfile<<"\n Area Energy: "<< EnergyArea<<std::endl;
            logfile<<"\n Volume Energy: "<< EnergyVolume<<std::endl;
            logfile<<"\n Adhesion Energy: "<< EnergyAdhesion<<std::endl;
            logfile<<"\n Biased Energy: "<< E_bias<<std::endl;
            logfile<<"\n Total Energy: "<< EnergyTotal_new<<std::endl;
            logfile<<"\n number of vertices: "<< numV<<std::endl;
            logfile.close();
          }
          igl::writeOFF(resfilename,V_new,F);
        }


        if (EnergyChange<tolerance && tolerance_flag != 0){//1 means con
          std::cout<<"Change of Energy is very small \n Reached Equilibrioum Shape"<<std::endl;
          std::cout<<"\n EnergyChange:  "<<EnergyChange<<std::endl;
          std::cout<<"\n Bending Energy: "<< EnergyBending<<std::endl;
          std::cout<<"\n Area Energy: "<< EnergyArea<<std::endl;
          std::cout<<"\n Volume Energy: "<< EnergyVolume<<std::endl;
          std::cout<<"\n Total Energy: "<< EnergyTotal<<std::endl;
          igl::writeOFF(outfilename,V_new,F);
          break;
        }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << duration.count() << std::endl;
    igl::writeOFF(resfilename,V_new,F);
}
void readParameter(){
    std::string line;
    std::ifstream runfile;
    runfile.open("run_file.txt");
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
    runfile >> parameter.gamma;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.particle_radious;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.adhesion_strength;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.potential_range;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.wrapping_fraction;
    getline(runfile, line);
    getline(runfile, line);
    getline(runfile, parameter.meshFile);
    getline(runfile, line);
    getline(runfile, parameter.outFile);
    getline(runfile, line);
    getline(runfile, parameter.resFile);
    getline(runfile, line);
}
