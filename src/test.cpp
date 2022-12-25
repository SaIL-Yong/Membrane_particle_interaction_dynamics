#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <igl/edges.h>
#include <igl/grad.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Geometry>
#include <vector>
#include <igl/intrinsic_delaunay_triangulation.h>
#include <igl/writeOFF.h>
#include "meshops.h"
#include "energy.h"
#include "parameters.h"
using namespace std::chrono;
int numV;                                               // number of vertices
int numF;                                               // number of faces
double area_avg;                                        // average area of each triangle mesh
Eigen::MatrixXd V,V_new,F_normals;                                      // matrix storing vertice coordinates
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
    double rp=parameter.particle_radious;

    Energy EN;
    Eigen::MatrixXd ForceArea,ForceVolume,ForceBending,velocity,ForceTotal;
    double EnergyVolume,EnergyArea,EnergyBending,EnergyTotal,EnergyTotal_new, EnergyChange;
    Eigen::MatrixXd Force_Adhesion;
    double adhesion_energy;
    EN.compute_bendingenergy_force(V,F,ForceBending,EnergyBending);
    EN.compute_areaenergy_force(V,F,ForceArea,EnergyArea);
    EN.compute_adhesion_energy_force(V,F,Force_Adhesion,adhesion_energy);
    //std::cout<< V(0,2)<<std::endl;
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
    getline(runfile, parameter.meshFile);
    getline(runfile, line);
    getline(runfile, parameter.outFile);
    getline(runfile, line);
    getline(runfile, parameter.resFile);
    getline(runfile, line);
}

// #include "meshops.h"
// #include "energy.h"
// #include "parameters.h"
// using namespace std::chrono;
// int numV;                                               // number of vertices
// int numF;                                               // number of faces
// double area_avg;                                        // average area of each triangle mesh
// Eigen::MatrixXd V,V_new,area_voronoi;                                      // matrix storing vertice coordinates
// Eigen::MatrixXi F;
// Eigen::VectorXd dblA;
// Parameter parameter;
//
// int main(){
//     readParameter();
//     std::string filename = parameter.meshFile;
//     igl::readOFF(filename, V, F);
//     numF = F.rows();
//     numV = V.rows();
//     //igl::doublearea(V,F,dblA);
//     int iterations=parameter.iterations;
//     int logfrequency=100;
//     Eigen::MatrixXd ForceArea,ForceVolume,ForceBending,velocity,ForceTotal; //forces
//
//     double dt=parameter.dt; //time step
//     double gamma=parameter.gamma;
//     double tolerance=parameter.tolerance;
//     float reduced_volume=parameter.reduced_volume;
//     double EnergyVolume,EnergyArea,EnergyBending,EnergyTotal,EnergyTotal_new, EnergyChange;  //energies
//
//     Mesh M1;
//     Energy E1;
//     E1.compute_bendingenergy_force(V,F,ForceBending,EnergyBending);
//     E1.compute_areaenergy_force(V,F,ForceArea,EnergyArea);
//     E1.compute_volumeenergy_force(V,F,reduced_volume,ForceVolume,EnergyVolume);
//     EnergyTotal=EnergyBending+EnergyArea+EnergyVolume;
//     ForceTotal=ForceBending+ForceArea+ForceVolume;
//     std::cout<<"Bending Energy \n"<< EnergyBending<<std::endl;
//     std::cout<<"\n reduced_volume "<< reduced_volume<<std::endl;
//     //std::cout<<"\n force_area "<<ForceArea<<std::endl;
//
//
// }
//
// void readParameter(){
//     std::string line;
//     std::ifstream runfile;
//     runfile.open("run_file.txt");
//     getline(runfile, line);
//     runfile >> parameter.iterations;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.dt;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.Kb;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.Ka;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.Kv;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.reduced_volume;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.tolerance;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.gamma;
//     getline(runfile, line);
//     getline(runfile, line);
//     getline(runfile, parameter.meshFile);
//     getline(runfile, line);
//     getline(runfile, parameter.outFile);
//     getline(runfile, line);
//     getline(runfile, parameter.resFile);
//     getline(runfile, line);
// }
