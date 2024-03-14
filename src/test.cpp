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
////

    float bending_modulus=parameter.Kb;
    float rp=parameter.particle_radious;
    float u=parameter.adhesion_strength;
    float rho=(parameter.potential_range)*rp;
    float U=(bending_modulus*u)/(pow(rp,2)) ;
    float rc=5*rho;
    float X=V.maxCoeff()+rp+rho*1,Y=0.0,Z=0.0;
    float chi=parameter.wrapping_fraction;
    float Area_w_t=chi*4.0*PI*(pow(rp,2));
    float Ew_t=-U*Area_w_t;
    float K_bias=10.0;

/////
    float reduced_volume=parameter.reduced_volume;
    Energy E1;
    Eigen::MatrixXd Force_Area,Force_Volume,Force_Bending,Force_Adhesion,velocity,ForceTotal; //forces

    double EnergyVolume,EnergyArea,EnergyBending,EnergyTotal,EnergyAdhesion,EnergyTotal_new, EnergyChange;  //energies
    //
    E1.compute_bendingenergy_force(V,F,Force_Bending,EnergyBending);
    E1.compute_areaenergy_force(V,F,Force_Area,EnergyArea);
    E1.compute_volumeenergy_force(V,F,reduced_volume,Force_Volume,EnergyVolume);
    E1.compute_adhesion_energy_force(V,F,X,Y,Z,rp,rho,u,U,rc,Ew_t,K_bias,Force_Adhesion,EnergyAdhesion);

    EnergyTotal=EnergyBending+EnergyArea+EnergyVolume+EnergyAdhesion;
    //
    // float chi=parameter.wrapping_fraction;
    // float Area_w_t=chi*4.0*PI*(pow(rp,2));
    // double Ew_t=-U*Area_w_t;
    // std::cout<<"Biased-Energy::"<<Area_w_t<<std::endl;
    // std::cout<<"Biased-Energy::"<<Ew_t<<std::endl;
    // std::vector<std::vector<double> > VF;
    // std::vector<std::vector<double> > VFi;
    // igl::vertex_triangle_adjacency(V,F,VF, VFi);
    // for (int i = 0; i < 2; i++){
    // for (int j = 0; j < VF[i].size(); j++){
    //     std::cout << "Adjacent Faces::"<<VF[i][j]<< " " ;}
    // std::cout <<std::endl;}
    // for (int i = 0; i < 2; i++){
    // for (int j = 0; j < VFi[i].size(); j++){
    //     std::cout << VFi[i][j]<< " " ;}
    // std::cout <<std::endl;}
    // std::vector<std::vector<double> > A;
    // igl::adjacency_list(F,A);
    // for (int i = 0; i < 2; i++){
    // for (int j = 0; j < A[i].size(); j++){
    //     std::cout << "Adjacent Vertices:"<<A[i][j]<< " " ;}
    // std::cout <<std::endl;}
    //Eigen::MatrixXd Mod_Bias= Eigen::MatrixXd::Zero(V.rows(),3);
    //std::cout<<"Matrix::"<<Mod_Bias<<std::endl;

    //std::cout<<"Vertices::"<<VFi[0][1]<<std::endl;
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
