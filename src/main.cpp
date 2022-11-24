#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <chrono>
#include <cmath>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>
#include <igl/unique_edge_map.h>
#include <igl/doublearea.h>
#include <igl/oriented_facets.h>
#include <igl/edge_flaps.h>
#include <igl/edges.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Geometry>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <vector>
#include <igl/adjacency_matrix.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/gaussian_curvature.h>
#include <igl/intrinsic_delaunay_triangulation.h>
#include <igl/barycenter.h>
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
    igl::readOFF(filename, V, F);
    numF = F.rows();
    numV = V.rows();
    int iterations=parameter.iteration;
    int logfrequency=100;
    Eigen::MatrixXd ForceArea,ForceVolume,ForceBending,velocity,ForceTotal; //forces

    double dt=parameter.dt; //time step
    double gamma=parameter.gamma;
    double tolerance=parameter.tolerance;
    float reduced_volume=parameter.reduced_volume;
    double EnergyVolume,EnergyArea,EnergyBending,EnergyTotal,EnergyTotal_new, EnergyChange;  //energies

    Mesh M1;
    Energy E1;
    E1.compute_bendingenergy_force(V,F,ForceBending,EnergyBending);
    E1.compute_areaenergy_force(V,F,ForceArea,EnergyArea);
    E1.compute_volumeenergy_force(V,F,reduced_volume,ForceVolume,EnergyVolume);
    EnergyTotal=EnergyBending+EnergyArea+EnergyVolume;
    ForceTotal=ForceBending+ForceArea+ForceVolume;
    std::cout<<"Bending Energy \n"<< EnergyBending<<std::endl;
    std::cout<<"\n reduced_volume "<< reduced_volume<<std::endl;
    V_new=V;
  auto start = high_resolution_clock::now();
  for (int i=0; i< iterations; i++){
         velocity=ForceTotal/gamma;
         V_new += velocity*dt;
         Eigen::MatrixXd l;
         igl::edge_lengths(V_new,F,l);
         V_new=M1.vertex_smoothing(V_new,F);
         igl::intrinsic_delaunay_triangulation(l,F,l,F);
         E1.compute_bendingenergy_force(V_new,F,ForceBending,EnergyBending);
         E1.compute_areaenergy_force(V_new,F,ForceArea,EnergyArea);
         E1.compute_volumeenergy_force(V_new,F,reduced_volume,ForceVolume,EnergyVolume);
         ForceTotal=ForceBending+ForceArea+ForceVolume;
         EnergyTotal_new=EnergyBending+EnergyArea+EnergyVolume;
         EnergyChange=abs(EnergyTotal_new-EnergyTotal);
         EnergyTotal=EnergyTotal_new;

         time +=dt;
         //std::cout<<"time \n "<< time<<"\n iteration \n"<<i<<std::endl;
         //std::cout<<"EnergyChange \n "<<EnergyChange<<std::endl;

        if(i%logfrequency ==0){
          std::cout<<"time "<< time<<"\n iteration "<<i<<std::endl;
          std::cout<<"\n EnergyChange:  "<<EnergyChange<<std::endl;
          std::cout<<"\n Bending Energy: "<< EnergyBending<<std::endl;
          std::cout<<"\n Area Energy: "<< EnergyArea<<std::endl;
          std::cout<<"\n Volume Energy: "<< EnergyVolume<<std::endl;
          std::cout<<"\n Total Energy: "<< EnergyTotal<<std::endl;
          std::cout<<"\n number of vertices: "<< numV<<std::endl;

        }


        if (EnergyChange<tolerance){
          std::cout<<"Change of Energy is very small \n Reached Equilibrioum Shape"<<std::endl;
          std::cout<<"\n EnergyChange:  "<<EnergyChange<<std::endl;
          std::cout<<"\n Bending Energy: "<< EnergyBending<<std::endl;
          std::cout<<"\n Area Energy: "<< EnergyArea<<std::endl;
          std::cout<<"\n Volume Energy: "<< EnergyVolume<<std::endl;
          std::cout<<"\n Total Energy: "<< EnergyTotal<<std::endl;
          igl::writeOFF("final_icosahedron5.off",V_new,F);
          break;
        }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    //igl::writeOFF("restart.off",V_new,F);
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
    getline(runfile, parameter.meshFile);
    getline(runfile, line);
    getline(runfile, parameter.outFile);
    getline(runfile, line);
    getline(runfile, parameter.resFile);
    getline(runfile, line);
}
