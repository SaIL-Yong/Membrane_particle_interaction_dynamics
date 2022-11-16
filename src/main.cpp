#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <cmath>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/writeOBJ.h>
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
#include "meshops.h"
#include "energy.h"
int numV;                                               // number of vertices
int numF;                                               // number of faces
double area_avg;                                        // average area of each triangle mesh
Eigen::MatrixXd V,V_new;                                      // matrix storing vertice coordinates
Eigen::MatrixXi F;

int main(){
    igl::readOFF("oblate.off", V, F);
    numF = F.rows();
    numV = V.rows();
    int iterations=10000;

    Eigen::MatrixXd ForceArea,ForceVolume,ForceBending,velocity,ForceTotal; //forces

    double dt=0.001,time; //time step
    double gamma=2.0;
    double tolerance=1.0e-7;
    float reduced_volume=0.65;
    double EnergyVolume,EnergyArea,EnergyBending,EnergyTotal,EnergyTotal_new, EnergyChange;  //energies


    Energy E1;
    E1.compute_bendingenergy_force(V,F,ForceBending,EnergyBending);
    E1.compute_areaenergy_force(V,F,ForceArea,EnergyArea);
    E1.compute_volumeenergy_force(V,F,reduced_volume,ForceVolume,EnergyVolume);
    EnergyTotal=EnergyBending+EnergyArea+EnergyVolume;
    ForceTotal=ForceBending+ForceArea+ForceVolume;
    std::cout<<"Bending EnergyE \n"<< EnergyBending<<std::endl;
//    velocity= ForceTotal/gamma;//testing
    V_new=V;
//    V_new += velocity*dt;
  for (int i=0; i< iterations; i++){
         velocity=ForceTotal/gamma;
         V_new += velocity*dt;
         Eigen::MatrixXd l;
         //igl::edge_lengths(V,F,l);
         //igl::intrinsic_delaunay_triangulation(l,F,l,F);
         E1.compute_bendingenergy_force(V_new,F,ForceBending,EnergyBending);
         E1.compute_areaenergy_force(V_new,F,ForceArea,EnergyArea);
         E1.compute_volumeenergy_force(V_new,F,reduced_volume,ForceVolume,EnergyVolume);
         ForceTotal=ForceBending+ForceArea+ForceVolume;
         EnergyTotal_new=EnergyBending+EnergyArea+EnergyVolume;
         EnergyChange=abs(EnergyTotal_new-EnergyTotal);
         EnergyTotal=EnergyTotal_new;

         time +=dt;
         std::cout<<"time \n "<< time<<"\n iteration \n"<<i<<std::endl;
         std::cout<<"EnergyChange \n "<<EnergyChange<<std::endl;


        if (EnergyChange<tolerance){
          std::cout<<"Bending EnergyE \n"<< EnergyBending<<std::endl;
          std::cout<<"Area Energy \n"<< EnergyArea<<std::endl;
          std::cout<<"Volume Energy \n"<< EnergyVolume<<std::endl;
          std::cout<<"Total Energy \n"<< EnergyTotal<<std::endl;
          igl::writeOFF("cube.off",V_new,F);
          break;
        }
    }
}
