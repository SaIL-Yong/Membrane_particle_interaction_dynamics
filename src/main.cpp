//#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <cmath>
#include <igl/readOFF.h>
#include <igl/per_face_normals.h>
#include <igl/writeOBJ.h>
#include <igl/adjacency_list.h>
#include "igl/unique_edge_map.h"
#include "igl/doublearea.h"
#include <iostream>
#include "igl/oriented_facets.h"
#include "igl/edge_flaps.h"
#include "igl/edges.h"
#include "igl/triangle_triangle_adjacency.h"
#include "Eigen/Geometry"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

#include "energy.h" 
#include "meshops.h"

int numV;                                               // number of vertices
int numF;                                               // number of faces
double area_avg;                                        // average area of each triangle mesh
Eigen::MatrixXd V;                                      // matrix storing vertice coordinates
Eigen::MatrixXi F;                                      // matrix storing face information, every face is one row with three integers
Eigen::MatrixXi TT, TTi;                                // face adjacency information
Eigen::MatrixXd F_normals;                              // face normals
Eigen::VectorXd dblA;                                   // store face edge information

Eigen::MatrixXd HN,H;
Eigen::SparseMatrix<double> L,M,Minv;
Eigen::VectorXd EB;


//double total_bending_energy;

//double totalEB;

int main(){
  // Load a mesh in OFF format
  igl::readOFF("oblate.off", V, F);

  //Print the vertices and faces matrices
  //std::cout << "Vertices: " << std::endl << V << std::endl;
  //std::cout << "Faces:    " << std::endl << F << std::endl;
//   numF = F.rows();
    float reduced_volume=0.7;
    numV = V.rows();
    igl::doublearea(V,F,dblA);
    double area_total = dblA.sum();
    area_avg = dblA.mean();
    std::cout << "average double area " << dblA.mean() << std::endl;
    std::cout << "total area " << dblA.sum()/2 << std::endl;
    Energy E1;
    double total_bending_energy=E1.BendingEnergy(V,F);
     //Model m;
    std::cout << "total_bending_energy" <<total_bending_energy<< std::endl;
    double area_energy=E1.AreaEnergy(V,F);
    std::cout << "area_energy" <<area_energy<< std::endl;
    double volume_energy=E1.VolumeEnergy(V,F,reduced_volume);
    std::cout << "volume_energy" <<volume_energy<< std::endl;
    
    Mesh M1;
    double volume= M1.cal_volume(V,F);
    std::cout << "total_volume" <<volume<< std::endl;
     
  // Save the mesh in OBJ format
  //igl::writeOBJ("cube.obj",V,F);
}

     



