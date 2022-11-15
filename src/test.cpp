#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <cmath>
#include <igl/readOFF.h>
#include <igl/per_vertex_normals.h>
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
#include <vector>
#include <igl/adjacency_matrix.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/gaussian_curvature.h>
#include "meshops.h"
#include "energy.h"
int numV;                                               // number of vertices
int numF;                                               // number of faces
double area_avg;                                        // average area of each triangle mesh
Eigen::MatrixXd V;                                      // matrix storing vertice coordinates
Eigen::MatrixXi F;

int main(){
    igl::readOFF("oblate.off", V, F);
    numF = F.rows();
    numV = V.rows();

    Eigen::MatrixXd FV;
    float reduced_volume=0.7;
    double EV1;
    //Energy EN;
    //EN.compute_bendingenergy_force(V,F,FB,EB1);
    Energy E1;
    E1.compute_volumeenergy_force(V,F,reduced_volume,FV,EV1);

  // Mesh m;
   //double volume=m.cal_volume(V,F);
  //  Eigen::MatrixXd VG=m.volume_grad(V,F);

    std::cout<<"total_energy_area"<< EV1<<std::endl;
    std::cout<<"force_area"<<FV<<std::endl;
}

/*
Eigen::VectorXd Lap_H;
Eigen::MatrixXd V_normals,Force_Bending;
Eigen::VectorXd dblA;
igl::per_vertex_normals(V,F,V_normals);
igl::doublearea(V,F,dblA);
Eigen::MatrixXd H,HN;
Eigen::SparseMatrix<double> L,M,Minv;
Eigen::VectorXd area_voronoi,K,H_squared;
Eigen::VectorXd EB,force_density;
igl::cotmatrix(V,F,L);
igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
igl::invert_diag(M,Minv);
igl::gaussian_curvature(V,F,K);
K = (Minv*K).eval();
HN= -Minv*(L*V)/2.0;
H = HN.rowwise().norm(); //up to sign
Eigen::MatrixXd H_X_N= (HN.array()*V_normals.array());
Eigen::VectorXd abc=H_X_N.rowwise().sum();
Eigen::VectorXd sign_of_H= abc.array().sign();
Eigen::VectorXd H_signed=H.array()*sign_of_H.array();
H_squared=H.array().square();
area_voronoi=M.diagonal();
Lap_H=Minv*(L*H_signed);
force_density=(2.0*H_signed.array()*(H_squared-K).array() )+ Lap_H.array();
Eigen::VectorXd vector_term=force_density.array()*area_voronoi.array();


Force_Bending=(2*0.01)*(V_normals.array().colwise()*vector_term.array());//*area_voronoi.array();

Eigen::MatrixXd V_normals,Force_Bending;
igl::per_vertex_normals(V,F,V_normals);
Eigen::MatrixXd H,HN,H_squared;
Eigen::SparseMatrix<double> L,M,Minv;
Eigen::VectorXd area_voronoi,Lap_H,force_density,K;
Eigen::VectorXd EB;
igl::cotmatrix(V,F,L);
igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
igl::invert_diag(M,Minv);
igl::gaussian_curvature(V,F,K);
K = (Minv*K).eval();
HN= -Minv*(L*V)/2.0;
H = HN.rowwise().norm(); //up to sign
Eigen::MatrixXd H_X_N= (HN.array()*V_normals.array());
Eigen::VectorXd abc=H_X_N.rowwise().sum();
Eigen::VectorXd sign_of_H= abc.array().sign();
Eigen::VectorXd H_signed=H.array()*sign_of_H.array();
H_squared=H.array().square();
area_voronoi=M.diagonal();
EB = 2.0*0.01*(H_squared.transpose() * area_voronoi).diagonal();
double total_EB=EB.sum();
Lap_H=Minv*(L*H_signed);
force_density=(2.0*H_signed.array()*(H_squared-K).array() )+ Lap_H.array();
Eigen::VectorXd vector_term=force_density.array()*area_voronoi.array();
Force_Bending=(2*0.01)*(V_normals.array().colwise()*vector_term.array());
std::cout<<"force_density"<<Force_Bending<<std::endl;
*/
