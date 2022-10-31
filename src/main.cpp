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

int numV;                                               // number of vertices
int numF;                                               // number of faces
double area_total;                                      // total area of triangle mesh
double area_avg;                                        // average area of each triangle mesh
Eigen::MatrixXd V;                                      // matrix storing vertice coordinates
Eigen::MatrixXi F;                                      // matrix storing face information, every face is one row with three integers
Eigen::MatrixXi TT, TTi;                                // face adjacency information
Eigen::MatrixXd F_normals;                              // face normals
Eigen::VectorXd dblA;                                   // store face edge information

Eigen::MatrixXd HN,H;
Eigen::SparseMatrix<double> L,M,Minv;
Eigen::VectorXd EB;


double total_bending_energy;

//double totalEB;

int main(){
  // Load a mesh in OFF format
  igl::readOFF("oblate.off", V, F);

  //Print the vertices and faces matrices
  //std::cout << "Vertices: " << std::endl << V << std::endl;
  //std::cout << "Faces:    " << std::endl << F << std::endl;
    numF = F.rows();
    numV = V.rows();
    igl::doublearea(V,F,dblA);
    std::cout << "average double area " << dblA.mean() << std::endl;
    std::cout << "total area " << dblA.sum() << std::endl;
    area_total = dblA.sum();
    area_avg = dblA.mean();
/*    Eigen::VectorXd area_voronoi;
    Eigen::VectorXd EB;
    double total_EB;
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    igl::invert_diag(M,Minv);
    HN= -Minv*(L*V);
    H = HN.rowwise().norm(); //up to sign
    area_voronoi=M.diagonal();
    EB = 0.01*(H.transpose() * area_voronoi).diagonal();
    //EB=2*0.01*H.dot(area_voronoi);
    total_EB=EB.sum();
*/
    Energy E1;
    total_bending_energy=E1.BendingEnergy(V,F);
     //Model m;
    std::cout << "total_bending_energy" <<total_bending_energy<< std::endl;
    //std::cout << "total_bending_energy" <<area_voronoi<< std::endl;
     
  // Save the mesh in OBJ format
  //igl::writeOBJ("cube.obj",V,F);
}
/*double BendingEnergy();

double BendingEnergy(Eigen::MatrixXd V1,Eigen::MatrixXi F1){
	Eigen::VectorXd area_voronoi;
        Eigen::VectorXd EB;
        double total_EB;
	igl::cotmatrix(V1,F1,L);
	igl::massmatrix(V1,F1,igl::MASSMATRIX_TYPE_VORONOI,M);
	igl::invert_diag(M,Minv);
	HN= -Minv*(L*V1);
	H = HN.rowwise().norm(); //up to sign
	area_voronoi=M.diagonal();
	EB=2*0.01*H*area_voronoi;
	total_EB=EB.sum();
	return total_EB;
}
*/

     
/*Eigen::VectorXd area_voronoi;
Eigen::VectorXd EB;
double total_EB;
//BE=BendingEnergy->BE;
/*void BE(Eigen::MatrixXd V,Eigen::MatrixXi F){
//	V= BendingEnergy->V;
//	F=BendingEnergy->F;
//	L= BendingEnergy->L;
//	M= BendingEnergy->M;
//	Minv=BendingEnergy->Minv;
//	EB=BendingEnergy->EB;
//	total_EB=BendingEnergy->totalEB;
//	area_voronoi=BendignEnergy->area-voronoi;
//	
	igl::cotmatrix(V,F,L);
	igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
	igl::invert_diag(M,Minv);
	HN = -Minv*(L*V);
	H = HN.rowwise().norm(); //up to sign
	area_voronoi=M.diagonal();
	EB=2*0.01*H;  //*area_voronoi
	total_EB=EB.sum();
	
	 std::cout << "TEB " <<total_EB<< std::endl;
}
*/


