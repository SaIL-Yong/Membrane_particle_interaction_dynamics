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
#include <vector>
#include <igl/adjacency_matrix.h>
#include <igl/vertex_triangle_adjacency.h>
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
    Eigen::MatrixXd F_normals;
    Eigen::VectorXd dblA;
    int numV=V.rows();
    igl::per_face_normals(V,F,F_normals);
    igl::doublearea(V,F,dblA);
    Energy EN;
    double bendingenergy=EN.BendingEnergy(V,F);
    std:: cout<< "Mean_Curvature"<<bendingenergy <<std::endl;
}
