#pragma once
#include <iostream>
#include <cmath>
#include <igl/readOFF.h>
#include <igl/per_face_normals.h>
#include "igl/doublearea.h"
#include "igl/edges.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


//void BendingEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F);




class Energy {
           public:
           /*Eigen::MatrixXd V;
           Eigen::MatrixXi F;
           Eigen::MatrixXd HN;
           Eigen::SparseMatrix<double> L,M,Minv;
           Eigen::VectorXd area_voronoi;
           Eigen::VectorXd EB;
           double total_EB;*/
           double BendingEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F);
 
};   


