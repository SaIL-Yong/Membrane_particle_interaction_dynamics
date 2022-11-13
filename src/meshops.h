#pragma once
#include <iostream>
#include <cmath>
#include <igl/readOFF.h>
#include "igl/doublearea.h"
#include "igl/edges.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_face_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/cotmatrix.h>


class Mesh {
           public:

           double cal_volume(Eigen::MatrixXd V,Eigen::MatrixXi F);
           //double AreaEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F);
           //Eigen::Eigen::MatrixXd V
           Eigen::MatrixXd area_grad(Eigen::MatrixXd V,Eigen::MatrixXi F);
           Eigen::MatrixXd volume_grad(Eigen::MatrixXd V,Eigen::MatrixXi F);

};
