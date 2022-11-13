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
#include <igl/vertex_triangle_adjacency.h>
#include <igl/cotmatrix.h>
#include "energy.h"
#include "meshops.h"

class Force {
           public:
           Eigen::MatrixXd force_bending(Eigen::MatrixXd V,Eigen::MatrixXi F);
           Eigen::MatrixXd force_area(Eigen::MatrixXd V,Eigen::MatrixXi F);
           Eigen::MatrixXd force_volume(Eigen::MatrixXd V,Eigen::MatrixXi F);

};
