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


class Mesh {
           public:

           double cal_volume(Eigen::MatrixXd V,Eigen::MatrixXi F);
           //double AreaEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F);
 
};   



