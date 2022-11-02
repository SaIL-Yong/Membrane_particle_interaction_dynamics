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
#ifndef PI
#define PI 3.14159265358979323846
#endif



class Energy {
           public:

           double BendingEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F);
           double AreaEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F);
           double VolumeEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F,double reduced_volume);
 
};   


