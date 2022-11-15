#pragma once
#include <iostream>
#include <cmath>
#include <igl/readOFF.h>
#include <igl/per_vertex_normals.h>
#include "igl/doublearea.h"
#include "igl/edges.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/gaussian_curvature.h>
#ifndef PI
#define PI 3.14159265358979323846
#endif

class Energy {
public:
           //double BendingEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F);
           //double AreaEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F);
           //double VolumeEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F,double reduced_volume);
           void compute_bendingenergy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,Eigen::MatrixXd& Force_Bending,double& total_EB);
           void compute_areaenergy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,Eigen::MatrixXd& Force_Area,double& area_energy);
           void compute_volumeenergy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,float reduced_volume,Eigen::MatrixXd& Force_Volume,double& volume_energy);
           /*struct calculated_values{
             Eigen::MatrixXd HN,H,H_squared;
             Eigen::SparseMatrix<double> L,M,Minv;
             Eigen::VectorXd area_voronoi;
             Eigen::VectorXd EB;
             double total_EB;

           }; */

};
