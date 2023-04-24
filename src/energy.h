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
#include <Eigen/Dense>
#include "meshops.h"
#ifndef PI
#define PI 3.14159265358979323846
#endif
#ifndef EPS
#define EPS 1.0e-10
#endif

class Energy {
 public:
  // float rp=parameter.particle_radious;
  // float u=parameter.adhesion_strength;
  // float rho=(parameter.potential_range)*rp;
  // float U=(bending_modulus*u)/(pow(rp,2)) ;
      //float X,Y,Z,rp,u,rho,U;
   void compute_bendingenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Kb,double C_0, Eigen::MatrixXd& Force_Bending, double& bending_energy, Mesh m);
  void compute_areaenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Ka, double area_target, Eigen::MatrixXd& Force_Area, double& area_energy, Mesh m);
  void compute_volumeenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Kv, double volume_target, Eigen::MatrixXd& Force_Volume, double& volume_energy, Mesh m);
  //Eigen::VectorXd area_voronoi;
  //void compute_localareaenergy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,Eigen::MatrixXd& Force_Local_Area,double& local_area_energy);
  //void dosomething();
  void compute_adhesion_energy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::RowVector3d COM,
                                     double rp, double rho, double U, double rc, int angle_flag, int particle_position, double Ew_t, double Kw,
                                     Eigen::MatrixXd& Force_Adhesion, double& EnergyAdhesion, double& E_bias, Eigen::RowVector3d& PF, Mesh m);

 private:
  // variables for bending energy force calculation
  Eigen::VectorXd Lap_H, force_density;
  Eigen::VectorXd EB;
  Eigen::VectorXd vector_term;

  // variables for area and volume energy force calculation
  Eigen::MatrixXd AG;
  Eigen::MatrixXd VG;
  double da, dv, scalar_term;

  // variables for adhesion energy force calculation
  Eigen::VectorXd coefficient_derivative_x, coefficient_derivative_y, coefficient_derivative_z, coefficient, distance, dc;
  //Eigen::VectorXd Mod_Bias;
  Eigen::MatrixXd coefficient_of_derivative;
  Eigen::MatrixXd First_Term, Second_Term, Sum, Ead;
  Eigen::MatrixXd Force_Biased;
  Eigen::MatrixXd comvec;
  double angle, dEw;
};
