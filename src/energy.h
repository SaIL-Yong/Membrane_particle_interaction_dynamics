#pragma once
#include <iostream>
#include <cmath>
#include <random>
#include <igl/readOFF.h>
#include<igl/centroid.h>
#include <igl/per_vertex_normals.h>
#include "igl/doublearea.h"
#include "igl/edges.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/gaussian_curvature.h>
#include <Eigen/Dense>
#include "meshops.h"
#include "meshedparticle.h"
#ifndef PI
#define PI 3.14159265358979323846
#endif
#ifndef EPS
const double EPS=1.0e-10;
#endif

class Energy {
 public:
 
  void compute_bendingenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Kb, Eigen::MatrixXd& Force_Bending, double& bending_energy, Mesh m);
  void compute_areaenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Ka, double area_target, Eigen::MatrixXd& Force_Area, double& area_energy, Mesh m);
  void compute_volumeenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Kv, double volume_target, Eigen::MatrixXd& Force_Volume, double& volume_energy, Mesh m);
  
  void compute_adhesion_energy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd V_particle,double rho, double U,double r_equilibrium,double rc,
                                          int angle_flag,int particle_position,double Ew_t, double Kw,Eigen::MatrixXd& Force_Adhesion,
                                          Eigen::MatrixXd& Force_Repulsion,Eigen::MatrixXd signed_distance, double& EnergyAdhesion,double& EnergyBias,  Mesh m);

 void redistributeAdhesionForce(Eigen::MatrixXd V2, Eigen::MatrixXi F2, Eigen::MatrixXd closest_points, Eigen::MatrixXd Force_Adhesion, Eigen::VectorXi facet_index, Eigen::MatrixXd& ForcesOnVertices) ; 


 void compute_drag_force(Eigen::MatrixXd velocity, double gamma, double mass, Eigen::MatrixXd& Force_Drag);
 std::mt19937 gen;  // Mersenne Twister random number generator

 Energy() : gen(std::random_device{}()) {}

 void compute_random_force(double gamma, double kbT, double mass, double dt, Eigen::MatrixXd& Force_Random);

 private:
  // variables for bending energy force calculation
  int f_idx;
  Eigen::Vector3i face;
  Eigen::Matrix<double, 3, 3> triangle;
  Eigen::VectorXd Lap_H, force_density;
  Eigen::VectorXd EB;
  Eigen::VectorXd vector_term;
  Eigen::Vector3d force;

  // variables for area and volume energy force calculation
  Eigen::MatrixXd AG;
  Eigen::MatrixXd VG;
  double da, dv, scalar_term;
  double dx,dy,dz,r2,dc;

  // variables for adhesion energy force calculation
  Eigen::VectorXd coefficient_derivative_x, coefficient_derivative_y, coefficient_derivative_z, coefficient, distance;//, dc;
  Eigen::VectorXd lennard_jones_potential;
  //Eigen::VectorXd Mod_Bias;
  Eigen::MatrixXd coefficient_of_derivative;
  Eigen::MatrixXd First_Term, Second_Term, Sum, Ead,lennard_jones_force;
  Eigen::MatrixXd Force_Biased;
  Eigen::MatrixXd comvec;
  Eigen::RowVector3d particle_center;
  Eigen::Vector3d r_ij;
  Eigen::RowVector3d r_ij_transpose;
  Eigen::Vector3d particle_centroid;
  double angle, dEw;
  double s_by_r6,dc_mag;
};
