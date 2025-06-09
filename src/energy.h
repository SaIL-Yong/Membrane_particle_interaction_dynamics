#pragma once

#include <iostream>
#include <cmath>
#include <random>
#include <igl/readOFF.h>
#include <igl/centroid.h>
#include <igl/per_vertex_normals.h>
#include <igl/doublearea.h>
#include <igl/edges.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/gaussian_curvature.h>
#include <Eigen/Dense>
#include "meshops.h"
#include "meshedvesicle.h"

// Remove the PI macro to avoid conflicts with igl::PI
// If you need π in implementation files, use igl::PI or ::PI (from simulation.h)

#ifndef EPS
static constexpr double EPS = 1.0e-10;
#endif

class Energy {
 public:
  void compute_bendingenergy_force(
      const Eigen::MatrixXd& V,
      const Eigen::MatrixXi& F,
      double Kb,
      Eigen::MatrixXd&       Force_Bending,
      double&                bending_energy,
      Mesh                   m
  );

  void compute_areaenergy_force(
      const Eigen::MatrixXd& V,
      const Eigen::MatrixXi& F,
      double                 Ka,
      double                 area_target,
      Eigen::MatrixXd&       Force_Area,
      double&                area_energy,
      Mesh                   m
  );

  void compute_volumeenergy_force(
      const Eigen::MatrixXd& V,
      const Eigen::MatrixXi& F,
      double                 Kv,
      double                 volume_target,
      Eigen::MatrixXd&       Force_Volume,
      double&                volume_energy,
      Mesh                   m
  );

void compute_adhesion_energy_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& closest_points,
    double                 rho,
    double                 U,
    double                 r_equilibrium,
    double                 rc,
    int                    angle_flag,
    int                    vesicle_position,
    double                 Ew_t,
    double                 Kw,
    Eigen::MatrixXd&       Force_Adhesion,
    Eigen::MatrixXd&       Force_Repulsion,
    const Eigen::MatrixXd& signed_distance,
    double&                EnergyAdhesion,
    double&                EnergyBias,
    Mesh                   m
);

  void redistributeAdhesionForce(
      const Eigen::MatrixXd& V2,
      const Eigen::MatrixXi& F2,
      const Eigen::MatrixXd& closest_points,
      const Eigen::MatrixXd& Force_Adhesion,
      const Eigen::VectorXi& facet_index,
      Eigen::MatrixXd&       ForcesOnVertices
  );

  void compute_drag_force(
      const Eigen::MatrixXd& velocity,
      double                 gamma,
      double                 mass,
      Eigen::MatrixXd&       Force_Drag
  );

  void compute_drag_force_vesicle(
      const Eigen::Vector3d& velocity,
      double                 gamma,
      double                 mass,
      Eigen::Vector3d&       Force_Drag_vesicle
  );

  std::mt19937 gen;  // Mersenne Twister random number generator

  Energy() : gen(std::random_device{}()) {}

  void compute_random_force(
      double         gamma,
      double         kbT,
      double         mass,
      double         dt,
      Eigen::MatrixXd& Force_Random
  );

 private:
  // Variables for bending energy force calculation
  int                     f_idx;
  Eigen::Vector3i         face;
  Eigen::Matrix<double,3,3> triangle;
  Eigen::VectorXd         Lap_H, force_density;
  Eigen::VectorXd         EB;
  Eigen::VectorXd         vector_term;
  Eigen::Vector3d         force;

  // Variables for area and volume energy force calculation
  Eigen::MatrixXd         AG;
  Eigen::MatrixXd         VG;
  double                  da, dv, scalar_term;
  double                  dx, dy, dz, r2, dc;

  // Variables for adhesion energy force calculation
  Eigen::VectorXd         coefficient_derivative_x;
  Eigen::VectorXd         coefficient_derivative_y;
  Eigen::VectorXd         coefficient_derivative_z;
  Eigen::VectorXd         coefficient;
  Eigen::VectorXd         distance;
  Eigen::VectorXd         lennard_jones_potential;
  Eigen::MatrixXd         coefficient_of_derivative;
  Eigen::MatrixXd         First_Term;
  Eigen::MatrixXd         Second_Term;
  Eigen::MatrixXd         Sum;
  Eigen::MatrixXd         Ead;
  Eigen::MatrixXd         lennard_jones_force;
  Eigen::MatrixXd         Force_Biased;
  Eigen::MatrixXd         comvec;
  Eigen::RowVector3d      vesicle_center;
  Eigen::Vector3d         r_ij;
  Eigen::RowVector3d      r_ij_transpose;
  Eigen::Vector3d         vesicle_centroid;
  double                  angle, dEw;
  double                  s_by_r6, dc_mag;
};
// Compare this snippet from energy.cpp