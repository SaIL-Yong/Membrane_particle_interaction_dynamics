#include <iostream>
#include <cmath>
#include "energy.h"

void Energy::compute_bendingenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Kb, Eigen::MatrixXd& Force_Bending, double& bending_energy)
{
  bending_energy = 0.0;
  Force_Bending.setZero();
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
  igl::invert_diag(M, Minv);
  igl::gaussian_curvature(V, F, K);
  K = (Minv * K).eval();
  HN = -Minv * (L * V) / 2.0;
  H = HN.rowwise().norm(); //up to sign
  igl::per_vertex_normals(V, F, V_normals);
  H_X_N = (HN.array() * V_normals.array());
  abc = H_X_N.rowwise().sum();
  sign_of_H = abc.array().sign();
  H_signed = H.array() * sign_of_H.array();
  H_squared = H.array().square();
  area_voronoi = M.diagonal();
  EB = 2.0 * Kb * (H_squared.transpose() * area_voronoi).diagonal();
  bending_energy = EB.sum();

  Lap_H = Minv * (L * H_signed);
  force_density = (2.0 * H_signed.array() * (H_squared - K).array()) + Lap_H.array();
  vector_term = force_density.array() * area_voronoi.array();
  Force_Bending = (2.0 *Kb)*(V_normals.array().colwise() * vector_term.array());
}


void Energy::compute_areaenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Ka, double area_target, Eigen::MatrixXd& Force_Area, double& area_energy)
{
  area_energy = 0.0;
  Force_Area.setZero();
  igl::doublearea(V, F, dblA);
  area_current = dblA.sum() / 2;
  da = area_current - area_target;
  area_energy = Ka * da * da / area_target;
  scalar_term = -2.0 * Ka * da / area_target;
  AG = m.area_grad(V, F);
  Force_Area = scalar_term * AG;
}


void Energy::compute_volumeenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Kv, double volume_target, Eigen::MatrixXd& Force_Volume, double& volume_energy)
{
  volume_energy = 0.0;
  volume_current = m.cal_volume(V, F);
  Force_Volume.setZero();
  if (fabs(Kv) > EPS) {
    VG = m.volume_grad(V, F);
    dv = volume_current - volume_target;
    volume_energy = Kv * dv * dv / volume_target;
    scalar_term = -2.0 * Kv * dv / volume_target;
    Force_Volume = scalar_term * VG;
  }
}


void Energy::compute_adhesion_energy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double X, double Y, double Z,
                                           double Rp, double rho, double U, double rc, double Ew_t, double Kw,
                                           Eigen::MatrixXd& Force_Adhesion, double& EnergyAdhesion, double& EnergyBias)
{
  Force_Adhesion.setZero();
  coefficient.resize(V.rows());
  coefficient_derivative_x.resize(V.rows());
  coefficient_derivative_y.resize(V.rows());
  coefficient_derivative_z.resize(V.rows());
  distance.resize(V.rows());
  dc.resize(V.rows());
  Mod_Bias.resize(V.rows());
  Mod_Bias.setZero();
  EnergyAdhesion = 0.0;
  EnergyBias = 0.0;

  // can this loop be replaced by the eigen operation?
  for (int i = 0; i < V.rows(); i++) {
    distance(i) = sqrt((V(i,0)-X)*(V(i,0)-X)+(V(i,1)-Y)*(V(i,1)-Y)+(V(i,2)-Z)*(V(i,2)-Z));
    dc(i) = distance(i) - Rp;

    if (abs(dc(i)) > rc) continue;

    coefficient(i) = U * (exp(-(2.0*dc(i))/rho) - 2.0*exp(-dc(i)/rho));
    coefficient_derivative_x(i) = (U/(distance(i)*rho))
                                *(-exp(-(2.0*dc(i))/rho) + exp(-dc(i)/rho)) * 2.0 * (V(i,0)-X);
    coefficient_derivative_y(i) = (U/(distance(i)*rho))
                                *(-exp(-(2.0*dc(i))/rho) + exp(-dc(i)/rho)) * 2.0 * (V(i,1)-Y);
    coefficient_derivative_z(i) = (U/(distance(i)*rho))
                                *(-exp(-(2.0*dc(i))/rho) + exp(-dc(i)/rho)) * 2.0 * (V(i,2)-Z);

    if (dc(i) > EPS && fabs(Kw) > EPS) Mod_Bias(i) = 1.0;
  }

  coefficient_of_derivative.resize(V.rows(), 3);
  coefficient_of_derivative.col(0)=coefficient_derivative_x.transpose();
  coefficient_of_derivative.col(1)=coefficient_derivative_y.transpose();
  coefficient_of_derivative.col(2)=coefficient_derivative_z.transpose();
  
  First_Term = -(AG.array().colwise()*coefficient.array());
  Second_Term = -(coefficient_of_derivative.array().colwise()*area_voronoi.array());

  Sum = First_Term + Second_Term;
  Ead = coefficient.array()*area_voronoi.array();

  EnergyAdhesion = Ead.sum();

  if (fabs(Kw) > EPS) {
    dEw = EnergyAdhesion - Ew_t;
    EnergyBias = 0.5 * Kw * dEw * dEw;
    Force_Biased = Kw * dEw * (Sum.array().colwise() * Mod_Bias.array());
    Force_Adhesion = First_Term + Second_Term + Force_Biased;
  } else Force_Adhesion = First_Term + Second_Term; 
}