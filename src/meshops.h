#pragma once
#include <iostream>
#include <cmath>
#include <igl/readOFF.h>
#include <igl/doublearea.h>
#include <igl/edges.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_face_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>
#include <igl/unique_edge_map.h>
#include <igl/doublearea.h>
#include <igl/oriented_facets.h>
#include <igl/edge_flaps.h>
#include <igl/edges.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Geometry>
#include <vector>
#include <igl/adjacency_matrix.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/gaussian_curvature.h>
#include <igl/intrinsic_delaunay_triangulation.h>
#include <igl/barycenter.h>


class Mesh {
 public:
  void mesh_cal(Eigen::MatrixXd V, Eigen::MatrixXi F);
  double cal_volume2(Eigen::MatrixXd V, Eigen::MatrixXi F);
  double cal_volume(Eigen::MatrixXd V, Eigen::MatrixXi F);
  //double AreaEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F);
  //Eigen::Eigen::MatrixXd V
  Eigen::MatrixXd area_grad(Eigen::MatrixXd V, Eigen::MatrixXi F);
  Eigen::MatrixXd volume_grad(Eigen::MatrixXd V, Eigen::MatrixXi F);
  Eigen::MatrixXd vertex_smoothing(Eigen::MatrixXd V, Eigen::MatrixXi F);

  int numF, numV;
  double volume_total, area_total;
  Eigen::SparseMatrix<double> L, M, Minv;
  Eigen::MatrixXd H, HN, H_squared, H_X_N;
  Eigen::VectorXd abc, sign_of_H;
  Eigen::VectorXd H_signed;
  Eigen::VectorXd K;
  Eigen::MatrixXd F_normals, V_normals;
  Eigen::VectorXd dblA;
  std::vector<std::vector<double>> VF;
  std::vector<std::vector<double>> VFi;
  Eigen::VectorXd area_voronoi;

  // //Moment of Inertia
  // Eigen::Vector3d COM,r;
  // Eigen::Matrix3d I;
  // //Eigen::Vector3d r = Eigen::Vector3d::Zero();


 private:
  double sum;

  double p0x, p0y, p0z, p1x, p1y, p1z, p2x, p2y, p2z;
  double v321, v231, v312, v132, v213, v123;

  Eigen::MatrixXd AG, VG;
  Eigen::RowVector3d vol_ij;

  Eigen::MatrixXd V_b, V_new;
  double area_avg, lambda, face_area_sum;
  Eigen::MatrixXd f_norm, V_centroid;
  Eigen::VectorXd face_area;
  Eigen::VectorXd sum_of_area_centroid, V_avg, f_norm_sum;
};
