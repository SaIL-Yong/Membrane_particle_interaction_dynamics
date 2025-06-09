//meshops.cpp
#include "energy.h"
#include "meshops.h"
#include "parameters.h"
#include <omp.h>


void Mesh::mesh_cal(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
  numV = V.rows();
  numF = F.rows();
  igl::cotmatrix(V, F, L);
  igl::vertex_triangle_adjacency(numV, F, VF, VFi);
  igl::per_face_normals(V, F, F_normals);
  igl::per_vertex_normals(V, F, V_normals);
  igl::doublearea(V, F, dblA);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
  igl::invert_diag(M, Minv);
  igl::gaussian_curvature(V, F, K);
  

  K = (Minv * K).eval();
  HN = -Minv * (L * V) / 2.0;
  H = HN.rowwise().norm(); //up to sign

  H_X_N = (HN.array() * V_normals.array());
  abc = H_X_N.rowwise().sum();
  sign_of_H = abc.array().sign();
  H_signed = H.array() * sign_of_H.array();
  H_squared = H.array().square();
  area_voronoi.resize(numV);
  area_voronoi = M.diagonal();

  volume_total = cal_volume2(V, F);
  area_total = dblA.sum() / 2;

}

double Mesh::cal_volume2(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
  Eigen::MatrixXd V2(V.rows() + 1, V.cols());
  V2.topRows(V.rows()) = V;
  V2.bottomRows(1).setZero();
  Eigen::MatrixXi T(F.rows(), 4);
  T.leftCols(3) = F;
  T.rightCols(1).setConstant(V.rows());
  Eigen::VectorXd vol;
  igl::volume(V2, T, vol);
  return std::abs(vol.sum());
}

double Mesh::cal_volume(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
  volume_total = 0.0;
  sum = 0.0;

  //#pragma omp parallel for reduction(+:volume_total) private(p0x, p0y, p0z, p1x, p1y, p1z, p2x, p2y, p2z, v321, v231, v312, v132, v213, v123, sum)
  for (int i = 0; i < numF; i++) {

    p0x = V(F(i,0),0);
    p0y = V(F(i,0),1);
    p0z = V(F(i,0),2);
    p1x = V(F(i,1),0);
    p1y = V(F(i,1),1);
    p1z = V(F(i,1),2);
    p2x = V(F(i,2),0);
    p2y = V(F(i,2),1);
    p2z = V(F(i,2),2);
    v321 = p2x*p1y*p0z;
    v231 = p1x*p2y*p0z;
    v312 = p2x*p0y*p1z;
    v132 = p0x*p2y*p1z;
    v213 = p1x*p0y*p2z;
    v123 = p0x*p1y*p2z;
    sum = (-v321+ v231+ v312-v132-v213+ v123) / 6.0;
    volume_total += sum;
  }

  return volume_total;
}


Eigen::MatrixXd Mesh::area_grad(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
  AG = -L * V;
  return AG;
}


Eigen::MatrixXd Mesh::volume_grad(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
  //Eigen::MatrixXd vol_grad;
  VG.resize(numV, 3);
  VG.setZero();

  int i, j, k;

  //#pragma omp parallel for private(j, k) shared(VG)
  for (i = 0; i < numV; i++) {
    Eigen::RowVector3d vol_ij;
    vol_ij.setZero();

    for (j = 0; j < VF[i].size(); j++) {
      k = VF[i][j];
      vol_ij += (1.0/3.0)*(dblA(k)/2.0)*F_normals.row(k);
    }

    //#pragma omp critical
    VG.row(i) = vol_ij;
  }

  //std::cout<<"vol_grad \n"<<vol_grad<<std::endl;
  return VG;
}


Eigen::MatrixXd Mesh::vertex_smoothing(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
  V_new.resize(numV, 3);
  V_new.setZero();
  igl::barycenter(V, F, V_b);

  int i, j, k;

  //#pragma omp parallel for private(j, k)
  for (i = 0; i < numV; i++) {
    Eigen::MatrixXd f_norm(VF[i].size(), 3);
    Eigen::MatrixXd V_centroid(VF[i].size(), 3);
    Eigen::VectorXd face_area(VF[i].size());
          
    double face_area_sum = 0.0;
    
    for (j = 0; j < VF[i].size(); j++) {
      k = VF[i][j];
      face_area_sum += dblA(k);
      face_area(j) = dblA(k);
      f_norm.row(j) = F_normals.row(k);
      V_centroid.row(j) = V_b.row(k);
    }

    Eigen::VectorXd sum_of_area_centroid = face_area.transpose() * V_centroid;
    Eigen::RowVector3d V_avg = sum_of_area_centroid / face_area_sum;
    Eigen::RowVector3d f_norm_sum = f_norm.colwise().sum();
    double lambda = ((V_avg.dot(f_norm_sum)) - (V.row(i).dot(f_norm_sum))) / (f_norm_sum.dot(f_norm_sum));

    V_new.row(i) = V_avg - lambda*f_norm_sum;
  }

  return V_new;
}
