#include <iostream>
#include <cmath>
#include "energy.h"
#include "meshops.h"
#include "parameters.h"

extern Parameter parameter;

void Energy::compute_bendingenergy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,Eigen::MatrixXd& Force_Bending,double& total_EB){
  bending_modulus=parameter.Kb;
  Eigen::MatrixXd H,HN,H_squared;
  Eigen::SparseMatrix<double> L,M,Minv;
  Eigen::VectorXd Lap_H,force_density,K;
  Eigen::VectorXd EB;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  igl::invert_diag(M,Minv);
  igl::gaussian_curvature(V,F,K);
  K = (Minv*K).eval();
  HN= -Minv*(L*V)/2.0;
  H = HN.rowwise().norm(); //up to sign
  Eigen::MatrixXd V_normals;
  igl::per_vertex_normals(V,F,V_normals);
  Eigen::MatrixXd H_X_N= (HN.array()*V_normals.array());
  Eigen::VectorXd abc=H_X_N.rowwise().sum();
  Eigen::VectorXd sign_of_H= abc.array().sign();
  Eigen::VectorXd H_signed=H.array()*sign_of_H.array();
  H_squared=H.array().square();
  area_voronoi=M.diagonal();

  EB = 2.0*bending_modulus*(H_squared.transpose() * area_voronoi).diagonal();
  total_EB=EB.sum();
  Lap_H=Minv*(L*H_signed);
  force_density=(2.0*H_signed.array()*(H_squared-K).array() )+ Lap_H.array();
  Eigen::VectorXd vector_term=force_density.array()*area_voronoi.array();
  Force_Bending=(2*bending_modulus)*(V_normals.array().colwise()*vector_term.array());
}

void Energy::compute_areaenergy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,Eigen::MatrixXd& Force_Area,double& area_energy){
  double area_modulus=parameter.Ka;
  float radious=1.0;
  double area_target=4.0*PI*radious*radious;
  Eigen::VectorXd dblA;                                   // store face edge information
  igl::doublearea(V,F,dblA);
  double area_current = dblA.sum()/2;
  area_energy=area_modulus*(pow((area_current-area_target),2)/area_target);
  double scalar_term=-2.0*area_modulus*((area_current-area_target)/area_target);
  Mesh m;
  AG =m.area_grad(V,F);
  Force_Area=(scalar_term*AG);


}
void Energy::compute_volumeenergy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,float reduced_volume,Eigen::MatrixXd& Force_Volume,double& volume_energy){
  double volume_modulus=parameter.Kv;
  float radious=1.0;
  double volume_target=reduced_volume*(4.0/3.0)*PI*pow(radious,3);
  Mesh m;
  double volume_current= m.cal_volume(V,F);
  Eigen::MatrixXd VG =m.volume_grad(V,F);
  volume_energy=volume_modulus*(pow((volume_current-volume_target),2)/volume_target);
  double scalar_term=-2.0*volume_modulus*((volume_current-volume_target)/volume_target);
  Force_Volume=scalar_term*VG;
}

void Energy::compute_adhesion_energy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,float X,float Y,float Z,float rp,float rho,float u,float U,float rc,Eigen::MatrixXd& Force_Adhesion,double& EnergyAdhesion){
    // float rp=parameter.particle_radious;
    // float u=parameter.adhesion_strength;
    // float rho=(parameter.potential_range)*rp;
    // float U=(bending_modulus*u)/(pow(rp,2)) ;
    //float X=1.0+rp+rho*1,Y=0.0,Z=0.0;
    Eigen::VectorXd coefficient_derivative_x(V.rows()),coefficient_derivative_y(V.rows()),
    coefficient_derivative_z(V.rows()),coefficient(V.rows()),distance(V.rows()),dc(V.rows());
    float tol=1e-10;
    for (int i=0; i<V.rows(); i++){
      distance(i)=sqrt((V(i,0)-X)*(V(i,0)-X)+(V(i,1)-Y)*(V(i,1)-Y)+(V(i,2)-Z)*(V(i,2)-Z));
      dc(i)=distance(i)-rp;
      coefficient(i)=U*(exp(-(2.0*dc(i))/rho) - 2.0*exp(-dc(i)/rho));
      coefficient_derivative_x(i)  = (U/(distance(i)*rho))
      *(-exp(-(2.0*dc(i))/rho) +  exp(-dc(i)/rho)) * 2.0 * (V(i,0)-X);
      coefficient_derivative_y(i)  = (U/(distance(i)*rho))
      *(-exp(-(2.0*dc(i))/rho) +  exp(-dc(i)/rho)) * 2.0 * (V(i,1)-Y);
      coefficient_derivative_z(i)  = (U/(distance(i)*rho))
      *(-exp(-(2.0*dc(i))/rho) +  exp(-dc(i)/rho)) * 2.0 * (V(i,2)-Z);
      if (abs(coefficient(i))<tol){
        coefficient(i)=0;
      }
      if (abs(dc(i))>rc || abs(coefficient_derivative_x(i))<tol){
        coefficient_derivative_x(i)  = 0;
      }
      if (abs(dc(i))>rc || abs(coefficient_derivative_y(i))<tol){
        coefficient_derivative_y(i)  = 0;
      }
      if (abs(dc(i))>rc || abs(coefficient_derivative_z(i))<tol){
        coefficient_derivative_z(i)  = 0;
      }
}
    Eigen::MatrixXd coefficient_of_derivative(V.rows(),3);
    coefficient_of_derivative.col(0)=coefficient_derivative_x.transpose();
    coefficient_of_derivative.col(1)=coefficient_derivative_y.transpose();
    coefficient_of_derivative.col(2)=coefficient_derivative_z.transpose();
    Eigen::MatrixXd First_Term = -(AG.array().colwise()*coefficient.array());
    Eigen::MatrixXd Second_Term= -(coefficient_of_derivative.array().colwise()*area_voronoi.array());

    Force_Adhesion= First_Term +Second_Term;
    //std::cout<<"coefficient :"<<Force_Adhesion<<std::endl;
    Eigen::VectorXd Ead= coefficient.array()*area_voronoi.array();
    EnergyAdhesion=Ead.sum();
    // std::cout<<"Adhesion_Energy::"<<EnergyAdhesion<<std::endl;
    // std::cout<<"Adhesion_Energy::"<<Force_Adhesion<<std::endl;

}
