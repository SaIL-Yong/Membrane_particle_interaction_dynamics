#include <iostream>
#include <cmath>
#include "energy.h"
#include "meshops.h"
//#include "constant.h"

void Energy::compute_bendingenergy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,Eigen::MatrixXd& Force_Bending,double& total_EB){
  Eigen::MatrixXd H,HN,H_squared;
  Eigen::SparseMatrix<double> L,M,Minv;
  Eigen::VectorXd area_voronoi,Lap_H,force_density,K;
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
  EB = 2.0*0.01*(H_squared.transpose() * area_voronoi).diagonal();
  total_EB=EB.sum();
  Lap_H=Minv*(L*H_signed);
  force_density=(2.0*H_signed.array()*(H_squared-K).array() )+ Lap_H.array();
  Eigen::VectorXd vector_term=force_density.array()*area_voronoi.array();
  Force_Bending=(2*0.01)*(V_normals.array().colwise()*vector_term.array());
}

void Energy::compute_areaenergy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,Eigen::MatrixXd& Force_Area,double& area_energy){
  float radious=1.0;
  double area_target=4.0*PI*radious*radious;
  Eigen::VectorXd dblA;                                   // store face edge information
  igl::doublearea(V,F,dblA);
  double area_current = dblA.sum()/2;
  area_energy=2.0*(pow((area_current-area_target),2)/area_target);
  double scalar_term=-2.0*2.0*((area_current-area_target)/area_target);
  Mesh m;
  Eigen::MatrixXd AG =m.area_grad(V,F);
  Force_Area=scalar_term*AG;
}
void Energy::compute_volumeenergy_force(Eigen::MatrixXd V,Eigen::MatrixXi F,float reduced_volume,Eigen::MatrixXd& Force_Volume,double& volume_energy){
  float radious=1.0;
  double volume_target=reduced_volume*(4.0/3.0)*PI*pow(radious,3);
  Mesh m;
  double volume_current= m.cal_volume(V,F);
  Eigen::MatrixXd VG =m.volume_grad(V,F);
  volume_energy=5.0*(pow((volume_current-volume_target),2)/volume_target);
  double scalar_term=-2.0*5.0*((volume_current-volume_target)/volume_target);
  Force_Volume=scalar_term*VG;
}

/*double Energy::BendingEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F){
           Eigen::MatrixXd HN,H,H_squared;
           Eigen::SparseMatrix<double> L,M,Minv;
           Eigen::VectorXd area_voronoi;
           Eigen::VectorXd EB;
           double total_EB;

          igl::cotmatrix(V,F,L);
          igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
          igl::invert_diag(M,Minv);
          HN= -Minv*(L*V);
          H = HN.rowwise().norm()/2.0; //up to sign
          H_squared=H.array().square();
          area_voronoi=M.diagonal();
          EB = 2.0*0.01*(H_squared.transpose() * area_voronoi).diagonal();
          total_EB=EB.sum();
          return total_EB;
}
*/
/*double Energy::AreaEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F){
        float radious=1.0;
        double area_target=4.0*PI*radious*radious;
        Eigen::VectorXd dblA;                                   // store face edge information
        igl::doublearea(V,F,dblA);

        double area_total = dblA.sum()/2;
       // std::cout << "total area " << area_total << std::endl;
        double area_energy=2.0*(pow((area_total-area_target),2)/area_target);
        return area_energy;

}
double Energy::VolumeEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F,double reduced_volume){
        float radious=1.0;
        double volume_target=reduced_volume*(4/3)*PI*pow(radious,3);
        Mesh M1;
        double volume_total= M1.cal_volume(V,F);
        double volume_energy=1*(pow((volume_total-volume_target),2)/volume_target);
        return volume_energy;

}




/*void BendingEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F){
	igl::cotmatrix(V,F,L);
	igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
	igl::invert_diag(M,Minv);
	HN= -Minv*(L*V);
	H = HN.rowwise().norm(); //up to sign
	area_voronoi=M1.diagonal();
	EB=2*0.01*H*area_voronoi;
	total_EB1=EB.sum();

	return total_EB1;
}
//BE=BendingEnergy->BE;

/*Energy En1;
void BE(Eigen::MatrixXd V,Eigen::MatrixXi F){
	V1= En1.V;
	F1= En1.F;
	L1=En1.L;
	M1=En1.M;
	Minv1=En1.Minv;
	area_voronoi1=En1.area_voronoi;
	EB1=En1.EB;
	total_EB1=En1.total_EB;


//	V= BendingEnergy->V;
//	F=BendingEnergy->F;
//	L= BendingEnergy->L;
//	M= BendingEnergy->M;
//	Minv=BendingEnergy->Minv;
//	EB=BendingEnergy->EB;
//	total_EB=BendingEnergy->totalEB;
//	area_voronoi=BendignEnergy->area-voronoi;
//

	igl::cotmatrix(V1,F1,L1);
	igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M1);
	igl::invert_diag(M1,Minv1);
	HN1 = -Minv1*(L1*V1);
	H1 = HN1.rowwise().norm(); //up to sign
	area_voronoi1=M1.diagonal();
	EB1=2*0.01*H1*area_voronoi1;
	total_EB1=EB1.sum();

	return total_EB1;
}
*/
