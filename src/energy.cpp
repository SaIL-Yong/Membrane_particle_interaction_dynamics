
#include <iostream>
#include <cmath>
#include <igl/readOFF.h>
#include <igl/per_face_normals.h>
#include "igl/doublearea.h"
#include "igl/edges.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include "energy.h"
//#include "constant.h"

double Energy::BendingEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F){ 
           Eigen::MatrixXd HN,H,H_squared;
           Eigen::SparseMatrix<double> L,M,Minv;
           Eigen::VectorXd area_voronoi;
           Eigen::VectorXd EB;
           double total_EB;
           
           
          igl::cotmatrix(V,F,L);
          igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
          igl::invert_diag(M,Minv);
          HN= -Minv*(L*V);
          H = HN.rowwise().norm(); //up to sign
          H_squared=H.array().square();
          area_voronoi=M.diagonal();
          EB = 0.01*(H_squared.transpose() * area_voronoi).diagonal();
          total_EB=EB.sum();
          return total_EB;


}

double Energy::AreaEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F){ 
        double radious=1.0; 
        double area_target=4*PI*radious*radious;
        Eigen::VectorXd dblA;                                   // store face edge information
        igl::doublearea(V,F,dblA);                  
    
        double area_total = dblA.sum()/2;
       // std::cout << "total area " << area_total << std::endl;        
        double area_energy=2*(pow((area_total-area_target),2)/area_target);
        return area_energy;

}



/*double Energy::BendingEnergy(){
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    igl::invert_diag(M,Minv);
    HN= -Minv*(L*V);
    H = HN.rowwise().norm(); //up to sign
    area_voronoi=M.diagonal();
    EB = 0.01*(H.transpose() * area_voronoi).diagonal();
    total_EB=EB.sum();
    
    return total_Eb;
    //std::cout << "total_bending_energy" <<total_EB<< std::endl;


}


/*
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd HN,H;
Eigen::SparseMatrix<double> L,M,Minv;
Eigen::VectorXd area_voronoi;
Eigen::VectorXd EB;
double total_EB;
*/

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

