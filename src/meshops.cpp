#include "energy.h"
#include "meshops.h"
#include "parameters.h"

double Mesh::cal_volume(Eigen::MatrixXd V,Eigen::MatrixXi F){
          double Volume_total;
	  double sum=0;
          int numF = F.rows();
          for (int i = 0; i < numF; i++){

          double p0x=V(F(i,0),0);
          double p0y=V(F(i,0),1);
          double p0z=V(F(i,0),2);
          double p1x=V(F(i,1),0);
          double p1y=V(F(i,1),1);
          double p1z=V(F(i,1),2);
          double p2x=V(F(i,2),0);
          double p2y=V(F(i,2),1);
          double p2z=V(F(i,2),2);
          double v321= p2x*p1y*p0z;
          double v231= p1x*p2y*p0z;
          double v312= p2x*p0y*p1z;
          double v132= p0x*p2y*p1z;
          double v213= p1x*p0y*p2z;
          double v123= p0x*p1y*p2z;
          sum=(-v321+ v231+ v312-v132-v213+ v123) / 6.0;
          Volume_total+=sum;

          }

          return Volume_total;
}


Eigen::MatrixXd Mesh::area_grad(Eigen::MatrixXd V,Eigen::MatrixXi F){
           Eigen::SparseMatrix<double> L;
           Eigen::MatrixXd AG;
           igl::cotmatrix(V,F,L);
           AG=-L*V;
           return AG;
}
Eigen::MatrixXd Mesh::volume_grad(Eigen::MatrixXd V,Eigen::MatrixXi F){
           //Eigen::MatrixXd vol_grad;
           Eigen::MatrixXd F_normals;
           Eigen::VectorXd dblA;
           int numV=V.rows();
           std::vector<std::vector<double> > VF;
           std::vector<std::vector<double> > VFi;
           igl::vertex_triangle_adjacency(numV,F,VF, VFi);
           igl::per_face_normals(V,F,F_normals);

           igl::doublearea(V,F,dblA);
           Eigen::MatrixXd vol_grad(numV,3);
           for (int i=0; i<numV; i++){
             Eigen::RowVectorXd vol_ij {{0.0, 0.0, 0.0}};

             for (int j=0;j<VF[i].size();j++){
               int k=VF[i][j];
               vol_ij += (1.0/3.0)*(dblA(k)/2.0)*F_normals.row(k);
             }
           vol_grad.row(i)=vol_ij;
           }
           //std::cout<<"vol_grad \n"<<vol_grad<<std::endl;
          return vol_grad;
}
Eigen::MatrixXd Mesh::vertex_smoothing(Eigen::MatrixXd V,Eigen::MatrixXi F){
        int numV = V.rows();
        Eigen::MatrixXd V_b,F_normals;
        Eigen::MatrixXd V_new(numV,3);
        Eigen::VectorXd dblA;
        std::vector<std::vector<double> > VF;
        std::vector<std::vector<double> > VFi;
        igl::vertex_triangle_adjacency(V,F,VF, VFi);
        igl::doublearea(V,F,dblA);
        igl::barycenter(V,F,V_b);
        igl::per_face_normals(V,F,F_normals);  //Eigen::ArrayXf V_b = Eigen::ArrayXf::Zero(VF[0].size());
        double area_avg;
        for (int i=0; i<numV; i++){
          Eigen::MatrixXd f_norm(VF[i].size(),3);
          Eigen::MatrixXd V_centroid(VF[i].size(),3);
          Eigen::VectorXd face_area(VF[i].size());
          double face_area_sum=0;
          for (int j=0;j<VF[i].size();j++){
            int k=VF[i][j];
            //std::cout<<"facearea \n "<< dblA(k)<<std::endl;
            face_area_sum += (dblA(k));
            face_area(j)= dblA(k);
            f_norm.row(j)=F_normals.row(k);
            V_centroid.row(j)=V_b.row(k);
          }
          Eigen::VectorXd sum_of_area_centroid=face_area.transpose()*V_centroid;
          Eigen::VectorXd V_avg= sum_of_area_centroid/face_area_sum;
          Eigen::VectorXd f_norm_sum=f_norm.colwise().sum();//fsum in python
          //std::cout<<"f_sum \n "<< V.row(i).dot(f_norm_sum)<<std::endl;
          double lamda=((V_avg.dot(f_norm_sum))-(V.row(i).dot(f_norm_sum)))/(f_norm_sum.dot(f_norm_sum));

          V_new.row(i)=V_avg-lamda*f_norm_sum;

        }
        return V_new;

}
