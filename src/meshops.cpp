
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
#include "meshops.h"



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
