#include <iostream>
#include <cmath>
#include "energy.h"

void Energy::compute_bendingenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Kb, Eigen::MatrixXd& Force_Bending, double& bending_energy, Mesh m)
{
  bending_energy = 0.0;
  Force_Bending.setZero();
  EB = 2.0 * Kb * (m.H_squared.transpose() * m.area_voronoi).diagonal();
  bending_energy = EB.sum();

  Lap_H = m.Minv * (m.L * m.H_signed);
  force_density = (2.0 * m.H_signed.array() * (m.H_squared - m.K).array()) + Lap_H.array();
  vector_term = force_density.array() * m.area_voronoi.array();
  Force_Bending = (2.0 *Kb)*(m.V_normals.array().colwise() * vector_term.array());
}


void Energy::compute_areaenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Ka, double area_target, Eigen::MatrixXd& Force_Area, double& area_energy, Mesh m)
{
  area_energy = 0.0;
  Force_Area.setZero();
  da = m.area_total - area_target;
  area_energy = Ka * da * da / area_target;
  scalar_term = -2.0 * Ka * da / area_target;
  AG = m.area_grad(V, F);
  Force_Area = scalar_term * AG;
}


void Energy::compute_volumeenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Kv, double volume_target, Eigen::MatrixXd& Force_Volume, double& volume_energy, Mesh m)
{
  volume_energy = 0.0;
  Force_Volume.setZero();
  if (std::abs(Kv) > EPS) {
    VG = m.volume_grad(V, F);
    dv = m.volume_total - volume_target;
    volume_energy = Kv * dv * dv / volume_target;
    scalar_term = -2.0 * Kv * dv / volume_target;
    Force_Volume = scalar_term * VG;
  }
}
void Energy::compute_adhesion_energy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd V_particle,double rho, double U,double r_equilibrium,double rc,
                                           int angle_flag,int particle_position,double sigma,double Ew_t, double Kw,
                                           Eigen::MatrixXd& Force_Adhesion,Eigen::MatrixXd signed_distance, double& EnergyAdhesion,double& EnergyBias,  Mesh m)
{
  Force_Adhesion.setZero();
  coefficient.resize(V.rows());
  coefficient_derivative_x.resize(V.rows());
  coefficient_derivative_y.resize(V.rows());
  coefficient_derivative_z.resize(V.rows());
  distance.resize(V.rows());
  //comvec.resize(V.rows(),3);
  // lennard_jones_force.resize(V.rows(),3);
  // lennard_jones_potential.resize(V.rows());
  //dc.resize(V.rows());
  //Mod_Bias.resize(V.rows());
  coefficient.setZero();
  coefficient_derivative_x.setZero();
  coefficient_derivative_y.setZero();
  coefficient_derivative_z.setZero();
  distance.setZero();
  lennard_jones_force.setZero();
  lennard_jones_potential.setZero();
  
  //comvec = V.rowwise() - particle_center;
  comvec = V-V_particle;


  // can this loop be replaced by the eigen operation

  for (int i = 0; i < V.rows(); i++) {
    dc = signed_distance(i);
    //std::cout<<"signed distance"<<dc<<std::endl;

    

    if (std::abs(dc) > rc) continue;

    // // angle between connecting vector and vertex normal
    angle = acos((comvec.row(i)).dot(m.V_normals.row(i)) / (comvec.row(i).norm() * m.V_normals.row(i).norm()));
    

    if (angle_flag) {
       if (dc > 0 && angle <= 0.5*PI) continue;
      //  if (dc< 0 && angle >= 0.5*PI) continue;
      }


    //distance(bonds[i].first) = sqrt((V(i,0)-particle_centroid.x())*(V(i,0)-particle_centroid.x())+(V(i,1)-particle_centroid.y())*(V(i,1)-particle_centroid.y())
    //+(V(i,2)-particle_centroid.z())*(V(i,2)-particle_centroid.z()));
    dx=(V(i,0)-(V_particle(i,0)));
    dy=(V(i,1)-(V_particle(i,1)));
    dz=(V(i,2)-(V_particle(i,2)));

    //dc(bonds[i].first) =  sqrt(dx*dx + dy*dy +dz*dz) - r_equilibrium;
    // r_ij=(V.row((bonds[i].first)) - V_particle.row((bonds[i].second)));
    // dc = r_ij.norm();
    // dc_mag=r_ij.squaredNorm();
    //r_ij_transpose= r_ij.transpose();
    ////std::cout<<"r_ij vector form:"<<r_ij_transpose<<std::endl;

    coefficient(i) = U * (exp(-(2.0*dc)/rho) - 2.0*(exp(-dc/rho)));
    coefficient_derivative_x(i) = (U/(dc*rho))
                                *(-exp(-(2.0*dc)/rho) + exp(-dc/rho)) * 2.0 * dx;
    coefficient_derivative_y(i) = (U/(dc*rho))
                                *(-exp(-(2.0*dc)/rho) + exp(-dc/rho)) * 2.0 * dy;
    coefficient_derivative_z(i) = (U/(dc*rho))
                                *(-exp(-(2.0*dc)/rho) + exp(-dc/rho)) * 2.0 * dz;

    // s_by_r6 = pow(sigma / dc, 6.0);
    // if (dc>= (pow(2,(1.0/6.0))*sigma)){
    //   lennard_jones_potential(i) = 0.0;
    //   lennard_jones_force.row(i).setZero();     
    // }else{
    // lennard_jones_potential.row(bonds[i].first) << 4.0 * epsilon * (pow(s_by_r6, 2.0) - s_by_r6);
    // lennard_jones_force.row(bonds[i].first) += 24.0 * (epsilon / dc_mag) * (2.0 * pow(s_by_r6, 2.0) - s_by_r6) * (r_ij_transpose/dc_mag );                    
    // } 


  }
  //std::cout<<coefficient<<std::endl;
  coefficient_of_derivative.resize(V.rows(), 3);
  coefficient_of_derivative.col(0)=coefficient_derivative_x.transpose();
  coefficient_of_derivative.col(1)=coefficient_derivative_y.transpose();
  coefficient_of_derivative.col(2)=coefficient_derivative_z.transpose();

  // lennard_jones_force.resize(V.rows(), 3);
  // lennard_jones_potential.resize(V.rows());
 

  
  First_Term = -(AG.array().colwise()*coefficient.array());
  Second_Term = -(coefficient_of_derivative.array().colwise()*m.area_voronoi.array());

  Sum = First_Term + Second_Term;
  Ead = coefficient.array()*m.area_voronoi.array();

  EnergyAdhesion = Ead.sum();

  // if (epsilon > 0.0) {

  // Force_Adhesion = Sum+lennard_jones_force;
  // EnergyAdhesion = Ead.sum()+lennard_jones_potential.sum();
  // }else {
  // EnergyAdhesion = Ead.sum();
  // Force_Adhesion = Sum;
  // }


  if (std::abs(Kw) > EPS) {
    dEw = EnergyAdhesion - Ew_t;
    EnergyBias = 0.5 * Kw * dEw * dEw;
    //Force_Biased = Kw * dEw * (Sum.array().colwise() * Mod_Bias.array());
    Force_Biased = Kw * dEw * Sum;
    Force_Adhesion = Sum + Force_Biased;
  } else Force_Adhesion = Sum; 



}


//   // can this loop be replaced by the eigen operation?

//   for (int i = 0; i <bonds.size(); i++) {
//     //distance(bonds[i].first) = sqrt((V(i,0)-particle_centroid.x())*(V(i,0)-particle_centroid.x())+(V(i,1)-particle_centroid.y())*(V(i,1)-particle_centroid.y())
//     //+(V(i,2)-particle_centroid.z())*(V(i,2)-particle_centroid.z()));
//     dx=(V(bonds[i].first,0)-(V_particle(bonds[i].second,0)));
//     dy=(V(bonds[i].first,1)-(V_particle(bonds[i].second,1)));
//     dz=(V(bonds[i].first,2)-(V_particle(bonds[i].second,2)));

//     //dc(bonds[i].first) =  sqrt(dx*dx + dy*dy +dz*dz) - r_equilibrium;
//     r_ij=(V.row((bonds[i].first)) - V_particle.row((bonds[i].second)));
//     dc = r_ij.norm();
//     dc_mag=r_ij.squaredNorm();
//     //r_ij_transpose= r_ij.transpose();
//     ////std::cout<<"r_ij vector form:"<<r_ij_transpose<<std::endl;

//     coefficient.row(bonds[i].first) << U * (exp(-(2.0*dc)/rho) - 2.0*(exp(-dc/rho)));
//     coefficient_derivative_x.row(bonds[i].first) << (U/(dc*rho))
//                                 *(-exp(-(2.0*dc)/rho) + exp(-dc/rho)) * 2.0 * dx;
//     coefficient_derivative_y.row(bonds[i].first) << (U/(dc*rho))
//                                 *(-exp(-(2.0*dc)/rho) + exp(-dc/rho)) * 2.0 * dy;
//     coefficient_derivative_z.row(bonds[i].first) << (U/(dc*rho))
//                                 *(-exp(-(2.0*dc)/rho) + exp(-dc/rho)) * 2.0 * dz;

    
//     s_by_r6 = pow(sigma / dc, 6.0);
//     //double pot = 4.0 * eps * (pow(s_by_r6, 2.0) - s_by_r6);
//     if (dc>= (pow(2,(1.0/6.0))*sigma)){
//       lennard_jones_potential.row(bonds[i].first) << 0.0;
//       lennard_jones_force.row(bonds[i].first).setZero();     
//     }else{
//     lennard_jones_potential.row(bonds[i].first) << 4.0 * epsilon * (pow(s_by_r6, 2.0) - s_by_r6);
//     lennard_jones_force.row(bonds[i].first) += 24.0 * (epsilon / dc_mag) * (2.0 * pow(s_by_r6, 2.0) - s_by_r6) * (r_ij_transpose/dc_mag );                    
//     }           

//   }
//   //std::cout<<coefficient<<std::endl;
//   coefficient_of_derivative.resize(V.rows(), 3);
//   coefficient_of_derivative.col(0)=coefficient_derivative_x.transpose();
//   coefficient_of_derivative.col(1)=coefficient_derivative_y.transpose();
//   coefficient_of_derivative.col(2)=coefficient_derivative_z.transpose();
//   lennard_jones_force.resize(V.rows(), 3);
//   lennard_jones_potential.resize(V.rows());
  
//   First_Term = -(AG.array().colwise()*coefficient.array());
//   Second_Term = -(coefficient_of_derivative.array().colwise()*m.area_voronoi.array());

//   Sum = First_Term + Second_Term;
//   Ead = coefficient.array()*m.area_voronoi.array();


//   //std::cout<<"lennard jone potential"<<lennard_jones_potential<<std::endl;
//   //Force_Adhesion = Sum;

//   if (epsilon > 0.0) {

//   Force_Adhesion = Sum+lennard_jones_force;
//   EnergyAdhesion = Ead.sum()+lennard_jones_potential.sum();
//   }else {
//   EnergyAdhesion = Ead.sum();
//   Force_Adhesion = Sum;
//   }

//   if (std::abs(Kw) > EPS) {
//     dEw = EnergyAdhesion - Ew_t;
//     EnergyBias = 0.5 * Kw * dEw * dEw;
//     //Force_Biased = Kw * dEw * (Sum.array().colwise() * Mod_Bias.array());
//     Force_Biased = Kw * dEw * Sum;
//     Force_Adhesion = Sum + Force_Biased;
//   } else Force_Adhesion = Sum; 

// }

// void Energy::compute_adhesion_energy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double X, double Y, double Z,
//                                            double Rp, double rho, double U, double rc, int angle_flag, int particle_position, double Ew_t, double Kw,
//                                            Eigen::MatrixXd& Force_Adhesion, double& EnergyAdhesion, double& EnergyBias, Mesh m)
// {
//   Force_Adhesion.setZero();
//   coefficient.resize(V.rows());
//   coefficient_derivative_x.resize(V.rows());
//   coefficient_derivative_y.resize(V.rows());
//   coefficient_derivative_z.resize(V.rows());
//   distance.resize(V.rows());
//   dc.resize(V.rows());
//   //Mod_Bias.resize(V.rows());
//   coefficient.setZero();
//   coefficient_derivative_x.setZero();
//   coefficient_derivative_y.setZero();
//   coefficient_derivative_z.setZero();
//   distance.setZero();
//   dc.setZero();
//   //Mod_Bias.setZero();
//   EnergyAdhesion = 0.0;
//   EnergyBias = 0.0;

//   particle_center<<X, Y, Z;
//   // vector connecting cetner of the particle and to the vertices;
//   comvec = V.rowwise() - particle_center;

//   // can this loop be replaced by the eigen operation?
//   for (int i = 0; i < V.rows(); i++) {
//     distance(i) = sqrt((V(i,0)-X)*(V(i,0)-X)+(V(i,1)-Y)*(V(i,1)-Y)+(V(i,2)-Z)*(V(i,2)-Z));
//     dc(i) = distance(i) - Rp;

//     // angle between connecting vector and vertex normal
//     angle = acos((comvec.row(i)).dot(m.V_normals.row(i)) / (comvec.row(i).norm() * m.V_normals.row(i).norm()));
    
//     if (std::abs(dc(i)) > rc) continue;
//     if (angle_flag) {
//       if (particle_position > 0 && angle <= 0.5*PI) continue;
//       if (particle_position < 0 && angle >= 0.5*PI) continue;
//     }

//     coefficient(i) = U * (exp(-(2.0*dc(i))/rho) - 2.0*exp(-dc(i)/rho));
//     coefficient_derivative_x(i) = (U/(distance(i)*rho))
//                                 *(-exp(-(2.0*dc(i))/rho) + exp(-dc(i)/rho)) * 2.0 * (V(i,0)-X);
//     coefficient_derivative_y(i) = (U/(distance(i)*rho))
//                                 *(-exp(-(2.0*dc(i))/rho) + exp(-dc(i)/rho)) * 2.0 * (V(i,1)-Y);
//     coefficient_derivative_z(i) = (U/(distance(i)*rho))
//                                 *(-exp(-(2.0*dc(i))/rho) + exp(-dc(i)/rho)) * 2.0 * (V(i,2)-Z);

//     // if (dc(i) > EPS && std::abs(Kw) > EPS) Mod_Bias(i) = 1.0;
//   }

//   coefficient_of_derivative.resize(V.rows(), 3);
//   coefficient_of_derivative.col(0)=coefficient_derivative_x.transpose();
//   coefficient_of_derivative.col(1)=coefficient_derivative_y.transpose();
//   coefficient_of_derivative.col(2)=coefficient_derivative_z.transpose();
  
//   First_Term = -(AG.array().colwise()*coefficient.array());
//   Second_Term = -(coefficient_of_derivative.array().colwise()*m.area_voronoi.array());

//   Sum = First_Term + Second_Term;
//   Ead = coefficient.array()*m.area_voronoi.array();

//   EnergyAdhesion = Ead.sum();

//   if (std::abs(Kw) > EPS) {
//     dEw = EnergyAdhesion - Ew_t;
//     EnergyBias = 0.5 * Kw * dEw * dEw;
//     //Force_Biased = Kw * dEw * (Sum.array().colwise() * Mod_Bias.array());
//     Force_Biased = Kw * dEw * Sum;
//     Force_Adhesion = Sum + Force_Biased;
//   } else Force_Adhesion = Sum; 
// }