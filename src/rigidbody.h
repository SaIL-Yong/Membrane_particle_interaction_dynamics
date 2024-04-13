#pragma once
#ifndef RIGIDBODY_H
#define RIGIDBODY_H
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// Assuming "meshops.h" is needed elsewhere

class RigidBody {
private:
   
    Eigen::Matrix3d moment_of_inertia = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d inverse_moment_of_inertia = Eigen::Matrix3d::Zero();
    //bool propertiesCalculated = false;
     Eigen::Quaterniond orientation;
public:
    void calculateProperties(const Eigen::MatrixXd points, double mass/*,Eigen::Matrix3d& moment_of_inertia, Eigen::Matrix3d& inverse_moment_of_inertia*/);
     // Accessor methods for moment of inertia
    Eigen::Matrix3d getMomentOfInertia() const { return moment_of_inertia; }
    Eigen::Matrix3d getInverseMomentOfInertia() const { return inverse_moment_of_inertia; }

    // Declare torqueFromForce to take force, point_of_application, and centerOfMass by reference

    void calculate_center_of_mass(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::Vector3d& center_of_mass);
    void calculate_torque(Eigen::MatrixXd force, Eigen::MatrixXd point_of_application,Eigen::Vector3d center_of_mass,Eigen::Vector3d& torque);
    void printTorque(Eigen::MatrixXd force,  Eigen::MatrixXd point_of_application, Eigen::Vector3d center_of_mass);
};
#endif // RIGIDBODY_H
/*   // M1.mesh_cal(V1, F1);
    // E1.compute_bendingenergy_force(V1, F1, Kb, Force_Bending, EnergyBending, M1);
    // E1.compute_areaenergy_force(V1, F1, Ka, area_target, Force_Area, EnergyArea, M1);
    // E1.compute_volumeenergy_force(V1, F1, Kv, volume_target, Force_Volume, EnergyVolume, M1);
    // if (parameter.random_force_flag)E1.compute_random_force(V1, gamma, kbT, mass, dt, Force_Random);
    // if(i % bondfrequency == 0 && particle_flag )igl::signed_distance(V1, V2, F2, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, signed_distance, facet_index, closest_points, normals_closest_points);


    // if (particle_flag) E1.compute_adhesion_energy_force(V1, F1, closest_points, rho, U,r_equilibrium,rc,angle_flag,
    //                                 particle_position,sigma,  Ew_t, Kw,Force_Adhesion,signed_distance, EnergyAdhesion,EnergyBias, M1);
    // //if (particle_flag) E1.compute_adhesion_energy_force(V1, F1, X0, Y0, Z0, Rp, rho, U, rc, angle_flag, particle_position, Ew_t, Kw, Force_Adhesion, EnergyAdhesion, EnergyBias, M1);
    
    // EnergyTotal = EnergyBending + EnergyArea + EnergyVolume + EnergyAdhesion + EnergyBias;
    // Force_Total = Force_Bending + Force_Area + Force_Volume + Force_Adhesion;// + Force_Random;

    // //acceleration = Force_Total/mass;
    // acceleration_half_step = Force_Total / mass;

    // //velocity_half_step = velocity + 0.5 *dt* (acceleration_half_step- gamma*velocity) + Force_Random ;

    // V1 += velocity * dt + 0.5 * acceleration_half_step * (dt * dt);

    //V1 += (velocity * dt);// + 0.5 * acceleration_half_step * (dt * dt) ;
    //Repeat the force calucaltion here
    M1.mesh_cal(V1, F1);
    E1.compute_bendingenergy_force(V1, F1, Kb, Force_Bending, EnergyBending, M1);
    E1.compute_areaenergy_force(V1, F1, Ka, area_target, Force_Area, EnergyArea, M1);
    E1.compute_volumeenergy_force(V1, F1, Kv, volume_target, Force_Volume, EnergyVolume, M1);
    if (parameter.random_force_flag)E1.compute_random_force(V1, gamma, kbT, mass, dt, Force_Random);
    if(i % bondfrequency == 0 && particle_flag){
    igl::signed_distance(V1, V2, F2, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, signed_distance, facet_index, closest_points, normals_closest_points);
   //std::cout << "Facet Index" << facet_index<<std::endl;
    std::ofstream outfile("signed_distance.txt");
   // Check if the file was successfully opened
    if (outfile.is_open()) {
    outfile << signed_distance << std::endl;
    outfile << facet_index << std::endl;
    outfile.close();
    }
    else {
    std::cout << "Error: cannot open adhesion force file." <<std::endl;
        }
    }

    if (particle_flag) E1.compute_adhesion_energy_force(V1, F1, closest_points, rho, U,r_equilibrium,rc,angle_flag,
                                    particle_position,sigma,  Ew_t, Kw,Force_Adhesion,signed_distance, EnergyAdhesion,EnergyBias, M1);
    //if (particle_flag) E1.compute_adhesion_energy_force(V1, F1, X0, Y0, Z0, Rp, rho, U, rc, angle_flag, particle_position, Ew_t, Kw, Force_Adhesion, EnergyAdhesion, EnergyBias, M1);
    EnergyTotal = EnergyBending + EnergyArea + EnergyVolume + EnergyAdhesion + EnergyBias;
    Force_Total = Force_Bending + Force_Area + Force_Volume + Force_Adhesion;// + Force_Random;
    force_residual = Force_Total.norm();


    // acceleration = Force_Total / mass;

    // velocity = 0.5 * (acceleration + acceleration_half_step) * dt;
    // Update velocities with average acceleration
    //velocity = velocity_half_step + 0.5 *dt*(acceleration - (gamma * velocity_half_step)) + Force_Random ;

    //ForcesonParticleVertices
    //if (particle_flag) E1.redistributeAdhesionForce(V2,F2,closest_points, Force_Adhesion, facet_index,ForcesOnVertices); 

    
    //body.calculate_center_of_mass(V2,F2,COM);
    //body.printTorque(ForcesOnVertices, closest_points, COM);
    //std::cout << "Force_Total" << Force_Total << std::endl;
    

    /*  Rigid Body Calculations 
    // Calculate the properties
    //RigidBody::MeshProperties props = RigidBody::calculate_properties(V2, mass);
    Eigen::Matrix3d moment_of_inertia;
    Eigen::Matrix3d inverse_moment_of_inertia;
    R1.calculateProperties(V2,mass,moment_of_inertia, inverse_moment_of_inertia);
    // Use the calculated properties
    std::cout << "Inertia Tensor: \n" << moment_of_inertia << std::endl;   
    R1.calculate_center_of_mass(V2,F2,center_of_mass);
    R1.calculate_torque(V2,ForcesOnVertices,center_of_mass, torque);
    std::cout << "Center of Mass: \n" << center_of_mass.transpose() << std::endl;
    std::cout << "Torque: \n" << torque.transpose() << std::endl;
    
    */