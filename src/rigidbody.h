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
    Eigen::Vector3d center_of_mass;
    Eigen::Vector3d particle_acceleration;
    // //bool propertiesCalculated = false;
    Eigen::Quaterniond orientation;
    Eigen::Matrix3d rotation_matrix;
public:
    void calculate_properties(Eigen::MatrixXd points, double mass/*,Eigen::Matrix3d& moment_of_inertia, Eigen::Matrix3d& inverse_moment_of_inertia*/);
     // Accessor methods for moment of inertia
    Eigen::Vector3d getCenterOfMass() const { return center_of_mass; }
    Eigen::Matrix3d getMomentOfInertia() const { return moment_of_inertia; }
    Eigen::Matrix3d getInverseMomentOfInertia() const { return inverse_moment_of_inertia; }

    // Declare torqueFromForce to take force, point_of_application, and centerOfMass by reference

    void calculate_center_of_mass(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::Vector3d& center_of_mass);
    void calculate_torque(Eigen::MatrixXd force, Eigen::MatrixXd point_of_application,Eigen::Vector3d center_of_mass,Eigen::Vector3d& torque);
    void angular_momentum (Eigen::Vector3d torque,double dt, Eigen::Vector3d& angular_momentum);
    void calculate_omega(Eigen::Vector3d angular_momentum, Eigen::Matrix3d rot_mat,Eigen::Vector3d idiag,Eigen::Vector3d& angular_velocity);
    void update_quaternion(Eigen::Quaterniond current_quaternion, Eigen::Vector3d angular_velocity, double dt, Eigen::Quaterniond& new_quaternion);
    void update_vertex_position(Eigen::MatrixXd& V, Eigen::MatrixXd forces,  Eigen::Quaterniond quaternion, double dt, Eigen::Vector3d& particle_velocity);

    void diagonalize_inertia_tensor(Eigen::Matrix3d inertia_tensor, Eigen::Matrix3d& principal_axes,Eigen::Matrix3d& principal_moments);
    void exyz_to_q(Eigen::Matrix3d R ,Eigen::Quaterniond& quat);
    //void printTorque(Eigen::MatrixXd force,  Eigen::MatrixXd point_of_application, Eigen::Vector3d center_of_mass);
    // Declare a method to calculate the orientation of the rigid body
    Eigen::Vector3d  ang_mom,omega; 
    

};
#endif // RIGIDBODY_H
