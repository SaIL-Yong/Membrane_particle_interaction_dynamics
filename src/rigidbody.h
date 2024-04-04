#pragma once
#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include "meshops.h"
// Assuming "meshops.h" is needed elsewhere

class RigidBody {
private:
   
    Eigen::Matrix3d moment_of_inertia = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d inverse_moment_of_inertia = Eigen::Matrix3d::Zero();
    //bool propertiesCalculated = false;
public:
    void calculateProperties(const Eigen::MatrixXd points, double mass,Eigen::Matrix3d& moment_of_inertia, Eigen::Matrix3d& inverse_moment_of_inertia);
    //Eigen::Matrix3d getMomentOfInertia() const { return momentOfInertia; }
    //Eigen::Matrix3d getInverseMomentOfInertia();

    // Declare torqueFromForce to take force, point_of_application, and centerOfMass by reference

    void calculate_center_of_mass(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::Vector3d& center_of_mass);
    void calculate_torque(Eigen::MatrixXd force, Eigen::MatrixXd point_of_application,Eigen::Vector3d center_of_mass,Eigen::Vector3d& torque);
};
#endif // RIGIDBODY_H
