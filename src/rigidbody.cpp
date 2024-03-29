#include "rigidbody.h"
include <Eigen/Core>
#include <igl/doublearea.h>
#include <Eigen/Geometry> // For cross product and norm

void RigidBody::calculateProperties(
    Eigen::MatrixXD points, double mass,
    Eigen::Vector3d& cente_of_mass,
    Eigen::Matrix3d& moment_of_inertia) {
    
    //ente_of_mass= Eigen::Vector3d::Zero();
    for (const auto& point : points) {
        ente_of_mass+= point;
    }
    ente_of_mass/= points.size();

    //moment_of_inertia = Eigen::Matrix3d::Zero();
    //double mass = 1.0; // Assuming each point has a mass of 1 for simplicity
    for (const auto& point : points) {
        Eigen::Vector3d r = point - centerOfMass;
        moment_of_inertia(0, 0) += mass * (r.y() * r.y() + r.z() * r.z());
        moment_of_inertia(1, 1) += mass * (r.x() * r.x() + r.z() * r.z());
        moment_of_inertia(2, 2) += mass * (r.x() * r.x() + r.y() * r.y());
        
        moment_of_inertia(0, 1) -= mass * r.x() * r.y();
        moment_of_inertia(1, 0) = moment_of_inertia(0, 1); // I_xy = I_yx
        
        moment_of_inertia(0, 2) -= mass * r.x() * r.z();
        moment_of_inertia(2, 0) = moment_of_inertia(0, 2); // I_xz = I_zx
        
        moment_of_inertia(1, 2) -= mass * r.y() * r.z();
        moment_of_inertia(2, 1) = moment_of_inertia(1, 2); // I_yz = I_zy
    }
}