#include "rigidbody.h"
#include "meshops.h"
#include <Eigen/Core>
#include <igl/doublearea.h>
#include <Eigen/Geometry> // For cross product and norm

void RigidBody::calculateProperties(const Eigen::MatrixXd points, double mass/*,Eigen::Matrix3d& moment_of_inertia, Eigen::Matrix3d& inverse_moment_of_inertia*/) {  
    Eigen::Vector3d center_of_mass = Eigen::Vector3d::Zero();
    //int numPoints = points.rows(); 
// Assuming 'points' is of type Eigen::MatrixXd (or Eigen::MatrixXf) and each row is a point
    center_of_mass = points.colwise().mean();

    // Calculate the moment of inertia
    
    for (int i = 0; i < points.rows(); ++i) {
        Eigen::Vector3d r = points.row(i).transpose() - center_of_mass; // Corrected vector subtraction
        moment_of_inertia(0, 0) += mass * (r.y() * r.y() + r.z() * r.z());
        moment_of_inertia(1, 1) += mass * (r.x() * r.x() + r.z() * r.z());
        moment_of_inertia(2, 2) += mass * (r.x() * r.x() + r.y() * r.y());

        moment_of_inertia(0, 1) -= mass * r.x() * r.y();
        moment_of_inertia(1, 0) = moment_of_inertia(0, 1); // Symmetry in the inertia tensor
        moment_of_inertia(0, 2) -= mass * r.x() * r.z();
        moment_of_inertia(2, 0) = moment_of_inertia(0, 2); // Symmetry in the inertia tensor
        moment_of_inertia(1, 2) -= mass * r.y() * r.z();
        moment_of_inertia(2, 1) = moment_of_inertia(1, 2); // Symmetry in the inertia tensor
    }


    if (moment_of_inertia.determinant() != 0) {
            inverse_moment_of_inertia = moment_of_inertia.inverse();
        } else {
            std::cerr << "Warning: Moment of inertia matrix is not invertible." << std::endl;
    }
 // Return the calculated properties
}

void RigidBody::calculate_center_of_mass(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::Vector3d& center_of_mass) {
     center_of_mass = V.colwise().mean();

}
void RigidBody::calculate_torque(Eigen::MatrixXd force, Eigen::MatrixXd point_of_application,Eigen::Vector3d center_of_mass,Eigen::Vector3d& torque) {
        torque.setZero();
        for (int i = 0; i < force.rows(); ++i) {
        Eigen::Vector3d forceVec = force.row(i).transpose();
        // Correctly extracting a point of application vector
        Eigen::Vector3d pointVec = point_of_application.row(i).transpose();

        // Calculate the vector from the center of mass to the point of application
        Eigen::Vector3d r = pointVec - center_of_mass;

        // Calculate the torque for this force and add it to the total
        torque += r.cross(forceVec);
        }
}

void RigidBody::printTorque( Eigen::MatrixXd force,  Eigen::MatrixXd point_of_application,  Eigen::Vector3d center_of_mass) {
    Eigen::Vector3d torque;
    calculate_torque(force, point_of_application, center_of_mass, torque);
    std::cout << "Torque: " << torque.transpose() << std::endl;
}
// Eigen::Quaterniond RigidBody::calculateOrientation(Eigen::MatrixXd forces, Eigen::MatrixXd point_of_application,Eigen::Vector3d torque,  double dt,  Eigen::Quaterniond& current_orientation) {
//     // Calculate the net torque acting on the rigid body
//     Eigen::Vector3d net_torque = torque + calculate_torque(forces, point_of_application, center_of_mass);

//     // Calculate the change in orientation using the angular velocity
//     Eigen::Quaterniond delta_orientation;
//     delta_orientation.setIdentity();
//     delta_orientation.coeffs() = 0.5 * dt * Eigen::Vector3d(net_torque.x(), net_torque.y(), net_torque.z());

//     // Update the orientation of the rigid body
//     Eigen::Quaterniond new_orientation = current_orientation * delta_orientation;
//     new_orientation.normalize();

//     return new_orientation;
// }

