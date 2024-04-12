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
/*
void RigidBody::update_position_orientation(const Eigen::MatrixXd& forces, const Eigen::MatrixXd& point_of_application, const Eigen::Vector3d& torque, const double dt, Eigen::Vector3d& new_position) {
    // Calculate the net force acting on the center of mass
    Eigen::Vector3d net_force = forces.colwise().sum();

    // Calculate the linear acceleration of the center of mass
    Eigen::Vector3d linear_acceleration = net_force / mass;

    // Calculate the change in position using the linear acceleration
    Eigen::Vector3d delta_position = linear_acceleration * dt;

    // Calculate the change in orientation using the angular velocity
    Eigen::Quaterniond delta_orientation;
    delta_orientation.setIdentity();
    delta_orientation.coeffs() = 0.5 * dt * Eigen::Vector3d(angular_velocity.x(), angular_velocity.y(), angular_velocity.z());

    // Update the position of the center of mass
    new_position += delta_position;

    // Update the orientation of the center of mass
    Eigen::Quaterniond current_orientation;
    current_orientation.setIdentity();
    current_orientation.coeffs() = orientation;

    Eigen::Quaterniond new_orientation = current_orientation * delta_orientation;
    new_orientation.normalize();

    // Update the angular velocity
    Eigen::Vector3d new_angular_velocity = angular_velocity + torque * dt / moment_of_inertia;

    // Set the new position, orientation, and angular velocity
    position = new_position;
    orientation = new_orientation.coeffs();
    angular_velocity = new_angular_velocity;
}
*/

