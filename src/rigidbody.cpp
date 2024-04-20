#include "rigidbody.h"
#include "meshops.h"
#include <Eigen/Core>
#include <igl/doublearea.h>
#include <Eigen/Geometry> // For cross product and norm

void RigidBody::calculate_properties(Eigen::MatrixXd points, double mass/*,Eigen::Matrix3d& moment_of_inertia, Eigen::Matrix3d& inverse_moment_of_inertia*/) {  
    Eigen::Vector3d center_of_mass = Eigen::Vector3d::Zero();
    //int numPoints = points.rows(); 
// Assuming 'points' is of type Eigen::MatrixXd (or Eigen::MatrixXf) and each row is a point
    center_of_mass = points.colwise().mean();
    std::cout << "Center of mass: " << center_of_mass.transpose()<< std::endl;

    // Calculate the moment of inertia
    
    for (int i = 0; i < points.rows(); ++i) {
        Eigen::Vector3d r = points.row(i).transpose() - center_of_mass; // Corrected vector subtraction
        //std::cout << "r: " << r.y() << std::endl;
        moment_of_inertia(0, 0) += mass * (r.y() * r.y() + r.z() * r.z());
        //std::cout<<"moment_of_inertia(0, 0): "<<moment_of_inertia(0, 0)<<std::endl;
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
void RigidBody::angular_momentum (Eigen::Vector3d torque,double dt, Eigen::Vector3d& ang_mom) {
    ang_mom += torque * dt;  // Update angular momentum using ΔL = torque * Δt
}
void RigidBody::calculate_omega(Eigen::Vector3d angular_momentum, Eigen::Matrix3d inverse_moment_of_inertia,Eigen::Quaterniond quat,Eigen::Vector3d& angular_velocity) {
    Eigen::Matrix3d Iinv = quat.toRotationMatrix()*inverse_moment_of_inertia*quat.toRotationMatrix().transpose();
    angular_velocity = Iinv* angular_momentum; // ω = I^-1 * L
}


//Update quaternion based on angular velocity and time step
void RigidBody::update_quaternion(Eigen::Quaterniond current_quaternion, Eigen::Vector3d angular_velocity, double dt, Eigen::Quaterniond& new_quaternion) {
        // Convert angular velocity to a pure quaternion (zero real part)
        Eigen::Quaterniond omega_q(0, angular_velocity.x(), angular_velocity.y(), angular_velocity.z());

        // // Calculate the quaternion derivative using quaternion multiplication
         Eigen::Quaterniond product = current_quaternion * omega_q;

        // // Correctly applying scalar multiplication to the quaternion components
        // Eigen::Quaterniond quaternion_derivative(product.w() * 0.5, product.x() * 0.5, product.y() * 0.5, product.z() * 0.5);
            // Quaternion derivative should be half of current_quaternion * omega_q
            // Multiply quaternion coefficients by 0.5
        Eigen::Vector4d coeffs = product.coeffs() * 0.5;

        // Reconstruct the quaternion from scaled coefficients
        Eigen::Quaterniond quaternion_derivative(coeffs);

        // Update the quaternion using the derivative and the time step
           // Compute the magnitude of the angular velocity vector
        double omega_magnitude = angular_velocity.norm();
    
        // Avoid division by zero in case of very small angular velocities
        if(omega_magnitude > std::numeric_limits<double>::epsilon()) {
        // Normalized angular velocity vector
        Eigen::Vector3d omega_normalized = angular_velocity.normalized();
        // Compute the scaled angle for the quaternion exponential
        double theta = omega_magnitude * dt * 0.5; 
        // Compute the quaternion exponential
        Eigen::Quaterniond delta_q(std::cos(theta), 
                                   omega_normalized.x() * std::sin(theta),
                                   omega_normalized.y() * std::sin(theta),
                                   omega_normalized.z() * std::sin(theta));
        
        // Update the quaternion by multiplying with the exponential quaternion
        new_quaternion = (delta_q * current_quaternion).normalized();
        } else {
        // For very small angular velocities, you might simply normalize the quaternion
        // as the change is negligible
        new_quaternion.normalize();
    }
    //  Output the rotation matrix for debugging/verification
    //std::cout << "Updated Rotation Matrix:\n" << new_quaternion.toRotationMatrix() << std::endl;
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

