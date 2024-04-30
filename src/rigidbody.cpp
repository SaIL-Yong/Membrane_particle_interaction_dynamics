#include "rigidbody.h"
#include "meshops.h"
#include <Eigen/Core>
#include <igl/doublearea.h>
#include <Eigen/Geometry> // For cross product and norm


////This function calculate momemt of inertia and center of mass in space frame
void RigidBody::calculate_properties(Eigen::MatrixXd points, double mass,Eigen::Matrix3d& principal_axes, Eigen::Matrix3d& principal_moments,Eigen::MatrixXd& displace) { 
   center_of_mass = Eigen::Vector3d::Zero();
   //int numPoints = points.rows();
 // Assuming 'points' is of type Eigen::MatrixXd (or Eigen::MatrixXf) and each row is a point
   center_of_mass = points.colwise().mean();
   std::cout << "Center of mass: " << center_of_mass.transpose()<< std::endl;
   displace = points.rowwise() - center_of_mass.transpose();
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


   /* moment of inertia inverse
   if (moment_of_inertia.determinant() != 0) {
           inverse_moment_of_inertia = moment_of_inertia.inverse();
       } else {
           std::cerr << "Warning: Moment of inertia matrix is not invertible." << std::endl;
   }
   */
  moment_of_inertia.normalize();// normalizing inertia tensor following LAMMPS
   // Diagonalize the moment of inertia tensor
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(moment_of_inertia);
   if (eigensolver.info() != Eigen::Success) abort();
  
   principal_axes = eigensolver.eigenvectors();
   principal_moments = eigensolver.eigenvalues().asDiagonal();


   std::cout << "Principal axes:\n" << principal_axes << std::endl;
   std::cout << "Principal moments of inertia (Diagonal):\n" << principal_moments << std::endl;
   // Optionally: Ensure the rotation matrix follows the right-hand rule
   // rotation_matrix.col(2) = rotation_matrix.col(0).cross(rotation_matrix.col(1));
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
        //Eigen::Vector3d pointVec = point_of_application.row(i).transpose();

        // Calculate the vector from the center of mass to the point of application
        Eigen::Vector3d r = point_of_application.row(i).transpose() - center_of_mass;

        // Calculate the torque for this force and add it to the total
        torque += r.cross(forceVec);
        }
}
void RigidBody::angular_momentum (Eigen::Vector3d torque,double dt, Eigen::Vector3d& ang_mom) {
    ang_mom += torque * dt;  // Update angular momentum using ΔL = torque * Δt
}


void RigidBody::calculate_omega(Eigen::Vector3d angular_momentum, Eigen::Matrix3d rot_mat, Eigen::Matrix3d idiag, Eigen::Vector3d& angular_velocity) {
        // Transform angular momentum from space to body frame
        mom_body = rot_mat.transpose() * angular_momentum;  // P^T * Angular_Momentum_space

        // Compute angular velocity in the body frame
        omega_body;//space_frame
        for (int i = 0; i < 3; ++i) {
            // Check if diagonal inertia (inverse inertia component) is zero to prevent division by zero
            omega_body(i) = idiag(i, i) == 0.0 ? 0.0 : mom_body(i) / idiag(i, i);
        }

        // Transform angular velocity back to the space frame
        angular_velocity = rot_mat * omega_body;  // P * ω_body(body_frame) 
}


//simple euler update
//RichardsonIterativeMethod for quaternion update would be implemented sooon
//Update quaternion based on angular velocity and time step
void RigidBody::update_quaternion(Eigen::Quaterniond current_quaternion, Eigen::Vector3d angular_velocity, double dt, Eigen::Quaterniond& new_quaternion) {
        // Convert angular velocity to a pure quaternion (zero real part)
        Eigen::Quaterniond omega_q(0, angular_velocity.x(), angular_velocity.y(), angular_velocity.z());

        // // Calculate the quaternion derivative using quaternion multiplication
         Eigen::Quaterniond product = omega_q*current_quaternion;
        // Update the quaternion using the derivative and time step
        new_quaternion.w() = current_quaternion.w() + 0.5 *dt* product.w(),
        new_quaternion.x() = current_quaternion.x() + 0.5 *dt* product.x(),
        new_quaternion.y() = current_quaternion.y() + 0.5 *dt* product.y(),
        new_quaternion.z() = current_quaternion.z() + 0.5 *dt* product.z();
        // normalizing the quaternion
        new_quaternion.normalize();//space_frame
}

// void RigidBody::diagonalize_inertia_tensor(Eigen::Matrix3d inertia_tensor, Eigen::Matrix3d& principal_axes, Eigen::Matrix3d& principal_moments) {
//     Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(inertia_tensor);
//     if (eigensolver.info() != Eigen::Success) abort();
    
//     principal_axes = eigensolver.eigenvectors();
//     principal_moments = eigensolver.eigenvalues().asDiagonal();
    
//     std::cout << "Principal axes:\n" << principal_axes << std::endl;
//     std::cout << "Principal moments of inertia (Diagonal):\n" << principal_moments << std::endl;
//     // Optionally: Ensure the rotation matrix follows the right-hand rule
//     // rotation_matrix.col(2) = rotation_matrix.col(0).cross(rotation_matrix.col(1)); 
// }
void RigidBody::exyz_to_q(Eigen::Matrix3d R ,Eigen::Quaterniond& quat) {
    // Convert the rotation matrix to a quaternion
    quat = Eigen::Quaterniond(R);
    quat.normalize();//space_frame
}
void RigidBody::q_to_exyz(Eigen::Quaterniond quat, Eigen::Matrix3d& R)
{
    // Convert the quaternion to a rotation matrix
    R = quat.toRotationMatrix();//space_frame
}

//// V= vcm + omega x r  (r= V - com)
void RigidBody::update_vertex_velocities_positions(Eigen::MatrixXd& V,Eigen::Matrix3d rot_mat,Eigen::Vector3d vcm,Eigen::Vector3d omega,Eigen::MatrixXd displace,double dt, Eigen::MatrixXd& node_velocities) {
        // velocities is a matrix with the same dimensions as V
        // Calculate vertex velocities
        // Calculate vertex velocities
        // for (int i = 0; i < V.rows(); i++) {
        //     Eigen::Vector3d r = rot_mat * displace.row(i).transpose(); // Apply rotation to displacement
        //     Eigen::Vector3d rotational_velocity = omega.cross(r);      // ω_space×r_space
        //     node_velocities.row(i) = vcm+ rotational_velocity; // v=vcm_space +ω_space×r_space
        // }
        // Update vertex positions
        //Rotate the vertices

        V = V * rot_mat.transpose();

}

void RigidBody::rotate_vertices(Eigen::MatrixXd& vertices,Eigen::Vector3d center_of_mass ,Eigen::Matrix3d rot_mat) {
    // Step 1: Translate vertices to center the COM at the origin
    vertices.rowwise() -= center_of_mass.transpose();

    // Step 2: Apply the rotation matrix
    vertices = (vertices * rot_mat.transpose());// V_rotated = Rotation_Matrix * V

    // Step 3: Translate vertices back to the original center of mass
    vertices.rowwise() += center_of_mass.transpose();
}
/*
// Function to update vertex positions using matrix operations
void RigidBody::update_vertex_position(Eigen::MatrixXd& V, Eigen::MatrixXd forces,  Eigen::Quaterniond quaternion, double dt, Eigen::Vector3d& particle_velocity_com) {
        // Calculate the net force (assuming each column in forces is a force vector for the corresponding vertex)
        //Eigen::Vector3d net_force = forces.rowwise().sum();
        //std::cout << "Net Force: " << net_force.transpose() << std::endl;

        // Calculate the acceleration of the center of mass based on the net force
        particle_acceleration_com =  forces.colwise().sum() / V.rows();

        // Update the velocity of the center of mass based on the acceleration
        particle_velocity_com += particle_acceleration * dt;
        //std::cout << "Particle Velocity: " << particle_velocity.transpose() << std::endl;

        // Update all vertex positions by translating with the velocity
        V.rowwise() += (particle_velocity * dt).transpose();

        // Calculate the rotation matrix from the quaternion
        rotation_matrix = quaternion.toRotationMatrix();
        //Affine3d T = quaternion * Translation3d(V);

               // Apply the rotation to all vertex positions
        V = (rotation_matrix * V.transpose()).transpose();
        //V= T * V;
}
*/