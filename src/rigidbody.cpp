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
   //std::cout << "Center of mass: " << center_of_mass.transpose()<< std::endl;
   // Calculate the moment of inertia
   Eigen::MatrixXd delta = points.rowwise() - center_of_mass.transpose(); // Vertex displacement vector in space frame (in the initial orientation)
  
   for (int i = 0; i < points.rows(); ++i) {
       r = delta.row(i); // Corrected vector subtraction
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
   // Diagonalize the moment of inertia tensor
   Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(moment_of_inertia);
   if (eigensolver.info() != Eigen::Success) abort();
  
   principal_axes = eigensolver.eigenvectors();
   principal_moments = eigensolver.eigenvalues().asDiagonal();
   Eigen::Vector3d cross= principal_axes.col(0).cross(principal_axes.col(1));
   if (cross.dot(principal_axes.col(2))<0)
   {
       principal_axes.col(2)=-principal_axes.col(2);
   }



   std::cout << "Principal axes:\n" << principal_axes << std::endl;
   //std::cout << "Principal moments of inertia (Diagonal):\n" << principal_moments << std::endl;
   // Optionally: Ensure the rotation matrix follows the right-hand rule
   // rotation_matrix.col(2) = rotation_matrix.col(0).cross(rotation_matrix.col(1));
// Return the calculated properties

    // Calculate the displacement vector for each point
   displace = delta * principal_axes; // Displacement vector in the principal axes frame/ body frame

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
        Eigen::Vector3d r = point_of_application.row(i) - center_of_mass.transpose();

        // Calculate the torque for this force and add it to the total
        torque += r.cross(forceVec);
        }
        //std::cout << "Torque: " << torque.transpose() << std::endl;
}
void RigidBody::angular_momentum (Eigen::Vector3d torque,double dt, Eigen::Vector3d& ang_mom) {
    ang_mom += torque * dt;  // Update angular momentum using ΔL = torque * Δt
}


void RigidBody::calculate_omega(Eigen::Vector3d angular_momentum, Eigen::Matrix3d rot_mat, Eigen::Matrix3d idiag, Eigen::Vector3d& angular_velocity) {
        // Transform angular momentum from space to body frame
        mom_body = rot_mat.transpose() * angular_momentum;  // P^T * Angular_Momentum_space

        // Compute angular velocity in the body frame
        // omega_body;//space_frame
        // for (int i = 0; i < 3; ++i) {
        //     // Check if diagonal inertia (inverse inertia component) is zero to prevent division by zero
        //     omega_body(i) = idiag(i, i) != 0 ? mom_body(i) / idiag(i, i) : 0.0;
        // }
        omega_body = mom_body.array() / idiag.diagonal().array(); // Safe division
        omega_body = omega_body.unaryExpr([](double v) { return std::isfinite(v) ? v : 0.0; }); // Handling division by zero

        // Transform angular velocity back to the space frame
        angular_velocity = rot_mat * omega_body;  // P * ω_body(body_frame) 
}


//simple euler update
//RichardsonIterativeMethod for quaternion update would be implemented sooon
//Update quaternion based on angular velocity and time step
void RigidBody::update_quaternion(Eigen::Quaterniond current_quaternion, Eigen::Vector3d angular_velocity, double dt, Eigen::Quaterniond& new_quaternion) {
        // Convert angular velocity to a pure quaternion (zero real part)
        Eigen::Quaterniond omega_q(0, angular_velocity.x(), angular_velocity.y(), angular_velocity.z());

        //omega_q.normalize();//space_frame

        Eigen::Quaterniond product =(omega_q * current_quaternion);

        //dq= (1/2)*q⊗Ω

        // // // Calculate the quaternion derivative using quaternion multiplication
        Eigen::Quaterniond quaternion_derivative(0.5 * product.w(),
                              0.5 * product.x(), 
                              0.5 * product.y(),
                              0.5 * product.z());
        //Eigen::Vector4d q_derivative = quaternion_derivative.coeffs();
        
        //std::cout<<"quaternion_derivative: ";
        //std::cout<<quaternion_derivative.w()<<" "<<quaternion_derivative.x()<<" "<<quaternion_derivative.y()<<" "<<quaternion_derivative.z()<<std::endl;
        // Update the quaternion using the derivative and time step
        new_quaternion.w() = current_quaternion.w() + quaternion_derivative.w()*dt;
        new_quaternion.x() = current_quaternion.x() + quaternion_derivative.x()*dt;
        new_quaternion.y() = current_quaternion.y() + quaternion_derivative.y()*dt;
        new_quaternion.z() = current_quaternion.z() + quaternion_derivative.z()*dt;
        // we need to multiply each component of the quaternion by 'dt' manually
        //Eigen::Vector4d q_derivative = quaternion_derivative.coeffs() * dt;

       // Update the quaternion by adding the scaled derivative to the current quaternion
        //new_quaternion.coeffs() = current_quaternion.coeffs() + 0.5*q_derivative*dt;

        // normalizing the quaternion
        new_quaternion.normalize();//space_frame
        //std::cout << "New Quaternion: " << new_quaternion.w() << " " << new_quaternion.x() << " " << new_quaternion.y() << " " << new_quaternion.z() << std::endl;
}


void RigidBody::exyz_to_q(Eigen::Matrix3d R ,Eigen::Quaterniond& quat) {
    // Convert the rotation matrix to a quaternion
    quat = Eigen::Quaterniond(R);
    quat.normalize();//space_frame
}
void RigidBody::q_to_exyz(Eigen::Quaterniond quat, Eigen::Matrix3d& R)
{   
    //quat.normalize();//space_frame
    // Convert the quaternion to a rotation matrix
    R = quat.toRotationMatrix();//space_frame
}

void RigidBody::rotate_vertices(Eigen::MatrixXd& vertices, Eigen::Vector3d center_of_mass, Eigen::MatrixXd displace, Eigen::Matrix3d rot_mat) {

    // Step 1: Apply the rotation matrix (convert vertice displacement vector from body frame to space frame)
    vertices =  displace * rot_mat.transpose();// V_rotated = Rotation_Matrix * V

    // Step 2: Translate vertices back to the original center of mass
    vertices.rowwise() += center_of_mass.transpose();
}
// Normalize the quaternion
void RigidBody::Quaternion::normalize() {
    double norm = std::sqrt(w * w + x * x + y * y + z * z);
    w /= norm;
    x /= norm;
    y /= norm;
    z /= norm;
}
// Quaternion multiplication using the definition
RigidBody::Quaternion RigidBody::Quaternion::operator*(const Quaternion& other) const {
    return Quaternion(
        w * other.w - x * other.x - y * other.y - z * other.z,
        w * other.x + x * other.w + y * other.z - z * other.y,
        w * other.y - x * other.z + y * other.w + z * other.x,
        w * other.z + x * other.y - y * other.x + z * other.w
    );
}
// Update the quaternion based on angular velocity
void RigidBody::update_Quaternion(const Quaternion& currentQuaternion, const Eigen::Vector3d& angularVelocity, double dt, Quaternion& newQuaternion) {
    // Convert angular velocity to a pure quaternion (zero real part)
    Quaternion omega_q(0, angularVelocity.x(), angularVelocity.y(), angularVelocity.z());

    // Quaternion derivative calculation
    Quaternion product = omega_q * currentQuaternion;
    Quaternion quaternion_derivative(
        0.5 * product.w, 0.5 * product.x, 0.5 * product.y, 0.5 * product.z);

    // Update the quaternion using the derivative and time step
    newQuaternion = Quaternion(
        currentQuaternion.w + quaternion_derivative.w * dt,
        currentQuaternion.x + quaternion_derivative.x * dt,
        currentQuaternion.y + quaternion_derivative.y * dt,
        currentQuaternion.z + quaternion_derivative.z * dt
    );

    // Normalize the quaternion to maintain unit length
    newQuaternion.normalize();
}
//// V= vcm + omega x r  (r= V - com)
