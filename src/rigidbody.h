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
    Eigen::Vector3d r;
    //ang_mom_ to omega calculation
    Eigen::Vector3d mom_body;
    
public:
    Eigen::Vector3d omega_body;
    void calculate_properties(Eigen::MatrixXd points, double mass,Eigen::Matrix3d& principal_axes, Eigen::Matrix3d& principal_moments,Eigen::MatrixXd& displace);
     // Accessor methods for moment of inertia
    Eigen::Vector3d getCenterOfMass() const { return center_of_mass; }

    void calculate_center_of_mass(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::Vector3d& center_of_mass);
    //Torque calculation
    void calculate_torque(Eigen::MatrixXd force, Eigen::MatrixXd point_of_application,Eigen::Vector3d center_of_mass,Eigen::Vector3d& torque);
    void calculate_drag_torque(Eigen::Vector3d omega_body, Eigen::Matrix3d rot_mat,Eigen::Matrix3d idiag, double gamma, Eigen::Vector3d& drag_torque);

    void angular_momentum (Eigen::Vector3d torque,double dt, Eigen::Vector3d& angular_momentum);
    void calculate_omega(Eigen::Vector3d angular_momentum, Eigen::Matrix3d rot_mat, Eigen::Matrix3d idiag, Eigen::Vector3d& angular_velocity, Eigen::Vector3d& omega_body);
    void update_quaternion(Eigen::Quaterniond current_quaternion, Eigen::Vector3d angular_velocity, double dt, Eigen::Quaterniond& new_quaternion);
    //void update_vertex_position(Eigen::MatrixXd& V, Eigen::MatrixXd forces,  Eigen::Quaterniond quaternion, double dt, Eigen::Vector3d& particle_velocity_com);
    void update_vertex_velocities_positions(Eigen::MatrixXd& V,Eigen::Matrix3d rot_mat,Eigen::Vector3d vcm,Eigen::Vector3d omega,Eigen::MatrixXd displace,double dt, Eigen::MatrixXd& node_velocities);
    //void diagonalize_inertia_tensor(Eigen::Matrix3d inertia_tensor, Eigen::Matrix3d& principal_axes,Eigen::Matrix3d& principal_moments);
    void exyz_to_q(Eigen::Matrix3d R ,Eigen::Quaterniond& quat);
    void q_to_exyz(Eigen::Quaterniond quat, Eigen::Matrix3d& R);
    //void rotate_vertices(Eigen::MatrixXd& vertices,Eigen::Vector3d center_of_mass, Eigen::Matrix3d rot_mat);
    void rotate_vertices(Eigen::MatrixXd& vertices, Eigen::Vector3d center_of_mass, Eigen::MatrixXd displace, Eigen::Matrix3d rot_mat);
    //void printTorque(Eigen::MatrixXd force,  Eigen::MatrixXd point_of_application, Eigen::Vector3d center_of_mass);
    // Declare a method to calculate the orientation of the rigid body
    //Eigen::Vector3d  ang_mom,omega; 



    //quaternion opretations
    struct Quaternion {
    double w, x, y, z;
        
    // Default and parameterized constructor
    Quaternion(double w = 1.0, double x = 0.0, double y = 0.0, double z = 0.0)
            : w(w), x(x), y(y), z(z) {}

    // Conversion from Eigen::Quaterniond
        // Conversion constructor from Eigen::Quaterniond
    Quaternion(const Eigen::Quaterniond& eq)
            : w(eq.w()), x(eq.x()), y(eq.y()), z(eq.z()) {}

    // Normalize the quaternion
    void normalize(); 
    // Multiply this quaternion by another and return the result
    Quaternion operator*(const Quaternion& other) const;

    };

    void update_Quaternion(const Quaternion& currentQuaternion, const Eigen::Vector3d& angularVelocity, double dt, Quaternion& newQuaternion);
    

    

};
#endif // RIGIDBODY_H