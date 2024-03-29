#pragma once
#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <cmath>
#include "meshops.h"



#include <Eigen/Dense>
#include <vector>

class RigidBody {
public:
    struct MeshProperties {
        Eigen::Vector3d center_of_mass;
        Eigen::Matrix3d inertia_tensor;
    };

    void calculateProperties(
        Eigen::MatrixXD points, double mass,
        Eigen::Vector3d& cente_of_mass,
        Eigen::Matrix3d& moment_of_inertia);
};

#endif // RIGIDBODY_H