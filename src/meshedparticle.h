#pragma once
#include <iostream>
#include <cmath>
#include "meshedparticle.h"
#include <igl/point_mesh_squared_distance.h>
#include "igl/signed_distance.h"
#include <vector>
#include <iostream>
#include "meshops.h"


class ParticleAdhesion {
public:
void find_pairs(Eigen::MatrixXd V,Eigen::MatrixXi F, Eigen::MatrixXd V_particle, Eigen::VectorXd signed_distance,
                const double distance_threshold,Eigen::Vector3d center_of_mass, std::vector<std::pair<int, int>>& bonds);

void remove_long_bonds(std::vector<std::pair<int, int>>& bonds, Eigen::MatrixXd& V, Eigen::MatrixXd& V_particle,
                 double max_bond_length);                
private:
Eigen::VectorXi nearest;
};
