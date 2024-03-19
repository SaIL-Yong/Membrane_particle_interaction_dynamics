#pragma once
#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#endif
#include <iostream>
#include <fstream>
#include <string>

void readParameter();
struct Parameter {
    int iterations, tolerance_flag, forced_wrapping_flag;
    int logfrequency, dumpfrequency, resfrequency,bondfrequency;
    double tolfrequency;
    int particle_flag, particle_position;
    double dt, Kb, Kv, Ka, gamma, tolerance, wrapping_fraction, wrapping_bias_strength;
    double reduced_volume, particle_radius, adhesion_strength, potential_range,r_equilibrium;
    std::string meshFile, outFile, resFile,particleFile;
    int mesh_reg_frequency, vertex_smoothing_flag, delaunay_triangulation_flag,lennard_jones_flag;
    int angle_condition_flag;
    double X0, Y0, Z0;
    int particle_coord_flag;
    double epsilon,sigma;
};
