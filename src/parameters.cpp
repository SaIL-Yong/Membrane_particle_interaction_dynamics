// parameters.cpp

#include "parameters.h"
#include <iostream>
#include <fstream>
#include <string>

void Parameter::readParameter()
{
    std::ifstream runfile("run_parameters.txt");
    if (!runfile.is_open())
    {
        std::cerr << "ERROR: cannot open run_parameters.txt\n";
        return;
    }

    std::string key;

    // 1) iterations
    runfile >> key >> iterations;
    // 2) dt
    runfile >> key >> dt;

    // --- Mesh 1 (membrane) energetics ---
    runfile >> key >> Kb1;
    runfile >> key >> Ka1;
    runfile >> key >> Kv1;
    runfile >> key >> reduced_volume1;

    runfile >> key >> tolerance;
    runfile >> key >> tolerance_flag;
    runfile >> key >> tolfrequency;

    runfile >> key >> gamma;
    runfile >> key >> mass;
    runfile >> key >> kbT;

    // --- Mesh 2 (vesicle) energetics ---
    runfile >> key >> Kb2;
    runfile >> key >> Ka2;
    runfile >> key >> Kv2;
    runfile >> key >> reduced_volume2;

    // --- Vesicle-adhesion flags ---
    runfile >> key >> vesicle_flag;
    // runfile >> key >> vesicle_position;
    // runfile >> key >> vesicle_coord_flag;
    // runfile >> key >> vesicle_radius;

    runfile >> key >> adhesion_strength;
    runfile >> key >> mass_vesicle_ratio;
    runfile >> key >> potential_range;
    runfile >> key >> r_equilibrium;
    runfile >> key >> angle_condition_flag;
    runfile >> key >> forced_wrapping_flag;
    runfile >> key >> wrapping_fraction;
    runfile >> key >> wrapping_bias_strength;

    // --- File paths for mesh1 & mesh2 ---
    runfile >> key >> meshFile1;
    runfile >> key >> meshFile2;
    runfile >> key >> outFile1;
    runfile >> key >> outFile2;
    runfile >> key >> resFile1;
    runfile >> key >> resFile2;

    // logfrequency
    runfile >> key >> logfrequency;
    // dumpfrequency
    runfile >> key >> dumpfrequency;
    // resfrequency
    runfile >> key >> resfrequency;
    // mesh_reg_frequency
    runfile >> key >> mesh_reg_frequency;
    // bondfrequency
    runfile >> key >> bondfrequency;

    runfile >> key >> vertex_smoothing_flag;
    runfile >> key >> delaunay_triangulation_flag;
    runfile >> key >> random_force_flag;

    runfile.close();
}
