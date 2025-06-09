#pragma once
#include <string>

// -----------------------------------------------------------------------------
// Parameter struct: holds all settings read from run_parameters.txt
// -----------------------------------------------------------------------------
struct Parameter {
    // ───── Simulation control ─────────────────────────────────────────────────
    int    iterations;
    double dt;

    double tolerance;
    int    tolerance_flag;   // 1 = ON, 0 = OFF
    double tolfrequency;     // check energy‐change every this many time units

    int    logfrequency;
    int    dumpfrequency;
    int    resfrequency;
    int    mesh_reg_frequency;
    int    bondfrequency;

    int    vertex_smoothing_flag;       // 0 = OFF, 1 = ON
    int    delaunay_triangulation_flag; // 0 = OFF, 1 = ON
    int    random_force_flag;           // 0 = OFF, 1 = ON

    double gamma;  // drag coefficient (same for both meshes)
    double mass;   // mass coefficient (same for both meshes)
    double kbT;    // thermal energy (for random‐force flag)

    // ───── Primary Mesh (mesh1: “membrane”) parameters ────────────────────────
    double Kb1;             // bending modulus for mesh 1
    double Ka1;             // stretching modulus for mesh 1
    double Kv1;             // osmotic modulus for mesh 1
    double reduced_volume1; // target reduced volume for mesh 1

    std::string meshFile1;  // OFF file path for mesh 1
    std::string outFile1;   // where to write final OFF for mesh 1
    std::string resFile1;   // restart file path (mesh 1)

    // ───── Secondary Mesh (mesh2: “vesicle”) parameters ──────────────────────
    double Kb2;             // bending modulus for mesh 2
    double Ka2;             // stretching modulus for mesh 2
    double Kv2;             // osmotic modulus for mesh 2
    double reduced_volume2; // target reduced volume for mesh 2

    std::string meshFile2;  // OFF file path for mesh 2
    std::string outFile2;   // where to write final OFF for mesh 2
    std::string resFile2;   // restart file path (mesh 2)

    // ───── Vesicle‐adhesion flags ────────────────────────────────────────────
    int    vesicle_flag;         // 1 = use mesh 2 as a vesicle, 0 = do not
    int    vesicle_position;     // –1 = inside, +1 = outside (legacy usage)
    int    vesicle_coord_flag;   // 1 = initial coords provided, 0 otherwise
    double vesicle_radius;       // legacy “vesicle radius” (unused if no vesicle)

    double adhesion_strength;    // “u” between mesh 1↔mesh 2
    double mass_vesicle_ratio;   // mass(mesh 2) / mass(mesh 1)
    double potential_range;      // ρ for adhesion
    double r_equilibrium;        // equilibrium separation
    int    angle_condition_flag; // 1 = enforce angle criterion, 0 = ignore

    int    forced_wrapping_flag;    // 1 = ON, 0 = OFF
    double wrapping_fraction;       // fraction of vesicle area to forcibly wrap
    double wrapping_bias_strength;  // Kw

    // ───── Member function to read from run_parameters.txt ────────────────
    void readParameter();
};
// -----------------------------------------------------------------------------