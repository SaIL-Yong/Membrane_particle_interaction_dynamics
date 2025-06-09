#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Dense>
#include <string>

// Replace the macro‐style PI definition with a constexpr to avoid conflicts
static constexpr double PI = 3.14159265358979323846;

struct SimulationData
{
    // ─── Mesh 1 (membrane) ───────────────────────────────────────────────────
    Eigen::MatrixXd V1;  // (#V1 × 3) vertex positions
    Eigen::MatrixXi F1;  // (#F1 × 3) faces
    int numV1, numF1;

    // Per‐vertex forces on Mesh 1:
    Eigen::MatrixXd Force_Bending1;   // (#V1 × 3)
    Eigen::MatrixXd Force_Area1;      // (#V1 × 3)
    Eigen::MatrixXd Force_Volume1;    // (#V1 × 3)
    Eigen::MatrixXd Force_Random1;    // (#V1 × 3)
    Eigen::MatrixXd Force_Drag1;      // (#V1 × 3)
    Eigen::MatrixXd Force_Adhesion;   // (#V1 × 3)
    Eigen::MatrixXd Force_Total1;     // (#V1 × 3)

    // Velocities (mesh 1)
    Eigen::MatrixXd velocity1;        // (#V1 × 3)

    // Energies (mesh 1)
    double EnergyBending1;
    double EnergyArea1;
    double EnergyVolume1;
    double EnergyPotential1;  // = bending + area + volume + adhesion + bias
    double EnergyKinetic1;    // = ½·mass1·∑‖v_i‖²

    // Mesh‐regularization helpers (mesh 1)
    Eigen::MatrixXd l1;              // edge lengths for intrinsic Delaunay
    double area_target1;             // target area (mesh 1)
    double volume_target1;           // target volume (mesh 1)
    double rVol1;                    // current reduced‐volume = 6√π·V1·A1⁻³/²
    double rVol_t1;                  // target reduced‐volume (mesh 1)
    double Kb1, Ka1, Kv1;            // bending/area/osmotic moduli (mesh 1)
    double Rv1;                      // initial “radius” for mesh 1
    double gamma1, mass1, kbT1;      // drag, mass, thermal (mesh 1)
    int    v_smooth_flag1;           // 1 = enable vertex smoothing
    int    delaunay_tri_flag1;       // 1 = enable intrinsic Delaunay
    int    random_force_flag1;       // 1 = add random‐force
    int    mesh_reg_frequency1;      // how often to run mesh‐reg (in steps)

    // ─── Mesh 2 (vesicle) ────────────────────────────────────────────────────
    Eigen::MatrixXd V2;  // (#V2 × 3) vertex positions
    Eigen::MatrixXi F2;  // (#F2 × 3) faces
    int numV2, numF2;

    // Per‐vertex forces on Mesh 2:
    Eigen::MatrixXd Force_Bending2;   // (#V2 × 3)
    Eigen::MatrixXd Force_Area2;      // (#V2 × 3)
    Eigen::MatrixXd Force_Volume2;    // (#V2 × 3)
    Eigen::MatrixXd Force_Random2;    // (#V2 × 3)
    Eigen::MatrixXd Force_Drag2;      // (#V2 × 3)
    Eigen::MatrixXd Force_Repulsion;  // (#V2 × 3)
    Eigen::MatrixXd Force_Total2;     // (#V2 × 3)

    // Velocities (mesh 2)
    Eigen::MatrixXd velocity2;        // (#V2 × 3)

    // Energies (mesh 2)
    double EnergyBending2;
    double EnergyArea2;
    double EnergyVolume2;
    double EnergyPotential2;  // = bending + area + volume
    double EnergyKinetic2;    // = ½·mass2·∑‖v_i‖²

    // Mesh‐regularization helpers (mesh 2)
    Eigen::MatrixXd l2;              // edge lengths for intrinsic Delaunay
    double area_target2;             // target area (mesh 2)
    double volume_target2;           // target volume (mesh 2)
    double rVol2;                    // current reduced‐volume (mesh 2)
    double rVol_t2;                  // target reduced‐volume (mesh 2)
    double Kb2, Ka2, Kv2;            // bending/area/osmotic moduli (mesh 2)
    double Rv2;                      // initial “radius” for mesh 2
    double gamma2, mass2, kbT2;      // drag, mass, thermal (mesh 2)
    int    v_smooth_flag2;           // 1 = enable vertex smoothing (mesh 2)
    int    delaunay_tri_flag2;       // 1 = enable intrinsic Delaunay (mesh 2)
    int    random_force_flag2;       // 1 = add random‐force (mesh 2)
    int    mesh_reg_frequency2;      // how often to run mesh‐reg (mesh 2)

    // ─── Shared adhesion fields ─────────────────────────────────────────────
    Eigen::MatrixXd Force_Total;          // Combined force field (mesh 1 + mesh 2)
    double force_residual;                // ‖Force_Total‖

    Eigen::MatrixXd signed_distance;      // (#V1 × 1) signed‐distance to mesh 2
    Eigen::VectorXi facet_index;          // (#V1 × 1) nearest‐face index on mesh 2
    Eigen::MatrixXd closest_points;       // (#V1 × 3) closest points on mesh 2
    Eigen::MatrixXd normals_closest_points; // (#V1 × 3) normals at closest points

    double EnergyAdhesion;         // ∑ᵢ u(dᵢ)·(Voronoi areaᵢ)
    double EnergyBias;             // = 0.5·Kw·(EAd – Ew_t)²
    double Ew_t;                   // target wrapping energy = –u·Area_w_t
    double Kw;                     // wrapping‐bias modulus
    double rc;                     // cutoff radius = 10·ρ
    double potential_range;        // ρ (interaction range)
    double adhesion_strength;      // u (adhesion strength)
    double r_equilibrium;          // equilibrium distance parameter (for bias)
    double u;  // “raw” adhesion strength from parameter file
    double U;  // scaled adhesion strength = (Kb * u)/(Rp²)

    // Vesicle‐flag parameters
    int    vesicle_flag;           // 1 = use mesh 2 as vesicle, 0 = skip mesh 2
    int    vesicle_position;       // –1 = inside, +1 = outside
    int    vesicle_coord_flag;     // 1 = initial coords provided
    double vesicle_radius;         // initial guess for vesicle radius
    double mass_vesicle_ratio;     // = mass2 / mass1
    int    angle_flag;             // 1 = enforce angle‐condition, 0 = ignore
    int    forced_wrapping_flag;   // 1 = ON, 0 = OFF
    double wrapping_fraction;      // fraction of vesicle area to force‐wrap

    // ─── Logging & convergence ●────────────────────────────────────────────
    double EnergyPotential;        // mesh 1 + mesh 2 potentials
    double EnergyKinetic;          // mesh 1 + mesh 2 kinetics
    double EnergyTotal;            // total = potential + kinetic
    double EnergyTotalOld_log;     // for computing “log‐interval” change
    double EnergyChangeRate_log;   // = (EnergyTotal – EnergyTotalOld_log) / (dt · logfreq)
    double EnergyChangeRate_avg;   // running average over tolmean_steps
    Eigen::VectorXd etol;          // stores recent change‐rates (length ≈ iterations/logfreq)
    int tolsteps;                  // how many steps between tolerance checks
    int tolmean_steps;             // how many log-intervals to average
    double tolerance;              // threshold for convergence
    int    tolerance_flag;         // 1 = check convergence, 0 = ignore
    double tolfrequency;           // how often (in simulation-time) to check

    // ─── Time-stepping control ●──────────────────────────────────────────────
    int    iterations;             // total timesteps
    int    logfrequency;           // how often (in steps) to write to logfile
    int    dumpfrequency;          // how often (in steps) to write OFF dumps
    int    resfrequency;           // restart-file frequency
    int    mesh_reg_frequency;     // common default; per-mesh we have 1 and 2
    int    bondfrequency;          // (unused but read from file)

    double dt;                     // timestep
    double time;                   // current simulation time = iteration · dt
    double dtf;                    // = dt/2 (half-step)

    // ─── File-name fields (used by integrate.cpp) ───────────────────────────
    std::string outFile1;          // "./out_membrane.off"
    std::string outFile2;          // "./out_second.off"
    std::string resFile1;          // restart-file for mesh 1
    std::string resFile2;          // restart-file for mesh 2

    // ─── “Constructor”-style helpers to zero-initialize buffers ──────────────
    void initialize(int vertexCount)
    {
        numV1 = vertexCount;
        Force_Bending1.setZero(numV1, 3);
        Force_Area1.setZero(numV1, 3);
        Force_Volume1.setZero(numV1, 3);
        Force_Random1.setZero(numV1, 3);
        Force_Drag1.setZero(numV1, 3);
        Force_Adhesion.setZero(numV1, 3);
        Force_Total1.setZero(numV1, 3);

        velocity1.setZero(numV1, 3);

        signed_distance.setZero(numV1, 1);
        facet_index.setZero(numV1);
        closest_points.setZero(numV1, 3);
        normals_closest_points.setZero(numV1, 3);
    }

    void initialize_mesh2(int vertexCount)
    {
        numV2 = vertexCount;
        Force_Bending2.setZero(numV2, 3);
        Force_Area2.setZero(numV2, 3);
        Force_Volume2.setZero(numV2, 3);
        Force_Random2.setZero(numV2, 3);
        Force_Drag2.setZero(numV2, 3);
        Force_Repulsion.setZero(numV2, 3);
        Force_Total2.setZero(numV2, 3);

        velocity2.setZero(numV2, 3);
    }
};

#endif
