// simulation.h
#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Dense>
#include <Eigen/Geometry>  // For Eigen::Quaterniond

struct SimulationData {
    // Mesh data
    Eigen::MatrixXd V1, V2;
    Eigen::MatrixXi F1, F2;
    int numV, numF, numFp, numVp;
    Eigen::MatrixXd l;

    //File input output
    std::string outFile,meshFile,particleFile;

    // Dynamics and forces
    Eigen::MatrixXd Force_Area, Force_Volume, Force_Bending, Force_Adhesion,ForcesOnParticle;
    Eigen::MatrixXd Force_Random, Force_Drag, Force_Repulsion, Force_Total;
    Eigen::MatrixXd acceleration, acceleration_half_step, velocity, velocity_half_step;
    Eigen::MatrixXd particle_velocities, displace;

    // Collision detection data
    //Eigen::Vector3d center_of_mass;
    Eigen::MatrixXd signed_distance;
    Eigen::MatrixXi facet_index;
    Eigen::MatrixXd closest_points, normals_closest_points;

    // Energies
    double EnergyVolume, EnergyArea, EnergyBending, EnergyAdhesion, EnergyBias;
    double EnergyPotential, EnergyKinetic, EnergyTotal, EnergyTotalold_log, EnergyChangeRate_log, EnergyChangeRate_avg;
    double EnergyParticleKinetic, EnergyParticleKineticTranslation, EnergyParticleKineticRotation;

    // Physical Membrane properties
    double gamma, mass, kbT, Kb, Kv, Ka, Rv, area_target, volume_target, rVol;
    double Rp, u, U, rho, rc, r_equilibrium;
    double mass_particle, total_mass_particle;
    double force_residual, distance_threshold;
    double rVol_t; // true reduced volume
    // Rigid body dynamics properties
    Eigen::Matrix3d rotation_matrix, idiag;
    Eigen::Vector3d torque, ang_momentum, ang_velocity, particle_velocity_com, particle_acceleration_com;
    Eigen::Quaterniond current_quaternion, new_quaternion;
    Eigen::Vector3d COM,center_of_mass;
    // Simulation control parameters
    int iterations, logfrequency, dumpfrequency, bondfrequency, resfrequency;
    double dt, time, dtf;
    double tolerance;
    int tolerance_flag;
    double tolfrequency;
    int tolsteps, tolmean_steps;
    Eigen::VectorXd etol;
    // Particle adhesion and wrapping properties
    int particle_position;
    double Ew_t, Kw;  // Energy and strength for wrapping
    double forced_wrapping_fraction;

    // Mesh regularization properties
    int v_smooth_flag, delaunay_tri_flag, mesh_reg_frequency;
    //flags
    int angle_flag, random_force_flag, particle_flag,forced_wrapping_flag;

    // Constructor
    // SimulationData() : numV(1000), numF(0), numFp(0), numVp(0), iterations(1000), logfrequency(10), 
    //                    dumpfrequency(100), bondfrequency(50), resfrequency(500), dt(0.01), time(0.0), 
    //                    dtf(0.005), tolerance(0.01), tolerance_flag(1), tolfrequency(100.0), tolsteps(10), 
    //                    tolmean_steps(1), gamma(1.0), mass(1.0), kbT(1.0), Kb(0.01), Kv(0.01), Ka(0.01), 
    //                    Rv(1.0), area_target(12.5664), volume_target(4.18879), rVol(0.95), Rp(0.5), u(0.5), 
    //                    U(0.01), rho(2.0), rc(20.0), r_equilibrium(1.0), mass_particle(0.1), total_mass_particle(10.0), 
    //                    force_residual(0.0), distance_threshold(0.1), forced_wrapping_fraction(0.05), v_smooth_flag(1), 
    //                    delaunay_tri_flag(1), mesh_reg_frequency(100) {
    //     setVertexCount(numV);
    // }

    void initialize(int vertexCount) {
        numV = vertexCount;
        // Resize and initialize matrices
        Force_Area.resize(numV, 3); Force_Volume.resize(numV, 3);
        Force_Bending.resize(numV, 3); Force_Adhesion.resize(numV, 3);
        Force_Random.resize(numV, 3); Force_Drag.resize(numV, 3);
        Force_Repulsion.resize(numV, 3); Force_Total.resize(numV, 3);
        acceleration.resize(numV, 3); acceleration_half_step.resize(numV, 3);
        velocity.resize(numV, 3); velocity_half_step.resize(numV, 3);
        // particle_velocities.resize(numV, 3); displace.resize(numV, 3);
        // Set matrices to zero
        Force_Area.setZero(); Force_Volume.setZero();
        Force_Bending.setZero(); Force_Adhesion.setZero();
        Force_Random.setZero(); Force_Drag.setZero();
        Force_Repulsion.setZero(); Force_Total.setZero();
        acceleration.setZero(); acceleration_half_step.setZero();
        velocity.setZero(); velocity_half_step.setZero();
        particle_velocities.setZero(); displace.setZero();

        rotation_matrix.setIdentity();
        idiag.setZero();
        torque.setZero();
        ang_momentum.setZero();
        ang_velocity.setZero();
        particle_velocity_com.setZero();
        particle_acceleration_com.setZero();
        center_of_mass.setZero();
        current_quaternion.setIdentity();

    }
};

#endif
