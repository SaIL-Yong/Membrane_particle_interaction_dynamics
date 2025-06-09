#include "simulation.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <Eigen/Dense> // For matrix operations
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/edge_lengths.h>
#include <igl/intrinsic_delaunay_triangulation.h>
#include <igl/signed_distance.h>
#include "energy.h"
#include "parameters.h"
#include "integrate.h"
#include "meshops.h"

#define PI 3.14159265358979323846
#define EPS 1.0e-10

//──────────────────────────────────────────────────────────────────────────────
// Verlet‐integration routine (membrane + optional vesicle + adhesion energetics)
//──────────────────────────────────────────────────────────────────────────────
void verlet_integration(SimulationData& sim_data, std::fstream& logfile) {
    // ─── Logfile Header ───────────────────────────────────────────────────────
    logfile << "Iteration  Time  "
               "Area1  Volume1  ReducedVolume1  BendingEnergy1  AreaEnergy1  VolumeEnergy1  "
               "AdhesionEnergy  PotentialEnergy1  KineticEnergy1  "
               "Area2  Volume2  ReducedVolume2  BendingEnergy2  AreaEnergy2  VolumeEnergy2  "
               "PotentialEnergy2  KineticEnergy2  "
               "TotalEnergy  EnergyChangeRate  ForceResidual"
            << std::endl;

    // ─── Console Header ───────────────────────────────────────────────────────
    std::cout << "Iteration  ReducedVolume1  BendingEnergy1  AdhesionEnergy  TotalEnergy  "
                 "EnergyChangeRate  ForceResidual"
              << std::endl;

    Mesh M1, M2;
    Energy E1, E2;
    int toln = 0;

    // Initial forces/energies at iteration 0
    calculate_forces(sim_data, M1, E1, M2, E2, /*current_iteration=*/0);

    for (int i = 0; i < sim_data.iterations; i++) {
        // ── First half‐step for Mesh 1 ────────────────────────────────────────
        sim_data.velocity1 += (sim_data.Force_Total1 / sim_data.mass1) * 0.5 * sim_data.dt;
        sim_data.V1       += sim_data.velocity1 * sim_data.dt;

        // ── First half‐step for Mesh 2 (if vesicle_flag=1) ────────────────────
        if (sim_data.vesicle_flag) {
            sim_data.velocity2 += (sim_data.Force_Total2 / sim_data.mass2) * 0.5 * sim_data.dt;
            sim_data.V2       += sim_data.velocity2 * sim_data.dt;
        }

        // ── Mesh Regularization for Mesh 1 ───────────────────────────────────
        if (sim_data.v_smooth_flag1 || sim_data.delaunay_tri_flag1) {
            if ((i + 1) % sim_data.mesh_reg_frequency1 == 0) {
                if (sim_data.v_smooth_flag1) {
                    sim_data.V1 = M1.vertex_smoothing(sim_data.V1, sim_data.F1);
                }
                if (sim_data.delaunay_tri_flag1) {
                    igl::edge_lengths(sim_data.V1, sim_data.F1, sim_data.l1);
                    igl::intrinsic_delaunay_triangulation(sim_data.l1, sim_data.F1, sim_data.l1, sim_data.F1);
                }
            }
        }

        // ── Mesh Regularization for Mesh 2 ───────────────────────────────────
        if (sim_data.vesicle_flag && (sim_data.v_smooth_flag2 || sim_data.delaunay_tri_flag2)) {
            if ((i + 1) % sim_data.mesh_reg_frequency2 == 0) {
                if (sim_data.v_smooth_flag2) {
                    sim_data.V2 = M2.vertex_smoothing(sim_data.V2, sim_data.F2);
                }
                if (sim_data.delaunay_tri_flag2) {
                    igl::edge_lengths(sim_data.V2, sim_data.F2, sim_data.l2);
                    igl::intrinsic_delaunay_triangulation(sim_data.l2, sim_data.F2, sim_data.l2, sim_data.F2);
                }
            }
        }

        // ── Recompute forces & energies for both meshes + adhesion ──────────
        calculate_forces(sim_data, M1, E1, M2, E2, i);

        // ── Second half‐step for Mesh 1 ──────────────────────────────────────
        sim_data.velocity1 += (sim_data.Force_Total1 / sim_data.mass1) * 0.5 * sim_data.dt;

        // ── Second half‐step for Mesh 2 ──────────────────────────────────────
        if (sim_data.vesicle_flag) {
            sim_data.velocity2 += (sim_data.Force_Total2 / sim_data.mass2) * 0.5 * sim_data.dt;
        }

        // ── Compute reduced volumes ─────────────────────────────────────────
        sim_data.rVol1 = 6.0 * sqrt(PI) * M1.volume_total * pow(M1.area_total, -1.5);
        if (sim_data.vesicle_flag) {
            sim_data.rVol2 = 6.0 * sqrt(PI) * M2.volume_total * pow(M2.area_total, -1.5);
        } else {
            sim_data.rVol2 = 0.0;
        }

        // ── Logging & Console Output (every logfrequency steps) ─────────────
        if (i % sim_data.logfrequency == 0) {
            // Energy change rate for Mesh 1 + Mesh 2 combined
            sim_data.EnergyChangeRate_log =
                (sim_data.EnergyTotal - sim_data.EnergyTotalOld_log)
                / (sim_data.logfrequency * sim_data.dt);
            sim_data.EnergyTotalOld_log = sim_data.EnergyTotal;
            sim_data.etol(toln++) = sim_data.EnergyChangeRate_log;

            // Console output (mesh 1-centric)
            std::cout << i << "  "
                      << sim_data.rVol1 << "  "
                      << sim_data.EnergyBending1 << "  "
                      << sim_data.EnergyAdhesion << "  "
                      << sim_data.EnergyTotal << "  "
                      << sim_data.EnergyChangeRate_log << "  "
                      << sim_data.force_residual << std::endl;

            // Logfile output (all quantities)
            logfile << i << "  "
                    << sim_data.time << "  "
                    // Mesh 1 geometry & energies
                    << M1.area_total << "  "
                    << M1.volume_total << "  "
                    << sim_data.rVol1 << "  "
                    << sim_data.EnergyBending1 << "  "
                    << sim_data.EnergyArea1 << "  "
                    << sim_data.EnergyVolume1 << "  "
                    // Adhesion energy
                    << sim_data.EnergyAdhesion << "  "
                    // Mesh 1 potential & kinetic
                    << sim_data.EnergyPotential1 << "  "
                    << sim_data.EnergyKinetic1 << "  ";

            if (sim_data.vesicle_flag) {
                // Mesh 2 geometry & energies
                logfile << M2.area_total << "  "
                        << M2.volume_total << "  "
                        << sim_data.rVol2 << "  "
                        << sim_data.EnergyBending2 << "  "
                        << sim_data.EnergyArea2 << "  "
                        << sim_data.EnergyVolume2 << "  "
                        // Mesh 2 potential & kinetic
                        << sim_data.EnergyPotential2 << "  "
                        << sim_data.EnergyKinetic2 << "  ";
            } else {
                // if no vesicle, pad with zeros
                logfile << 0.0 << "  "  // Area2
                        << 0.0 << "  "  // Volume2
                        << 0.0 << "  "  // ReducedVolume2
                        << 0.0 << "  "  // BendingEnergy2
                        << 0.0 << "  "  // AreaEnergy2
                        << 0.0 << "  "  // VolumeEnergy2
                        << 0.0 << "  "  // PotentialEnergy2
                        << 0.0 << "  "; // KineticEnergy2
            }

            // Total combined energy, change rate, residual
            logfile << sim_data.EnergyTotal << "  "
                    << sim_data.EnergyChangeRate_log << "  "
                    << sim_data.force_residual << std::endl;
        }

        // ── Convergence Check ────────────────────────────────────────────────
        if (i % sim_data.tolsteps == 0) {
            if (i != 0) {
                sim_data.EnergyChangeRate_avg =
                    sim_data.etol(Eigen::seq(toln - 1 - sim_data.tolmean_steps, toln - 1)).mean();
                if (std::abs(sim_data.EnergyChangeRate_avg) < sim_data.tolerance &&
                    sim_data.tolerance_flag) {
                    std::cout << "Energy change rate reaches threshold.\n"
                              << "Simulation at equilibrium.\n";
                    break;
                }
            }
        }

        // ── OFF dumps every dumpfrequency steps ─────────────────────────────
        if (i % sim_data.dumpfrequency == 0) {
            char dump_m1[128], dump_m2[128];
            sprintf(dump_m1, "dump%08d_mesh1.off", i);
            igl::writeOFF(dump_m1, sim_data.V1, sim_data.F1);

            if (sim_data.vesicle_flag) {
                sprintf(dump_m2, "dump%08d_mesh2.off", i);
                igl::writeOFF(dump_m2, sim_data.V2, sim_data.F2);
            }
        }

        sim_data.time += sim_data.dt;

        if ((i + 1) == sim_data.iterations) {
            std::cout << "Simulation reached max iterations.\n";
        }

        // ── Always overwrite final mesh outputs ───────────────────────────
        igl::writeOFF(sim_data.outFile1, sim_data.V1, sim_data.F1);
        if (sim_data.vesicle_flag) {
            igl::writeOFF(sim_data.outFile2, sim_data.V2, sim_data.F2);
        }
    }
}


//──────────────────────────────────────────────────────────────────────────────
// calculate_forces: compute bending/area/volume for Mesh 1 & Mesh 2 + adhesion
//──────────────────────────────────────────────────────────────────────────────
void calculate_forces(
    SimulationData& sim_data,
    Mesh&           M1,
    Energy&         E1,
    Mesh&           M2,
    Energy&         E2,
    int             current_iteration
) {
    // ─── MESH 1 (MEMBRANE) CALCULATIONS ──────────────────────────────────────
    M1.mesh_cal(sim_data.V1, sim_data.F1);

    // Bending force & energy (mesh 1)
    E1.compute_bendingenergy_force(
        sim_data.V1,
        sim_data.F1,
        sim_data.Kb1,
        sim_data.Force_Bending1,
        sim_data.EnergyBending1,
        M1
    );

    // Area force & energy (mesh 1)
    E1.compute_areaenergy_force(
        sim_data.V1,
        sim_data.F1,
        sim_data.Ka1,
        sim_data.area_target1,
        sim_data.Force_Area1,
        sim_data.EnergyArea1,
        M1
    );

    // Volume force & energy (mesh 1)
    E1.compute_volumeenergy_force(
        sim_data.V1,
        sim_data.F1,
        sim_data.Kv1,
        sim_data.volume_target1,
        sim_data.Force_Volume1,
        sim_data.EnergyVolume1,
        M1
    );

    // Random force on mesh 1 (if enabled)
    if (sim_data.random_force_flag1) {
        E1.compute_random_force(
            sim_data.gamma1,
            sim_data.kbT1,
            sim_data.mass1,
            sim_data.dt,
            sim_data.Force_Random1
        );
    }

    // Drag force on mesh 1
    E1.compute_drag_force(
        sim_data.velocity1,
        sim_data.gamma1,
        sim_data.mass1,
        sim_data.Force_Drag1
    );

    // ─── MESH 2 (VESICLE) CALCULATIONS (only if vesicle_flag) ───────────────
    if (sim_data.vesicle_flag) {
        M2.mesh_cal(sim_data.V2, sim_data.F2);

        // Bending force & energy (mesh 2)
        E2.compute_bendingenergy_force(
            sim_data.V2,
            sim_data.F2,
            sim_data.Kb2,
            sim_data.Force_Bending2,
            sim_data.EnergyBending2,
            M2
        );

        // Area force & energy (mesh 2)
        E2.compute_areaenergy_force(
            sim_data.V2,
            sim_data.F2,
            sim_data.Ka2,
            sim_data.area_target2,
            sim_data.Force_Area2,
            sim_data.EnergyArea2,
            M2
        );

        // Volume force & energy (mesh 2)
        E2.compute_volumeenergy_force(
            sim_data.V2,
            sim_data.F2,
            sim_data.Kv2,
            sim_data.volume_target2,
            sim_data.Force_Volume2,
            sim_data.EnergyVolume2,
            M2
        );

        // Random force on mesh 2 (if enabled)
        if (sim_data.random_force_flag2) {
            E2.compute_random_force(
                sim_data.gamma2,
                sim_data.kbT2,
                sim_data.mass2,
                sim_data.dt,
                sim_data.Force_Random2
            );
        }

        // Drag force on mesh 2
        E2.compute_drag_force(
            sim_data.velocity2,
            sim_data.gamma2,
            sim_data.mass2,
            sim_data.Force_Drag2
        );
    }

    // ─── ADHESION BETWEEN MESH 1 & MESH 2 (only if vesicle_flag) ────────────
    if (sim_data.vesicle_flag) {
        igl::signed_distance(
            sim_data.V1,
            sim_data.V2,
            sim_data.F2,
            igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,
            sim_data.signed_distance,
            sim_data.facet_index,
            sim_data.closest_points,
            sim_data.normals_closest_points
        );
        E1.compute_adhesion_energy_force(
            sim_data.V1,
            sim_data.F1,
            //sim_data.V2,
            sim_data.closest_points,
            sim_data.potential_range,
            sim_data.U,//adhesion_strength,
            sim_data.r_equilibrium,
            sim_data.rc,
            sim_data.angle_flag,
            sim_data.vesicle_position,
            sim_data.Ew_t,
            sim_data.Kw,
            sim_data.Force_Adhesion,
            sim_data.Force_Repulsion,
            sim_data.signed_distance,
            sim_data.EnergyAdhesion,
            sim_data.EnergyBias,
            M1
        );
    } else {
        // If no vesicle, zero out adhesion
        sim_data.Force_Adhesion.setZero(sim_data.numV1, 3);
        sim_data.EnergyAdhesion = 0.0;
        sim_data.EnergyBias = 0.0;
    }

    // ─── ASSEMBLE ENERGIES FOR MESH 1 ────────────────────────────────────────
    sim_data.EnergyPotential1 =
        sim_data.EnergyBending1 +
        sim_data.EnergyArea1 +
        sim_data.EnergyVolume1 +
        sim_data.EnergyAdhesion +
        sim_data.EnergyBias;
    sim_data.EnergyKinetic1 = 0.5 * sim_data.mass1
                            * (sim_data.velocity1.rowwise().squaredNorm().sum());

    // ─── ASSEMBLE ENERGIES FOR MESH 2 ────────────────────────────────────────
    if (sim_data.vesicle_flag) {
        sim_data.EnergyPotential2 =
            sim_data.EnergyBending2 +
            sim_data.EnergyArea2 +
            sim_data.EnergyVolume2;
        sim_data.EnergyKinetic2 = 0.5 * sim_data.mass2
                                * (sim_data.velocity2.rowwise().squaredNorm().sum());
    } else {
        sim_data.EnergyPotential2 = 0.0;
        sim_data.EnergyKinetic2   = 0.0;
    }

    // ─── GRAND TOTAL ENERGY ─────────────────────────────────────────────────
    sim_data.EnergyPotential = sim_data.EnergyPotential1 + sim_data.EnergyPotential2;
    sim_data.EnergyKinetic   = sim_data.EnergyKinetic1   + sim_data.EnergyKinetic2;
    sim_data.EnergyTotal     = sim_data.EnergyPotential + sim_data.EnergyKinetic;

    // ─── ASSEMBLE FORCES FOR MESH 1 ─────────────────────────────────────────
    sim_data.Force_Total1 =
        sim_data.Force_Bending1 +
        sim_data.Force_Area1 +
        sim_data.Force_Volume1 +
        sim_data.Force_Adhesion +
        sim_data.Force_Random1 +
        sim_data.Force_Drag1;

    // ─── ASSEMBLE FORCES FOR MESH 2 ─────────────────────────────────────────
    if (sim_data.vesicle_flag) {
        sim_data.Force_Total2 =
            sim_data.Force_Bending2 +
            sim_data.Force_Area2 +
            sim_data.Force_Volume2 +
            sim_data.Force_Random2 +
            sim_data.Force_Drag2;
    } else {
        sim_data.Force_Total2.setZero(sim_data.numV2, 3);
    }

    // ─── COMBINED RESIDUAL NORM (no need to literally add two different-sized matrices)
    sim_data.force_residual =
      std::sqrt(
        sim_data.Force_Total1.squaredNorm() +
        sim_data.Force_Total2.squaredNorm()
      );
}


void initialize_simulation(
    SimulationData& sim_data,
    const Parameter&      parameter,
    std::fstream&   logfile
) {
    // ─── Read Mesh 1 (membrane) and verify success ─────────────────────────
    bool ok1 = igl::readOFF(parameter.meshFile1, sim_data.V1, sim_data.F1);
    if (!ok1 || sim_data.V1.rows() == 0 || sim_data.F1.cols() != 3) {
        std::cerr << "ERROR: failed to load meshFile1 = '"
                  << parameter.meshFile1 << "'\n";
        std::exit(1);
    }

    // ─── Read Mesh 2 (vesicle) if requested and verify ────────────────────
    if (parameter.vesicle_flag) {
        bool ok2 = igl::readOFF(parameter.meshFile2, sim_data.V2, sim_data.F2);
        if (!ok2 || sim_data.V2.rows() == 0 || sim_data.F2.cols() != 3) {
            std::cerr << "ERROR: failed to load meshFile2 = '"
                      << parameter.meshFile2 << "'\n";
            std::exit(1);
        }
    } else {
        sim_data.numV2 = 0;
        sim_data.numF2 = 0;
    }

    // ─── Number of vertices/faces for mesh 1 ──────────────────────────────
    sim_data.numV1 = sim_data.V1.rows();
    sim_data.numF1 = sim_data.F1.rows();

    // ─── Number of vertices/faces for mesh 2 ──────────────────────────────
    if (parameter.vesicle_flag) {
        sim_data.numV2 = sim_data.V2.rows();
        sim_data.numF2 = sim_data.F2.rows();
    }

    // ─── Initialize per‐vertex arrays for mesh 1 & mesh 2 ─────────────────
    sim_data.initialize(sim_data.numV1);
    if (parameter.vesicle_flag) {
        sim_data.initialize_mesh2(sim_data.numV2);
    }

    // ─── Log file header ───────────────────────────────────────────────────
    logfile << "=== Simulation Log ===\n";
    logfile << "Mesh 1: " << parameter.meshFile1;
    if (parameter.vesicle_flag) {
        logfile << "  |  Mesh 2: " << parameter.meshFile2;
    }
    logfile << "\n\n";

    // ─── Copy common settings into sim_data ────────────────────────────────
    sim_data.iterations         = parameter.iterations;
    sim_data.logfrequency       = parameter.logfrequency;
    sim_data.dumpfrequency      = parameter.dumpfrequency;
    sim_data.resfrequency       = parameter.resfrequency;
    sim_data.mesh_reg_frequency1 = parameter.mesh_reg_frequency;
    sim_data.mesh_reg_frequency2 = parameter.mesh_reg_frequency;
    sim_data.bondfrequency      = parameter.bondfrequency;
    sim_data.dt                 = parameter.dt;
    sim_data.time               = 0.0;
    sim_data.dtf                = sim_data.dt / 2.0;
    sim_data.tolerance          = parameter.tolerance;
    sim_data.tolerance_flag     = parameter.tolerance_flag;
    sim_data.tolfrequency       = parameter.tolfrequency;
    sim_data.tolsteps           = static_cast<int>(std::floor(sim_data.tolfrequency / sim_data.dt));
    sim_data.tolmean_steps      = static_cast<int>(std::floor(sim_data.tolsteps / sim_data.logfrequency));
    sim_data.etol.resize(static_cast<int>(std::floor(sim_data.iterations / sim_data.logfrequency)));
    sim_data.etol.setZero();

    // ─── Print Mesh 1 info ─────────────────────────────────────────────────
    std::cout << "=== Mesh 1 (Membrane) Info ===\n";
    std::cout << "Vertices: " << sim_data.numV1
              << "   Faces: "  << sim_data.numF1 << "\n";
    std::cout << "Iterations: "   << sim_data.iterations << "\n";
    std::cout << "Logfreq: "      << sim_data.logfrequency << "   Dumpfreq: "
              << sim_data.dumpfrequency << "\n";
    std::cout << "Tol steps: "    << sim_data.tolsteps << "\n";
    std::cout << "Time step: "    << sim_data.dt << "\n\n";

    logfile << "--- Mesh 1 Info ---\n";
    logfile << "Vertices: " << sim_data.numV1 << "   Faces: " << sim_data.numF1 << "\n";
    logfile << "Iterations: " << sim_data.iterations << "\n";
    logfile << "Logfreq: " << sim_data.logfrequency
            << "   Dumpfreq: " << sim_data.dumpfrequency << "\n";
    logfile << "Tol steps: " << sim_data.tolsteps << "\n";
    logfile << "Time step: " << sim_data.dt << "\n\n";

    // ─── Print Mesh 2 info (if applicable) ────────────────────────────────
    if (parameter.vesicle_flag) {
        std::cout << "=== Mesh 2 (Vesicle) Info ===\n";
        std::cout << "Vertices: " << sim_data.numV2
                  << "   Faces: "  << sim_data.numF2 << "\n\n";

        logfile << "--- Mesh 2 Info ---\n";
        logfile << "Vertices: " << sim_data.numV2 << "   Faces: " << sim_data.numF2 << "\n\n";
    }

    // ─── Copy “mesh1” physical properties into sim_data ────────────────────
    sim_data.gamma1   = parameter.gamma;
    sim_data.mass1    = parameter.mass;
    sim_data.kbT1     = parameter.kbT;
    sim_data.Kb1      = parameter.Kb1;
    sim_data.Ka1      = parameter.Ka1;
    sim_data.Kv1      = 0.0;           // default
    sim_data.Rv1      = 1.0;           // will update if Kv1 > EPS
    sim_data.area_target1   = 4.0 * PI * sim_data.Rv1 * sim_data.Rv1;
    sim_data.volume_target1 = 0.0;

    if (std::abs(parameter.Kv1) > EPS) {
        sim_data.rVol_t1       = parameter.reduced_volume1;
        sim_data.Kv1           = parameter.Kv1;
        sim_data.volume_target1 =
            sim_data.rVol_t1 * (4.0 / 3.0) * PI * std::pow(sim_data.Rv1, 3);
    }

    // Initialize mesh 1 velocity
    sim_data.velocity1 = Eigen::MatrixXd::Zero(sim_data.numV1, 3);

    // ─── Copy “mesh2” (vesicle) physical properties into sim_data ─────────
    if (parameter.vesicle_flag) {
        sim_data.gamma2   = parameter.gamma*parameter.mass_vesicle_ratio;  // same drag as mesh 1
        sim_data.mass2    = parameter.mass*parameter.mass_vesicle_ratio;
        sim_data.kbT2     = parameter.kbT;
        sim_data.Kb2      = parameter.Kb2;
        sim_data.Ka2      = parameter.Ka2;
        sim_data.Kv2      = 0.0;
        //sim_data.Rv2      = 1.0;
        // ─── Compute actual surface area of mesh 2 ─────────────────────────
        Eigen::VectorXd doubleA2;
        igl::doublearea(sim_data.V2, sim_data.F2, doubleA2);
        double area2 = doubleA2.sum() * 0.5; // sum(doubleA2)/2

        // Compute Rv2 from that surface area:  A = 4π R²  ⇒  Rv2 = sqrt(A/(4π))
        sim_data.Rv2 = std::sqrt(area2 / (4.0 * PI));
        sim_data.area_target2   = 4.0 * PI * sim_data.Rv2 * sim_data.Rv2;
        sim_data.volume_target2 = 0.0;

        if (std::abs(parameter.Kv2) > EPS) {
            sim_data.rVol_t2       = parameter.reduced_volume2;
            sim_data.Kv2           = parameter.Kv2;
            sim_data.volume_target2 =
                sim_data.rVol_t2 * (4.0 / 3.0) * PI * std::pow(sim_data.Rv2, 3);
        }

        // Initialize mesh 2 velocity
        sim_data.velocity2 = Eigen::MatrixXd::Zero(sim_data.numV2, 3);
    }

    // ─── Flags for mesh 1 and mesh 2 ───────────────────────────────────────
    sim_data.random_force_flag1   = parameter.random_force_flag;
    sim_data.v_smooth_flag1       = parameter.vertex_smoothing_flag;
    sim_data.delaunay_tri_flag1   = parameter.delaunay_triangulation_flag;
    sim_data.mesh_reg_frequency1  = parameter.mesh_reg_frequency;

    if (parameter.vesicle_flag) {
        sim_data.random_force_flag2   = parameter.random_force_flag;
        sim_data.v_smooth_flag2       = parameter.vertex_smoothing_flag;
        sim_data.delaunay_tri_flag2   = parameter.delaunay_triangulation_flag;
        sim_data.mesh_reg_frequency2  = parameter.mesh_reg_frequency;
    } else {
        sim_data.random_force_flag2   = 0;
        sim_data.v_smooth_flag2       = 0;
        sim_data.delaunay_tri_flag2   = 0;
        sim_data.mesh_reg_frequency2  = 0;
    }


    // ─── Write mesh 1 properties to console & log ───────────────────────────
    std::cout << "=== Mesh 1 (Membrane) Properties ===\n";
    std::cout << "Drag coeff: " << sim_data.gamma1
              << "   Mass: " << sim_data.mass1 << "\n";
    std::cout << "Kb1: " << sim_data.Kb1 << "   Ka1: " << sim_data.Ka1 
              << "   Kv1: " << sim_data.Kv1 << "\n";
    std::cout << "Area target1: " << sim_data.area_target1
              << "   Volume target1: " << sim_data.volume_target1 << "\n\n";

    logfile << "--- Mesh 1 (Membrane) Properties ---\n";
    logfile << "Drag coeff: " << sim_data.gamma1 << "   Mass: " << sim_data.mass1 << "\n";
    logfile << "Kb1: " << sim_data.Kb1 << "   Ka1: " << sim_data.Ka1 
            << "   Kv1: " << sim_data.Kv1 << "\n";
    logfile << "Area target1: " << sim_data.area_target1
            << "   Volume target1: " << sim_data.volume_target1 << "\n\n";

    // ─── Write mesh 2 properties to console & log (if vesicle_flag) ─────────
    if (parameter.vesicle_flag) {
        std::cout << "=== Mesh 2 (Vesicle) Properties ===\n";
        std::cout << "Drag coeff: " << sim_data.gamma2 
                  << "   Mass: " << sim_data.mass2 << "\n";
        std::cout << "Kb2: " << sim_data.Kb2 << "   Ka2: " << sim_data.Ka2 
                  << "   Kv2: " << sim_data.Kv2 << "\n";
        std::cout << "Area target2: " << sim_data.area_target2
                  << "   Volume target2: " << sim_data.volume_target2 << "\n\n";

        logfile << "--- Mesh 2 (Vesicle) Properties ---\n";
        logfile << "Drag coeff: " << sim_data.gamma2 << "   Mass: " << sim_data.mass2 << "\n";
        logfile << "Kb2: " << sim_data.Kb2 << "   Ka2: " << sim_data.Ka2 
                << "   Kv2: " << sim_data.Kv2 << "\n";
        logfile << "Area target2: " << sim_data.area_target2
                << "   Volume target2: " << sim_data.volume_target2 << "\n\n";
    }

    // ─── Copy vesicle‐adhesion info into sim_data ────────────────────────────
    sim_data.vesicle_flag         = parameter.vesicle_flag;
    sim_data.vesicle_position     = parameter.vesicle_position;
    sim_data.vesicle_coord_flag   = parameter.vesicle_coord_flag;
    sim_data.vesicle_radius       = parameter.vesicle_radius;

    // Log vesicle-adhesion info
    logfile << "--- Vesicle-Adhesion Info ---\n";
    logfile << "vesicle_flag: " << sim_data.vesicle_flag << "\n";
    //logfile << "vesicle_position: " << sim_data.vesicle_position.transpose() << "\n";
    //logfile << "vesicle_coord_flag: " << sim_data.vesicle_coord_flag << "\n";
    //logfile << "vesicle_radius: " << sim_data.vesicle_radius << "\n\n";
    //sim_data.adhesion_strength    = parameter.adhesion_strength;
    sim_data.u    = parameter.adhesion_strength;
    sim_data.U   = (sim_data.Kb1 * sim_data.u) / (sim_data.Rv2 * sim_data.Rv2);
    sim_data.mass_vesicle_ratio   = parameter.mass_vesicle_ratio;
    sim_data.potential_range      = parameter.potential_range;
    sim_data.r_equilibrium        = parameter.r_equilibrium;
    sim_data.angle_flag           = parameter.angle_condition_flag;
    sim_data.forced_wrapping_flag = parameter.forced_wrapping_flag;
    sim_data.wrapping_fraction    = parameter.wrapping_fraction;
    sim_data.Kw                   = parameter.wrapping_bias_strength;

    // Log vesicle-adhesion parameters
    logfile << "adhesion_strength: " << sim_data.u << "\n";
    logfile << "U: " << sim_data.U << "\n";
    logfile << "mass_vesicle_ratio: " << sim_data.mass_vesicle_ratio << "\n";
    logfile << "potential_range: " << sim_data.potential_range << "\n";
    logfile << "r_equilibrium: " << sim_data.r_equilibrium << "\n";
    logfile << "angle_flag: " << sim_data.angle_flag << "\n";
    logfile << "forced_wrapping_flag: " << sim_data.forced_wrapping_flag << "\n";
    logfile << "wrapping_fraction: " << sim_data.wrapping_fraction << "\n";
    logfile << "wrapping_bias_strength (Kw): " << sim_data.Kw << "\n";

    // Precompute rc (10×ρ), Ew_t (–U×area) if forced_wrapping
    sim_data.rc = 10.0 * sim_data.potential_range;
    if (sim_data.forced_wrapping_flag && parameter.vesicle_flag) {
        Mesh tmpM;
        tmpM.mesh_cal(sim_data.V2, sim_data.F2);
        double vesicle_area = tmpM.area_total;
        double Area_w_t = sim_data.wrapping_fraction * vesicle_area;
        sim_data.Ew_t = - sim_data.U * Area_w_t;
    } else {
        sim_data.Ew_t = 0.0;
    }
    // ─── Copy output‐file names into sim_data ─────────────────────────────
    sim_data.outFile1 = parameter.outFile1;
    sim_data.outFile2 = parameter.outFile2;
    sim_data.resFile1 = parameter.resFile1;
    sim_data.resFile2 = parameter.resFile2;
    // ─── Final initialization message ───────────────────────────────────────
    logfile << "Initialization complete. Beginning integration...\n\n";
    std::cout << "Initialization complete. Beginning integration...\n\n";
}
