// energy.cpp

#include <iostream>
#include <cmath>
#include <random>
#include <igl/barycentric_coordinates.h>
#include "energy.h"
#include "meshops.h"

// ─── compute_bendingenergy_force ───────────────────────────────────────────────────
void Energy::compute_bendingenergy_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    double                 Kb,
    Eigen::MatrixXd&       Force_Bending,
    double&                bending_energy,
    Mesh                   m
) {
    int nV = V.rows();
    bending_energy = 0.0;
    Force_Bending.setZero(nV, 3);

    // EB = 2·Kb·(H²)·area_voronoi → scalar per vertex
    EB = 2.0 * Kb * (m.H_squared.transpose() * m.area_voronoi).diagonal();
    bending_energy = EB.sum();

    // Lap_H = Minv · (L · H_signed)
    Lap_H = m.Minv * (m.L * m.H_signed);

    // force_density = 2·H_signed·(H_squared – K) + Lap_H
    force_density = (2.0 * m.H_signed.array() * (m.H_squared - m.K).array())
                    + Lap_H.array();

    // vector_term = force_density · area_voronoi
    vector_term = force_density.array() * m.area_voronoi.array();

    // Force_Bending = 2·Kb · (V_normals ∘ vector_term)
    //   (i.e. multiply each normal by the scalar "vector_term[i]")
    Force_Bending = (2.0 * Kb)
                  * (m.V_normals.array().colwise() * vector_term.array());
}


// ─── compute_areaenergy_force ─────────────────────────────────────────────────────
void Energy::compute_areaenergy_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    double                 Ka,
    double                 area_target,
    Eigen::MatrixXd&       Force_Area,
    double&                area_energy,
    Mesh                   m
) {
    int nV = V.rows();
    area_energy = 0.0;
    Force_Area.setZero(nV, 3);

    // da = current_area − target_area
    da = m.area_total - area_target;

    // area_energy = Ka·(da² / area_target)
    area_energy = Ka * da * da / area_target;

    // scalar_term = −2·Ka·da / area_target
    scalar_term = -2.0 * Ka * da / area_target;

    // AG = ∇(area)  (nV×3 matrix)
    AG = m.area_grad(V, F);

    // Force_Area = scalar_term · AG
    Force_Area = scalar_term * AG;
}


// ─── compute_volumeenergy_force ───────────────────────────────────────────────────
void Energy::compute_volumeenergy_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    double                 Kv,
    double                 volume_target,
    Eigen::MatrixXd&       Force_Volume,
    double&                volume_energy,
    Mesh                   m
) {
    int nV = V.rows();
    volume_energy = 0.0;
    Force_Volume.setZero(nV, 3);

    if (std::abs(Kv) > EPS) {
        // VG = ∇(volume)  (nV×3)
        VG = m.volume_grad(V, F);

        // dv = current_volume − target_volume
        dv = m.volume_total - volume_target;

        // volume_energy = Kv·(dv² / volume_target)
        volume_energy = Kv * dv * dv / volume_target;

        // scalar_term = −2·Kv·dv / volume_target
        scalar_term = -2.0 * Kv * dv / volume_target;

        // Force_Volume = scalar_term · VG
        Force_Volume = scalar_term * VG;
    }
}


// ─── compute_adhesion_energy_force ────────────────────────────────────────────────
// Note: we now take `closest_points` (nV×3), not V_vesicle.  “signed_distance” is also nV×1.
void Energy::compute_adhesion_energy_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& closest_points,
    double                 rho,
    double                 U,
    double                 r_equilibrium,
    double                 rc,
    int                    angle_flag,
    int                    vesicle_position,   // (unused in this version)
    double                 Ew_t,
    double                 Kw,
    Eigen::MatrixXd&       Force_Adhesion,
    Eigen::MatrixXd&       Force_Repulsion,
    const Eigen::MatrixXd& signed_distance,
    double&                EnergyAdhesion,
    double&                EnergyBias,
    Mesh                   m
) {
    int nV = V.rows();
    int nC = closest_points.rows();
    if (nV != nC) {
        std::cerr << "ERROR in compute_adhesion_energy_force():\n"
                  << "    V.rows() = " << nV << "  but closest_points.rows() = " << nC << "\n"
                  << "    That means signed_distance was called with mismatched arguments.\n";
        std::exit(1);
    }

    Force_Adhesion.setZero(nV, 3);

    coefficient.resize(nV);
    coefficient_derivative_x.resize(nV);
    coefficient_derivative_y.resize(nV);
    coefficient_derivative_z.resize(nV);
    distance.resize(nV);

    coefficient.setZero();
    coefficient_derivative_x.setZero();
    coefficient_derivative_y.setZero();
    coefficient_derivative_z.setZero();
    distance.setZero();
    // lennard_jones_force.setZero(nV, 3);
    // lennard_jones_potential.setZero(nV);

   // 1) Compute a Morse-like coefficient at each vertex (only if |dc| ≤ rc)
    for (int i = 0; i < nV; i++) {
        double dc = signed_distance(i);
        if (std::abs(dc) > rc) {
            continue; // no adhesion beyond cutoff
        }

        // Vector from the vesicle SURFACE (closest_points[i]) to membrane vertex V[i]
        Eigen::RowVector3d diff = V.row(i) - closest_points.row(i);

        // // Angle between “diff” and the membrane normal at i:
        // double numerator = diff.dot(m.V_normals.row(i));

        // double denom     = diff.norm() * m.V_normals.row(i).norm();
        // double angle     = 0.0;
        // if (denom > 0.0) {
        //     angle = std::acos(numerator / denom);
        // }

        // if (angle_flag) {
        //     // Only allow adhesion if normal is roughly “facing” the vesicle
        //     if (dc > 0 && angle <= 0.5 * igl::PI) {
        //         continue;
        //     }
        // }

        // Components of “diff”
        dx = diff.x();
        dy = diff.y();
        dz = diff.z();

        // Morse‐like coefficient:
        coefficient(i) = U * (std::exp(-(2.0 * dc) / rho)
                            - 2.0 * std::exp(-dc / rho));

        // Derivatives w.r.t. x,y,z (guard against dc=0)
        if (std::abs(dc) > EPS) {
            double base   = -std::exp(-(2.0 * dc) / rho)
                            + std::exp(-dc / rho);
            double factor = U / (dc * rho);
            coefficient_derivative_x(i) = factor * base * 2.0 * dx;
            coefficient_derivative_y(i) = factor * base * 2.0 * dy;
            coefficient_derivative_z(i) = factor * base * 2.0 * dz;
        } else {
            coefficient_derivative_x(i) = 0.0;
            coefficient_derivative_y(i) = 0.0;
            coefficient_derivative_z(i) = 0.0;
        }
    }

    // 2) Assemble coefficient_of_derivative (nV×3)
    coefficient_of_derivative.resize(nV, 3);
    coefficient_of_derivative.col(0) = coefficient_derivative_x.transpose();
    coefficient_of_derivative.col(1) = coefficient_derivative_y.transpose();
    coefficient_of_derivative.col(2) = coefficient_derivative_z.transpose();

    // 3) Get area gradient once (nV×3)
    AG = m.area_grad(V, F);

    // 4) First_Term = –(AG ∘ coefficient)   => (nV×3)
    Eigen::MatrixXd tmp1 = -(AG.array().colwise() * coefficient.array());

    // 5) Second_Term = –(coefficient_of_derivative ∘ area_voronoi)  => (nV×3)
    Eigen::MatrixXd tmp2 = -(coefficient_of_derivative.array()
                             .colwise() * m.area_voronoi.array());

    // 6) Sum = tmp1 + tmp2  (both are nV×3, so no dimension mismatch)
    Sum = tmp1 + tmp2;

    // 7) Adhesion energy: Ead = coefficient ∘ area_voronoi  (length nV)
    Ead = coefficient.array() * m.area_voronoi.array();
    EnergyAdhesion = Ead.sum();

    if (std::abs(Kw) > EPS) {
        dEw = EnergyAdhesion - Ew_t;
        EnergyBias = 0.5 * Kw * dEw * dEw;
        Force_Biased = Kw * dEw * Sum;
        Force_Adhesion = Sum + Force_Biased;
    } else {
        EnergyBias    = 0.0;
        Force_Adhesion = Sum;
    }

    // 8) The purely‐repulsive part is tmp2
    Force_Repulsion = tmp2;
}


// ─── redistributeAdhesionForce ───────────────────────────────────────────────────
void Energy::redistributeAdhesionForce(
    const Eigen::MatrixXd& V2,
    const Eigen::MatrixXi& F2,
    const Eigen::MatrixXd& closest_points,
    const Eigen::MatrixXd& Force_Adhesion,
    const Eigen::VectorXi& facet_index,
    Eigen::MatrixXd&       ForcesOnVertices
) {
    int nV2 = V2.rows();
    ForcesOnVertices.setZero(nV2, 3);

    for (int i = 0; i < closest_points.rows(); ++i) {
        // Negate for redistribution
        Eigen::RowVector3d f = -Force_Adhesion.row(i);
        if (f.norm() <= 1e-10) continue;

        int f_idx = facet_index(i);
        Eigen::Vector3i face = F2.row(f_idx);

        // Build the 3×3 triangle
        Eigen::Matrix<double, 3, 3> triangle;
        triangle.row(0) = V2.row(face(0));
        triangle.row(1) = V2.row(face(1));
        triangle.row(2) = V2.row(face(2));

        Eigen::RowVector3d p = closest_points.row(i);
        Eigen::Matrix<double, 1, 3> barycentricCoords;
        igl::barycentric_coordinates(
            p,
            triangle.row(0),
            triangle.row(1),
            triangle.row(2),
            barycentricCoords
        );

        for (int j = 0; j < 3; ++j) {
            ForcesOnVertices.row(face(j))
              += f * barycentricCoords(0, j);
        }
    }
}


// ─── compute_random_force ───────────────────────────────────────────────────────
void Energy::compute_random_force(
    double           gamma,
    double           kbT,
    double           mass,
    double           dt,
    Eigen::MatrixXd& Force_Random
) {
    std::normal_distribution<double> dist(0.0, 1.0);
    double scale = std::sqrt(2 * gamma * kbT * mass / dt);
    Force_Random = Eigen::MatrixXd::NullaryExpr(
        Force_Random.rows(),
        Force_Random.cols(),
        [&]() { return scale * dist(gen); }
    );
}


// ─── compute_drag_force ─────────────────────────────────────────────────────────
void Energy::compute_drag_force(
    const Eigen::MatrixXd& velocity,
    double                 gamma,
    double                 mass,
    Eigen::MatrixXd&       Force_Drag
) {
    Force_Drag = -gamma * velocity;
}


// ─── compute_drag_force_vesicle ─────────────────────────────────────────────────
void Energy::compute_drag_force_vesicle(
    const Eigen::Vector3d& velocity,
    double                 gamma,
    double                 mass,
    Eigen::Vector3d&       Force_Drag_vesicle
) {
    Force_Drag_vesicle = -gamma * velocity;
}
