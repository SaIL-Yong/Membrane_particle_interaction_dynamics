#include "meshedparticle.h"
#include "meshops.h"
#include "energy.h"
#include <cstdlib>
#include <ctime>





//Temporary Bond Creation based on surface distance
void ParticleAdhesion::find_pairs(Eigen::MatrixXd V,Eigen::MatrixXi F, Eigen::MatrixXd V_particle, Eigen::VectorXd signed_distance,
                /**
                 * @brief The maximum distance between two particles for them to form a bond.
                 * 
                 * If the distance between two particles is less than or equal to this threshold,
                 * they will form a bond.
                 */
                const double distance_threshold, Eigen::Vector3d center_of_mass,std::vector<std::pair<int, int>>& bonds)
{
    bonds.clear(); // Clear any existing bonds since you don't want to keep them
    //double radius = 0.2;  // Adjust the radius as needed

    // Create a matrix to store vertices within the specified radius
    // Eigen::MatrixXd vertices_within_radius(0, 3);  // Initialize as an empty matrix

    // // Iterate through the vertices and check if they are within the radius
    // for (int i = 0; i < V.rows(); i++) {
    //     Eigen::Vector3d vertex = V.row(i);
    //     if (distance(vertex, center_point) <= radius) {
    //         // Append the vertex to the matrix
    //         vertices_within_radius.conservativeResize(vertices_within_radius.rows() + 1, 3);
    //         vertices_within_radius.row(vertices_within_radius.rows() - 1) = vertex;
    //     }
    // }



    for (int i = 0; i < V.rows(); i++) {
        double min_distance = std::numeric_limits<double>::max();
        int nearest_index = -1;

        for (int j = 0; j < V_particle.rows(); j++) {
            double distance = sqrt(pow(V(i, 0) - V_particle(j, 0), 2)
                + pow(V(i, 1) - V_particle(j, 1), 2)
                + pow(V(i, 2) - V_particle(j, 2), 2));

            if (distance < min_distance) {
                min_distance = distance;
                nearest_index = j;
            }
        }

        if (nearest_index != -1 && min_distance < distance_threshold) {
            bonds.emplace_back(i, nearest_index); // Generate a new bond based on distance
        }
    }

    std::ofstream outputFile("bonds.txt");

    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    for (const auto& bond : bonds) {
        int vertex1_idx = bond.first;
        int vertex2_idx = bond.second;

        // Write the bond information to the file
        double distance = sqrt(pow(V(vertex1_idx, 0) - V_particle(vertex2_idx, 0), 2)
            + pow(V(vertex1_idx, 1) - V_particle(vertex2_idx, 1), 2)
            + pow(V(vertex1_idx, 2) - V_particle(vertex2_idx, 2), 2));

        outputFile << "Bond: (" << vertex1_idx << ", " << vertex2_idx << "), Distance: " << distance << std::endl;
    }

    outputFile.close();

    //Eigen::Vector3d center_of_mass;
    //igl::centroid(V, F, center_of_mass);
    //std::cout << "COM: " << center_of_mass << std::endl;


    /* random point generator
    srand(static_cast<unsigned>(time(0))); // Seed the random number generator

    int num_points = 100;
    Eigen::MatrixXd random_points(num_points, 3);

    for (int i = 0; i < num_points; i++) {
        // Select a random vertex
        int random_vertex = rand() % V.rows(); // Assuming V is a matrix of vertex coordinates

        // Generate random offsets for x, y, and z
        double x_offset = (rand() % 2000 - 1000) / 10000.0; // Adjust the range as needed
        double y_offset = (rand() % 2000 - 1000) / 10000.0;
        double z_offset = (rand() % 2000 - 1000) / 10000.0;

        // Calculate the random point
        random_points.row(i) = V.row(random_vertex) + Eigen::RowVector3d(x_offset, y_offset, z_offset);
    }
    */


    //saving in a text file
    // std::ofstream outputFile1("random_points.txt");

    // if (outputFile1.is_open()) {
    //     for (int i = 0; i < num_points; i++) {
    //         outputFile1 << random_points(i, 0) << " " << random_points(i, 1) << " " << random_points(i, 2) << std::endl;
    //     }

    //     outputFile1.close();
    //     std::cout << "Random points saved to random_points.txt." << std::endl;
    // } else {
    //     std::cerr << "Unable to open the output file." << std::endl;
        
    // }




}





//Temporary Bond Creation
/*void ParticleAdhesion::find_pairs(Eigen::MatrixXd V,Eigen::MatrixXi F, Eigen::MatrixXd V_particle, double distance_threshold, std::vector<std::pair<int, int>>& bonds)
{
    bonds.clear(); // Clear any existing bonds since you don't want to keep them

    for (int i = 0; i < V.rows(); i++) {
        double min_distance = std::numeric_limits<double>::max();
        int nearest_index = -1;

        for (int j = 0; j < V_particle.rows(); j++) {
            double distance = sqrt(pow(V(i, 0) - V_particle(j, 0), 2)
                + pow(V(i, 1) - V_particle(j, 1), 2)
                + pow(V(i, 2) - V_particle(j, 2), 2));

            if (distance < min_distance) {
                min_distance = distance;
                nearest_index = j;
            }
        }

        if (nearest_index != -1 && min_distance < distance_threshold) {
            bonds.emplace_back(i, nearest_index); // Generate a new bond based on distance
        }
    }

    std::ofstream outputFile("bonds.txt");

    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    for (const auto& bond : bonds) {
        int vertex1_idx = bond.first;
        int vertex2_idx = bond.second;

        // Write the bond information to the file
        double distance = sqrt(pow(V(vertex1_idx, 0) - V_particle(vertex2_idx, 0), 2)
            + pow(V(vertex1_idx, 1) - V_particle(vertex2_idx, 1), 2)
            + pow(V(vertex1_idx, 2) - V_particle(vertex2_idx, 2), 2));

        outputFile << "Bond: (" << vertex1_idx << ", " << vertex2_idx << "), Distance: " << distance << std::endl;
    }

    outputFile.close();
}

*/

void ParticleAdhesion::remove_long_bonds(std::vector<std::pair<int, int>>& bonds, Eigen::MatrixXd& V, Eigen::MatrixXd& V_particle, double max_bond_length)
{
    // Create a new vector to store the updated bonds
    std::vector<std::pair<int, int>> updated_bonds;

    for (const auto& bond : bonds) {
        int vertex1_idx = bond.first;
        int vertex2_idx = bond.second;

        // Calculate the bond length between the two vertices
        double bond_length = (V.row(vertex1_idx) - V_particle.row(vertex2_idx)).norm();

        // Check if the bond length is less than or equal to the maximum allowed length
        if (bond_length <= max_bond_length) {
            // If the bond length is within the threshold, keep this bond
            updated_bonds.push_back(bond);
        }
        // If the bond length is greater than the threshold, skip this bond (delete it).
    }

    // Replace the original bonds vector with the updated one
    bonds = updated_bonds;
}

