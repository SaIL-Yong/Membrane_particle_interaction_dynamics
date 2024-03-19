// #include "RigidBody.h"
// include <Eigen/Core>
// #include <igl/massmatrix.h>
// #include <igl/doublearea.h>




// // // using namespace Eigen;

// // // MatrixXd V; // Vertex positions
// // // MatrixXi F; // Triangle indices

// // // double density = 1.0; // Density of the material
// // // double m; // Mass of the object
// // // Vector3d center_of_mass; // Center of mass of the object
// // // Matrix3d I; // Inertia tensor of the object

// // // // Read in the mesh
// // // igl::readOBJ("mesh.obj", V, F);

// // // // Compute the mass properties
// // // double volume = igl::volume(V, F) / 3.0; // Volume of the object
// // // m = density * volume; // Compute the mass of the object
// // // igl::centroid(V, F, center_of_mass); // Compute the center of mass of the object

// // // // Compute the inertia tensor
// // // I.setZero(); // Initialize the inertia tensor to zero
// // // for (int i = 0; i < F.rows(); i++) {
// // //     Vector3d v0 = V.row(F(i, 0));
// // //     Vector3d v1 = V.row(F(i, 1));
// // //     Vector3d v2 = V.row(F(i, 2));
// // //     Vector3d c = (v0 + v1 + v2) / 3.0 - center_of_mass;
// // //     double area = igl::doublearea(v0, v1, v2);
// // //     I(0,0) += area * (c(1)*c(1) + c(2)*c(2));
// // //     I(1,1) += area * (c(0)*c(0) + c(2)*c(2));
// // //     I(2,2) += area * (c(0)*c(0) + c(1)*c(1));
// // //     I(0,1) -= area * c(0) * c(1);
// // //     I(0,2) -= area * c(0) * c(2);
// // //     I(1,2) -= area * c(1) * c(2);
// // // }
// // // I(1,0) = I(0,1);
// // // I(2,0) = I(0,2);
// // // I(2,1) = I(1,2);
// // // I *= density / 60.0;

// // // // Output the inertia tensor
// // // std::cout << "Inertia tensor:" << std::endl << I << std::endl;

// // // Eigen::Quaterniond RigidBody::q()
// // // {
// // //     return Eigen::Quaterniond();
// // // }
