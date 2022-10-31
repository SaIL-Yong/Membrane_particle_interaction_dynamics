#pragma once
#include<Eigen/Dense>
#include<iostream>
#include<string>
#include<vector>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



/*
Mesh class
A mesh is a collection of triangles. 
The mesh class contains the following information:
1) a vector of vertices 
2) a vector of normal vectors on the faces
3) a vector of coordinate transformation operators, including local to global, global to local

4) a vector of face adjacency information. Each face is guaranteed to have 3 neighbors.
5) a vector of inverse face adjacency information

6) Edge face information and convention

                 **
                 * *
     f:1, edge 1 *  *  f: 2, edge: 1
                 *   *
                 ******
                 f:0, edge:2
*/


struct Mesh {
    int numV;                                               // number of vertices
    int numF;                                               // number of faces
    double area_total;                                      // total area of triangle mesh
    double area_avg;                                        // average area of each triangle mesh
    Eigen::MatrixXd V;                                      // matrix storing vertice coordinates
    Eigen::MatrixXi F;                                      // matrix storing face information, every face is one row with three integers
    Eigen::MatrixXi TT, TTi;                                // face adjacency information
    Eigen::MatrixXd F_normals;                              // face normals
    Eigen::VectorXd dblA;                                   // store face edge information

    std::vector<Eigen::Matrix3d> Face_Edges;                // edge information for each faces

    
    
    
   
    
    

    // read mesh file
    void readMeshFile(std::string filename);                                
    
    // initialization
    void initialize();

    

};
