#pragma once
#include <iostream>
#include <cmath>
#include <igl/readOFF.h>
#include <igl/doublearea.h>
#include <igl/edges.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_face_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>
#include <igl/unique_edge_map.h>
#include <igl/doublearea.h>
#include <igl/oriented_facets.h>
#include <igl/edge_flaps.h>
#include <igl/edges.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Geometry>
#include <vector>
#include <igl/adjacency_matrix.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/gaussian_curvature.h>
#include <igl/intrinsic_delaunay_triangulation.h>
#include <igl/barycenter.h>



class Mesh {
           public:

           double cal_volume(Eigen::MatrixXd V,Eigen::MatrixXi F);
           //double AreaEnergy(Eigen::MatrixXd V,Eigen::MatrixXi F);
           //Eigen::Eigen::MatrixXd V
           Eigen::MatrixXd area_grad(Eigen::MatrixXd V,Eigen::MatrixXi F);
           Eigen::MatrixXd volume_grad(Eigen::MatrixXd V,Eigen::MatrixXi F);
           Eigen::MatrixXd vertex_smoothing(Eigen::MatrixXd V,Eigen::MatrixXi F);

};
