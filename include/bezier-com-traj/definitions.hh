/*
 * Copyright 2018, LAAS-CNRS
 * Author: Steve Tonneau
 */

#ifndef BEZIER_COM_TRAJ_DEFINITIONS_H
#define BEZIER_COM_TRAJ_DEFINITIONS_H

#include <centroidal-dynamics-lib/centroidal_dynamics.hh>
#include <spline/bezier_curve.h>
#include <Eigen/Dense>

namespace bezier_com_traj
{

typedef double value_type;
typedef Eigen::Matrix <value_type, 3, 3>                           Matrix3;
typedef Eigen::Matrix <value_type, 6, 3>                           Matrix63;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 3>              MatrixX3;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, Eigen::Dynamic> MatrixXX;
typedef centroidal_dynamics::Vector3 Vector3;
typedef centroidal_dynamics::Vector6 Vector6;
typedef centroidal_dynamics::VectorX VectorX;

typedef Eigen::Ref<Vector3>     Ref_vector3;
typedef Eigen::Ref<VectorX>     Ref_vectorX;
typedef Eigen::Ref<MatrixX3>    Ref_matrixX3;
typedef Eigen::Ref<MatrixXX>    Ref_matrixXX;

typedef const Eigen::Ref<const Vector3>     & Cref_vector3;
typedef const Eigen::Ref<const Vector6>     & Cref_vector6;
typedef const Eigen::Ref<const VectorX>     & Cref_vectorX;
typedef const Eigen::Ref<const MatrixXX>    & Cref_matrixXX;
typedef const Eigen::Ref<const MatrixX3>    & Cref_matrixX3;

typedef Matrix63 matrix6_t;
typedef Vector6 point6_t;
typedef Matrix3 matrix3_t;
typedef Vector3 point3_t;

typedef Eigen::Vector3d point_t;
typedef const Eigen::Ref<const point_t>& point_t_tC;

typedef spline::bezier_curve  <double, double, 3, true, point_t > bezier_t;
typedef spline::bezier_curve  <double, double, 6, true, point6_t> bezier6_t;

typedef std::vector< std::pair<double, int> > T_time;
typedef T_time::const_iterator CIT_time;

/**
* @brief waypoint_t a waypoint is composed of a  6*3 matrix that depend
* on the variable x, and of a 6d vector independent of x, such that
* each control point of the target bezier curve is given by pi = wix * x + wis
*/
typedef std::pair<matrix6_t, point6_t> waypoint6_t;
typedef std::pair<matrix3_t, point3_t> waypoint3_t;
typedef std::pair<double,point3_t> coefs_t;

} // end namespace bezier_com_traj

#endif
