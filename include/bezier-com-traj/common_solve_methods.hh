/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef BEZIER_COM_TRAJ_LIB_COMMON_SOLVE_H
#define BEZIER_COM_TRAJ_LIB_COMMON_SOLVE_H

#include <Eigen/Dense>
#include <bezier-com-traj/config.hh>
#include <bezier-com-traj/data.hh>


namespace centroidal_dynamics
{

typedef Matrix63 matrix6_t;
typedef Vector6 point6_t;
typedef std::pair<matrix6_t, point6_t> waypoint_t;

typedef const Eigen::Ref<const point_t>& point_t_tC;

typedef spline::bezier_curve  <double, double, 6, true, point6_t> bezier6_t;

Matrix3 skew(point_t_tC x);
waypoint_t initwp();
std::vector<waypoint_t> ComputeDiscretizedWaypoints(const std::vector<waypoint_t>& wps, const std::vector<spline::Bern<double> >& bernstein, int numSteps);
std::pair<MatrixXX, VectorX> compute6dControlPointInequalities(const ContactData& cData, const std::vector<waypoint_t>& wps, const std::vector<waypoint_t>& wpL, const bool useAngMomentum, double T, double timeStep);

/**
 * @brief solve Qp problem 0.5*||D*x - d||^2, subject to A*x <= b
 * @param A Inequality matrix
 * @param b Inequality vecto
 * @param D Cost matrix
 * @param d cost Vector
 * @return
 */
ResultData solve(Cref_matrixXX A, Cref_vectorX b, Cref_matrixXX D, Cref_vectorX d);


} // end namespace centroidal_dynamics

#endif
