/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef BEZIER_COM_TRAJ_LIB_COMMON_SOLVE_H
#define BEZIER_COM_TRAJ_LIB_COMMON_SOLVE_H

#include <Eigen/Dense>
#include <bezier-com-traj/config.hh>
#include <bezier-com-traj/data.hh>


namespace bezier_com_traj
{

typedef Matrix63 matrix6_t;
typedef Vector6 point6_t;
typedef std::pair<matrix6_t, point6_t> waypoint_t;

typedef const Eigen::Ref<const point_t>& point_t_tC;

typedef spline::bezier_curve  <double, double, 6, true, point6_t> bezier6_t;

Matrix3 skew(point_t_tC x);
waypoint_t initwp();
std::vector<waypoint_t> ComputeDiscretizedWaypoints(const std::vector<waypoint_t>& wps, const std::vector<spline::Bern<double> >& bernstein, int numSteps);
std::pair<MatrixXX, VectorX> compute6dControlPointInequalities(const ContactData& cData, const std::vector<waypoint_t>& wps, const std::vector<waypoint_t>& wpL, const bool useAngMomentum);


/**
 * @brief solve x' h x + 2 g' x, subject to A*x <= b using quadprog
 * @param A Inequality matrix
 * @param b Inequality vector
 * @param H Cost matrix
 * @param g cost Vector
 * @return
 */
ResultData solve(Cref_matrixXX A, Cref_vectorX b, Cref_matrixXX H, Cref_vectorX g, Cref_vectorX initGuess);


} // end namespace bezier_com_traj

#endif
