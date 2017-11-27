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
/**
 * @brief waypoint_t a waypoint is composed of a  6*3 matrix that depend
 * on the variable x, and of a 6d vector independent of x, such that
 * each control point of the target bezier curve is given by pi = wix * x + wis
 */
typedef std::pair<matrix6_t, point6_t> waypoint_t;
typedef const Eigen::Ref<const point_t>& point_t_tC;
typedef spline::bezier_curve  <double, double, 6, true, point6_t> bezier6_t;

BEZIER_COM_TRAJ_DLLAPI Matrix3 skew(point_t_tC x);
BEZIER_COM_TRAJ_DLLAPI waypoint_t initwp();
/**
 * @brief ComputeDiscretizedWaypoints Given the waypoints defining a bezier curve,
 * computes a discretization of the curve
 * @param wps original waypoints
 * @param bernstein berstein polynoms for
 * @param numSteps desired number of wayoints
 * @return a vector of waypoint representing the discretization of the curve
 */
BEZIER_COM_TRAJ_DLLAPI  std::vector<waypoint_t> ComputeDiscretizedWaypoints(const std::vector<waypoint_t>& wps, const std::vector<spline::Bern<double> >& bernstein, int numSteps);

/**
 * @brief compute6dControlPointInequalities Given linear and angular control waypoints,
 * compute the inequality matrices A and b, A x <= b that constrain the desired control point x.
 * @param cData data for the current contact phase
 * @param wps waypoints or the linear part of the trajectory
 * @param wpL waypoints or the angular part of the trajectory
 * @param useAngMomentum whether the angular momentum is consider or equal to 0
 * @param fail set to true if problem is found infeasible
 * @return
 */
BEZIER_COM_TRAJ_DLLAPI  std::pair<MatrixXX, VectorX> compute6dControlPointInequalities(const ContactData& cData, const std::vector<waypoint_t>& wps, const std::vector<waypoint_t>& wpL, const bool useAngMomentum, bool& fail);


/**
 * @brief solve x' h x + 2 g' x, subject to A*x <= b using quadprog
 * @param A Inequality matrix
 * @param b Inequality vector
 * @param H Cost matrix
 * @param g cost Vector
 * @return
 */
BEZIER_COM_TRAJ_DLLAPI ResultData solve(Cref_matrixXX A, Cref_vectorX b, Cref_matrixXX H, Cref_vectorX g, Cref_vectorX initGuess);


} // end namespace bezier_com_traj

#endif
