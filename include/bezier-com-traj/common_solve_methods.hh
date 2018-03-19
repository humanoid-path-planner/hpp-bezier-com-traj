/*
 * Copyright 2018, LAAS-CNRS
 * Author: Steve Tonneau
 */

#ifndef BEZIER_COM_TRAJ_LIB_COMMON_SOLVE_H
#define BEZIER_COM_TRAJ_LIB_COMMON_SOLVE_H

#include <Eigen/Dense>
#include <bezier-com-traj/config.hh>
#include <bezier-com-traj/data.hh>


namespace bezier_com_traj
{


BEZIER_COM_TRAJ_DLLAPI Matrix3 skew(point_t_tC x);
template<typename T> T initwp();
int Normalize(Ref_matrixXX A, Ref_vectorX b);


/**
 * @brief ComputeDiscretizedWaypoints Given the waypoints defining a bezier curve,
 * computes a discretization of the curve
 * @param wps original waypoints
 * @param bernstein berstein polynoms for
 * @param numSteps desired number of wayoints
 * @return a vector of waypoint representing the discretization of the curve
 */
BEZIER_COM_TRAJ_DLLAPI  std::vector<waypoint6_t> ComputeDiscretizedWaypoints(const std::vector<waypoint6_t>& wps, const std::vector<spline::Bern<double> >& bernstein, int numSteps);

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
BEZIER_COM_TRAJ_DLLAPI  std::pair<MatrixXX, VectorX> compute6dControlPointInequalities(const ContactData& cData, const std::vector<waypoint6_t>& wps, const std::vector<waypoint6_t>& wpL, const bool useAngMomentum, bool& fail);


/**
 * @brief solve x' h x + 2 g' x, subject to A*x <= b using quadprog
 * @param A Inequality matrix
 * @param b Inequality vector
 * @param H Cost matrix
 * @param g cost Vector
 * @return
 */
BEZIER_COM_TRAJ_DLLAPI ResultData solve(Cref_matrixXX A, Cref_vectorX b, Cref_matrixXX H, Cref_vectorX g, Cref_vectorX initGuess);


/**
 * @brief solve x' h x + 2 g' x, subject to A*x <= b using quadprog, with x of fixed dimension 3
 * @param Ab Inequality matrix and vector
 * @param Hg Cost matrix and vector
 * @return
 */
BEZIER_COM_TRAJ_DLLAPI ResultData solve(const std::pair<MatrixXX, VectorX>& Ab,const std::pair<MatrixXX, VectorX>& Hg,  const Vector3& init);

/**
 * @brief Compute the Bernstein polynoms for a given degree
 * @param degree required degree
 * @return
 */
std::vector<spline::Bern<double> > ComputeBersteinPolynoms(const unsigned int degree);

} // end namespace bezier_com_traj

#endif
