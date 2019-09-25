/*
 * Copyright 2018, LAAS-CNRS
 * Author: Steve Tonneau
 */

#ifndef BEZIER_COM_TRAJ_LIB_COMMON_SOLVE_H
#define BEZIER_COM_TRAJ_LIB_COMMON_SOLVE_H

#include <hpp/bezier-com-traj/local_config.hh>
#include <hpp/bezier-com-traj/data.hh>
#include <hpp/bezier-com-traj/waypoints/waypoints_definition.hh>
#include <hpp/bezier-com-traj/solver/solver-abstract.hpp>

#include <Eigen/Dense>

namespace bezier_com_traj {

/**
 * @brief ComputeDiscretizedWaypoints Given the waypoints defining a bezier curve,
 * computes a discretization of the curve
 * @param wps original waypoints
 * @param bernstein berstein polynoms for
 * @param numSteps desired number of wayoints
 * @return a vector of waypoint representing the discretization of the curve
 */
BEZIER_COM_TRAJ_DLLAPI std::vector<waypoint6_t> ComputeDiscretizedWaypoints(
    const std::vector<waypoint6_t>& wps, const std::vector<spline::Bern<double> >& bernstein, int numSteps);

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
BEZIER_COM_TRAJ_DLLAPI std::pair<MatrixXX, VectorX> compute6dControlPointInequalities(
    const ContactData& cData, const std::vector<waypoint6_t>& wps, const std::vector<waypoint6_t>& wpL,
    const bool useAngMomentum, bool& fail);

/**
 * @brief compute6dControlPointEqualities Given linear and angular control waypoints,
 * compute the equality matrices D and d, D [x; Beta]' = d that constrain the desired control point x and contact
 * forces Beta.
 * @param cData data for the current contact phase
 * @param wps waypoints or the linear part of the trajectory
 * @param wpL waypoints or the angular part of the trajectory
 * @param useAngMomentum whether the angular momentum is consider or equal to 0
 * @param fail set to true if problem is found infeasible
 * @return
 */
BEZIER_COM_TRAJ_DLLAPI std::pair<MatrixXX, VectorX> compute6dControlPointEqualities(
    const ContactData& cData, const std::vector<waypoint6_t>& wps, const std::vector<waypoint6_t>& wpL,
    const bool useAngMomentum, bool& fail);

/**
 * @brief solve x' h x + 2 g' x, subject to A*x <= b using quadprog
 * @param A Inequality matrix
 * @param b Inequality vector
 * @param H Cost matrix
 * @param g cost Vector
 * @param x initGuess initial guess
 * @param minBounds lower bounds on x values. Can be of size 0 if all elements of x are unbounded in that direction, or
 * a size equal to x. Unbounded elements should be lesser or equal to solvers::UNBOUNDED_UP;
 * @param maxBounds upper bounds on x values. Can be of size 0 if all elements of x are unbounded in that direction, or
 * a size equal to x Unbounded elements should be higher or lower than solvers::UNBOUNDED_DOWN;
 * @param solver solver used to solve QP or LP. If LGPK is used, Hessian is not considered as an lp is solved
 * @return
 */
BEZIER_COM_TRAJ_DLLAPI ResultData solve(Cref_matrixXX A, Cref_vectorX b, Cref_matrixXX H, Cref_vectorX g,
                                        Cref_vectorX initGuess, Cref_vectorX minBounds, Cref_vectorX maxBounds,
                                        const solvers::SolverType solver = solvers::SOLVER_QUADPROG);
/**
 * @brief solve x' h x + 2 g' x, subject to A*x <= b and D*x = c using quadprog
 * @param A Inequality matrix
 * @param b Inequality vector
 * @param D Equality matrix
 * @param d Equality vector
 * @param H Cost matrix
 * @param g cost Vector
 * @return
 */
BEZIER_COM_TRAJ_DLLAPI ResultData solve(Cref_matrixXX A, Cref_vectorX b, Cref_matrixXX D, Cref_vectorX d,
                                        Cref_matrixXX H, Cref_vectorX g, Cref_vectorX initGuess,
                                        const solvers::SolverType solver = solvers::SOLVER_QUADPROG);

/**
 * @brief solve x' h x + 2 g' x, subject to A*x <= b using quadprog, with x of fixed dimension 3
 * @param Ab Inequality matrix and vector
 * @param Hg Cost matrix and vector
 * @return
 */
BEZIER_COM_TRAJ_DLLAPI ResultData solve(const std::pair<MatrixXX, VectorX>& Ab, const std::pair<MatrixXX, VectorX>& Hg,
                                        const VectorX& init,
                                        const solvers::SolverType solver = solvers::SOLVER_QUADPROG);

/**
 * @brief solve x' h x + 2 g' x, subject to A*x <= b  and D*x = c using quadprog, with x of fixed dimension 3
 * @param Ab Inequality matrix and vector
 * @param Dd Equality matrix and vector
 * @param Hg Cost matrix and vector
 * @param minBounds lower bounds on x values. Can be of size 0 if all elements of x are unbounded in that direction, or
 * a size equal to x. Unbounded elements should be equal to -std::numeric_limits<double>::infinity();
 * @param maxBounds upper bounds on x values. Can be of size 0 if all elements of x are unbounded in that direction, or
 * a size equal to x Unbounded elements should be equal to  std::numeric_limits<double>::infinity();
 * @param solver solver used to solve QP or LP. If LGPK is used, Hessian is not considered as an lp is solved
 * @return
 */
BEZIER_COM_TRAJ_DLLAPI ResultData solve(const std::pair<MatrixXX, VectorX>& Ab, const std::pair<MatrixXX, VectorX>& Dd,
                                        const std::pair<MatrixXX, VectorX>& Hg, Cref_vectorX minBounds,
                                        Cref_vectorX maxBounds, const VectorX& init,
                                        const solvers::SolverType solver = solvers::SOLVER_QUADPROG);

template <typename Point>
BEZIER_COM_TRAJ_DLLAPI std::vector<std::pair<double, Point> > computeDiscretizedWaypoints(const ProblemData& pData,
                                                                                          double T,
                                                                                          const T_time& timeArray);

template <typename Point>
BEZIER_COM_TRAJ_DLLAPI std::vector<std::pair<double, Point> > computeDiscretizedAccelerationWaypoints(
    const ProblemData& pData, double T, const T_time& timeArray);

}  // end namespace bezier_com_traj

#include "common_solve_methods.inl"

#endif
