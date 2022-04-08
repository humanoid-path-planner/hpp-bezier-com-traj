/*
 * Copyright 2018, LAAS-CNRS
 * Author: Steve Tonneau
 */

#ifndef BEZIER_COM_TRAJ_LIB_SOLVE_H
#define BEZIER_COM_TRAJ_LIB_SOLVE_H

#include <Eigen/Dense>
#include <hpp/bezier-com-traj/data.hh>
#include <hpp/bezier-com-traj/local_config.hh>

namespace bezier_com_traj {
/**
 * @brief solve0step Tries to solve the 0-step capturability problem. Given the
 * current contact phase, a COM position, and an initial velocity, tries to
 * compute a feasible COM trajectory that stops the character without falling.
 * In this specific implementation, the considered constraints are:
 * init position and velocity, 0 velocity constraints (acceleration constraints
 * are ignored)
 * @param pData problem Data. Should contain only one contact phase.
 * @param Ts timelength of each contact phase. Should only contain one value
 * @param timeStep time that the solver has to stop.
 * @return ResultData a struct containing the resulting trajectory, if success
 * is true.
 */
BEZIER_COM_TRAJ_DLLAPI ResultDataCOMTraj
solve0step(const ProblemData& pData, const std::vector<double>& Ts,
           const double timeStep = -1);

/**
 * @brief computeCOMTraj Tries to solve the one step problem :  Given two or
 * three contact phases, an initial and final com position and velocity, try to
 * compute the CoM trajectory (as a Bezier curve) that connect them
 * @param pData problem Data.
 * @param Ts timelength of each contact phase. Should be the same legnth as
 * pData.contacts
 * @param timeStep time step used by the discretization
 * @return ResultData a struct containing the resulting trajectory, if success
 * is true.
 */
BEZIER_COM_TRAJ_DLLAPI ResultDataCOMTraj
computeCOMTrajFixedSize(const ProblemData& pData, const VectorX& Ts,
                        const unsigned int pointsPerPhase = 3);

/**
 * @brief computeCOMTraj Tries to solve the one step problem :  Given two or
 * three contact phases, an initial and final com position and velocity, try to
 * compute the CoM trajectory (as a Bezier curve) that connect them
 * @param pData problem Data.
 * @param Ts timelength of each contact phase. Should be the same length as
 * pData.contacts
 * @param timeStep time step used by the discretization, if -1 : use the
 * continuous fomulation
 * @param solver solver used to perform optimization. WARNING: if the continuous
 * force formulation is used, it is highly recommended to use the SOLVER_GLPK
 * solver if available and a quadratic cost is not necessary, as these solvers
 * are increasely more computationnaly efficient for the problem
 * @return ResultData a struct containing the resulting trajectory, if success
 * is true.
 */
BEZIER_COM_TRAJ_DLLAPI ResultDataCOMTraj computeCOMTraj(
    const ProblemData& pData, const VectorX& Ts, const double timeStep = -1,
    const solvers::SolverType solver = solvers::SOLVER_QUADPROG);

}  // end namespace bezier_com_traj

#endif
