/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef BEZIER_COM_TRAJ_LIB_SOLVE_H
#define BEZIER_COM_TRAJ_LIB_SOLVE_H

#include <Eigen/Dense>
#include <bezier-com-traj/config.hh>
#include <bezier-com-traj/data.hh>

namespace bezier_com_traj
{
      /**
      * @brief solve0step Tries to solve the 0-step capturability problem. Given the current contact phase,
      * a COM position, and an initial velocity, tries to compute a feasible COM trajectory that
      * stops the character without falling.
      * @param pData problem Data. Should contain only one contact phase.
      * @param Ts timelength of each contact phase. Should only contain one value
      * @param timeStep time that the solver has to stop.
      * @return ResultData a struct containing the resulting trajectory, if success is true.
      */
     BEZIER_COM_TRAJ_DLLAPI ResultDataCOMTraj solve0step(const ProblemData& pData, const std::vector<double>& Ts, const double timeStep = -1);

     /**
     * @brief solveEndEffector Tries to produce a trajectory represented as a bezier curve
     * that satisfy position, velocity and acceleration constraint for the initial and final point
     * and that follow as close as possible the input trajectory
     * @param pData problem Data.
     * @param path the path to follow, the class Path must implement the operator (double t) , t \in [0,1] return a Vector3
     * that give the position on the path for a given time
     * @param T time lenght of the trajectory
     * @param timeStep time that the solver has to stop
     * @return ResultData a struct containing the resulting trajectory, if success is true.
     */
     template<typename Path> BEZIER_COM_TRAJ_DLLAPI ResultDataCOMTraj solveEndEffector(const ProblemData& pData,Path path, const double T, const double timeStep);

} // end namespace bezier_com_traj

#endif
