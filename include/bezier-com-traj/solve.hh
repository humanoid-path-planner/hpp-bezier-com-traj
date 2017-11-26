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
     BEZIER_COM_TRAJ_DLLAPI ResultData solve0step(const ProblemData& pData, const std::vector<double>& Ts, const double timeStep = -1);
} // end namespace bezier_com_traj

#endif
