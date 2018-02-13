/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef BEZIER_COM_TRAJ_LIB_SOLVE_H
#define BEZIER_COM_TRAJ_LIB_SOLVE_H

#include <bezier-com-traj/config.hh>
#include <bezier-com-traj/data.hh>
#include <Eigen/Dense>

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



     /// Methods for transition test :

     /**
      * @brief solveIntersection Solve the QP problem, that find a point inside the constraints Ab that minimise the cost Hg
      * @param Ab s.t. Ax <= b
      * @param Hg min  x'Hx  + g'x
      * @param init x_init
      * @return ResultData
      */
     BEZIER_COM_TRAJ_DLLAPI ResultData solveIntersection(const std::pair<MatrixXX, VectorX>& Ab,const std::pair<MatrixXX, VectorX>& Hg,  const Vector3& init);

     /**
     * @brief solveOnestep Tries to solve the one step problem :  Given two contact phases, an initial and final com position and velocity,
     *  try to compute the CoM trajectory (as a Bezier curve) that connect them
     * @param pData problem Data. Should contain only two contact phase.
     * @param Ts timelength of each contact phase. Should only contain two value
     * @param timeStep time step used by the discretization
     * @return ResultData a struct containing the resulting trajectory, if success is true.
     */
    BEZIER_COM_TRAJ_DLLAPI ResultDataCOMTraj solveOnestep(const ProblemData& pData, const VectorX& Ts, const Vector3& init_guess,const int pointsPerPhase = 3 ,const double feasability_treshold = 0.);




} // end namespace bezier_com_traj

#endif
