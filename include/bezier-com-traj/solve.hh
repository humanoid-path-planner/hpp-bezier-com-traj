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
    ResultData solve0step(const ProblemData& pData, const std::vector<double> Ts, const double timeStep = -1);
} // end namespace bezier_com_traj

#endif
