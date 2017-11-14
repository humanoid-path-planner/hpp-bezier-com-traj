/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef BEZIER_COM_TRAJ_LIB_SOLVE_H
#define BEZIER_COM_TRAJ_LIB_SOLVE_H

#include <Eigen/Dense>
#include <bezier-com-traj/config.hh>
#include <bezier-com-traj/data.hh>

namespace centroidal_dynamics
{
    ResultData solve0step(const ProblemData& pData, const double T);
} // end namespace centroidal_dynamics

#endif
