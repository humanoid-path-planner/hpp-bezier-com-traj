/*
 * Copyright 2018, LAAS-CNRS
 * Author: Steve Tonneau
 */

#ifndef BEZIER_COM_COST_WP_DEF_H
#define BEZIER_COM_COST_WP_DEF_H

#include <bezier-com-traj/data.hh>
#include "boost/assign.hpp"



namespace bezier_com_traj{
/**
* This file contains definitions for the different cost functions used in qp minimization
*/
namespace cost {


/** @brief genCostFunction generate a cost function according to the constraints
 * of the problem, and the flag selected in ProblemData.
 * The cost has the form x' H x + 2 g' x.
 * @param Ts times per phase
 * @param pointsPerPhase TODO replace
 * @param H hessian cost matrix to be filled
 * @param g vector matrix
 */
void genCostFunction(const ProblemData& pData,const VectorX& Ts, const double T,
                     const int pointsPerPhase, MatrixXX& H, VectorX& g);

} // namespace cost
} // namespace bezier_com_traj

#endif
