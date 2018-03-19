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

void computeCostMinAcceleration(const ProblemData& pData,const VectorX& Ts,
                                const int pointsPerPhase, MatrixXX& H, VectorX& g);


typedef void (*costCompute) (const ProblemData&, const VectorX&, const int, MatrixXX&, VectorX&);
typedef std::map<CostFunction,costCompute > T_costCompute;
static const T_costCompute costs = boost::assign::map_list_of
        (ACCELERATION, computeCostMinAcceleration);

/** @brief genCostFunction generate a cost function according to the constraints
 * of the problem, and the flag selected in ProblemData.
 * The cost has the form x' H x + 2 g' x.
 * @param Ts times per phase
 * @param pointsPerPhase TODO replace
 * @param H hessian cost matrix to be filled
 * @param g vector matrix
 */
void genCostFunction(const ProblemData& pData,const VectorX& Ts, const int pointsPerPhase, MatrixXX& H, VectorX& g)
{
    T_costCompute::const_iterator cit = costs.find(pData.costFunction_);
    if(cit != costs.end())
        return cit->second(pData, Ts, pointsPerPhase, H, g);
    else
    {
        std::cout<<"Unknown cost function"<<std::endl;
        throw std::runtime_error("Unknown cost function");
    }
}


} // namespace cost
} // namespace bezier_com_traj

#endif
