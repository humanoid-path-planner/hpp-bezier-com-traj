/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#ifndef BEZIER_COM_TRAJ_WP_DEF_H
#define BEZIER_COM_TRAJ_WP_DEF_H

#include <bezier-com-traj/data.hh>



namespace bezier_com_traj{
/**
* This file is used to choose the correct expressions of the curves waypoints,
* depending on the options set in ProblemData.constraints
*/

/** @brief evaluateCurveAtTime compute the expression of the point on the curve c at t,
 * defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
coefs_t evaluateCurveAtTime(const ProblemData& pData, const std::vector<point_t>& pi,double t);

/** @brief evaluateAccelerationCurveAtTime compute the expression of the point on the curve ddc at t,
* defined by the waypoint pi and one free waypoint (x)
* @param pi constant waypoints of the curve
* @param t param (normalized !)
* @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
*/
coefs_t evaluateAccelerationCurveAtTime(const ProblemData& pData, const std::vector<point_t>& pi,double T,double t);

/**
 * @brief computeConstantWaypoints compute the constant waypoints of c(t)
 * defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T);

/**
 * @brief computeWwaypoints compute the constant waypoints of w(t)
 * defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
std::vector<waypoint6_t> computeWwaypoints(const ProblemData& pData,double T);

coefs_t computeFinalAccelerationPoint(const ProblemData& pData,double T);

coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T);

}

#endif
