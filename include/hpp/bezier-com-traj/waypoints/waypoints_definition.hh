/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#ifndef BEZIER_COM_TRAJ_WP_DEF_H
#define BEZIER_COM_TRAJ_WP_DEF_H

#include <hpp/bezier-com-traj/data.hh>

namespace bezier_com_traj {
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
coefs_t evaluateCurveAtTime(const ProblemData& pData, const std::vector<point_t>& pi, double t);

/** @brief evaluateVelocityCurveAtTime compute the expression of the point on the curve dc at t,
 * defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
coefs_t evaluateVelocityCurveAtTime(const ProblemData& pData, const std::vector<point_t>& pi, double T, double t);

/** @brief evaluateAccelerationCurveAtTime compute the expression of the point on the curve ddc at t,
 * defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
coefs_t evaluateAccelerationCurveAtTime(const ProblemData& pData, const std::vector<point_t>& pi, double T, double t);

/** @brief evaluateAccelerationCurveAtTime compute the expression of the point on the curve ddc at t,
 * defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
coefs_t evaluateJerkCurveAtTime(const ProblemData& pData, const std::vector<point_t>& pi, double T, double t);

/**
 * @brief computeConstantWaypoints compute the constant waypoints of c(t)
 * defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
std::vector<point_t> computeConstantWaypoints(const ProblemData& pData, double T);

/**
 * @brief computeConstantWaypointsSymbolic compute the constant waypoints of c(t)
 * defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return the waypoints expressed as a polynom of the free waypoint
 */
bezier_wp_t::t_point_t computeConstantWaypointsSymbolic(const ProblemData& pData, double T);

/**
 * @brief computeWwaypoints compute the constant waypoints of dc(t)
 * defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
std::vector<waypoint_t> computeVelocityWaypoints(const ProblemData& pData, const double T,
                                                 std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>());

/**
 * @brief computeWwaypoints compute the constant waypoints of ddc(t)
 * defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
std::vector<waypoint_t> computeAccelerationWaypoints(
    const ProblemData& pData, const double T, std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>());

/**
 * @brief computeWwaypoints compute the constant waypoints of dddc(t)
 * defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
std::vector<waypoint_t> computeJerkWaypoints(const ProblemData& pData, const double T,
                                             std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>());

waypoint_t evaluateCurveWaypointAtTime(const ProblemData& pData, const std::vector<point_t>& pi, double t);
waypoint_t evaluateVelocityCurveWaypointAtTime(const ProblemData& pData, const double T,
                                               const std::vector<point_t>& pi, double t);
waypoint_t evaluateAccelerationCurveWaypointAtTime(const ProblemData& pData, const double T,
                                                   const std::vector<point_t>& pi, double t);
waypoint_t evaluateJerkCurveWaypointAtTime(const ProblemData& pData, const double T, const std::vector<point_t>& pi,
                                           double t);

/**
 * @brief computeConstantWaypoints compute the constant waypoints of w(t)
 * defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
bezier_wp_t::t_point_t computeWwaypoints(const ProblemData& pData, double T);

coefs_t computeFinalVelocityPoint(const ProblemData& pData, double T);

size_t dimVar(const ProblemData& pData);

std::pair<MatrixXX, VectorX> computeVelocityCost(const ProblemData& pData, double T,
                                                 std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>());

}  // namespace bezier_com_traj

#endif
