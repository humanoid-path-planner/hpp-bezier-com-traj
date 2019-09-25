/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#ifndef BEZIER_COM_TRAJ_c0_dc0_ddc0_H
#define BEZIER_COM_TRAJ_c0_dc0_ddc0_H

#include <hpp/bezier-com-traj/data.hh>

namespace bezier_com_traj {
namespace c0_dc0_ddc0 {

static const ConstraintFlag flag = INIT_POS | INIT_VEL | INIT_ACC;

/// ### EQUATION FOR CONSTRAINts on initial position, velocity and acceleration, and only final position (degree = 4)
/**
 * @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one
 * free waypoint (x)
 * @param pi constant waypoints of the curve, assume pi[0] pi[1] x pi[2] p3
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
inline coefs_t evaluateCurveAtTime(const std::vector<point_t>& pi, double t) {
  coefs_t wp;
  double t2 = t * t;
  double t3 = t2 * t;
  // equation found with sympy
  // (-1.0*pi[0] + 3.0*pi[1] - 3.0*pi[2] + 1.0*x)*t**3 + (3.0*pi[0] - 6.0*pi[1] + 3.0*pi[2])*T2 +
  // (-3.0*pi[0] + 3.0*pi[1])*t + 1.0*pi[0],
  wp.first = t3;
  wp.second =
      t3 * (3 * (pi[1] - pi[2]) - pi[0]) + t2 * (3 * (pi[0] + pi[2]) - 6 * pi[1]) + 3 * t * (pi[1] - pi[0]) + pi[0];
  return wp;
}

inline coefs_t evaluateAccelerationCurveAtTime(const std::vector<point_t>& pi, double T, double t) {
  coefs_t wp;
  double alpha = 1. / (T * T);
  // equation found with sympy
  // 6.0*t*alpha*x + (-6.0*pi[0] + 18.0*pi[1] - 18.0*pi[2])/T2*t + (6.0*pi[0] - 12.0*pi[1] + 6.0*pi[2])/T2
  wp.first = 6.0 * t * alpha;
  wp.second = (18. * (pi[1] - pi[2]) - 6. * pi[0]) * alpha * t + (6. * (pi[0] + pi[2]) - 12.0 * pi[1]) * alpha;
  return wp;
}

inline std::vector<point_t> computeConstantWaypoints(const ProblemData& pData, double T) {
  // equation for constraint on initial position, velocity and acceleration, and only final position (degree =
  // 4)(degree 4, 4 constant waypoint and one free (p3)) first, compute the constant waypoints that only depend on
  // pData :
  double n = 3.;
  std::vector<point_t> pi;
  pi.push_back(pData.c0_);                                                                      // pi[0]
  pi.push_back((pData.dc0_ * T / n) + pData.c0_);                                               // pi[1]
  pi.push_back((pData.ddc0_ * T * T / (n * (n - 1))) + (2. * pData.dc0_ * T / n) + pData.c0_);  // pi[2]
  pi.push_back(point_t::Zero());                                                                // x
  return pi;
}

inline bezier_wp_t::t_point_t computeWwaypoints(const ProblemData& pData, double T) {
  bezier_wp_t::t_point_t wps;
  const int DIM_POINT = 6;
  const int DIM_VAR = 3;
  std::vector<point_t> pi = computeConstantWaypoints(pData, T);
  std::vector<Matrix3> Cpi;
  for (std::size_t i = 0; i < pi.size(); ++i) {
    Cpi.push_back(skew(pi[i]));
  }
  const Vector3 g = pData.contacts_.front().contactPhase_->m_gravity;
  const Matrix3 Cg = skew(g), Id = Matrix3::Identity();
  const double T2 = T * T;
  const double alpha = 1 / (T2);

  // equation of waypoints for curve w found with sympy
  waypoint_t w0 = initwp(DIM_POINT, DIM_VAR);
  w0.second.head<3>() = (6 * pi[0] - 12 * pi[1] + 6 * pi[2]) * alpha;
  w0.second.tail<3>() = 1.0 * (1.0 * Cg * T2 * pi[0] - 12.0 * Cpi[0] * pi[1] + 6.0 * Cpi[0] * pi[2]) * alpha;
  wps.push_back(w0);
  waypoint_t w1 = initwp(DIM_POINT, DIM_VAR);
  w1.first.block<3, 3>(0, 0) = 2.0 * alpha * Id;
  w1.first.block<3, 3>(3, 0) = 2.0 * Cpi[0] * alpha;
  w1.second.head<3>() = 1.0 * (4.0 * pi[0] - 6.0 * pi[1]) * alpha;
  w1.second.tail<3>() = 1.0 * (1.0 * Cg * T2 * pi[1] - 6.0 * Cpi[0] * pi[2] + 6.0 * Cpi[1] * pi[2]) * alpha;
  wps.push_back(w1);
  waypoint_t w2 = initwp(DIM_POINT, DIM_VAR);
  w2.first.block<3, 3>(0, 0) = 4.0 * alpha * Id;
  w2.first.block<3, 3>(3, 0) = 1.0 * (-2.0 * Cpi[0] + 6.0 * Cpi[1]) * alpha;
  w2.second.head<3>() = 1.0 * (2.0 * pi[0] - 6.0 * pi[2]) * alpha;
  w2.second.tail<3>() = 1.0 * (1.0 * Cg * T2 * pi[2] - 6.0 * Cpi[1] * pi[2]) * alpha;
  wps.push_back(w2);
  waypoint_t w3 = initwp(DIM_POINT, DIM_VAR);
  w3.first.block<3, 3>(0, 0) = 6 * alpha * Id;
  w3.first.block<3, 3>(3, 0) = 1.0 * (1.0 * Cg * T2 - 6.0 * Cpi[1] + 12.0 * Cpi[2]) * alpha;
  w3.second.head<3>() = (6 * pi[1] - 12 * pi[2]) * alpha;
  // w3.second.head<3>() = 0;
  wps.push_back(w3);
  return wps;
}

inline coefs_t computeFinalVelocityPoint(const ProblemData& pData, double T) {
  coefs_t v;
  std::vector<point_t> pi = computeConstantWaypoints(pData, T);
  // equation found with sympy
  // 3.0*(-pi[2] + x)/T
  v.first = 3. / T;
  v.second = -3. * pi[2] / T;
  return v;
}

}  // namespace c0_dc0_ddc0
}  // namespace bezier_com_traj

#endif
