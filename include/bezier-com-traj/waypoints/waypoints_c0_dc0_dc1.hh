/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#ifndef BEZIER_COM_TRAJ_C0DC0D1_H
#define BEZIER_COM_TRAJ_C0DC0D1_H

#include <bezier-com-traj/data.hh>

namespace bezier_com_traj{
namespace c0_dc0_dc1{

static const ConstraintFlag flag =INIT_POS | INIT_VEL | END_VEL;

/// ### EQUATION FOR CONSTRAINTS ON INIT AND FINAL POSITION AND VELOCITY (DEGREE = 4)
/** @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve, assume p0 p1 x p2 p3
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
inline coefs_t evaluateCurveAtTime(const std::vector<point_t>& pi,double t){
    coefs_t wp;
    double t2 = t*t;
    double t3 = t2*t;
    // equation found with sympy
    wp.first  = -2.0*t3 + 3.0*t2;
    wp.second = -1.0*pi[0]*t3 + 3.0*pi[0]*t2 - 3.0*pi[0]*t + 1.0*pi[0] + 3.0*pi[1]*t3 - 6.0*pi[1]*t2 + 3.0*pi[1]*t;
    return wp;
}

inline coefs_t evaluateAccelerationCurveAtTime(const std::vector<point_t>& pi,double T,double t){
    coefs_t wp;
    double alpha = 1./(T*T);
    // equation found with sympy
    wp.first  = (-12.0*t + 6.0)*alpha;
    wp.second = (-6.0*pi[0]*t + 6.0*pi[0] + 18.0*pi[1]*t - 12.0*pi[1])*alpha;
    return wp;
}


inline std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
    // equation for constraint on initial and final position and velocity (degree 3, 2 constant waypoints and two free (p2 = p3))
    // first, compute the constant waypoints that only depend on pData :
    if(pData.dc1_.norm() != 0.)
        throw std::runtime_error("Capturability not implemented for spped different than 0");
    std::vector<point_t> pi;
    pi.push_back(pData.c0_); //p0
    pi.push_back((pData.dc0_ * T / 3. )+  pData.c0_); // p1
    pi.push_back(point_t::Zero()); // p2 = x
    pi.push_back(point_t::Zero()); // p3 = x
    return pi;
}

inline std::vector<waypoint6_t> computeWwaypoints(const ProblemData& pData,double T){
    std::vector<waypoint6_t> wps;
    std::vector<point_t> pi = c0_dc0_dc1::computeConstantWaypoints(pData,T);
    std::vector<Matrix3> Cpi;
    for(std::size_t i = 0 ; i < pi.size() ; ++i){
        Cpi.push_back(skew(pi[i]));
    }
    const Vector3 g = pData.contacts_.front().contactPhase_->m_gravity;
    const Matrix3  Cg = skew(g);
    const double T2 = T*T;
    const double alpha = 1./(T2);
    // equation of waypoints for curve w found with sympy
    /*waypoint6_t w0 = initwp<waypoint6_t>();
    w0.first.block<3,3>(0,0) = 6*alpha*Matrix3::Identity();
    w0.first.block<3,3>(3,0) = 6.0*Cpi[0]*alpha;
    w0.second.head<3>() = (6*pi[0] - 12*pi[1])*alpha;
    w0.second.tail<3>() = 1.0*(1.0*Cg*T2*pi[0] - 12.0*Cpi[0]*pi[1])*alpha;
    wps.push_back(w0);
    waypoint6_t w1 = initwp<waypoint6_t>();
    w1.first.block<3,3>(0,0) = 2.0*alpha*Matrix3::Identity();
    w1.first.block<3,3>(3,0) = 1.0*(-4.0*Cpi[0] + 6.0*Cpi[1])*alpha;
    w1.second.head<3>() = 1.0*(4.0*pi[0] - 6.0*pi[1])*alpha;
    w1.second.tail<3>() = 1.0*Cg*pi[1];
    wps.push_back(w1);
    waypoint6_t w2 = initwp<waypoint6_t>();
    w2.first.block<3,3>(0,0) = -2.0*alpha*Matrix3::Identity();
    w2.first.block<3,3>(3,0) = 1.0*(1.0*Cg* - 2.0*Cpi[0])*alpha;
    w2.second.head<3>() = 2.0*pi[0]*alpha;
    wps.push_back(w2);
    waypoint6_t w3 = initwp<waypoint6_t>();
    w3.first.block<3,3>(0,0) = -6*alpha*Matrix3::Identity();
    w3.first.block<3,3>(3,0) = 1.0*(1.0*Cg* - 6.0*Cpi[1])*alpha;
    w3.second.head<3>() = 6*pi[1]*alpha;
    wps.push_back(w3);
    return wps;*/

    waypoint6_t w0 = initwp<waypoint6_t>();
    w0.first.block<3,3>(0,0) = 6*alpha*Matrix3::Identity();
    w0.first.block<3,3>(3,0) = 6.0*Cpi[0]*alpha;
    w0.second.head<3>() = (6*pi[0] - 12*pi[1])*alpha;
    w0.second.tail<3>() = (-Cpi[0])*(12.0*pi[1]*alpha + g);
    wps.push_back(w0);
    waypoint6_t w1 = initwp<waypoint6_t>();
    w1.first.block<3,3>(0,0) =  3*alpha*Matrix3::Identity();
    w1.first.block<3,3>(3,0) = skew(1.5 * (3*pi[1] - pi[0]))*alpha;
    w1.second.head<3>() = 1.5 *alpha* (3* pi[0] - 5*pi[1]);
    w1.second.tail<3>() = (3*alpha*pi[0]).cross(-pi[1]) + 0.25 * (Cg * (3*pi[1] + pi[0]));
    wps.push_back(w1);
    waypoint6_t w2 = initwp<waypoint6_t>();
    w2.first.block<3,3>(0,0) = 0*alpha*Matrix3::Identity();
    w2.first.block<3,3>(3,0) = skew(0.5*g - 3*alpha* pi[0] + 3*alpha*pi[1]);
    w2.second.head<3>() = 3*alpha*(pi[0] - pi[1]);
    w2.second.tail<3>() = 0.5 * Cg*pi[1];
    wps.push_back(w2);
    waypoint6_t w3 = initwp<waypoint6_t>();
    w3.first.block<3,3>(0,0) = -3*alpha*Matrix3::Identity();
    w3.first.block<3,3>(3,0) = skew(g - 1.5 *alpha* (pi[1] + pi[0]));
    w3.second.head<3>() = 1.5*alpha * (pi[1] + pi[0]);
    wps.push_back(w3);
    waypoint6_t w4 = initwp<waypoint6_t>();
    w4.first.block<3,3>(0,0) = -6*alpha * Matrix3::Identity();
    w4.first.block<3,3>(3,0) = skew(g - 6*alpha* pi[1]);
    w4.second.head<3>() = 6*pi[1]*alpha;
    wps.push_back(w4);
    return wps;
}

inline coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T){
    coefs_t v;
    std::vector<point_t> pi = c0_dc0_dc1::computeConstantWaypoints(pData,T);
    // equation found with sympy
    v.first = 0.;
    v.second = point3_t::Zero();
    return v;
}


}
}

#endif
