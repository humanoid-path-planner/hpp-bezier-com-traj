/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */


#ifndef BEZIER_COM_TRAJ_c0_dc0_ddc0_dc1_c1_H
#define BEZIER_COM_TRAJ_c0_dc0_ddc0_dc1_c1_H

#include <bezier-com-traj/data.hh>

namespace bezier_com_traj{
namespace c0_dc0_ddc0_dc1_c1{

static const ConstraintFlag flag = INIT_POS | INIT_VEL | INIT_ACC | END_VEL | END_POS;

/// ### EQUATION FOR CONSTRAINTS ON INIT AND FINAL POSITION AND VELOCITY AND INIT ACCELERATION (DEGREE = 5)
///
/**
 * @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve, assume p0 p1 x p2 p3
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
inline coefs_t evaluateCurveAtTime(const std::vector<point_t>& pi, double t){
    coefs_t wp;
    double t2 = t*t;
    double t3 = t2*t;
    double t4 = t3*t;
    double t5 = t4*t;
    // equation found with sympy
    wp.first = 10.0*t5 - 20.0*t4 + 10.0*t3;
    wp.second = -1.0*pi[0]*t5 + 5.0*pi[0]*t4 - 10.0*pi[0]*t3 + 10.0*pi[0]*t2 - 5.0*pi[0]*t + 1.0*pi[0] + 5.0*pi[1]*t5 - 20.0*pi[1]*t4 + 30.0*pi[1]*t3 - 20.0*pi[1]*t2 + 5.0*pi[1]*t - 10.0*pi[2]*t5 + 30.0*pi[2]*t4 - 30.0*pi[2]*t3 + 10.0*pi[2]*t2 - 5.0*pi[4]*t5 + 5.0*pi[4]*t4 + 1.0*pi[5]*t5;
    return wp;
}

inline coefs_t evaluateAccelerationCurveAtTime(const std::vector<point_t>& pi,double T,double t){
    coefs_t wp;
    double alpha = 1./(T*T);
    double t2 = t*t;
    double t3 = t2*t;
    // equation found with sympy
    wp.first = (200.0*t3 - 240.0*t2 + 60.0*t)*alpha;
    wp.second = 1.0*(-20.0*pi[0]*t3 + 60.0*pi[0]*t2 - 60.0*pi[0]*t + 20.0*pi[0] + 100.0*pi[1]*t3 - 240.0*pi[1]*t2 + 180.0*pi[1]*t - 40.0*pi[1] - 200.0*pi[2]*t3 + 360.0*pi[2]*t2 - 180.0*pi[2]*t + 20.0*pi[2] - 100.0*pi[4]*t3 + 60.0*pi[4]*t2 + 20.0*pi[5]*t3)*alpha;
    return wp;
}


inline std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
    // equation for constraint on initial and final position and velocity and initial acceleration(degree 5, 5 constant waypoint and one free (p3))
    // first, compute the constant waypoints that only depend on pData :
    double n = 5.;
    std::vector<point_t> pi;
    pi.push_back(pData.c0_); //p0
    pi.push_back((pData.dc0_ * T / n )+  pData.c0_); // p1
    pi.push_back((pData.ddc0_*T*T/(n*(n-1))) + (2.*pData.dc0_ *T / n) + pData.c0_); // p2
    pi.push_back(point_t::Zero()); // x
    pi.push_back((-pData.dc1_ * T / n) + pData.c1_); // p4
    pi.push_back(pData.c1_); // p5
    return pi;
}

inline std::vector<waypoint6_t> computeWwaypoints(const ProblemData& pData,double T){
    std::vector<waypoint6_t> wps;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    std::vector<Matrix3> Cpi;
    for(std::size_t i = 0 ; i < pi.size() ; ++i){
        Cpi.push_back(skew(pi[i]));
    }
    const Vector3 g = pData.contacts_.front().contactPhase_->m_gravity;
    const Matrix3  Cg = skew( g);
    const double T2 = T*T;
    const double alpha = 1/(T2);

    // equation of waypoints for curve w found with sympy
    waypoint6_t w0 = initwp<waypoint6_t>();
    w0.second.head<3>() = (20*pi[0] - 40*pi[1] + 20*pi[2])*alpha;
    w0.second.tail<3>() = 1.0*(1.0*Cg*T2*pi[0] - 40.0*Cpi[0]*pi[1] + 20.0*Cpi[0]*pi[2])*alpha;
    wps.push_back(w0);
    waypoint6_t w1 = initwp<waypoint6_t>();
    w1.first.block<3,3>(0,0) = 8.57142857142857*alpha*Matrix3::Identity();
    w1.first.block<3,3>(3,0) = 8.57142857142857*Cpi[0]*alpha;
    w1.second.head<3>() = 1.0*(11.4285714285714*pi[0] - 14.2857142857143*pi[1] - 5.71428571428572*pi[2])*alpha;
    w1.second.tail<3>() = 1.0*(0.285714285714286*Cg*T2*pi[0] + 0.714285714285714*Cg*T2*pi[1] - 20.0*Cpi[0]*pi[2] + 14.2857142857143*Cpi[1]*pi[2])*alpha;
    wps.push_back(w1);
    waypoint6_t w2 = initwp<waypoint6_t>();
    w2.first.block<3,3>(0,0) = 5.71428571428571*alpha*Matrix3::Identity();
    w2.first.block<3,3>(3,0) = 1.0*(-8.57142857142857*Cpi[0] + 14.2857142857143*Cpi[1])*alpha;
    w2.second.head<3>() = 1.0*(5.71428571428571*pi[0] - 14.2857142857143*pi[2] + 2.85714285714286*pi[4])*alpha;
    w2.second.tail<3>() = 1.0*(0.0476190476190479*Cg*T2*pi[0] + 0.476190476190476*Cg*T2*pi[1] + 0.476190476190476*Cg*T2*pi[2] + 2.85714285714286*Cpi[0]*pi[4] - 14.2857142857143*Cpi[1]*pi[2])*alpha;
    wps.push_back(w2);
    waypoint6_t w3 = initwp<waypoint6_t>();
    w3.first.block<3,3>(0,0) = -2.85714285714286*alpha*Matrix3::Identity();
    w3.first.block<3,3>(3,0) = 1.0*(0.285714285714286*Cg*T2 - 14.2857142857143*Cpi[1] + 11.4285714285714*Cpi[2])*alpha;
    w3.second.head<3>() = 1.0*(2.28571428571429*pi[0] + 5.71428571428571*pi[1] - 11.4285714285714*pi[2] + 5.71428571428571*pi[4] + 0.571428571428571*pi[5])*alpha;
    w3.second.tail<3>() = 1.0*(0.142857142857143*Cg*T2*pi[1] + 0.571428571428571*Cg*T2*pi[2] - 2.85714285714286*Cpi[0]*pi[4] + 0.571428571428571*Cpi[0]*pi[5] + 8.57142857142857*Cpi[1]*pi[4])*alpha;
    wps.push_back(w3);
    waypoint6_t w4 = initwp<waypoint6_t>();
    w4.first.block<3,3>(0,0) = -11.4285714285714*alpha*Matrix3::Identity();
    w4.first.block<3,3>(3,0) = 1.0*(0.571428571428571*Cg*T2 - 11.4285714285714*Cpi[2])*alpha;
    w4.second.head<3>() = 1.0*(0.571428571428571*pi[0] + 5.71428571428571*pi[1] - 2.85714285714286*pi[2] + 5.71428571428571*pi[4] + 2.28571428571429*pi[5])*alpha;
    w4.second.tail<3>() = 1.0*(0.285714285714286*Cg*T2*pi[2] + 0.142857142857143*Cg*T2*pi[4] - 0.571428571428572*Cpi[0]*pi[5] - 8.57142857142857*Cpi[1]*pi[4] + 2.85714285714286*Cpi[1]*pi[5] + 14.2857142857143*Cpi[2]*pi[4])*alpha;
    wps.push_back(w4);
    waypoint6_t w5 = initwp<waypoint6_t>();
    w5.first.block<3,3>(0,0) = -14.2857142857143*alpha*Matrix3::Identity();
    w5.first.block<3,3>(3,0) = 1.0*(0.476190476190476*Cg*T2 - 14.2857142857143*Cpi[4])*alpha;
    w5.second.head<3>() = 1.0*(2.85714285714286*pi[1] + 5.71428571428571*pi[2] + 5.71428571428571*pi[5])*alpha;
    w5.second.tail<3>() = 1.0*(0.476190476190476*Cg*T2*pi[4] + 0.0476190476190476*Cg*T2*pi[5] - 2.85714285714286*Cpi[1]*pi[5] - 14.2857142857143*Cpi[2]*pi[4] + 8.57142857142857*Cpi[2]*pi[5])*alpha;
    wps.push_back(w5);
    waypoint6_t w6 = initwp<waypoint6_t>();
    w6.first.block<3,3>(0,0) = -5.71428571428572*alpha*Matrix3::Identity();
    w6.first.block<3,3>(3,0) = 1.0*(14.2857142857143*Cpi[4] - 20.0*Cpi[5])*alpha;
    w6.second.head<3>() = 1.0*(8.57142857142857*pi[2] - 14.2857142857143*pi[4] + 11.4285714285714*pi[5])*alpha;
    w6.second.tail<3>() = 1.0*(0.714285714285714*Cg*T2*pi[4] + 0.285714285714286*Cg*T2*pi[5] - 8.57142857142858*Cpi[2]*pi[5])*alpha;
    wps.push_back(w6);
    waypoint6_t w7 = initwp<waypoint6_t>();
    w7.first.block<3,3>(0,0) = 20*alpha*Matrix3::Identity();
    w7.first.block<3,3>(3,0) = 1.0*(20.0*Cpi[5])*alpha;
    w7.second.head<3>() = (-40*pi[4] + 20*pi[5])*alpha;
    w7.second.tail<3>() = 1.0*(1.0*Cg*T2*pi[5]  + 40.0*Cpi[4]*pi[5])*alpha;
    wps.push_back(w7);
    return wps;
}

inline coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T){
    coefs_t v;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    // equation found with sympy
    v.first = 0.;
    v.second = (-5.0*pi[4] + 5.0*pi[5])/ T;
    return v;
}

}
}

#endif
