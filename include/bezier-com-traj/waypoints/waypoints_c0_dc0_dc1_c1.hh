/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#ifndef BEZIER_COM_TRAJ_C0DC0D1C1_H
#define BEZIER_COM_TRAJ_C0DC0D1C1_H

#include <bezier-com-traj/data.hh>

namespace bezier_com_traj{
namespace c0_dc0_dc1_c1{

static const ConstraintFlag flag =INIT_POS | INIT_VEL | END_POS | END_VEL;

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
    double t4 = t3*t;
    // equation found with sympy
    wp.first = ( 6.0*t4 - 12.0*t3 + 6.0*t2);
    wp.second = 1.0*pi[0]*t4 - 4.0*pi[0]*t3 + 6.0*pi[0]*t2 - 4.0*pi[0]*t + 1.0*pi[0] - 4.0*pi[1]*t4 + 12.0*pi[1]*t3 - 12.0*pi[1]*t2 + 4.0*pi[1]*t - 4.0*pi[3]*t4 + 4.0*pi[3]*t3 + 1.0*pi[4]*t4;
    // std::cout<<"wp at t = "<<t<<std::endl;
    // std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
    return wp;
}

inline coefs_t evaluateAccelerationCurveAtTime(const std::vector<point_t>& pi,double T,double t){
    coefs_t wp;
    double alpha = 1./(T*T);
    // equation found with sympy
    wp.first = (72.0*t*t - 72.0*t + 12.0)*alpha;
    wp.second = (12.0*pi[0]*t*t - 24.0*pi[0]*t + 12.0*pi[0] - 48.0*pi[1]*t*t + 72.0*pi[1]*t - 24.0*pi[1] - 48.0*pi[3]*t*t + 24.0*pi[3]*t + 12.0*pi[4]*t*t)*alpha;
    // std::cout<<"acc_wp at t = "<<t<<std::endl;
    // std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
    return wp;
}


inline std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
    // equation for constraint on initial and final position and velocity (degree 4, 4 constant waypoint and one free (p2))
    // first, compute the constant waypoints that only depend on pData :
    int n = 4;
    std::vector<point_t> pi;
    pi.push_back(pData.c0_); //p0
    pi.push_back((pData.dc0_ * T / n )+  pData.c0_); // p1
    pi.push_back(point_t::Zero()); // p2 = x
    pi.push_back((-pData.dc1_ * T / n) + pData.c1_); // p3
    pi.push_back(pData.c1_); // p4
    /* for(int i = 0 ; i < pi.size() ; ++i){
        std::cout<<" p"<<i<<" = "<<pi[i].transpose()<<std::endl;
    }*/
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
    w0.first.block<3,3>(0,0) = 12.*alpha*Matrix3::Identity();
    w0.first.block<3,3>(3,0) = 12.*alpha*Cpi[0];
    w0.second.head<3>() = (12.*pi[0] - 24.*pi[1])*alpha;
    w0.second.tail<3>() = 1.0*Cg*pi[0] - (24.0*Cpi[0]*pi[1])*alpha;
    wps.push_back(w0);
    waypoint6_t w1 = initwp<waypoint6_t>();
    w1.first.block<3,3>(0,0) =  -2.4*alpha*Matrix3::Identity();
    w1.first.block<3,3>(3,0) =(-12.0*Cpi[0] + 9.6*Cpi[1])*alpha;
    w1.second.head<3>() = (7.2*pi[0] - 9.6*pi[1] + 4.8*pi[3])*alpha;
    w1.second.tail<3>() = (0.2*Cg*T2*pi[0] + 0.8*Cg*T2*pi[1] + 4.8*Cpi[0]*pi[3])*alpha;
    wps.push_back(w1);
    waypoint6_t w2 = initwp<waypoint6_t>();
    w2.first.block<3,3>(0,0) =  -9.6*alpha*Matrix3::Identity();
    w2.first.block<3,3>(3,0) = (0.6*Cg*T2 - 9.6*Cpi[1])*alpha;
    w2.second.head<3>() = (3.6*pi[0]  + 4.8*pi[3] + 1.2*pi[4])*alpha;
    w2.second.tail<3>() = (0.4*Cg*T2*pi[1] - 4.8*Cpi[0]*pi[3] + 1.2*Cpi[0]*pi[4] + 9.6*Cpi[1]*pi[3])*alpha;
    wps.push_back(w2);
    waypoint6_t w3 = initwp<waypoint6_t>();
    w3.first.block<3,3>(0,0) = -9.6*alpha*Matrix3::Identity();
    w3.first.block<3,3>(3,0) = (0.6*Cg*T2  - 9.6*Cpi[3])*alpha;
    w3.second.head<3>() = (1.2*pi[0] + 4.8*pi[1]  + 3.6*pi[4])*alpha;
    w3.second.tail<3>() = (0.4*Cg*T2*pi[3]  - 1.2*Cpi[0]*pi[4] - 9.6*Cpi[1]*pi[3] + 4.8*Cpi[1]*pi[4])*alpha;
    wps.push_back(w3);
    waypoint6_t w4 = initwp<waypoint6_t>();
    w4.first.block<3,3>(0,0) = -2.4*alpha*Matrix3::Identity();
    w4.first.block<3,3>(3,0) =(9.6*Cpi[3] - 12.0*Cpi[4])*alpha;
    w4.second.head<3>() = (4.8*pi[1] - 9.6*pi[3] + 7.2*pi[4])*alpha;
    w4.second.tail<3>() = (0.8*Cg*T2*pi[3] + 0.2*Cg*T2*pi[4] - 4.8*Cpi[1]*pi[4])*alpha;
    wps.push_back(w4);
    waypoint6_t w5 = initwp<waypoint6_t>();
    w5.first.block<3,3>(0,0) =12*alpha*Matrix3::Identity();
    w5.first.block<3,3>(3,0) =12.0*Cpi[4]*alpha;
    w5.second.head<3>() = (-24*pi[3] + 12*pi[4])*alpha;
    w5.second.tail<3>() =(Cg*T2*pi[4]  + 24.0*Cpi[3]*pi[4])*alpha;
    wps.push_back(w5);
    return wps;
}

inline coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T){
    coefs_t v;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    // equation found with sympy
    v.first = 0.;
    v.second = (-4.0*pi[3] + 4.0*pi[4])/ T;
    return v;
}

inline coefs_t computeFinalAccelerationPoint(const ProblemData& pData,double T){
    coefs_t v;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    // equation found with sympy
    v.first = 12./(T*T);
    v.second = (-24.0*pi[3] + 12.*pi[4])/ (T*T);
    return v;
}


}
}

#endif
