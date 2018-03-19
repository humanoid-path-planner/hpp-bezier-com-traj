/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */


#ifndef BEZIER_COM_TRAJ_c0_dc0_ddc0_c1_H
#define BEZIER_COM_TRAJ_c0_dc0_ddc0_c1_H

#include <bezier-com-traj/data.hh>

namespace bezier_com_traj{
namespace c0_dc0_ddc0_c1{

static const ConstraintFlag flag =INIT_POS | INIT_VEL | INIT_ACC | END_POS;

/// ### EQUATION FOR CONSTRAINts on initial position, velocity and acceleration, and only final position (degree = 4)
/**
 * @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve, assume p0 p1 x p2 p3
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
coefs_t evaluateCurveAtTime(std::vector<point_t> pi,double t){
    coefs_t wp;
    double t2 = t*t;
    double t3 = t2*t;
    double t4 = t3*t;
    // equation found with sympy
    wp.first = -4.0*t4 + 4.0*t3;
    wp.second =1.0*pi[0]*t4 - 4.0*pi[0]*t3 + 6.0*pi[0]*t2 - 4.0*pi[0]*t + 1.0*pi[0] - 4.0*pi[1]*t4 + 12.0*pi[1]*t3 - 12.0*pi[1]*t2 + 4.0*pi[1]*t + 6.0*pi[2]*t4 - 12.0*pi[2]*t3 + 6.0*pi[2]*t2 + 1.0*pi[4]*t4;
    //std::cout<<"wp at t = "<<t<<std::endl;
    //std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
    return wp;
}

coefs_t evaluateAccelerationCurveAtTime(std::vector<point_t> pi,double T,double t){
    coefs_t wp;
    double alpha = 1./(T*T);
    double t2 = t*t;
    // equation found with sympy
    wp.first = (-48.0*t2 + 24.0*t)*alpha;
    wp.second = (12.0*pi[0]*t2 - 24.0*pi[0]*t + 12.0*pi[0] - 48.0*pi[1]*t2 + 72.0*pi[1]*t - 24.0*pi[1] + 72.0*pi[2]*t2 - 72.0*pi[2]*t + 12.0*pi[2] + 12.0*pi[4]*t2)*alpha;
    //std::cout<<"acc_wp at t = "<<t<<std::endl;
    //std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
    return wp;
}


std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
    // equation for constraint on initial position, velocity and acceleration, and only final position (degree = 4)(degree 4, 4 constant waypoint and one free (p3))
    // first, compute the constant waypoints that only depend on pData :
    double n = 4.;
    std::vector<point_t> pi;
    pi.push_back(pData.c0_); //p0
    pi.push_back((pData.dc0_ * T / n )+  pData.c0_); // p1
    pi.push_back((pData.ddc0_*T*T/(n*(n-1))) + (2.*pData.dc0_ *T / n) + pData.c0_); // p2
    pi.push_back(point_t::Zero()); // x
    pi.push_back(pData.c1_); // p4
    /*std::cout<<"fixed waypoints : "<<std::endl;
    for(std::vector<point_t>::const_iterator pit = pi.begin() ; pit != pi.end() ; ++pit){
        std::cout<<" pi = "<<*pit<<std::endl;
    }*/
    return pi;
}

std::vector<waypoint6_t> computeWwaypoints(const ProblemData& pData,double T){
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
    w0.second.head<3>() = (12*pi[0] - 24*pi[1] + 12*pi[2])*alpha;
    w0.second.tail<3>() = 1.0*(1.0*Cg*T2*pi[0] - 24.0*Cpi[0]*pi[1] + 12.0*Cpi[0]*pi[2])*alpha;
    wps.push_back(w0);
    waypoint6_t w1 = initwp<waypoint6_t>();
    w1.first.block<3,3>(0,0) = 4.8*alpha*Matrix3::Identity();
    w1.first.block<3,3>(3,0) = 4.8*Cpi[0]*alpha;
    w1.second.head<3>() = 1.0*(7.2*pi[0] - 9.6*pi[1] - 2.4*pi[2])*alpha;
    w1.second.tail<3>() = 1.0*(0.2*Cg*T2*pi[0] + 0.8*Cg*T2*pi[1] - 12.0*Cpi[0]*pi[2] + 9.6*Cpi[1]*pi[2])*alpha;
    wps.push_back(w1);
    waypoint6_t w2 = initwp<waypoint6_t>();
    w2.first.block<3,3>(0,0) = 4.8*alpha*Matrix3::Identity();
    w2.first.block<3,3>(3,0) = 1.0*(-4.8*Cpi[0] + 9.6*Cpi[1])*alpha;
    w2.second.head<3>() = 1.0*(3.6*pi[0] - 9.6*pi[2] + 1.2*pi[4])*alpha;
    w2.second.tail<3>() = 1.0*(0.4*Cg*T2*pi[1] + 0.6*Cg*T2*pi[2] + 1.2*Cpi[0]*pi[4] - 9.6*Cpi[1]*pi[2])*alpha;
    wps.push_back(w2);
    waypoint6_t w3 = initwp<waypoint6_t>();
    w3.first.block<3,3>(3,0) = 1.0*(0.4*Cg*T2 - 9.6*Cpi[1] + 9.6*Cpi[2])*alpha;
    w3.second.head<3>() = 1.0*(1.2*pi[0] + 4.8*pi[1] - 9.6*pi[2] + 3.6*pi[4])*alpha;
    w3.second.tail<3>() = 1.0*(0.6*Cg*T2*pi[2] - 1.2*Cpi[0]*pi[4]  + 4.8*Cpi[1]*pi[4])*alpha;
    wps.push_back(w3);
    waypoint6_t w4 = initwp<waypoint6_t>();
    w4.first.block<3,3>(0,0) = -9.6*alpha*Matrix3::Identity();
    w4.first.block<3,3>(3,0) = 1.0*(0.8*Cg*T2 - 9.6*Cpi[2])*alpha;
    w4.second.head<3>() = 1.0*(4.8*pi[1] - 2.4*pi[2] + 7.2*pi[4])*alpha;
    w4.second.tail<3>() = 1.0*(0.2*Cg*T2*pi[4] - 4.8*Cpi[1]*pi[4] + 12.0*Cpi[2]*pi[4])*alpha;
    wps.push_back(w4);
    waypoint6_t w5 = initwp<waypoint6_t>();
    w5.first.block<3,3>(0,0) = -24*alpha*Matrix3::Identity();
    w5.first.block<3,3>(3,0) = 1.0*(- 24.0*Cpi[4])*alpha;
    w5.second.head<3>() = (12*pi[2] + 12*pi[4])*alpha;
    w5.second.tail<3>() = 1.0*(1.0*Cg*T2*pi[4] - 12.0*Cpi[2]*pi[4])*alpha;
    wps.push_back(w5);
    return wps;
}

coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T){
    coefs_t v;
    // equation found with sympy
    v.first = -4./T;
    v.second = 4.* pData.c1_ / T;
    return v;
}

coefs_t computeFinalAccelerationPoint(const ProblemData& pData,double T){
    coefs_t v;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    // equation found with sympy
    v.first = -24./(T*T);
    v.second = (12.0*pi[2] + 12.*pi[4])/ (T*T);
    return v;
}


}
}

#endif
