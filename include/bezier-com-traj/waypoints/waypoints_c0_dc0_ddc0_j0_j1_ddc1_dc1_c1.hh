#ifndef BEZIER_COM_TRAJ_C0_DC0_DDC0_J0_J1_DDC1_DC1_C1_HH
#define BEZIER_COM_TRAJ_C0_DC0_DDC0_J0_J1_DDC1_DC1_C1_HH

/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/data.hh>

namespace bezier_com_traj{
namespace c0_dc0_ddc0_j0_j1_ddc1_dc1_c1{

static const ConstraintFlag flag = INIT_POS | INIT_VEL | INIT_ACC | END_ACC | END_VEL | END_POS | INIT_JERK | END_JERK;
static const size_t DIM_VAR = 3;
static const size_t DIM_POINT = 3;
/// ### EQUATION FOR CONSTRAINTS ON INIT AND FINAL POSITION AND VELOCITY AND ACCELERATION AND JERK (DEGREE = 8)
///
/**
 * @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve, assume pi[8] pi[1] pi[2] pi[3] x pi[4] pi[5] pi[6] pi[7]
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
inline coefs_t evaluateCurveAtTime(const std::vector<point_t>& pi,double t){
    coefs_t wp;
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    const double t6 = t5*t;
    const double t7 = t6*t;
    const double t8 = t7*t;
    // equation found with sympy
    wp.first = 70.0*t8 - 280.0*t7 + 420.0*t6 - 280.0*t5 + 70.0*t4;
    wp.second = 1.0*pi[8]*t8 - 8.0*pi[8]*t7 + 28.0*pi[8]*t6 - 56.0*pi[8]*t5 + 70.0*pi[8]*t4 - 56.0*pi[8]*t3 + 28.0*pi[8]*t2 - 8.0*pi[8]*t + 1.0*pi[8] - 8.0*pi[1]*t8 + 56.0*pi[1]*t7 - 168.0*pi[1]*t6 + 280.0*pi[1]*t5 - 280.0*pi[1]*t4 + 168.0*pi[1]*t3 - 56.0*pi[1]*t2 + 8.0*pi[1]*t + 28.0*pi[2]*t8 - 168.0*pi[2]*t7 + 420.0*pi[2]*t6 - 560.0*pi[2]*t5 + 420.0*pi[2]*t4 - 168.0*pi[2]*t3 + 28.0*pi[2]*t2 - 56.0*pi[3]*pow(t,8 )+ 280.0*pi[3]*t7 - 560.0*pi[3]*t6 + 560.0*pi[3]*t5 - 280.0*pi[3]*t4 + 56.0*pi[3]*t3 - 56.0*pi[5]*t8 + 168.0*pi[5]*t7 - 168.0*pi[5]*t6 + 56.0*pi[5]*pow(t,5 )+ 28.0*pi[6]*t8 - 56.0*pi[6]*t7 + 28.0*pi[6]*t6 - 8.0*pi[7]*t8 + 8.0*pi[7]*t7 + 1.0*pi[8]*t8;
    return wp;
}

inline coefs_t evaluateVelocityCurveAtTime(const std::vector<point_t>& pi,double T,double t){
    coefs_t wp;
    const double alpha = 1./(T);
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    const double t6 = t5*t;
    const double t7 = t6*t;
    // equation found with sympy
    wp.first = (560.0*t7 - 1960.0*t6 + 2520.0*t5 - 1400.0*t4 + 280.0*t3)*alpha;
    wp.second = (8.0*pi[8]*t7 - 56.0*pi[8]*t6 + 168.0*pi[8]*t5 - 280.0*pi[8]*t4 + 280.0*pi[8]*t3 - 168.0*pi[8]*t2 + 56.0*pi[8]*t - 8.0*pi[8] - 64.0*pi[1]*t7 + 392.0*pi[1]*t6 - 1008.0*pi[1]*t5 + 1400.0*pi[1]*t4 - 1120.0*pi[1]*t3 + 504.0*pi[1]*t2 - 112.0*pi[1]*t + 8.0*pi[1] + 224.0*pi[2]*t7 - 1176.0*pi[2]*t6 + 2520.0*pi[2]*t5 - 2800.0*pi[2]*t4 + 1680.0*pi[2]*t3 - 504.0*pi[2]*t2 + 56.0*pi[2]*t - 448.0*pi[3]*t7 + 1960.0*pi[3]*t6 - 3360.0*pi[3]*t5 + 2800.0*pi[3]*t4 - 1120.0*pi[3]*t3 + 168.0*pi[3]*t2 - 448.0*pi[5]*t7 + 1176.0*pi[5]*t6 - 1008.0*pi[5]*t5 + 280.0*pi[5]*t4 + 224.0*pi[6]*t7 - 392.0*pi[6]*t6 + 168.0*pi[6]*t5 - 64.0*pi[7]*t7 + 56.0*pi[7]*t6 + 8.0*pi[8]*t7)*alpha;
    return wp;
}



inline coefs_t evaluateAccelerationCurveAtTime(const std::vector<point_t>& pi,double T,double t){
    coefs_t wp;
    const double alpha = 1./(T*T);
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    const double t6 = t5*t;
    // equation found with sympy
    wp.first = ((3920.0*t6 - 11760.0*t5 + 12600.0*t4 - 5600.0*t3 + 840.0*t2))*alpha;
    wp.second = (56.0*pi[8]*t6 - 336.0*pi[8]*t5 + 840.0*pi[8]*t4 - 1120.0*pi[8]*t3 + 840.0*pi[8]*t2 - 336.0*pi[8]*t + 56.0*pi[8] - 448.0*pi[1]*t6 + 2352.0*pi[1]*t5 - 5040.0*pi[1]*t4 + 5600.0*pi[1]*t3 - 3360.0*pi[1]*t2 + 1008.0*pi[1]*t - 112.0*pi[1] + 1568.0*pi[2]*t6 - 7056.0*pi[2]*t5 + 12600.0*pi[2]*t4 - 11200.0*pi[2]*t3 + 5040.0*pi[2]*t2 - 1008.0*pi[2]*t + 56.0*pi[2] - 3136.0*pi[3]*t6 + 11760.0*pi[3]*t5 - 16800.0*pi[3]*t4 + 11200.0*pi[3]*t3 - 3360.0*pi[3]*t2+ 336.0*pi[3]*t - 3136.0*pi[5]*t6 + 7056.0*pi[5]*t5 - 5040.0*pi[5]*t4 + 1120.0*pi[5]*t3 + 1568.0*pi[6]*t6 - 2352.0*pi[6]*t5 + 840.0*pi[6]*t4 - 448.0*pi[7]*t6 + 336.0*pi[7]*t5 + 56.0*pi[8]*t6)*alpha;
    return wp;
}

inline coefs_t evaluateJerkCurveAtTime(const std::vector<point_t>& pi,double T,double t){
    coefs_t wp;
    const double alpha = 1./(T*T*T);
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    // equation found with sympy
    wp.first = (23520.0*t5 - 58800.0*t4 + 50400.0*t3 - 16800.0*t2 + 1680.0*t)*alpha;
    wp.second = 1.0*(336.0*pi[0]*t5 - 1680.0*pi[0]*t4 + 3360.0*pi[0]*t3 - 3360.0*pi[0]*t2 + 1680.0*pi[0]*t - 336.0*pi[0] - 2688.0*pi[1]*t5 + 11760.0*pi[1]*t4 - 20160.0*pi[1]*t3 + 16800.0*pi[1]*t2 - 6720.0*pi[1]*t + 1008.0*pi[1] + 9408.0*pi[2]*t5 - 35280.0*pi[2]*t4 + 50400.0*pi[2]*t3 - 33600.0*pi[2]*t2 + 10080.0*pi[2]*t - 1008.0*pi[2] - 18816.0*pi[3]*t5 + 58800.0*pi[3]*t4 - 67200.0*pi[3]*t3 + 33600.0*pi[3]*t2 - 6720.0*pi[3]*t + 336.0*pi[3] - 18816.0*pi[5]*t5 + 35280.0*pi[5]*t4 - 20160.0*pi[5]*t3 + 3360.0*pi[5]*t2 + 9408.0*pi[6]*t5 - 11760.0*pi[6]*t4 + 3360.0*pi[6]*t3 - 2688.0*pi[7]*t5 + 1680.0*pi[7]*t4 + 336.0*pi[8]*t5)*alpha;
    return wp;
}
inline waypoint_t evaluateCurveWaypointAtTime(const std::vector<point_t>& pi,double t){
    coefs_t coef = evaluateCurveAtTime(pi,t);
    waypoint_t wp;
    wp.first = Matrix3::Identity()*coef.first;
    wp.second = coef.second;
    return wp;

}
inline waypoint_t evaluateVelocityCurveWaypointAtTime(const std::vector<point_t>& pi,double T,double t){
    coefs_t coef = evaluateVelocityCurveAtTime(pi,T,t);
    waypoint_t wp;
    wp.first = Matrix3::Identity()*coef.first;
    wp.second = coef.second;
    return wp;

}
inline waypoint_t evaluateAccelerationCurveWaypointAtTime(const std::vector<point_t>& pi,double T,double t){
    coefs_t coef = evaluateAccelerationCurveAtTime(pi,T,t);
    waypoint_t wp;
    wp.first = Matrix3::Identity()*coef.first;
    wp.second = coef.second;
    return wp;

}

inline waypoint_t evaluateJerkCurveWaypointAtTime(const std::vector<point_t>& pi,double T,double t){
    coefs_t coef = evaluateJerkCurveAtTime(pi,T,t);
    waypoint_t wp;
    wp.first = Matrix3::Identity()*coef.first;
    wp.second = coef.second;
    return wp;
}



inline std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
    // equation for constraint on initial and final position and velocity and initial acceleration(degree 5, 5 constant waypoint and one free (pi[3]))
    // first, compute the constant waypoints that only depend on pData :
    double n = 8.;
    std::vector<point_t> pi;
    pi.push_back(pData.c0_);
    pi.push_back((pData.dc0_ * T / n )+  pData.c0_);
    pi.push_back((pData.ddc0_*T*T/(n*(n-1))) + (2*pData.dc0_ *T / n) + pData.c0_); // * T because derivation make a T appear
    pi.push_back((pData.j0_*T*T*T/(n*(n-1)*(n-2)))+ (3*pData.ddc0_*T*T/(n*(n-1))) + (3*pData.dc0_ *T / n) + pData.c0_);
    pi.push_back(point_t::Zero());
    pi.push_back((-pData.j1_*T*T*T/(n*(n-1)*(n-2))) + (3*pData.ddc1_ *T*T / (n*(n-1))) - (3 * pData.dc1_ *T / n) + pData.c1_ ); // * T ??
    pi.push_back((pData.ddc1_ *T*T / (n*(n-1))) - (2 * pData.dc1_ *T / n) + pData.c1_ ); // * T ??
    pi.push_back((-pData.dc1_ * T / n) + pData.c1_); // * T ?
    pi.push_back(pData.c1_);
    return pi;
}

inline bezier_wp_t::t_point_t computeWwaypoints(const ProblemData& pData,double T){
    bezier_wp_t::t_point_t wps;
    const int DIM_POINT = 6;
    const int DIM_VAR = 3;
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
    waypoint_t w0 = initwp(DIM_POINT,DIM_VAR);
    w0.second.head<3>() = (30*pi[0] - 60*pi[1] + 30*pi[2])*alpha;
    w0.second.tail<3>() = 1.0*(1.0*Cg*T2*pi[0] - 60.0*Cpi[0]*pi[1] + 30.0*Cpi[0]*pi[2])*alpha;
    wps.push_back(w0);
    waypoint_t w1 = initwp(DIM_POINT,DIM_VAR);
    w1.first.block<3,3>(0,0) = 13.3333333333333*alpha*Matrix3::Identity();
    w1.first.block<3,3>(3,0) = 13.3333333333333*Cpi[0]*alpha;
    w1.second.head<3>() = 1.0*(16.6666666666667*pi[0] - 20.0*pi[1] - 10.0*pi[2])*alpha;
    w1.second.tail<3>() = 1.0*(0.333333333333333*Cg*T2*pi[0] + 0.666666666666667*Cg*T2*pi[1] - 30.0*Cpi[0]*pi[2] + 20.0*Cpi[1]*pi[2])*alpha;
    wps.push_back(w1);
    waypoint_t w2 = initwp(DIM_POINT,DIM_VAR);
    w2.first.block<3,3>(0,0) = 6.66666666666667*alpha*Matrix3::Identity();
    w2.first.block<3,3>(3,0) = 1.0*(-13.3333333333333*Cpi[0] + 20.0*Cpi[1])*alpha;
    w2.second.head<3>() = 1.0*(8.33333333333333*pi[0] - 20.0*pi[2] + 5.0*pi[4])*alpha;
    w2.second.tail<3>() = 1.0*(0.0833333333333334*Cg*T2*pi[0] + 0.5*Cg*T2*pi[1] + 0.416666666666667*Cg*T2*pi[2] + 5.0*Cpi[0]*pi[4] - 20.0*Cpi[1]*pi[2])*alpha;
    wps.push_back(w2);
    waypoint_t w3 = initwp(DIM_POINT,DIM_VAR);
    w3.first.block<3,3>(0,0) = -5.71428571428572*alpha*Matrix3::Identity();
    w3.first.block<3,3>(3,0) = 1.0*(0.238095238095238*Cg*T2 - 20.0*Cpi[1] + 14.2857142857143*Cpi[2])*alpha;
    w3.second.head<3>() = 1.0*(3.57142857142857*pi[0] + 7.14285714285714*pi[1] - 14.2857142857143*pi[2] + 7.85714285714286*pi[4] + 1.42857142857143*pi[5])*alpha;
    w3.second.tail<3>() = 1.0*(0.0119047619047619*Cg*T2*pi[0] + 0.214285714285714*Cg*T2*pi[1] + 0.535714285714286*Cg*T2*pi[2] - 5.0*Cpi[0]*pi[4] + 1.42857142857143*Cpi[0]*pi[5] + 12.8571428571429*Cpi[1]*pi[4])*alpha;
    wps.push_back(w3);
    waypoint_t w4 = initwp(DIM_POINT,DIM_VAR);
    w4.first.block<3,3>(0,0) = -14.2857142857143*alpha*Matrix3::Identity();
    w4.first.block<3,3>(3,0) = 1.0*(0.476190476190476*Cg*T2 - 14.2857142857143*Cpi[2])*alpha;
    w4.second.head<3>() = 1.0*(1.19047619047619*pi[0] + 7.14285714285714*pi[1] - 3.57142857142857*pi[2] + 5.0*pi[4] + 4.28571428571429*pi[5] + 0.238095238095238*pi[6])*alpha;
    w4.second.tail<3>() = 1.0*( 0.0476190476190471*Cg*T2*pi[1] + 0.357142857142857*Cg*T2*pi[2] + 0.119047619047619*Cg*T2*pi[4] - 1.42857142857143*Cpi[0]*pi[5] + 0.238095238095238*Cpi[0]*pi[6] - 12.8571428571429*Cpi[1]*pi[4] + 5.71428571428571*Cpi[1]*pi[5] + 17.8571428571429*Cpi[2]*pi[4])*alpha;
    wps.push_back(w4);
    waypoint_t w5 = initwp(DIM_POINT,DIM_VAR);
    w5.first.block<3,3>(0,0) = -14.2857142857143*alpha*Matrix3::Identity();
    w5.first.block<3,3>(3,0) = 1.0*(0.476190476190476*Cg*T2  - 14.2857142857143*Cpi[4])*alpha;
    w5.second.head<3>() = 1.0*(0.238095238095238*pi[0] + 4.28571428571429*pi[1] + 5.0*pi[2] - 3.57142857142857*pi[4] + 7.14285714285714*pi[5] + 1.19047619047619*pi[6])*alpha;
    w5.second.tail<3>() = 1.0*( + 0.11904761904762*Cg*T2*pi[2] + 0.357142857142857*Cg*T2*pi[4] + 0.0476190476190476*Cg*T2*pi[5]  - 0.238095238095238*Cpi[0]*pi[6] - 5.71428571428572*Cpi[1]*pi[5] + 1.42857142857143*Cpi[1]*pi[6] - 17.8571428571429*Cpi[2]*pi[4] + 12.8571428571429*Cpi[2]*pi[5])*alpha;
    wps.push_back(w5);
    waypoint_t w6 = initwp(DIM_POINT,DIM_VAR);
    w6.first.block<3,3>(0,0) = -5.71428571428571*alpha*Matrix3::Identity();
    w6.first.block<3,3>(3,0) = 1.0*(0.238095238095238*Cg*T2 + 14.2857142857143*Cpi[4] - 20.0*Cpi[5])*alpha;
    w6.second.head<3>() = 1.0*(1.42857142857143*pi[1] + 7.85714285714286*pi[2] - 14.2857142857143*pi[4] + 7.14285714285715*pi[5] + 3.57142857142857*pi[6])*alpha;
    w6.second.tail<3>() = 1.0*(0.535714285714286*Cg*T2*pi[4] + 0.214285714285714*Cg*T2*pi[5] + 0.0119047619047619*Cg*T2*pi[6] - 1.42857142857143*Cpi[1]*pi[6]  - 12.8571428571429*Cpi[2]*pi[5] + 5.0*Cpi[2]*pi[6])*alpha;
    wps.push_back(w6);
    waypoint_t w7 = initwp(DIM_POINT,DIM_VAR);
    w7.first.block<3,3>(0,0) = 6.66666666666667*alpha*Matrix3::Identity();
    w7.first.block<3,3>(3,0) = 1.0*( 20.0*Cpi[5] - 13.3333333333333*Cpi[6])*alpha;
    w7.second.head<3>() = 1.0*(5.0*pi[2] - 20.0*pi[4]  + 8.33333333333333*pi[6])*alpha;
    w7.second.tail<3>() = 1.0*( 0.416666666666667*Cg*T2*pi[4] + 0.5*Cg*T2*pi[5] + 0.0833333333333333*Cg*T2*pi[6]  - 5.0*Cpi[2]*pi[6] + 20.0*Cpi[4]*pi[5])*alpha;
    wps.push_back(w7);
    waypoint_t w8 = initwp(DIM_POINT,DIM_VAR);
    w8.first.block<3,3>(0,0) = 13.3333333333333*alpha*Matrix3::Identity();
    w8.first.block<3,3>(3,0) = 1.0*( 13.3333333333333*Cpi[6])*alpha;
    w8.second.head<3>() = 1.0*(-9.99999999999999*pi[4] - 20.0*pi[5] + 16.6666666666667*pi[6])*alpha;
    w8.second.tail<3>() = 1.0*( 0.666666666666667*Cg*T2*pi[5] + 0.333333333333333*Cg*T2*pi[6]  - 20.0*Cpi[4]*pi[5] + 30.0*Cpi[4]*pi[6])*alpha;
    wps.push_back(w8);
    waypoint_t w9 = initwp(DIM_POINT,DIM_VAR);
    w9.second.head<3>() = (30*pi[4] - 60*pi[5] + 30*pi[6])*alpha;
    w9.second.tail<3>() = 1.0*(1.0*Cg*T2*pi[6] - 30.0*Cpi[4]*pi[6] + 60.0*Cpi[5]*pi[6])*alpha;
    wps.push_back(w9);
    return wps;
}

std::vector<waypoint_t> computeVelocityWaypoints(const ProblemData& pData,const double T, std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    std::vector<waypoint_t> wps;
    assert(pi.size() == 9);

    double alpha = 1. / (T);
    waypoint_t w = initwp(DIM_POINT,DIM_VAR);
    // assign w0:
    w.second= alpha*8*(-pi[0]+pi[1]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w1:
    w.second = alpha*8*(-pi[1]+pi[2]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w2:
    w.second = alpha*8*(-pi[2]+pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w3:
    w.first = 8*alpha*Matrix3::Identity();
    w.second = alpha*-8*pi[3];
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w4:
    w.first = -8*alpha*Matrix3::Identity();
    w.second = alpha*8*pi[5];
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.second=alpha*8*(-pi[5]+pi[6]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w6:
    w.second=alpha*8*(-pi[6]+pi[7]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w7:
    w.second=alpha*8*(-pi[7]+pi[8]);
    wps.push_back(w);
    return wps;
}

std::vector<waypoint_t> computeAccelerationWaypoints(const ProblemData& pData,const double T, std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    std::vector<waypoint_t> wps;
    assert(pi.size() == 9);
    double alpha = 1. / (T*T);

    waypoint_t w = initwp(DIM_POINT,DIM_VAR);

    // assign w0:
    w.second= 56*alpha*(pi[0] - 2*pi[1] + pi[2]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w1:
    w.second= 56*alpha*(pi[1] - 2*pi[2] + pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w2:
    w.first = 56*alpha*Matrix3::Identity();
    w.second = (56*pi[2] - 112*pi[3])*alpha;
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w3:
    w.first = -112*alpha*Matrix3::Identity();
    w.second = (56*pi[3] +56*pi[8])*alpha;
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w4:
    w.first = 56*alpha*Matrix3::Identity();
    w.second = (-112*pi[5] + 56*pi[6])*alpha;
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.second=56*alpha*(pi[5]-2*pi[6]+pi[7]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.second=56*alpha*(pi[6]-2*pi[7]+pi[8]);
    wps.push_back(w);
    return wps;
}


std::vector<waypoint_t> computeJerkWaypoints(const ProblemData& pData,const double T, std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    std::vector<waypoint_t> wps;
    assert(pi.size() == 9);

    double alpha = 1. / (T*T*T);

    waypoint_t w = initwp(DIM_POINT,DIM_VAR);

    // assign w0:
    w.second= 336*(-pi[0] +3*pi[1] - 3*pi[2] + pi[3])*alpha;
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w1:
    w.first = 336*alpha*Matrix3::Identity();
    w.second= 336*(-pi[1] + 3*pi[2] - 3*pi[3])*alpha;
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w2:
    w.first = -3*336*alpha*Matrix3::Identity();
    w.second = 336*(-pi[2] + 3*pi[3] + pi[5])*alpha;
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w3:
    w.first = 3*336*alpha*Matrix3::Identity();
    w.second = 336*(-pi[3] - 3*pi[5] + pi[6])*alpha;
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w4:
    w.first = -336*alpha*Matrix3::Identity();
    w.second =  336*(3*pi[5] - 3*pi[6] + pi[7])*alpha;
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.second= 336*(-pi[5] + 3*pi[6] - 3*pi[7] + pi[8])*alpha;
    wps.push_back(w);
    return wps;
}

inline coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T){
    coefs_t v;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    // equation found with sympy
    v.first = 0.;
    v.second = (-6.0*pi[5] + 6.0*pi[6])/ T;
    return v;
}


inline std::pair<MatrixXX,VectorX> computeVelocityCost(const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    MatrixXX H = MatrixXX::Zero(3,3);
    VectorX g  = VectorX::Zero(3);
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    g =  (-7.8321678321748*pi[0] - 7.83216783237586*pi[1] + 9.13752913728184*pi[3] + 9.13752913758454*pi[5]  - 7.83216783216697*pi[7] - 7.83216783216777*pi[8])/(2*T);
    H = Matrix3::Identity() *  6.52680652684107 / (T);

    double norm=H.norm();
    H /= norm;
    g /= norm;


    return std::make_pair(H,g);
}

}

}


#endif // WAYPOINTS_C0_DC0_DDC0_J0_J1_DDC1_DC1_C1_HH
