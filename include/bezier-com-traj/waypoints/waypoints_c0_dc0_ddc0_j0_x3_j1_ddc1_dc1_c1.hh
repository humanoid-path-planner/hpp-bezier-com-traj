#ifndef BEZIER_COM_TRAJ_C0_DC0_DDC0_J0_X3_J1_DDC1_DC1_C1_HH
#define BEZIER_COM_TRAJ_C0_DC0_DDC0_J0_X3_J1_DDC1_DC1_C1_HH

/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/data.hh>

namespace bezier_com_traj{
namespace c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1{

static const ConstraintFlag flag = INIT_POS | INIT_VEL | INIT_ACC | END_ACC | END_VEL | END_POS | INIT_JERK | END_JERK | THREE_FREE_VAR;
static const size_t DIM_VAR = 9;
static const size_t DIM_POINT = 3;
/// ### EQUATION FOR CONSTRAINTS ON INIT AND FINAL POSITION AND VELOCITY AND ACCELERATION AND JERK AND 3 variables in the middle (DEGREE = 10)
///
/**
 * @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve, assume pi[8] pi[1] pi[2] pi[3] x0 x1 x2 pi[4] pi[5] pi[6] pi[7]
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
//TODO
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

//TODO
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


//TODO
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

//TODO
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


inline std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
    // equation for constraint on initial and final position and velocity and initial acceleration(degree 5, 5 constant waypoint and one free (pi[3]))
    // first, compute the constant waypoints that only depend on pData :
    double n = 10.;
    std::vector<point_t> pi;
    pi.push_back(pData.c0_);
    pi.push_back((pData.dc0_ * T / n )+  pData.c0_);
    pi.push_back((pData.ddc0_*T*T/(n*(n-1))) + (2*pData.dc0_ *T / n) + pData.c0_); // * T because derivation make a T appear
    pi.push_back((pData.j0_*T*T*T/(n*(n-1)*(n-2)))+ (3*pData.ddc0_*T*T/(n*(n-1))) + (3*pData.dc0_ *T / n) + pData.c0_);
    pi.push_back(point_t::Zero());
    pi.push_back(point_t::Zero());
    pi.push_back(point_t::Zero());
    pi.push_back((-pData.j1_*T*T*T/(n*(n-1)*(n-2))) + (3*pData.ddc1_ *T*T / (n*(n-1))) - (3 * pData.dc1_ *T / n) + pData.c1_ ); // * T ??
    pi.push_back((pData.ddc1_ *T*T / (n*(n-1))) - (2 * pData.dc1_ *T / n) + pData.c1_ ); // * T ??
    pi.push_back((-pData.dc1_ * T / n) + pData.c1_); // * T ?
    pi.push_back(pData.c1_);
    return pi;
}

//TODO
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
    std::cout<<"NOT IMPLEMENTED YET"<<std::endl;
    return wps;
}

std::vector<waypoint_t> computeVelocityWaypoints(const ProblemData& pData,const double T, std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    std::vector<waypoint_t> wps;
    assert(pi.size() == 11);

    double alpha = 1. / (T);
    waypoint_t w = initwp(DIM_POINT,DIM_VAR);

    // assign w0:
    w.second= alpha*10*(-pi[0] + pi[1]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w1:
    w.second = alpha*10*(-pi[1] + pi[2]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w2:
    w.second = alpha*10*(-pi[2] + pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w3:
    w.first.block<3,3>(0,0) = 10*alpha*Matrix3::Identity(); // x0
    w.second = alpha*10*(-pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w4:
    w.first.block<3,3>(0,0) = -10*alpha*Matrix3::Identity(); // x0
    w.first.block<3,3>(0,3) = 10*alpha*Matrix3::Identity(); // x1
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.first.block<3,3>(0,3) = -10*alpha*Matrix3::Identity(); // x1
    w.first.block<3,3>(0,6) = 10*alpha*Matrix3::Identity(); // x2
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w6:
    w.first.block<3,3>(0,6) = -10*alpha*Matrix3::Identity(); // x2
    w.second=alpha*10*pi[7];
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w7:
    w.second=alpha*10*(-pi[7] + pi[8]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w8:
    w.second=alpha*10*(-pi[8] + pi[9]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w9:
    w.second=alpha*10*(pi[10] - pi[9]);
    wps.push_back(w);
    return wps;
}

std::vector<waypoint_t> computeAccelerationWaypoints(const ProblemData& pData,const double T, std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    std::vector<waypoint_t> wps;
    assert(pi.size() == 11);
    double alpha = 1. / (T*T);

    waypoint_t w = initwp(DIM_POINT,DIM_VAR);

    // assign w0:
    w.second= alpha*90*(pi[0] - 2*pi[1] + pi[2]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w1:
    w.second= alpha*90*(pi[1] - 2*pi[2] + pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w2:
    w.first.block<3,3>(0,0) = 90*alpha*Matrix3::Identity();//x0
    w.second = alpha*90*(pi[2]-pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w3:
    w.first.block<3,3>(0,0) = -180*alpha*Matrix3::Identity(); // x0
    w.first.block<3,3>(0,3) = 90*alpha*Matrix3::Identity(); // x1
    w.second = alpha*90*pi[3];
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w4:
    w.first.block<3,3>(0,0) = 90*alpha*Matrix3::Identity(); // x0
    w.first.block<3,3>(0,3) = -180*alpha*Matrix3::Identity(); // x1
    w.first.block<3,3>(0,6) = 90*alpha*Matrix3::Identity(); // x2
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.first.block<3,3>(0,3) = 90*alpha*Matrix3::Identity(); // x1
    w.first.block<3,3>(0,6) = -180*alpha*Matrix3::Identity(); // x2
    w.second=alpha*90*pi[7];
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w6:
    w.first.block<3,3>(0,6) = 90*alpha*Matrix3::Identity(); // x2
    w.second=alpha*90*(-2*pi[7] + pi[8]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w7:
    w.second=alpha*90*(pi[7] - 2*pi[8] + pi[9]);
    wps.push_back(w);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w8:
    w.second=alpha*90*(pi[10] + pi[8] - 2*pi[9]);
    wps.push_back(w);
    return wps;
}

std::vector<waypoint_t> computeJerkWaypoints(const ProblemData& pData,const double T, std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    std::vector<waypoint_t> wps;
    assert(pi.size() == 11);

    double alpha = 1. / (T*T*T);

    waypoint_t  w = initwp(DIM_POINT,DIM_VAR);

    // assign w0:
    w.second= alpha*720*(-pi[0] + 3*pi[1] - 3*pi[2] + pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w1:
    w.first.block<3,3>(0,0) = 720*alpha*Matrix3::Identity();    //x0
    w.second= alpha*720*(-pi[1] + 3*pi[2] - 3*pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w2:
    w.first.block<3,3>(0,0) = 720*-3*alpha*Matrix3::Identity();  // x0
    w.first.block<3,3>(0,3) = 720*alpha*Matrix3::Identity(); //x1
    w.second = alpha*720*(-pi[2] + 3*pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w3:
    w.first.block<3,3>(0,0) = 720*3*alpha*Matrix3::Identity();  // x0
    w.first.block<3,3>(0,3) = 720*-3*alpha*Matrix3::Identity(); //x1
    w.first.block<3,3>(0,6) = 720*alpha*Matrix3::Identity(); // x2
    w.second = alpha*720*(-pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w4:
    w.first.block<3,3>(0,0) = -720*alpha*Matrix3::Identity();  // x0
    w.first.block<3,3>(0,3) = 720*3*alpha*Matrix3::Identity(); //x1
    w.first.block<3,3>(0,6) = 720*-3*alpha*Matrix3::Identity(); // x2
    w.second = alpha*720*pi[7];
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.first.block<3,3>(0,3) = -720*alpha*Matrix3::Identity(); //x1
    w.first.block<3,3>(0,6) = 720*3*alpha*Matrix3::Identity(); // x2
    w.second=alpha* 720*(-3*pi[7] + pi[8]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.first.block<3,3>(0,6) = -720*alpha*Matrix3::Identity(); // x2
    w.second=alpha*720*(3*pi[7] - 3*pi[8] + pi[9]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w6:
    w.second=alpha*720*(pi[10] - pi[7] + 3*pi[8] - 3*pi[9]);
    wps.push_back(w);
    return wps;
}

//TODO
inline coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T){
    coefs_t v;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    // equation found with sympy
    v.first = 0.;
    v.second = (-6.0*pi[5] + 6.0*pi[6])/ T;
    return v;
}

}

}


#endif // WAYPOINTS_C0_DC0_DDC0_J0_J1_DDC1_DC1_C1_HH
