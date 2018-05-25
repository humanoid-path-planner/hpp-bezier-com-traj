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
inline waypoint_t evaluateCurveWaypointAtTime(const std::vector<point_t>& pi,double t){
    waypoint_t wp = initwp(DIM_POINT,DIM_VAR);
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    const double t6 = t5*t;
    const double t7 = t6*t;
    const double t8 = t7*t;
    const double t9 = t8*t;
    const double t10 = t9*t;

    // equation found with sympy
    wp.first.block<3,3>(0,0) = Matrix3::Identity()*(t4*210.0 -t5*1260.0 +t6*3150.0 -4200.0*t7 + 3150.0*t8 -1260.0*t9 + 210.0*t10); // x0
    wp.first.block<3,3>(0,3) = Matrix3::Identity()*(252.0*t5 - 1260.0*t6 + 2520.0*t7 - 2520.0*t8 + 1260.0*t9 - 252.0*t10); //x1
    wp.first.block<3,3>(0,6) = Matrix3::Identity()*(210.0*t6 - 840.0*t7 + 1260.0*t8 - 840.0*t9 + 210.0*t10); // x2
    wp.second = 1.0*pi[0] +
            t*(-10.0*pi[0] + 10.0*pi[1]) +
            t2*( 45.0*pi[0] - 90.0*pi[1] + 45.0*pi[2]) +
            t3*( -120.0*pi[0] + 360.0*pi[1] - 360.0*pi[2] + 120.0*pi[3]) +
            t4*( 210.0*pi[0] - 840.0*pi[1] + 1260.0*pi[2] - 840.0*pi[3] ) +
            t5*( -252.0*pi[0] + 1260.0*pi[1] - 2520.0*pi[2] + 2520.0*pi[3]) +
            t6*( 210.0*pi[0] - 1260.0*pi[1] + 3150.0*pi[2] - 4200.0*pi[3] ) +
            t7*( -120.0*pi[0] + 840.0*pi[1] - 2520.0*pi[2] + 4200.0*pi[3] + 120.0*pi[7]) +
            t8*( 45.0*pi[0] - 360.0*pi[1] + 1260.0*pi[2] - 2520.0*pi[3] - 360.0*pi[7] + 45.0*pi[8]) +
            t9*( -10.0*pi[0] + 90.0*pi[1] - 360.0*pi[2] + 840.0*pi[3] + 360.0*pi[7] - 90.0*pi[8] + 10.0*pi[9]) +
            t10*( 1.0*pi[0] + 1.0*pi[10] - 10.0*pi[1] + 45.0*pi[2] - 120.0*pi[3] - 120.0*pi[7] + 45.0*pi[8] - 10.0*pi[9]);
    return wp;
}

//TODO
inline waypoint_t evaluateVelocityCurveWaypointAtTime(const std::vector<point_t>& pi,double T,double t){
    waypoint_t wp = initwp(DIM_POINT,DIM_VAR);
    const double alpha = 1./(T);
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    const double t6 = t5*t;
    const double t7 = t6*t;
    const double t8 = t7*t;
    const double t9 = t8*t;
    // equation found with sympy
    wp.first.block<3,3>(0,0) = Matrix3::Identity()*alpha*(1.0*(2100.0*t9 - 11340.0*t8 + 25200.0*t7 - 29400.0*t6 + 18900.0*t5 - 6300.0*t4 + 840.0*t3)); // x0
    wp.first.block<3,3>(0,3) = Matrix3::Identity()*alpha*( 1.0*(-2520.0*t9 + 11340.0*t8 - 20160.0*t7 + 17640.0*t6 - 7560.0*t5 + 1260.0*t4)); //x1
    wp.first.block<3,3>(0,6) = Matrix3::Identity()*alpha*(1.0*(2100.0*t9 - 7560.0*t8 + 10080.0*t7 - 5880.0*t6 + 1260.0*t5)); // x2
    wp.second = (1.0*(-10.0*pi[0] + 10.0*pi[1])
            + t*(1.0*(90.0*pi[0] - 180.0*pi[1] + 90.0*pi[2]))
            + t2*(1.0*(-360.0*pi[0] + 1080.0*pi[1] - 1080.0*pi[2] + 360.0*pi[3]))
            + t3*(1.0*(-90.0*pi[0] + 810.0*pi[1] - 3240.0*pi[2] + 7560.0*pi[3] + 3240.0*pi[7] - 810.0*pi[8] + 90.0*pi[9]))
            + t4*(1.0*(840.0*pi[0] - 3360.0*pi[1] + 5040.0*pi[2] - 3360.0*pi[3]))
            + t5*(1.0*(10.0*pi[0] + 10.0*pi[10] - 100.0*pi[1] + 450.0*pi[2] - 1200.0*pi[3] - 1200.0*pi[7] + 450.0*pi[8] - 100.0*pi[9]))
            + t6*(1.0*(-1260.0*pi[0] + 6300.0*pi[1] - 12600.0*pi[2] + 12600.0*pi[3]))
            + t7*(1.0*(1260.0*pi[0] - 7560.0*pi[1] + 18900.0*pi[2] - 25200.0*pi[3]))
            + t8*(1.0*(-840.0*pi[0] + 5880.0*pi[1] - 17640.0*pi[2] + 29400.0*pi[3] + 840.0*pi[7]))
            + t9*(1.0*(360.0*pi[0] - 2880.0*pi[1] + 10080.0*pi[2] - 20160.0*pi[3] - 2880.0*pi[7] + 360.0*pi[8])))*alpha;
    return wp;
}


//TODO
inline waypoint_t evaluateAccelerationCurveWaypointAtTime(const std::vector<point_t>& pi,double T,double t){
    waypoint_t wp = initwp(DIM_POINT,DIM_VAR);
    const double alpha = 1./(T*T);
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    const double t6 = t5*t;
    const double t7 = t6*t;
    const double t8 = t7*t;
    // equation found with sympy
    wp.first.block<3,3>(0,0) = Matrix3::Identity()*alpha*( 1.0*(18900.0*t8 - 90720.0*t7 + 176400.0*t6 - 176400.0*t5 + 94500.0*t4 - 25200.0*t3 + 2520.0*t2)); // x0
    wp.first.block<3,3>(0,3) = Matrix3::Identity()*alpha*(1.0*(-22680.0*t8 + 90720.0*t7 - 141120.0*t6 + 105840.0*t5 - 37800.0*t4 + 5040.0*t3)); //x1
    wp.first.block<3,3>(0,6) = Matrix3::Identity()*alpha*(1.0*(18900.0*t8 - 60480.0*t7 + 70560.0*t6 - 35280.0*t5 + 6300.0*t4)); // x2
    wp.second = (1.0*(90.0*pi[0] - 180.0*pi[1] + 90.0*pi[2])
            + t*(1.0*(-720.0*pi[0] + 2160.0*pi[1] - 2160.0*pi[2] + 720.0*pi[3]))
            + t2*(1.0*(2520.0*pi[0] - 10080.0*pi[1] + 15120.0*pi[2] - 10080.0*pi[3]))
            + t3*(1.0*(90.0*pi[0] + 90.0*pi[10] - 900.0*pi[1] + 4050.0*pi[2] - 10800.0*pi[3] - 10800.0*pi[7] + 4050.0*pi[8] - 900.0*pi[9]))
            + t4*(1.0*(-5040.0*pi[0] + 25200.0*pi[1] - 50400.0*pi[2] + 50400.0*pi[3]))
            + t5*(1.0*(6300.0*pi[0] - 37800.0*pi[1] + 94500.0*pi[2] - 126000.0*pi[3]))
            + t6*(1.0*(-5040.0*pi[0] + 35280.0*pi[1] - 105840.0*pi[2] + 176400.0*pi[3] + 5040.0*pi[7]))
            + t7*(1.0*(2520.0*pi[0] - 20160.0*pi[1] + 70560.0*pi[2] - 141120.0*pi[3] - 20160.0*pi[7] + 2520.0*pi[8]))
            + t8*(1.0*(-720.0*pi[0] + 6480.0*pi[1] - 25920.0*pi[2] + 60480.0*pi[3] + 25920.0*pi[7] - 6480.0*pi[8] + 720.0*pi[9])))*alpha;
    return wp;
}

//TODO
inline waypoint_t evaluateJerkCurveWaypointAtTime(const std::vector<point_t>& pi,double T,double t){
    waypoint_t wp = initwp(DIM_POINT,DIM_VAR);
    const double alpha = 1./(T*T*T);
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    const double t6 = t5*t;
    const double t7 = t6*t;
    // equation found with sympy
    wp.first.block<3,3>(0,0) = Matrix3::Identity()*alpha*( 1.0*(151200.0*t7 - 635040.0*t6 + 1058400.0*t5 - 882000.0*t4 + 378000.0*t3 - 75600.0*t2 + 5040.0*t)); // x0
    wp.first.block<3,3>(0,3) = Matrix3::Identity()*alpha*(1.0*(-181440.0*t7 + 635040.0*t6 - 846720.0*t5 + 529200.0*t4 - 151200.0*t3 + 15120.0*t2)); //x1
    wp.first.block<3,3>(0,6) = Matrix3::Identity()*alpha*(1.0*(151200.0*t7 - 423360.0*t6 + 423360.0*t5 - 176400.0*t4 + 25200.0*t3)); // x2
    wp.second = (1.0*(-720.0*pi[0] + 2160.0*pi[1] - 2160.0*pi[2] + 720.0*pi[3])
            + t*(1.0*(5040.0*pi[0] - 20160.0*pi[1] + 30240.0*pi[2] - 20160.0*pi[3]))
            + t2*(1.0*(-15120.0*pi[0] + 75600.0*pi[1] - 151200.0*pi[2] + 151200.0*pi[3]))
            + t3*(1.0*(25200.0*pi[0] - 151200.0*pi[1] + 378000.0*pi[2] - 504000.0*pi[3]))
            + t4*(1.0*(-25200.0*pi[0] + 176400.0*pi[1] - 529200.0*pi[2] + 882000.0*pi[3] + 25200.0*pi[7]))
            + t5*(1.0*(15120.0*pi[0] - 120960.0*pi[1] + 423360.0*pi[2] - 846720.0*pi[3] - 120960.0*pi[7] + 15120.0*pi[8]))
            + t6*(1.0*(-5040.0*pi[0] + 45360.0*pi[1] - 181440.0*pi[2] + 423360.0*pi[3] + 181440.0*pi[7] - 45360.0*pi[8] + 5040.0*pi[9]))
            + t7*(1.0*(720.0*pi[0] + 720.0*pi[10] - 7200.0*pi[1] + 32400.0*pi[2] - 86400.0*pi[3] - 86400.0*pi[7] + 32400.0*pi[8] - 7200.0*pi[9])))*alpha;
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


inline std::pair<MatrixXX,VectorX> computeVelocityCost(const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    MatrixXX H = MatrixXX::Zero(9,9);
    VectorX g  = VectorX::Zero(9);
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    g.segment<3>(0) =  (-12.352941184069*pi[0] - 2.03619909502433*pi[10] - 10.5882353430148*pi[1] + 1.2217194516605*pi[2] + 12.2171947000329*pi[3] - 4.66474701697538*pi[7] - 7.21925133730399*pi[8] - 5.42986425333795*pi[9])/(2*T); // x0
    g.segment<3>(3) = (-5.29411764601331*pi[0] - 5.29411764705762*pi[10] - 8.95927605247282*pi[1] - 6.10859723220821*pi[2] + 2.2213080007358*pi[3] + 2.22130810120924*pi[7] - 6.10859728485633*pi[8] - 8.95927601808432*pi[9] )/(2*T); // x1
    g.segment<3>(6) = (-2.03619909557297*pi[0] - 12.3529411764706*pi[10] - 5.42986425052241*pi[1] - 7.21925133714926*pi[2] - 4.66474700749421*pi[3] + 12.2171945706055*pi[7] + 1.22171945695766*pi[8] - 10.5882352941172*pi[9] )/(2*T); // x2

    H.block<3,3>(0,0) = Matrix3::Identity() *  7.77457833646806 / (T); // x0^2
    H.block<3,3>(3,3) = Matrix3::Identity() *  7.25627312788583 / (T); // x1^2
    H.block<3,3>(6,6) = Matrix3::Identity() *  7.77457836216558 / (T); // x2^2
    H.block<3,3>(0,3) = Matrix3::Identity() *   10.8844097406652/ (2*T); // x0*x1 / 2
    H.block<3,3>(3,0) = Matrix3::Identity() *   10.8844097406652/ (2*T); // x0*x1 / 2
    H.block<3,3>(0,6) = Matrix3::Identity() *   2.41875768460934/ (2*T); // x0*x2 / 2
    H.block<3,3>(6,0) = Matrix3::Identity() *   2.41875768460934/ (2*T); // x0*x2 / 2
    H.block<3,3>(3,6) = Matrix3::Identity() *   10.8844097036619/ (2*T); // x1*x2 / 2
    H.block<3,3>(6,3) = Matrix3::Identity() *   10.8844097036619/ (2*T); // x1*x2 / 2

    double norm=H.norm();
    H /= norm;
    g /= norm;


    return std::make_pair(H,g);
}



}

}


#endif // WAYPOINTS_C0_DC0_DDC0_J0_J1_DDC1_DC1_C1_HH
