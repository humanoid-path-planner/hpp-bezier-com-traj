#ifndef BEZIER_COM_TRAJ_C0_DC0_DDC0_J0_X5_J1_DDC1_DC1_C1_HH
#define BEZIER_COM_TRAJ_C0_DC0_DDC0_J0_X5_J1_DDC1_DC1_C1_HH

/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/data.hh>

namespace bezier_com_traj{
namespace c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1{

static const ConstraintFlag flag = INIT_POS | INIT_VEL | INIT_ACC | END_ACC | END_VEL | END_POS | INIT_JERK | END_JERK | FIVE_FREE_VAR;
static const size_t DIM_VAR = 3*5;
static const size_t DIM_POINT = 3;
/// ### EQUATION FOR CONSTRAINTS ON INIT AND FINAL POSITION AND VELOCITY AND ACCELERATION AND JERK AND 5 variables in the middle (DEGREE = 10)
///
/**
 * @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve, assume pi[8] pi[1] pi[2] pi[3] x0 x1 x2 x3 x4 pi[4] pi[5] pi[6] pi[7]
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
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
    const double t11 = t10*t;
    const double t12 = t11*t;

    // equation found with sympy
    wp.first.block<3,3>(0,0) = Matrix3::Identity()*(+ 495.0*t4 - 3960.0*t5 + 13860.0*t6 - 27720.0*t7 + 34650.0*t8 - 27720.0*t9 + 13860.0*t10 - 3960.0*t11 + 495.0*t12); // x0
    wp.first.block<3,3>(0,3) = Matrix3::Identity()*(+ 792.0*t5 - 5544.0*t6 + 16632.0*t7 - 27720.0*t8 + 27720.0*t9 - 16632.0*t10 + 5544.0*t11 - 792.0*t12); //x1
    wp.first.block<3,3>(0,6) = Matrix3::Identity()*(+ 924.0*t6 - 5544.0*t7 + 13860.0*t8 - 18480.0*t9 + 13860.0*t10 - 5544.0*t11 + 924.0*t12); // x2
    wp.first.block<3,3>(0,9) = Matrix3::Identity()*(+ 792.0*t7 - 3960.0*t8 + 7920.0*t9 - 7920.0*t10 + 3960.0*t11 - 792.0*t12); // x3
    wp.first.block<3,3>(0,12) = Matrix3::Identity()*(+ 495.0*t8 - 1980.0*t9 + 2970.0*t10  - 1980.0*t11 + 495.0*t12); // x4


    wp.second =  1.0*pi[0] +
             t*( -12.0*pi[0] + 12.0*pi[1]) +
             t2*( 66.0*pi[0] - 132.0*pi[1] + 66.0*pi[2]) +
             t3*( -220.0*pi[0] + 660.0*pi[1] - 660.0*pi[2] + 220.0*pi[3]) +
             t4*( 495.0*pi[0] - 1980.0*pi[1] + 2970.0*pi[2] - 1980.0*pi[3]) +
             t5*( -792.0*pi[0] + 3960.0*pi[1] - 7920.0*pi[2] + 7920.0*pi[3] ) +
             t6*( 924.0*pi[0] - 5544.0*pi[1] + 13860.0*pi[2] - 18480.0*pi[3]) +
             t7*( -792.0*pi[0] + 5544.0*pi[1] - 16632.0*pi[2] + 27720.0*pi[3]) +
             t8*( 495.0*pi[0] - 3960.0*pi[1] + 13860.0*pi[2] - 27720.0*pi[3]) +
             t9*( -220.0*pi[0] + 1980.0*pi[1] - 7920.0*pi[2] + 18480.0*pi[3] + 220.0*pi[9]) +
             t10*( 66.0*pi[0] + 66.0*pi[10] - 660.0*pi[1] + 2970.0*pi[2] - 7920.0*pi[3] - 660.0*pi[9]) +
             t11*( -12.0*pi[0] - 132.0*pi[10] + 12.0*pi[11] + 132.0*pi[1] - 660.0*pi[2] + 1980.0*pi[3] + 660.0*pi[9]) +
             t12*( 1.0*pi[0] + 66.0*pi[10] - 12.0*pi[11] + 1.0*pi[12] - 12.0*pi[1] + 66.0*pi[2] - 220.0*pi[3] - 220.0*pi[9]);
    return wp;
}

//TODO
inline waypoint_t evaluateVelocityCurveWaypointAtTime(const std::vector<point_t>& pi,double T,double t){
    waypoint_t wp = initwp(DIM_POINT,DIM_VAR);
    std::cout<<"NOT IMPLEMENTED YET"<<std::endl;

    const double alpha = 1./(T);
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    const double t6 = t5*t;
    const double t7 = t6*t;
    const double t8 = t7*t;
    const double t9 = t8*t;
    const double t10 = t9*t;
    const double t11 = t10*t;

    // equation found with sympy
    wp.first.block<3,3>(0,0) = Matrix3::Identity()*alpha; // x0
    wp.first.block<3,3>(0,3) = Matrix3::Identity()*alpha; //x1
    wp.first.block<3,3>(0,6) = Matrix3::Identity()*alpha; // x2
    wp.first.block<3,3>(0,9) = Matrix3::Identity()*alpha; // x3
    wp.first.block<3,3>(0,12) = Matrix3::Identity()*alpha; // x4

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
    std::cout<<"NOT IMPLEMENTED YET"<<std::endl;

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
    std::cout<<"NOT IMPLEMENTED YET"<<std::endl;

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
    double n = 12.;
    std::vector<point_t> pi;
    pi.push_back(pData.c0_);
    pi.push_back((pData.dc0_ * T / n )+  pData.c0_);
    pi.push_back((pData.ddc0_*T*T/(n*(n-1))) + (2*pData.dc0_ *T / n) + pData.c0_); // * T because derivation make a T appear
    pi.push_back((pData.j0_*T*T*T/(n*(n-1)*(n-2)))+ (3*pData.ddc0_*T*T/(n*(n-1))) + (3*pData.dc0_ *T / n) + pData.c0_);
    pi.push_back(point_t::Zero());
    pi.push_back(point_t::Zero());
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
inline bezier_wp_t::t_point_t computeWwaypoints(const ProblemData& pData,double T){
    bezier_wp_t::t_point_t wps;
    const int DIM_POINT = 6;
    const int DIM_VAR = 15;
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
    assert(pi.size() == 13);

    double alpha = 1. / (T);
    waypoint_t w = initwp(DIM_POINT,DIM_VAR);

    // assign w0:
    w.second= alpha*12*(-pi[0] + pi[1]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w1:
    w.second = alpha*12*(-pi[1] + pi[2]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w2:
    w.second = alpha*12*(-pi[2] + pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w3:
    w.first.block<3,3>(0,0) = 12*alpha*Matrix3::Identity(); // x0
    w.second = alpha*12*(-pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w4:
    w.first.block<3,3>(0,0) = -12*alpha*Matrix3::Identity(); // x0
    w.first.block<3,3>(0,3) = 12*alpha*Matrix3::Identity(); // x1
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.first.block<3,3>(0,3) = -12*alpha*Matrix3::Identity(); // x1
    w.first.block<3,3>(0,6) = 12*alpha*Matrix3::Identity(); // x2
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w6:
    w.first.block<3,3>(0,6) = -12*alpha*Matrix3::Identity(); // x2
    w.first.block<3,3>(0,9) = 12*alpha*Matrix3::Identity(); // x3
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w7:
    w.first.block<3,3>(0,9) = -12*alpha*Matrix3::Identity(); // x3
    w.first.block<3,3>(0,12) = 12*alpha*Matrix3::Identity(); // x4
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w8:
    w.first.block<3,3>(0,12) = -12*alpha*Matrix3::Identity(); // x4
    w.second=alpha*12*pi[9];
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w9:
    w.second=alpha*12*(-pi[9] + pi[10]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w10:
    w.second=alpha*12*(-pi[10] + pi[11]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w11:
    w.second=alpha*12*(-pi[11] + pi[12]);
    wps.push_back(w);
    return wps;
}

std::vector<waypoint_t> computeAccelerationWaypoints(const ProblemData& pData,const double T, std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    std::vector<waypoint_t> wps;
    assert(pi.size() == 13);
    double alpha = 1. / (T*T);

    waypoint_t w = initwp(DIM_POINT,DIM_VAR);

    // assign w0:
    w.second= alpha*132*(pi[0] - 2*pi[1] + pi[2]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w1:
    w.second= alpha*132*(pi[1] - 2*pi[2] + pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w2:
    w.first.block<3,3>(0,0) = 132*alpha*Matrix3::Identity();//x0
    w.second = alpha*132*(pi[2]-pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w3:
    w.first.block<3,3>(0,0) = -264*alpha*Matrix3::Identity(); // x0
    w.first.block<3,3>(0,3) = 132*alpha*Matrix3::Identity(); // x1
    w.second = alpha*132*pi[3];
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w4:
    w.first.block<3,3>(0,0) = 132*alpha*Matrix3::Identity(); // x0
    w.first.block<3,3>(0,3) = -264*alpha*Matrix3::Identity(); // x1
    w.first.block<3,3>(0,6) = 132*alpha*Matrix3::Identity(); // x2
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.first.block<3,3>(0,3) = 132*alpha*Matrix3::Identity(); // x1
    w.first.block<3,3>(0,6) = -264*alpha*Matrix3::Identity(); // x2
    w.first.block<3,3>(0,9) = 132*alpha*Matrix3::Identity(); // x3
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w6:
    w.first.block<3,3>(0,6) = 132*alpha*Matrix3::Identity(); // x2
    w.first.block<3,3>(0,9) = -264*alpha*Matrix3::Identity(); // x3
    w.first.block<3,3>(0,12) = 132*alpha*Matrix3::Identity(); // x4
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w7:
    w.first.block<3,3>(0,9) = 132*alpha*Matrix3::Identity(); // x3
    w.first.block<3,3>(0,12) = -264*alpha*Matrix3::Identity(); // x4
    w.second=alpha*132*pi[7];
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w8:
    w.first.block<3,3>(0,12) = 132*alpha*Matrix3::Identity(); // x4
    w.second=alpha*132*(-2*pi[9] + pi[10]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w9:
    w.second=alpha*132*(pi[9] - 2*pi[10] + pi[11]);
    wps.push_back(w);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w10:
    w.second=alpha*132*(pi[12] + pi[10] - 2*pi[11]);
    wps.push_back(w);
    return wps;
}

std::vector<waypoint_t> computeJerkWaypoints(const ProblemData& pData,const double T, std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    std::vector<waypoint_t> wps;
    assert(pi.size() == 13);

    double alpha = 1. / (T*T*T);

    waypoint_t  w = initwp(DIM_POINT,DIM_VAR);

    // assign w0:
    w.second= alpha*1320*(-pi[0] + 3*pi[1] - 3*pi[2] + pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w1:
    w.first.block<3,3>(0,0) = 1320*alpha*Matrix3::Identity();    //x0
    w.second= alpha*1320*(-pi[1] + 3*pi[2] - 3*pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w2:
    w.first.block<3,3>(0,0) = 1320*-3*alpha*Matrix3::Identity();  // x0
    w.first.block<3,3>(0,3) = 1320*alpha*Matrix3::Identity(); //x1
    w.second = alpha*1320*(-pi[2] + 3*pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w3:
    w.first.block<3,3>(0,0) = 1320*3*alpha*Matrix3::Identity();  // x0
    w.first.block<3,3>(0,3) = 1320*-3*alpha*Matrix3::Identity(); //x1
    w.first.block<3,3>(0,6) = 1320*alpha*Matrix3::Identity(); // x2
    w.second = alpha*1320*(-pi[3]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w4:
    w.first.block<3,3>(0,0) = -1320*alpha*Matrix3::Identity();  // x0
    w.first.block<3,3>(0,3) = 1320*3*alpha*Matrix3::Identity(); //x1
    w.first.block<3,3>(0,6) = 1320*-3*alpha*Matrix3::Identity(); // x2
    w.first.block<3,3>(0,9) = 1320*alpha*Matrix3::Identity(); // x3
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.first.block<3,3>(0,3) = -1320*alpha*Matrix3::Identity();  // x1
    w.first.block<3,3>(0,6) = 1320*3*alpha*Matrix3::Identity(); //x2
    w.first.block<3,3>(0,9) = 1320*-3*alpha*Matrix3::Identity(); // x3
    w.first.block<3,3>(0,12) = 1320*alpha*Matrix3::Identity(); // x4
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w4:
    w.first.block<3,3>(0,6) = -1320*alpha*Matrix3::Identity();  // x2
    w.first.block<3,3>(0,9) = 1320*3*alpha*Matrix3::Identity(); //x3
    w.first.block<3,3>(0,12) = 1320*-3*alpha*Matrix3::Identity(); // x4
    w.second = alpha*1320*pi[9];
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.first.block<3,3>(0,9) = -1320*alpha*Matrix3::Identity(); //x3
    w.first.block<3,3>(0,12) = 1320*3*alpha*Matrix3::Identity(); // x4
    w.second=alpha* 1320*(-3*pi[9] + pi[10]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w5:
    w.first.block<3,3>(0,12) = -1320*alpha*Matrix3::Identity(); // x6
    w.second=alpha*1320*(3*pi[9] - 3*pi[10] + pi[11]);
    wps.push_back(w);
    w = initwp(DIM_POINT,DIM_VAR);
    // assign w6:
    w.second=alpha*1320*(pi[12] - pi[9] + 3*pi[10] - 3*pi[11]);
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

//TODO
inline std::pair<MatrixXX,VectorX> computeVelocityCost(const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi = std::vector<bezier_t::point_t>()){
    MatrixXX H = MatrixXX::Zero(DIM_VAR,DIM_VAR);
    VectorX g  = VectorX::Zero(DIM_VAR);
    if(pi.size() == 0)
        pi = computeConstantWaypoints(pData,T);

    g.segment<3>(0) =  ((-17.8646615739593*pi[0] - 4.24835843773412*pi[10] - 1.80981866649436*pi[11] - 0.408668730537654*pi[12] - 13.8947369836412*pi[1] + 2.56878303943036*pi[2] + 16.0548432515434*pi[3] - 6.66893486967885*pi[9]))/(2*T); // x0
    g.segment<3>(3) = ((-7.93984965058761*pi[0] - 7.0641309185535*pi[10] - 4.08668730511085*pi[11] - 1.22600619206665*pi[12] - 12.1432972410894*pi[1] - 7.06413670827152*pi[2] + 3.85315840674136*pi[3] - 7.50872663647158*pi[9] ))/(2*T); // x1
    g.segment<3>(6) = (-3.26934980716442*pi[0] - 8.9907120599184*pi[10] - 7.76470588258269*pi[11] - 3.26934984520124*pi[12] - 7.76470597007188*pi[1] - 8.99071211730055*pi[2] - 4.49535589801788*pi[3] - 4.49535607858364*pi[9]  )/(2*T); // x2
    g.segment<3>(9) = (-1.22600620726636*pi[0] - 7.06413092270385*pi[10] - 12.1432994250704*pi[11] - 7.93984962409094*pi[12] - 4.08668774398579*pi[1] - 7.0641311269266*pi[2] - 7.50872489092664*pi[3] + 3.85316232209763*pi[9] )/(2*T); // x3
    g.segment<3>(12) = (-0.408668732514974*pi[0] + 2.56877487851457*pi[10] - 13.8947368423667*pi[11] - 17.8646616541281*pi[12] - 1.80981880873492*pi[1] - 4.2483587965255*pi[2] - 6.66893350792178*pi[3] + 16.0548429731073*pi[9])/(2*T); // x4

    H.block<3,3>(0,0) = Matrix3::Identity() * 9.63290527229048  / (T); // x0^2
    H.block<3,3>(3,3) = Matrix3::Identity() * 8.29911962311903  / (T); // x1^2
    H.block<3,3>(6,6) = Matrix3::Identity() * 7.92188615942945  / (T); // x2^2
    H.block<3,3>(9,9) = Matrix3::Identity() * 8.29911871865983  / (T); // x3^2
    H.block<3,3>(12,12) = Matrix3::Identity() *9.63290582796267   / (T); // x4^2

    H.block<3,3>(0,3) = Matrix3::Identity() * 13.4860690009623  / (2*T); // x0*x1 /2
    H.block<3,3>(3,0) = Matrix3::Identity() * 13.4860690009623  / (2*T); // x0*x1 /2
    H.block<3,3>(0,6) = Matrix3::Identity() * 4.14955180440231  / (2*T); // x0*x2 /2
    H.block<3,3>(6,0) = Matrix3::Identity() * 4.14955180440231  / (2*T); // x0*x2 /2
    H.block<3,3>(0,9) = Matrix3::Identity() * - 3.55676093144659  / (2*T); // x0*x3 /2
    H.block<3,3>(9,0) = Matrix3::Identity() * - 3.55676093144659  / (2*T); // x0*x3 /2
    H.block<3,3>(0,12) = Matrix3::Identity() *  - 7.07311260219052 / (2*T); // x0*x4 /2
    H.block<3,3>(12,0) = Matrix3::Identity() *  - 7.07311260219052 / (2*T); // x0*x4 /2

    H.block<3,3>(3,6) = Matrix3::Identity() * 12.4486856197374  / (2*T); // x1*x2 /2
    H.block<3,3>(6,3) = Matrix3::Identity() * 12.4486856197374  / (2*T); // x1*x2 /2
    H.block<3,3>(3,9) = Matrix3::Identity() * 4.20345048607838  / (2*T); // x1*x3 /2
    H.block<3,3>(9,3) = Matrix3::Identity() * 4.20345048607838  / (2*T); // x1*x3 /2
    H.block<3,3>(3,12) = Matrix3::Identity() * - 3.55676456195318  / (2*T); // x1*x4 /2
    H.block<3,3>(12,3) = Matrix3::Identity() * - 3.55676456195318  / (2*T); // x1*x4 /2

    H.block<3,3>(6,9) = Matrix3::Identity() * 12.448679688301  / (2*T); // x2*x3 /2
    H.block<3,3>(9,6) = Matrix3::Identity() * 12.448679688301  / (2*T); // x2*x3 /2
    H.block<3,3>(6,12) = Matrix3::Identity() * 4.149559164651  / (2*T); // x2*x4 /2
    H.block<3,3>(12,6) = Matrix3::Identity() * 4.149559164651  / (2*T); // x2*x4 /2

    H.block<3,3>(9,12) = Matrix3::Identity() * 13.4860680294621  / (2*T); // x3*x4 /2
    H.block<3,3>(12,9) = Matrix3::Identity() * 13.4860680294621  / (2*T); // x3*x4 /2


    double norm=H.norm();
    H /= norm;
    g /= norm;


    return std::make_pair(H,g);
}



}

}


#endif // WAYPOINTS_C0_DC0_DDC0_J0_J1_DDC1_DC1_C1_HH
