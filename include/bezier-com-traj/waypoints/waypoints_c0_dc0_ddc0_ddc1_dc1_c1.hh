/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/data.hh>

namespace bezier_com_traj{
namespace c0_dc0_ddc0_ddc1_dc1_c1{

typedef waypoint3_t waypoint_t;
typedef std::pair<double,point3_t> coefs_t;

static const ConstraintFlag flag = INIT_POS | INIT_VEL | INIT_ACC | END_ACC | END_VEL | END_POS;

/*bool useThisConstraints(Constraints constraints){
    if(constraints.c0_ && constraints.dc0_ && constraints.ddc0_ && constraints.ddc1_ && constraints.dc1_ && constraints.c1_)
        return true;
    else
        return false;
}*/

/// ### EQUATION FOR CONSTRAINTS ON INIT AND FINAL POSITION AND VELOCITY AND ACCELERATION (DEGREE = 6)
///
/**
 * @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve, assume p0 p1 p2 x p3 p4 p5
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
coefs_t evaluateCurveAtTime(std::vector<point_t> pi,double t){
    coefs_t wp;
    double t2 = t*t;
    double t3 = t2*t;
    double t4 = t3*t;
    double t5 = t4*t;
    double t6 = t5*t;
    // equation found with sympy
    wp.first = -20.0*t6 + 60.0*t5 - 60.0*t4 + 20.0*t3;
    wp.second = 1.0*pi[0]*t6 - 6.0*pi[0]*t5 + 15.0*pi[0]*t4 - 20.0*pi[0]*t3 + 15.0*pi[0]*t2 - 6.0*pi[0]*t + 1.0*pi[0] - 6.0*pi[1]*t6 + 30.0*pi[1]*t5 - 60.0*pi[1]*t4 + 60.0*pi[1]*t3 - 30.0*pi[1]*t2 + 6.0*pi[1]*t + 15.0*pi[2]*t6 - 60.0*pi[2]*t5 + 90.0*pi[2]*t4 - 60.0*pi[2]*t3 + 15.0*pi[2]*t2 + 15.0*pi[4]*t6 - 30.0*pi[4]*t5 + 15.0*pi[4]*t4 - 6.0*pi[5]*t6 + 6.0*pi[5]*t5 + 1.0*pi[6]*t6;
    // std::cout<<"wp at t = "<<t<<std::endl;
    // std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
    return wp;
}

coefs_t evaluateAccelerationCurveAtTime(std::vector<point_t> pi,double T,double t){
    coefs_t wp;
    double alpha = 1./(T*T);
    double t2 = t*t;
    double t3 = t2*t;
    double t4 = t3*t;
    // equation found with sympy
    wp.first = 1.0*(-600.0*t4 + 1200.0*t3 - 720.0*t2 + 120.0*t)*alpha;
    wp.second = 1.0*(30.0*pi[0]*t4 - 120.0*pi[0]*t3 + 180.0*pi[0]*t2 - 120.0*pi[0]*t + 30.0*pi[0] - 180.0*pi[1]*t4 + 600.0*pi[1]*t3 - 720.0*pi[1]*t2 + 360.0*pi[1]*t - 60.0*pi[1] + 450.0*pi[2]*t4 - 1200.0*pi[2]*t3 + 1080.0*pi[2]*t2 - 360.0*pi[2]*t + 30.0*pi[2] + 450.0*pi[4]*t4 - 600.0*pi[4]*t3 + 180.0*pi[4]*t2 - 180.0*pi[5]*t4 + 120.0*pi[5]*t3 + 30.0*pi[6]*t4)*alpha;
    // std::cout<<"acc_wp at t = "<<t<<std::endl;
    // std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
    return wp;
}


std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
    // equation for constraint on initial and final position and velocity and initial acceleration(degree 5, 5 constant waypoint and one free (p3))
    // first, compute the constant waypoints that only depend on pData :
    double n = 6.;
    std::vector<point_t> pi;
    pi.push_back(pData.c0_); //p0
    pi.push_back((pData.dc0_ * T / n )+  pData.c0_); // p1
    pi.push_back((pData.ddc0_*T*T/(n*(n-1))) + (2.*pData.dc0_ *T / n) + pData.c0_); // p2
    pi.push_back(point_t::Zero()); // x
    pi.push_back((pData.ddc1_*T*T/(n*(n-1))) - (2*pData.dc1_*T/n) + pData.c1_) ;
    pi.push_back((-pData.dc1_ * T / n) + pData.c1_); // p4
    pi.push_back(pData.c1_); // p5
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
    for(int i = 0 ; i < pi.size() ; ++i){
        Cpi.push_back(skew(pi[i]));
    }
    const Vector3 g = pData.contacts_.front().contactPhase_->m_gravity;
    const Matrix3  Cg = skew( g);
    const double T2 = T*T;
    const double alpha = 1/(T2);

    // equation of waypoints for curve w found with sympy
    waypoint6_t w0 = initwp<waypoint6_t>();
    w0.second.head<3>() = (30*pi[0] - 60*pi[1] + 30*pi[2])*alpha;
    w0.second.tail<3>() = 1.0*(1.0*Cg*T2*pi[0] - 60.0*Cpi[0]*pi[1] + 30.0*Cpi[0]*pi[2])*alpha;
    wps.push_back(w0);
    waypoint6_t w1 = initwp<waypoint6_t>();
    w1.first.block<3,3>(0,0) = 13.3333333333333*alpha*Matrix3::Identity();
    w1.first.block<3,3>(3,0) = 13.3333333333333*Cpi[0]*alpha;
    w1.second.head<3>() = 1.0*(16.6666666666667*pi[0] - 20.0*pi[1] - 10.0*pi[2])*alpha;
    w1.second.tail<3>() = 1.0*(0.333333333333333*Cg*T2*pi[0] + 0.666666666666667*Cg*T2*pi[1] - 30.0*Cpi[0]*pi[2] + 20.0*Cpi[1]*pi[2])*alpha;
    wps.push_back(w1);
    waypoint6_t w2 = initwp<waypoint6_t>();
    w2.first.block<3,3>(0,0) = 6.66666666666667*alpha*Matrix3::Identity();
    w2.first.block<3,3>(3,0) = 1.0*(-13.3333333333333*Cpi[0] + 20.0*Cpi[1])*alpha;
    w2.second.head<3>() = 1.0*(8.33333333333333*pi[0] - 20.0*pi[2] + 5.0*pi[4])*alpha;
    w2.second.tail<3>() = 1.0*(0.0833333333333334*Cg*T2*pi[0] + 0.5*Cg*T2*pi[1] + 0.416666666666667*Cg*T2*pi[2] + 5.0*Cpi[0]*pi[4] - 20.0*Cpi[1]*pi[2])*alpha;
    wps.push_back(w2);
    waypoint6_t w3 = initwp<waypoint6_t>();
    w3.first.block<3,3>(0,0) = -5.71428571428572*alpha*Matrix3::Identity();
    w3.first.block<3,3>(3,0) = 1.0*(0.238095238095238*Cg*T2 - 20.0*Cpi[1] + 14.2857142857143*Cpi[2])*alpha;
    w3.second.head<3>() = 1.0*(3.57142857142857*pi[0] + 7.14285714285714*pi[1] - 14.2857142857143*pi[2] + 7.85714285714286*pi[4] + 1.42857142857143*pi[5])*alpha;
    w3.second.tail<3>() = 1.0*(0.0119047619047619*Cg*T2*pi[0] + 0.214285714285714*Cg*T2*pi[1] + 0.535714285714286*Cg*T2*pi[2] - 5.0*Cpi[0]*pi[4] + 1.42857142857143*Cpi[0]*pi[5] + 12.8571428571429*Cpi[1]*pi[4])*alpha;
    wps.push_back(w3);
    waypoint6_t w4 = initwp<waypoint6_t>();
    w4.first.block<3,3>(0,0) = -14.2857142857143*alpha*Matrix3::Identity();
    w4.first.block<3,3>(3,0) = 1.0*(0.476190476190476*Cg*T2 - 14.2857142857143*Cpi[2])*alpha;
    w4.second.head<3>() = 1.0*(1.19047619047619*pi[0] + 7.14285714285714*pi[1] - 3.57142857142857*pi[2] + 5.0*pi[4] + 4.28571428571429*pi[5] + 0.238095238095238*pi[6])*alpha;
    w4.second.tail<3>() = 1.0*( 0.0476190476190471*Cg*T2*pi[1] + 0.357142857142857*Cg*T2*pi[2] + 0.119047619047619*Cg*T2*pi[4] - 1.42857142857143*Cpi[0]*pi[5] + 0.238095238095238*Cpi[0]*pi[6] - 12.8571428571429*Cpi[1]*pi[4] + 5.71428571428571*Cpi[1]*pi[5] + 17.8571428571429*Cpi[2]*pi[4])*alpha;
    wps.push_back(w4);
    waypoint6_t w5 = initwp<waypoint6_t>();
    w5.first.block<3,3>(0,0) = -14.2857142857143*alpha*Matrix3::Identity();
    w5.first.block<3,3>(3,0) = 1.0*(0.476190476190476*Cg*T2  - 14.2857142857143*Cpi[4])*alpha;
    w5.second.head<3>() = 1.0*(0.238095238095238*pi[0] + 4.28571428571429*pi[1] + 5.0*pi[2] - 3.57142857142857*pi[4] + 7.14285714285714*pi[5] + 1.19047619047619*pi[6])*alpha;
    w5.second.tail<3>() = 1.0*( + 0.11904761904762*Cg*T2*pi[2] + 0.357142857142857*Cg*T2*pi[4] + 0.0476190476190476*Cg*T2*pi[5]  - 0.238095238095238*Cpi[0]*pi[6] - 5.71428571428572*Cpi[1]*pi[5] + 1.42857142857143*Cpi[1]*pi[6] - 17.8571428571429*Cpi[2]*pi[4] + 12.8571428571429*Cpi[2]*pi[5])*alpha;
    wps.push_back(w5);
    waypoint6_t w6 = initwp<waypoint6_t>();
    w6.first.block<3,3>(0,0) = -5.71428571428571*alpha*Matrix3::Identity();
    w6.first.block<3,3>(3,0) = 1.0*(0.238095238095238*Cg*T2 + 14.2857142857143*Cpi[4] - 20.0*Cpi[5])*alpha;
    w6.second.head<3>() = 1.0*(1.42857142857143*pi[1] + 7.85714285714286*pi[2] - 14.2857142857143*pi[4] + 7.14285714285715*pi[5] + 3.57142857142857*pi[6])*alpha;
    w6.second.tail<3>() = 1.0*(0.535714285714286*Cg*T2*pi[4] + 0.214285714285714*Cg*T2*pi[5] + 0.0119047619047619*Cg*T2*pi[6] - 1.42857142857143*Cpi[1]*pi[6]  - 12.8571428571429*Cpi[2]*pi[5] + 5.0*Cpi[2]*pi[6])*alpha;
    wps.push_back(w6);
    waypoint6_t w7 = initwp<waypoint6_t>();
    w7.first.block<3,3>(0,0) = 6.66666666666667*alpha*Matrix3::Identity();
    w7.first.block<3,3>(3,0) = 1.0*( 20.0*Cpi[5] - 13.3333333333333*Cpi[6])*alpha;
    w7.second.head<3>() = 1.0*(5.0*pi[2] - 20.0*pi[4]  + 8.33333333333333*pi[6])*alpha;
    w7.second.tail<3>() = 1.0*( 0.416666666666667*Cg*T2*pi[4] + 0.5*Cg*T2*pi[5] + 0.0833333333333333*Cg*T2*pi[6]  - 5.0*Cpi[2]*pi[6] + 20.0*Cpi[4]*pi[5])*alpha;
    wps.push_back(w7);
    waypoint6_t w8 = initwp<waypoint6_t>();
    w8.first.block<3,3>(0,0) = 13.3333333333333*alpha*Matrix3::Identity();
    w8.first.block<3,3>(3,0) = 1.0*( 13.3333333333333*Cpi[6])*alpha;
    w8.second.head<3>() = 1.0*(-9.99999999999999*pi[4] - 20.0*pi[5] + 16.6666666666667*pi[6])*alpha;
    w8.second.tail<3>() = 1.0*( 0.666666666666667*Cg*T2*pi[5] + 0.333333333333333*Cg*T2*pi[6]  - 20.0*Cpi[4]*pi[5] + 30.0*Cpi[4]*pi[6])*alpha;
    wps.push_back(w8);
    waypoint6_t w9 = initwp<waypoint6_t>();
    w9.second.head<3>() = (30*pi[4] - 60*pi[5] + 30*pi[6])*alpha;
    w9.second.tail<3>() = 1.0*(1.0*Cg*T2*pi[6] - 30.0*Cpi[4]*pi[6] + 60.0*Cpi[5]*pi[6])*alpha;
    wps.push_back(w9);
    return wps;
}

coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T){
    coefs_t v;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    // equation found with sympy
    v.first = 0.;
    v.second = (-6.0*pi[5] + 6.0*pi[6])/ T;
    return v;
}

coefs_t computeFinalAccelerationPoint(const ProblemData& pData,double T){
    coefs_t v;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    // equation found with sympy
    v.first = 0.;
    v.second = (-30.0*pi[4] - 60.*pi[5] + 30.*pi[6])/ (T*T);
    return v;
}
}

}
