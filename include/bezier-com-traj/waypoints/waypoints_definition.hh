/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/waypoints/waypoints_c0_dc0_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_dc1_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_dc1_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh>



namespace bezier_com_traj{
typedef waypoint3_t waypoint_t;
typedef std::pair<double,point3_t> coefs_t;

/**
  * This file is used to choose the correct expressions of the curves waypoints, depending on the options set in ProblemData.constraints
  */


/** @brief evaluateCurveAtTime compute the expression of the point on the curve c at t, defined by the waypoint pi and one free waypoint (x)
     * @param pi constant waypoints of the curve
     * @param t param (normalized !)
     * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
     */
coefs_t evaluateCurveAtTime(const ProblemData& pData, std::vector<point_t> pi,double t){
    if(c0_dc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_c1::evaluateCurveAtTime(pi,t);
    else if(c0_dc0_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_dc1_c1::evaluateCurveAtTime(pi,t);
    else if(c0_dc0_ddc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_c1::evaluateCurveAtTime(pi,t);
    else if(c0_dc0_ddc0_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_dc1_c1::evaluateCurveAtTime(pi,t);
    else if(c0_dc0_ddc0_ddc1_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_ddc1_dc1_c1::evaluateCurveAtTime(pi,t);
    else{
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}

/** @brief evaluateAccelerationCurveAtTime compute the expression of the point on the curve ddc at t, defined by the waypoint pi and one free waypoint (x)
     * @param pi constant waypoints of the curve
     * @param t param (normalized !)
     * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
     */
coefs_t evaluateAccelerationCurveAtTime(const ProblemData& pData, std::vector<point_t> pi,double T,double t){
    if(c0_dc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_c1::evaluateAccelerationCurveAtTime(pi,T,t);
    else if(c0_dc0_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_dc1_c1::evaluateAccelerationCurveAtTime(pi,T,t);
    else if(c0_dc0_ddc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_c1::evaluateAccelerationCurveAtTime(pi,T,t);
    else if(c0_dc0_ddc0_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_dc1_c1::evaluateAccelerationCurveAtTime(pi,T,t);
    else if(c0_dc0_ddc0_ddc1_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_ddc1_dc1_c1::evaluateAccelerationCurveAtTime(pi,T,t);
    else{
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}

/**
 * @brief computeConstantWaypoints compute the constant waypoints of c(t) defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
    if(c0_dc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_c1::computeConstantWaypoints(pData,T);
    else if(c0_dc0_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_dc1_c1::computeConstantWaypoints(pData,T);
    else if(c0_dc0_ddc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_c1::computeConstantWaypoints(pData,T);
    else if(c0_dc0_ddc0_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_dc1_c1::computeConstantWaypoints(pData,T);
    else if(c0_dc0_ddc0_ddc1_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_ddc1_dc1_c1::computeConstantWaypoints(pData,T);
    else{
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}

/**
 * @brief computeWwaypoints compute the constant waypoints of w(t) defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
std::vector<waypoint6_t> computeWwaypoints(const ProblemData& pData,double T){
    if(c0_dc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_c1::computeWwaypoints(pData,T);
    else if(c0_dc0_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_dc1_c1::computeWwaypoints(pData,T);
    else if(c0_dc0_ddc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_c1::computeWwaypoints(pData,T);
    else if(c0_dc0_ddc0_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_dc1_c1::computeWwaypoints(pData,T);
    else if(c0_dc0_ddc0_ddc1_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_ddc1_dc1_c1::computeWwaypoints(pData,T);
    else{
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}

coefs_t computeFinalAccelerationPoint(const ProblemData& pData,double T){
    if(c0_dc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_c1::computeFinalAccelerationPoint(pData,T);
    else if(c0_dc0_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_dc1_c1::computeFinalAccelerationPoint(pData,T);
    else if(c0_dc0_ddc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_c1::computeFinalAccelerationPoint(pData,T);
    else if(c0_dc0_ddc0_dc1_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_dc1_c1::computeFinalAccelerationPoint(pData,T);
    else{
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}

coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T){
    if(c0_dc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_c1::computeFinalVelocityPoint(pData,T);
    else if(c0_dc0_ddc0_c1::useThisConstraints(pData.constraints_))
        return c0_dc0_ddc0_c1::computeFinalVelocityPoint(pData,T);
    else{
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}

void computeFinalAcceleration(const ProblemData& pData,double T,ResultDataCOMTraj& res){
    if(pData.constraints_.ddc1_){
        res.ddc1_ = pData.ddc1_;
    }else{
        coefs_t a = computeFinalAccelerationPoint(pData,T);
        res.ddc1_ = a.first*res.x + a.second;
    }
}

void computeFinalVelocity(const ProblemData& pData,double T,ResultDataCOMTraj& res){
    if(pData.constraints_.dc1_)
        res.dc1_ = pData.dc1_;
    else{
        coefs_t v = computeFinalVelocityPoint(pData,T);
        res.dc1_ = v.first*res.x + v.second;
    }
}


}
