/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>

namespace bezier_com_traj
{
typedef waypoint3_t waypoint_t;

ResultData solveIntersection(const std::pair<MatrixXX, VectorX>& Ab,const std::pair<MatrixXX, VectorX>& Hg,  const Vector3& init)
{
    return solve(Ab.first,Ab.second,Hg.first,Hg.second, init);
}



/**
 * @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve, assume p0 p1 x p2 p3
 * @param t param
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
waypoint_t evaluateCurveAtTime(std::vector<point_t> pi,double t){
    waypoint_t wp;
    double t2 = t*t;
    double t3 = t2*t;
    double t4 = t3*t;
    // equation found with sympy
    wp.first = Matrix3::Identity() *( 6.0*t4 - 12.0*t3 + 6.0*t2);
    wp.second = 1.0*pi[0]*t4 - 4.0*pi[0]*t3 + 6.0*pi[0]*t2 - 4.0*pi[0]*t + 1.0*pi[0] - 4.0*pi[1]*t4 + 12.0*pi[1]*t3 - 12.0*pi[1]*t2 + 4.0*pi[1]*t - 4.0*pi[2]*t4 + 4.0*pi[2]*t3 + 1.0*pi[3]*t4;
    return wp;
}


std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
    // equation for constraint on initial and final position and velocity (degree 4, 4 constant waypoint and one free (p2))
    // first, compute the constant waypoints that only depend on pData :
    int n = 4;
    std::vector<point_t> pi;
    pi.push_back(pData.c0_); //p0
    pi.push_back((pData.dc0_ * T / n )+  pData.c0_); // p1
    pi.push_back((-pData.dc1_ * T / n) + pData.c1_); // p3
    pi.push_back(pData.c1_); // p4
    return pi;
}

std::vector<waypoint_t> computeDiscretizedWaypoints(const ProblemData& pData,double T,double timeStep){
    int numStep = int(T / timeStep);
    std::vector<waypoint_t> wps;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    double t = 0;
    while(t<T){
        wps.push_back(evaluateCurveAtTime(pi,t));
        t+= timeStep;
    }
    return wps;
}






std::pair<MatrixX3, VectorX> computeConstraintsOneStep(const ProblemData& pData,const std::vector<double>& Ts,const double timeStep){
    // compute the list of discretized waypoint :
    double t_total = 0.;
    for(int i = 0 ; i < Ts.size() ; ++i)
        t_total+=Ts[i];

    // Compute all the discretized wayPoint
}


std::pair<MatrixX3, VectorX> computeCostFunctionOneStep(const ProblemData&pData){

}


void computeBezierCurve(const ProblemData& pData, const std::vector<double>& Ts, ResultDataCOMTraj& res)
{
    std::vector<Vector3> wps;
    double T = 0;
    for(int i = 0 ; i < Ts.size() ; ++i)
        T+=Ts[i];

    std::vector<Vector3> pi = computeConstantWaypoints(pData,T);
    wps.push_back(pi[0]);
    wps.push_back(pi[1]);
    wps.push_back(res.x);
    wps.push_back(pi[2]);
    wps.push_back(pi[3]);
    res.c_of_t_ = bezier_t (wps.begin(), wps.end(),T);
}

ResultDataCOMTraj solveOnestep(const ProblemData& pData, const std::vector<double>& Ts, const double timeStep){
    assert(pData.contacts_.size() ==2);
    assert(Ts.size() == pData.contacts_.size());
    bool fail = true;
    ResultDataCOMTraj res;
    std::pair<MatrixX3, VectorX> Ab = computeConstraintsOneStep(pData,Ts,timeStep);
    std::pair<MatrixX3, VectorX> Hg = computeCostFunctionOneStep(pData);
    Vector3 midPoint = (pData.c0_ + pData.c1_)/2.;
    // rewriting 0.5 || Dx -d ||^2 as x'Hx  + g'x
    ResultData resQp = solve(Ab.first,Ab.second,Hg.first,Hg.second, midPoint);
    if(resQp.success_)
    {
        res.success_ = true;
        res.x = resQp.x;
        computeBezierCurve (pData,Ts,res);
    }
    return res;
}


} // namespace
