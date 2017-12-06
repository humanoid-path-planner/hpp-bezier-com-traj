/*
 * Copyright 2017, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>
using namespace bezier_com_traj;
typedef waypoint3_t waypoint_t;
const int DIM_POINT=3;
namespace bezier_com_traj
{

std::vector<waypoint_t> createEndEffectorWaypoints(double T,const ProblemData& pData){
    // create the waypoint from the analytical expressions :
    std::vector<waypoint_t> wps;
    double n = 6; // degree
    point_t p0 = pData.c0_;
    point_t p1 = pData.dc0_ * T / n +  pData.c0_;
    point_t p2 = pData.ddc0_*T*T/(n*(n-1)) + 2*pData.dc0_ *T / n + pData.c0_; // * T because derivation make a T appear
    point_t p6 = pData.c1_;
    point_t p5 = -pData.dc1_ * T / n + pData.c1_; // * T ?
    point_t p4 = pData.ddc1_ *T*T / (n*(n-1)) - 2 * pData.dc1_ *T / n + pData.c1_ ; // * T ??
    double alpha = 1. / (T*T);

    waypoint_t w = initwp<waypoint_t>();
    // assign w0:
    w.second= 30*alpha*(p0 - 2*p1 + p2);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w1:
    w.first = 30*alpha*Matrix3::Identity();
    w.second = (30*p1 - 60*p2)*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w2:
    w.first = -60*alpha*Matrix3::Identity();
    w.second = (30*p2 +30*p4)*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w3:
    w.first = 30*alpha*Matrix3::Identity();
    w.second = (-60*p4 + 30*p5)*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w4:
    w.second=30*alpha*(p4-2*p5+p6);
    wps.push_back(w);
    return wps;
}

template <typename Path>
ResultDataCOMTraj solveEndEffector(const ProblemData& pData,Path path, const double T, const double timeStep){
    std::vector<waypoint_t> wps=createEndEffectorWaypoints(T,pData);
    // stack the constraint for each waypoint :
    MatrixXX A = MatrixXX::Zero(DIM_POINT*wps.size(),DIM_POINT);
    VectorX b = VectorX::Zero(DIM_POINT);
    int i = 0;
    for (std::vector<waypoint_t>::const_iterator wpcit = wps.begin(); wpcit != wps.end(); ++wpcit)
    {
        A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = wpcit->first;
        b.segment<DIM_POINT>(i*DIM_POINT)   = -wpcit->second;
        ++i;
    }

}


}//namespace bezier_com_traj
