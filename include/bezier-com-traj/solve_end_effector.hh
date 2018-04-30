/*
 * Copyright 2017, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>
#include <limits>

using namespace bezier_com_traj;

namespace bezier_com_traj
{
typedef waypoint3_t waypoint_t;
typedef std::pair<double,point3_t> coefs_t;
const int DIM_POINT=3;
const int NUM_DISCRETIZATION = 11;
const bool verbose = false;

/**
* @brief solveEndEffector Tries to produce a trajectory represented as a bezier curve
* that satisfy position, velocity and acceleration constraint for the initial and final point
* and that follow as close as possible the input trajectory
* @param pData problem Data.
* @param path the path to follow, the class Path must implement the operator (double t) , t \in [0,1] return a Vector3
* that give the position on the path for a given time
* @param T time lenght of the trajectory
* @param timeStep time that the solver has to stop
* @return ResultData a struct containing the resulting trajectory, if success is true.
*/
template<typename Path>
ResultDataCOMTraj solveEndEffector(const ProblemData& pData,const Path& path, const double T, const double weightDistance, bool useVelCost = true);


coefs_t initCoefs(){
    coefs_t c;
    c.first=0;
    c.second=point3_t::Zero();
    return c;
}


void computeConstantWaypoints(const ProblemData& pData,double T,double n,point_t& p0,point_t& p1, point_t& p2, point_t& p4, point_t& p5, point_t& p6){
    p0 = pData.c0_;
    p1 = (pData.dc0_ * T / n )+  pData.c0_;
    p2 = (pData.ddc0_*T*T/(n*(n-1))) + (2*pData.dc0_ *T / n) + pData.c0_; // * T because derivation make a T appear
    p6 = pData.c1_;
    p5 = (-pData.dc1_ * T / n) + pData.c1_; // * T ?
    p4 = (pData.ddc1_ *T*T / (n*(n-1))) - (2 * pData.dc1_ *T / n) + pData.c1_ ; // * T ??
}

void computeConstantWaypoints(const ProblemData& pData,double T,double n,point_t& p0,point_t& p1, point_t& p2, point_t& p3, point_t& p5, point_t& p6,point_t& p7,point_t& p8){
    p0 = pData.c0_;
    p1 = (pData.dc0_ * T / n )+  pData.c0_;
    p2 = (pData.ddc0_*T*T/(n*(n-1))) + (2*pData.dc0_ *T / n) + pData.c0_; // * T because derivation make a T appear
    p3 = (pData.j0_*T*T*T/(n*(n-1)*(n-2)))+ (3*pData.ddc0_*T*T/(n*(n-1))) + (3*pData.dc0_ *T / n) + pData.c0_;

    p8 = pData.c1_;
    p7 = (-pData.dc1_ * T / n) + pData.c1_; // * T ?
    p6 = (pData.ddc1_ *T*T / (n*(n-1))) - (2 * pData.dc1_ *T / n) + pData.c1_ ; // * T ??
    p5 = (-pData.j1_*T*T*T/(n*(n-1)*(n-2))) + (3*pData.ddc1_ *T*T / (n*(n-1))) - (3 * pData.dc1_ *T / n) + pData.c1_ ; // * T ??
}

// with jerk and jerk derivative constrained to 0
void computeConstantWaypoints(const ProblemData& pData,double T,double n,point_t& p0,point_t& p1, point_t& p2, point_t& p3, point_t& p4, point_t& p5,point_t& p6,point_t& p7,point_t& p8,point_t& p9){
    p0 = pData.c0_;
    p1 = (pData.dc0_ * T / n )+  pData.c0_;
    p2 = (n*n*pData.c0_ - n*pData.c0_ + 2*n*pData.dc0_*T - 2*pData.dc0_*T + pData.ddc0_*T*T)/(n*(n - 1)); // * T because derivation make a T appear
    p3 = (n*n*pData.c0_ - n*pData.c0_ + 3*n*pData.dc0_*T - 3*pData.dc0_*T + 3*pData.ddc0_*T*T)/(n*(n - 1));
    p4 = (n*n*pData.c0_ - n*pData.c0_ + 4*n*pData.dc0_*T - 4*pData.dc0_ *T+ 6*pData.ddc0_*T*T)/(n*(n - 1)) ;

    p9 = pData.c1_;
    p8 = (-pData.dc1_ * T / n) + pData.c1_; // * T ?
    p7 = (n*n*pData.c1_ - n*pData.c1_ - 2*n*pData.dc1_*T + 2*pData.dc1_*T + pData.ddc1_*T*T)/(n*(n - 1)) ; // * T ??
    p6 = (n*n*pData.c1_ - n*pData.c1_ - 3*n*pData.dc1_*T + 3*pData.dc1_*T + 3*pData.ddc1_*T*T)/(n*(n - 1)) ; // * T ??
    p5 = (n*n*pData.c1_ - n*pData.c1_ - 4*n*pData.dc1_*T + 4*pData.dc1_*T + 6*pData.ddc1_*T*T)/(n*(n - 1)) ;
}


// with jerk and jerk second derivative constrained to 0
void computeConstantWaypoints(const ProblemData& pData,double T,double n,point_t& p0,point_t& p1, point_t& p2, point_t& p3, point_t& p4, point_t& p5,point_t& p6,point_t& p7,point_t& p8,point_t& p9,point_t& p10,point_t& p11){
    p0 = pData.c0_;
    p1 = (pData.dc0_ * T / n )+  pData.c0_;
    p2 = (n*n*pData.c0_ - n*pData.c0_ + 2*n*pData.dc0_*T - 2*pData.dc0_*T + pData.ddc0_*T*T)/(n*(n - 1)); // * T because derivation make a T appear
    p3 = (n*n*pData.c0_ - n*pData.c0_ + 3*n*pData.dc0_*T - 3*pData.dc0_*T + 3*pData.ddc0_*T*T)/(n*(n - 1));
    p4 = (n*n*pData.c0_ - n*pData.c0_ + 4*n*pData.dc0_*T - 4*pData.dc0_ *T+ 6*pData.ddc0_*T*T)/(n*(n - 1)) ;
    p5 = (n*n*pData.c0_ - n*pData.c0_ + 5*n*pData.dc0_*T - 5*pData.dc0_ *T+ 10*pData.ddc0_*T*T)/(n*(n - 1)) ;

    p11 = pData.c1_;
    p10 = (-pData.dc1_ * T / n) + pData.c1_; // * T ?
    p9 = (n*n*pData.c1_ - n*pData.c1_ - 2*n*pData.dc1_*T + 2*pData.dc1_*T + pData.ddc1_*T*T)/(n*(n - 1)) ; // * T ??
    p8 = (n*n*pData.c1_ - n*pData.c1_ - 3*n*pData.dc1_*T + 3*pData.dc1_*T + 3*pData.ddc1_*T*T)/(n*(n - 1)) ; // * T ??
    p7 = (n*n*pData.c1_ - n*pData.c1_ - 4*n*pData.dc1_*T + 4*pData.dc1_*T + 6*pData.ddc1_*T*T)/(n*(n - 1)) ;
    p6 = (n*n*pData.c1_ - n*pData.c1_ - 5*n*pData.dc1_*T + 5*pData.dc1_*T + 10*pData.ddc1_*T*T)/(n*(n - 1)) ;
}

// up to jerk second derivativ constraints for init, pos vel and acc constraint for goal
std::vector<bezier_t::point_t> computeConstantWaypointsInitPredef(const ProblemData& pData,double T){
    const double n = 8;
    std::vector<bezier_t::point_t> pts;
    pts.push_back(pData.c0_); // c0
    pts.push_back((pData.dc0_ * T / n )+  pData.c0_); //dc0
    pts.push_back((n*n*pData.c0_ - n*pData.c0_ + 2*n*pData.dc0_*T - 2*pData.dc0_*T + pData.ddc0_*T*T)/(n*(n - 1)));//ddc0 // * T because derivation make a T appear
    pts.push_back((n*n*pData.c0_ - n*pData.c0_ + 3*n*pData.dc0_*T - 3*pData.dc0_*T + 3*pData.ddc0_*T*T)/(n*(n - 1))); //j0
    pts.push_back((n*n*pData.c0_ - n*pData.c0_ + 4*n*pData.dc0_*T - 4*pData.dc0_ *T+ 6*pData.ddc0_*T*T)/(n*(n - 1))) ; //dj0
    pts.push_back((n*n*pData.c0_ - n*pData.c0_ + 5*n*pData.dc0_*T - 5*pData.dc0_ *T+ 10*pData.ddc0_*T*T)/(n*(n - 1))) ; //ddj0

    pts.push_back((n*n*pData.c1_ - n*pData.c1_ - 2*n*pData.dc1_*T + 2*pData.dc1_*T + pData.ddc1_*T*T)/(n*(n - 1))) ; //ddc1 // * T ??
    pts.push_back((-pData.dc1_ * T / n) + pData.c1_); // dc1
    pts.push_back(pData.c1_); //c1

    return pts;
}


// up to jerk second derivativ constraints for goal, pos vel and acc constraint for init
std::vector<bezier_t::point_t> computeConstantWaypointsGoalPredef(const ProblemData& pData,double T){
    const double n = 8;
    std::vector<bezier_t::point_t> pts;
    pts.push_back(pData.c0_); //c0
    pts.push_back((pData.dc0_ * T / n )+  pData.c0_); //dc0
    pts.push_back((n*n*pData.c0_ - n*pData.c0_ + 2*n*pData.dc0_*T - 2*pData.dc0_*T + pData.ddc0_*T*T)/(n*(n - 1))); //ddc0 // * T because derivation make a T appear

    pts.push_back((n*n*pData.c1_ - n*pData.c1_ - 5*n*pData.dc1_*T + 5*pData.dc1_*T + 10*pData.ddc1_*T*T)/(n*(n - 1))) ; //ddj1
    pts.push_back((n*n*pData.c1_ - n*pData.c1_ - 4*n*pData.dc1_*T + 4*pData.dc1_*T + 6*pData.ddc1_*T*T)/(n*(n - 1))) ; //dj1
    pts.push_back((n*n*pData.c1_ - n*pData.c1_ - 3*n*pData.dc1_*T + 3*pData.dc1_*T + 3*pData.ddc1_*T*T)/(n*(n - 1))) ; // j1
    pts.push_back((n*n*pData.c1_ - n*pData.c1_ - 2*n*pData.dc1_*T + 2*pData.dc1_*T + pData.ddc1_*T*T)/(n*(n - 1))) ; //ddc1 * T ??
    pts.push_back((-pData.dc1_ * T / n) + pData.c1_); // dc1
    pts.push_back(pData.c1_); //c1

    return pts;
}

std::vector<bezier_t::point_t> computeConstantWaypoints(const ProblemData& pData,double T,double n){
    std::vector<bezier_t::point_t> pts;

    if(n<=6){
        point_t p0,p1,p2,p3,p4,p5;
        computeConstantWaypoints(pData,T,n,p0,p1,p2,p3,p4,p5);
        pts.push_back(p0);
        pts.push_back(p1);
        pts.push_back(p2);
        pts.push_back(bezier_t::point_t());
        pts.push_back(p3);
        pts.push_back(p4);
        pts.push_back(p5);
    }else if(n<=8){
        point_t p0,p1,p2,p3,p4,p5,p6,p7;
        computeConstantWaypoints(pData,T,n,p0,p1,p2,p3,p4,p5,p6,p7);
        pts.push_back(p0);
        pts.push_back(p1);
        pts.push_back(p2);
        pts.push_back(p3);
        pts.push_back(bezier_t::point_t());
        pts.push_back(p4);
        pts.push_back(p5);
        pts.push_back(p6);
        pts.push_back(p7);
    }else if (n<=10){
        point_t p0,p1,p2,p3,p4,p5,p6,p7,p8,p9;
        computeConstantWaypoints(pData,T,n,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9);
        pts.push_back(p0);
        pts.push_back(p1);
        pts.push_back(p2);
        pts.push_back(p3);
        pts.push_back(p4);
        pts.push_back(bezier_t::point_t());
        pts.push_back(p5);
        pts.push_back(p6);
        pts.push_back(p7);
        pts.push_back(p8);
        pts.push_back(p9);
    }else if (n<=12){
        point_t p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11;
        computeConstantWaypoints(pData,T,n,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11);
        pts.push_back(p0);
        pts.push_back(p1);
        pts.push_back(p2);
        pts.push_back(p3);
        pts.push_back(p4);
        pts.push_back(p5);
        pts.push_back(bezier_t::point_t());
        pts.push_back(p6);
        pts.push_back(p7);
        pts.push_back(p8);
        pts.push_back(p9);
        pts.push_back(p10);
        pts.push_back(p11);
    }
    return pts;
}


std::vector<waypoint_t> createEndEffectorJerkWaypoints(double T,const ProblemData& pData,const std::vector<bezier_t::point_t> pi){
    // create the waypoint from the analytical expressions :
    std::vector<waypoint_t> wps;
    assert(pi.size() == 9);
    double alpha = 1. / (T*T*T);

    waypoint_t w = initwp<waypoint_t>();
    // assign w0:
    w.second= 336*(-pi[0] +3*pi[1] - 3*pi[2] + pi[3])*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w1:
    w.first = 336*alpha*Matrix3::Identity();
    w.second= 336*(-pi[1] + 3*pi[2] - 3*pi[3])*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w2:
    w.first = -3*336*alpha*Matrix3::Identity();
    w.second = 336*(-pi[2] + 3*pi[3] + pi[5])*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w3:
    w.first = 3*336*alpha*Matrix3::Identity();
    w.second = 336*(-pi[3] - 3*pi[5] + pi[6])*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w4:
    w.first = -336*alpha*Matrix3::Identity();
    w.second =  336*(3*pi[5] - 3*pi[6] + pi[7])*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w5:
    w.second= 336*(-pi[5] + 3*pi[6] - 3*pi[7] + pi[8])*alpha;
    wps.push_back(w);
    return wps;
}

std::vector<waypoint_t> createEndEffectorAccelerationWaypoints(double T,const ProblemData& pData,const std::vector<bezier_t::point_t> pi){
    // create the waypoint from the analytical expressions :
    std::vector<waypoint_t> wps;
    assert(pi.size() == 9);

    if(verbose){
      std::cout<<"Create end eff waypoints, constant waypoints = :"<<std::endl<<
                 "p0 = "<<pi[0].transpose()<<std::endl<<"p1 = "<<pi[1].transpose()<<std::endl<<"p2 = "<<pi[2].transpose()<<std::endl<<
                 "p3 = "<<pi[3].transpose()<<std::endl<<"pi[5] = "<<pi[5].transpose()<<std::endl<<"pi[6] = "<<pi[6].transpose()<<std::endl<<"pi[7] = "<<pi[7].transpose()<<"pi[8] = "<<pi[8].transpose()<<std::endl;
    }
    double alpha = 1. / (T*T);

    waypoint_t w = initwp<waypoint_t>();
    // assign w0:
    w.second= 56*alpha*(pi[0] - 2*pi[1] + pi[2]);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w1:
    w.second= 56*alpha*(pi[1] - 2*pi[2] + pi[3]);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w2:
    w.first = 56*alpha*Matrix3::Identity();
    w.second = (56*pi[2] - 112*pi[3])*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w3:
    w.first = -112*alpha*Matrix3::Identity();
    w.second = (56*pi[3] +56*pi[8])*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w4:
    w.first = 56*alpha*Matrix3::Identity();
    w.second = (-112*pi[5] + 56*pi[6])*alpha;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w5:
    w.second=56*alpha*(pi[5]-2*pi[6]+pi[7]);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w5:
    w.second=56*alpha*(pi[6]-2*pi[7]+pi[8]);
    wps.push_back(w);
    return wps;
}


std::vector<waypoint_t> createEndEffectorVelocityWaypoints(double T,const ProblemData& pData,const std::vector<bezier_t::point_t> pi){
    // create the waypoint from the analytical expressions :
    std::vector<waypoint_t> wps;
    assert(pi.size() == 9);

   /* std::cout<<"Create end eff waypoints, constant waypoints = :"<<std::endl<<
               "p0 = "<<p0.transpose()<<std::endl<<"p1 = "<<p1.transpose()<<std::endl<<"p2 = "<<p2.transpose()<<std::endl<<
               "p4 = "<<p4.transpose()<<std::endl<<"p5 = "<<p5.transpose()<<std::endl<<"pi[6] = "<<pi[6].transpose()<<std::endl;*/
    double alpha = 1. / (T);

    waypoint_t w = initwp<waypoint_t>();
    // assign w0:
    w.second= alpha*8*(-pi[0]+pi[1]);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w1:
    w.second = alpha*8*(-pi[1]+pi[2]);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w2:
    w.second = alpha*8*(-pi[2]+pi[3]);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w3:
    w.first = 8*alpha*Matrix3::Identity();
    w.second = alpha*-8*pi[3];
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w4:
    w.first = -8*alpha*Matrix3::Identity();
    w.second = alpha*8*pi[5];
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w5:
    w.second=alpha*8*(-pi[5]+pi[6]);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w6:
    w.second=alpha*8*(-pi[6]+pi[7]);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w7:
    w.second=alpha*8*(-pi[7]+pi[8]);
    wps.push_back(w);
    return wps;
}


void computeConstraintsMatrix(const ProblemData& pData,const std::vector<waypoint_t>& wps_acc,const std::vector<waypoint_t>& wps_vel,const VectorX& acc_bounds,const VectorX& vel_bounds,MatrixXX& A,VectorX& b,const std::vector<waypoint_t>& wps_jerk = std::vector<waypoint_t>(),const VectorX& jerk_bounds=VectorX(DIM_POINT)  ){
    assert(acc_bounds.size() == DIM_POINT && "Acceleration bounds should have the same dimension as the points");
    assert(vel_bounds.size() == DIM_POINT && "Velocity bounds should have the same dimension as the points");
    assert(jerk_bounds.size() == DIM_POINT && "Jerk bounds should have the same dimension as the points");

    int empty_acc=0;
    int empty_vel=0;
    int empty_jerk=0;
    for (std::vector<waypoint_t>::const_iterator wpcit = wps_acc.begin(); wpcit != wps_acc.end(); ++wpcit)
    {
        if(wpcit->first.isZero(std::numeric_limits<double>::epsilon()))
            empty_acc++;
    }
    for (std::vector<waypoint_t>::const_iterator wpcit = wps_vel.begin(); wpcit != wps_vel.end(); ++wpcit)
    {
        if(wpcit->first.isZero(std::numeric_limits<double>::epsilon()))
            empty_vel++;
    }
    for (std::vector<waypoint_t>::const_iterator wpcit = wps_jerk.begin(); wpcit != wps_jerk.end(); ++wpcit)
    {
        if(wpcit->first.isZero(std::numeric_limits<double>::epsilon()))
            empty_jerk++;
    }

    A = MatrixXX::Zero(2*DIM_POINT*(wps_acc.size()-empty_acc+wps_vel.size()-empty_vel+wps_jerk.size() - empty_jerk)+DIM_POINT,DIM_POINT); // *2 because we have to put the lower and upper bound for each one, +DIM_POINT for the constraint on x[z]
    b = VectorX::Zero(A.rows());
    int i = 0;
    //upper acc bounds
    for (std::vector<waypoint_t>::const_iterator wpcit = wps_acc.begin(); wpcit != wps_acc.end(); ++wpcit)
    {
        if(! wpcit->first.isZero(std::numeric_limits<double>::epsilon())){
            A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = wpcit->first;
            b.segment<DIM_POINT>(i*DIM_POINT)   = acc_bounds - wpcit->second;
            ++i;
        }
    }
    //lower acc bounds
    for (std::vector<waypoint_t>::const_iterator wpcit = wps_acc.begin(); wpcit != wps_acc.end(); ++wpcit)
    {
        if(! wpcit->first.isZero(std::numeric_limits<double>::epsilon())){
            A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = -wpcit->first;
            b.segment<DIM_POINT>(i*DIM_POINT)   = acc_bounds + wpcit->second;
            ++i;
        }
    }
    //upper velocity bounds
    for (std::vector<waypoint_t>::const_iterator wpcit = wps_vel.begin(); wpcit != wps_vel.end(); ++wpcit)
    {
        if(! wpcit->first.isZero(std::numeric_limits<double>::epsilon())){
            A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = wpcit->first;
            b.segment<DIM_POINT>(i*DIM_POINT)   = vel_bounds - wpcit->second;
            ++i;
        }
    }
    //lower velocity bounds
    for (std::vector<waypoint_t>::const_iterator wpcit = wps_vel.begin(); wpcit != wps_vel.end(); ++wpcit)
    {
        if(! wpcit->first.isZero(std::numeric_limits<double>::epsilon())){
            A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = -wpcit->first;
            b.segment<DIM_POINT>(i*DIM_POINT)   = vel_bounds + wpcit->second;
            ++i;
        }
    }

    //upper jerk bounds
    for (std::vector<waypoint_t>::const_iterator wpcit = wps_vel.begin(); wpcit != wps_vel.end(); ++wpcit)
    {
        if(! wpcit->first.isZero(std::numeric_limits<double>::epsilon())){
            A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = wpcit->first;
            b.segment<DIM_POINT>(i*DIM_POINT)   = vel_bounds - wpcit->second;
            ++i;
        }
    }
    //lower jerk bounds
    for (std::vector<waypoint_t>::const_iterator wpcit = wps_jerk.begin(); wpcit != wps_jerk.end(); ++wpcit)
    {
        if(! wpcit->first.isZero(std::numeric_limits<double>::epsilon())){
            A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = -wpcit->first;
            b.segment<DIM_POINT>(i*DIM_POINT)   = jerk_bounds + wpcit->second;
            ++i;
        }
    }

    // test : constraint x[z] to be always higher than init[z] and goal[z].
    // TODO replace z with the direction of the contact normal ... need to change the API
    MatrixXX mxz = MatrixXX::Zero(DIM_POINT,DIM_POINT);
    mxz(DIM_POINT-1,DIM_POINT-1) = -1;
    VectorX nxz = VectorX::Zero(DIM_POINT);
    nxz[2]= - std::min(pData.c0_[2],pData.c1_[2]);
    A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = mxz;
    b.segment<DIM_POINT>(i*DIM_POINT)   = nxz;


    //TEST :
  /*  A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = Matrix3::Identity();
    b.segment<DIM_POINT>(i*DIM_POINT)   = Vector3(10,10,10);
    i++;
    A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = -Matrix3::Identity();
    b.segment<DIM_POINT>(i*DIM_POINT)   =  Vector3(10,10,10);*/
}

/*std::vector<coefs_t> createDiscretizationPoints(const ProblemData& pData){
    // equation found with sympy (for 11 points)
    std::vector<coefs_t> cks; // one element for each discretization points :
    //.first is the term that depend on x and . second is the constant term
    coefs_t ck = initCoefs();
    // 0
    ck.first = 0.;
    ck.second =  pData.c0_;
    cks.push_back(ck);
    ck = initCoefs();
    // 1
    ck.first =  729./50000.;
    ck.second =  19683.*pData.c0_/20000. + 127.*pData.c1_/100000. + 45927.*pData.dc0_/500000. - 207.*pData.dc1_/500000. + 6561.*pData.ddc0_/2000000. + 81.*pData.ddc1_/2000000.;
    cks.push_back(ck);
    ck = initCoefs();
    // 2
    ck.first =  256./3125.;
    ck.second =  2816.*pData.c0_/3125. + 53.*pData.c1_/3125. + 2304.*pData.dc0_/15625. - 84.*pData.dc1_/15625. + 128.*pData.ddc0_/15625. + 8.*pData.ddc1_/15625. ;
    cks.push_back(ck);
    ck = initCoefs();
    // 3
    ck.first = 9261./50000. ;
    ck.second =  74431.*pData.c0_/100000. + 7047.*pData.c1_/100000. + 79233.*pData.dc0_/500000. - 10773.*pData.dc1_/500000. + 21609.*pData.ddc0_/2000000. + 3969.*pData.ddc1_/2000000.;
    cks.push_back(ck);
    ck = initCoefs();
    // 4
    ck.first =  864./3125.;
    ck.second =  1701.*pData.c0_/3125. + 112.*pData.c1_/625.+ 2106.*pData.dc0_/15625. - 816.*pData.dc1_/15625. + 162.*pData.ddc0_/15625. + 72.*pData.ddc1_/15625.;
    cks.push_back(ck);
    ck = initCoefs();
    // 5
    ck.first =  5./16.;
    ck.second = 11.*pData.c0_/32. + 11.*pData.c1_/32. + 3.*pData.dc0_/32. - 3.*pData.dc1_/32. + pData.ddc0_/128. + pData.ddc1_/128. ;
    cks.push_back(ck);
    ck = initCoefs();
    // 6
    ck.first =  864./3125.;
    ck.second = 112.*pData.c0_/625. + 1701.*pData.c1_/3125. + 816.*pData.dc0_/15625. - 2106.*pData.dc1_/15625. + 72.*pData.ddc0_/15625. + 162.*pData.ddc1_/15625. ;
    cks.push_back(ck);
    ck = initCoefs();
    // 7
    ck.first =  9261./50000.;
    ck.second =  704699999999999.*pData.c0_/10000000000000000. + 74431.*pData.c1_/100000. + 10773.*pData.dc0_/500000. - 79233.*pData.dc1_/500000. + 3969.*pData.ddc0_/2000000. + 21609.*pData.ddc1_/2000000.;
    cks.push_back(ck);
    ck = initCoefs();
    // 8
    ck.first =  256./3125.;
    ck.second =   53.*pData.c0_/3125. + 2816.*pData.c1_/3125. + 84.*pData.dc0_/15625. - 2304.*pData.dc1_/15625. + 8.*pData.ddc0_/15625. + 128.*pData.ddc1_/15625.;
    cks.push_back(ck);
    ck = initCoefs();
    // 9
    ck.first =  729./50000.;
    ck.second =  127.*pData.c0_/100000. + 19683.*pData.c1_/20000. + 207.*pData.dc0_/500000. - 45927.*pData.dc1_/500000. + 81.*pData.ddc0_/2000000. + 6561.*pData.ddc1_/2000000. ;
    cks.push_back(ck);
    ck = initCoefs();
    // 10
    ck.first =  0.;
    ck.second = pData.c1_ ;
    cks.push_back(ck);
    return cks;
}*/


/**
 * @brief evaluateCurve Evaluate the curve at a given parameter
 * @param pData
 * @param T
 * @param t Normalized : between 0 and 1
 * @return
 */
coefs_t evaluateCurve(const ProblemData& pData,double T, double t){
    point_t p0,p1,p2,p3,p5,p6,p7,p8;
    computeConstantWaypoints(pData,T,8,p0,p1,p2,p3,p5,p6,p7,p8);
    coefs_t coefs;
    coefs.first = 70.0*pow(t,8) - 280.0*pow(t,7) + 420.0*pow(t,6) - 280.0*pow(t,5) + 70.0*pow(t,4);
    coefs.second = 1.0*p0*pow(t,8) - 8.0*p0*pow(t,7) + 28.0*p0*pow(t,6) - 56.0*p0*pow(t,5) + 70.0*p0*pow(t,4) - 56.0*p0*pow(t,3) + 28.0*p0*pow(t,2) - 8.0*p0*t + 1.0*p0 - 8.0*p1*pow(t,8) + 56.0*p1*pow(t,7) - 168.0*p1*pow(t,6) + 280.0*p1*pow(t,5) - 280.0*p1*pow(t,4) + 168.0*p1*pow(t,3) - 56.0*p1*pow(t,2) + 8.0*p1*t + 28.0*p2*pow(t,8) - 168.0*p2*pow(t,7) + 420.0*p2*pow(t,6) - 560.0*p2*pow(t,5) + 420.0*p2*pow(t,4) - 168.0*p2*pow(t,3) + 28.0*p2*pow(t,2) - 56.0*p3*pow(t,8 )+ 280.0*p3*pow(t,7) - 560.0*p3*pow(t,6) + 560.0*p3*pow(t,5) - 280.0*p3*pow(t,4) + 56.0*p3*pow(t,3) - 56.0*p5*pow(t,8) + 168.0*p5*pow(t,7) - 168.0*p5*pow(t,6) + 56.0*p5*pow(t,5 )+ 28.0*p6*pow(t,8) - 56.0*p6*pow(t,7) + 28.0*p6*pow(t,6) - 8.0*p7*pow(t,8) + 8.0*p7*pow(t,7) + 1.0*p8*pow(t,8);
    return coefs;
}


/**
 * @brief evaluateAccCurve Evaluate the acceleration at a given parameter
 * @param pData
 * @param T
 * @param param Normalized : between 0 and 1
 * @return
 */
coefs_t evaluateVelCurve(const ProblemData& pData, double T, double t){
    point_t p0,p1,p2,p3,p5,p6,p7,p8;
    computeConstantWaypoints(pData,T,8,p0,p1,p2,p3,p5,p6,p7,p8);
    coefs_t coefs;
    const double alpha = 1./(T);
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    const double t6 = t5*t;
    const double t7 = t6*t;
    //equations found with sympy
    coefs.first= (560.0*t7 - 1960.0*t6 + 2520.0*t5 - 1400.0*t4 + 280.0*t3)*alpha;
    coefs.second=(8.0*p0*t7 - 56.0*p0*t6 + 168.0*p0*t5 - 280.0*p0*t4 + 280.0*p0*t3 - 168.0*p0*t2 + 56.0*p0*t - 8.0*p0 - 64.0*p1*t7 + 392.0*p1*t6 - 1008.0*p1*t5 + 1400.0*p1*t4 - 1120.0*p1*t3 + 504.0*p1*t2 - 112.0*p1*t + 8.0*p1 + 224.0*p2*t7 - 1176.0*p2*t6 + 2520.0*p2*t5 - 2800.0*p2*t4 + 1680.0*p2*t3 - 504.0*p2*t2 + 56.0*p2*t - 448.0*p3*t7 + 1960.0*p3*t6 - 3360.0*p3*t5 + 2800.0*p3*t4 - 1120.0*p3*t3 + 168.0*p3*t2 - 448.0*p5*t7 + 1176.0*p5*t6 - 1008.0*p5*t5 + 280.0*p5*t4 + 224.0*p6*t7 - 392.0*p6*t6 + 168.0*p6*t5 - 64.0*p7*t7 + 56.0*p7*t6 + 8.0*p8*t7)*alpha;
    return coefs;
}


/**
 * @brief evaluateAccCurve Evaluate the acceleration at a given parameter
 * @param pData
 * @param T
 * @param param Normalized : between 0 and 1
 * @return
 */
coefs_t evaluateAccCurve(const ProblemData& pData, double T, double t){
    point_t p0,p1,p2,p3,p5,p6,p7,p8;
    computeConstantWaypoints(pData,T,8,p0,p1,p2,p3,p5,p6,p7,p8);
    coefs_t coefs;
    const double alpha = 1./(T*T);
    //equations found with sympy
    coefs.first= ((3920.0*pow(t,6) - 11760.0*pow(t,5) + 12600.0*pow(t,4) - 5600.0*pow(t,3) + 840.0*pow(t,2)))*alpha;
    coefs.second=(56.0*p0*pow(t,6) - 336.0*p0*pow(t,5) + 840.0*p0*pow(t,4) - 1120.0*p0*pow(t,3) + 840.0*p0*pow(t,2) - 336.0*p0*t + 56.0*p0 - 448.0*p1*pow(t,6) + 2352.0*p1*pow(t,5) - 5040.0*p1*pow(t,4) + 5600.0*p1*pow(t,3) - 3360.0*p1*pow(t,2) + 1008.0*p1*t - 112.0*p1 + 1568.0*p2*pow(t,6) - 7056.0*p2*pow(t,5) + 12600.0*p2*pow(t,4) - 11200.0*p2*pow(t,3) + 5040.0*p2*pow(t,2) - 1008.0*p2*t + 56.0*p2 - 3136.0*p3*pow(t,6) + 11760.0*p3*pow(t,5) - 16800.0*p3*pow(t,4) + 11200.0*p3*pow(t,3) - 3360.0*p3*pow(t,2)+ 336.0*p3*t - 3136.0*p5*pow(t,6) + 7056.0*p5*pow(t,5) - 5040.0*p5*pow(t,4) + 1120.0*p5*pow(t,3) + 1568.0*p6*pow(t,6) - 2352.0*p6*pow(t,5) + 840.0*p6*pow(t,4) - 448.0*p7*pow(t,6) + 336.0*p7*pow(t,5) + 56.0*p8*pow(t,6))*alpha;
    return coefs;
}


/**
 * @brief evaluateAccCurve Evaluate the acceleration at a given parameter
 * @param pData
 * @param T
 * @param param Normalized : between 0 and 1
 * @return
 */
coefs_t evaluateJerkCurve(const ProblemData& pData, double T, double t){
    point_t p0,p1,p2,p3,p5,p6,p7,p8;
    computeConstantWaypoints(pData,T,8,p0,p1,p2,p3,p5,p6,p7,p8);
    coefs_t coefs;
    const double t2 = t*t;
    const double t3 = t2*t;
    const double t4 = t3*t;
    const double t5 = t4*t;
    const double alpha = 1./(T*T*T);

    //equations found with sympy
    coefs.first= (23520.0*t5 - 58800.0*t4 + 50400.0*t3 - 16800.0*t2 + 1680.0*t)*alpha;
    coefs.second= 1.0*(336.0*p0*t5 - 1680.0*p0*t4 + 3360.0*p0*t3 - 3360.0*p0*t2 + 1680.0*p0*t - 336.0*p0 - 2688.0*p1*t5 + 11760.0*p1*t4 - 20160.0*p1*t3 + 16800.0*p1*t2 - 6720.0*p1*t + 1008.0*p1 + 9408.0*p2*t5 - 35280.0*p2*t4 + 50400.0*p2*t3 - 33600.0*p2*t2 + 10080.0*p2*t - 1008.0*p2 - 18816.0*p3*t5 + 58800.0*p3*t4 - 67200.0*p3*t3 + 33600.0*p3*t2 - 6720.0*p3*t + 336.0*p3 - 18816.0*p5*t5 + 35280.0*p5*t4 - 20160.0*p5*t3 + 3360.0*p5*t2 + 9408.0*p6*t5 - 11760.0*p6*t4 + 3360.0*p6*t3 - 2688.0*p7*t5 + 1680.0*p7*t4 + 336.0*p8*t5)*alpha;
    return coefs;
}

template <typename Path>
void computeDistanceCostFunction(int numPoints,const ProblemData& pData, double T,const Path& path, MatrixXX& H,VectorX& g){
    H = MatrixXX::Zero(DIM_POINT,DIM_POINT);
    g  = VectorX::Zero(DIM_POINT);
    double step = 1./(numPoints-1);
    std::vector<coefs_t> cks;
    for(size_t i = 0 ; i < numPoints ; ++i){
        cks.push_back(evaluateCurve(pData,T,i*step));
    }
    point3_t pk;
    size_t i = 0;
    for (std::vector<coefs_t>::const_iterator ckcit = cks.begin(); ckcit != cks.end(); ++ckcit){
        pk=path(i*step);
      //  std::cout<<"pk = "<<pk.transpose()<<std::endl;
      //  std::cout<<"coef First : "<<ckcit->first<<std::endl;
      //  std::cout<<"coef second : "<<ckcit->second.transpose()<<std::endl;
        H += (ckcit->first * ckcit->first * Matrix3::Identity());
        g += (ckcit->first * ckcit->second) - (pk * ckcit->first);
        i++;
    }
    double norm=H(0,0); // because H is always diagonal.
    H /= norm;
    g /= norm;
}

void computeC_of_T (const ProblemData& pData,double T, ResultDataCOMTraj& res){
    std::vector<Vector3> wps;
    point_t p0,p1,p2,p3,p5,p6,p7,p8;
    computeConstantWaypoints(pData,T,8,p0,p1,p2,p3,p5,p6,p7,p8);
    wps.push_back(p0);
    wps.push_back(p1);
    wps.push_back(p2);
    wps.push_back(p3);
    wps.push_back(res.x);
    wps.push_back(p5);
    wps.push_back(p6);
    wps.push_back(p7);
    wps.push_back(p8);
    res.c_of_t_ = bezier_t (wps.begin(), wps.end(),T);
    if(verbose)
      std::cout<<"bezier curve created, size = "<<res.c_of_t_.size_<<std::endl;
}

void computeVelCostFunction(int numPoints,const ProblemData& pData,double T, MatrixXX& H,VectorX& g){
    double step = 1./(numPoints-1);
    H = MatrixXX::Zero(DIM_POINT,DIM_POINT);
    g  = VectorX::Zero(DIM_POINT);
    std::vector<coefs_t> cks;
    for(int i = 0 ; i < numPoints ; ++i){
        cks.push_back(evaluateVelCurve(pData,T,i*step));
    }
    for (std::vector<coefs_t>::const_iterator ckcit = cks.begin(); ckcit != cks.end(); ++ckcit){
        H+=(ckcit->first * ckcit->first * Matrix3::Identity());
        g+=ckcit->first*ckcit->second;
    }
    //TEST : don't consider z axis for minimum acceleration cost
    //H(2,2) = 1e-6;
    //g[2] = 1e-6 ;
    //normalize :
    double norm=H(0,0); // because H is always diagonal
    H /= norm;
    g /= norm;
}

void computeAccelerationCostFunction(int numPoints,const ProblemData& pData,double T, MatrixXX& H,VectorX& g){
    double step = 1./(numPoints-1);
    H = MatrixXX::Zero(DIM_POINT,DIM_POINT);
    g  = VectorX::Zero(DIM_POINT);
    std::vector<coefs_t> cks;
    for(int i = 0 ; i < numPoints ; ++i){
        cks.push_back(evaluateAccCurve(pData,T,i*step));
    }
    for (std::vector<coefs_t>::const_iterator ckcit = cks.begin(); ckcit != cks.end(); ++ckcit){
        H+=(ckcit->first * ckcit->first * Matrix3::Identity());
        g+=ckcit->first*ckcit->second;
    }
    //TEST : don't consider z axis for minimum acceleration cost
    //H(2,2) = 1e-6;
    //g[2] = 1e-6 ;
    //normalize :
    double norm=H(0,0); // because H is always diagonal
    H /= norm;
    g /= norm;
}

void computeJerkCostFunction(int numPoints,const ProblemData& pData,double T, MatrixXX& H,VectorX& g){
    double step = 1./(numPoints-1);
    H = MatrixXX::Zero(DIM_POINT,DIM_POINT);
    g  = VectorX::Zero(DIM_POINT);
    std::vector<coefs_t> cks;
    for(int i = 0 ; i < numPoints ; ++i){
        cks.push_back(evaluateJerkCurve(pData,T,i*step));
    }
    for (std::vector<coefs_t>::const_iterator ckcit = cks.begin(); ckcit != cks.end(); ++ckcit){
        H+=(ckcit->first * ckcit->first * Matrix3::Identity());
        g+=ckcit->first*ckcit->second;
    }
    //TEST : don't consider z axis for minimum acceleration cost
    //H(2,2) = 1e-6;
    //g[2] = 1e-6 ;
    //normalize :
    double norm=H(0,0); // because H is always diagonal
    H /= norm;
    g /= norm;
}


template <typename Path>
ResultDataCOMTraj solveEndEffector(const ProblemData& pData,const Path& path, const double T, const double weightDistance, bool useVelCost){

    if(verbose)
      std::cout<<"solve end effector, T = "<<T<<std::endl;
    assert (weightDistance>=0. && weightDistance<=1. && "WeightDistance must be between 0 and 1");
    double weightSmooth = 1. - weightDistance;
    std::vector<bezier_t::point_t> pi = computeConstantWaypoints(pData,T,8);
    std::vector<waypoint_t> wps_jerk=createEndEffectorJerkWaypoints(T,pData,pi);
    std::vector<waypoint_t> wps_acc=createEndEffectorAccelerationWaypoints(T,pData,pi);
    std::vector<waypoint_t> wps_vel=createEndEffectorVelocityWaypoints(T,pData,pi);
    // stack the constraint for each waypoint :
    MatrixXX A;
    VectorX b;
    Vector3 jerk_bounds(5000,5000,5000);
    Vector3 acc_bounds(5000,5000,5000);
    Vector3 vel_bounds(5000,5000,5000);
    computeConstraintsMatrix(pData,wps_acc,wps_vel,acc_bounds,vel_bounds,A,b,wps_jerk,jerk_bounds);
  //  std::cout<<"End eff A = "<<std::endl<<A<<std::endl;
 //   std::cout<<"End eff b = "<<std::endl<<b<<std::endl;
    // compute cost function (discrete integral under the curve defined by 'path')
    MatrixXX H_rrt=MatrixXX::Zero(DIM_POINT,DIM_POINT),H_acc,H_jerk,H_smooth,H;
    VectorX g_rrt=VectorX::Zero(DIM_POINT),g_acc,g_jerk,g_smooth,g;
    if(weightDistance>0)
        computeDistanceCostFunction<Path>(50,pData,T,path,H_rrt,g_rrt);
    if(useVelCost)
        computeVelCostFunction(50,pData,T,H_smooth,g_smooth);
    else
        computeJerkCostFunction(50,pData,T,H_smooth,g_smooth);
  /*  std::cout<<"End eff H_rrt = "<<std::endl<<H_rrt<<std::endl;
    std::cout<<"End eff g_rrt = "<<std::endl<<g_rrt<<std::endl;
    std::cout<<"End eff H_acc = "<<std::endl<<H_acc<<std::endl;
    std::cout<<"End eff g_acc = "<<std::endl<<g_acc<<std::endl;
*/

    // add the costs :
    H = MatrixXX::Zero(DIM_POINT,DIM_POINT);
    g  = VectorX::Zero(DIM_POINT);
    H = weightSmooth*(H_smooth) + weightDistance*H_rrt;
    g = weightSmooth*(g_smooth) + weightDistance*g_rrt;
    if(verbose){
      std::cout<<"End eff H = "<<std::endl<<H<<std::endl;
      std::cout<<"End eff g = "<<std::endl<<g<<std::endl;
    }
    // call the solver
    VectorX init = VectorX(DIM_POINT);
    init = (pData.c0_ + pData.c1_)/2.;
   // init =pData.c0_;
    if(verbose)
      std::cout<<"Init = "<<std::endl<<init.transpose()<<std::endl;
    ResultData resQp = solve(A,b,H,g, init);

    ResultDataCOMTraj res;
    if(resQp.success_)
    {
        res.success_ = true;
        res.x = resQp.x;
       // computeRealCost(pData, res);
        computeC_of_T (pData,T,res);
       // computedL_of_T(pData,Ts,res);
    }
   if(verbose){
     std::cout<<"Solved, success = "<<res.success_<<" x = "<<res.x.transpose()<<std::endl;
     std::cout<<"Final cost : "<<resQp.cost_<<std::endl;
    }
   return res;
}


}//namespace bezier_com_traj
