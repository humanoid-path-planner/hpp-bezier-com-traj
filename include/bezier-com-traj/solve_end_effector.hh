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
ResultDataCOMTraj solveEndEffector(const ProblemData& pData,const Path& path, const double T, const double timeStep);


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

std::vector<bezier_t::point_t> computeConstantWaypoints(const ProblemData& pData,double T,double n){
    point_t p0,p1,p2,p3,p4,p5;
    computeConstantWaypoints(pData,T,n,p0,p1,p2,p3,p4,p5);
    std::vector<bezier_t::point_t> pts;
    pts.push_back(p0);
    pts.push_back(p1);
    pts.push_back(p2);
    pts.push_back(p3);
    pts.push_back(p4);
    pts.push_back(p5);
    return pts;
}


std::vector<waypoint_t> createEndEffectorAccelerationWaypoints(double T,const ProblemData& pData){
    // create the waypoint from the analytical expressions :
    std::vector<waypoint_t> wps;
    point_t p0,p1,p2,p4,p5,p6;
    computeConstantWaypoints(pData,T,6,p0,p1,p2,p4,p5,p6);
    std::cout<<"Create end eff waypoints, constant waypoints = :"<<std::endl<<
               "p0 = "<<p0.transpose()<<std::endl<<"p1 = "<<p1.transpose()<<std::endl<<"p2 = "<<p2.transpose()<<std::endl<<
               "p4 = "<<p4.transpose()<<std::endl<<"p5 = "<<p5.transpose()<<std::endl<<"p6 = "<<p6.transpose()<<std::endl;
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


std::vector<waypoint_t> createEndEffectorVelocityWaypoints(double T,const ProblemData& pData){
    // create the waypoint from the analytical expressions :
    std::vector<waypoint_t> wps;
    point_t p0,p1,p2,p4,p5,p6;
    computeConstantWaypoints(pData,T,6,p0,p1,p2,p4,p5,p6);
   /* std::cout<<"Create end eff waypoints, constant waypoints = :"<<std::endl<<
               "p0 = "<<p0.transpose()<<std::endl<<"p1 = "<<p1.transpose()<<std::endl<<"p2 = "<<p2.transpose()<<std::endl<<
               "p4 = "<<p4.transpose()<<std::endl<<"p5 = "<<p5.transpose()<<std::endl<<"p6 = "<<p6.transpose()<<std::endl;*/
    double alpha = 1. / (T);

    waypoint_t w = initwp<waypoint_t>();
    // assign w0:
    w.second= alpha*6*(-p0+p1);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w1:
    w.second = alpha*6*(-p1+p2);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w2:
    w.first = 6*alpha*Matrix3::Identity();
    w.second = alpha*-6*p2;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w3:
    w.first = -6*alpha*Matrix3::Identity();
    w.second = alpha*6*p4;
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w4:
    w.second=alpha*6*(-p4+p5);
    wps.push_back(w);
    w = initwp<waypoint_t>();
    // assign w5:
    w.second=alpha*6*(-p5+p6);
    wps.push_back(w);
    return wps;
}


void computeConstraintsMatrix(const std::vector<waypoint_t>& wps_acc,const std::vector<waypoint_t>& wps_vel,const VectorX& acc_bounds,const VectorX& vel_bounds,MatrixXX& A,VectorX& b){
    assert(acc_bounds.length() == DIM_POINT && "Acceleration bounds should have the same dimension as the points");
    assert(vel_bounds.length() == DIM_POINT && "Velocity bounds should have the same dimension as the points");
    int empty_acc=0;
    int empty_vel=0;
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

    A = MatrixXX::Zero(2*DIM_POINT*(wps_acc.size()-empty_acc+wps_vel.size()-empty_vel),DIM_POINT); // *2 because we have to put the lower and upper bound for each one
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

    //TEST :
  /*  A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = Matrix3::Identity();
    b.segment<DIM_POINT>(i*DIM_POINT)   = Vector3(10,10,10);
    i++;
    A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = -Matrix3::Identity();
    b.segment<DIM_POINT>(i*DIM_POINT)   =  Vector3(10,10,10);*/
}

std::vector<coefs_t> createDiscretizationPoints(const ProblemData& pData){
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
}

/**
 * @brief evaluateAccCurve Evaluate the acceleration at a given parameter
 * @param pData
 * @param T
 * @param param Normalized : between 0 and 1
 * @return
 */
coefs_t evaluateAccCurve(const ProblemData& pData, double T, double param){
    point_t p0,p1,p2,p4,p5,p6;
    computeConstantWaypoints(pData,T,6,p0,p1,p2,p4,p5,p6);
    coefs_t coefs;
    double alpha = 1./(T*T);
    //equations found with sympy
    coefs.first= (-600.0*param*param*param*param + 1200.0*param*param*param - 720.0*param*param + 120.0*param)*alpha;
    coefs.second=(30.0*p0 - 180.0*p1 + 450.0*p2 + 450.0*p4 - 180.0*p5 + 30.0*p6)*alpha*param*param*param*param + 1.0*(-120.0*p0 + 600.0*p1 - 1200.0*p2 - 600.0*p4 + 120.0*p5)*alpha*param*param*param + 1.0*(180.0*p0 - 720.0*p1 + 1080.0*p2 + 180.0*p4)*alpha*param*param + 1.0*(-120.0*p0 + 360.0*p1 - 360.0*p2)*alpha*param + 1.0*(30.0*p0 - 60.0*p1 + 30.0*p2)*alpha;
    return coefs;
}

template <typename Path>
void computeDistanceCostFunction(const std::vector<coefs_t>& cks , const Path& path, MatrixXX& H,VectorX& g){
    H = MatrixXX::Zero(DIM_POINT,DIM_POINT);
    g  = VectorX::Zero(DIM_POINT);
    int numPoints = cks.size();
    double step = 1./(numPoints-1);
  //  std::cout<<"compute cost; step = "<<step<<std::endl;
    int i = 0;
    point3_t pk;
    for (std::vector<coefs_t>::const_iterator ckcit = cks.begin(); ckcit != cks.end(); ++ckcit){
        pk=path(i*step);
   //     std::cout<<"path ( "<<i*step<<" ) = "<<pk.transpose()<<std::endl;
        H += (ckcit->first * ckcit->first * Matrix3::Identity());
        g += (ckcit->first* (2*ckcit->second - 2*pk ) );
        i++;
    }
    //normalize :
   /* double norm=H(0,0); // because H is always diagonal
    H /= norm;
    g /= sqrt(norm);*/
}

void computeC_of_T (const ProblemData& pData,double T, ResultDataCOMTraj& res){
    std::vector<Vector3> wps;
    point_t p0,p1,p2,p4,p5,p6;
    computeConstantWaypoints(pData,T,6,p0,p1,p2,p4,p5,p6);
    wps.push_back(p0);
    wps.push_back(p1);
    wps.push_back(p2);
    wps.push_back(res.x);
    wps.push_back(p4);
    wps.push_back(p5);
    wps.push_back(p6);
    res.c_of_t_ = bezier_t (wps.begin(), wps.end(),T);
    std::cout<<"bezier curve created, size = "<<res.c_of_t_.size_<<std::endl;
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
    //normalize :
    /*double norm=H(0,0); // because H is always diagonal
    H /= norm;
    g /= sqrt(norm);*/
   // H(1,1)=1e-6;
   // g(1)=0;
}


template <typename Path>
ResultDataCOMTraj solveEndEffector(const ProblemData& pData,const Path& path, const double T, const double weightDistance){
    std::cout<<"solve end effector, T = "<<T<<std::endl;
    assert (weightDistance>=0. && weightDistance<=1. && "WeightDistance must be between 0 and 1");
    double weightAcc = 1. - weightDistance;
    std::vector<waypoint_t> wps_acc=createEndEffectorAccelerationWaypoints(T,pData);
    std::vector<waypoint_t> wps_vel=createEndEffectorVelocityWaypoints(T,pData);
    // stack the constraint for each waypoint :
    MatrixXX A;
    VectorX b;
    Vector3 acc_bounds(10,10,10);
    Vector3 vel_bounds(5,5,5);
    computeConstraintsMatrix(wps_acc,wps_vel,acc_bounds,vel_bounds,A,b);
 //   std::cout<<"End eff A = "<<std::endl<<A<<std::endl;
 //   std::cout<<"End eff b = "<<std::endl<<b<<std::endl;
    // compute cost function (discrete integral under the curve defined by 'path')
    std::vector<coefs_t> cks = createDiscretizationPoints(pData);
    MatrixXX H_rrt,H_acc,H;
    VectorX g_rrt,g_acc,g;
    computeDistanceCostFunction<Path>(cks,path,H_rrt,g_rrt);
    computeAccelerationCostFunction(50,pData,T,H_acc,g_acc);
  /*  std::cout<<"End eff H_rrt = "<<std::endl<<H_rrt<<std::endl;
    std::cout<<"End eff g_rrt = "<<std::endl<<g_rrt<<std::endl;
    std::cout<<"End eff H_acc = "<<std::endl<<H_acc<<std::endl;
    std::cout<<"End eff g_acc = "<<std::endl<<g_acc<<std::endl;
*/
    // add the costs :
    H = MatrixXX::Zero(DIM_POINT,DIM_POINT);
    g  = VectorX::Zero(DIM_POINT);
    H = weightAcc*H_acc + weightDistance*H_rrt;
    g = weightAcc*g_acc + weightDistance*g_rrt;
    std::cout<<"End eff H = "<<std::endl<<H<<std::endl;
    std::cout<<"End eff g = "<<std::endl<<g<<std::endl;

    // call the solver
    VectorX init = VectorX(DIM_POINT);
    init = (pData.c0_ + pData.c1_)/2.;
   // init =pData.c0_;
    std::cout<<"Init = "<<std::endl<<init<<std::endl;
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
   std::cout<<"Solved, success = "<<res.success_<<" x = "<<res.x.transpose()<<std::endl;
   std::cout<<"Final cost : "<<resQp.cost_<<std::endl;
   return res;
}


}//namespace bezier_com_traj
