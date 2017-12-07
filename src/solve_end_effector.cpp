/*
 * Copyright 2017, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>
using namespace bezier_com_traj;
typedef waypoint3_t waypoint_t;
typedef std::pair<double,point3_t> coefs_t;
const int DIM_POINT=3;
const int NUM_DISCRETIZATION = 11;
namespace bezier_com_traj
{

coefs_t initCoefs(){
    coefs_t c;
    c.first=0;
    c.second=point3_t::Zero();
    return c;
}

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

void computeConstraintsMatrix(std::vector<waypoint_t> wps,MatrixXX& A,VectorX& b){
    A = MatrixXX::Zero(DIM_POINT*wps.size(),DIM_POINT);
    b = VectorX::Zero(DIM_POINT);
    int i = 0;
    for (std::vector<waypoint_t>::const_iterator wpcit = wps.begin(); wpcit != wps.end(); ++wpcit)
    {
        A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = wpcit->first;
        b.segment<DIM_POINT>(i*DIM_POINT)   = -wpcit->second;
        ++i;
    }
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

template <typename Path>
void computeCostFunction(std::vector<coefs_t>cks ,MatrixXX& H,VectorX& g,Path path){
    H = MatrixXX::Zero(DIM_POINT,DIM_POINT);
    g  = VectorX::Zero(DIM_POINT);
    int numPoints = cks.size();
    double step = 1./(numPoints-1);
    std::cout<<"compute cost; step = "<<step;
    int i = 0;
    point3_t pk;
    for (std::vector<coefs_t>::const_iterator ckcit = cks.begin(); ckcit != cks.end(); ++ckcit){
        pk=path(i*step);
        std::cout<<"path ( "<<i*step<<" ) = "<<pk<<std::endl;
        H += (ckcit->first * ckcit->first * Matrix3::Identity());
        g += (ckcit->first* (2*ckcit->second - 2*pk ) );
    }

}


template <typename Path>
ResultDataCOMTraj solveEndEffector(const ProblemData& pData,Path path, const double T, const double timeStep){
    std::vector<waypoint_t> wps=createEndEffectorWaypoints(T,pData);
    // stack the constraint for each waypoint :
    MatrixXX A;
    VectorX b;
    computeConstraintsMatrix(wps,A,b);
    std::cout<<"End eff A = "<<A<<std::endl;
    std::cout<<"End eff b = "<<b<<std::endl;
    // compute cost function (discrete integral under the curve defined by 'path')
    std::vector<coefs_t> cks = createDiscretizationPoints(pData);
    MatrixXX H;
    VectorX g;
    computeCostFunction<Path>(cks,path,H,g);
    std::cout<<"End eff H = "<<H<<std::endl;
    std::cout<<"End eff g = "<<g<<std::endl;


    // call the solver
    VectorX init = VectorX(DIM_POINT);
    init = pData.c0_;
    ResultData resQp = solve(A,b,H,g, init);


}


}//namespace bezier_com_traj
