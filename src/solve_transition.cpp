/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>
#include  <limits>


#ifndef QHULL
#define QHULL 0
#endif


namespace bezier_com_traj
{
typedef waypoint3_t waypoint_t;
typedef std::pair<double,point3_t> coefs_t;

ResultData solveIntersection(const std::pair<MatrixXX, VectorX>& Ab,const std::pair<MatrixXX, VectorX>& Hg,  const Vector3& init)
{
    return solve(Ab.first,Ab.second,Hg.first,Hg.second, init);
}

void printQHullFile(const std::pair<MatrixXX, VectorX>& Ab,VectorX intPoint,const std::string& fileName,bool clipZ = false){
     std::ofstream file;
     using std::endl;
     std::string path("/home/pfernbac/Documents/com_ineq_test/");
     path.append(fileName);
     file.open(path.c_str(),std::ios::out | std::ios::trunc);
     file<<"3 1"<<endl;
     file<<"\t "<<intPoint[0]<<"\t"<<intPoint[1]<<"\t"<<intPoint[2]<<endl;
     file<<"4"<<endl;
     clipZ ? file<<Ab.first.rows()+2<<endl : file<<Ab.first.rows()<<endl;
     for(size_t i = 0 ; i < Ab.first.rows() ; ++i){
         file<<"\t"<<Ab.first(i,0)<<"\t"<<Ab.first(i,1)<<"\t"<<Ab.first(i,2)<<"\t"<<-Ab.second[i]-0.001<<endl;
     }
     if(clipZ){
         file<<"\t"<<0<<"\t"<<0<<"\t"<<1.<<"\t"<<-3.<<endl;
         file<<"\t"<<0<<"\t"<<0<<"\t"<<-1.<<"\t"<<-1.<<endl;
     }
     file.close();
}


/// ### EQUATION FOR CONSTRAINTS ON INIT AND FINAL POSITION AND VELOCITY (DEGREE = 4)
/** @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve, assume p0 p1 x p2 p3
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
coefs_t evaluateCurveAtTime(std::vector<point_t> pi,double t){
    coefs_t wp;
    double t2 = t*t;
    double t3 = t2*t;
    double t4 = t3*t;
    // equation found with sympy
    wp.first = ( 6.0*t4 - 12.0*t3 + 6.0*t2);
    wp.second = 1.0*pi[0]*t4 - 4.0*pi[0]*t3 + 6.0*pi[0]*t2 - 4.0*pi[0]*t + 1.0*pi[0] - 4.0*pi[1]*t4 + 12.0*pi[1]*t3 - 12.0*pi[1]*t2 + 4.0*pi[1]*t - 4.0*pi[3]*t4 + 4.0*pi[3]*t3 + 1.0*pi[4]*t4;
   // std::cout<<"wp at t = "<<t<<std::endl;
   // std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
    return wp;
}

coefs_t evaluateAccelerationCurveAtTime(std::vector<point_t> pi,double T,double t){
    coefs_t wp;
    double alpha = 1./(T*T);
    // equation found with sympy
    wp.first = (72.0*t*t - 72.0*t + 12.0)*alpha;
    wp.second = (12.0*pi[0]*t*t - 24.0*pi[0]*t + 12.0*pi[0] - 48.0*pi[1]*t*t + 72.0*pi[1]*t - 24.0*pi[1] - 48.0*pi[3]*t*t + 24.0*pi[3]*t + 12.0*pi[4]*t*t)*alpha;
   // std::cout<<"acc_wp at t = "<<t<<std::endl;
   // std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
    return wp;
}


std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
    // equation for constraint on initial and final position and velocity (degree 4, 4 constant waypoint and one free (p2))
    // first, compute the constant waypoints that only depend on pData :
    int n = 4;
    std::vector<point_t> pi;
    pi.push_back(pData.c0_); //p0
    pi.push_back((pData.dc0_ * T / n )+  pData.c0_); // p1
    pi.push_back(point_t::Zero()); // p2 = x
    pi.push_back((-pData.dc1_ * T / n) + pData.c1_); // p3
    pi.push_back(pData.c1_); // p4
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
    w0.first.block<3,3>(0,0) = 12.*alpha*Matrix3::Identity();
    w0.first.block<3,3>(3,0) = 12.*alpha*Cpi[0];
    w0.second.head<3>() = (12.*pi[0] - 24.*pi[1])*alpha;
    w0.second.tail<3>() = 1.0*Cg*pi[0] - (24.0*Cpi[0]*pi[1])*alpha;
    wps.push_back(w0);
    waypoint6_t w1 = initwp<waypoint6_t>();
    w1.first.block<3,3>(0,0) =  -2.4*alpha*Matrix3::Identity();
    w1.first.block<3,3>(3,0) =(-12.0*Cpi[0] + 9.6*Cpi[1])*alpha;
    w1.second.head<3>() = (7.2*pi[0] - 9.6*pi[1] + 4.8*pi[3])*alpha;
    w1.second.tail<3>() = (0.2*Cg*T2*pi[0] + 0.8*Cg*T2*pi[1] + 4.8*Cpi[0]*pi[3])*alpha;
    wps.push_back(w1);
    waypoint6_t w2 = initwp<waypoint6_t>();
    w2.first.block<3,3>(0,0) =  -9.6*alpha*Matrix3::Identity();
    w2.first.block<3,3>(3,0) = (0.6*Cg*T2 - 9.6*Cpi[1])*alpha;
    w2.second.head<3>() = (3.6*pi[0]  + 4.8*pi[3] + 1.2*pi[4])*alpha;
    w2.second.tail<3>() = (0.4*Cg*T2*pi[1] - 4.8*Cpi[0]*pi[3] + 1.2*Cpi[0]*pi[4] + 9.6*Cpi[1]*pi[3])*alpha;
    wps.push_back(w2);
    waypoint6_t w3 = initwp<waypoint6_t>();
    w3.first.block<3,3>(0,0) = -9.6*alpha*Matrix3::Identity();
    w3.first.block<3,3>(3,0) = (0.6*Cg*T2  - 9.6*Cpi[3])*alpha;
    w3.second.head<3>() = (1.2*pi[0] + 4.8*pi[1]  + 3.6*pi[4])*alpha;
    w3.second.tail<3>() = (0.4*Cg*T2*pi[3]  - 1.2*Cpi[0]*pi[4] - 9.6*Cpi[1]*pi[3] + 4.8*Cpi[1]*pi[4])*alpha;
    wps.push_back(w3);
    waypoint6_t w4 = initwp<waypoint6_t>();
    w4.first.block<3,3>(0,0) = -2.4*alpha*Matrix3::Identity();
    w4.first.block<3,3>(3,0) =(9.6*Cpi[3] - 12.0*Cpi[4])*alpha;
    w4.second.head<3>() = (4.8*pi[1] - 9.6*pi[3] + 7.2*pi[4])*alpha;
    w4.second.tail<3>() = (0.8*Cg*T2*pi[3] + 0.2*Cg*T2*pi[4] - 4.8*Cpi[1]*pi[4])*alpha;
    wps.push_back(w4);
    waypoint6_t w5 = initwp<waypoint6_t>();
    w5.first.block<3,3>(0,0) =12*alpha*Matrix3::Identity();
    w5.first.block<3,3>(3,0) =12.0*Cpi[4]*alpha;
    w5.second.head<3>() = (-24*pi[3] + 12*pi[4])*alpha;
    w5.second.tail<3>() =(Cg*T2*pi[4]  + 24.0*Cpi[3]*pi[4])*alpha;
    wps.push_back(w5);
    return wps;
}



/// ### EQUATION FOR CONSTRAINTS ON INIT AND FINAL POSITION AND VELOCITY AND INIT ACCELERATION (DEGREE = 5)
///

///**
// * @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
// * @param pi constant waypoints of the curve, assume p0 p1 x p2 p3
// * @param t param (normalized !)
// * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
// */
//coefs_t evaluateCurveAtTime(std::vector<point_t> pi,double t){
//    coefs_t wp;
//    double t2 = t*t;
//    double t3 = t2*t;
//    double t4 = t3*t;
//    double t5 = t4*t;
//    // equation found with sympy
//    wp.first = 10.0*t5 - 20.0*t4 + 10.0*t3;
//    wp.second = -1.0*pi[0]*t5 + 5.0*pi[0]*t4 - 10.0*pi[0]*t3 + 10.0*pi[0]*t2 - 5.0*pi[0]*t + 1.0*pi[0] + 5.0*pi[1]*t5 - 20.0*pi[1]*t4 + 30.0*pi[1]*t3 - 20.0*pi[1]*t2 + 5.0*pi[1]*t - 10.0*pi[2]*t5 + 30.0*pi[2]*t4 - 30.0*pi[2]*t3 + 10.0*pi[2]*t2 - 5.0*pi[4]*t5 + 5.0*pi[4]*t4 + 1.0*pi[5]*t5;
//   // std::cout<<"wp at t = "<<t<<std::endl;
//   // std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
//    return wp;
//}

//coefs_t evaluateAccelerationCurveAtTime(std::vector<point_t> pi,double T,double t){
//    coefs_t wp;
//    double alpha = 1./(T*T);
//    double t2 = t*t;
//    double t3 = t2*t;
//    // equation found with sympy
//    wp.first = (200.0*t3 - 240.0*t2 + 60.0*t)*alpha;
//    wp.second = 1.0*(-20.0*pi[0]*t3 + 60.0*pi[0]*t2 - 60.0*pi[0]*t + 20.0*pi[0] + 100.0*pi[1]*t3 - 240.0*pi[1]*t2 + 180.0*pi[1]*t - 40.0*pi[1] - 200.0*pi[2]*t3 + 360.0*pi[2]*t2 - 180.0*pi[2]*t + 20.0*pi[2] - 100.0*pi[4]*t3 + 60.0*pi[4]*t2 + 20.0*pi[5]*t3)*alpha;
//   // std::cout<<"acc_wp at t = "<<t<<std::endl;
//   // std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
//    return wp;
//}


//std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
//    // equation for constraint on initial and final position and velocity and initial acceleration(degree 5, 5 constant waypoint and one free (p3))
//    // first, compute the constant waypoints that only depend on pData :
//    double n = 5.;
//    std::vector<point_t> pi;
//    pi.push_back(pData.c0_); //p0
//    pi.push_back((pData.dc0_ * T / n )+  pData.c0_); // p1
//    pi.push_back((pData.ddc0_*T*T/(n*(n-1))) + (2.*pData.dc0_ *T / n) + pData.c0_); // p2
//    pi.push_back(point_t::Zero()); // x
//    pi.push_back((-pData.dc1_ * T / n) + pData.c1_); // p4
//    pi.push_back(pData.c1_); // p5
//    /*std::cout<<"fixed waypoints : "<<std::endl;
//    for(std::vector<point_t>::const_iterator pit = pi.begin() ; pit != pi.end() ; ++pit){
//        std::cout<<" pi = "<<*pit<<std::endl;
//    }*/
//    return pi;
//}

/// ### EQUATION FOR CONSTRAINts on initial position, velocity and acceleration, and only final position (degree = 4)
///
///**
// * @brief evaluateCurveAtTime compute the expression of the point on the curve at t, defined by the waypoint pi and one free waypoint (x)
// * @param pi constant waypoints of the curve, assume p0 p1 x p2 p3
// * @param t param (normalized !)
// * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
// */
//coefs_t evaluateCurveAtTime(std::vector<point_t> pi,double t){
//    coefs_t wp;
//    double t2 = t*t;
//    double t3 = t2*t;
//    double t4 = t3*t;
//    // equation found with sympy
//    wp.first = -4.0*t4 + 4.0*t3;
//    wp.second =1.0*pi[0]*t4 - 4.0*pi[0]*t3 + 6.0*pi[0]*t2 - 4.0*pi[0]*t + 1.0*pi[0] - 4.0*pi[1]*t4 + 12.0*pi[1]*t3 - 12.0*pi[1]*t2 + 4.0*pi[1]*t + 6.0*pi[2]*t4 - 12.0*pi[2]*t3 + 6.0*pi[2]*t2 + 1.0*pi[4]*t4;
//    std::cout<<"wp at t = "<<t<<std::endl;
//    std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
//    return wp;
//}

//coefs_t evaluateAccelerationCurveAtTime(std::vector<point_t> pi,double T,double t){
//    coefs_t wp;
//    double alpha = 1./(T*T);
//    double t2 = t*t;
//    // equation found with sympy
//    wp.first = (-48.0*t2 + 24.0*t)*alpha;
//    wp.second = (12.0*pi[0]*t2 - 24.0*pi[0]*t + 12.0*pi[0] - 48.0*pi[1]*t2 + 72.0*pi[1]*t - 24.0*pi[1] + 72.0*pi[2]*t2 - 72.0*pi[2]*t + 12.0*pi[2] + 12.0*pi[4]*t2)*alpha;
//    std::cout<<"acc_wp at t = "<<t<<std::endl;
//    std::cout<<" first : "<<wp.first<<" ; second : "<<wp.second.transpose()<<std::endl;
//    return wp;
//}


//std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T){
//    // equation for constraint on initial position, velocity and acceleration, and only final position (degree = 4)(degree 4, 4 constant waypoint and one free (p3))
//    // first, compute the constant waypoints that only depend on pData :
//    double n = 4.;
//    std::vector<point_t> pi;
//    pi.push_back(pData.c0_); //p0
//    pi.push_back((pData.dc0_ * T / n )+  pData.c0_); // p1
//    pi.push_back((pData.ddc0_*T*T/(n*(n-1))) + (2.*pData.dc0_ *T / n) + pData.c0_); // p2
//    pi.push_back(point_t::Zero()); // x
//    pi.push_back(pData.c1_); // p4
//    std::cout<<"fixed waypoints : "<<std::endl;
//    for(std::vector<point_t>::const_iterator pit = pi.begin() ; pit != pi.end() ; ++pit){
//        std::cout<<" pi = "<<*pit<<std::endl;
//    }
//    return pi;
//}


//coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T){
//     coefs_t v4;
//     // equation found with sympy
//     v4.first = -4./T;
//     v4.second = 4.* pData.c1_ / T;
//     return v4;
//}

//void computeFinalVelocity(const ProblemData& pData,double T,ResultDataCOMTraj& res){
//    coefs_t v = computeFinalVelocityPoint(pData,T);
//    res.dc1_ = v.first*res.x + v.second;
//}


//void computeFinalAcceleration(const ProblemData& pData,double T,ResultDataCOMTraj& res){
//    point_t p2,p4;
//    p2 = (pData.ddc0_*T*T/(4*(3))) + (2.*pData.dc0_ *T / 4) + pData.c0_;
//    p4 = pData.c1_;
//    coefs_t a;
//    a.first = -24/(T*T);
//    a.second = 12*(p2 + p4)/(T*T);

//    res.ddc1_ = a.first * res.x + a.second;
//}

/**
 * @brief computeDiscretizedTime build an array of discretized points in time, such that there is the same number of point in each phase. Doesn't contain t=0, is of size pointsPerPhase*phaseTimings.size()
 * @param phaseTimings
 * @param pointsPerPhase
 * @return
 */
std::vector<double> computeDiscretizedTime(const VectorX& phaseTimings,const int pointsPerPhase ){
    std::vector<double> timeArray;
    double t = 0;
    for(int i = 0 ; i < phaseTimings.size() ; ++i){
        double step = (double) phaseTimings[i] / pointsPerPhase;
        for(int j = 0 ; j < pointsPerPhase ; ++j){
            t += step;
            timeArray.push_back(t);
        }
    }
    return timeArray;
}

std::vector<coefs_t> computeDiscretizedWaypoints(const ProblemData& pData,double T,const std::vector<double>& timeArray){
    //int numStep = int(T / timeStep);
    std::vector<coefs_t> wps;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    // evaluate curve work with normalized time !
    double t;
    for (int i = 0 ; i<timeArray.size() ; ++i ){
        t = timeArray[i] / T;
        if(t>1)
            t=1.;
        wps.push_back(evaluateCurveAtTime(pi,t));
    }
    return wps;
}


std::vector<waypoint6_t> computeDiscretizedWwaypoints(const ProblemData& pData,double T,const std::vector<double>& timeArray){
    std::vector<waypoint6_t> wps = computeWwaypoints(pData,T);
    std::vector<waypoint6_t> res;
    std::vector<spline::Bern<double> > berns = ComputeBersteinPolynoms(wps.size()-1);

    for(int i = 0 ; i < timeArray.size() ; ++i){
        waypoint6_t w = initwp<waypoint6_t>();
        for (int j = 0 ; j < wps.size() ; ++j){
            double b = berns[j](timeArray[i]);
            w.first +=b*(wps[j].first );
            w.second+=b*(wps[j].second);
        }
        res.push_back(w);
    }
    return res;
}



std::vector<coefs_t> computeDiscretizedAccelerationWaypoints(const ProblemData& pData,double T,const std::vector<double>& timeArray){
    //int numStep = int(T / timeStep);
    std::vector<coefs_t> wps;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    double t;
    for (int i = 0 ; i<timeArray.size() ; ++i ){
        t = timeArray[i] / T;
        if(t>1)
            t=1.;
        wps.push_back(evaluateAccelerationCurveAtTime(pi,T,t));
    }
    return wps;
}

 std::pair<MatrixXX,VectorX> dynamicStabilityConstraints_cross(const MatrixXX& mH,const VectorX& h,const Vector3& g,const coefs_t& c,const coefs_t& ddc){
     Matrix3 S_hat;
     int dimH = (int)(mH.rows());
     MatrixXX A(dimH,4);
     VectorX b(dimH);
     S_hat = skew(c.second*ddc.first - ddc.second*c.first + g*c.first);
     A.block(0,0,dimH,3) = mH.block(0,3,dimH,3) * S_hat + mH.block(0,0,dimH,3) * ddc.first;
     b = h + mH.block(0,0,dimH,3)*(g - ddc.second) + mH.block(0,3,dimH,3)*(c.second.cross(g) - c.second.cross(ddc.second));
    // Normalize(A,b);
     // add 1 for the slack variable :
     A.block(0,3,dimH,1) = VectorX::Ones(dimH);
     return std::make_pair<MatrixXX,VectorX>(A,b);
}

std::pair<MatrixXX,VectorX> dynamicStabilityConstraints(const MatrixXX& mH,const VectorX& h,const Vector3& g,const waypoint6_t& w){
    int dimH = (int)(mH.rows());
    MatrixXX A(dimH,4);
    VectorX b(dimH);
    VectorX g_(6);
    g_.head<3>() = g;
    A.block(0,0,dimH,3) = mH*w.first;
    b = h + mH*(g_ - w.second);
    //Normalize(A,b);
    A.block(0,3,dimH,1) = VectorX::Ones(dimH);
    return std::make_pair<MatrixXX,VectorX>(A,b);
}

std::pair<MatrixXX,VectorX> staticStabilityConstraints(const MatrixXX& mH,const VectorX& h, const Vector3& g,const coefs_t& c){
     int dimH = (int)(mH.rows());
     MatrixXX A(dimH,4);
     VectorX b(dimH);
     A.block(0,0,dimH,3) = mH.block(0,3,dimH,3) * c.first * skew(g);
     b = h + mH.block(0,0,dimH,3)*g - mH.block(0,3,dimH,3)*g.cross(c.second);
     // add 1 for the slack variable :
     A.block(0,3,dimH,1) = VectorX::Ones(dimH);
    // Normalize(A,b);
     return std::make_pair<MatrixXX,VectorX>(A,b);
}

void compareStabilityMethods(const MatrixXX& mH,const VectorX& h,const Vector3& g,const coefs_t& c,const coefs_t& ddc,const waypoint6_t& w){
    std::pair<MatrixXX,VectorX> Ab_cross,Ab_w;
    /*
    Vector3 wd;
    wd = c.cross(ddc) + g.cross(c);
    std::cout<<"wu cross : "<<ddc.first<<std::endl;
    std::cout<<"wu       : "<<w.first.block<3,3>(0,0)<<std::endl;
    error_wu = ddc.first - w.first(0,0);
    */

    Ab_cross = dynamicStabilityConstraints_cross(mH,h,g,c,ddc);
    Ab_w = dynamicStabilityConstraints(mH,h,g,w);
    Normalize(Ab_cross.first,Ab_cross.second);
    Normalize(Ab_w.first,Ab_w.second);


    MatrixXX A_error = Ab_cross.first - Ab_w.first;
    VectorX b_error = Ab_cross.second - Ab_w.second;
    double A_error_norm = A_error.lpNorm<Eigen::Infinity>();
    double b_error_norm = b_error.lpNorm<Eigen::Infinity>();
    std::cout<<" max a error : "<<A_error_norm<<" ; b : "<<b_error_norm<<std::endl;
    std::cout<<"A error : "<<std::endl<<A_error<<std::endl;
    std::cout<<"b error : "<<std::endl<<b_error<<std::endl;

    assert(A_error_norm < 1e-4 && b_error_norm < 1e-4 && "Both method didn't find the same results.");
}


std::pair<MatrixXX, VectorX> computeConstraintsOneStep(const ProblemData& pData,const VectorX& Ts,const int pointsPerPhase,VectorX& constraints_equivalence){
    // compute the list of discretized waypoint :
    double t_total = 0.;
    for(int i = 0 ; i < Ts.size() ; ++i)
        t_total+=Ts[i];
    // Compute all the discretized wayPoint
    std::cout<<"total time : "<<t_total<<std::endl;
    std::vector<double> timeArray = computeDiscretizedTime(Ts,pointsPerPhase);
    std::vector<coefs_t> wps_c = computeDiscretizedWaypoints(pData,t_total,timeArray);
    //std::vector<coefs_t> wps_ddc = computeDiscretizedAccelerationWaypoints(pData,t_total,timeArray);
    std::vector<waypoint6_t> wps_w = computeDiscretizedWwaypoints(pData,t_total,timeArray);
    std::cout<<" number of discretized waypoints c: "<<wps_c.size()<<std::endl;
    std::cout<<" number of discretized waypoints w: "<<wps_w.size()<<std::endl;
    assert(/*wps_c.size() == wps_ddc.size() &&*/  wps_w.size() == wps_c.size());
    std::vector<int> stepIdForPhase; // stepIdForPhase[i] is the id of the last step of phase i / first step of phase i+1 (overlap)
    for(int i = 0 ; i < Ts.size() ; ++i)
        stepIdForPhase.push_back(pointsPerPhase*(i+1)-1);

    assert(stepIdForPhase.back() == (wps_c.size()-1)); // -1 because the first one is the index (start at 0) and the second is the size
    // compute the total number of inequalities (to initialise A and b)
    int num_ineq = 0;
    int num_stab_ineq = 0;
    int num_kin_ineq = 0;
    int numStepForPhase;
    centroidal_dynamics::MatrixXX Hrow;
    VectorX h;
    MatrixXX H,mH;
    for(int i = 0 ; i < Ts.size() ; ++i){
        pData.contacts_[i].contactPhase_->getPolytopeInequalities(Hrow,h);
        numStepForPhase = pointsPerPhase;
        if(i > 0 )
            ++numStepForPhase; // because at the switch point between phases we add the constraints of both phases.
        std::cout<<"constraint size : Kin = "<<pData.contacts_[i].kin_.rows()<<" ; stab : "<<Hrow.rows()<<" times "<<numStepForPhase<<" steps"<<std::endl;
        num_stab_ineq += Hrow.rows() * numStepForPhase;
        if(i == Ts.size()-1)
            --numStepForPhase; // we don't consider kinematics constraints for the last point (because it's independant of x)
        num_kin_ineq += pData.contacts_[i].kin_.rows() * numStepForPhase;
    }
    num_ineq = num_stab_ineq + num_kin_ineq;
    std::cout<<"total of inequalities : "<<num_ineq<<std::endl;
    // init constraints matrix :
    MatrixXX A = MatrixXX::Zero(num_ineq,4); // 3 + 1 :  because of the slack constraints
    VectorX b(num_ineq);
    std::pair<MatrixXX,VectorX> Ab_stab;

    int id_rows = 0;
    int current_size;

    int id_phase = 0;
    ContactData phase = pData.contacts_[id_phase];
    // compute some constant matrice for the current phase :
    const Vector3& g = phase.contactPhase_->m_gravity;
    //std::cout<<"g = "<<g.transpose()<<std::endl;
    //std::cout<<"mass = "<<phase.contactPhase_->m_mass<<std::endl;
    //const Matrix3 gSkew = bezier_com_traj::skew(g);
    phase.contactPhase_->getPolytopeInequalities(Hrow,h);
    H = -Hrow;
    H.rowwise().normalize();
    int dimH = (int)(H.rows());
    mH = phase.contactPhase_->m_mass * H;

    // assign the Stability constraints  for each discretized waypoints :
    // we don't consider the first point, because it's independant of x.
    for(int id_step = 0 ; id_step <  timeArray.size() ; ++id_step ){
        // add stability constraints :

        Ab_stab = dynamicStabilityConstraints(mH,h,g,wps_w[id_step]);
        //compareStabilityMethods(mH,h,g,wps_c[id_step],wps_ddc[id_step],wps_w[id_step]);
        A.block(id_rows,0,dimH,4) = Ab_stab.first;
        b.segment(id_rows,dimH) = Ab_stab.second;
        id_rows += dimH ;

        // check if we are going to switch phases :
        for(int i = 0 ; i < (stepIdForPhase.size()-1) ; ++i){
            if(id_step == stepIdForPhase[i]){
                // switch to phase i
                id_phase=i+1;
                phase = pData.contacts_[id_phase];
                phase.contactPhase_->getPolytopeInequalities(Hrow,h);
                H = -Hrow;
                H.rowwise().normalize();
                dimH = (int)(H.rows());
                mH = phase.contactPhase_->m_mass * H;
                // the current waypoint must have the constraints of both phases. So we add it again :
                // TODO : filter for redunbdant constraints ...
                // add stability constraints :
                Ab_stab = dynamicStabilityConstraints(mH,h,g,wps_w[id_step]);
                //compareStabilityMethods(mH,h,g,wps_c[id_step],wps_ddc[id_step],wps_w[id_step]);
                A.block(id_rows,0,dimH,4) = Ab_stab.first;
                b.segment(id_rows,dimH) = Ab_stab.second;
                id_rows += dimH ;
            }
        }
    }

    id_phase = 0;
    phase = pData.contacts_[id_phase];
    // assign the Kinematics constraints  for each discretized waypoints :
    // we don't consider the first point, because it's independant of x.
    for(int id_step = 0 ; id_step <  timeArray.size() ; ++id_step ){
        // add constraints for wp id_step, on current phase :
        // add kinematics constraints :
        // constraint are of the shape A c <= b . But here c(x) = Fx + s so : AFx <= b - As
        if(id_step != timeArray.size()-1){ // we don't consider kinematics constraints for the last point (because it's independant of x)
            current_size = phase.kin_.rows();
            A.block(id_rows,0,current_size,3) = (phase.Kin_ * wps_c[id_step].first);
            b.segment(id_rows,current_size) = phase.kin_ - (phase.Kin_*wps_c[id_step].second);
            id_rows += current_size;
        }

        // check if we are going to switch phases :
        for(int i = 0 ; i < (stepIdForPhase.size()-1) ; ++i){
            if(id_step == stepIdForPhase[i]){
                // switch to phase i
                id_phase=i+1;
                phase = pData.contacts_[id_phase];
                // the current waypoint must have the constraints of both phases. So we add it again :
                // TODO : filter for redunbdant constraints ...
                current_size = phase.kin_.rows();
                A.block(id_rows,0,current_size,3) = (phase.Kin_ * wps_c[id_step].first);
                b.segment(id_rows,current_size) = phase.kin_ - (phase.Kin_*wps_c[id_step].second);
                id_rows += current_size;
            }
        }
    }


    std::cout<<"id rows : "<<id_rows<<" ; total rows : "<<A.rows()<<std::endl;
    assert(id_rows == (A.rows()) && "The constraints matrices were not fully filled.");
    return std::make_pair(A,b);
}

//cost : min distance between x and midPoint :
std::pair<MatrixXX, VectorX> computeCostMidPoint(const ProblemData& pData, size_t size){
    Vector3 midPoint = (pData.c0_ + pData.c1_)/2.; // todo : replace it with point found by planning ??
    MatrixXX H = MatrixXX::Identity(size,size);
    VectorX g = VectorX::Zero(size);
    g.head<3>() = -midPoint;
    return std::make_pair(H,g);
}

//cost : min distance between end velocity and the one computed by planning
//std::pair<MatrixX3, VectorX> computeCostEndVelocity(const ProblemData& pData,const double T){
//    coefs_t v = computeFinalVelocityPoint(pData,T);
//    Matrix3 H = Matrix3::Identity() * v.first * v.first;
//    Vector3 g = v.first*(v.second - pData.dc1_);
//    return std::make_pair(H,g);
//}


void computeBezierCurve(const ProblemData& pData, const double T, ResultDataCOMTraj& res)
{
    std::vector<Vector3> wps;

    std::vector<Vector3> pi = computeConstantWaypoints(pData,T);
    wps.push_back(pi[0]);
    wps.push_back(pi[1]);
    //wps.push_back(pi[2]);
    wps.push_back(res.x);
    wps.push_back(pi[3]);
    wps.push_back(pi[4]);
    res.c_of_t_ = bezier_t (wps.begin(), wps.end(),T);
}

double analyseSlack(const VectorX& slack,const VectorX& constraint_equivalence ){
    //TODO
    assert(slack.size() == constraint_equivalence.size() && "slack variables and constraints equivalence should have the same size." );
    //std::cout<<"slack : "<<slack<<std::endl;
    std::cout<<"list of violated constraints : "<<std::endl;
    double previous_id = -1;
    for(size_t i = 0 ; i < slack.size() ; ++i){
        if((slack[i]*slack[i]) > std::numeric_limits<double>::epsilon()){
            if(constraint_equivalence[i] != previous_id){
                std::cout<<"step "<<constraint_equivalence[i]<<std::endl;
                previous_id = constraint_equivalence[i];
            }
        }
    }
    //return (slack.squaredNorm())/(slack.size());
    return slack.lpNorm<Eigen::Infinity>();
}

ResultDataCOMTraj solveOnestep(const ProblemData& pData, const VectorX& Ts,const Vector3& init_guess,const int pointsPerPhase, const double feasability_treshold){
    assert(pData.contacts_.size() ==2 || pData.contacts_.size() ==3);
    assert(Ts.size() == pData.contacts_.size());
    double T = 0;
    for(int i = 0 ; i < Ts.size() ; ++i)
        T+=Ts[i];
   // bool fail = true;
    ResultDataCOMTraj res;
    VectorX constraint_equivalence;
    std::pair<MatrixXX, VectorX> Ab = computeConstraintsOneStep(pData,Ts,pointsPerPhase,constraint_equivalence);
   // std::pair<MatrixX3, VectorX> Hg = computeCostEndVelocity(pData,T);
    std::pair<MatrixXX, VectorX> Hg = computeCostMidPoint(pData,constraint_equivalence.size()+3);

    std::cout<<"Init = "<<std::endl<<init_guess.transpose()<<std::endl;
    VectorX x = VectorX::Zero(3 + constraint_equivalence.size());
    x.head<3>() = init_guess;

    // rewriting 0.5 || Dx -d ||^2 as x'Hx  + g'x
    ResultData resQp = solve(Ab.first,Ab.second,Hg.first,Hg.second, x);


    double feasability = analyseSlack(resQp.x.tail(constraint_equivalence.size()),constraint_equivalence);
    std::cout<<"feasability : "<<feasability<<std::endl;

    if(feasability<=feasability_treshold)
    {
        res.success_ = true;
        res.x = resQp.x.head<3>();
        computeBezierCurve (pData,T,res);
//        computeFinalVelocity(pData,T,res);
//        computeFinalAcceleration(pData,T,res);
        std::cout<<"Solved, success "<<" x = ["<<res.x[0]<<","<<res.x[1]<<","<<res.x[2]<<"]"<<std::endl;
    }
    #if QHULL
    printQHullFile(Ab,resQp.x,"bezier_wp.txt");
    #endif
    std::cout<<"Final cost : "<<resQp.cost_<<std::endl;
    return res;
}


} // namespace
