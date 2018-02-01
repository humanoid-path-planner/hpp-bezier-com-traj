/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>

namespace bezier_com_traj
{
typedef waypoint3_t waypoint_t;
typedef std::pair<double,point3_t> coefs_t;

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
coefs_t evaluateCurveAtTime(std::vector<point_t> pi,double t){
    coefs_t wp;
    double t2 = t*t;
    double t3 = t2*t;
    double t4 = t3*t;
    // equation found with sympy
    wp.first = ( 6.0*t4 - 12.0*t3 + 6.0*t2);
    wp.second = 1.0*pi[0]*t4 - 4.0*pi[0]*t3 + 6.0*pi[0]*t2 - 4.0*pi[0]*t + 1.0*pi[0] - 4.0*pi[1]*t4 + 12.0*pi[1]*t3 - 12.0*pi[1]*t2 + 4.0*pi[1]*t - 4.0*pi[2]*t4 + 4.0*pi[2]*t3 + 1.0*pi[3]*t4;
    return wp;
}

coefs_t evaluateAccelerationCurveAtTime(std::vector<point_t> pi,double T,double t){
    coefs_t wp;
    double alpha = 1./(T*T);
    // equation found with sympy
    wp.first = (72.0*t*t - 72.0*t + 12.0)*alpha;
    wp.second = (12.0*pi[0]*t*t - 24.0*pi[0]*t + 12.0*pi[0] - 48.0*pi[1]*t*t + 72.0*pi[1]*t - 24.0*pi[1] - 48.0*pi[2]*t*t + 24.0*pi[2]*t + 12.0*pi[3]*t*t)*alpha;
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

std::vector<coefs_t> computeDiscretizedWaypoints(const ProblemData& pData,double T,double timeStep){
    //int numStep = int(T / timeStep);
    std::vector<coefs_t> wps;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    double t = 0;
    while(t<T){
        wps.push_back(evaluateCurveAtTime(pi,t));
        t+= timeStep;
    }
    return wps;
}


std::vector<coefs_t> computeDiscretizedAccelerationWaypoints(const ProblemData& pData,double T,double timeStep){
    //int numStep = int(T / timeStep);
    std::vector<coefs_t> wps;
    std::vector<point_t> pi = computeConstantWaypoints(pData,T);
    double t = 0;
    while(t<T){
        wps.push_back(evaluateAccelerationCurveAtTime(pi,T,t));
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
    std::vector<coefs_t> wps = computeDiscretizedWaypoints(pData,t_total,timeStep);
    std::vector<coefs_t> acc_wps = computeDiscretizedAccelerationWaypoints(pData,t_total,timeStep);
    int numStep = int(t_total / timeStep);
    assert(wps.size() == acc_wps.size());
    assert(numStep == wps.size());
    std::vector<int> stepPerPhases;
    for(int i = 0 ; i < Ts.size() ; ++i){
        if(i == 0)
            stepPerPhases.push_back(int(Ts[i] / timeStep));
        else
            stepPerPhases.push_back((int(Ts[i] / timeStep))+stepPerPhases.back());

    }
    // compute the total number of inequalities (to initialise A and b)
    int num_ineq = 0;
    int numStepForPhase;
    centroidal_dynamics::MatrixXX Hrow;
    VectorX h;
    MatrixX3 H,mH;
    for(int i = 0 ; i < Ts.size() ; ++i){
        if(i==0)
            numStepForPhase = stepPerPhases[i];
        else
            numStepForPhase = stepPerPhases[i] -stepPerPhases[i-1];
        num_ineq += pData.contacts_[i].kin_.rows() * numStepForPhase;
        pData.contacts_[i].contactPhase_->getPolytopeInequalities(Hrow,h);
        num_ineq += Hrow.rows() * numStepForPhase;
    }
    std::cout<<"total of inequalities : "<<num_ineq<<std::endl;

    // assign the constraints (kinematics and stability for each waypoints :
    MatrixX3 A(num_ineq,3);
    VectorX b(num_ineq);
    Matrix3 S_hat;

    int id_rows = 0;
    double t = 0;
    int id_phase = 0;
    int current_size;

    ContactData phase = pData.contacts_[id_phase];
    // compute some constant matrice for the current phase :
    const Vector3& g = phase.contactPhase_->m_gravity;
    //const Matrix3 gSkew = bezier_com_traj::skew(g);
    phase.contactPhase_->getPolytopeInequalities(Hrow,h);
    H = -Hrow;
    H.rowwise().normalize();
    int dimH = (int)(H.rows());
    mH = phase.contactPhase_->m_mass * H;

    for(int id_step = 0 ; id_step < numStep ; ++id_step , t+=timeStep){
        // add constraints for wp id_step, on current phase :
        // add kinematics constraints :
        // constraint are of the forme A c <= b . But here c(x) = Fx + s so => AFx <= b - As
        current_size = phase.kin_.rows();
        A.block(0,id_rows,current_size,3) = (phase.Kin_ * wps[id_step].first);
        b.segment(id_rows,current_size) = phase.kin_ - (phase.Kin_*wps[id_step].second);
        id_rows += current_size;
        // add stability constraints :
        // compute skew matrix :
        S_hat = skew(g*wps[id_step].first + wps[id_step].second*acc_wps[id_step].first - acc_wps[id_step].second*wps[id_step].first);
        A.block(0,id_rows,dimH,3) = mH.block(0,3,dimH,3)*S_hat + mH.block(0,0,dimH,3)*acc_wps[id_step].first;
        b.segment(id_rows,dimH) = h + mH.block(0,0,dimH,3)*(g - acc_wps[id_step].second) + (acc_wps[id_step].second - g).cross(wps[id_step].second);
        id_rows += dimH ;

        // check if we are going to switch phases :
        for(int i = 0 ; i < stepPerPhases.size() ; ++i){
            if(id_step == stepPerPhases[i]){
                // switch to phase i
                phase = pData.contacts_[id_phase];
                phase.contactPhase_->getPolytopeInequalities(Hrow,h);
                H = -Hrow;
                H.rowwise().normalize();
                dimH = (int)(H.rows());
                mH = phase.contactPhase_->m_mass * H;
                // the current waypoint must have the constraints of both phases. So we add it again :
                // TODO : filter for redunbdant constraints ...
                current_size = phase.kin_.rows();
                A.block(0,id_rows,current_size,3) = (phase.Kin_ * wps[id_step].first);
                b.segment(id_rows,current_size) = phase.kin_ - (phase.Kin_*wps[id_step].second);
                id_rows += current_size;
                // add stability constraints :
                // compute skew matrix :
                S_hat = skew(g*wps[id_step].first + wps[id_step].second*acc_wps[id_step].first - acc_wps[id_step].second*wps[id_step].first);
                A.block(0,id_rows,dimH,3) = mH.block(0,3,dimH,3)*S_hat + mH.block(0,0,dimH,3)*acc_wps[id_step].first;
                b.segment(id_rows,dimH) = h + mH.block(0,0,dimH,3)*(g - acc_wps[id_step].second) + (acc_wps[id_step].second - g).cross(wps[id_step].second);
                id_rows += dimH ;
            }
        }
    }

    return std::make_pair(A,b);
}


std::pair<MatrixX3, VectorX> computeCostFunctionOneStep(const ProblemData& pData){
    Vector3 midPoint = (pData.c0_ + pData.c1_)/2.; // todo : replace it with point found by planning ??
    //cost : min distance between x and midPoint :
    MatrixXX H = Matrix3::Identity();
    VectorX g = -midPoint;
    return std::make_pair(H,g);
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
    assert(pData.contacts_.size() ==2 || pData.contacts_.size() ==3);
    assert(Ts.size() == pData.contacts_.size());
   // bool fail = true;
    ResultDataCOMTraj res;
    std::pair<MatrixX3, VectorX> Ab = computeConstraintsOneStep(pData,Ts,timeStep);
    std::pair<MatrixX3, VectorX> Hg = computeCostFunctionOneStep(pData);
    Vector3 midPoint = (pData.c0_ + pData.c1_)/2.; // todo : replace it with point found by planning ??
    std::cout<<"Init = "<<std::endl<<midPoint.transpose()<<std::endl;

    // rewriting 0.5 || Dx -d ||^2 as x'Hx  + g'x
    ResultData resQp = solve(Ab.first,Ab.second,Hg.first,Hg.second, midPoint);
    if(resQp.success_)
    {
        res.success_ = true;
        res.x = resQp.x;
        computeBezierCurve (pData,Ts,res);
    }
    std::cout<<"Solved, success = "<<res.success_<<" x = "<<res.x.transpose()<<std::endl;
    std::cout<<"Final cost : "<<resQp.cost_<<std::endl;
    return res;
}


} // namespace
