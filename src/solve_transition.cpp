/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>
#include <bezier-com-traj/waypoints/waypoints_definition.hh>
#include <bezier-com-traj/cost/costfunction_definition.hh>

#include <centroidal-dynamics-lib/centroidal_dynamics.hh>

#include  <limits>
#include  <algorithm>

#ifndef QHULL
#define QHULL 1
#endif

const int numCol = 3;


namespace bezier_com_traj
{
typedef waypoint3_t waypoint_t;
typedef std::pair<double,point3_t> coefs_t;

std::vector<waypoint6_t> computeDiscretizedWwaypoints(const ProblemData& pData,double T, const T_time& timeArray)
{
    std::vector<waypoint6_t> wps = computeWwaypoints(pData,T);
    std::vector<waypoint6_t> res;
    std::vector<spline::Bern<double> > berns = ComputeBersteinPolynoms((int)wps.size()-1);
    double t, b;
    for (CIT_time cit = timeArray.begin(); cit != timeArray.end(); ++cit)
    {
        waypoint6_t w = initwp<waypoint6_t>();
        t = std::min(cit->first / T, 1.);
        for (std::size_t j = 0 ; j < wps.size() ; ++j)
        {
            b = berns[j](t);
            w.first +=b*(wps[j].first );
            w.second+=b*(wps[j].second);
        }
        res.push_back(w);
    }
    return res;
}

 std::pair<MatrixXX,VectorX> dynamicStabilityConstraints_cross(const MatrixXX& mH,const VectorX& h,const Vector3& g,
                                                               const coefs_t& c,const coefs_t& ddc)
{
     Matrix3 S_hat;
     int dimH = (int)(mH.rows());
     MatrixXX A(dimH,4);
     VectorX b(dimH);
     S_hat = skew(c.second*ddc.first - ddc.second*c.first + g*c.first);
     A.block(0,0,dimH,3) = mH.block(0,3,dimH,3) * S_hat + mH.block(0,0,dimH,3) * ddc.first;
     b = h + mH.block(0,0,dimH,3)*(g - ddc.second) + mH.block(0,3,dimH,3)*(c.second.cross(g) - c.second.cross(ddc.second));
     Normalize(A,b);
     // add 1 for the slack variable :
     A.block(0,3,dimH,1) = VectorX::Ones(dimH);
     return std::make_pair<MatrixXX,VectorX>(A,b);
}

std::pair<MatrixXX,VectorX> dynamicStabilityConstraints(const MatrixXX& mH,const VectorX& h,const Vector3& g,const waypoint6_t& w){
    int dimH = (int)(mH.rows());
    MatrixXX A(dimH,numCol);
    VectorX b(dimH);
    VectorX g_= VectorX::Zero(6);
    g_.head<3>() = g;
    A.block(0,0,dimH,3) = mH*w.first;
    b = h + mH*(g_ - w.second);
    Normalize(A,b);
    return std::make_pair<MatrixXX,VectorX>(A,b);
}

std::vector<int> stepIdPerPhase(const T_time& timeArray) // const int pointsPerPhase)
{
    std::vector<int> res;
    int i = 0;
    CIT_time cit = timeArray.begin();
    for (CIT_time cit2 = timeArray.begin()+1; cit2 != timeArray.end(); ++cit, ++cit2, ++i)
    {
        if (cit2->second != cit->second)
        {
            res.push_back(i);
        }
    }
    res.push_back(i);
    return res;
}

long int computeNumIneq(const ProblemData& pData, const VectorX& Ts, const std::vector<int>& phaseSwitch)
{
    const size_t numPoints = phaseSwitch.back() +1;
    long int num_stab_ineq = 0;
    long int num_kin_ineq = 0;
    int numStepForPhase;
    int numStepsCumulated= 0;
    centroidal_dynamics::MatrixXX Hrow; VectorX h;
    for(int i = 0 ; i < Ts.size() ; ++i)
    {
        pData.contacts_[i].contactPhase_->getPolytopeInequalities(Hrow,h);
        numStepForPhase = phaseSwitch[i]+1 - numStepsCumulated; // pointsPerPhase;
        numStepsCumulated= phaseSwitch[i]+1;
        if(i > 0 )
            ++numStepForPhase; // because at the switch point between phases we add the constraints of both phases.
        num_stab_ineq += Hrow.rows() * numStepForPhase;
        if(i == Ts.size()-1)
            --numStepForPhase; // we don't consider kinematics constraints for the last point (because it's independant of x)
        num_kin_ineq += pData.contacts_[i].kin_.rows() * numStepForPhase;
    }
    long int res = num_stab_ineq + num_kin_ineq;
    if(pData.constraints_.constraintAcceleration_)
        res += 2*3 *(numPoints) ; // upper and lower bound on acceleration for each discretized waypoint (exept the first one)
    return res;
}

void updateH(const ProblemData& pData, const ContactData& phase, MatrixXX& mH, VectorX& h, int& dimH)
{
    VectorX hrow;
    centroidal_dynamics::MatrixXX Hrow;
    phase.contactPhase_->getPolytopeInequalities(Hrow,hrow);
    mH = -Hrow * phase.contactPhase_->m_mass;
    mH.rowwise().normalize();
    h = hrow;
    dimH = (int)(mH.rows());
    if(pData.constraints_.reduce_h_ > 0 )
        h -= VectorX::Ones(h.rows())*pData.constraints_.reduce_h_;
}

void assignStabilityConstraintsForTimeStep(MatrixXX& mH, VectorX& h, const waypoint6_t& wp_w, const int dimH, long int& id_rows,
              MatrixXX& A, VectorX& b, const Vector3& g)
{
    std::pair<MatrixXX,VectorX> Ab_stab = dynamicStabilityConstraints(mH,h,g,wp_w);
    A.block(id_rows,0,dimH,numCol) = Ab_stab.first;
    b.segment(id_rows,dimH) = Ab_stab.second;    
    id_rows += dimH ;
}

void switchContactPhase(const ProblemData& pData,
                        MatrixXX& A, VectorX& b,
                        MatrixXX& mH, VectorX& h,
                        const waypoint6_t& wp_w,  ContactData& phase,
                        const long int id_phase, long int& id_rows, int& dimH)
{
    phase = pData.contacts_[id_phase];
    updateH(pData, phase, mH, h, dimH);
    // the current waypoint must have the constraints of both phases. So we add it again :
    // TODO : filter for redunbdant constraints ...
    // add stability constraints :
    assignStabilityConstraintsForTimeStep(mH, h, wp_w, dimH, id_rows, A, b, phase.contactPhase_->m_gravity);
}

long int assignStabilityConstraints(const ProblemData& pData, MatrixXX& A, VectorX& b, const T_time& timeArray,
                                    const double t_total, const std::vector<int>& stepIdForPhase)
{
    long int id_rows = 0;
    std::vector<waypoint6_t> wps_w = computeDiscretizedWwaypoints(pData,t_total,timeArray);

    std::size_t id_phase = 0;
    ContactData phase = pData.contacts_[id_phase];
    const Vector3& g = phase.contactPhase_->m_gravity;

    //reference to current stability matrix
    MatrixXX mH; VectorX h; int dimH;
    updateH(pData, phase, mH, h, dimH);

    // assign the Stability constraints  for each discretized waypoints :
    // we don't consider the first point, because it's independant of x.
    for(std::size_t id_step = 0 ; id_step <  timeArray.size() ; ++id_step)
    {
        // add stability constraints :
        assignStabilityConstraintsForTimeStep(mH, h, wps_w[id_step], dimH, id_rows, A, b, g);
        // check if we are going to switch phases :
        for(std::vector<int>::const_iterator it_switch = stepIdForPhase.begin() ; it_switch != (stepIdForPhase.end()-1) ; ++it_switch)
        {
            if((int)id_step == (*it_switch))
            {
                id_phase++;
                switchContactPhase(pData, A,b, mH, h,
                            wps_w[id_step], phase, id_phase, id_rows, dimH);
            }
        }
    }
    return id_rows;
}

void assignKinematicConstraints(const ProblemData& pData, MatrixXX& A, VectorX& b, const T_time& timeArray,
                                const double t_total, const std::vector<int>& stepIdForPhase, long int& id_rows)
{
    std::size_t id_phase = 0;
    std::vector<coefs_t> wps_c = computeDiscretizedWaypoints<point_t>(pData,t_total,timeArray);
    ContactData phase = pData.contacts_[id_phase];
    long int current_size;
    id_phase = 0;
    phase = pData.contacts_[id_phase];
    // assign the Kinematics constraints  for each discretized waypoints :
    // we don't consider the first point, because it's independant of x.
    for(std::size_t id_step = 0 ; id_step <  timeArray.size() ; ++id_step )
    {
        // add constraints for wp id_step, on current phase :
        // add kinematics constraints :
        // constraint are of the shape A c <= b . But here c(x) = Fx + s so : AFx <= b - As
        if(id_step != timeArray.size()-1)
        { // we don't consider kinematics constraints for the last point (because it's independant of x)
            current_size = (int)phase.kin_.rows();
            A.block(id_rows,0,current_size,3) = (phase.Kin_ * wps_c[id_step].first);
            b.segment(id_rows,current_size) = phase.kin_ - (phase.Kin_*wps_c[id_step].second);
            id_rows += current_size;
        }

        // check if we are going to switch phases :
        for(std::size_t i = 0 ; i < (stepIdForPhase.size()-1) ; ++i)
        {
            if(id_step == (std::size_t)stepIdForPhase[i]){
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
}


void assignAccelerationConstraints(const ProblemData& pData, MatrixXX& A, VectorX& b, const T_time& timeArray,
                                const double t_total, long int& id_rows)
{
    if(pData.constraints_.constraintAcceleration_)
    {   // assign the acceleration constraints  for each discretized waypoints :
        std::vector<coefs_t> wps_ddc = computeDiscretizedAccelerationWaypoints<point_t>(pData,t_total,timeArray);
        Vector3 acc_bounds = Vector3::Ones()*pData.constraints_.maxAcceleration_;
        for(std::size_t id_step = 0 ; id_step <  timeArray.size() ; ++id_step )
        {
            A.block(id_rows,0,3,3) = Matrix3::Identity() * wps_ddc[id_step].first; // upper
            b.segment(id_rows,3) = acc_bounds - wps_ddc[id_step].second;
            A.block(id_rows+3,0,3,3) = -Matrix3::Identity() * wps_ddc[id_step].first; // lower
            b.segment(id_rows+3,3) = acc_bounds + wps_ddc[id_step].second;
            id_rows += 6;
        }
    }
}

std::pair<MatrixXX, VectorX> computeConstraintsOneStep(const ProblemData& pData, const VectorX& Ts,
                                                       const double t_total, const T_time& timeArray)
{
    // Compute all the discretized wayPoint
    std::vector<int> stepIdForPhase = stepIdPerPhase(timeArray);
    // stepIdForPhase[i] is the id of the last step of phase i / first step of phase i+1 (overlap)

    // init constraints matrix :
    long int num_ineq = computeNumIneq(pData, Ts, stepIdForPhase);
    MatrixXX A = MatrixXX::Zero(num_ineq,numCol); VectorX b(num_ineq);

    // assign dynamic and kinematic constraints per timestep, including additional acceleration
    // constraints.
    long int id_rows = assignStabilityConstraints(pData, A, b, timeArray, t_total, stepIdForPhase);
    assignKinematicConstraints(pData, A, b, timeArray, t_total, stepIdForPhase, id_rows);
    assignAccelerationConstraints(pData, A, b, timeArray, t_total, id_rows);

    assert(id_rows == (A.rows()) && "The constraints matrices were not fully filled.");
    return std::make_pair(A,b);
}

void computeFinalAcceleration(ResultDataCOMTraj& res){
    res.ddc1_ = res.c_of_t_.derivate(res.c_of_t_.max(), 2);
}

void computeFinalVelocity(ResultDataCOMTraj& res){
    res.dc1_ = res.c_of_t_.derivate(res.c_of_t_.max(), 1);
}

std::pair<MatrixXX, VectorX> genCostFunction(const ProblemData& pData,const VectorX& Ts,
                                             const double T, const T_time& timeArray)
{
    MatrixXX H(numCol,numCol);
    VectorX g(numCol);
    cost::genCostFunction(pData,Ts,T,timeArray,H,g);
    return std::make_pair(H,g);
}

ResultDataCOMTraj genTraj(ResultData resQp, const ProblemData& pData, const double T )
{
    ResultDataCOMTraj res;
    if(resQp.success_)
    {
        res.success_ = true;
        res.x = resQp.x.head<3>();
        std::vector<Vector3> pis = computeConstantWaypoints(pData,T);
        res.c_of_t_ = computeBezierCurve<bezier_t, point_t> (pData.constraints_.flag_,T,pis,res.x);
        computeFinalVelocity(res);
        computeFinalAcceleration(res);
    }
    return res;
}

double computeTotalTime(const VectorX& Ts)
{
    double T = 0;
    for(int i = 0 ; i < Ts.size() ; ++i)
        T+=Ts[i];
    return T;
}

ResultDataCOMTraj solveOnestep(const ProblemData& pData, const VectorX& Ts,const Vector3& init_guess,
                               const int pointsPerPhase, const double /*feasability_treshold*/){
    assert(Ts.size() == pData.contacts_.size());
    double T = Ts.sum();
    T_time timeArray = computeDiscretizedTime(Ts,pointsPerPhase);
    std::pair<MatrixXX, VectorX> Ab = computeConstraintsOneStep(pData,Ts,T,timeArray);
    std::pair<MatrixXX, VectorX> Hg = genCostFunction(pData,Ts,T,timeArray);
    VectorX x = VectorX::Zero(numCol); x.head<3>() = init_guess;
    ResultData resQp = solve(Ab,Hg, x);
#if QHULL
    if (resQp.success_) printQHullFile(Ab,resQp.x,"bezier_wp.txt");
#endif
    return genTraj(resQp, pData, T);
}


} // namespace
