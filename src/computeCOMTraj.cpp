/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <hpp/bezier-com-traj/solve.hh>
#include <hpp/bezier-com-traj/common_solve_methods.hh>
#include <hpp/bezier-com-traj/waypoints/waypoints_definition.hh>
#include <hpp/bezier-com-traj/cost/costfunction_definition.hh>

#include <hpp/bezier-com-traj/solver/solver-abstract.hpp>

#include <hpp/centroidal-dynamics/centroidal_dynamics.hh>

#include  <limits>
#include  <algorithm>

static const int dimVarX = 3;
static const int numRowsForce = 6;


namespace bezier_com_traj
{
typedef std::pair<double,point3_t> coefs_t;

bezier_wp_t::t_point_t computeDiscretizedWwaypoints(const ProblemData& pData,double T, const T_time& timeArray)
{
    bezier_wp_t::t_point_t wps = computeWwaypoints(pData,T);
    bezier_wp_t::t_point_t res;
    const int DIM_VAR = (int)dimVar(pData);
    std::vector<spline::Bern<double> > berns = ComputeBersteinPolynoms((int)wps.size()-1);
    double t, b;
    for (CIT_time cit = timeArray.begin(); cit != timeArray.end(); ++cit)
    {
        waypoint_t w = initwp(6,DIM_VAR);
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
     return std::make_pair(A,b);
}

std::pair<MatrixXX,VectorX> dynamicStabilityConstraints(const MatrixXX& mH,const VectorX& h,const Vector3& g,const waypoint_t& w){
    int dimH = (int)(mH.rows());
    MatrixXX A(dimH,dimVarX);
    VectorX b(dimH);
    VectorX g_= VectorX::Zero(6);
    g_.head<3>() = g;
    A.block(0,0,dimH,3) = mH*w.first;
    b = h + mH*(g_ - w.second);
    Normalize(A,b);
    return std::make_pair(A,b);
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
        numStepsCumulated = phaseSwitch[i]+1;
        if(i > 0 )
            ++numStepForPhase; // because at the switch point between phases we add the constraints of both phases.
        num_stab_ineq += Hrow.rows() * numStepForPhase;
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
    centroidal_dynamics::LP_status status =  phase.contactPhase_->getPolytopeInequalities(Hrow,hrow);
    assert(status == centroidal_dynamics::LP_STATUS_OPTIMAL && "Error in centroidal dynamics lib while computing inequalities");
    mH = -Hrow * phase.contactPhase_->m_mass;
    mH.rowwise().normalize();
    h = hrow;
    dimH = (int)(mH.rows());
    if(pData.constraints_.reduce_h_ > 0 )
        h -= VectorX::Ones(h.rows())*pData.constraints_.reduce_h_;
}

void assignStabilityConstraintsForTimeStep(MatrixXX& mH, VectorX& h, const waypoint_t& wp_w, const int dimH, long int& id_rows,
              MatrixXX& A, VectorX& b, const Vector3& g)
{
    std::pair<MatrixXX,VectorX> Ab_stab = dynamicStabilityConstraints(mH,h,g,wp_w);
    A.block(id_rows,0,dimH,dimVarX) = Ab_stab.first;
    b.segment(id_rows,dimH) = Ab_stab.second;
    id_rows += dimH ;
}

// mG is -G
void assignStabilityConstraintsForTimeStepForce(const waypoint_t& wp_w, const long int rowIdx, long int& id_cols, const centroidal_dynamics::Matrix6X mG,
              MatrixXX& D, VectorX& d, const double mass, const Vector3& g)
{
    D.block(rowIdx,id_cols,numRowsForce,mG.cols()) = mG;
    id_cols += mG.cols() ;
    D.block(rowIdx,0,numRowsForce,dimVarX) = mass * wp_w.first;
    VectorX g_= VectorX::Zero(6);
    g_.head<3>() = g;
    d.segment(rowIdx,numRowsForce) = mass * (g_ - wp_w.second);
}

void switchContactPhase(const ProblemData& pData,
                        MatrixXX& A, VectorX& b,
                        MatrixXX& mH, VectorX& h,
                        const waypoint_t& wp_w, const ContactData& phase, long int& id_rows, int& dimH)
{
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
    bezier_wp_t::t_point_t wps_w = computeDiscretizedWwaypoints(pData,t_total,timeArray);

    std::size_t id_phase = 0;
    const Vector3& g = pData.contacts_[id_phase].contactPhase_->m_gravity;

    //reference to current stability matrix
    MatrixXX mH; VectorX h; int dimH;
    updateH(pData, pData.contacts_[id_phase], mH, h, dimH);

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
                            wps_w[id_step], pData.contacts_[id_phase], id_rows, dimH);
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
    long int current_size;
    id_phase = 0;
    // assign the Kinematics constraints  for each discretized waypoints :
    // we don't consider the first point, because it's independant of x.
    for(std::size_t id_step = 0 ; id_step <  timeArray.size() ; ++id_step )
    {
        // add constraints for wp id_step, on current phase :
        // add kinematics constraints :
        // constraint are of the shape A c <= b . But here c(x) = Fx + s so : AFx <= b - As
        current_size = (int)pData.contacts_[id_phase].kin_.rows();
        if(current_size > 0)
        {
            A.block(id_rows,0,current_size,3) = (pData.contacts_[id_phase].Kin_ * wps_c[id_step].first);
            b.segment(id_rows,current_size) = pData.contacts_[id_phase].kin_ - (pData.contacts_[id_phase].Kin_*wps_c[id_step].second);
            id_rows += current_size;
        }

        // check if we are going to switch phases :
        for(std::size_t i = 0 ; i < (stepIdForPhase.size()-1) ; ++i)
        {
            if(id_step == (std::size_t)stepIdForPhase[i])
            {
                // switch to phase i
                id_phase=i+1;
                // the current waypoint must have the constraints of both phases. So we add it again :
                // TODO : filter for redunbdant constraints ...
                current_size = pData.contacts_[id_phase].kin_.rows();
                if(current_size > 0)
                {
                    A.block(id_rows,0,current_size,3) = (pData.contacts_[id_phase].Kin_ * wps_c[id_step].first);
                    b.segment(id_rows,current_size) = pData.contacts_[id_phase].kin_ - (pData.contacts_[id_phase].Kin_*wps_c[id_step].second);
                    id_rows += current_size;
                }
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
    MatrixXX A = MatrixXX::Zero(num_ineq,dimVarX); VectorX b(num_ineq);

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
                                             const double T, const T_time& timeArray, const long int& dim)
{
    MatrixXX H =  MatrixXX::Identity(dim,dim)*1e-6;
    VectorX g = VectorX::Zero(dim);
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
        res.cost_ = resQp.cost_;
        std::vector<Vector3> pis = computeConstantWaypoints(pData,T);
        res.c_of_t_ = computeBezierCurve<bezier_t, point_t> (pData.constraints_.flag_,T,pis,res.x);
        computeFinalVelocity(res);
        computeFinalAcceleration(res);
        res.dL_of_t_ = bezier_t::zero(T);
    }
    return res;
}


/**
 * @brief computeNumIneqContinuous compute the number of inequalitie required by all the phases
 * @param pData
 * @param Ts
 * @param degree
 * @param w_degree //FIXME : cannot use 2n+3 because for capturability the degree doesn't correspond (cf waypoints_c0_dc0_dc1 )
 * @param useDD whether double description or force formulation is used to compute dynamic constraints)
 * @return
 */
long int computeNumIneqContinuous(const ProblemData& pData, const VectorX& Ts,const int degree,const int w_degree, const bool useDD){
    long int num_ineq = 0;
    centroidal_dynamics::MatrixXX Hrow; VectorX h;
    for(std::vector<ContactData>::const_iterator it = pData.contacts_.begin() ; it != pData.contacts_.end() ; ++it){
        //kinematics :
        num_ineq += it->kin_.rows()*(degree+1);
        //stability :
        if(useDD)
        {
            it->contactPhase_->getPolytopeInequalities(Hrow,h);
            num_ineq += Hrow.rows()*(w_degree+1);
        }
    }
    // acceleration constraints : 6 per points
    num_ineq += (degree-1)*6*Ts.size();
    return num_ineq;
}

/**
 * @brief computeNumEqContinuous compute the number of equalities required by all the phases
 * @param pData
 * @param w_degree //FIXME : cannot use 2n+3 because for capturability the degree doesn't correspond (cf waypoints_c0_dc0_dc1 )
 * @param useForce whether double description or force formulation is used to compute dynamic constraints)
 * @param forceVarDim reference value that stores the total force variable size
 * @return
 */
long int computeNumEqContinuous(const ProblemData& pData, const int w_degree, const bool useForce, long int& forceVarDim){
    long int numEq = 0;
    if(useForce)
    {
        long int cSize = 0;
        int currentDegree; //remove redundant equalities depending on more constraining problem
        std::vector<ContactData>::const_iterator it2 = pData.contacts_.begin(); ++it2;
        for(std::vector<ContactData>::const_iterator it = pData.contacts_.begin() ; it != pData.contacts_.end() ; ++it, ++it2)
        {
            currentDegree = w_degree +1;
            if(cSize != 0 && it->contactPhase_->m_G_centr.cols() > cSize)
                currentDegree -= 1;
            cSize =it->contactPhase_->m_G_centr.cols();
            if(it2 != pData.contacts_.end() && it2->contactPhase_->m_G_centr.cols() <= cSize)
                currentDegree -= 1;
            forceVarDim += cSize * (currentDegree);
            numEq += (numRowsForce ) * (currentDegree);
        }
    }
    return numEq;
}

std::pair<MatrixXX, VectorX> computeConstraintsContinuous(const ProblemData& pData, const VectorX& Ts,
                                                          std::pair<MatrixXX, VectorX>& Dd,
                                                          VectorX& minBounds,
                                                          VectorX& /*maxBounds*/){

    // determine whether to use force or double description
    // based on equilibrium value
    bool useDD = pData.representation_ == DOUBLE_DESCRIPTION;//!(pData.contacts_.begin()->contactPhase_->getAlgorithm() == centroidal_dynamics::EQUILIBRIUM_ALGORITHM_PP);

    double T = Ts.sum();

    // create the curves for c and w with symbolic waypoints (ie. depend on y)
    bezier_wp_t::t_point_t wps_c = computeConstantWaypointsSymbolic(pData,T);
    bezier_wp_t::t_point_t wps_w = computeWwaypoints(pData,T);
    bezier_wp_t c(wps_c.begin(),wps_c.end(),T);
    bezier_wp_t ddc = c.compute_derivate(2);
    bezier_wp_t w(wps_w.begin(),wps_w.end(),T);


    // for each splitted curves : add the constraints for each waypoints
    const long int num_ineq = computeNumIneqContinuous(pData,Ts,(int)c.degree_,(int)w.degree_,  useDD);
    long int forceVarDim = 0;
    const long int num_eq   = computeNumEqContinuous  (pData,(int)w.degree_, !useDD,forceVarDim);
    long int totalVarDim = dimVarX + forceVarDim ;
    long int id_rows = 0;
    long int id_cols = dimVarX; // start after x variable
    //MatrixXX A = MatrixXX::Zero(num_ineq+forceVarDim,totalVarDim);
    MatrixXX A = MatrixXX::Zero(num_ineq,totalVarDim);
    VectorX b = VectorX::Zero(num_ineq);
    MatrixXX& D = Dd.first;
    VectorX& d =  Dd.second;
    D = MatrixXX::Zero(num_eq,totalVarDim);
    d = VectorX::Zero(num_eq);

    double current_t = 0.;
    ContactData phase;
    int size_kin;
    MatrixXX mH; VectorX h; int dimH;
    Vector3 acc_bounds = Vector3::Ones()*pData.constraints_.maxAcceleration_;

    long int cSize = 0;
    long int rowIdx = 0;
    // for each phases, split the curve and compute the waypoints of the splitted curves
    for(size_t id_phase = 0 ; id_phase < (size_t)Ts.size() ; ++id_phase){
        bezier_wp_t cs = c.extract(current_t, Ts[id_phase]+current_t);
        bezier_wp_t ddcs = ddc.extract(current_t, Ts[id_phase]+current_t);
        bezier_wp_t ws = w.extract(current_t , Ts[id_phase]+current_t);
        current_t += Ts[id_phase];

        phase = pData.contacts_[id_phase];

        // add kinematics constraints :
        size_kin = (int)phase.kin_.rows();
        for(bezier_wp_t::cit_point_t wpit = cs.waypoints().begin() ; wpit != cs.waypoints().end() ; ++wpit){
            A.block(id_rows,0,size_kin,3) = (phase.Kin_ * wpit->first);
            b.segment(id_rows,size_kin) = phase.kin_ - (phase.Kin_*wpit->second);
            id_rows += size_kin;
        }
        // add stability constraints
        if(useDD)
        {
            updateH(pData, phase, mH, h, dimH);
            for(bezier_wp_t::cit_point_t wpit = ws.waypoints().begin() ; wpit != ws.waypoints().end() ; ++wpit)
                assignStabilityConstraintsForTimeStep(mH, h, *wpit, dimH, id_rows, A, b, phase.contactPhase_->m_gravity);
        }
        else
        {
            bezier_wp_t::cit_point_t start = ws.waypoints().begin() ;
            bezier_wp_t::cit_point_t stop = ws.waypoints().end() ;
            if (id_phase +1 < (size_t)Ts.size() && pData.contacts_[id_phase+1].contactPhase_->m_G_centr.cols() <= phase.contactPhase_->m_G_centr.cols())
                --stop;
            if (cSize > 0 && cSize < phase.contactPhase_->m_G_centr.cols())
                ++start;
            cSize = phase.contactPhase_->m_G_centr.cols();
            for(bezier_wp_t::cit_point_t wpit = start ; wpit != stop ; ++wpit, rowIdx+=6)
                assignStabilityConstraintsForTimeStepForce(*wpit, rowIdx, id_cols, phase.contactPhase_->m_G_centr, D, d, phase.contactPhase_->m_mass, phase.contactPhase_->m_gravity);
        }
        // add acceleration constraints :
        for(bezier_wp_t::cit_point_t wpit = ddcs.waypoints().begin() ; wpit != ddcs.waypoints().end() ; ++wpit){
            A.block(id_rows,0,3,3) =  wpit->first; // upper
            b.segment(id_rows,3) = acc_bounds - wpit->second;
            A.block(id_rows+3,0,3,3) = - wpit->first; // lower
            b.segment(id_rows+3,3) = acc_bounds + wpit->second;
            id_rows += 6;
        }
    }
    if(!useDD)
    {
        // add positive constraints on forces
        //A.block(id_rows,3,forceVarDim,forceVarDim)=MatrixXX::Identity(forceVarDim,forceVarDim)*(-1);
        //id_rows += forceVarDim;
        minBounds = VectorX::Zero(A.cols()); // * (-std::numeric_limits<double>::infinity());
        minBounds.head<3>() = VectorX::Ones(3) * (solvers::UNBOUNDED_DOWN);
        minBounds.tail(forceVarDim) = VectorX::Zero(forceVarDim) ;
    }
    assert(id_rows == (A.rows()) && "The inequality constraints matrices were not fully filled.");
    assert(id_rows == (b.rows()) && "The inequality constraints matrices were not fully filled.");
    assert(id_cols == (D.cols()) && "The equality constraints matrices were not fully filled.");
    assert(rowIdx  == (d.rows()) && "The equality constraints matrices were not fully filled.");
    //int out = Normalize(D,d);
    return std::make_pair(A,b);
}


ResultDataCOMTraj computeCOMTrajFixedSize(const ProblemData& pData, const VectorX& Ts,
                               const unsigned int pointsPerPhase)
{
    assert (pData.representation_ == DOUBLE_DESCRIPTION);
    if(Ts.size() != (int) pData.contacts_.size())
        throw std::runtime_error("Time phase vector has different size than the number of contact phases");
    double T = Ts.sum();
    T_time timeArray = computeDiscretizedTimeFixed(Ts,pointsPerPhase);
    std::pair<MatrixXX, VectorX> Ab = computeConstraintsOneStep(pData,Ts,T,timeArray);
    std::pair<MatrixXX, VectorX> Hg = genCostFunction(pData,Ts,T,timeArray, dimVarX);
    VectorX x = VectorX::Zero(dimVarX);
    ResultData resQp = solve(Ab,Hg, x);
#if PRINT_QHULL_INEQ
    if (resQp.success_) printQHullFile(Ab,resQp.x, "bezier_wp.txt");
#endif
    return genTraj(resQp, pData, T);
}

ResultDataCOMTraj   computeCOMTraj(const ProblemData& pData, const VectorX& Ts,
                               const double timeStep, const solvers::SolverType solver)
{
    if(Ts.size() != (int) pData.contacts_.size())
        throw std::runtime_error("Time phase vector has different size than the number of contact phases");
    double T = Ts.sum();
    T_time timeArray;
    std::pair<MatrixXX, VectorX> Ab;
    std::pair<MatrixXX, VectorX> Dd;
    VectorX minBounds, maxBounds;
    if(timeStep > 0 ){ // discretized
        assert (pData.representation_ == DOUBLE_DESCRIPTION);
        timeArray = computeDiscretizedTime(Ts,timeStep);
        Ab = computeConstraintsOneStep(pData,Ts,T,timeArray);
        Dd = std::make_pair(MatrixXX::Zero(0,Ab.first.cols()),VectorX::Zero(0));
        minBounds =VectorX::Zero(0);
        maxBounds =VectorX::Zero(0);
    }else{ // continuous
        Ab = computeConstraintsContinuous(pData,Ts, Dd, minBounds, maxBounds);
        timeArray = computeDiscretizedTimeFixed(Ts,7); // FIXME : hardcoded value for discretization for cost function in case of continuous formulation for the constraints
    }
    std::pair<MatrixXX, VectorX> Hg = genCostFunction(pData,Ts,T,timeArray,Ab.first.cols());
    VectorX x = VectorX::Zero(Ab.first.cols());
    ResultData resQp = solve(Ab,Dd,Hg, minBounds, maxBounds, x, solver);
#if PRINT_QHULL_INEQ
    if (resQp.success_) printQHullFile(Ab,resQp.x, "bezier_wp.txt");
#endif
    return genTraj(resQp, pData, T);
}


} // namespace
