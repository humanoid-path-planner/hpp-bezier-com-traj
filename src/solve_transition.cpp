/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>
#include  <limits>
#include <bezier-com-traj/waypoints/waypoints_definition.hh>

#ifndef QHULL
#define QHULL 1
#endif

#ifndef USE_SLACK
#define USE_SLACK 0
#endif


namespace bezier_com_traj
{
typedef waypoint3_t waypoint_t;
typedef std::pair<double,point3_t> coefs_t;

ResultData solveIntersection(const std::pair<MatrixXX, VectorX>& Ab,const std::pair<MatrixXX, VectorX>& Hg,  const Vector3& init)
{
    return solve(Ab.first,Ab.second,Hg.first,Hg.second, init);
}

void printQHullFile(const std::pair<MatrixXX, VectorX>& Ab,VectorX intPoint,const std::string& fileName,bool clipZ){
     std::ofstream file;
     using std::endl;
     std::string path("/local/fernbac/bench_iros18/constraints_obj/");
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







/**
 * @brief computeDiscretizedTime build an array of discretized points in time, such that there is the same number of point in each phase. Doesn't contain t=0, is of size pointsPerPhase*phaseTimings.size()
 * @param phaseTimings
 * @param pointsPerPhase
 * @return
 */
std::vector<double> computeDiscretizedTime(const VectorX& phaseTimings,const int pointsPerPhase ){
    std::vector<double> timeArray;
    double t = 0;
    double t_total = 0;
    for(size_t i = 0 ; i < phaseTimings.size() ; ++i)
        t_total += phaseTimings[i];

    for(int i = 0 ; i < phaseTimings.size() ; ++i){
        double step = (double) phaseTimings[i] / pointsPerPhase;
        for(int j = 0 ; j < pointsPerPhase ; ++j){
            t += step;
            timeArray.push_back(t);
        }
    }
    timeArray.pop_back();
    timeArray.push_back(t_total); // avoid numerical errors
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
        wps.push_back(evaluateCurveAtTime(pData,pi,t));
    }
    return wps;
}


std::vector<waypoint6_t> computeDiscretizedWwaypoints(const ProblemData& pData,double T,const std::vector<double>& timeArray){
    std::vector<waypoint6_t> wps = computeWwaypoints(pData,T);
    std::vector<waypoint6_t> res;
    std::vector<spline::Bern<double> > berns = ComputeBersteinPolynoms(wps.size()-1);
    double t;
    double b;
    for(int i = 0 ; i < timeArray.size() ; ++i){
        waypoint6_t w = initwp<waypoint6_t>();
        for (int j = 0 ; j < wps.size() ; ++j){
            t = timeArray[i]/T;
            if(t>1.)
                t=1.;
            b = berns[j](t);
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
        wps.push_back(evaluateAccelerationCurveAtTime(pData,pi,T,t));
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
     Normalize(A,b);
     // add 1 for the slack variable :
     A.block(0,3,dimH,1) = VectorX::Ones(dimH);
     return std::make_pair<MatrixXX,VectorX>(A,b);
}

std::pair<MatrixXX,VectorX> dynamicStabilityConstraints(const MatrixXX& mH,const VectorX& h,const Vector3& g,const waypoint6_t& w){
    int dimH = (int)(mH.rows());
    int numCol;
    #if USE_SLACK
    numCol = 4;
    #else
    numCol = 3;
    #endif
    MatrixXX A(dimH,numCol);
    VectorX b(dimH);
    VectorX g_= VectorX::Zero(6);
    g_.head<3>() = g;
    A.block(0,0,dimH,3) = mH*w.first;
    b = h + mH*(g_ - w.second);
    Normalize(A,b);
    #if USE_SLACK
    A.block(0,3,dimH,1) = VectorX::Ones(dimH); // slack variable, after normalization !!
    #endif
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
    //std::cout<<"total time : "<<t_total<<std::endl;
    std::vector<double> timeArray = computeDiscretizedTime(Ts,pointsPerPhase);
    std::vector<coefs_t> wps_c = computeDiscretizedWaypoints(pData,t_total,timeArray);
    std::vector<coefs_t> wps_ddc;
    Vector3 acc_bounds;
    if(pData.constraints_.constraintAcceleration_){
        wps_ddc = computeDiscretizedAccelerationWaypoints(pData,t_total,timeArray);
        acc_bounds = Vector3::Ones()*pData.constraints_.maxAcceleration_;
    }
    std::vector<waypoint6_t> wps_w = computeDiscretizedWwaypoints(pData,t_total,timeArray);
    //std::cout<<" number of discretized waypoints c: "<<wps_c.size()<<std::endl;
    //std::cout<<" number of discretized waypoints w: "<<wps_w.size()<<std::endl;
    assert(/*wps_c.size() == wps_ddc.size() &&*/  wps_w.size() == wps_c.size());
    std::vector<int> stepIdForPhase; // stepIdForPhase[i] is the id of the last step of phase i / first step of phase i+1 (overlap)
    for(int i = 0 ; i < Ts.size() ; ++i)
        stepIdForPhase.push_back(pointsPerPhase*(i+1)-1);

    assert(stepIdForPhase.back() == (wps_c.size()-1)); // -1 because the first one is the index (start at 0) and the second is the size
    // compute the total number of inequalities (to initialise A and b)
    int numCol;
    #if USE_SLACK
    numCol = 4;
    #else
    numCol = 3;
    #endif
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
        //std::cout<<"constraint size : Kin = "<<pData.contacts_[i].kin_.rows()<<" ; stab : "<<Hrow.rows()<<" times "<<numStepForPhase<<" steps"<<std::endl;
        num_stab_ineq += Hrow.rows() * numStepForPhase;
        if(i == Ts.size()-1)
            --numStepForPhase; // we don't consider kinematics constraints for the last point (because it's independant of x)
        num_kin_ineq += pData.contacts_[i].kin_.rows() * numStepForPhase;
    }
    num_ineq = num_stab_ineq + num_kin_ineq;
    if(pData.constraints_.constraintAcceleration_){
        num_ineq += 2*3 *(wps_c.size()) ; // upper and lower bound on acceleration for each discretized waypoint (exept the first one)
    }
    //std::cout<<"total of inequalities : "<<num_ineq<<std::endl;
    // init constraints matrix :
    MatrixXX A = MatrixXX::Zero(num_ineq,numCol); // 3 + 1 :  because of the slack constraints
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
    if(pData.constraints_.reduce_h_ > 0 )
        h -= VectorX::Ones(h.rows())*pData.constraints_.reduce_h_;

    // assign the Stability constraints  for each discretized waypoints :
    // we don't consider the first point, because it's independant of x.
    for(int id_step = 0 ; id_step <  timeArray.size() ; ++id_step ){
        // add stability constraints :

        Ab_stab = dynamicStabilityConstraints(mH,h,g,wps_w[id_step]);
        //compareStabilityMethods(mH,h,g,wps_c[id_step],wps_ddc[id_step],wps_w[id_step]);
        A.block(id_rows,0,dimH,numCol) = Ab_stab.first;
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
                if(pData.constraints_.reduce_h_ > 0 )
                    h -= VectorX::Ones(h.rows())*pData.constraints_.reduce_h_;
                // the current waypoint must have the constraints of both phases. So we add it again :
                // TODO : filter for redunbdant constraints ...
                // add stability constraints :
                Ab_stab = dynamicStabilityConstraints(mH,h,g,wps_w[id_step]);
                //compareStabilityMethods(mH,h,g,wps_c[id_step],wps_ddc[id_step],wps_w[id_step]);
                A.block(id_rows,0,dimH,numCol) = Ab_stab.first;
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

    if(pData.constraints_.constraintAcceleration_) {   // assign the acceleration constraints  for each discretized waypoints :
        for(int id_step = 0 ; id_step <  timeArray.size() ; ++id_step ){
            A.block(id_rows,0,3,3) = Matrix3::Identity() * wps_ddc[id_step].first; // upper
            b.segment(id_rows,3) = acc_bounds - wps_ddc[id_step].second;
            A.block(id_rows+3,0,3,3) = -Matrix3::Identity() * wps_ddc[id_step].first; // lower
            b.segment(id_rows+3,3) = acc_bounds + wps_ddc[id_step].second;
            id_rows += 6;
        }
    }
    //std::cout<<"id rows : "<<id_rows<<" ; total rows : "<<A.rows()<<std::endl;
    assert(id_rows == (A.rows()) && "The constraints matrices were not fully filled.");
    return std::make_pair(A,b);
}

void addSlackInCost( MatrixXX& H, VectorX& g){
    H(3,3) = 1e9;
    g[3] = 0;
}

//cost : min distance between x and midPoint :
void computeCostMidPoint(const ProblemData& pData, MatrixXX& H, VectorX& g){
    // cost : x' H x + 2 x g'
    Vector3 midPoint = (pData.c0_ + pData.c1_)/2.; // todo : replace it with point found by planning ??
    H.block<3,3>(0,0) = Matrix3::Identity();
    g.head<3>() = -midPoint;
}

void computeCostMinAcceleration(const ProblemData& pData,const VectorX& Ts, const int pointsPerPhase, MatrixXX& H, VectorX& g){
    double t_total = 0.;
    for(int i = 0 ; i < Ts.size() ; ++i)
        t_total+=Ts[i];
    std::vector<double> timeArray = computeDiscretizedTime(Ts,pointsPerPhase);
    std::vector<coefs_t> wps_ddc = computeDiscretizedAccelerationWaypoints(pData,t_total,timeArray);
    // cost : x' H x + 2 x g'
    H.block<3,3>(0,0) = Matrix3::Zero();
    g.head<3>() = Vector3::Zero();
    for(size_t i = 0 ; i < wps_ddc.size() ; ++i){
        H.block<3,3>(0,0) += Matrix3::Identity() * wps_ddc[i].first * wps_ddc[i].first;
        g.head<3>() += wps_ddc[i].first*wps_ddc[i].second;
    }
}

//cost : min distance between end velocity and the one computed by planning
void computeCostEndVelocity(const ProblemData& pData,const double T, MatrixXX& H, VectorX& g){
    coefs_t v = computeFinalVelocityPoint(pData,T);
    H.block<3,3>(0,0) = Matrix3::Identity() * v.first * v.first;
    g.head<3> () = v.first*(v.second - pData.dc1_);
}


void computeBezierCurve(const ProblemData& pData, const double T, ResultDataCOMTraj& res)
{
    std::vector<Vector3> wps;

    std::vector<Vector3> pi = computeConstantWaypoints(pData,T);
    size_t i = 0;
    if(pData.constraints_.c0_){
        wps.push_back(pi[i]);
        i++;
        if(pData.constraints_.dc0_){
            wps.push_back(pi[i]);
            i++;
            if(pData.constraints_.ddc0_){
                wps.push_back(pi[i]);
                i++;
            }
        }
    }
    wps.push_back(res.x);
    i++;
    if(pData.constraints_.ddc1_){
        assert(pData.constraints_.dc1_ && "You cannot constraint final acceleration if final velocity is not constrained.");
        wps.push_back(pi[i]);
        i++;
    }
    if(pData.constraints_.dc1_){
        assert(pData.constraints_.c1_ && "You cannot constraint final velocity if final position is not constrained.");
        wps.push_back(pi[i]);
        i++;
    }
    if(pData.constraints_.c1_){
        wps.push_back(pi[i]);
        i++;
    }

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
    int sizeX;
    #if USE_SLACK
    sizeX = 4;
    #else
    sizeX = 3;
    #endif
    MatrixXX H(sizeX,sizeX);
    VectorX g(sizeX);
    ResultDataCOMTraj res;
    VectorX constraint_equivalence;
    std::pair<MatrixXX, VectorX> Ab = computeConstraintsOneStep(pData,Ts,pointsPerPhase,constraint_equivalence);
    //computeCostMidPoint(pData,H,g);
    if(pData.constraints_.dc1_)
        computeCostMinAcceleration(pData,Ts,pointsPerPhase,H,g);
    else
        computeCostEndVelocity(pData,T,H,g);

    #if USE_SLACK
    addSlackInCost(H,g);
    #endif

    //std::cout<<"Init = "<<std::endl<<init_guess.transpose()<<std::endl;

    VectorX x = VectorX::Zero(sizeX); // 3 + slack
    x.head<3>() = init_guess;

    // rewriting 0.5 || Dx -d ||^2 as x'Hx  + g'x
    ResultData resQp = solve(Ab.first,Ab.second,H,g, x);
    bool success;
     #if USE_SLACK
    double feasability = fabs(resQp.x[3]);
    //std::cout<<"feasability : "<<feasability<<"     treshold = "<<feasability_treshold<<std::endl;
    success = feasability<=feasability_treshold;
    #else
    success = resQp.success_;
    #endif

    if(success)
    {
        res.success_ = true;
        res.x = resQp.x.head<3>();
        computeBezierCurve (pData,T,res);
        computeFinalVelocity(pData,T,res);
        computeFinalAcceleration(pData,T,res);
        //std::cout<<"Solved, success "<<" x = ["<<res.x[0]<<","<<res.x[1]<<","<<res.x[2]<<"]"<<std::endl;
        #if QHULL
        printQHullFile(Ab,resQp.x,"bezier_wp.txt");
        #endif
    }else{
        //std::cout<<"Over treshold,  x = ["<<resQp.x[0]<<","<<resQp.x[1]<<","<<resQp.x[2]<<"]"<<std::endl;
    }

    //std::cout<<"Final cost : "<<resQp.cost_<<std::endl;
    return res;
}


} // namespace
