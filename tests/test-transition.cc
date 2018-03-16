//
// Copyright (c) 2018 CNRS
// Authors: Pierre Fernbach
//
// This file is part of bezier_COM_traj
// hpp-core is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-core is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-core  If not, see
// <http://www.gnu.org/licenses/>.


#define BOOST_TEST_MODULE transition
#include <boost/test/included/unit_test.hpp>
#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>
#include <bezier-com-traj/data.hh>
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>
#include <test_helper.hh>

#define NOMINAL_COM_HEIGHT 0.795


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

bool check_constraints(const bezier_com_traj::ContactData& contactPhase, Vector3 c,Vector3 dc, Vector3 ddc){

    BOOST_CHECK(verifyKinematicConstraints(std::make_pair(contactPhase.Kin_,contactPhase.kin_),c));
    BOOST_CHECK(verifyStabilityConstraintsDLP(*contactPhase.contactPhase_,c,dc,ddc));
    BOOST_CHECK(verifyStabilityConstraintsPP(*contactPhase.contactPhase_,c,dc,ddc));

}


void check_transition(bezier_com_traj::ProblemData& pData, VectorX Ts,bool shouldFail=false){
    BOOST_CHECK_EQUAL(pData.contacts_.size(),Ts.size());

    double t_total = 0;
    for(size_t i = 0 ; i < Ts.size() ; ++i)
        t_total += Ts[i];

    Vector3 init = (pData.c1_ - pData.c0_)/2.;
    int pointsPerPhase = 5;

    // check if transition is feasible (should be)
    bezier_com_traj::ResultDataCOMTraj res = bezier_com_traj::solveOnestep(pData,Ts,init,pointsPerPhase);
    if(shouldFail){
        BOOST_CHECK(!res.success_);
        return;
    }

    BOOST_CHECK(res.success_);

    if(res.success_){
         // check if timing is respected
        std::vector<double> timings = computeDiscretizedTime(Ts,pointsPerPhase);
        BOOST_CHECK_EQUAL(timings.back(),t_total);

        Vector3 c,dc,ddc;
        bezier_com_traj::bezier_t ct = res.c_of_t_;
        bezier_com_traj::bezier_t dct = ct.compute_derivate(1);
        bezier_com_traj::bezier_t ddct = dct.compute_derivate(1);

        BOOST_CHECK_EQUAL(ct.min(),0);
        BOOST_CHECK_EQUAL(dct.min(),0);
        BOOST_CHECK_EQUAL(ddct.min(),0);
        BOOST_CHECK_EQUAL(ct.max(),t_total);
        BOOST_CHECK_EQUAL(dct.max(),t_total);
        BOOST_CHECK_EQUAL(ddct.max(),t_total);

        // check if each discretized point is feasible :

        std::vector<int> stepIdForPhase; // stepIdForPhase[i] is the id of the last step of phase i / first step of phase i+1 (overlap)
        for(int i = 0 ; i < Ts.size() ; ++i)
            stepIdForPhase.push_back(pointsPerPhase*(i+1)-1);
        int id_phase = 0;
        bezier_com_traj::ContactData phase = pData.contacts_[id_phase];

        for(size_t id_step = 0 ; id_step < timings.size() ; ++id_step){
            c = ct(timings[id_step]);
            dc = dct(timings[id_step]);
            ddc = ddct(timings[id_step]);

            // check i (c,dc,ddc) verify the constraints of current phase
            check_constraints(phase,c,dc,ddc);

            // check if we switch phases
            for(int i = 0 ; i < (stepIdForPhase.size()-1) ; ++i){
                if(id_step == stepIdForPhase[i]){
                    id_phase=i+1;
                    phase = pData.contacts_[id_phase];
                    check_constraints(phase,c,dc,ddc);
                }
            }
        }
    }
    for(size_t i = 0 ; i < pData.contacts_.size() ; ++i){
        delete pData.contacts_[i].contactPhase_;
    }
}



BOOST_AUTO_TEST_SUITE( flat_ground )

// one step (left foot) : init pos of end effector : (0,0.1)(0,-0.1)
//                                         end pos : (0.3,0.1)(0,-0.1)
// three phases : both feet in contact, only right feet, both feet

bezier_com_traj::ContactData phase0_flat(){
    bezier_com_traj::ContactData cData;
    MatrixX3 normals(2,3),positions(2,3);
    normals.block<1,3>(0,0)=Vector3(0,0,1);
    positions.block<1,3>(0,0)=Vector3(0,0.1,0);
    normals.block<1,3>(1,0)=Vector3(0,0,1);
    positions.block<1,3>(1,0)=Vector3(0,-0.1,0);
    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    cData.contactPhase_  = new centroidal_dynamics::Equilibrium(ComputeContactCone(contacts.first,contacts.second));

    ConstraintsPair kin1 = generateKinematicsConstraints(Matrix3::Identity(),Vector3(0,-0.1,0));
    ConstraintsPair kin = stackConstraints(kin1,generateKinematicsConstraints(Matrix3::Identity(),Vector3(0,0.1,0)));
    cData.Kin_ = kin.first;
    cData.kin_ = kin.second;

    return cData;
}

bezier_com_traj::ContactData phase1_flat(){
    bezier_com_traj::ContactData cData;
    MatrixX3 normals(1,3),positions(1,3);
    normals.block<1,3>(0,0)=Vector3(0,0,1);
    positions.block<1,3>(0,0)=Vector3(0,-0.1,0);
    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    cData.contactPhase_  = new centroidal_dynamics::Equilibrium(ComputeContactCone(contacts.first,contacts.second));

    ConstraintsPair kin = generateKinematicsConstraints(Matrix3::Identity(),Vector3(0,-0.1,0));
    cData.Kin_ = kin.first;
    cData.kin_ = kin.second;

    return cData;
}

bezier_com_traj::ContactData phase2_flat(){
    bezier_com_traj::ContactData cData;
    MatrixX3 normals(2,3),positions(2,3);
    normals.block<1,3>(0,0)=Vector3(0,0,1);
    positions.block<1,3>(0,0)=Vector3(0.3,0.1,0);
    normals.block<1,3>(1,0)=Vector3(0,0,1);
    positions.block<1,3>(1,0)=Vector3(0,-0.1,0);
    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    cData.contactPhase_  = new centroidal_dynamics::Equilibrium(ComputeContactCone(contacts.first,contacts.second));

    ConstraintsPair kin1 = generateKinematicsConstraints(Matrix3::Identity(),Vector3(0,-0.1,0));
    ConstraintsPair kin = stackConstraints(kin1,generateKinematicsConstraints(Matrix3::Identity(),Vector3(0.3,0.1,0)));
    cData.Kin_ = kin.first;
    cData.kin_ = kin.second;

    return cData;
}


bezier_com_traj::ProblemData gen_problem_data_flat(){
    bezier_com_traj::ProblemData pData;
    pData.c0_ = Vector3(0,0,NOMINAL_COM_HEIGHT);
    pData.c1_ = Vector3(0.15,0,NOMINAL_COM_HEIGHT);
    pData.dc0_ = Vector3::Zero();
    pData.dc1_ = Vector3::Zero();
    pData.ddc0_ = Vector3::Zero();
    pData.ddc1_ = Vector3::Zero();

    pData.contacts_.push_back(phase0_flat());
    pData.contacts_.push_back(phase1_flat());
    pData.contacts_.push_back(phase2_flat());

    return pData;
}

BOOST_AUTO_TEST_CASE(quasi_static){

// compute kinematic constraints for the right foot :
    ConstraintsPair kin1 = generateKinematicsConstraints(Matrix3::Identity(),Vector3(0,-0.1,0));
    ConstraintsPair kin0 = stackConstraints(kin1,generateKinematicsConstraints(Matrix3::Identity(),Vector3(0,0.1,0)));
    ConstraintsPair kin2 = stackConstraints(kin1,generateKinematicsConstraints(Matrix3::Identity(),Vector3(0.3,0.1,0)));

    bezier_com_traj::ContactData cDataMid = phase1_flat();
    ConstraintsPair stab = generateStabilityConstraints(*cDataMid.contactPhase_);

    std::pair<Matrix3,Vector3> Hg = computeCost();
    Vector3 init = Vector3::Zero();

    ConstraintsPair Ab_first = stackConstraints(kin0,stab);
    bezier_com_traj::ResultData res_first = bezier_com_traj::solveIntersection(Ab_first,Hg,init);
    BOOST_CHECK(res_first.success_);

    ConstraintsPair Ab_second = stackConstraints(kin2,stab);
    bezier_com_traj::ResultData res_second = bezier_com_traj::solveIntersection(Ab_second,Hg,init);
    BOOST_CHECK(res_second.success_);
    delete cDataMid.contactPhase_;
}


BOOST_AUTO_TEST_CASE(transition){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}


BOOST_AUTO_TEST_CASE(transition_noDc1){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ ^= bezier_com_traj::END_VEL;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}

BOOST_AUTO_TEST_CASE(transition_ddc0){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ |= bezier_com_traj::INIT_ACC;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}

BOOST_AUTO_TEST_CASE(transition_ddc0_ddc1){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ |= bezier_com_traj::INIT_ACC;
    pData.constraints_.flag_ |= bezier_com_traj::END_ACC;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}

BOOST_AUTO_TEST_CASE(transition_noAcc){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.constraintAcceleration_ = false;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}

BOOST_AUTO_TEST_CASE(transition_noDc1_noAcc){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ ^= bezier_com_traj::END_VEL;
    pData.constraints_.constraintAcceleration_ = false;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}
BOOST_AUTO_TEST_CASE(transition_ddc0_noAcc){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ |= bezier_com_traj::INIT_ACC;
    pData.constraints_.constraintAcceleration_ = false;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}

BOOST_AUTO_TEST_CASE(transition_ddc0_ddc1_noAcc){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ |= bezier_com_traj::INIT_ACC | bezier_com_traj::END_ACC ;
    pData.constraints_.constraintAcceleration_ = false;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}


BOOST_AUTO_TEST_CASE(transition_Acc1){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.maxAcceleration_=1.;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}

BOOST_AUTO_TEST_CASE(transition_noDc1_Acc1){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ ^= bezier_com_traj::END_VEL;
    pData.constraints_.maxAcceleration_=1.;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}
BOOST_AUTO_TEST_CASE(transition_ddc0_Acc2){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ |= bezier_com_traj::INIT_ACC;
    pData.constraints_.maxAcceleration_=2.;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}

BOOST_AUTO_TEST_CASE(transition_ddc0_ddc1_Acc2){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ |= bezier_com_traj::INIT_ACC | bezier_com_traj::END_ACC ;
    pData.constraints_.maxAcceleration_=2.;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}

BOOST_AUTO_TEST_CASE(transition_ddc0_ddc1_Acc05){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ |= bezier_com_traj::INIT_ACC | bezier_com_traj::END_ACC ;
    pData.constraints_.maxAcceleration_=0.5;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}

BOOST_AUTO_TEST_CASE(transition_Acc05){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.maxAcceleration_=0.5;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts);
}

// constraints that should fails :

BOOST_AUTO_TEST_CASE(transition_Acc02){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.maxAcceleration_=0.2;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts,true);
}

BOOST_AUTO_TEST_CASE(transition_noDc1_Acc05){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ ^= bezier_com_traj::END_VEL;
    pData.constraints_.maxAcceleration_=0.5;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts,true);
}


BOOST_AUTO_TEST_CASE(transition_ddc0_Acc1){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ |= bezier_com_traj::INIT_ACC;
    pData.constraints_.maxAcceleration_=1.;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts,true);
}

BOOST_AUTO_TEST_CASE(transition_ddc0_ddc1_Acc02){
    bezier_com_traj::ProblemData pData = gen_problem_data_flat();
    pData.constraints_.flag_ |= bezier_com_traj::INIT_ACC | bezier_com_traj::END_ACC ;
    pData.constraints_.maxAcceleration_=0.2;
    VectorX Ts(3);
    Ts<<0.6,0.6,0.6;
    check_transition(pData,Ts,true);
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE( platform )

// platform : first step :
// (0.35,0.1,0) ; (0.35,-0.1,0) -> (0.775, 0.23, -0.02);(0.35,-0.1,0)  (normal : 0.0, -0.423, 0.906)

// second step :
// (0.775, 0.23, -0.02);(0.35,-0.1,0) -> (0.775, 0.23, -0.02);(1.15,-0.1,0)
// unfeasible in quasi-static

//third step :
//(0.775, 0.23, -0.02);(1.15,-0.1,0) -> (1.15,0.1,0);(1.15,-0.1,0)

bezier_com_traj::ContactData phase0_platform(){
    bezier_com_traj::ContactData cData;
    MatrixX3 normals(2,3),positions(2,3);
    normals.block<1,3>(0,0)=Vector3(0,0,1);
    positions.block<1,3>(0,0)=Vector3(0.35,0.1,0);
    normals.block<1,3>(1,0)=Vector3(0,0,1);
    positions.block<1,3>(1,0)=Vector3(0.35,-0.1,0);
    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    cData.contactPhase_  = new centroidal_dynamics::Equilibrium(ComputeContactCone(contacts.first,contacts.second));

    ConstraintsPair kin1 = generateKinematicsConstraints(Matrix3::Identity(),Vector3(0.35,-0.1,0));
    ConstraintsPair kin = stackConstraints(kin1,generateKinematicsConstraints(Matrix3::Identity(),Vector3(0.35,0.1,0)));
    cData.Kin_ = kin.first;
    cData.kin_ = kin.second;

    return cData;
}

bezier_com_traj::ContactData phase1_platform(){
    bezier_com_traj::ContactData cData;
    MatrixX3 normals(1,3),positions(1,3);
    normals.block<1,3>(0,0)=Vector3(0,0,1);
    positions.block<1,3>(0,0)=Vector3(0.35,-0.1,0);
    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    cData.contactPhase_  = new centroidal_dynamics::Equilibrium(ComputeContactCone(contacts.first,contacts.second));

    ConstraintsPair kin = generateKinematicsConstraints(Matrix3::Identity(),Vector3(0.35,-0.1,0));
    cData.Kin_ = kin.first;
    cData.kin_ = kin.second;

    return cData;
}

bezier_com_traj::ContactData phase2_platform(){
    bezier_com_traj::ContactData cData;
    MatrixX3 normals(2,3),positions(2,3);
    normals.block<1,3>(0,0)=Vector3(0.0, -0.423, 0.906).normalized();
    positions.block<1,3>(0,0)=Vector3(0.775, 0.23, -0.02);
    normals.block<1,3>(1,0)=Vector3(0,0,1);
    positions.block<1,3>(1,0)=Vector3(0.35,-0.1,0);
    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    cData.contactPhase_  = new centroidal_dynamics::Equilibrium(ComputeContactCone(contacts.first,contacts.second));

    ConstraintsPair kin1 = generateKinematicsConstraints(Matrix3::Identity(),Vector3(0.35,-0.1,0));
    Eigen::Quaterniond quat = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(),Eigen::Vector3d(0.0, -0.423, 0.906));
    Matrix3 rot = quat.normalized().toRotationMatrix();
    ConstraintsPair kin = stackConstraints(kin1,generateKinematicsConstraints(rot,Vector3(0.775, 0.23, -0.02)));
    cData.Kin_ = kin.first;
    cData.kin_ = kin.second;

    return cData;
}

bezier_com_traj::ContactData phase3_platform(){
    bezier_com_traj::ContactData cData;
    MatrixX3 normals(1,3),positions(1,3);
    normals.block<1,3>(0,0)=Vector3(0.0, -0.423, 0.906).normalized();
    positions.block<1,3>(0,0)=Vector3(0.775, 0.23, -0.02);

    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    cData.contactPhase_  = new centroidal_dynamics::Equilibrium(ComputeContactCone(contacts.first,contacts.second));

    Eigen::Quaterniond quat = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(),Eigen::Vector3d(0.0, -0.423, 0.906));
    Matrix3 rot = quat.normalized().toRotationMatrix();
    ConstraintsPair kin = generateKinematicsConstraints(rot,Vector3(0.775, 0.23, -0.02));
    cData.Kin_ = kin.first;
    cData.kin_ = kin.second;

    return cData;
}

bezier_com_traj::ContactData phase4_platform(){
    bezier_com_traj::ContactData cData;
    MatrixX3 normals(2,3),positions(2,3);
    normals.block<1,3>(0,0)=Vector3(0.0, -0.423, 0.906).normalized();
    positions.block<1,3>(0,0)=Vector3(0.775, 0.23, -0.02);
    normals.block<1,3>(1,0)=Vector3(0,0,1);
    positions.block<1,3>(1,0)=Vector3(1.15,-0.1,0);
    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    cData.contactPhase_  = new centroidal_dynamics::Equilibrium(ComputeContactCone(contacts.first,contacts.second));

    ConstraintsPair kin1 = generateKinematicsConstraints(Matrix3::Identity(),Vector3(1.15,-0.1,0));
    Eigen::Quaterniond quat = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(),Eigen::Vector3d(0.0, -0.423, 0.906));
    Matrix3 rot = quat.normalized().toRotationMatrix();
    ConstraintsPair kin = stackConstraints(kin1,generateKinematicsConstraints(rot,Vector3(0.775, 0.23, -0.02)));
    cData.Kin_ = kin.first;
    cData.kin_ = kin.second;

    return cData;
}

bezier_com_traj::ContactData phase5_platform(){
    bezier_com_traj::ContactData cData;
    MatrixX3 normals(1,3),positions(1,3);
    normals.block<1,3>(0,0)=Vector3(0,0,1);
    positions.block<1,3>(0,0)=Vector3(1.15,-0.1,0);

    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    cData.contactPhase_  = new centroidal_dynamics::Equilibrium(ComputeContactCone(contacts.first,contacts.second));

    ConstraintsPair kin = generateKinematicsConstraints(Matrix3::Identity(),Vector3(1.15,-0.1,0));
    cData.Kin_ = kin.first;
    cData.kin_ = kin.second;

    return cData;
}

bezier_com_traj::ContactData phase6_platform(){
    bezier_com_traj::ContactData cData;
    MatrixX3 normals(2,3),positions(2,3);
    normals.block<1,3>(0,0)=Vector3(0,0,1);
    positions.block<1,3>(0,0)=Vector3(1.15,-0.1,0);
    normals.block<1,3>(1,0)=Vector3(0,0,1);
    positions.block<1,3>(1,0)=Vector3(1.15,0.1,0);

    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    cData.contactPhase_  = new centroidal_dynamics::Equilibrium(ComputeContactCone(contacts.first,contacts.second));

    ConstraintsPair kin1 = generateKinematicsConstraints(Matrix3::Identity(),Vector3(1.15,-0.1,0));
    ConstraintsPair kin = stackConstraints(kin1,generateKinematicsConstraints(Matrix3::Identity(),Vector3(1.15,0.1,0)));
    cData.Kin_ = kin.first;
    cData.kin_ = kin.second;

    return cData;
}



bezier_com_traj::ProblemData gen_problem_data_Step0(){
    bezier_com_traj::ProblemData pData;
    pData.c0_ = Vector3(0.35,0,NOMINAL_COM_HEIGHT);
    pData.c1_ = Vector3(0.56,0.04,NOMINAL_COM_HEIGHT);
    pData.dc0_ = Vector3::Zero();
    pData.dc1_ = Vector3::Zero();
    pData.ddc0_ = Vector3::Zero();
    pData.ddc1_ = Vector3::Zero();

    pData.contacts_.push_back(phase0_platform());
    pData.contacts_.push_back(phase1_platform());
    pData.contacts_.push_back(phase2_platform());

    return pData;
}


bezier_com_traj::ProblemData gen_problem_data_Step1(){
    bezier_com_traj::ProblemData pData;
    pData.c0_ = Vector3(0.56,0.04,0.765);
    pData.c1_ = Vector3(0.98,0.02,0.77);
    pData.dc0_ = Vector3::Zero();
    pData.dc1_ = Vector3::Zero();
    pData.ddc0_ = Vector3::Zero();
    pData.ddc1_ = Vector3::Zero();

    pData.contacts_.push_back(phase2_platform());
    pData.contacts_.push_back(phase3_platform());
    pData.contacts_.push_back(phase4_platform());

    return pData;
}


bezier_com_traj::ProblemData gen_problem_data_Step2(){
    bezier_com_traj::ProblemData pData;
    pData.c0_ = Vector3(0.98,0.04,NOMINAL_COM_HEIGHT);
    pData.c1_ = Vector3(1.15,0,NOMINAL_COM_HEIGHT);
    pData.dc0_ = Vector3::Zero();
    pData.dc1_ = Vector3::Zero();
    pData.ddc0_ = Vector3::Zero();
    pData.ddc1_ = Vector3::Zero();

    pData.contacts_.push_back(phase4_platform());
    pData.contacts_.push_back(phase5_platform());
    pData.contacts_.push_back(phase6_platform());

    return pData;
}

BOOST_AUTO_TEST_CASE(quasi_static_0){
    // should be successfull in quasiStatic
// compute kinematic constraints for the right foot :
    bezier_com_traj::ContactData cDataFirst = phase0_platform();
    bezier_com_traj::ContactData cDataMid = phase1_platform();
    bezier_com_traj::ContactData cDataSecond = phase2_platform();

    ConstraintsPair kin_first = std::make_pair(cDataFirst.Kin_,cDataFirst.kin_);
    ConstraintsPair kin_second = std::make_pair(cDataSecond.Kin_,cDataSecond.kin_);
    ConstraintsPair stab = generateStabilityConstraints(*cDataMid.contactPhase_);
    std::pair<Matrix3,Vector3> Hg = computeCost();
    Vector3 init = Vector3::Zero();

    ConstraintsPair Ab_first = stackConstraints(kin_first,stab);
    bezier_com_traj::ResultData res_first = bezier_com_traj::solveIntersection(Ab_first,Hg,init);
    BOOST_CHECK(res_first.success_);

    ConstraintsPair Ab_second = stackConstraints(kin_second,stab);
    bezier_com_traj::ResultData res_second = bezier_com_traj::solveIntersection(Ab_second,Hg,init);
    BOOST_CHECK(res_second.success_);

    delete cDataFirst.contactPhase_;
    delete cDataMid.contactPhase_;
    delete cDataSecond.contactPhase_;
}

BOOST_AUTO_TEST_CASE(quasi_static_1){
    // should NOT be successfull in quasiStatic
// compute kinematic constraints for the right foot :
    bezier_com_traj::ContactData cDataFirst = phase2_platform();
    bezier_com_traj::ContactData cDataMid = phase3_platform();
    bezier_com_traj::ContactData cDataSecond = phase4_platform();

    ConstraintsPair kin_first = std::make_pair(cDataFirst.Kin_,cDataFirst.kin_);
    ConstraintsPair kin_second = std::make_pair(cDataSecond.Kin_,cDataSecond.kin_);
    ConstraintsPair stab = generateStabilityConstraints(*cDataMid.contactPhase_);
    std::pair<Matrix3,Vector3> Hg = computeCost();
    Vector3 init = Vector3::Zero();

    ConstraintsPair Ab_first = stackConstraints(kin_first,stab);
    bezier_com_traj::ResultData res_first = bezier_com_traj::solveIntersection(Ab_first,Hg,init);
    BOOST_CHECK(! res_first.success_);

    ConstraintsPair Ab_second = stackConstraints(kin_second,stab);
    bezier_com_traj::ResultData res_second = bezier_com_traj::solveIntersection(Ab_second,Hg,init);
    BOOST_CHECK(! res_second.success_);

    delete cDataFirst.contactPhase_;
    delete cDataMid.contactPhase_;
    delete cDataSecond.contactPhase_;
}

BOOST_AUTO_TEST_CASE(quasi_static_2){
    // should be successfull in quasiStatic
// compute kinematic constraints for the right foot :
    bezier_com_traj::ContactData cDataFirst = phase4_platform();
    bezier_com_traj::ContactData cDataMid = phase5_platform();
    bezier_com_traj::ContactData cDataSecond = phase6_platform();

    ConstraintsPair kin_first = std::make_pair(cDataFirst.Kin_,cDataFirst.kin_);
    ConstraintsPair kin_second = std::make_pair(cDataSecond.Kin_,cDataSecond.kin_);
    ConstraintsPair stab = generateStabilityConstraints(*cDataMid.contactPhase_);
    std::pair<Matrix3,Vector3> Hg = computeCost();
    Vector3 init = Vector3::Zero();

    ConstraintsPair Ab_first = stackConstraints(kin_first,stab);
    bezier_com_traj::ResultData res_first = bezier_com_traj::solveIntersection(Ab_first,Hg,init);
    BOOST_CHECK(res_first.success_);

    ConstraintsPair Ab_second = stackConstraints(kin_second,stab);
    bezier_com_traj::ResultData res_second = bezier_com_traj::solveIntersection(Ab_second,Hg,init);
    BOOST_CHECK(res_second.success_);


    delete cDataFirst.contactPhase_;
    delete cDataMid.contactPhase_;
    delete cDataSecond.contactPhase_;
}


BOOST_AUTO_TEST_CASE(transition_0){
    bezier_com_traj::ProblemData pData = gen_problem_data_Step0();
    VectorX Ts(3);
    Ts<<0.8,0.6,0.8;
    check_transition(pData,Ts);

}

BOOST_AUTO_TEST_CASE(transition_1){
    bezier_com_traj::ProblemData pData = gen_problem_data_Step1();
    VectorX Ts(3);
    Ts<<0.4,0.2,0.4;

    check_transition(pData,Ts);

}

BOOST_AUTO_TEST_CASE(transition_2){
    bezier_com_traj::ProblemData pData = gen_problem_data_Step2();
    VectorX Ts(3);
    Ts<<0.8,0.6,0.8;

    check_transition(pData,Ts);

}

BOOST_AUTO_TEST_SUITE_END()
