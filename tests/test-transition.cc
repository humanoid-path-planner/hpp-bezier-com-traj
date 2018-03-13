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
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>
#include <test_helper.hh>

#define NOMINAL_COM_HEIGHT 0.795

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


BOOST_AUTO_TEST_SUITE_END()

