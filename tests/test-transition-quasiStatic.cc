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


#define BOOST_TEST_MODULE transition-quasiStatic
#include <boost/test/included/unit_test.hpp>
#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>
#include <test_helper.hh>




BOOST_AUTO_TEST_SUITE( quasiStatic )

BOOST_AUTO_TEST_CASE(single_support){
    // check if state valid with only one contact :
    MatrixX3 normal(1,3);
    normal << 0,0,1;
    MatrixX3 position(1,3);
    position<< 0,0,0;
    std::pair<MatrixX3,VectorX> Ab = generateConstraints(normal,position,Matrix3::Identity(),Vector3::Zero());
    std::pair<Matrix3,Vector3> Hg = computeCost();
    Vector3 init = Vector3::Zero();
    bezier_com_traj::ResultData res = bezier_com_traj::solveIntersection(Ab,Hg,init);
    BOOST_CHECK(res.success_);

    // sample positions for second foot inside kinematics constraints and check for feasibility :

}

BOOST_AUTO_TEST_CASE(quasiStatic_exist){
    // sample positions for second foot inside kinematics constraints and check for feasibility :
    //Vector3 firstLeg_n(0,0,1);
    //Vector3 firstLeg_p(0,0,0);

    for(size_t i = 0 ; i < 500 ; ++i){
        MatrixX3 normal(1,3);
        MatrixX3 position(1,3);
        //normal.block<1,3>(0,0) = firstLeg_n;
        //position.block<1,3>(0,0) = firstLeg_p;
        double x = fRandom(KIN_X_MIN,KIN_X_MAX);
        double y = fRandom(KIN_Y_MIN,KIN_Y_MAX);
        normal.block<1,3>(0,0) = Vector3(0,0,1);
        position.block<1,3>(0,0) = Vector3(x,y,0.);

        std::pair<MatrixX3,VectorX> Ab = generateConstraints(normal,position,Matrix3::Identity(),Vector3::Zero());
        std::pair<Matrix3,Vector3> Hg = computeCost();
        Vector3 init = Vector3::Zero();
        bezier_com_traj::ResultData res = bezier_com_traj::solveIntersection(Ab,Hg,init);
        BOOST_CHECK(res.success_);
    }
}


BOOST_AUTO_TEST_CASE(quasiStatic_empty_upX){
    for(size_t i = 0 ; i < 500 ; ++i){
        MatrixX3 normal(1,3);
        MatrixX3 position(1,3);
        //normal.block<1,3>(0,0) = firstLeg_n;
        //position.block<1,3>(0,0) = firstLeg_p;
        double x = fRandom(KIN_X_MAX+(LX/2.)+0.001,10);
        double y = fRandom(KIN_Y_MIN,KIN_Y_MAX);
        normal.block<1,3>(0,0) = Vector3(0,0,1);
        position.block<1,3>(0,0) = Vector3(x,y,0.);

        std::pair<MatrixX3,VectorX> Ab = generateConstraints(normal,position,Matrix3::Identity(),Vector3::Zero());
        std::pair<Matrix3,Vector3> Hg = computeCost();
        Vector3 init = Vector3::Zero();
        bezier_com_traj::ResultData res = bezier_com_traj::solveIntersection(Ab,Hg,init);
        BOOST_CHECK(! res.success_);
    }
}

BOOST_AUTO_TEST_CASE(quasiStatic_empty_downX){
    for(size_t i = 0 ; i < 500 ; ++i){
        MatrixX3 normal(1,3);
        MatrixX3 position(1,3);
        //normal.block<1,3>(0,0) = firstLeg_n;
        //position.block<1,3>(0,0) = firstLeg_p;
        double x = fRandom(-10,KIN_X_MIN-(LX/2.)-0.001);
        double y = fRandom(KIN_Y_MIN,KIN_Y_MAX);
        normal.block<1,3>(0,0) = Vector3(0,0,1);
        position.block<1,3>(0,0) = Vector3(x,y,0.);

        std::pair<MatrixX3,VectorX> Ab = generateConstraints(normal,position,Matrix3::Identity(),Vector3::Zero());
        std::pair<Matrix3,Vector3> Hg = computeCost();
        Vector3 init = Vector3::Zero();
        bezier_com_traj::ResultData res = bezier_com_traj::solveIntersection(Ab,Hg,init);
        BOOST_CHECK(! res.success_);
    }
}

BOOST_AUTO_TEST_CASE(quasiStatic_empty_downY){
    for(size_t i = 0 ; i < 500 ; ++i){
        MatrixX3 normal(1,3);
        MatrixX3 position(1,3);
        //normal.block<1,3>(0,0) = firstLeg_n;
        //position.block<1,3>(0,0) = firstLeg_p;
        double x = fRandom(KIN_X_MIN,KIN_X_MAX);
        double y = fRandom(-10,KIN_Y_MIN-(LY/2.)-0.001);
        normal.block<1,3>(0,0) = Vector3(0,0,1);
        position.block<1,3>(0,0) = Vector3(x,y,0.);

        std::pair<MatrixX3,VectorX> Ab = generateConstraints(normal,position,Matrix3::Identity(),Vector3::Zero());
        std::pair<Matrix3,Vector3> Hg = computeCost();
        Vector3 init = Vector3::Zero();
        bezier_com_traj::ResultData res = bezier_com_traj::solveIntersection(Ab,Hg,init);
        BOOST_CHECK(! res.success_);
    }
}

BOOST_AUTO_TEST_CASE(quasiStatic_empty_upY){
    for(size_t i = 0 ; i < 500 ; ++i){
        MatrixX3 normal(1,3);
        MatrixX3 position(1,3);
        //normal.block<1,3>(0,0) = firstLeg_n;
        //position.block<1,3>(0,0) = firstLeg_p;
        double x = fRandom(KIN_X_MIN,KIN_X_MAX);
        double y = fRandom(KIN_Y_MAX+(LY/2.)+0.001,10);
        normal.block<1,3>(0,0) = Vector3(0,0,1);
        position.block<1,3>(0,0) = Vector3(x,y,0.);

        std::pair<MatrixX3,VectorX> Ab = generateConstraints(normal,position,Matrix3::Identity(),Vector3::Zero());
        std::pair<Matrix3,Vector3> Hg = computeCost();
        Vector3 init = Vector3::Zero();
        bezier_com_traj::ResultData res = bezier_com_traj::solveIntersection(Ab,Hg,init);
        BOOST_CHECK(! res.success_);
    }
}


BOOST_AUTO_TEST_SUITE_END()
