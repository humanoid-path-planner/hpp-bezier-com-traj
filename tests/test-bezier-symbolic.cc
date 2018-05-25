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
#include <spline/bezier_curve.h>

using namespace bezier_com_traj;

std::vector<point_t> generate_wps(){
    ProblemData pData;
    pData.c0_ = Vector3(0,0.5,5.);
    pData.c1_ = Vector3(2,-0.5,5.);
    pData.dc0_ = Vector3::Zero();
    pData.dc1_ = Vector3::Zero();
    pData.ddc0_ = Vector3::Zero();
    pData.ddc1_ = Vector3::Zero();
    pData.constraints_.flag_ = INIT_POS | INIT_VEL | END_VEL | END_POS;
    return computeConstantWaypoints(pData,1.);
}


bezier_wp_t::t_point_t generate_wps_symbolic(){

    const int DIM_VAR = 3;
    const int DIM_POINT = 3;

    std::vector<point_t> pts = generate_wps();
    bezier_wp_t::t_point_t wps;
    for(std::vector<point_t>::const_iterator pit = pts.begin() ; pit != pts.end() ; ++pit ){
        waypoint_t w = initwp(DIM_POINT,DIM_VAR);
        if(*pit == bezier_t::point_t::Zero()){
            w.first = MatrixXX::Identity(DIM_POINT,DIM_VAR);
        }else{
            w.second = *pit;
        }
        wps.push_back(w);
    }
    return wps;
}


point_t eval(const waypoint_t& w, const point_t& x){
    return w.first*x + w.second;
}

void vectorEqual(const point_t& a , const point_t& b, const double EPS = 1e-14){
    BOOST_CHECK((a-b).norm() < EPS);
}

BOOST_AUTO_TEST_SUITE( symbolic )

BOOST_AUTO_TEST_CASE(symbolic_eval_c){
    std::vector<point_t> pts = generate_wps();
    bezier_wp_t::t_point_t wps =  generate_wps_symbolic();
    point_t y(1,0.2,4.5);
    pts[2] = y;

    bezier_t c (pts.begin(),pts.end());
    bezier_wp_t c_sym (wps.begin(),wps.end());

    double t = 0.;
    while(t<1.){
        vectorEqual(c(t),eval(c_sym(t),y));
        t += 0.01;
    }
}

BOOST_AUTO_TEST_CASE(symbolic_eval_dc){
    std::vector<point_t> pts = generate_wps();
    bezier_wp_t::t_point_t wps =  generate_wps_symbolic();
    point_t y(1,0.2,4.5);
    pts[2] = y;

    bezier_t c (pts.begin(),pts.end());
    bezier_t dc = c.compute_derivate(1);
    bezier_wp_t c_sym (wps.begin(),wps.end());
    bezier_wp_t dc_sym = c_sym.compute_derivate(1);

    double t = 0.;
    while(t<1.){
        vectorEqual(dc(t),eval(dc_sym(t),y));
        t += 0.01;
    }
}

BOOST_AUTO_TEST_CASE(symbolic_eval_ddc){
    std::vector<point_t> pts = generate_wps();
    bezier_wp_t::t_point_t wps =  generate_wps_symbolic();
    point_t y(1,0.2,4.5);
    pts[2] = y;

    bezier_t c (pts.begin(),pts.end());
    bezier_t ddc = c.compute_derivate(2);
    bezier_wp_t c_sym (wps.begin(),wps.end());
    bezier_wp_t ddc_sym = c_sym.compute_derivate(2);

    double t = 0.;
    while(t<1.){
        vectorEqual(ddc(t),eval(ddc_sym(t),y),1e-10);
        t += 0.01;
    }
}


BOOST_AUTO_TEST_CASE(symbolic_eval_jc){
    std::vector<point_t> pts = generate_wps();
    bezier_wp_t::t_point_t wps =  generate_wps_symbolic();
    point_t y(1,0.2,4.5);
    pts[2] = y;

    bezier_t c (pts.begin(),pts.end());
    bezier_t jc = c.compute_derivate(3);
    bezier_wp_t c_sym (wps.begin(),wps.end());
    bezier_wp_t jc_sym = c_sym.compute_derivate(3);


    double t = 0.;
    while(t<1.){
        vectorEqual(jc(t),eval(jc_sym(t),y),1e-10);
        t += 0.01;
    }
}



BOOST_AUTO_TEST_SUITE_END()



