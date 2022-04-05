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
#include <ndcurves/bezier_curve.h>

#include <boost/test/included/unit_test.hpp>
#include <hpp/bezier-com-traj/common_solve_methods.hh>
#include <hpp/bezier-com-traj/data.hh>
#include <hpp/bezier-com-traj/solve.hh>
#include <hpp/centroidal-dynamics/centroidal_dynamics.hh>

#include "test_helper.hh"

using namespace bezier_com_traj;
const double T = 1.5;

ProblemData buildPData(const centroidal_dynamics::EquilibriumAlgorithm algo =
                           centroidal_dynamics::EQUILIBRIUM_ALGORITHM_PP) {
  ProblemData pData;
  pData.c0_ = Vector3(0, 0.5, 5.);
  pData.c1_ = Vector3(2, -0.5, 5.);
  pData.dc0_ = Vector3::Zero();
  pData.dc1_ = Vector3::Zero();
  pData.ddc0_ = Vector3::Zero();
  pData.ddc1_ = Vector3::Zero();
  pData.constraints_.flag_ = INIT_POS | INIT_VEL | END_VEL | END_POS;

  MatrixX3 normals(2, 3), positions(2, 3);
  normals.block<1, 3>(0, 0) = Vector3(0, 0, 1);
  positions.block<1, 3>(0, 0) = Vector3(0, 0.1, 0);
  normals.block<1, 3>(1, 0) = Vector3(0, 0, 1);
  positions.block<1, 3>(1, 0) = Vector3(0, -0.1, 0);
  std::pair<MatrixX3, MatrixX3> contacts =
      computeRectangularContacts(normals, positions, LX, LY);
  pData.contacts_.push_back(new centroidal_dynamics::Equilibrium(
      ComputeContactCone(contacts.first, contacts.second, algo)));

  return pData;
}

std::vector<point_t> generate_wps() {
  return computeConstantWaypoints(buildPData(), T);
}

bezier_wp_t::t_point_t generate_wps_symbolic() {
  return computeConstantWaypointsSymbolic(buildPData(), T);
}

VectorX eval(const waypoint_t& w, const point_t& x) {
  return w.first * x + w.second;
}

void vectorEqual(const VectorX& a, const VectorX& b, const double EPS = 1e-14) {
  BOOST_CHECK_EQUAL(a.size(), b.size());
  BOOST_CHECK((a - b).norm() < EPS);
}

BOOST_AUTO_TEST_SUITE(symbolic)

BOOST_AUTO_TEST_CASE(symbolic_eval_c) {
  std::vector<point_t> pts = generate_wps();
  bezier_wp_t::t_point_t wps = generate_wps_symbolic();
  point_t y(1, 0.2, 4.5);
  pts[2] = y;

  bezier_t c(pts.begin(), pts.end(), 0., T);
  bezier_wp_t c_sym(wps.begin(), wps.end(), 0., T);

  double t = 0.;
  while (t < T) {
    vectorEqual(c(t), eval(c_sym(t), y));
    t += 0.01;
  }
}

BOOST_AUTO_TEST_CASE(symbolic_eval_dc) {
  std::vector<point_t> pts = generate_wps();
  bezier_wp_t::t_point_t wps = generate_wps_symbolic();
  point_t y(1, 0.2, 4.5);
  pts[2] = y;

  bezier_t c(pts.begin(), pts.end(), 0., T);
  bezier_t dc = c.compute_derivate(1);
  bezier_wp_t c_sym(wps.begin(), wps.end(), 0., T);
  bezier_wp_t dc_sym = c_sym.compute_derivate(1);

  double t = 0.;
  while (t < T) {
    vectorEqual(dc(t), eval(dc_sym(t), y));
    t += 0.01;
  }
}

BOOST_AUTO_TEST_CASE(symbolic_eval_ddc) {
  std::vector<point_t> pts = generate_wps();
  bezier_wp_t::t_point_t wps = generate_wps_symbolic();
  point_t y(1, 0.2, 4.5);
  pts[2] = y;

  bezier_t c(pts.begin(), pts.end(), 0., T);
  bezier_t ddc = c.compute_derivate(2);
  bezier_wp_t c_sym(wps.begin(), wps.end(), 0., T);
  bezier_wp_t ddc_sym = c_sym.compute_derivate(2);

  double t = 0.;
  while (t < T) {
    vectorEqual(ddc(t), eval(ddc_sym(t), y), 1e-10);
    t += 0.01;
  }
}

BOOST_AUTO_TEST_CASE(symbolic_eval_jc) {
  std::vector<point_t> pts = generate_wps();
  bezier_wp_t::t_point_t wps = generate_wps_symbolic();
  point_t y(1, 0.2, 4.5);
  pts[2] = y;

  bezier_t c(pts.begin(), pts.end(), 0., T);
  bezier_t jc = c.compute_derivate(3);
  bezier_wp_t c_sym(wps.begin(), wps.end(), 0., T);
  bezier_wp_t jc_sym = c_sym.compute_derivate(3);

  double t = 0.;
  while (t < T) {
    vectorEqual(jc(t), eval(jc_sym(t), y), 1e-10);
    t += 0.01;
  }
}

BOOST_AUTO_TEST_CASE(symbolic_split_c) {
  std::vector<point_t> pts = generate_wps();
  bezier_wp_t::t_point_t wps = generate_wps_symbolic();
  point_t y(1, 0.2, 4.5);
  pts[2] = y;

  bezier_t c(pts.begin(), pts.end(), 0., T);
  bezier_wp_t c_sym(wps.begin(), wps.end(), 0., T);

  double a, b, t, t1, t2;
  for (size_t i = 0; i < 100; ++i) {
    a = (rand() / (double)RAND_MAX) * T;
    b = (rand() / (double)RAND_MAX) * T;
    t1 = std::min(a, b);
    t2 = std::max(a, b);
    // std::cout<<"try extract between : ["<<t1<<";"<<t2<<"] "<<std::endl;
    bezier_t c_e = c.extract(t1, t2);
    bezier_wp_t c_sym_e = c_sym.extract(t1, t2);
    t = t1;
    while (t < t2) {
      vectorEqual(c_e(t), eval(c_sym_e(t), y));
      vectorEqual(c(t), eval(c_sym_e(t), y));
      t += 0.01;
    }
  }
}

BOOST_AUTO_TEST_CASE(symbolic_split_c_bench) {
  using namespace std;

  std::vector<point_t> pts = generate_wps();
  bezier_wp_t::t_point_t wps = generate_wps_symbolic();
  point_t y(1, 0.2, 4.5);
  pts[2] = y;

  bezier_wp_t c_sym(wps.begin(), wps.end(), 0., T);

  std::vector<double> values;
  for (int i = 0; i < 100000; ++i) values.push_back((double)rand() / RAND_MAX);

  clock_t s0, e0;
  std::pair<bezier_wp_t, bezier_wp_t> splitted = c_sym.split(0.5);
  s0 = clock();
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    splitted = c_sym.split(*cit);
  }
  e0 = clock();

  std::cout << "Time required to split a c curve : "
            << ((double)(e0 - s0) / CLOCKS_PER_SEC) / 100. << " ms "
            << std::endl;
}

BOOST_AUTO_TEST_CASE(symbolic_split_w) {
  bezier_wp_t::t_point_t wps = computeWwaypoints(buildPData(), T);
  point_t y(1, 0.2, 4.5);

  bezier_wp_t w(wps.begin(), wps.end(), 0., T);

  double a, b, t, t1, t2;
  for (size_t i = 0; i < 100; ++i) {
    a = (rand() / (double)RAND_MAX) * T;
    b = (rand() / (double)RAND_MAX) * T;
    t1 = std::min(a, b);
    t2 = std::max(a, b);
    // std::cout<<"try extract between : ["<<t1<<";"<<t2<<"] "<<std::endl;
    bezier_wp_t w_e = w.extract(t1, t2);
    t = t1;
    while (t < t2) {
      vectorEqual(eval(w(t), y), eval(w_e(t), y), 1e-12);
      t += 0.01;
    }
  }
}

BOOST_AUTO_TEST_CASE(symbolic_split_w_bench) {
  bezier_wp_t::t_point_t wps = computeWwaypoints(buildPData(), T);
  point_t y(1, 0.2, 4.5);

  bezier_wp_t w(wps.begin(), wps.end(), 0., T);

  std::vector<double> values;
  for (int i = 0; i < 100000; ++i) values.push_back((double)rand() / RAND_MAX);

  clock_t s0, e0;
  std::pair<bezier_wp_t, bezier_wp_t> splitted = w.split(0.5);
  s0 = clock();
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    splitted = w.split(*cit);
  }
  e0 = clock();

  std::cout << "Time required to split a w curve : "
            << ((double)(e0 - s0) / CLOCKS_PER_SEC) / 100. << " ms "
            << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
