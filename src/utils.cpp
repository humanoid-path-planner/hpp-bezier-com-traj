/*
 * Copyright 2017, LAAS-CNRS
 * Author: Steve Tonneau
 */

#include <hpp/bezier-com-traj/utils.hh>
namespace bezier_com_traj {

waypoint_t initwp(const size_t rows, const size_t cols) {
  waypoint_t w;
  w.first = MatrixXX::Zero(rows, cols);
  w.second = VectorX::Zero(rows);
  return w;
}

waypoint_t operator+(const waypoint_t& w1, const waypoint_t& w2) {
  if (w1.second.rows() != w2.second.rows() || w1.first.rows() != w2.first.rows() || w1.first.cols() != w2.first.cols())
    throw std::runtime_error("You cannot add waypoint_t of different size.");
  return waypoint_t(w1.first + w2.first, w1.second + w2.second);
}

waypoint_t operator-(const waypoint_t& w1, const waypoint_t& w2) {
  if (w1.second.rows() != w2.second.rows() || w1.first.rows() != w2.first.rows() || w1.first.cols() != w2.first.cols())
    throw std::runtime_error("You cannot add waypoint_t of different size.");
  return waypoint_t(w1.first - w2.first, w1.second - w2.second);
}

waypoint_t operator*(const double k, const waypoint_t& w) { return waypoint_t(k * w.first, k * w.second); }

waypoint_t operator*(const waypoint_t& w, const double k) { return waypoint_t(k * w.first, k * w.second); }

template <>
waypoint9_t initwp<waypoint9_t>() {
  waypoint9_t w;
  w.first = Matrix39::Zero();
  w.second = point3_t::Zero();
  return w;
}

template <>
waypoint6_t initwp<waypoint6_t>() {
  waypoint6_t w;
  w.first = matrix6_t::Zero();
  w.second = point6_t::Zero();
  return w;
}

template <>
waypoint3_t initwp<waypoint3_t>() {
  waypoint3_t w;
  w.first = matrix3_t::Zero();
  w.second = point3_t::Zero();
  return w;
}

Matrix3 skew(point_t_tC x) {
  Matrix3 res = Matrix3::Zero();
  res(0, 1) = -x(2);
  res(0, 2) = x(1);
  res(1, 0) = x(2);
  res(1, 2) = -x(0);
  res(2, 0) = -x(1);
  res(2, 1) = x(0);
  return res;
}

std::vector<ndcurves::Bern<double> > ComputeBersteinPolynoms(const unsigned int degree) {
  std::vector<ndcurves::Bern<double> > res;
  for (unsigned int i = 0; i <= (unsigned int)degree; ++i) res.push_back(ndcurves::Bern<double>(degree, i));
  return res;
}

T_time computeDiscretizedTimeFixed(const VectorX& phaseTimings, const unsigned int pointsPerPhase) {
  T_time timeArray;
  double t = 0;
  double t_total = phaseTimings.sum();
  timeArray.push_back(std::make_pair(0., 0));
  for (int i = 0; i < phaseTimings.size(); ++i) {
    double step = (double)phaseTimings[i] / pointsPerPhase;
    for (size_t j = 0; j < pointsPerPhase; ++j) {
      t += step;
      timeArray.push_back(std::make_pair(t, i));
    }
  }
  timeArray.pop_back();
  timeArray.push_back(std::make_pair(t_total, phaseTimings.size() - 1));  // avoid numerical errors
  return timeArray;
}

T_time computeDiscretizedTime(const VectorX& phaseTimings, const double timeStep) {
  T_time timeArray;
  double t = 0;
  double currentTiming = 0.;
  for (int i = 0; i < phaseTimings.size(); ++i) {
    assert(timeStep * 2 <= phaseTimings[i] &&
           "Time step too high: should allow to contain at least 2 points per phase");
    t = currentTiming;
    currentTiming += phaseTimings[i];
    while (t < currentTiming) {
      timeArray.push_back(std::make_pair(t, i));
      t += timeStep;
    }
    timeArray.push_back(std::make_pair(currentTiming, i));
  }
  return timeArray;
}

void printQHullFile(const std::pair<MatrixXX, VectorX>& Ab, VectorX intPoint, const std::string& fileName,
                    bool clipZ) {
  std::ofstream file;
  using std::endl;
  std::string path("/local/fernbac/bench_iros18/constraints_obj/");
  path.append(fileName);
  file.open(path.c_str(), std::ios::out | std::ios::trunc);
  file << "3 1" << endl;
  file << "\t " << intPoint[0] << "\t" << intPoint[1] << "\t" << intPoint[2] << endl;
  file << "4" << endl;
  clipZ ? file << Ab.first.rows() + 2 << endl : file << Ab.first.rows() << endl;
  for (int i = 0; i < Ab.first.rows(); ++i) {
    file << "\t" << Ab.first(i, 0) << "\t" << Ab.first(i, 1) << "\t" << Ab.first(i, 2) << "\t" << -Ab.second[i] - 0.001
         << endl;
  }
  if (clipZ) {
    file << "\t" << 0 << "\t" << 0 << "\t" << 1. << "\t" << -3. << endl;
    file << "\t" << 0 << "\t" << 0 << "\t" << -1. << "\t" << -1. << endl;
  }
  file.close();
}

}  // namespace bezier_com_traj
