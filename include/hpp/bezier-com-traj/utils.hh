/*
 * Copyright 2018, LAAS-CNRS
 * Author: Steve Tonneau
 */

#ifndef BEZIER_COM_TRAJ_LIB_UTILS_H
#define BEZIER_COM_TRAJ_LIB_UTILS_H

#include <hpp/bezier-com-traj/local_config.hh>
#include <hpp/bezier-com-traj/definitions.hh>
#include <hpp/bezier-com-traj/flags.hh>

#include <Eigen/Dense>

#include <vector>

namespace bezier_com_traj {

template <typename T>
T initwp();
waypoint_t initwp(const size_t rows, const size_t cols);
waypoint_t operator+(const waypoint_t& w1, const waypoint_t& w2);
waypoint_t operator-(const waypoint_t& w1, const waypoint_t& w2);
waypoint_t operator*(const double k, const waypoint_t& w);
waypoint_t operator*(const waypoint_t& w, const double k);

struct waypoint_t {
  MatrixXX first;
  VectorX second;

  waypoint_t() : first(MatrixXX()), second(VectorX()) {}

  waypoint_t(MatrixXX A, VectorX b) : first(A), second(b) {}

  static waypoint_t Zero(size_t dim) { return initwp(dim, dim); }

  size_t size() const{return second.size();}

  bool isApprox(const waypoint_t& other, const value_type prec = Eigen::NumTraits<value_type>::dummy_precision()) const{
    return first.isApprox(other.first,prec) && second.isApprox(other.second,prec);
  }

  bool operator==(const waypoint_t& other) const{ return isApprox(other); }

  bool operator!=(const waypoint_t& other) const{ return !(*this == other); }
};

/**
 * @brief Compute the Bernstein polynoms for a given degree
 * @param degree required degree
 * @return the bernstein polynoms
 */
BEZIER_COM_TRAJ_DLLAPI std::vector<curves::Bern<double> > ComputeBersteinPolynoms(const unsigned int degree);

/**
 * @brief given the constraints of the problem, and a set of waypoints, return
 * the bezier curve corresponding
 * @param pData problem data
 * @param T total trajectory time
 * @param pis list of waypoints
 * @return the bezier curve
 */
template <typename Bezier, typename Point>
BEZIER_COM_TRAJ_DLLAPI Bezier computeBezierCurve(const ConstraintFlag& flag, const double T,
                                                 const std::vector<Point>& pi, const Point& x);

/**
 * @brief computeDiscretizedTime build an array of discretized points in time,
 * such that there is the same number of point in each phase. Doesn't contain t=0,
 * is of size pointsPerPhase*phaseTimings.size()
 * @param phaseTimings
 * @param pointsPerPhase
 * @return
 */
T_time computeDiscretizedTimeFixed(const VectorX& phaseTimings, const unsigned int pointsPerPhase);

/**
 * @brief computeDiscretizedTime build an array of discretized points in time,
 * given the timestep. Doesn't contain t=0,
 * is of size pointsPerPhase*phaseTimings.size()
 * @param phaseTimings
 * @param timeStep
 * @return */
T_time computeDiscretizedTime(const VectorX& phaseTimings, const double timeStep);

/**
 * @brief write a polytope describe by A x <= b linear constraints in
 * a given filename
 * @return the bernstein polynoms
 */
void printQHullFile(const std::pair<MatrixXX, VectorX>& Ab, VectorX intPoint, const std::string& fileName,
                    bool clipZ = false);

/**
 * @brief skew symmetric matrix
 */
BEZIER_COM_TRAJ_DLLAPI Matrix3 skew(point_t_tC x);

/**
 * @brief normalize inequality constraints
 */
int Normalize(Ref_matrixXX A, Ref_vectorX b);

}  // end namespace bezier_com_traj

template <typename Bezier, typename Point>
Bezier bezier_com_traj::computeBezierCurve(const ConstraintFlag& flag, const double T, const std::vector<Point>& pi,
                                           const Point& x) {
  std::vector<Point> wps;
  size_t i = 0;
  if (flag & INIT_POS) {
    wps.push_back(pi[i]);
    i++;
    if (flag & INIT_VEL) {
      wps.push_back(pi[i]);
      i++;
      if (flag & INIT_ACC) {
        wps.push_back(pi[i]);
        i++;
      }
    }
  }
  wps.push_back(x);
  i++;
  if (flag & (END_VEL) && !(flag & (END_POS))) {
    wps.push_back(x);
    i++;
  } else {
    if (flag & END_ACC) {
      assert(flag & END_VEL && "You cannot constrain final acceleration if final velocity is not constrained.");
      wps.push_back(pi[i]);
      i++;
    }
    if (flag & END_VEL) {
      assert(flag & END_POS && "You cannot constrain final velocity if final position is not constrained.");
      wps.push_back(pi[i]);
      i++;
    }
    if (flag & END_POS) {
      wps.push_back(pi[i]);
      i++;
    }
  }
  return Bezier(wps.begin(), wps.end(), 0.,T);
}

#endif
