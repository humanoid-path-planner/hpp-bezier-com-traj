/*
 * Copyright 2017, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#include <hpp/bezier-com-traj/solve.hh>
#include <hpp/bezier-com-traj/common_solve_methods.hh>
#include <limits>
#include <hpp/bezier-com-traj/waypoints/waypoints_definition.hh>

using namespace bezier_com_traj;

namespace bezier_com_traj {
typedef std::pair<double, point3_t> coefs_t;
const int DIM_POINT = 3;
// const int NUM_DISCRETIZATION = 11;
const bool verbose = false;

/**
 * @brief solveEndEffector Tries to produce a trajectory represented as a bezier curve
 * that satisfy position, velocity and acceleration constraint for the initial and final point
 * and that follow as close as possible the input trajectory
 * @param pData problem Data.
 * @param path the path to follow, the class Path must implement the operator (double t) , t \in [0,1] return a Vector3
 * that give the position on the path for a given time
 * @param T time lenght of the trajectory
 * @param timeStep time that the solver has to stop
 * @return ResultData a struct containing the resulting trajectory, if success is true.
 */
template <typename Path>
ResultDataCOMTraj solveEndEffector(const ProblemData& pData, const Path& path, const double T,
                                   const double weightDistance, bool useVelCost = true);

coefs_t initCoefs() {
  coefs_t c;
  c.first = 0;
  c.second = point3_t::Zero();
  return c;
}

// up to jerk second derivativ constraints for init, pos vel and acc constraint for goal
std::vector<bezier_t::point_t> computeConstantWaypointsInitPredef(const ProblemData& pData, double T) {
  const double n = 4;
  std::vector<bezier_t::point_t> pts;
  pts.push_back(pData.c0_);                         // c0
  pts.push_back((pData.dc0_ * T / n) + pData.c0_);  // dc0
  pts.push_back(
      (n * n * pData.c0_ - n * pData.c0_ + 2 * n * pData.dc0_ * T - 2 * pData.dc0_ * T + pData.ddc0_ * T * T) /
      (n * (n - 1)));  // ddc0 // * T because derivation make a T appear
  pts.push_back(
      (n * n * pData.c0_ - n * pData.c0_ + 3 * n * pData.dc0_ * T - 3 * pData.dc0_ * T + 3 * pData.ddc0_ * T * T) /
      (n * (n - 1)));  // j0
  //  pts.push_back((n*n*pData.c0_ - n*pData.c0_ + 4*n*pData.dc0_*T - 4*pData.dc0_ *T+ 6*pData.ddc0_*T*T)/(n*(n - 1)))
  //  ; //dj0
  //   pts.push_back((n*n*pData.c0_ - n*pData.c0_ + 5*n*pData.dc0_*T - 5*pData.dc0_ *T+ 10*pData.ddc0_*T*T)/(n*(n -
  //   1))) ; //ddj0

  // pts.push_back((n*n*pData.c1_ - n*pData.c1_ - 2*n*pData.dc1_*T + 2*pData.dc1_*T + pData.ddc1_*T*T)/(n*(n - 1))) ;
  // //ddc1 // * T ?? pts.push_back((-pData.dc1_ * T / n) + pData.c1_); // dc1
  pts.push_back(pData.c1_);  // c1
  return pts;
}

// up to jerk second derivativ constraints for goal, pos vel and acc constraint for init
std::vector<bezier_t::point_t> computeConstantWaypointsGoalPredef(const ProblemData& pData, double T) {
  const double n = 4;
  std::vector<bezier_t::point_t> pts;
  pts.push_back(pData.c0_);  // c0
  // pts.push_back((pData.dc0_ * T / n )+  pData.c0_); //dc0
  // pts.push_back((n*n*pData.c0_ - n*pData.c0_ + 2*n*pData.dc0_*T - 2*pData.dc0_*T + pData.ddc0_*T*T)/(n*(n - 1)));
  // //ddc0 // * T because derivation make a T appear

  //  pts.push_back((n*n*pData.c1_ - n*pData.c1_ - 5*n*pData.dc1_*T + 5*pData.dc1_*T + 10*pData.ddc1_*T*T)/(n*(n - 1)))
  //  ; //ddj1 pts.push_back((n*n*pData.c1_ - n*pData.c1_ - 4*n*pData.dc1_*T + 4*pData.dc1_*T +
  //  6*pData.ddc1_*T*T)/(n*(n - 1))) ; //dj1
  pts.push_back(
      (n * n * pData.c1_ - n * pData.c1_ - 3 * n * pData.dc1_ * T + 3 * pData.dc1_ * T + 3 * pData.ddc1_ * T * T) /
      (n * (n - 1)));  // j1
  pts.push_back(
      (n * n * pData.c1_ - n * pData.c1_ - 2 * n * pData.dc1_ * T + 2 * pData.dc1_ * T + pData.ddc1_ * T * T) /
      (n * (n - 1)));                                // ddc1 * T ??
  pts.push_back((-pData.dc1_ * T / n) + pData.c1_);  // dc1
  pts.push_back(pData.c1_);                          // c1
  return pts;
}

void computeConstraintsMatrix(const ProblemData& pData, const std::vector<waypoint_t>& wps_acc,
                              const std::vector<waypoint_t>& wps_vel, const VectorX& acc_bounds,
                              const VectorX& vel_bounds, MatrixXX& A, VectorX& b,
                              const std::vector<waypoint_t>& wps_jerk = std::vector<waypoint_t>(),
                              const VectorX& jerk_bounds = VectorX(DIM_POINT)) {
  assert(acc_bounds.size() == DIM_POINT && "Acceleration bounds should have the same dimension as the points");
  assert(vel_bounds.size() == DIM_POINT && "Velocity bounds should have the same dimension as the points");
  assert(jerk_bounds.size() == DIM_POINT && "Jerk bounds should have the same dimension as the points");
  const int DIM_VAR = dimVar(pData);
  int empty_acc = 0;
  int empty_vel = 0;
  int empty_jerk = 0;
  for (std::vector<waypoint_t>::const_iterator wpcit = wps_acc.begin(); wpcit != wps_acc.end(); ++wpcit) {
    if (wpcit->first.isZero(std::numeric_limits<double>::epsilon())) empty_acc++;
  }
  for (std::vector<waypoint_t>::const_iterator wpcit = wps_vel.begin(); wpcit != wps_vel.end(); ++wpcit) {
    if (wpcit->first.isZero(std::numeric_limits<double>::epsilon())) empty_vel++;
  }
  for (std::vector<waypoint_t>::const_iterator wpcit = wps_jerk.begin(); wpcit != wps_jerk.end(); ++wpcit) {
    if (wpcit->first.isZero(std::numeric_limits<double>::epsilon())) empty_jerk++;
  }

  A = MatrixXX::Zero(
      (2 * DIM_POINT * (wps_acc.size() - empty_acc + wps_vel.size() - empty_vel + wps_jerk.size() - empty_jerk)) +
          DIM_VAR,
      DIM_VAR);  // *2 because we have to put the lower and upper bound for each one, +DIM_VAR for the constraint on
                 // x[z]
  b = VectorX::Zero(A.rows());
  int i = 0;

  // upper acc bounds
  for (std::vector<waypoint_t>::const_iterator wpcit = wps_acc.begin(); wpcit != wps_acc.end(); ++wpcit) {
    if (!wpcit->first.isZero(std::numeric_limits<double>::epsilon())) {
      A.block(i * DIM_POINT, 0, DIM_POINT, DIM_VAR) = wpcit->first;
      b.segment<DIM_POINT>(i * DIM_POINT) = acc_bounds - wpcit->second;
      ++i;
    }
  }
  // lower acc bounds
  for (std::vector<waypoint_t>::const_iterator wpcit = wps_acc.begin(); wpcit != wps_acc.end(); ++wpcit) {
    if (!wpcit->first.isZero(std::numeric_limits<double>::epsilon())) {
      A.block(i * DIM_POINT, 0, DIM_POINT, DIM_VAR) = -wpcit->first;
      b.segment<DIM_POINT>(i * DIM_POINT) = acc_bounds + wpcit->second;
      ++i;
    }
  }

  // upper velocity bounds
  for (std::vector<waypoint_t>::const_iterator wpcit = wps_vel.begin(); wpcit != wps_vel.end(); ++wpcit) {
    if (!wpcit->first.isZero(std::numeric_limits<double>::epsilon())) {
      A.block(i * DIM_POINT, 0, DIM_POINT, DIM_VAR) = wpcit->first;
      b.segment<DIM_POINT>(i * DIM_POINT) = vel_bounds - wpcit->second;
      ++i;
    }
  }
  // lower velocity bounds
  for (std::vector<waypoint_t>::const_iterator wpcit = wps_vel.begin(); wpcit != wps_vel.end(); ++wpcit) {
    if (!wpcit->first.isZero(std::numeric_limits<double>::epsilon())) {
      A.block(i * DIM_POINT, 0, DIM_POINT, DIM_VAR) = -wpcit->first;
      b.segment<DIM_POINT>(i * DIM_POINT) = vel_bounds + wpcit->second;
      ++i;
    }
  }

  // upper jerk bounds
  for (std::vector<waypoint_t>::const_iterator wpcit = wps_jerk.begin(); wpcit != wps_jerk.end(); ++wpcit) {
    if (!wpcit->first.isZero(std::numeric_limits<double>::epsilon())) {
      A.block(i * DIM_POINT, 0, DIM_POINT, DIM_VAR) = wpcit->first;
      b.segment<DIM_POINT>(i * DIM_POINT) = vel_bounds - wpcit->second;
      ++i;
    }
  }
  // lower jerk bounds
  for (std::vector<waypoint_t>::const_iterator wpcit = wps_jerk.begin(); wpcit != wps_jerk.end(); ++wpcit) {
    if (!wpcit->first.isZero(std::numeric_limits<double>::epsilon())) {
      A.block(i * DIM_POINT, 0, DIM_POINT, DIM_VAR) = -wpcit->first;
      b.segment<DIM_POINT>(i * DIM_POINT) = jerk_bounds + wpcit->second;
      ++i;
    }
  }

  // test : constraint x[z] to be always higher than init[z] and goal[z].
  // TODO replace z with the direction of the contact normal ... need to change the API
  MatrixXX mxz = MatrixXX::Zero(DIM_VAR, DIM_VAR);
  size_t j = DIM_POINT - 1;
  VectorX nxz = VectorX::Zero(DIM_VAR);
  while (j < (DIM_VAR)) {
    mxz(j, j) = -1;
    nxz[j] = -std::min(pData.c0_[2], pData.c1_[2]);
    j += DIM_POINT;
  }
  A.block(i * DIM_POINT, 0, DIM_VAR, DIM_VAR) = mxz;
  b.segment(i * DIM_POINT, DIM_VAR) = nxz;

  //  std::cout<<"(i*DIM_POINT + DIM_VAR) = " << (i*DIM_POINT + DIM_VAR)<<std::endl;
  //  std::cout<<"A rows = "<<A.rows()<<std::endl;
  assert((i * DIM_POINT + DIM_VAR) == A.rows() && "Constraints matrix were not correctly initialized");
  // TEST :
  /*  A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = Matrix3::Identity();
    b.segment<DIM_POINT>(i*DIM_POINT)   = Vector3(10,10,10);
    i++;
    A.block<DIM_POINT,DIM_POINT>(i*DIM_POINT,0) = -Matrix3::Identity();
    b.segment<DIM_POINT>(i*DIM_POINT)   =  Vector3(10,10,10);*/
}

std::pair<MatrixXX, VectorX> computeDistanceCostFunction(int numPoints, const ProblemData& pData, double T,
                                                         std::vector<point3_t> pts_path) {
  assert(numPoints == pts_path.size() && "Pts_path size must be equal to numPoints");
  double step = 1. / (numPoints - 1);
  std::vector<point_t> pi = computeConstantWaypoints(pData, T);
  waypoint_t c_wp;
  MatrixXX H = MatrixXX::Zero(dimVar(pData), dimVar(pData));
  VectorX g = VectorX::Zero(dimVar(pData));
  point3_t pk;
  for (size_t i = 0; i < numPoints; ++i) {
    c_wp = evaluateCurveWaypointAtTime(pData, pi, i * step);
    pk = pts_path[i];
    //  std::cout<<"pk = "<<pk.transpose()<<std::endl;
    //  std::cout<<"coef First : "<<ckcit->first<<std::endl;
    //  std::cout<<"coef second : "<<ckcit->second.transpose()<<std::endl;
    H += (c_wp.first.transpose() * c_wp.first);
    g += ((c_wp.second - pk).transpose() * c_wp.first).transpose();
  }
  double norm = H.norm();
  H /= norm;
  g /= norm;
  return std::make_pair(H, g);
}

template <typename Path>
std::pair<MatrixXX, VectorX> computeDistanceCostFunction(int numPoints, const ProblemData& pData, double T,
                                                         const Path& path) {
  double step = 1. / (numPoints - 1);
  std::vector<point3_t> pts_path;
  for (size_t i = 0; i < numPoints; ++i) pts_path.push_back(path((double)(i * step)));
  return computeDistanceCostFunction(numPoints, pData, T, pts_path);
}

// TODO
void computeC_of_T(const ProblemData& pData, double T, ResultDataCOMTraj& res) {
  std::vector<Vector3> wps = computeConstantWaypoints(pData, T);
  if (dimVar(pData) == 3)
    wps[4] = res.x;  // FIXME : compute id from constraints
  else if (dimVar(pData) == 9) {
    wps[4] = res.x.segment<3>(0);
    wps[5] = res.x.segment<3>(3);
    wps[6] = res.x.segment<3>(6);
  } else if (dimVar(pData) == 15) {
    wps[4] = res.x.segment<3>(0);
    wps[5] = res.x.segment<3>(3);
    wps[6] = res.x.segment<3>(6);
    wps[7] = res.x.segment<3>(9);
    wps[8] = res.x.segment<3>(12);
  }
  res.c_of_t_ = bezier_t(wps.begin(), wps.end(), T);
  if (verbose) std::cout << "bezier curve created, size = " << res.c_of_t_.size_ << std::endl;
}

void computeVelCostFunctionDiscretized(int numPoints, const ProblemData& pData, double T, MatrixXX& H, VectorX& g) {
  double step = 1. / (numPoints - 1);
  std::vector<waypoint_t> cks;
  std::vector<point_t> pi = computeConstantWaypoints(pData, T);
  for (int i = 0; i < numPoints; ++i) {
    cks.push_back(evaluateVelocityCurveWaypointAtTime(pData, T, pi, i * step));
  }
  H = MatrixXX::Zero(dimVar(pData), dimVar(pData));
  g = VectorX::Zero(dimVar(pData));
  for (std::vector<waypoint_t>::const_iterator ckcit = cks.begin(); ckcit != cks.end(); ++ckcit) {
    // H+=(ckcit->first.transpose() * ckcit->first);
    // g+=ckcit->second.transpose() * ckcit->first;
    for (size_t i = 0; i < (dimVar(pData) / 3); ++i) {
      H.block<3, 3>(i * 3, i * 3) += Matrix3::Identity() * ckcit->first(0, i * 3) * ckcit->first(0, i * 3);
      g.segment<3>(i * 3) += ckcit->second.segment<3>(0) * ckcit->first(0, i * 3);
    }
  }
  // TEST : don't consider z axis for minimum acceleration cost
  // H(2,2) = 1e-6;
  // g[2] = 1e-6 ;
  // normalize :
  //  double norm=H.norm(); // because H is always diagonal
  //  H /= norm;
  //  g /= norm;
}

void computeAccelerationCostFunctionDiscretized(int numPoints, const ProblemData& pData, double T, MatrixXX& H,
                                                VectorX& g) {
  double step = 1. / (numPoints - 1);
  std::vector<waypoint_t> cks;
  std::vector<point_t> pi = computeConstantWaypoints(pData, T);
  for (int i = 0; i < numPoints; ++i) {
    cks.push_back(evaluateAccelerationCurveWaypointAtTime(pData, T, pi, i * step));
  }
  H = MatrixXX::Zero(dimVar(pData), dimVar(pData));
  g = VectorX::Zero(dimVar(pData));
  for (std::vector<waypoint_t>::const_iterator ckcit = cks.begin(); ckcit != cks.end(); ++ckcit) {
    H += (ckcit->first.transpose() * ckcit->first);
    g += ckcit->first.transpose() * ckcit->second;
  }
  // TEST : don't consider z axis for minimum acceleration cost
  // H(2,2) = 1e-6;
  // g[2] = 1e-6 ;
  // normalize :
  // double norm=H.norm(); // because H is always diagonal
  //  H /= norm;
  //  g /= norm;
}

void computeJerkCostFunctionDiscretized(int numPoints, const ProblemData& pData, double T, MatrixXX& H, VectorX& g) {
  double step = 1. / (numPoints - 1);

  std::vector<waypoint_t> cks;
  std::vector<point_t> pi = computeConstantWaypoints(pData, T);
  for (int i = 0; i < numPoints; ++i) {
    cks.push_back(evaluateJerkCurveWaypointAtTime(pData, T, pi, i * step));
  }
  H = MatrixXX::Zero(dimVar(pData), dimVar(pData));
  g = VectorX::Zero(dimVar(pData));
  for (std::vector<waypoint_t>::const_iterator ckcit = cks.begin(); ckcit != cks.end(); ++ckcit) {
    H += (ckcit->first.transpose() * ckcit->first);
    g += ckcit->first.transpose() * ckcit->second;
  }
  // TEST : don't consider z axis for minimum acceleration cost
  // H(2,2) = 1e-6;
  // g[2] = 1e-6 ;
  // normalize :
  // double norm=H.norm(); // because H is always diagonal
  // H /= norm;
  // g /= norm;
}

std::pair<MatrixXX, VectorX> computeEndEffectorConstraints(const ProblemData& pData, const double T,
                                                           std::vector<bezier_t::point_t> pi) {
  std::vector<waypoint_t> wps_jerk = computeJerkWaypoints(pData, T, pi);
  std::vector<waypoint_t> wps_acc = computeAccelerationWaypoints(pData, T, pi);
  std::vector<waypoint_t> wps_vel = computeVelocityWaypoints(pData, T, pi);
  // stack the constraint for each waypoint :
  Vector3 jerk_bounds(10000, 10000, 10000);  // TODO : read it from somewhere (ProblemData ?)
  Vector3 acc_bounds(10000, 10000, 10000);
  Vector3 vel_bounds(10000, 10000, 10000);
  MatrixXX A;
  VectorX b;
  computeConstraintsMatrix(pData, wps_acc, wps_vel, acc_bounds, vel_bounds, A, b, wps_jerk, jerk_bounds);
  return std::make_pair(A, b);
}

template <typename Path>
std::pair<MatrixXX, VectorX> computeEndEffectorCost(const ProblemData& pData, const Path& path, const double T,
                                                    const double weightDistance, bool /*useVelCost*/,
                                                    std::vector<bezier_t::point_t> pi) {
  assert(weightDistance >= 0. && weightDistance <= 1. && "WeightDistance must be between 0 and 1");
  double weightSmooth = 1. - weightDistance;
  const int DIM_VAR = dimVar(pData);
  // compute distance cost function (discrete integral under the curve defined by 'path')
  MatrixXX H;
  VectorX g;
  std::pair<MatrixXX, VectorX> Hg_smooth, Hg_rrt;

  if (weightDistance > 0)
    Hg_rrt = computeDistanceCostFunction<Path>(50, pData, T, path);
  else {
    Hg_rrt.first = MatrixXX::Zero(DIM_VAR, DIM_VAR);
    Hg_rrt.second = VectorX::Zero(DIM_VAR);
  }

  Hg_smooth = computeVelocityCost(pData, T, pi);

  /*  std::cout<<"End eff H_rrt = "<<std::endl<<H_rrt<<std::endl;
    std::cout<<"End eff g_rrt = "<<std::endl<<g_rrt<<std::endl;
    std::cout<<"End eff H_acc = "<<std::endl<<H_acc<<std::endl;
    std::cout<<"End eff g_acc = "<<std::endl<<g_acc<<std::endl;
*/
  // add the costs :
  H = MatrixXX::Zero(DIM_VAR, DIM_VAR);
  g = VectorX::Zero(DIM_VAR);
  H = weightSmooth * (Hg_smooth.first) + weightDistance * Hg_rrt.first;
  g = weightSmooth * (Hg_smooth.second) + weightDistance * Hg_rrt.second;
  // H = Hg_smooth.first;
  //  g = Hg_smooth.second;

  return std::make_pair(H, g);
}

template <typename Path>
ResultDataCOMTraj solveEndEffector(const ProblemData& pData, const Path& path, const double T,
                                   const double weightDistance, bool useVelCost) {
  if (verbose) std::cout << "solve end effector, T = " << T << std::endl;
  std::vector<bezier_t::point_t> pi = computeConstantWaypoints(pData, T);
  std::pair<MatrixXX, VectorX> Ab = computeEndEffectorConstraints(pData, T, pi);
  std::pair<MatrixXX, VectorX> Hg = computeEndEffectorCost(pData, path, T, weightDistance, useVelCost, pi);
  if (verbose) {
    std::cout << "End eff A = " << std::endl << Ab.first << std::endl;
    std::cout << "End eff b = " << std::endl << Ab.second << std::endl;
    std::cout << "End eff H = " << std::endl << Hg.first << std::endl;
    std::cout << "End eff g = " << std::endl << Hg.second << std::endl;
    std::cout << "Dim Var = " << dimVar(pData) << std::endl;
    std::cout << "Dim H   = " << Hg.first.rows() << " x " << Hg.first.cols() << std::endl;
    std::cout << "Dim g   = " << Hg.second.rows() << std::endl;
    std::cout << "Dim A   = " << Ab.first.rows() << " x " << Ab.first.cols() << std::endl;
    std::cout << "Dim b   = " << Ab.first.rows() << std::endl;
  }

  VectorX init = VectorX(dimVar(pData));
  // init = (pData.c0_ + pData.c1_)/2.;
  // init =pData.c0_;
  if (verbose) std::cout << "Init = " << std::endl << init.transpose() << std::endl;

  ResultData resQp = solve(Ab, Hg, init);

  ResultDataCOMTraj res;
  if (resQp.success_) {
    res.success_ = true;
    res.x = resQp.x;
    // computeRealCost(pData, res);
    computeC_of_T(pData, T, res);
    // computedL_of_T(pData,Ts,res);
  }
  if (verbose) {
    std::cout << "Solved, success = " << res.success_ << " x = " << res.x.transpose() << std::endl;
    std::cout << "Final cost : " << resQp.cost_ << std::endl;
  }
  return res;
}

}  // namespace bezier_com_traj
