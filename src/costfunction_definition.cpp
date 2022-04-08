/*
 * Copyright 2017, LAAS-CNRS
 * Author: Steve Tonneau
 */

#include <hpp/bezier-com-traj/common_solve_methods.hh>
#include <hpp/bezier-com-traj/cost/costfunction_definition.hh>
#include <hpp/bezier-com-traj/waypoints/waypoints_definition.hh>

namespace bezier_com_traj {
namespace cost {
// cost : min distance between x and midPoint :
void computeCostMidPoint(const ProblemData& pData, const VectorX& /*Ts*/,
                         const double /*T*/, const T_time& /*timeArray*/,
                         MatrixXX& H, VectorX& g) {
  // cost : x' H x + 2 x g'
  Vector3 midPoint = (pData.c0_ + pData.c1_) /
                     2.;  // todo : replace it with point found by planning ??
  H.block<3, 3>(0, 0) = Matrix3::Identity();
  g.head<3>() = -midPoint;
}

// cost : min distance between end velocity and the one computed by planning
void computeCostEndVelocity(const ProblemData& pData, const VectorX& /*Ts*/,
                            const double T, const T_time& /*timeArray*/,
                            MatrixXX& H, VectorX& g) {
  if (pData.constraints_.flag_ && END_VEL)
    throw std::runtime_error(
        "Can't use computeCostEndVelocity as cost function when end velocity "
        "is a contraint");
  coefs_t v = computeFinalVelocityPoint(pData, T);
  H.block<3, 3>(0, 0) = Matrix3::Identity() * v.first * v.first;
  g.head<3>() = v.first * (v.second - pData.dc1_);
}

// TODO this is temporary.The acceleration integral can be computed analitcally
void computeCostMinAcceleration(const ProblemData& pData, const VectorX& Ts,
                                const double /*T*/, const T_time& timeArray,
                                MatrixXX& H, VectorX& g) {
  double t_total = 0.;
  for (int i = 0; i < Ts.size(); ++i) t_total += Ts[i];
  std::vector<coefs_t> wps_ddc =
      computeDiscretizedAccelerationWaypoints<point3_t>(pData, t_total,
                                                        timeArray);
  // cost : x' H x + 2 x g'
  H.block<3, 3>(0, 0) = Matrix3::Zero();
  g.head<3>() = Vector3::Zero();
  for (size_t i = 0; i < wps_ddc.size(); ++i) {
    H.block<3, 3>(0, 0) +=
        Matrix3::Identity() * wps_ddc[i].first * wps_ddc[i].first;
    g.head<3>() += wps_ddc[i].first * wps_ddc[i].second;
  }
}

typedef void (*costCompute)(const ProblemData&, const VectorX&, const double,
                            const T_time&, MatrixXX&, VectorX&);
typedef std::map<CostFunction, costCompute> T_costCompute;
static const T_costCompute costs =
    boost::assign::map_list_of(ACCELERATION, computeCostMinAcceleration)(
        DISTANCE_TRAVELED, computeCostMidPoint)(TARGET_END_VELOCITY,
                                                computeCostEndVelocity);

void genCostFunction(const ProblemData& pData, const VectorX& Ts,
                     const double T, const T_time& timeArray, MatrixXX& H,
                     VectorX& g) {
  T_costCompute::const_iterator cit = costs.find(pData.costFunction_);
  if (cit != costs.end())
    return cit->second(pData, Ts, T, timeArray, H, g);
  else {
    std::cout << "Unknown cost function" << std::endl;
    throw std::runtime_error("Unknown cost function");
  }
}

}  // namespace cost
}  // namespace bezier_com_traj
