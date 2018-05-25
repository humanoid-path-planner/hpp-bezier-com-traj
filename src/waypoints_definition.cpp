/*
 * Copyright 2018, LAAS-CNRS
 * Author: Pierre Fernbach
 */

#ifndef BEZIER_COM_TRAJ_WP_DEF_H
#define BEZIER_COM_TRAJ_WP_DEF_H

#include <bezier-com-traj/data.hh>
#include <bezier-com-traj/waypoints/waypoints_definition.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_dc1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_dc1_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_dc1_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_j0_j1_ddc1_dc1_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1.hh>
#include <bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1.hh>

#include "boost/assign.hpp"



namespace bezier_com_traj{
/**
  * This file is used to choose the correct expressions of the curves waypoints, depending on the options set in ProblemData.constraints
  */


int dimVar(const ProblemData& pData){
    if(pData.constraints_.flag_ & FIVE_FREE_VAR)
        return 15;
    else if(pData.constraints_.flag_ & FOUR_FREE_VAR)
        return 12;
    else if(pData.constraints_.flag_ & THREE_FREE_VAR)
        return 9;
    else if(pData.constraints_.flag_ & TWO_FREE_VAR)
        return 6;
    else
        return 3;
}

typedef std::pair<double,point3_t> coefs_t;
typedef coefs_t (*evalCurveAtTime) (const std::vector<point_t>& pi,double t);
typedef std::map<ConstraintFlag,evalCurveAtTime > T_evalCurveAtTime;
typedef T_evalCurveAtTime::const_iterator         CIT_evalCurveAtTime;
static const T_evalCurveAtTime evalCurveAtTimes = boost::assign::map_list_of
        (c0_dc0_c1::flag                     , c0_dc0_c1::evaluateCurveAtTime)
        (c0_dc0_dc1::flag                    , c0_dc0_dc1::evaluateCurveAtTime)
        (c0_dc0_dc1_c1::flag                 , c0_dc0_dc1_c1::evaluateCurveAtTime)
        (c0_dc0_ddc0_c1::flag                , c0_dc0_ddc0_c1::evaluateCurveAtTime)
        (c0_dc0_ddc0_dc1_c1::flag            , c0_dc0_ddc0_dc1_c1::evaluateCurveAtTime)
        (c0_dc0_ddc0_ddc1_dc1_c1::flag       , c0_dc0_ddc0_ddc1_dc1_c1::evaluateCurveAtTime)
        (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::evaluateCurveAtTime);

/** @brief evaluateCurveAtTime compute the expression of the point on the curve c at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
// TODOin C++ 10, all these methods could be just one function :)
 coefs_t evaluateCurveAtTime(const ProblemData& pData, const std::vector<point_t>& pi,double t)
{
    CIT_evalCurveAtTime cit = evalCurveAtTimes.find(pData.constraints_.flag_);
    if(cit != evalCurveAtTimes.end())
        return cit->second(pi,t);
    else
    {
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}



typedef coefs_t (*evalAccCurveAtTime) (const std::vector<point_t>& pi,double T,double t);
typedef std::map<ConstraintFlag,evalAccCurveAtTime > T_evalAccCurveAtTime;
typedef T_evalAccCurveAtTime::const_iterator         CIT_evalAccCurveAtTime;
static const T_evalAccCurveAtTime evalAccCurveAtTimes = boost::assign::map_list_of
        (c0_dc0_c1::flag                , c0_dc0_c1::evaluateAccelerationCurveAtTime)
        (c0_dc0_dc1::flag               , c0_dc0_dc1::evaluateAccelerationCurveAtTime)
        (c0_dc0_dc1_c1::flag            , c0_dc0_dc1_c1::evaluateAccelerationCurveAtTime)
        (c0_dc0_ddc0_c1::flag           , c0_dc0_ddc0_c1::evaluateAccelerationCurveAtTime)
        (c0_dc0_ddc0_dc1_c1::flag       , c0_dc0_ddc0_dc1_c1::evaluateAccelerationCurveAtTime)
        (c0_dc0_ddc0_ddc1_dc1_c1::flag  , c0_dc0_ddc0_ddc1_dc1_c1::evaluateAccelerationCurveAtTime)
        (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::evaluateAccelerationCurveAtTime);


/** @brief evaluateAccelerationCurveAtTime compute the expression of the point on the curve ddc at t, defined by the waypoint pi and one free waypoint (x)
     * @param pi constant waypoints of the curve
     * @param t param (normalized !)
     * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
     */
 coefs_t evaluateAccelerationCurveAtTime(const ProblemData& pData, const std::vector<point_t>& pi,double T,double t)
{
    CIT_evalAccCurveAtTime cit = evalAccCurveAtTimes.find(pData.constraints_.flag_);
    if(cit != evalAccCurveAtTimes.end())
        return cit->second(pi,T,t);
    else
    {
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}


 typedef waypoint_t (*evalCurveWaypointAtTime) (const std::vector<point_t>& pi,double t);
 typedef std::map<ConstraintFlag,evalCurveWaypointAtTime > T_evalCurveWaypointAtTime;
 typedef T_evalCurveWaypointAtTime::const_iterator         CIT_evalCurveWaypointAtTime;
 static const T_evalCurveWaypointAtTime evalCurveWaypointAtTimes = boost::assign::map_list_of
         (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::evaluateCurveWaypointAtTime)
         (c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::evaluateCurveWaypointAtTime)
         (c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::evaluateCurveWaypointAtTime);


 /** @brief evaluateCurveAtTime compute the expression of the point on the curve c at t, defined by the waypoint pi and one free waypoint (x)
  * @param pi constant waypoints of the curve
  * @param t param (normalized !)
  * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
  */
 // TODOin C++ 10, all these methods could be just one function :)
  waypoint_t evaluateCurveWaypointAtTime(const ProblemData& pData, const std::vector<point_t>& pi, double t)
 {
     CIT_evalCurveWaypointAtTime cit = evalCurveWaypointAtTimes.find(pData.constraints_.flag_);
     if(cit != evalCurveWaypointAtTimes.end())
         return cit->second(pi,t);
     else
     {
         std::cout<<"Current constraints set are not implemented"<<std::endl;
         throw std::runtime_error("Current constraints set are not implemented");
     }
 }
  typedef waypoint_t (*evalVelCurveWaypointAtTime) (const std::vector<point_t>& pi, const double T,double t);
  typedef std::map<ConstraintFlag,evalVelCurveWaypointAtTime > T_evalVelCurveWaypointAtTime;
  typedef T_evalVelCurveWaypointAtTime::const_iterator         CIT_evalVelCurveWaypointAtTime;
  static const T_evalVelCurveWaypointAtTime evalVelCurveWaypointAtTimes = boost::assign::map_list_of
          (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::evaluateVelocityCurveWaypointAtTime)
          (c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::evaluateVelocityCurveWaypointAtTime)
          (c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::evaluateVelocityCurveWaypointAtTime);

  /** @brief evaluateCurveAtTime compute the expression of the point on the curve c at t, defined by the waypoint pi and one free waypoint (x)
   * @param pi constant waypoints of the curve
   * @param t param (normalized !)
   * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
   */
  // TODOin C++ 10, all these methods could be just one function :)
   waypoint_t evaluateVelocityCurveWaypointAtTime(const ProblemData& pData,  const double T,const std::vector<point_t>& pi,double t)
  {
      CIT_evalVelCurveWaypointAtTime cit = evalVelCurveWaypointAtTimes.find(pData.constraints_.flag_);
      if(cit != evalVelCurveWaypointAtTimes.end())
          return cit->second(pi,T,t);
      else
      {
          std::cout<<"Current constraints set are not implemented"<<std::endl;
          throw std::runtime_error("Current constraints set are not implemented");
      }
  }
   typedef waypoint_t (*evalAccCurveWaypointAtTime) (const std::vector<point_t>& pi, const double T,double t);
   typedef std::map<ConstraintFlag,evalAccCurveWaypointAtTime > T_evalAccCurveWaypointAtTime;
   typedef T_evalAccCurveWaypointAtTime::const_iterator         CIT_evalAccCurveWaypointAtTime;
   static const T_evalAccCurveWaypointAtTime evalAccCurveWaypointAtTimes = boost::assign::map_list_of
           (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::evaluateAccelerationCurveWaypointAtTime)
           (c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::evaluateAccelerationCurveWaypointAtTime)
           (c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::evaluateAccelerationCurveWaypointAtTime);

   /** @brief evaluateCurveAtTime compute the expression of the point on the curve c at t, defined by the waypoint pi and one free waypoint (x)
    * @param pi constant waypoints of the curve
    * @param t param (normalized !)
    * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
    */
   // TODOin C++ 10, all these methods could be just one function :)
    waypoint_t evaluateAccelerationCurveWaypointAtTime(const ProblemData& pData, const double T, const std::vector<point_t>& pi,double t)
   {
       CIT_evalAccCurveWaypointAtTime cit = evalAccCurveWaypointAtTimes.find(pData.constraints_.flag_);
       if(cit != evalAccCurveWaypointAtTimes.end())
           return cit->second(pi,T,t);
       else
       {
           std::cout<<"Current constraints set are not implemented"<<std::endl;
           throw std::runtime_error("Current constraints set are not implemented");
       }
   }
typedef waypoint_t (*evalJerkCurveWaypointAtTime) (const std::vector<point_t>& pi, const double T,double t);
typedef std::map<ConstraintFlag,evalJerkCurveWaypointAtTime > T_evalJerkCurveWaypointAtTime;
typedef T_evalJerkCurveWaypointAtTime::const_iterator         CIT_evalJerkCurveWaypointAtTime;
static const T_evalJerkCurveWaypointAtTime evalJerkCurveWaypointAtTimes = boost::assign::map_list_of
        (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::evaluateJerkCurveWaypointAtTime)
        (c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::evaluateJerkCurveWaypointAtTime)
        (c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::evaluateJerkCurveWaypointAtTime);


/** @brief evaluateCurveAtTime compute the expression of the point on the curve c at t, defined by the waypoint pi and one free waypoint (x)
 * @param pi constant waypoints of the curve
 * @param t param (normalized !)
 * @return the expression of the waypoint such that wp.first . x + wp.second = point on curve
 */
// TODOin C++ 10, all these methods could be just one function :)
 waypoint_t evaluateJerkCurveWaypointAtTime(const ProblemData& pData, const double T, const std::vector<point_t>& pi,double t)
{
    CIT_evalJerkCurveWaypointAtTime cit = evalJerkCurveWaypointAtTimes.find(pData.constraints_.flag_);
    if(cit != evalJerkCurveWaypointAtTimes.end())
        return cit->second(pi,T,t);
    else
    {
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}

typedef std::vector<point_t> (*compConsWp) (const ProblemData& pData,double T);
typedef std::map<ConstraintFlag,compConsWp > T_compConsWp;
typedef T_compConsWp::const_iterator         CIT_compConsWp;
static const T_compConsWp compConsWps = boost::assign::map_list_of
        (c0_dc0_c1::flag                , c0_dc0_c1::computeConstantWaypoints)
        (c0_dc0_dc1::flag               , c0_dc0_dc1::computeConstantWaypoints)
        (c0_dc0_dc1_c1::flag            , c0_dc0_dc1_c1::computeConstantWaypoints)
        (c0_dc0_ddc0_c1::flag           , c0_dc0_ddc0_c1::computeConstantWaypoints)
        (c0_dc0_ddc0_dc1_c1::flag       , c0_dc0_ddc0_dc1_c1::computeConstantWaypoints)
        (c0_dc0_ddc0_ddc1_dc1_c1::flag  , c0_dc0_ddc0_ddc1_dc1_c1::computeConstantWaypoints)
        (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::computeConstantWaypoints)
        (c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::computeConstantWaypoints)
        (c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::computeConstantWaypoints);

/**
 * @brief computeConstantWaypoints compute the constant waypoints of c(t) defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
 std::vector<point_t> computeConstantWaypoints(const ProblemData& pData,double T)
{
    CIT_compConsWp cit = compConsWps.find(pData.constraints_.flag_);
    if(cit != compConsWps.end())
        return cit->second(pData,T);
    else
    {
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}

bezier_wp_t::t_point_t computeConstantWaypointsSymbolic(const ProblemData& pData,double T)
{
    const int DIM_POINT = 3; // FIXME : always true ??
    const int DIM_VAR = dimVar(pData);
    std::vector<point_t> pts = computeConstantWaypoints(pData,T);
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


 typedef std::vector<waypoint_t> (*compVelWp) (const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi);
 typedef std::map<ConstraintFlag,compVelWp > T_compVelWp;
 typedef T_compVelWp::const_iterator         CIT_compVelWp;
 static const T_compVelWp compVelWps = boost::assign::map_list_of
         (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::computeVelocityWaypoints)
         (c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::computeVelocityWaypoints)
         (c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::computeVelocityWaypoints);


 /**
  * @brief computeConstantWaypoints compute the constant waypoints of c(t) defined by the constraints on initial and final states
  * @param pData
  * @param T
  * @return
  */
  std::vector<waypoint_t> computeVelocityWaypoints(const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi)
 {
     CIT_compVelWp cit = compVelWps.find(pData.constraints_.flag_);
     if(cit != compVelWps.end())
         return cit->second(pData,T,pi);
     else
     {
         std::cout<<"Current constraints set are not implemented"<<std::endl;
         throw std::runtime_error("Current constraints set are not implemented");
     }
 }


  typedef std::vector<waypoint_t> (*compAccWp) (const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi);
  typedef std::map<ConstraintFlag,compAccWp > T_compAccWp;
  typedef T_compAccWp::const_iterator         CIT_compAccWp;
  static const T_compAccWp compAccWps = boost::assign::map_list_of
          (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::computeAccelerationWaypoints)
          (c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::computeAccelerationWaypoints)
          (c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::computeAccelerationWaypoints);


  /**
   * @brief computeConstantWaypoints compute the constant waypoints of c(t) defined by the constraints on initial and final states
   * @param pData
   * @param T
   * @return
   */
   std::vector<waypoint_t> computeAccelerationWaypoints(const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi )
  {
      CIT_compAccWp cit = compAccWps.find(pData.constraints_.flag_);
      if(cit != compAccWps.end())
          return cit->second(pData,T,pi);
      else
      {
          std::cout<<"Current constraints set are not implemented"<<std::endl;
          throw std::runtime_error("Current constraints set are not implemented");
      }
  }


typedef std::vector<waypoint_t> (*compJerkWp) (const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi);
typedef std::map<ConstraintFlag,compJerkWp > T_compJerkWp;
typedef T_compJerkWp::const_iterator         CIT_compJerkWp;
static const T_compJerkWp compJerkWps = boost::assign::map_list_of
       (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::computeJerkWaypoints)
        (c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::computeJerkWaypoints)
        (c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::computeJerkWaypoints);



/**
* @brief computeConstantWaypoints compute the constant waypoints of c(t) defined by the constraints on initial and final states
* @param pData
* @param T
* @return
*/
std::vector<waypoint_t> computeJerkWaypoints(const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi )
{
   CIT_compJerkWp cit = compJerkWps.find(pData.constraints_.flag_);
   if(cit != compJerkWps.end())
       return cit->second(pData,T,pi);
   else
   {
       std::cout<<"Current constraints set are not implemented"<<std::endl;
       throw std::runtime_error("Current constraints set are not implemented");
   }
}


typedef  bezier_wp_t::t_point_t (*compWp) (const ProblemData& pData,double T);
typedef std::map<ConstraintFlag,compWp > T_compWp;
typedef T_compWp::const_iterator         CIT_compWp;
static const T_compWp compWps = boost::assign::map_list_of
        (c0_dc0_c1::flag                , c0_dc0_c1::computeWwaypoints)
        (c0_dc0_dc1::flag               , c0_dc0_dc1::computeWwaypoints)
        (c0_dc0_dc1_c1::flag            , c0_dc0_dc1_c1::computeWwaypoints)
        (c0_dc0_ddc0_c1::flag           , c0_dc0_ddc0_c1::computeWwaypoints)
        (c0_dc0_ddc0_dc1_c1::flag       , c0_dc0_ddc0_dc1_c1::computeWwaypoints)
        (c0_dc0_ddc0_ddc1_dc1_c1::flag  , c0_dc0_ddc0_ddc1_dc1_c1::computeWwaypoints)
        (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::computeWwaypoints);


/**
 * @brief computeWwaypoints compute the constant waypoints of w(t) defined by the constraints on initial and final states
 * @param pData
 * @param T
 * @return
 */
 bezier_wp_t::t_point_t computeWwaypoints(const ProblemData& pData,double T)
{
    CIT_compWp cit = compWps.find(pData.constraints_.flag_);
    if(cit != compWps.end())
        return cit->second(pData,T);
    else
    {
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}

typedef coefs_t (*compFinalVelP) (const ProblemData& pData,double T);
typedef std::map<ConstraintFlag,compFinalVelP > T_compFinalVelP;
typedef T_compFinalVelP::const_iterator         CIT_compFinalVelP;
static const T_compFinalVelP compFinalVelPs = boost::assign::map_list_of
        (c0_dc0_c1::flag                , c0_dc0_c1::computeFinalVelocityPoint)
        (c0_dc0_dc1::flag               , c0_dc0_dc1::computeFinalVelocityPoint)
        (c0_dc0_dc1_c1::flag            , c0_dc0_dc1_c1::computeFinalVelocityPoint)
        (c0_dc0_ddc0_c1::flag           , c0_dc0_ddc0_c1::computeFinalVelocityPoint)
        (c0_dc0_ddc0_dc1_c1::flag       , c0_dc0_ddc0_dc1_c1::computeFinalVelocityPoint)
        (c0_dc0_ddc0_ddc1_dc1_c1::flag  , c0_dc0_ddc0_ddc1_dc1_c1::computeFinalVelocityPoint);


 coefs_t computeFinalVelocityPoint(const ProblemData& pData,double T)
{
    CIT_compFinalVelP cit = compFinalVelPs.find(pData.constraints_.flag_);
    if(cit != compFinalVelPs.end())
        return cit->second(pData,T);
    else
    {
        std::cout<<"Current constraints set are not implemented"<<std::endl;
        throw std::runtime_error("Current constraints set are not implemented");
    }
}




typedef std::pair<MatrixXX,VectorX> (*compVelCost) (const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi);
typedef std::map<ConstraintFlag,compVelCost > T_compVelCost;
typedef T_compVelCost::const_iterator         CIT_compVelCost;
static const T_compVelCost compVelCosts = boost::assign::map_list_of
       (c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_j1_ddc1_dc1_c1::computeVelocityCost)
        (c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1::computeVelocityCost)
        (c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::flag , c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1::computeVelocityCost);



/**
* @brief computeVelocityCost the matrices H and g defining a cost that minimise the integral of the squared velocity
* @param pData
* @param T
* @return
*/
std::pair<MatrixXX,VectorX> computeVelocityCost(const ProblemData& pData,double T,std::vector<bezier_t::point_t> pi )
{
   CIT_compVelCost cit = compVelCosts.find(pData.constraints_.flag_);
   if(cit != compVelCosts.end())
       return cit->second(pData,T,pi);
   else
   {
       std::cout<<"Current constraints set are not implemented"<<std::endl;
       throw std::runtime_error("Current constraints set are not implemented");
   }
}


}

#endif
