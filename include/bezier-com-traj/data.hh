/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef BEZIER_COM_TRAJ_LIB_DATA_H
#define BEZIER_COM_TRAJ_LIB_DATA_H

#include <Eigen/Dense>
#include <bezier-com-traj/config.hh>
#include <spline/bezier_curve.h>
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>
#include <vector>

namespace centroidal_dynamics
{
    struct ContactData
    {
        ContactData()
            : contactPhase_(0)
            , Kin_(Eigen::Matrix3d::Zero())
            , kin_(Vector3::Zero()) {}
       ~ContactData(){}

       centroidal_dynamics::Equilibrium* contactPhase_;
       MatrixX3 Kin_;
       VectorX kin_;
    };

    struct ProblemData
    {
        ProblemData()
            : c0_(Vector3::Zero())
            ,dc0_(Vector3::Zero())
            ,useAngularMomentum_(false) {}

        std::vector<ContactData> contacts_;
        Vector3  c0_;
        Vector3 dc0_;
        Vector3  l0_;
        bool useAngularMomentum_;
    };

    typedef Eigen::Vector3d point_t;
    typedef spline::bezier_curve  <double, double, 3, true, point_t > bezier_t;
    struct ResultData
    {
        ResultData():
            success_(false)
          , c_of_t_(0)
          , dL_of_t_(0) {}
        ~ResultData(){}

        bool success_;
        bezier_t* c_of_t_;
        bezier_t* dL_of_t_;
    };

} // end namespace centroidal_dynamics

#endif
