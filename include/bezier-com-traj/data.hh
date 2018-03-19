/*
 * Copyright 2018, LAAS-CNRS
 * Author: Steve Tonneau
 */

#ifndef BEZIER_COM_TRAJ_LIB_DATA_H
#define BEZIER_COM_TRAJ_LIB_DATA_H

#include <bezier-com-traj/config.hh>
#include <bezier-com-traj/flags.hh>
#include <bezier-com-traj/definitions.hh>

#include <spline/bezier_curve.h>
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>
#include <Eigen/Dense>

#include <vector>

namespace bezier_com_traj
{
    struct BEZIER_COM_TRAJ_DLLAPI ContactData
    {
        ContactData()
            : contactPhase_(0)
            , Kin_(Eigen::Matrix3d::Zero())
            , kin_(VectorX::Zero(0))
            , Ang_(Eigen::Matrix3d::Zero())
            , ang_(VectorX::Zero(0)) {}
       ~ContactData(){}

       centroidal_dynamics::Equilibrium* contactPhase_;
       MatrixX3 Kin_;
       VectorX kin_;
       MatrixX3 Ang_;
       VectorX ang_;
    };

    struct BEZIER_COM_TRAJ_DLLAPI Constraints
    {
        Constraints()
            : flag_(INIT_POS | INIT_VEL | END_VEL | END_POS)
            , constraintAcceleration_(true)
            , maxAcceleration_(5.)
            , reduce_h_(1e-4) {}

        Constraints(ConstraintFlag flag)
            : flag_(flag)
            , constraintAcceleration_(true)
            , maxAcceleration_(5.)
            , reduce_h_(1e-4) {}

        ~Constraints(){}

        ConstraintFlag flag_;
        bool constraintAcceleration_;
        double maxAcceleration_;
        double reduce_h_;
    };


    struct BEZIER_COM_TRAJ_DLLAPI ProblemData
    {
        ProblemData()
            : c0_ (Vector3::Zero())
            ,dc0_ (Vector3::Zero())
            ,ddc0_(Vector3::Zero())
            , c1_ (Vector3::Zero())
            ,dc1_ (Vector3::Zero())
            ,ddc1_(Vector3::Zero())
            ,useAngularMomentum_(false)
            ,costFunction_(ACCELERATION) {}

        std::vector<ContactData> contacts_;
        Vector3  c0_,dc0_,ddc0_,c1_,dc1_,ddc1_;
        Vector3  l0_;
        bool useAngularMomentum_;
        Constraints constraints_;
        CostFunction costFunction_;
    };

    struct BEZIER_COM_TRAJ_DLLAPI ResultData
    {
        ResultData():
            success_(false)
          , cost_(-1.)
          , x(VectorX::Zero(0)){}

        ResultData(const bool success, const double cost, Cref_vectorX x ):
            success_(success)
          , cost_(cost)
          , x(x){}

        ResultData(const ResultData& other):
            success_(other.success_)
          , cost_(other.cost_)
          , x(other.x){}
        ~ResultData(){}

        ResultData& operator=(const ResultData& other)
        {
            success_= (other.success_);
            cost_ = (other.cost_);
            x = (other.x);
            return *this;
        }
        bool success_;
        double cost_;
        VectorX x;
    };

    struct BEZIER_COM_TRAJ_DLLAPI ResultDataCOMTraj : public ResultData
    {
        ResultDataCOMTraj():
            ResultData()
          , c_of_t_(bezier_t::zero())
          , dL_of_t_(bezier_t::zero())
          , dc1_(point_t::Zero())
          , ddc1_(point_t::Zero()) {}

        ~ResultDataCOMTraj(){}

        bezier_t c_of_t_;
        bezier_t dL_of_t_;
        point_t dc1_;
        point_t ddc1_;
    };
} // end namespace bezier_com_traj

#endif
