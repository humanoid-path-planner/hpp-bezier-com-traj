/*
 * Copyright 2018, LAAS-CNRS
 * Author: Steve Tonneau
 */

#ifndef BEZIER_COM_TRAJ_LIB_DATA_H
#define BEZIER_COM_TRAJ_LIB_DATA_H

#include <bezier-com-traj/config.hh>
#include <bezier-com-traj/flags.hh>
#include <bezier-com-traj/definitions.hh>
#include <bezier-com-traj/utils.hh>

#include <spline/bezier_curve.h>
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>
#include <Eigen/Dense>

#include <vector>

namespace bezier_com_traj
{
    /**
    * @brief Contact data contains all the contact information
    * relative to a contact phase: contact points and normals
    * (within Equilibrium object), as well as any additional
    * kinematic and angular constraints.
    */
    struct BEZIER_COM_TRAJ_DLLAPI ContactData
    {
        ContactData()
            : contactPhase_(0)
            , Kin_(Eigen::Matrix3d::Zero())
            , kin_(VectorX::Zero(0))
            , Ang_(Eigen::Matrix3d::Zero())
            , ang_(VectorX::Zero(0)) {}

        ContactData(const ContactData& other)
            : contactPhase_(new centroidal_dynamics::Equilibrium(*(other.contactPhase_)))
            , Kin_(other.Kin_)
            , kin_(other.kin_)
            , Ang_(other.Ang_)
            , ang_(other.ang_){}

        ContactData(centroidal_dynamics::Equilibrium* contactPhase)
            : contactPhase_(contactPhase)
            , Kin_(Eigen::Matrix3d::Zero())
            , kin_(VectorX::Zero(0))
            , Ang_(Eigen::Matrix3d::Zero())
            , ang_(VectorX::Zero(0)) {}
       ~ContactData(){}

       centroidal_dynamics::Equilibrium* contactPhase_;
       MatrixX3 Kin_; // inequality kinematic constraints
       VectorX  kin_;
       MatrixX3 Ang_; // inequality angular momentum constraints
       VectorX  ang_;
    };

    /**
    * @brief Used to define the constraints on the trajectory generation problem.
    * Flags are used to constrain initial and terminal com positions an derivatives.
    * Additionally, the maximum acceleration can be bounded.
    */
    struct BEZIER_COM_TRAJ_DLLAPI Constraints
    {
        Constraints()
            : flag_(INIT_POS | INIT_VEL | END_VEL | END_POS)
            , constraintAcceleration_(false)
            , maxAcceleration_(10.)
            , reduce_h_(1e-3) {}

        Constraints(ConstraintFlag flag)
            : flag_(flag)
            , constraintAcceleration_(false)
            , maxAcceleration_(10.)
            , reduce_h_(1e-3) {}

        ~Constraints(){}

        ConstraintFlag flag_;
        bool constraintAcceleration_;
        double maxAcceleration_;
        double reduce_h_;
    };


    /**
    * @brief Defines all the inputs of the problem:
    * Initial and terminal constraints, as well as selected
    * cost functions. Also,a list of ContactData defines the
    * different phases of the problem. While the method
    * can handle any phase greater than one, using more
    * than three phases is probably too constraining.
    */
    struct BEZIER_COM_TRAJ_DLLAPI ProblemData
    {
        ProblemData()
            : c0_ (point_t::Zero())
            ,dc0_ (point_t::Zero())
            ,ddc0_(point_t::Zero())
            ,j0_(point_t::Zero())
            , c1_ (point_t::Zero())
            ,dc1_ (point_t::Zero())
            ,ddc1_(point_t::Zero())
            ,j1_(point_t::Zero())
            ,useAngularMomentum_(false)
            ,costFunction_(ACCELERATION) {}

        std::vector<ContactData> contacts_;
        point_t  c0_,dc0_,ddc0_,j0_,c1_,dc1_,ddc1_,j1_;
        point_t  l0_;
        bool useAngularMomentum_;
        Constraints constraints_;
        CostFunction costFunction_;
    };


    /**
    * @brief Struct used to return the results of the trajectory generation
    * problem.
    */
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
        bool success_; // whether the optimization was successful
        double cost_; // cost evaluation for the solved control point
        VectorX x; //control point
    };


    /**
    * @brief Specialized ResultData that computes the Bezier curves
    * corresponding to the computed trajectory
    */
    struct BEZIER_COM_TRAJ_DLLAPI ResultDataCOMTraj : public ResultData
    {
        ResultDataCOMTraj():
            ResultData()
          , c_of_t_(bezier_t::zero())
          , dL_of_t_(bezier_t::zero())
          , dc1_(point_t::Zero())
          , ddc1_(point_t::Zero()) {}
        ~ResultDataCOMTraj(){}

        bezier_t c_of_t_; // center of mass trajectory
        bezier_t dL_of_t_; // angular momentum derivative trajectory
        point_t dc1_; // terminal velocity
        point_t ddc1_; //terminal acceleration
    };
} // end namespace bezier_com_traj

#endif
