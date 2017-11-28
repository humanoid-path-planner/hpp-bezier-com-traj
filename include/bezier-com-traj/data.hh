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

namespace bezier_com_traj
{
    typedef double value_type;
    typedef Eigen::Matrix <value_type, 3, 3>                           Matrix3;
    typedef Eigen::Matrix <value_type, 6, 3>                           Matrix63;
    typedef Eigen::Matrix <value_type, Eigen::Dynamic, 3>              MatrixX3;
    typedef Eigen::Matrix <value_type, Eigen::Dynamic, Eigen::Dynamic> MatrixXX;
    typedef centroidal_dynamics::Vector3 Vector3;
    typedef centroidal_dynamics::Vector6 Vector6;
    typedef centroidal_dynamics::VectorX VectorX;

    typedef Eigen::Ref<Vector3>     Ref_vector3;
    typedef Eigen::Ref<VectorX>     Ref_vectorX;
    typedef Eigen::Ref<MatrixX3>    Ref_matrixX3;
    typedef Eigen::Ref<MatrixXX>    Ref_matrixXX;

    typedef const Eigen::Ref<const Vector3>     & Cref_vector3;
    typedef const Eigen::Ref<const Vector6>     & Cref_vector6;
    typedef const Eigen::Ref<const VectorX>     & Cref_vectorX;
    typedef const Eigen::Ref<const MatrixXX>    & Cref_matrixXX;
    typedef const Eigen::Ref<const MatrixX3>    & Cref_matrixX3;

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

    struct BEZIER_COM_TRAJ_DLLAPI ProblemData
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
          , dL_of_t_(bezier_t::zero()) {}

        ~ResultDataCOMTraj(){}

        bezier_t c_of_t_;
        bezier_t dL_of_t_;
    };

} // end namespace bezier_com_traj

#endif
