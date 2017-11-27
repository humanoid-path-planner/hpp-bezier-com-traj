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
            , kin_(Vector3::Zero())
            , Ang_(Eigen::Matrix3d::Zero())
            , ang_(Vector3::Zero()) {}
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

        bool success_;
        double cost_;
        VectorX x;
    };

    struct BEZIER_COM_TRAJ_DLLAPI ResultDataCOMTraj : public ResultData
    {
        ResultDataCOMTraj():
            ResultData()
          , c_of_t_(0)
          , dL_of_t_(0) {}

        ResultDataCOMTraj(const ResultDataCOMTraj& other):
            ResultData(other.success_,other.cost_, other.x)
        {
            if(other.constC_of_t())
                c_of_t_   = new bezier_t(*(other.constC_of_t()));
            if(other.constDL_of_t())
                dL_of_t_  = new bezier_t(*(other.constDL_of_t()));
        }

        ~ResultDataCOMTraj()
        {
            if(c_of_t_)
                delete c_of_t_;
            if(dL_of_t_)
                delete dL_of_t_;
        }
        /**
         * @brief C_of_t return a copy of the trajectory curve.
         * Only valid if c_of_t_ was initialized(success_ is true.)
         * @return a copy of C_of_t
         */
        bezier_t C_of_t() const
        {
            assert(c_of_t_);
            return *c_of_t_;
        }
        /**
         * @brief DL_of_t return a copy of the trajectory curve.
         * Only valid if c_of_t_ was initialized(success_ is true.)
         * @return a copy of DL_of_t
         */
        bezier_t DL_of_t() const
        {
            assert(dL_of_t_);
            return *dL_of_t_;
        }

        const bezier_t* constC_of_t() const
        {
            return c_of_t_;
        }

        const bezier_t* constDL_of_t() const
        {
            return dL_of_t_;
        }

        void SetC_of_t(bezier_t* c_of_t)
        {
            c_of_t_ = c_of_t;
        }

        void SetDL_of_t(bezier_t* dL_of_t)
        {
            dL_of_t_ = dL_of_t;
        }


    private:
        bezier_t* c_of_t_;
        bezier_t* dL_of_t_;
    };

} // end namespace bezier_com_traj

#endif
