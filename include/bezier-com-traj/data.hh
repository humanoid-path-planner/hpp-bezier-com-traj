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


    typedef Matrix63 matrix6_t;
    typedef Vector6 point6_t;
    typedef Matrix3 matrix3_t;
    typedef Vector3 point3_t;
    /**
    * @brief waypoint_t a waypoint is composed of a  6*3 matrix that depend
    * on the variable x, and of a 6d vector independent of x, such that
    * each control point of the target bezier curve is given by pi = wix * x + wis
    */
    typedef std::pair<matrix6_t, point6_t> waypoint6_t;
    typedef std::pair<matrix3_t, point3_t> waypoint3_t;
    typedef std::pair<double,point3_t> coefs_t;


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

    enum BEZIER_COM_TRAJ_DLLAPI ConstraintFlag{
        INIT_POS = 0x00001,
        INIT_VEL = 0x00002,
        INIT_ACC = 0x00004,
        END_POS  = 0x00008,
        END_VEL  = 0x00010,
        END_ACC  = 0x00020
      };

    inline ConstraintFlag operator~(ConstraintFlag a)
    {return static_cast<ConstraintFlag>(~static_cast<const int>(a));}

    inline ConstraintFlag operator|(ConstraintFlag a, ConstraintFlag b)
    {return static_cast<ConstraintFlag>(static_cast<const int>(a) | static_cast<const int>(b));}

    inline ConstraintFlag operator&(ConstraintFlag a, ConstraintFlag b)
    {return static_cast<ConstraintFlag>(static_cast<const int>(a) & static_cast<const int>(b));}

    inline ConstraintFlag operator^(ConstraintFlag a, ConstraintFlag b)
    {return static_cast<ConstraintFlag>(static_cast<const int>(a) ^ static_cast<const int>(b));}

    inline ConstraintFlag& operator|=(ConstraintFlag& a, ConstraintFlag b)
    {return (ConstraintFlag&)((int&)(a) |= static_cast<const int>(b));}

    inline ConstraintFlag& operator&=(ConstraintFlag& a, ConstraintFlag b)
    {return (ConstraintFlag&)((int&)(a) &= static_cast<const int>(b));}

    inline ConstraintFlag& operator^=(ConstraintFlag& a, ConstraintFlag b)
    {return (ConstraintFlag&)((int&)(a) ^= static_cast<const int>(b));}

    struct BEZIER_COM_TRAJ_DLLAPI Constraints
    {

        Constraints()
            : flag_(INIT_POS | INIT_VEL | END_VEL | END_POS)
            , constraintAcceleration_(true)
            , maxAcceleration_(5.)
            ,reduce_h_(1e-4) {}

        Constraints(ConstraintFlag flag)
            : flag_(flag)
            , constraintAcceleration_(true)
            , maxAcceleration_(5.)
            ,reduce_h_(1e-4)
        {
            /*if(dc0_)
                assert(c0_ && "You cannot constraint init velocity if init position is not constrained.");
            if(ddc0_)
                assert(dc0_ && "You cannot constraint init acceleration if init velocity is not constrained.");
            if(dc1_)
                assert(c1_ && "You cannot constraint final velocity if final position is not constrained.");
            if(ddc1_)
                assert(dc1_ && "You cannot constraint final acceleration if final velocity is not constrained.");*/
        }

        ~Constraints(){}

        ConstraintFlag flag_;
        bool constraintAcceleration_;
        double maxAcceleration_;
        double reduce_h_;

    };


    struct BEZIER_COM_TRAJ_DLLAPI ProblemData
    {
        ProblemData()
            : c0_(Vector3::Zero())
            ,dc0_(Vector3::Zero())
            ,ddc0_(Vector3::Zero())
            , c1_(Vector3::Zero())
            ,dc1_(Vector3::Zero())
            ,ddc1_(Vector3::Zero())
            ,useAngularMomentum_(false) {}

        std::vector<ContactData> contacts_;
        Vector3  c0_,dc0_,ddc0_,c1_,dc1_,ddc1_;
        Vector3  l0_;
        bool useAngularMomentum_;
        Constraints constraints_;
    };

    typedef Eigen::Vector3d point_t;
    typedef const Eigen::Ref<const point_t>& point_t_tC;
    typedef spline::bezier_curve  <double, double, 3, true, point_t > bezier_t;
    typedef spline::bezier_curve  <double, double, 6, true, point6_t> bezier6_t;
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
