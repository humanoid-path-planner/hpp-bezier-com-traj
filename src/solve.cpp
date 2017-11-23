/*
 * Copyright 2017, LAAS-CNRS
 * Author: Steve Tonneau
 */

#include <bezier-com-traj/solve.hh>
#include <spline/bernstein.h>
//#include <qpOASES.hpp>


using namespace centroidal_dynamics;

typedef Matrix63 matrix6_t;
typedef Vector6 point6_t;
typedef std::pair<matrix6_t, point6_t> waypoint_t;

typedef const Eigen::Ref<const point_t>& point_t_tC;

typedef spline::bezier_curve  <double, double, 6, true, point6_t> bezier6_t;

Matrix3 skew(point_t_tC x)
{
    Matrix3 res = Matrix3::Zero();
    res(0,1) = - x(2); res(0,2) =   x(1);
    res(1,0) =   x(2); res(1,2) = - x(0);
    res(2,0) = - x(1); res(2,1) =   x(0);
    return res;
}

waypoint_t initwp()
{
    waypoint_t w;
    w.first  = matrix6_t::Zero();
    w.second = point6_t::Zero();
    return w;
}

waypoint_t w0(point_t_tC p0, point_t_tC p1, point_t_tC g, const Matrix3& p0X, const Matrix3& /*p1X*/, const Matrix3& /*gX*/, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = 6*alpha* Matrix3::Identity();
    w.first.block<3,3>(3,0) = 6*alpha*p0X;
    w.second.head(3) = 6*alpha*(p0 - 2*p1);
    w.second.tail(3) =-p0.cross(12*alpha*p1 + g);
    return w;
}

/*def w0(p0, p1, g, p0X, p1X, gX, alpha):
    wx, ws = __init_6D()
    wx[:3,:] = 6*alpha*identity(3);  wx[3:,:] = 6*alpha*p0X;
    ws[:3]   = 6*alpha*(p0 - 2*p1)
    ws[3:]   = X(-p0, 12*alpha*p1 + g )
    return  (wx, ws)*/

waypoint_t w1(point_t_tC p0, point_t_tC p1, point_t_tC /*g*/, const Matrix3& /*p0X*/, const Matrix3& /*p1X*/, const Matrix3& gX, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = 3*alpha*Matrix3::Identity();
    w.first.block<3,3>(3,0) = skew(1.5 * (3*p1 - p0))*alpha;
    w.second.head(3) = 1.5 *alpha* (3*p0 - 5*p1);
    w.second.tail(3) =(3*alpha*p0).cross(-p1) + 0.25 * (gX * (3*p1 + p0)) ;
    return w;
}

/*def w1(p0, p1, g, p0X, p1X, gX, alpha):
    wx, ws = __init_6D()
    wx[:3,:] = 3*alpha*identity(3);
    wx[3:,:] = skew(1.5 * (3*p1 - p0))*alpha
    ws[:3]   =  1.5 *alpha* (3*p0 - 5*p1);
    ws[3:]   = X(3*alpha*p0, -p1) + 0.25 * (gX.dot(3*p1 + p0))
    return  (wx, ws)
*/
waypoint_t w2(point_t_tC p0, point_t_tC p1, point_t_tC g, const Matrix3& /*p0X*/, const Matrix3& /*p1X*/, const Matrix3& gX, const double alpha)
{
    waypoint_t w = initwp();
    // w.first.block<3,3>(0,0) = 0;
    w.first.block<3,3>(3,0) = skew(0.5*g - 3*alpha* p0 + 3*alpha*p1);
    w.second.head(3) = 3*alpha*(p0 - p1);
    w.second.tail(3) = 0.5 * gX*p1;
    return w;
}

/*
def w2(p0, p1, g, p0X, p1X, gX, alpha):
    wx, ws = __init_6D()
    #~ wx[:3,:] = 0;
    wx[3:,:] = skew(0.5*g - 3*alpha* p0 + 3*alpha*p1)
    ws[:3]   =  3*alpha*(p0 - p1);
    ws[3:]   = 0.5 * gX.dot(p1)
    return  (wx, ws)
*/
waypoint_t w3 (point_t_tC p0, point_t_tC p1, point_t_tC g, const Matrix3& /*p0X*/, const Matrix3& /*p1X*/, const Matrix3& /*gX*/, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = -3*alpha*Matrix3::Identity();
    w.first.block<3,3>(3,0) = skew(g - 1.5 *alpha* (p1 + p0));
    w.second.head(3) = 1.5*alpha * (p1 + p0);
    //w.second.tail(3) = 0;
    return w;
}

/*
def w3(p0, p1, g, p0X, p1X, gX, alpha):
    wx, ws = __init_6D()
    wx[:3,:] = -3*alpha* identity(3);
    wx[3:,:] = skew(g - 1.5 *alpha* (p1 + p0))
    ws[:3]   = 1.5*alpha * (p1 + p0)
    #~ ws[3:]   = 0
    return  (wx, ws)
*/
waypoint_t w4 (point_t_tC /*p0*/, point_t_tC p1, point_t_tC g, const Matrix3& /*p0X*/, const Matrix3& /*p1X*/, const Matrix3& /*gX*/, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = -6*alpha * Matrix3::Identity();
    w.first.block<3,3>(3,0) = skew(g - 6*alpha* p1);
    w.second.head(3) = 6*alpha*p1;
    //w.second.tail(3) = 0;
    return w;
}

/*
def w4(p0, p1, g, p0X, p1X, gX, alpha):
    wx, ws = __init_6D()
    wx[:3,:] = -6*alpha *identity(3);
    wx[3:,:] = skew(g - 6*alpha* p1)
    ws[:3]   = 6*alpha*p1
    #~ ws[3:]   = 0
    return  (wx, ws)

#angular momentum waypoints
def u0(l0, alpha):
    ux, us = __init_6D()
    ux[3:] = identity(3)* 3 * alpha
    us[3:] = -3*alpha*l0[:]
    return  (ux, us)

def u1(l0, alpha):
    ux, us = __init_6D()
    us[3:] = -1.5*l0*alpha
    return  (ux, us)

def u2(l0, alpha):
    ux, us = __init_6D()
    ux[3:] = identity(3)* (-1.5) * alpha
    us[3:] = -l0 / 2. * alpha
    return  (ux, us)

def u3(l0, alpha):
    ux, us = __init_6D()
    ux[3:] = identity(3)*  (-1.5) * alpha
    return  (ux, us)

def u4(l0, alpha):
    ux, us = __init_6D()
    return  (ux, us)*/

std::vector<spline::Bern<double> > ComputeBersteinPolynoms()
{
    std::vector<spline::Bern<double> > res;
    for (unsigned int i =0; i <5; ++i)
        res.push_back(spline::Bern<double>(4,i));
    return res;
}

std::vector<waypoint_t> ComputeDiscretizedWaypoints(const std::vector<waypoint_t>& wps,  int numSteps)
{
    double dt = 1./numSteps;
    std::vector<spline::Bern<double> > berns = ComputeBersteinPolynoms();
    std::vector<waypoint_t> res;
    for (int i =0; i < numSteps + 1; ++i)
    {
        waypoint_t w = initwp();
        for (int j = 0; j <5; ++j)
        {
            double b = berns[j](i*dt);
            w.first +=b*(wps[j].first );
            w.second+=b*(wps[j].second);
        }
        res.push_back(w);
    }
    return res;
}


int computeNumSteps(const double T, const double timeStep)
{
    return timeStep > 0. ? int(T / timeStep) : -1;
}

std::vector<waypoint_t> ComputeAllWaypoints(point_t_tC p0, point_t_tC dc0, point_t_tC g, const double T, const double timeStep)
{
    int numSteps = computeNumSteps(T, timeStep);
    static const double  n = 3.; //degree
    point_t p1 = dc0 * T / n +  p0;
    Matrix3 p0X = skew(p0);
    Matrix3 p1X = skew(p1);
    Matrix3  gX = skew( g);
    double alpha = 1. / (T*T);
    std::vector<waypoint_t> wps;
    wps.push_back(w0(p0, p1, g, p0X, p1X, gX, alpha));
    wps.push_back(w1(p0, p1, g, p0X, p1X, gX, alpha));
    wps.push_back(w2(p0, p1, g, p0X, p1X, gX, alpha));
    wps.push_back(w3(p0, p1, g, p0X, p1X, gX, alpha));
    wps.push_back(w4(p0, p1, g, p0X, p1X, gX, alpha));
    if (numSteps > 0)
        wps = ComputeDiscretizedWaypoints(wps, numSteps);
    return wps;
}

MatrixXX initMatrixA(const int dimH, const std::vector<waypoint_t>& wps, Cref_vectorX kin)
{
    int dimKin = kin == point_t::Zero() ? 0 : (int)(kin.rows());
    return MatrixXX::Zero(dimH * wps.size() + dimKin, 3);
}


/*def __add_kinematic_and_normalize(self,A,b):
    if self._kinematic_constraints != None:
        dim_kin = self._kinematic_constraints[0].shape[0]
        A[-dim_kin:,:] = self._kinematic_constraints[0][:]
        b[-dim_kin:] =  self._kinematic_constraints[1][:]
    A, b = normalize(A,b)
    self.__Ain = A[:]; self.__Aub = b[:]*/

void addKinematicAndNormalize(Cref_matrixXX A, Cref_vectorX b, Cref_matrixX3 Kin, Cref_vectorX kin)
{

}

/* compute the inequality methods that determine the 6D bezier curve w(t)
as a function of a variable waypoint for the 3D COM trajectory.
The initial curve is of degree 3 (init pos and velocity, 0 velocity constraints + one free variable).
The 6d curve is of degree 2*n-2 = 4, thus 5 control points are to be computed.
Each control point produces a 6 * 3 inequality matrix wix, and a 6 *1 column right member wsi.
Premultiplying it by H gives mH w_xi * x <= mH_wsi where m is the mass
Stacking all of these results in a big inequality matrix A and a column vector x that determines the constraints
On the 6d curves, Ain x <= Aub
*/
std::pair<MatrixXX, VectorX> compute6dControlPointInequalities(const ContactData& cData, point_t_tC c0, point_t_tC dc0, point_t_tC l0, const bool useAngMomentum, double T, double timeStep)
{
    std::pair<MatrixXX, VectorX> res;
    MatrixXX& A = res.first;
    VectorX&  b = res.second;
    // gravity vector
    point_t g = point_t::Zero(); g(2) = -9.81;
    // compute waypoints
    std::vector<waypoint_t> wps = ComputeAllWaypoints(c0, dc0, g, T, timeStep);
    // compute GIWC
    MatrixXX H; VectorX h;
    cData.contactPhase_->getPolytopeInequalities(H,h);
    H = -H;
    int dimH = (int)(H.rows());
    MatrixXX mH = cData.contactPhase_->m_mass * H;
    // init and fill Ab matrix
    A = initMatrixA(dimH, wps, cData.kin_);
    b = VectorX::Zero(A.rows());
    point6_t bc = point6_t::Zero(); bc.head(3) = g; // constant part of Aub, Aubi = mH * (bc - wsi)
    int i = 0;
    for (std::vector<waypoint_t>::const_iterator wpcit = wps.begin(); wpcit != wps.end(); ++wpcit)
    {
        A.block(i*dimH,0, dimH, 6) = mH * wpcit->first;
        b.segment(i*dimH, dimH) = mH * (bc - wpcit->second);
        ++i;
    }
    addKinematicAndNormalize(A,b, cData.Kin_,cData.kin_);
    if (useAngMomentum)
    {

    }
    return res;

    /*
use_angular_momentum = l0 != None
        A,b = self.__add_kinematic_and_normalize(A,b, not use_angular_momentum)
        if use_angular_momentum:
            A,b = self.__add_angular_momentum(A,b, l0, T, num_steps)
        self.__Ain = A[:]; self.__Aub = b[:]
*/
}

// no angular momentum for now
ResultData solve0step(const ProblemData& pData,  const std::vector<double> Ts, const double timeStep)
{
    assert(pData.contacts_.size() ==1);
    assert(Ts.size() == pData.contacts_.size());
    compute6dControlPointInequalities(pData.contacts_.front(),pData.c0_, pData.dc0_, pData.l0_, false, Ts.front(),timeStep);

    ResultData res;
    return res;
}

