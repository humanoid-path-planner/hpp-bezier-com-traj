/*
 * Copyright 2017, LAAS-CNRS
 * Author: Steve Tonneau
 */

#include <bezier-com-traj/solve.hh>
#include <spline/bernstein.h>
//#include <qpOASES.hpp>


using namespace centroidal_dynamics;

typedef Eigen::Matrix<double, 6, 3, 0, 6, 3> matrix6_t;
typedef Eigen::Matrix<double, 6, 1, 0, 6, 1> point6_t;
typedef std::pair<matrix6_t, point6_t> waypoint_t;
typedef Eigen::Matrix3d matrix3_t;

typedef const Eigen::Ref<const matrix3_t>& matrix3_tC;
typedef const Eigen::Ref<const point_t>& point_t_tC;

typedef spline::bezier_curve  <double, double, 6, true, point6_t> bezier6_t;

matrix3_t skew(point_t_tC x)
{
    matrix3_t res = matrix3_t::Zero();
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

waypoint_t w0(point_t_tC p0, point_t_tC p1, point_t_tC g, matrix3_tC p0X, matrix3_tC /*p1X*/, matrix3_tC /*gX*/, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = 6*alpha*matrix3_t::Identity();
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

waypoint_t w1(point_t_tC p0, point_t_tC p1, point_t_tC /*g*/, matrix3_tC /*p0X*/, matrix3_tC /*p1X*/, matrix3_tC gX, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = 3*alpha*matrix3_t::Identity();
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
waypoint_t w2(point_t_tC p0, point_t_tC p1, point_t_tC g, matrix3_tC /*p0X*/, matrix3_tC /*p1X*/, matrix3_tC gX, const double alpha)
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
waypoint_t w3 (point_t_tC p0, point_t_tC p1, point_t_tC g, matrix3_tC /*p0X*/, matrix3_tC /*p1X*/, matrix3_tC /*gX*/, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = -3*alpha*matrix3_t::Identity();
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
waypoint_t w4 (point_t_tC /*p0*/, point_t_tC p1, point_t_tC g, matrix3_tC /*p0X*/, matrix3_tC /*p1X*/, matrix3_tC /*gX*/, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = -6*alpha *matrix3_t::Identity();
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

std::vector<waypoint_t> ComputeAllWaypoints(point_t_tC p0, point_t_tC p1, point_t_tC g, double T, int numSteps = -1 )
{
    matrix3_t p0X = skew(p0);
    matrix3_t p1X = skew(p1);
    matrix3_t  gX = skew( g);
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

// no angular momentum for now
ResultData solve0step(const ProblemData& pData, const double T, int numSteps = -1)
{
    assert(pData.contacts_.size() ==1);
    ResultData res;
    return res;
}

