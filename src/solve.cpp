/*
 * Copyright 2017, LAAS-CNRS
 * Author: Steve Tonneau
 */

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>

using namespace bezier_com_traj;

namespace bezier_com_traj
{
waypoint_t w0(point_t_tC p0, point_t_tC p1, point_t_tC g, const Matrix3& p0X, const Matrix3& /*p1X*/, const Matrix3& /*gX*/, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = 6*alpha* Matrix3::Identity();
    w.first.block<3,3>(3,0) = 6*alpha*p0X;
    w.second.head(3) = 6*alpha*(p0 - 2*p1);
    w.second.tail(3) =(-p0).cross(12*alpha*p1 + g);
    return w;
}

waypoint_t w1(point_t_tC p0, point_t_tC p1, point_t_tC /*g*/, const Matrix3& /*p0X*/, const Matrix3& /*p1X*/, const Matrix3& gX, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = 3*alpha*Matrix3::Identity();
    w.first.block<3,3>(3,0) = skew(1.5 * (3*p1 - p0))*alpha;
    w.second.head(3) = 1.5 *alpha* (3*p0 - 5*p1);
    w.second.tail(3) =(3*alpha*p0).cross(-p1) + 0.25 * (gX * (3*p1 + p0)) ;
    return w;
}

waypoint_t w2(point_t_tC p0, point_t_tC p1, point_t_tC g, const Matrix3& /*p0X*/, const Matrix3& /*p1X*/, const Matrix3& gX, const double alpha)
{
    waypoint_t w = initwp();
    // w.first.block<3,3>(0,0) = 0;
    w.first.block<3,3>(3,0) = skew(0.5*g - 3*alpha* p0 + 3*alpha*p1);
    w.second.head(3) = 3*alpha*(p0 - p1);
    w.second.tail(3) = 0.5 * gX*p1;
    return w;
}

waypoint_t w3 (point_t_tC p0, point_t_tC p1, point_t_tC g, const Matrix3& /*p0X*/, const Matrix3& /*p1X*/, const Matrix3& /*gX*/, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = -3*alpha*Matrix3::Identity();
    w.first.block<3,3>(3,0) = skew(g - 1.5 *alpha* (p1 + p0));
    w.second.head(3) = 1.5*alpha * (p1 + p0);
    //w.second.tail(3) = 0;
    return w;
}

waypoint_t w4 (point_t_tC /*p0*/, point_t_tC p1, point_t_tC g, const Matrix3& /*p0X*/, const Matrix3& /*p1X*/, const Matrix3& /*gX*/, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(0,0) = -6*alpha * Matrix3::Identity();
    w.first.block<3,3>(3,0) = skew(g - 6*alpha* p1);
    w.second.head(3) = 6*alpha*p1;
    //w.second.tail(3) = 0;
    return w;
}

waypoint_t u0 (point_t_tC l0, const double alpha)
{
    waypoint_t w = initwp();
    //w.first.block<3,3>(0,0) = 0;
    w.first.block<3,3>(3,0) = 3*alpha * Matrix3::Identity();
    //w.second.head(3) = 0;
    w.second.tail(3) = -3*alpha*l0;
    return w;
}

waypoint_t u1 (point_t_tC l0, const double alpha)
{
    waypoint_t w = initwp();
    //w.first.block<3,3>(0,0) = 0;
    //w.first.block<3,3>(3,0) = 0;
    //w.second.head(3) = 0;
    w.second.tail(3) = -1.5*alpha*l0;
    return w;
}

waypoint_t u2 (point_t_tC l0, const double alpha)
{
    waypoint_t w = initwp();
    //w.first.block<3,3>(0,0) = 0;
    w.first.block<3,3>(3,0) = -1.5*alpha * Matrix3::Identity();
    //w.second.head(3) = 0;
    w.second.tail(3) = -l0 / 2. * alpha;
    return w;
}

waypoint_t u3 (point_t_tC /*l0*/, const double alpha)
{
    waypoint_t w = initwp();
    w.first.block<3,3>(3,0) = -1.5*alpha * Matrix3::Identity();
    //w.second.head(3) = 0;
    //w.second.tail(3) = 0.;
    return w;
}


waypoint_t u4 (point_t_tC /*l0*/, const double /*alpha*/)
{
    waypoint_t w = initwp();
    //w.first.block<3,3>(0,0) = 0;
    //w.first.block<3,3>(3,0) = 0;
    //w.second.head(3) = 0;
    //w.second.tail(3) = 0.;
    return w;
}


std::vector<spline::Bern<double> > ComputeBersteinPolynoms()
{
    std::vector<spline::Bern<double> > res;
    for (unsigned int i =0; i <5; ++i)
        res.push_back(spline::Bern<double>(4,i));
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
    {
        std::vector<spline::Bern<double> > berns = ComputeBersteinPolynoms();
        wps = ComputeDiscretizedWaypoints(wps, berns, numSteps);
    }
    return wps;
}

std::vector<waypoint_t> ComputeAllWaypointsAngularMomentum(point_t_tC l0, const double T, const double timeStep)
{
    int numSteps = computeNumSteps(T, timeStep);
    double alpha = 1. / (T);
    std::vector<waypoint_t> wps;
    wps.push_back(u0(l0, alpha));
    wps.push_back(u1(l0, alpha));
    wps.push_back(u2(l0, alpha));
    wps.push_back(u3(l0, alpha));
    wps.push_back(u4(l0, alpha));
    if (numSteps > 0)
    {
        std::vector<spline::Bern<double> > berns = ComputeBersteinPolynoms();
        wps = ComputeDiscretizedWaypoints(wps, berns, numSteps);
    }
    return wps;
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
std::pair<MatrixXX, VectorX> compute6dControlPointInequalities(const ContactData& cData, point_t_tC c0, point_t_tC dc0, point_t_tC l0, const bool useAngMomentum, const double T, const double timeStep, bool& fail)
{
    std::vector<waypoint_t> wps, wpL;
    wps = ComputeAllWaypoints(c0, dc0, cData.contactPhase_->m_gravity, T, timeStep);
    if (useAngMomentum)
        wpL = ComputeAllWaypointsAngularMomentum(l0, T, timeStep);
    return compute6dControlPointInequalities(cData,wps,wpL, useAngMomentum, fail);
}

std::pair<MatrixXX, VectorX> computeCostFunction(point_t_tC p0, point_t_tC l0, const bool useAngMomentum)
{
    int dimPb = useAngMomentum ? 6 : 3;
    std::pair<MatrixXX, VectorX> res;
    res.first  = MatrixXX(dimPb,dimPb);
    res.second = VectorX (dimPb);
    Ref_matrixXX H = res.first;
    Ref_vectorX  g = res.second;

    //minimize distance to initial point
    double weightDist = useAngMomentum ? 0.1 : 1.;
    H.block<3,3>(0,0) = Matrix3::Identity() * weightDist;
    g.head(3) = - p0 * weightDist;

    // now angular momentum integral minimization
    if(useAngMomentum)
    {
        H.block<3,3>(3,3) = Matrix3::Identity() * 6./5.;
        g.tail(3) = 0.5 *(-(9.* l0) / 5.);
    }
    return res;
}

void computeRealCost(const ProblemData& pData, ResultData& resData)
{
    if(pData.useAngularMomentum_)
    {
        Vector3 xL = resData.x.tail(3);
        resData.cost_ = (1./5.)*(9.*pData.l0_.dot(pData.l0_) -  9.*pData.l0_.dot(xL) + 6.*xL.dot(xL));
    }
    else
        resData.cost_ = (pData.c0_ - resData.x).norm();
}

void computeC_of_T(const ProblemData& pData, const std::vector<double>& Ts, ResultDataCOMTraj& res)
{
    std::vector<Vector3> wps;
    wps.push_back(pData.c0_);
    wps.push_back(pData.dc0_ * Ts[0] / 3 + pData.c0_);
    wps.push_back(res.x.head(3));
    wps.push_back(res.x.head(3));
    res.c_of_t_ = bezier_t (wps.begin(), wps.end(),Ts[0]);
}

void computedL_of_T(const ProblemData& pData, const std::vector<double>& Ts, ResultDataCOMTraj& res)
{
    if(pData.useAngularMomentum_)
    {
        std::vector<Vector3> wps;
        wps.push_back(3*(res.x.tail(3)  - pData.l0_));
        wps.push_back(3*(-res.x.tail(3)));
        wps.push_back(Vector3::Zero());
        res.dL_of_t_ = bezier_t(wps.begin(), wps.end(),Ts[0], 1./Ts[0]);
    }
    else
        res.dL_of_t_ = bezier_t::zero(Ts[0]);
}

// no angular momentum for now
ResultDataCOMTraj solve0step(const ProblemData& pData,  const std::vector<double>& Ts, const double timeStep)
{
    assert(pData.contacts_.size() ==1);
    assert(Ts.size() == pData.contacts_.size());
    bool fail = true;
    std::pair<MatrixXX, VectorX> Ab = compute6dControlPointInequalities(pData.contacts_.front(),pData.c0_, pData.dc0_, pData.l0_, pData.useAngularMomentum_, Ts.front(),timeStep, fail);
    ResultDataCOMTraj res;
    if(fail)
        return res;
    std::pair<MatrixXX, VectorX> Hg = computeCostFunction(pData.c0_, pData.l0_, pData.useAngularMomentum_);
    int dimPb = pData.useAngularMomentum_ ? 6 : 3;
    VectorX init = VectorX(dimPb);
    init.head(3) = pData.c0_;
    if(dimPb > 3)
        init.tail(3) = pData.l0_;
    // rewriting 0.5 || Dx -d ||^2 as x'Hx  + g'x
    ResultData resQp = solve(Ab.first,Ab.second,Hg.first,Hg.second, init);
    if(resQp.success_)
    {
        res.success_ = true;
        res.x = resQp.x;
        computeRealCost(pData, res);
        computeC_of_T (pData,Ts,res);
        computedL_of_T(pData,Ts,res);
    }
    return res;
}
} // namespace bezier_com_traj
