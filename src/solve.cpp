/*
 * Copyright 2017, LAAS-CNRS
 * Author: Steve Tonneau
 */

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>

using namespace bezier_com_traj;

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
    return  (wx, ws)*/


waypoint_t u0 (point_t_tC l0, const double alpha)
{
    waypoint_t w = initwp();
    //w.first.block<3,3>(0,0) = 0;
    w.first.block<3,3>(3,0) = 3*alpha * Matrix3::Identity();
    //w.second.head(3) = 0;
    w.second.tail(3) = -3*alpha*l0;
    return w;
}

/*
#angular momentum waypoints
def u0(l0, alpha):
    ux, us = __init_6D()
    ux[3:] = identity(3)* 3 * alpha
    us[3:] = -3*alpha*l0[:]
    return  (ux, us)
*/


waypoint_t u1 (point_t_tC l0, const double alpha)
{
    waypoint_t w = initwp();
    //w.first.block<3,3>(0,0) = 0;
    //w.first.block<3,3>(3,0) = 0;
    //w.second.head(3) = 0;
    w.second.tail(3) = -1.5*alpha*l0;
    return w;
}


/*
def u1(l0, alpha):
    ux, us = __init_6D()
    us[3:] = -1.5*l0*alpha
    return  (ux, us)
*/

waypoint_t u2 (point_t_tC l0, const double alpha)
{
    waypoint_t w = initwp();
    //w.first.block<3,3>(0,0) = 0;
    w.first.block<3,3>(3,0) = -1.5*alpha * Matrix3::Identity();
    //w.second.head(3) = 0;
    w.second.tail(3) = -l0 / 2. * alpha;
    return w;
}


/*
def u2(l0, alpha):
    ux, us = __init_6D()
    ux[3:] = identity(3)* (-1.5) * alpha
    us[3:] = -l0 / 2. * alpha
*/

waypoint_t u3 (point_t_tC /*l0*/, const double alpha)
{
    waypoint_t w = initwp();
    //w.first.block<3,3>(0,0) = 0;
    w.first.block<3,3>(3,0) = -1.5*alpha * Matrix3::Identity();
    //w.second.head(3) = 0;
    //w.second.tail(3) = 0.;
    return w;
}


/*   return  (ux, us)

def u3(l0, alpha):
    ux, us = __init_6D()
    ux[3:] = identity(3)*  (-1.5) * alpha
    return  (ux, us)
*/

waypoint_t u4 (point_t_tC /*l0*/, const double /*alpha*/)
{
    waypoint_t w = initwp();
    //w.first.block<3,3>(0,0) = 0;
    //w.first.block<3,3>(3,0) = 0;
    //w.second.head(3) = 0;
    //w.second.tail(3) = 0.;
    return w;
}


/*
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
std::pair<MatrixXX, VectorX> compute6dControlPointInequalities(const ContactData& cData, point_t_tC c0, point_t_tC dc0, point_t_tC l0, const bool useAngMomentum, double T, double timeStep)
{
    std::vector<waypoint_t> wps, wpL;
    wps = ComputeAllWaypoints(c0, dc0, cData.contactPhase_->m_gravity, T, timeStep);
    if (useAngMomentum)
        wpL = ComputeAllWaypointsAngularMomentum(l0, T, timeStep);
    return compute6dControlPointInequalities(cData,wps,wpL, useAngMomentum,T, timeStep);
}

std::pair<MatrixXX, VectorX> computeCostFunction(point_t_tC p0, point_t_tC l0, const bool useAngMomentum)
{
    int dimPb = useAngMomentum ? 6 : 3;
    std::pair<MatrixXX, VectorX> res;
    res.first  = MatrixXX(dimPb,dimPb);
    res.second = VectorX (dimPb);
    Ref_matrixXX D = res.first;
    Ref_vectorX  d = res.second;

    //minimize distance to initial point
    double weightDist = useAngMomentum ? 0. : 1.;
    D.block<3,3>(0,0) = Matrix3::Identity() * weightDist;
    d.head(3) = p0 * weightDist;

    // now angular momentum integral minimization
    if(useAngMomentum)
    {
        double alpha = sqrt(12./5.);
        D.block<3,3>(3,3) = Matrix3::Identity() * alpha;
        d.tail(3) = (9.* l0) / (5. * alpha);
    }
    return res;
}
/*weight_dist_or = 1. if l0 == None else 0.
            #weight_dist_or = 0
            D = identity(dim_pb);
            alpha = sqrt(12./5.)
            for i in range(3):
                D[i,i] = weight_dist_or
            d = zeros(dim_pb);
            d[:3]= self._p0 * weight_dist_or
            if(l0 != None):
                # minimizing integral of angular momentum
                for i in range(3,6):
                    D[i,i] = alpha
                d[3:]= (9.* l0) / (5. * alpha)
            D = (D[:]); d = (d[:]); A = (self.__Ain[:]);
            lbA = (-100000.* ones(self.__Ain.shape[0]))[:]; ubA=(self.__Aub);
            lb = (-100. * ones(dim_pb))[:]; ub = (100. * ones(dim_pb))[:];
            self._qp_solver.setProblemData(D = D , d = d, A=A, lbA=lbA, ubA=ubA, lb = lb, ub = ub, x0=None)
            (x, imode) =  self._qp_solver.solve(D = D , d = d, A=A, lbA=lbA, ubA=ubA, lb = lb, ub = ub, x0=None)
            if l0 == None:
                cost = norm(self._p0 - x)
            else:
                cost = (1./5.)*(9.*l0.dot(l0) -  9.*l0.dot(x[3:]) + 6.*x[3:].dot(x[3:]))*/

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

// no angular momentum for now
ResultData solve0step(const ProblemData& pData,  const std::vector<double> Ts, const double timeStep)
{
    assert(pData.contacts_.size() ==1);
    assert(Ts.size() == pData.contacts_.size());    
    std::pair<MatrixXX, VectorX> Ab = compute6dControlPointInequalities(pData.contacts_.front(),pData.c0_, pData.dc0_, pData.l0_, pData.useAngularMomentum_, Ts.front(),timeStep);
    std::pair<MatrixXX, VectorX> Dd = computeCostFunction(pData.c0_, pData.l0_, pData.useAngularMomentum_);
    ResultData res = solve(Ab.first,Ab.second,Dd.first,Dd.second);
    if(res.success_)
        computeRealCost(pData, res);
    return res;
}

