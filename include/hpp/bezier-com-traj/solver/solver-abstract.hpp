//
// Copyright (c) 2017 CNRS
//
// This file is part of tsid
// tsid is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
// tsid is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// tsid If not, see
// <http://www.gnu.org/licenses/>.
//

#ifndef SOLVERABSTRACT_HH_
#define SOLVERABSTRACT_HH_

#include <hpp/bezier-com-traj/local_config.hh>
#include <Eigen/Dense>

namespace solvers
{

/**
* Possible states of the solver.
*/
enum optim_status
{
  OPTIM_OPTIMAL=0,
  OPTIM_INFEASIBLE=1
};

static const double UNBOUNDED_UP   =  100000.;
static const double UNBOUNDED_DOWN = -100000.;

typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::VectorXi VectorXi;
typedef const Eigen::Ref<const VectorXd>     & Cref_vectorX;

enum BEZIER_COM_TRAJ_DLLAPI SolverType
{
    SOLVER_QUADPROG         = 0x00001
    //SOLVER_QUADPROG_SPARSE  = 0x00002
#ifdef USE_GLPK_SOLVER
    ,SOLVER_GLPK            = 0x00002
#endif
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
      , x(VectorXd::Zero(0)){}

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
    VectorXd x; //control point
};

// min g'x
// st  CIx <= ci0
//     CEx  = ce0
/**
* @brief solve Solve a QP or LP given
* init position and velocity, 0 velocity constraints (acceleration constraints are ignored)
* @param pData problem Data. Should contain only one contact phase.
* @param Ts timelength of each contact phase. Should only contain one value
* @param timeStep time that the solver has to stop.
* @return ResultData a struct containing the resulting trajectory, if success is true.
*/
ResultData BEZIER_COM_TRAJ_DLLAPI solve(const MatrixXd & A,
                  const VectorXd & b,
                  const MatrixXd & D,
                  const VectorXd & d,
                  const MatrixXd & Hess,
                  const VectorXd & g,
                  const VectorXd & initGuess,
                  Cref_vectorX minBounds, Cref_vectorX maxBounds,
                  const SolverType solver);


} /* namespace solvers */

#endif /* SOLVERABSTRACT_HH_ */
