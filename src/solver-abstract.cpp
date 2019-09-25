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

#include "hpp/bezier-com-traj/solver/solver-abstract.hpp"
#ifdef USE_GLPK_SOLVER
#include <hpp/bezier-com-traj/solver/glpk-wrapper.hpp>
#include <glpk.h>
#endif
#include <hpp/bezier-com-traj/solver/eiquadprog-fast.hpp>

#include <Eigen/Sparse>
#include <stdexcept>

namespace solvers {

template <typename Derived>
inline bool is_nan(const Eigen::MatrixBase<Derived>& x) {
  bool isnan = !((x.array() == x.array()).all());
  return isnan;
}

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;
typedef Eigen::SparseVector<int> SpVeci;

namespace {
void addConstraintMinBoundQuadProg(solvers::Cref_vectorX minBounds, std::pair<MatrixXd, VectorXd>& data) {
  if (minBounds.size() == 0) return;
  MatrixXd& res = data.first;
  VectorXd& resv = data.second;
  MatrixXd D(res.rows() + res.cols(), res.cols());
  VectorXd d(resv.rows() + res.cols());
  D.block(0, 0, res.rows(), res.cols()) = res;
  D.block(res.rows(), 0, res.cols(), res.cols()) = (-1.) * MatrixXd::Identity(res.cols(), res.cols());
  d.head(resv.size()) = resv;
  d.tail(res.cols()) = -minBounds;
  data.first = D;
  data.second = d;
}

void addConstraintMaxBoundQuadProg(solvers::Cref_vectorX maxBounds, std::pair<MatrixXd, VectorXd>& data) {
  if (maxBounds.size() == 0) return;
  MatrixXd& res = data.first;
  VectorXd& resv = data.second;
  MatrixXd D(res.rows() + res.cols() - 3, res.cols());
  VectorXd d(resv.rows() + res.cols());
  D.block(0, 0, res.rows(), res.cols()) = res;
  D.block(res.rows(), 0, res.cols(), res.cols()) = MatrixXd::Identity(res.cols(), res.cols());
  d.head(resv.size()) = resv;
  d.tail(res.cols()) = maxBounds;
  data.first = D;
  data.second = d;
}

std::pair<MatrixXd, VectorXd> addBoundaryConstraintsQuadProg(solvers::Cref_vectorX minBounds,
                                                             solvers::Cref_vectorX maxBounds, const MatrixXd& CI,
                                                             const VectorXd& ci0) {
  std::pair<MatrixXd, VectorXd> data;
  data.first = CI;
  data.second = ci0;
  addConstraintMinBoundQuadProg(minBounds, data);
  addConstraintMaxBoundQuadProg(maxBounds, data);
  return data;
}
}  // namespace

ResultData solve(const MatrixXd& A, const VectorXd& b, const MatrixXd& D, const VectorXd& d, const MatrixXd& Hess,
                 const VectorXd& g, const VectorXd& initGuess, solvers::Cref_vectorX minBounds,
                 solvers::Cref_vectorX maxBounds, const SolverType solver) {
  assert(!(is_nan(A)));
  assert(!(is_nan(b)));
  assert(!(is_nan(D)));
  assert(!(is_nan(d)));
  assert(!(is_nan(initGuess)));
  ResultData res;
  res.x = initGuess;
  switch (solver) {
    /*
     * solves the problem
     * min. x' Hess x + 2 g0' x
     * s.t. CE x + ce0 = 0
     *      CI x + ci0 >= 0
     * Thus CI = -A; ci0 = b
     * CI = D; ce0 = -d
     */
    case SOLVER_QUADPROG:
      // case SOLVER_QUADPROG_SPARSE:
      {
        assert(!(is_nan(Hess)));
        std::pair<MatrixXd, VectorXd> CIp = addBoundaryConstraintsQuadProg(minBounds, maxBounds, A, b);
        VectorXd ce0 = -d;
        tsid::solvers::EiquadprogFast QPsolver = tsid::solvers::EiquadprogFast();
        tsid::solvers::EiquadprogFast_status status;
        // if(solver == SOLVER_QUADPROG)
        status = QPsolver.solve_quadprog(Hess, g, D, ce0, -CIp.first, CIp.second, res.x);
        /* else
         {
             SpMat Hsp = Hess.sparseView();
             status = QPsolver.solve_quadprog_sparse(Hsp,g,D,ce0,CI,b,res.x);
         }*/
        res.success_ = (status == tsid::solvers::EIQUADPROG_FAST_OPTIMAL);
        if (res.success_) res.cost_ = QPsolver.getObjValue();
        return res;
      }
#ifdef USE_GLPK_SOLVER
    case SOLVER_GLPK: {
      res.success_ = (solvers::solveglpk(g, D, d, A, b, minBounds, maxBounds, res.x, res.cost_) == GLP_OPT);
      return res;
    }
#endif
    default:
      throw std::runtime_error("Unknown solver type in solver-asbtract");
  }
}

} /* namespace solvers */
