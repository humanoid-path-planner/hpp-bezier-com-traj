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


#include "solver/solver-abstract.hpp"
#ifdef USE_GLPK_SOLVER
#include <solver/glpk-wrapper.hpp>
#endif
#include <solver/eiquadprog-fast.hpp>

#include <Eigen/Sparse>

namespace solvers
{

template<typename Derived>
inline bool is_nan(const Eigen::MatrixBase<Derived>& x)
{
    bool isnan = !((x.array()==x.array()).all());
    return isnan;
}


typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;
typedef Eigen::SparseVector<int>    SpVeci;

ResultData solve(const MatrixXd & A,
                  const VectorXd & b,
                  const MatrixXd & D,
                  const VectorXd & d,
                  const MatrixXd & Hess,
                  const VectorXd & g,
                  const VectorXd & initGuess,
                  const SOLVER_TYPE solver)
{
    assert (!(is_nan(A)));
    assert (!(is_nan(b)));
    assert (!(is_nan(D)));
    assert (!(is_nan(d)));
    assert (!(is_nan(initGuess)));
    assert (!(is_nan(Hess)));
    ResultData res;
    res.x = initGuess;
    switch(solver)
    {
        /*
       * solves the problem
       * min. x' Hess x + 2 g0' x
       * s.t. CE x + ce0 = 0
       *      CI x + ci0 >= 0
       * Thus CI = -A; ci0 = b
       * CI = D; ce0 = -d
       */
        case SOLVER_QUADPROG:
        case SOLVER_QUADPROG_SPARSE:
        {
            MatrixXd CI = -A;
            //MatrixXd CE = D;
            VectorXd ce0  = -d;
            tsid::solvers::EiquadprogFast QPsolver = tsid::solvers::EiquadprogFast();
            tsid::solvers::EiquadprogFast_status status;
            if(solver == SOLVER_QUADPROG)
                status = QPsolver.solve_quadprog(Hess,g,D,ce0,CI,b,res.x);
            else
                status = QPsolver.solve_quadprog_sparse(Hess.sparseView(),g,D,ce0,CI,b,res.x);
            res.success_ = (status == tsid::solvers::EIQUADPROG_FAST_OPTIMAL );
            if(res.success_)
                res.cost_ = QPsolver.getObjValue();
            return res;
        }
#ifdef USE_GLPK_SOLVER
        case SOLVER_GLPK:
    {
        res.success_ = (solvers::solve(g,D,d,A,b,res.x,res.cost_) == 0);
        return res;
    }
#endif
    default:
    throw std::runtime_error("Unknown solver type in solver-asbtract");
    }
}

} /* namespace solvers */
