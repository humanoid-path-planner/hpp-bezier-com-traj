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

#ifndef GLPKWRAPPER_HH_
#define GLPKWRAPPER_HH_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#define DEFAULT_MAX_ITER 1000


namespace solvers
{
    typedef Eigen::MatrixXd MatrixXd;
    typedef Eigen::VectorXd VectorXd;
    typedef Eigen::VectorXi VectorXi;
    enum glpk_status
    {
      glpk_OPTIMAL=0,
      glpk_INFEASIBLE=1,
      glpk_UNBOUNDED=2,
      glpk_MAX_ITER_REACHED=3,
    };

    // min g'x
    // st  CIx <= ci0
    //     CEx  = ce0
    glpk_status solve(const VectorXd & g0,
                      const MatrixXd & CE,
                      const VectorXd & ce0,
                      const MatrixXd & CI,
                      const VectorXd & ci0,
                      VectorXd& x);

} /* namespace solvers */

#endif /* GLPKWRAPPER_HH_ */
