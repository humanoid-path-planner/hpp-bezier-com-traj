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

#include <hpp/bezier-com-traj/solver/solver-abstract.hpp>

#include <Eigen/Dense>

namespace solvers {

// min g'x
// st  CIx <= ci0
//     CEx  = ce0
int solveglpk(const VectorXd& g0, const MatrixXd& CE, const VectorXd& ce0, const MatrixXd& CI, const VectorXd& ci0,
              solvers::Cref_vectorX minBounds, solvers::Cref_vectorX maxBounds, VectorXd& x, double& cost);

} /* namespace solvers */

#endif /* GLPKWRAPPER_HH_ */
