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

#include "hpp/bezier-com-traj/solver/eiquadprog-fast.hpp"

#include <iostream>
namespace tsid {
namespace solvers {

/// Compute sqrt(a^2 + b^2)
template <typename Scalar>
inline Scalar distance(Scalar a, Scalar b) {
  Scalar a1, b1, t;
  a1 = std::abs(a);
  b1 = std::abs(b);
  if (a1 > b1) {
    t = (b1 / a1);
    return a1 * std::sqrt(1.0 + t * t);
  } else if (b1 > a1) {
    t = (a1 / b1);
    return b1 * std::sqrt(1.0 + t * t);
  }
  return a1 * std::sqrt(2.0);
}

EiquadprogFast::EiquadprogFast() {
  m_maxIter = DEFAULT_MAX_ITER;
  q = 0;  // size of the active set A (containing the indices of the active
          // constraints)
  is_inverse_provided_ = false;
  m_nVars = 0;
  m_nEqCon = 0;
  m_nIneqCon = 0;
}

EiquadprogFast::~EiquadprogFast() {}

void EiquadprogFast::reset(int nVars, int nEqCon, int nIneqCon) {
  m_nVars = nVars;
  m_nEqCon = nEqCon;
  m_nIneqCon = nIneqCon;
  m_J.setZero(nVars, nVars);
  chol_.compute(m_J);
  R.resize(nVars, nVars);
  s.resize(nIneqCon);
  r.resize(nIneqCon + nEqCon);
  u.resize(nIneqCon + nEqCon);
  z.resize(nVars);
  d.resize(nVars);
  np.resize(nVars);
  A.resize(nIneqCon + nEqCon);
  iai.resize(nIneqCon);
  iaexcl.resize(nIneqCon);
  x_old.resize(nVars);
  u_old.resize(nIneqCon + nEqCon);
  A_old.resize(nIneqCon + nEqCon);

#ifdef OPTIMIZE_ADD_CONSTRAINT
  T1.resize(nVars);
#endif
}

bool EiquadprogFast::add_constraint(MatrixXd& R, MatrixXd& J, VectorXd& d,
                                    int& iq, double& R_norm) {
  long int nVars = J.rows();
#ifdef TRACE_SOLVER
  std::cout << "Add constraint " << iq << '/';
#endif
  long int j, k;
  double cc, ss, h, t1, t2, xny;

#ifdef OPTIMIZE_ADD_CONSTRAINT
  Eigen::Vector2d cc_ss;
#endif

  /* we have to find the Givens rotation which will reduce the element
           d(j) to zero.
           if it is already zero we don't have to do anything, except of
           decreasing j */
  for (j = d.size() - 1; j >= iq + 1; j--) {
    /* The Givens rotation is done with the matrix (cc cs, cs -cc).
                    If cc is one, then element (j) of d is zero compared with
       element (j - 1). Hence we don't have to do anything. If cc is zero, then
       we just have to switch column (j) and column (j - 1) of J. Since we only
       switch columns in J, we have to be careful how we update d depending on
       the sign of gs. Otherwise we have to apply the Givens rotation to these
       columns. The i - 1 element of d has to be updated to h. */
    cc = d(j - 1);
    ss = d(j);
    h = distance(cc, ss);
    if (h == 0.0) continue;
    d(j) = 0.0;
    ss = ss / h;
    cc = cc / h;
    if (cc < 0.0) {
      cc = -cc;
      ss = -ss;
      d(j - 1) = -h;
    } else
      d(j - 1) = h;
    xny = ss / (1.0 + cc);

// #define OPTIMIZE_ADD_CONSTRAINT
#ifdef OPTIMIZE_ADD_CONSTRAINT  // the optimized code is actually slower than
                                // the original
    T1 = J.col(j - 1);
    cc_ss(0) = cc;
    cc_ss(1) = ss;
    J.col(j - 1).noalias() = J.middleCols<2>(j - 1) * cc_ss;
    J.col(j) = xny * (T1 + J.col(j - 1)) - J.col(j);
#else
    // J.col(j-1) = J[:,j-1:j] * [cc; ss]
    for (k = 0; k < nVars; k++) {
      t1 = J(k, j - 1);
      t2 = J(k, j);
      J(k, j - 1) = t1 * cc + t2 * ss;
      J(k, j) = xny * (t1 + J(k, j - 1)) - t2;
    }
#endif
  }
  /* update the number of constraints added*/
  iq++;
  /* To update R we have to put the iq components of the d vector
into column iq - 1 of R
*/
  R.col(iq - 1).head(iq) = d.head(iq);
#ifdef TRACE_SOLVER
  std::cout << iq << std::endl;
#endif

  if (std::abs(d(iq - 1)) <= std::numeric_limits<double>::epsilon() * R_norm)
    // problem degenerate
    return false;
  R_norm = std::max<double>(R_norm, std::abs(d(iq - 1)));
  return true;
}

void EiquadprogFast::delete_constraint(MatrixXd& R, MatrixXd& J, VectorXi& A,
                                       VectorXd& u, int nEqCon, int& iq,
                                       int l) {
  const long int nVars = R.rows();
#ifdef TRACE_SOLVER
  std::cout << "Delete constraint " << l << ' ' << iq;
#endif
  int i, j, k;
  int qq = 0;
  double cc, ss, h, xny, t1, t2;

  /* Find the index qq for active constraint l to be removed */
  for (i = nEqCon; i < iq; i++)
    if (A(i) == l) {
      qq = i;
      break;
    }

  /* remove the constraint from the active set and the duals */
  for (i = qq; i < iq - 1; i++) {
    A(i) = A(i + 1);
    u(i) = u(i + 1);
    R.col(i) = R.col(i + 1);
  }

  A(iq - 1) = A(iq);
  u(iq - 1) = u(iq);
  A(iq) = 0;
  u(iq) = 0.0;
  for (j = 0; j < iq; j++) R(j, iq - 1) = 0.0;
  /* constraint has been fully removed */
  iq--;
#ifdef TRACE_SOLVER
  std::cout << '/' << iq << std::endl;
#endif

  if (iq == 0) return;

  for (j = qq; j < iq; j++) {
    cc = R(j, j);
    ss = R(j + 1, j);
    h = distance(cc, ss);
    if (h == 0.0) continue;
    cc = cc / h;
    ss = ss / h;
    R(j + 1, j) = 0.0;
    if (cc < 0.0) {
      R(j, j) = -h;
      cc = -cc;
      ss = -ss;
    } else
      R(j, j) = h;

    xny = ss / (1.0 + cc);
    for (k = j + 1; k < iq; k++) {
      t1 = R(j, k);
      t2 = R(j + 1, k);
      R(j, k) = t1 * cc + t2 * ss;
      R(j + 1, k) = xny * (t1 + R(j, k)) - t2;
    }
    for (k = 0; k < nVars; k++) {
      t1 = J(k, j);
      t2 = J(k, j + 1);
      J(k, j) = t1 * cc + t2 * ss;
      J(k, j + 1) = xny * (J(k, j) + t1) - t2;
    }
  }
}

template <class Derived>
void print_vector(const char* name, Eigen::MatrixBase<Derived>& x, int n) {
  if (x.size() < 10) std::cout << name << x.transpose() << std::endl;
}
template <class Derived>
void print_matrix(const char* name, Eigen::MatrixBase<Derived>& x, int n) {
  //  std::cout << name << std::endl << x << std::endl;
}

EiquadprogFast_status EiquadprogFast::solve_quadprog(
    const MatrixXd& Hess, const VectorXd& g0, const MatrixXd& CE,
    const VectorXd& ce0, const MatrixXd& CI, const VectorXd& ci0, VectorXd& x) {
  const int nVars = (int)g0.size();
  const int nEqCon = (int)ce0.size();
  const int nIneqCon = (int)ci0.size();

  if (nVars != m_nVars || nEqCon != m_nEqCon || nIneqCon != m_nIneqCon)
    reset(nVars, nEqCon, nIneqCon);

  assert(Hess.rows() == m_nVars && Hess.cols() == m_nVars);
  assert(g0.size() == m_nVars);
  assert(CE.rows() == m_nEqCon && CE.cols() == m_nVars);
  assert(ce0.size() == m_nEqCon);
  assert(CI.rows() == m_nIneqCon && CI.cols() == m_nVars);
  assert(ci0.size() == m_nIneqCon);

  int i, k, l;    // indices
  int ip;         // index of the chosen violated constraint
  int iq;         // current number of active constraints
  double psi;     // current sum of constraint violations
  double c1;      // Hessian trace
  double c2;      // Hessian Chowlesky factor trace
  double ss;      // largest constraint violation (negative for violation)
  double R_norm;  // norm of matrix R
  const double inf = std::numeric_limits<double>::infinity();
  double t, t1, t2;
  /* t is the step length, which is the minimum of the partial step length t1
   * and the full step length t2 */

  iter = 0;  // active-set iteration number

  /*
   * Preprocessing phase
   */
  /* compute the trace of the original matrix Hess */
  c1 = Hess.trace();

  /* decompose the matrix Hess in the form LL^T */
  if (!is_inverse_provided_) {
    START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_CHOWLESKY_DECOMPOSITION);
    chol_.compute(Hess);
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_CHOWLESKY_DECOMPOSITION);
  }

  /* initialize the matrix R */
  d.setZero(nVars);
  R.setZero(nVars, nVars);
  R_norm = 1.0;

  /* compute the inverse of the factorized matrix Hess^-1, this is the initial
   * value for H */
  // m_J = L^-T
  if (!is_inverse_provided_) {
    START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_CHOWLESKY_INVERSE);
    m_J.setIdentity(nVars, nVars);
#ifdef OPTIMIZE_HESSIAN_INVERSE
    chol_.matrixU().solveInPlace(m_J);
#else
    m_J = chol_.matrixU().solve(m_J);
#endif
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_CHOWLESKY_INVERSE);
  }

  c2 = m_J.trace();
#ifdef _SOLVER
  print_matrix("m_J", m_J, nVars);
#endif

  /* c1 * c2 is an estimate for cond(Hess) */

  /*
   * Find the unconstrained minimizer of the quadratic form 0.5 * x Hess x + g0
   * x this is a feasible point in the dual space x = Hess^-1 * g0
   */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1_UNCONSTR_MINIM);
  if (is_inverse_provided_) {
    x = m_J * (m_J.transpose() * g0);
  } else {
#ifdef OPTIMIZE_UNCONSTR_MINIM
    x = -g0;
    chol_.solveInPlace(x);
  }
#else
    x = chol_.solve(g0);
  }
  x = -x;
#endif
  /* and compute the current solution value */
  f_value = 0.5 * g0.dot(x);
#ifdef TRACE_SOLVER
  std::cout << "Unconstrained solution: " << f_value << std::endl;
  print_vector("x", x, nVars);
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1_UNCONSTR_MINIM);

  /* Add equality constraints to the working set A */

  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR);
  iq = 0;
  for (i = 0; i < nEqCon; i++) {
    START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR_1);
    np = CE.row(i);
    compute_d(d, m_J, np);
    update_z(z, m_J, d, iq);
    update_r(R, r, d, iq);

#ifdef TRACE_SOLVER
    print_matrix("R", R, iq);
    print_vector("z", z, nVars);
    print_vector("r", r, iq);
    print_vector("d", d, nVars);
#endif

    /* compute full step length t2: i.e., the minimum step in primal space s.t.
    the contraint becomes feasible */
    t2 = 0.0;
    if (std::abs(z.dot(z)) >
        std::numeric_limits<double>::epsilon())  // i.e. z != 0
      t2 = (-np.dot(x) - ce0(i)) / z.dot(np);

    x += t2 * z;

    /* set u = u+ */
    u(iq) = t2;
    u.head(iq) -= t2 * r.head(iq);

    /* compute the new solution value */
    f_value += 0.5 * (t2 * t2) * z.dot(np);
    A(i) = -i - 1;
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR_1);

    START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR_2);
    if (!add_constraint(R, m_J, d, iq, R_norm)) {
      // Equality constraints are linearly dependent
      return EIQUADPROG_FAST_REDUNDANT_EQUALITIES;
    }
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR_2);
  }
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR);

  /* set iai = K \ A */
  for (i = 0; i < nIneqCon; i++) iai(i) = i;

#ifdef USE_WARM_START
  //      DEBUG_STREAM("Gonna warm start using previous active
  //      set:\n"<<A.transpose()<<"\n")
  for (i = nEqCon; i < q; i++) {
    iai(i - nEqCon) = -1;
    ip = A(i);
    np = CI.row(ip);
    compute_d(d, m_J, np);
    update_z(z, m_J, d, iq);
    update_r(R, r, d, iq);

    /* compute full step length t2: i.e., the minimum step in primal space s.t.
    the contraint becomes feasible */
    t2 = 0.0;
    if (std::abs(z.dot(z)) >
        std::numeric_limits<double>::epsilon())  // i.e. z != 0
      t2 = (-np.dot(x) - ci0(ip)) / z.dot(np);
    else
      DEBUG_STREAM("[WARM START] z=0\n")

    x += t2 * z;

    /* set u = u+ */
    u(iq) = t2;
    u.head(iq) -= t2 * r.head(iq);

    /* compute the new solution value */
    f_value += 0.5 * (t2 * t2) * z.dot(np);

    if (!add_constraint(R, m_J, d, iq, R_norm)) {
      // constraints are linearly dependent
      std::cout << "[WARM START] Constraints are linearly dependent\n";
      return RT_EIQUADPROG_REDUNDANT_EQUALITIES;
    }
  }
#else

#endif

l1:
  iter++;
  if (iter >= m_maxIter) {
    q = iq;
    return EIQUADPROG_FAST_MAX_ITER_REACHED;
  }

  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1);

#ifdef TRACE_SOLVER
  print_vector("x", x, nVars);
#endif
  /* step 1: choose a violated constraint */
  for (i = nEqCon; i < iq; i++) {
    ip = A(i);
    iai(ip) = -1;
  }

  /* compute s(x) = ci^T * x + ci0 for all elements of K \ A */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1_2);
  ss = 0.0;
  ip = 0; /* ip will be the index of the chosen violated constraint */

#ifdef OPTIMIZE_STEP_1_2
  s = ci0;
  s.noalias() += CI * x;
  iaexcl.setOnes();
  psi = (s.cwiseMin(VectorXd::Zero(nIneqCon))).sum();
#else
  psi = 0.0; /* this value will contain the sum of all infeasibilities */
  for (i = 0; i < nIneqCon; i++) {
    iaexcl(i) = 1;
    s(i) = CI.row(i).dot(x) + ci0(i);
    psi += std::min(0.0, s(i));
  }
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1_2);
#ifdef TRACE_SOLVER
  print_vector("s", s, nIneqCon);
#endif

  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1);

  if (std::abs(psi) <=
      nIneqCon * std::numeric_limits<double>::epsilon() * c1 * c2 * 100.0) {
    /* numerically there are not infeasibilities anymore */
    q = iq;
    //        DEBUG_STREAM("Optimal active
    //        set:\n"<<A.head(iq).transpose()<<"\n\n")
    return EIQUADPROG_FAST_OPTIMAL;
  }

  /* save old values for u, x and A */
  u_old.head(iq) = u.head(iq);
  A_old.head(iq) = A.head(iq);
  x_old = x;

l2: /* Step 2: check for feasibility and determine a new S-pair */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2);
  // find constraint with highest violation (what about normalizing
  // constraints?)
  for (i = 0; i < nIneqCon; i++) {
    if (s(i) < ss && iai(i) != -1 && iaexcl(i)) {
      ss = s(i);
      ip = i;
    }
  }
  if (ss >= 0.0) {
    q = iq;
    //        DEBUG_STREAM("Optimal active set:\n"<<A.transpose()<<"\n\n")
    return EIQUADPROG_FAST_OPTIMAL;
  }

  /* set np = n(ip) */
  np = CI.row(ip);
  /* set u = (u 0)^T */
  u(iq) = 0.0;
  /* add ip to the active set A */
  A(iq) = ip;

  //      DEBUG_STREAM("Add constraint "<<ip<<" to active set\n")

#ifdef TRACE_SOLVER
  std::cout << "Trying with constraint " << ip << std::endl;
  print_vector("np", np, nVars);
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2);

l2a: /* Step 2a: determine step direction */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2A);
  /* compute z = H np: the step direction in the primal space (through m_J, see
   * the paper) */
  compute_d(d, m_J, np);
  //    update_z(z, m_J, d, iq);
  if (iq >= nVars) {
    //      throw std::runtime_error("iq >= m_J.cols()");
    z.setZero();
  } else {
    update_z(z, m_J, d, iq);
  }
  /* compute N* np (if q > 0): the negative of the step direction in the dual
   * space */
  update_r(R, r, d, iq);
#ifdef TRACE_SOLVER
  std::cout << "Step direction z" << std::endl;
  print_vector("z", z, nVars);
  print_vector("r", r, iq + 1);
  print_vector("u", u, iq + 1);
  print_vector("d", d, nVars);
  print_vector("A", A, iq + 1);
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2A);

  /* Step 2b: compute step length */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2B);
  l = 0;
  /* Compute t1: partial step length (maximum step in dual space without
   * violating dual feasibility */
  t1 = inf; /* +inf */
  /* find the index l s.t. it reaches the minimum of u+(x) / r */
  // l: index of constraint to drop (maybe)
  for (k = nEqCon; k < iq; k++) {
    double tmp;
    if (r(k) > 0.0 && ((tmp = u(k) / r(k)) < t1)) {
      t1 = tmp;
      l = A(k);
    }
  }
  /* Compute t2: full step length (minimum step in primal space such that the
   * constraint ip becomes feasible */
  if (std::abs(z.dot(z)) >
      std::numeric_limits<double>::epsilon())  // i.e. z != 0
    t2 = -s(ip) / z.dot(np);
  else
    t2 = inf; /* +inf */

  /* the step is chosen as the minimum of t1 and t2 */
  t = std::min(t1, t2);
#ifdef TRACE_SOLVER
  std::cout << "Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2
            << ") ";
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2B);

  /* Step 2c: determine new S-pair and take step: */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);
  /* case (i): no step in primal or dual space */
  if (t >= inf) {
    /* QPP is infeasible */
    q = iq;
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);
    return EIQUADPROG_FAST_UNBOUNDED;
  }
  /* case (ii): step in dual space */
  if (t2 >= inf) {
    /* set u = u +  t * [-r 1) and drop constraint l from the active set A */
    u.head(iq) -= t * r.head(iq);
    u(iq) += t;
    iai(l) = l;
    delete_constraint(R, m_J, A, u, nEqCon, iq, l);
#ifdef TRACE_SOLVER
    std::cout << " in dual space: " << f_value << std::endl;
    print_vector("x", x, nVars);
    print_vector("z", z, nVars);
    print_vector("A", A, iq + 1);
#endif
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);
    goto l2a;
  }

  /* case (iii): step in primal and dual space */
  x += t * z;
  /* update the solution value */
  f_value += t * z.dot(np) * (0.5 * t + u(iq));

  u.head(iq) -= t * r.head(iq);
  u(iq) += t;

#ifdef TRACE_SOLVER
  std::cout << " in both spaces: " << f_value << std::endl;
  print_vector("x", x, nVars);
  print_vector("u", u, iq + 1);
  print_vector("r", r, iq + 1);
  print_vector("A", A, iq + 1);
#endif

  if (t == t2) {
#ifdef TRACE_SOLVER
    std::cout << "Full step has taken " << t << std::endl;
    print_vector("x", x, nVars);
#endif
    /* full step has taken */
    /* add constraint ip to the active set*/
    if (!add_constraint(R, m_J, d, iq, R_norm)) {
      iaexcl(ip) = 0;
      delete_constraint(R, m_J, A, u, nEqCon, iq, ip);
#ifdef TRACE_SOLVER
      print_matrix("R", R, nVars);
      print_vector("A", A, iq);
#endif
      for (i = 0; i < nIneqCon; i++) iai(i) = i;
      for (i = 0; i < iq; i++) {
        A(i) = A_old(i);
        iai(A(i)) = -1;
        u(i) = u_old(i);
      }
      x = x_old;
      STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);
      goto l2; /* go to step 2 */
    } else
      iai(ip) = -1;
#ifdef TRACE_SOLVER
    print_matrix("R", R, nVars);
    print_vector("A", A, iq);
#endif
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);
    goto l1;
  }

  /* a partial step has been taken => drop constraint l */
  iai(l) = l;
  delete_constraint(R, m_J, A, u, nEqCon, iq, l);
  s(ip) = CI.row(ip).dot(x) + ci0(ip);

#ifdef TRACE_SOLVER
  std::cout << "Partial step has taken " << t << std::endl;
  print_vector("x", x, nVars);
  print_matrix("R", R, nVars);
  print_vector("A", A, iq);
  print_vector("s", s, nIneqCon);
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);

  goto l2a;
}

double trace_sparse(const EiquadprogFast::SpMat& Hess) {
  double sum = 0;
  for (int k = 0; k < Hess.outerSize(); ++k) sum += Hess.coeff(k, k);
  return sum;
}

EiquadprogFast_status EiquadprogFast::solve_quadprog_sparse(
    const EiquadprogFast::SpMat& Hess, const VectorXd& g0, const MatrixXd& CE,
    const VectorXd& ce0, const MatrixXd& CI, const VectorXd& ci0, VectorXd& x) {
  const int nVars = (int)g0.size();
  const int nEqCon = (int)ce0.size();
  const int nIneqCon = (int)ci0.size();

  if (nVars != m_nVars || nEqCon != m_nEqCon || nIneqCon != m_nIneqCon)
    reset(nVars, nEqCon, nIneqCon);

  assert(Hess.rows() == m_nVars && Hess.cols() == m_nVars);
  assert(g0.size() == m_nVars);
  assert(CE.rows() == m_nEqCon && CE.cols() == m_nVars);
  assert(ce0.size() == m_nEqCon);
  assert(CI.rows() == m_nIneqCon && CI.cols() == m_nVars);
  assert(ci0.size() == m_nIneqCon);

  int i, k, l;    // indices
  int ip;         // index of the chosen violated constraint
  int iq;         // current number of active constraints
  double psi;     // current sum of constraint violations
  double c1;      // Hessian trace
  double c2;      // Hessian Chowlesky factor trace
  double ss;      // largest constraint violation (negative for violation)
  double R_norm;  // norm of matrix R
  const double inf = std::numeric_limits<double>::infinity();
  double t, t1, t2;
  /* t is the step length, which is the minimum of the partial step length t1
   * and the full step length t2 */

  iter = 0;  // active-set iteration number

  /*
   * Preprocessing phase
   */
  /* compute the trace of the original matrix Hess */
  c1 = trace_sparse(Hess);

  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_CHOWLESKY_DECOMPOSITION);
  Eigen::SimplicialLLT<SpMat, Eigen::Lower> cholsp_(Hess);
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_CHOWLESKY_DECOMPOSITION);

  /* initialize the matrix R */
  d.setZero(nVars);
  R.setZero(nVars, nVars);
  R_norm = 1.0;

  /* compute the inverse of the factorized matrix Hess^-1, this is the initial
   * value for H */
  // m_J = L^-T
  if (!is_inverse_provided_) {
    START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_CHOWLESKY_INVERSE);
    m_J.setIdentity(nVars, nVars);
#ifdef OPTIMIZE_HESSIAN_INVERSE
    cholsp_.matrixU().solve(m_J);
#else
    m_J = cholsp_.matrixU().solve(m_J);
#endif
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_CHOWLESKY_INVERSE);
  }

  c2 = m_J.trace();
#ifdef _SOLVER
  print_matrix("m_J", m_J, nVars);
#endif

  /* c1 * c2 is an estimate for cond(Hess) */

  /*
   * Find the unconstrained minimizer of the quadratic form 0.5 * x Hess x + g0
   * x this is a feasible point in the dual space x = Hess^-1 * g0
   */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1_UNCONSTR_MINIM);
  if (is_inverse_provided_) {
    x = m_J * (m_J.transpose() * g0);
  } else {
#ifdef OPTIMIZE_UNCONSTR_MINIM
    // x = -g0;
    x = cholsp_.solve(-g0);
  }
#else
    x = cholsp_.solve(g0);
  }
  x = -x;
#endif
  /* and compute the current solution value */
  f_value = 0.5 * g0.dot(x);
#ifdef TRACE_SOLVER
  std::cout << "Unconstrained solution: " << f_value << std::endl;
  print_vector("x", x, nVars);
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1_UNCONSTR_MINIM);

  /* Add equality constraints to the working set A */

  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR);
  iq = 0;
  for (i = 0; i < nEqCon; i++) {
    START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR_1);
    np = CE.row(i);
    compute_d(d, m_J, np);
    update_z(z, m_J, d, iq);
    update_r(R, r, d, iq);

#ifdef TRACE_SOLVER
    print_matrix("R", R, iq);
    print_vector("z", z, nVars);
    print_vector("r", r, iq);
    print_vector("d", d, nVars);
#endif

    /* compute full step length t2: i.e., the minimum step in primal space s.t.
    the contraint becomes feasible */
    t2 = 0.0;
    if (std::abs(z.dot(z)) >
        std::numeric_limits<double>::epsilon())  // i.e. z != 0
      t2 = (-np.dot(x) - ce0(i)) / z.dot(np);

    x += t2 * z;

    /* set u = u+ */
    u(iq) = t2;
    u.head(iq) -= t2 * r.head(iq);

    /* compute the new solution value */
    f_value += 0.5 * (t2 * t2) * z.dot(np);
    A(i) = -i - 1;
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR_1);

    START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR_2);
    if (!add_constraint(R, m_J, d, iq, R_norm)) {
      // Equality constraints are linearly dependent
      return EIQUADPROG_FAST_REDUNDANT_EQUALITIES;
    }
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR_2);
  }
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_ADD_EQ_CONSTR);

  /* set iai = K \ A */
  for (i = 0; i < nIneqCon; i++) iai(i) = i;

#ifdef USE_WARM_START
  //      DEBUG_STREAM("Gonna warm start using previous active
  //      set:\n"<<A.transpose()<<"\n")
  for (i = nEqCon; i < q; i++) {
    iai(i - nEqCon) = -1;
    ip = A(i);
    np = CI.row(ip);
    compute_d(d, m_J, np);
    update_z(z, m_J, d, iq);
    update_r(R, r, d, iq);

    /* compute full step length t2: i.e., the minimum step in primal space s.t.
    the contraint becomes feasible */
    t2 = 0.0;
    if (std::abs(z.dot(z)) >
        std::numeric_limits<double>::epsilon())  // i.e. z != 0
      t2 = (-np.dot(x) - ci0(ip)) / z.dot(np);
    else
      DEBUG_STREAM("[WARM START] z=0\n")

    x += t2 * z;

    /* set u = u+ */
    u(iq) = t2;
    u.head(iq) -= t2 * r.head(iq);

    /* compute the new solution value */
    f_value += 0.5 * (t2 * t2) * z.dot(np);

    if (!add_constraint(R, m_J, d, iq, R_norm)) {
      // constraints are linearly dependent
      std::cout << "[WARM START] Constraints are linearly dependent\n";
      return RT_EIQUADPROG_REDUNDANT_EQUALITIES;
    }
  }
#else

#endif

l1:
  iter++;
  if (iter >= m_maxIter) {
    q = iq;
    return EIQUADPROG_FAST_MAX_ITER_REACHED;
  }

  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1);

#ifdef TRACE_SOLVER
  print_vector("x", x, nVars);
#endif
  /* step 1: choose a violated constraint */
  for (i = nEqCon; i < iq; i++) {
    ip = A(i);
    iai(ip) = -1;
  }

  /* compute s(x) = ci^T * x + ci0 for all elements of K \ A */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1_2);
  ss = 0.0;
  ip = 0; /* ip will be the index of the chosen violated constraint */

#ifdef OPTIMIZE_STEP_1_2
  s = ci0;
  s.noalias() += CI * x;
  iaexcl.setOnes();
  psi = (s.cwiseMin(VectorXd::Zero(nIneqCon))).sum();
#else
  psi = 0.0; /* this value will contain the sum of all infeasibilities */
  for (i = 0; i < nIneqCon; i++) {
    iaexcl(i) = 1;
    s(i) = CI.row(i).dot(x) + ci0(i);
    psi += std::min(0.0, s(i));
  }
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1_2);
#ifdef TRACE_SOLVER
  print_vector("s", s, nIneqCon);
#endif

  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_1);

  if (std::abs(psi) <=
      nIneqCon * std::numeric_limits<double>::epsilon() * c1 * c2 * 100.0) {
    /* numerically there are not infeasibilities anymore */
    q = iq;
    //        DEBUG_STREAM("Optimal active
    //        set:\n"<<A.head(iq).transpose()<<"\n\n")
    return EIQUADPROG_FAST_OPTIMAL;
  }

  /* save old values for u, x and A */
  u_old.head(iq) = u.head(iq);
  A_old.head(iq) = A.head(iq);
  x_old = x;

l2: /* Step 2: check for feasibility and determine a new S-pair */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2);
  // find constraint with highest violation (what about normalizing
  // constraints?)
  for (i = 0; i < nIneqCon; i++) {
    if (s(i) < ss && iai(i) != -1 && iaexcl(i)) {
      ss = s(i);
      ip = i;
    }
  }
  if (ss >= 0.0) {
    q = iq;
    //        DEBUG_STREAM("Optimal active set:\n"<<A.transpose()<<"\n\n")
    return EIQUADPROG_FAST_OPTIMAL;
  }

  /* set np = n(ip) */
  np = CI.row(ip);
  /* set u = (u 0)^T */
  u(iq) = 0.0;
  /* add ip to the active set A */
  A(iq) = ip;

  //      DEBUG_STREAM("Add constraint "<<ip<<" to active set\n")

#ifdef TRACE_SOLVER
  std::cout << "Trying with constraint " << ip << std::endl;
  print_vector("np", np, nVars);
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2);

l2a: /* Step 2a: determine step direction */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2A);
  /* compute z = H np: the step direction in the primal space (through m_J, see
   * the paper) */
  compute_d(d, m_J, np);
  //    update_z(z, m_J, d, iq);
  if (iq >= nVars) {
    //      throw std::runtime_error("iq >= m_J.cols()");
    z.setZero();
  } else {
    update_z(z, m_J, d, iq);
  }
  /* compute N* np (if q > 0): the negative of the step direction in the dual
   * space */
  update_r(R, r, d, iq);
#ifdef TRACE_SOLVER
  std::cout << "Step direction z" << std::endl;
  print_vector("z", z, nVars);
  print_vector("r", r, iq + 1);
  print_vector("u", u, iq + 1);
  print_vector("d", d, nVars);
  print_vector("A", A, iq + 1);
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2A);

  /* Step 2b: compute step length */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2B);
  l = 0;
  /* Compute t1: partial step length (maximum step in dual space without
   * violating dual feasibility */
  t1 = inf; /* +inf */
  /* find the index l s.t. it reaches the minimum of u+(x) / r */
  // l: index of constraint to drop (maybe)
  for (k = nEqCon; k < iq; k++) {
    double tmp;
    if (r(k) > 0.0 && ((tmp = u(k) / r(k)) < t1)) {
      t1 = tmp;
      l = A(k);
    }
  }
  /* Compute t2: full step length (minimum step in primal space such that the
   * constraint ip becomes feasible */
  if (std::abs(z.dot(z)) >
      std::numeric_limits<double>::epsilon())  // i.e. z != 0
    t2 = -s(ip) / z.dot(np);
  else
    t2 = inf; /* +inf */

  /* the step is chosen as the minimum of t1 and t2 */
  t = std::min(t1, t2);
#ifdef TRACE_SOLVER
  std::cout << "Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2
            << ") ";
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2B);

  /* Step 2c: determine new S-pair and take step: */
  START_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);
  /* case (i): no step in primal or dual space */
  if (t >= inf) {
    /* QPP is infeasible */
    q = iq;
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);
    return EIQUADPROG_FAST_UNBOUNDED;
  }
  /* case (ii): step in dual space */
  if (t2 >= inf) {
    /* set u = u +  t * [-r 1) and drop constraint l from the active set A */
    u.head(iq) -= t * r.head(iq);
    u(iq) += t;
    iai(l) = l;
    delete_constraint(R, m_J, A, u, nEqCon, iq, l);
#ifdef TRACE_SOLVER
    std::cout << " in dual space: " << f_value << std::endl;
    print_vector("x", x, nVars);
    print_vector("z", z, nVars);
    print_vector("A", A, iq + 1);
#endif
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);
    goto l2a;
  }

  /* case (iii): step in primal and dual space */
  x += t * z;
  /* update the solution value */
  f_value += t * z.dot(np) * (0.5 * t + u(iq));

  u.head(iq) -= t * r.head(iq);
  u(iq) += t;

#ifdef TRACE_SOLVER
  std::cout << " in both spaces: " << f_value << std::endl;
  print_vector("x", x, nVars);
  print_vector("u", u, iq + 1);
  print_vector("r", r, iq + 1);
  print_vector("A", A, iq + 1);
#endif

  if (t == t2) {
#ifdef TRACE_SOLVER
    std::cout << "Full step has taken " << t << std::endl;
    print_vector("x", x, nVars);
#endif
    /* full step has taken */
    /* add constraint ip to the active set*/
    if (!add_constraint(R, m_J, d, iq, R_norm)) {
      iaexcl(ip) = 0;
      delete_constraint(R, m_J, A, u, nEqCon, iq, ip);
#ifdef TRACE_SOLVER
      print_matrix("R", R, nVars);
      print_vector("A", A, iq);
#endif
      for (i = 0; i < nIneqCon; i++) iai(i) = i;
      for (i = 0; i < iq; i++) {
        A(i) = A_old(i);
        iai(A(i)) = -1;
        u(i) = u_old(i);
      }
      x = x_old;
      STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);
      goto l2; /* go to step 2 */
    } else
      iai(ip) = -1;
#ifdef TRACE_SOLVER
    print_matrix("R", R, nVars);
    print_vector("A", A, iq);
#endif
    STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);
    goto l1;
  }

  /* a partial step has been taken => drop constraint l */
  iai(l) = l;
  delete_constraint(R, m_J, A, u, nEqCon, iq, l);
  s(ip) = CI.row(ip).dot(x) + ci0(ip);

#ifdef TRACE_SOLVER
  std::cout << "Partial step has taken " << t << std::endl;
  print_vector("x", x, nVars);
  print_matrix("R", R, nVars);
  print_vector("A", A, iq);
  print_vector("s", s, nIneqCon);
#endif
  STOP_PROFILER_EIQUADPROG_FAST(EIQUADPROG_FAST_STEP_2C);

  goto l2a;
}

} /* namespace solvers */
} /* namespace tsid */
