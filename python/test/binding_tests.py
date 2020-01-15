import curves  # noqa - necessary to register curves::bezier_curve
import numpy as np
from hpp_centroidal_dynamics import Equilibrium, EquilibriumAlgorithm, SolverLP
from numpy import array, asmatrix, matrix

from hpp_bezier_com_traj import (SOLVER_QUADPROG, ConstraintFlag, Constraints, ContactData, ProblemData,
                                 computeCOMTraj, zeroStepCapturability)

# testing constructors
eq = Equilibrium("test", 54., 4)
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES)
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES)
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES, False)
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES, False, 1)
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES, True, 1, True)

# whether useWarmStart is enable (True by default)
previous = eq.useWarmStart()
# enable warm start in solver (only for QPOases)
eq.setUseWarmStart(False)
assert (previous != eq.useWarmStart())

# access solver name
assert (eq.getName() == "test")

z = array([0., 0., 1.])
P = asmatrix(array([array([x, y, 0]) for x in [-0.05, 0.05] for y in [-0.1, 0.1]]))
N = asmatrix(array([z for _ in range(4)]))

# setting contact positions and normals, as well as friction coefficients
eq.setNewContacts(asmatrix(P), asmatrix(N), 0.3, EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_PP)
# eq.setNewContacts(asmatrix(P),asmatrix(N),0.3,EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_LP)

# setting up optimization problem
c0 = matrix([0., 0., 1.]).T
# dc0 = matrix(np.random.uniform(-1, 1, size=3));
dc0 = matrix([0.1, 0., 0.]).T
l0 = matrix([0., 0., 0.]).T
T = 1.2
tstep = -1.

a = zeroStepCapturability(eq, c0, dc0, l0, False, T, tstep)

assert (a.success)
a.c_of_t(0)
a.dL_of_t(T)

Kin = matrix(np.identity(3))
kin = 10 * np.ones(3)
# TODO: Invalid sizes when resizing a matrix or array.
# a = zeroStepCapturability(eq, c0, dc0, l0, False, T, tstep, Kin, matrix(kin))
# assert (a.success)
#
# kin[2] = 0.5
# a = zeroStepCapturability(eq, c0, dc0, l0, False, T, tstep, Kin, matrix(kin))
# assert (np.asarray(a.x[2])[0][0] <= 0.5)
#
# a = zeroStepCapturability(eq, c0, dc0, l0, True, T, tstep, Kin, matrix(kin))

# testing contactData
cData = ContactData(Equilibrium("test", 54., 4))
ceq = cData.contactPhase_
ceq.setNewContacts(asmatrix(P), asmatrix(N), 0.3, EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_PP)
assert cData.contactPhase_.getAlgorithm(
) == EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_PP, "modifying ceq should modify cData.contactPhase_"

Id = matrix([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])

excep = False
try:
    cData.Kin_
except RuntimeError:
    excep = True
assert excep, "[ERROR] No kin assigned should have raised exception"
cData.setKinematicConstraints(Id, matrix([0., 0., 1.]).T)
cData.Kin_

excep = False
try:
    cData.setKinematicConstraints(Id, matrix([0., 0., 0., 1.]).T)
except RuntimeError:
    excep = True
assert excep, "[ERROR] Miss matching matrix and vector should raise an error"

excep = False
try:
    cData.Ang_
except RuntimeError:
    excep = True
assert excep, "[ERROR] No Ang_ assigned should have raised exception"
cData.setAngularConstraints(Id, matrix([0., 0., 1.]).T)
cData.Ang_

excep = False
try:
    cData.setAngularConstraints(Id, matrix([0., 0., 0., 1.]).T)
except RuntimeError:
    excep = True
assert excep, "[ERROR] Missmatching matrix and vector should raise an error"

# testing constraints
c = Constraints()
old = c.constrainAcceleration_
c.constrainAcceleration_ = not old
assert c.constrainAcceleration_ != old
old = c.flag_
assert c.flag_ == ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL | ConstraintFlag.END_VEL | ConstraintFlag.END_POS
c.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL
assert c.flag_ != old
old = c.maxAcceleration_
c.maxAcceleration_ = .235
assert c.maxAcceleration_ != old
old = c.reduce_h_
c.reduce_h_ = .235
assert c.reduce_h_ != old

# testing problem data
c = ProblemData()

nv = matrix([0., 0., 10.]).T
old = c.c0_
c.c0_ = nv
assert (c.c0_ != old).any()
old = c.dc0_
c.dc0_ = nv
assert (c.dc0_ != old).any()
old = c.ddc0_
c.ddc0_ = nv
assert (c.ddc0_ != old).any()
old = c.c1_
c.c1_ = nv
assert (c.c0_ != old).any()
old = c.dc1_
c.dc1_ = nv
assert (c.dc1_ != old).any()
old = c.ddc1_
c.ddc1_ = nv
assert (c.ddc1_ != old).any()
old = c.useAngularMomentum_
c.useAngularMomentum_ = not old
assert c.useAngularMomentum_ != old
pD = c
c = pD.constraints_
old = c.flag_
c.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL
assert pD.constraints_.flag_ != old

pD = ProblemData()
pD.constraints_.flag_ = ConstraintFlag.INIT_POS | ConstraintFlag.INIT_VEL | ConstraintFlag.END_VEL


def initContactData(pD):
    cData = ContactData(Equilibrium("test", 54., 4))
    cData.contactPhase_.setNewContacts(asmatrix(P), asmatrix(N), 0.3, EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_PP)
    pD.addContact(cData)


[initContactData(pD) for i in range(3)]

pD.c0_ = c0
pD.dc0_ = dc0
res = computeCOMTraj(pD, matrix([0.4, 0.4, 0.4]).T, 0.05)
# test glpk only if defined
try:
    test = SOLVER_GLPK
    res = computeCOMTraj(pD, matrix([0.4, 0.4, 0.4]).T, 0.05, SOLVER_GLPK)
except NameError:
    print("[WARNING] SOLVER_GLPK is not defined.")
    print("Consider installing GLPK if you are using CROC with a force formulation")

res = computeCOMTraj(pD, matrix([0.4, 0.4, 0.4]).T, 0.05, SOLVER_QUADPROG)
#  res = computeCOMTraj(pD,matrix([0.4,0.4,0.4]).T,0.05,SOLVER_QUADPROG_SPARSE)
assert np.linalg.norm(res.c_of_t.derivate(1.2, 1)) < 0.00000001

# non matching time step and contact phases
excep = False
try:
    res = computeCOMTraj(pD, matrix([0.4, 0.4]).T, 0.05)
except RuntimeError:
    excep = True
assert excep, "[ERROR] computeCOMTraj should have raised exception"

print("all tests passed")
