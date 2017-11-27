from centroidal_dynamics import *
from spline import *
from bezier_com_traj import *

#testing constructors
eq = Equilibrium("test", 54., 4) 
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES ) 
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES ) 
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES, False) 
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES, False, 1) 
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES, True, 1, True ) 

#whether useWarmStart is enable (True by default)
previous = eq.useWarmStart()
#enable warm start in solver (only for QPOases)
eq.setUseWarmStart(False)
assert(previous != eq.useWarmStart())

#access solver name
assert(eq.getName() == "test")


# creating contact points
from numpy import array, asmatrix, matrix
import numpy as np

z = array([0.,0.,1.])
P = asmatrix(array([array([x,y,0]) for x in [-0.05,0.05] for y in [-0.1,0.1]]))
N = asmatrix(array([z for _ in range(4)]))

#setting contact positions and normals, as well as friction coefficients 
eq.setNewContacts(asmatrix(P),asmatrix(N),0.3,EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_PP)
#~ eq.setNewContacts(asmatrix(P),asmatrix(N),0.3,EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_LP)

# setting up optimization problem
c0 = matrix([0.,0.,1.]) 
#~ dc0 = matrix(np.random.uniform(-1, 1, size=3)); 
dc0 =  matrix([0.1,0.,0.]) 
l0 = matrix([0.,0.,0.]) 
T = 1.2
tstep = -1.

a = zeroStepCapturability(eq,c0,dc0,l0,False,T,tstep)

assert(a.success)
a.c_of_t(0)
a.dL_of_t(T)

Kin = matrix(np.identity(3))
kin = 10*np.ones(3); 
a = zeroStepCapturability(eq,c0,dc0,l0,False,T,tstep,Kin,matrix(kin))
assert(a.success)

kin[2] = 0.5
a = zeroStepCapturability(eq,c0,dc0,l0,False,T,tstep,Kin,matrix(kin))
assert(np.asarray(a.x[2])[0][0] <=0.5)


a = zeroStepCapturability(eq,c0,dc0,l0,True,T,tstep,Kin,matrix(kin))
