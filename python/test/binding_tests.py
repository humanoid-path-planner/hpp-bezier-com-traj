from centroidal_dynamics import *

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

z = array([0.,0.,1.])
P = asmatrix(array([array([x,y,0]) for x in [-0.05,0.05] for y in [-0.1,0.1]]))
N = asmatrix(array([z for _ in range(4)]))

#setting contact positions and normals, as well as friction coefficients 
eq.setNewContacts(asmatrix(P),asmatrix(N),0.3,EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_LP)

c= asmatrix(array([0.,0.,1.]))

#computing robustness of a given configuration, first with no argument (0 acceleration, static equilibrium)
status, robustness = eq.computeEquilibriumRobustness(c)
assert (status == LP_STATUS_OPTIMAL), "LP should not fail"
assert (robustness > 0), "first test should be in equilibrirum"
	
#computing robustness of a given configuration with non zero acceleration
ddc= asmatrix(array([1000.,0.,0.]))
status, robustness = eq.computeEquilibriumRobustness(c,ddc)
assert (status == LP_STATUS_OPTIMAL), "LP should not fail"
assert (robustness < 0), "first test should NOT be in equilibrirum"

#now, use polytope projection algorithm
eq.setNewContacts(asmatrix(P),asmatrix(N),0.3,EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_PP)
H,h = eq.getPolytopeInequalities()

