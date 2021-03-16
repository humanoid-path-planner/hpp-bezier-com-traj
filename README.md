#  bezier_COM_Traj

[![Pipeline status](https://gepgitlab.laas.fr/humanoid-path-planner/hpp-bezier-com-traj/badges/master/pipeline.svg)](https://gepgitlab.laas.fr/humanoid-path-planner/hpp-bezier-com-traj/commits/master)
[![Coverage report](https://gepgitlab.laas.fr/humanoid-path-planner/hpp-bezier-com-traj/badges/master/coverage.svg?job=doc-coverage)](http://projects.laas.fr/gepetto/doc/humanoid-path-planner/hpp-bezier-com-traj/master/coverage/)


Copyright 2018-2020 LAAS-CNRS

Authors: Pierre Fernbach and Steve Tonneau

## Description
bezier_COM_Traj implements tools to compute Bezier trajectories given various sets of constraints: initial and terminal conditions (position, velocities, acceleration), additional linear constraints on the complete trajectory, and, most interestingly, constraints related to the center of mass dynamics.

The trajectories are genererated through the resolution of convex optimization (Quadratic Programms), and thus allow to specify a cost functional to minimize.

The library is implemented in C++, but also provides Python bindings.

Two types of applications can be used so far:
- First, zero step capturability: Given the centroidal state of a robot, determines whether it is possible for the robot to come to a stop without violating frictional constraints. In this formulation, the problem can be solved continuously, and angular momentum constraints can be used.

- Second, the general case (which encompasses zero step capturability):
Given a sequence of discrete contact configurations, and given the current state of the robot, and a desired target state, compute a kinematically and dynamically accurate trajectory for the center of mass of the robot. In this general case, the trajectory is checked at discrete intervals (the verification is not continuous). Furthermore, at this point angular momentum is not handled (this is a TODO and not a limitation of the approach).

More details can be found in the preprint paper:
CROC: Convex Resolution Of Centroidal dynamics trajectories to provide a feasibility criterion for the multi contact planning problem, by Fernbach et al.
https://hal.archives-ouvertes.fr/hal-01726155v1


## Dependencies
* [centroidal-dynamics-lib](https://github.com/stonneau/centroidal-dynamics-lib) Centroidal dynamics computation library
* [ndcurves](https://github.com/loco-3d/ndcurves) Bezier curves library
* [glpk](https://www.gnu.org/software/glpk/) GNU Linear Programming Kit

## Additional dependencies for python bindings
* [Boost.Python](http://www.boost.org/doc/libs/1_63_0/libs/python/doc/html/index.html)
* [eigenpy](https://github.com/stack-of-tasks/eigenpy)
* Additionally you will need to activate the python bindings for the above libraries

## Installation on ubuntu-14.04 64 bit

Once the required libraries are installed you can clone this repository using ssh:
```
git clone --recursive git@gitlab.com:stonneau/bezier_COM_traj.git $BEZIER_COM_DIR
```
or using http:
```
git clone --recursive https://gitlab.com/stonneau/bezier_COM_traj.git $BEZIER_COM_DIR
```
And you can build this library using CMake:
```
mkdir $BEZIER_COM_DIR/build
cd $BEZIER_COM_DIR/build
cmake -DCMAKE_INSTALL_PREFIX=${DEVEL_DIR}/install ..
make install
```

### Optional: Python bindings installation
To install the Python bindings, in the CMakeLists.txt file, first enable the BUILD_PYTHON_INTERFACE option:
```
OPTION (BUILD_PYTHON_INTERFACE "Build the python binding" ON)
```

Then rebuild the library:
```
cd $BEZIER_COM_DIR/build
cmake -DCMAKE_INSTALL_PREFIX=${DEVEL_DIR}/install ..
make install
```
The python bindings should then be accessible through the package bezier_com_traj.
To see all the possible uses, you can refer to the [test file](https://gitlab.com/stonneau/bezier_COM_traj/blob/master/python/test/binding_tests.py)

In spite of an exhaustive documentation, please refer to the C++ documentation, which mostly applies
to python.

## Python example : Zero step capturability

For the zero step capturability, we will first define a contact phase using the objects from centroidal_dynamics:
```
#importing the libraries of interest
import ndcurves  # noqa - necessary to register ndcurves::bezier_curve
import numpy as np
from numpy import array
from hpp_centroidal_dynamics import Equilibrium, EquilibriumAlgorithm, SolverLP
from hpp_bezier_com_traj import (SOLVER_QUADPROG, ConstraintFlag, Constraints, ContactData, ProblemData,
                                 computeCOMTraj, zeroStepCapturability)


# create an Equilibrium solver, for a robot of 54 kilos. We linearize the friction cone to four generating rays
eq = Equilibrium("test", 54., 4)

# Now define some contact points ...
P = array([[x, y, 0] for x in [-0.05, 0.05] for y in [-0.1, 0.1]])


#and normals
z = array([0., 0., 1.])
N = array([[0., 0., 1.]] * 4)


#setting contact positions and normals, as well as friction coefficient of 0.3
#EQUILIBRIUM_ALGORITHM_PP is the algorithm that will always be used for our problems
eq.setNewContacts(P, N, 0.3, EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_PP)

```

Then, we will define the initial state of our robot:
```
#c0 is the initial center of mass position
c0 = array([0., 0., 1.])


#we set the inital speed dc0 to a rather slow 10 cm / s along the x axis
dc0 = array([0.1, 0., 0.])
l0 = array([0., 0., 0.])
```

And finally, some optimization parameters:
The total duration of the trajectory, as well as
the discretization step. If the discretization step is < 0,
then the continuous formulation is used

```
#trajectory duration of 1.2 seconds
T = 1.2
#continuous resolution of the trajectory
tstep = -1.
```

We can now solve the problem:
```
# the boolean value indicates whether to use or not angular momentum
result = zeroStepCapturability(eq, c0, dc0, l0, False, T, tstep)
print(result.success)
#True the problem was feasible, and a trajectory was successfully computed

```

The found centroidal trajectory is accessible from the returned object, only if the problem
was feasible
```
result.c_of_t # a bezier curve object describing the com trajectory

#We can check that the end velocity is indeed zero:
dc_of_t = result.c_of_t.compute_derivate(1) # computing first derivative
print(np.linalg.norm(dc_of_t(dc_of_t.max())))
# 0.0

```

refer to the [test file](https://gitlab.com/stonneau/bezier_COM_traj/blob/master/python/test/binding_tests.py) for more advanced problems, including kinematic constraints,
mutiple contact phases handling and angular momentum
