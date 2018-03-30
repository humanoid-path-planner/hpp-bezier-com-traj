#  bezier_COM_Traj

Copyright 2018 LAAS-CNRS

Author: Steve Tonneau

##Description
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
* [centroidal-dynamics-lib] Centroidal dynamics computation library (https://github.com/stonneau/centroidal-dynamics-lib)
* [spline] Bezier curves library (https://github.com/stonneau/spline)

## Additional dependencies for python bindings
* [Boost.Python](http://www.boost.org/doc/libs/1_63_0/libs/python/doc/html/index.html)
* [eigenpy](https://github.com/stack-of-tasks/eigenpy)
* Additionally you will need to activate the python bindings for the above libraries

##Installation on ubuntu-14.04 64 bit

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
To see example of use, you can refer to the [test file](https://gitlab.com/stonneau/bezier_COM_traj/blob/master/python/test/binding_tests.py)

In spite of an exhaustive documentation, please refer to the C++ documentation, which mostly applies
to python.
