# libquadrotor
dependency-free C implementations of quadrotor control, state estimation, and dynamics.

The current main components are:

 - Nonlinear trajectory tracking controller
 - Nonlinear attitude tracking controller
 - Power distribution for "X" and "Plus" motor arrangements
 - First-order dynamics simulation
 - Extended Kalman Filter (EKF) fusion of IMU with motion capture

Possible future extensions include EKF support for more types of sensors, trajectory planning, more detailed physics model, power distribution for other kinds of multirotors, etc.

The goal of this project is to implement the core functions of quadrotor operation and simulation
in a form that can be used both on the onboard flight controller and on a PC for testing, simulation, machine learning, etc.
In pursuit of this goal, the code has the following characteristics:

- Written in C99 with no compiler extensions
- Uses none of the standard library except `<math.h>` and `<stdbool.h>` (not including test code)
- No external dependencies except the single header [`cmath3d`](https://github.com/jpreiss/cmath3d)
- No dynamic memory allocation
- Uses no global state

A CPU with hardware floating point operations is assumed.
Small structs are passed by value more often than in old-fashioned C code,
so the C compiler must be reasonably good at copy elision, destructuring, etc.
