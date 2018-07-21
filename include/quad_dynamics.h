#pragma once

#include "quad_common_types.h"

// simulate the dynamics of the quadrotor.
void quad_dynamics(
	struct quad_physical_params const *params,
	struct quad_state const *now,
	struct quad_accel const *acc,
	float const dt,
	struct quad_state *next);

// given the thrust in newtons produced by each of the 4 motors,
// compute the magnitude of thrust (normalized by mass)
// and the angular acceleration.
void quad_motor_forces(
	struct quad_physical_params const *params,
	float const motors[4],
	struct quad_accel *force);
