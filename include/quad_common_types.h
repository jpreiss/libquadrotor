#pragma once

#include "math3d.h"

static float const GRAV = 9.81;

// the mutable physical state of the quadrotor.
// this is used for both setpoints/goals and actual states.
struct quad_state
{
	struct vec pos;
	struct vec vel;
	struct vec acc;

	struct quat quat; // from body frame to world frame.
	struct vec omega; // radians per second in body frame.
	int ticks_since_qnormalize;
};
static inline void quad_zero_state(struct quad_state *s)
{
	s->pos = s->vel = s->acc = s->omega = vzero();
	s->quat = qeye();
	s->ticks_since_qnormalize = 0;
}

// the thrust magnitude in the body "up" direction
// and the angular acceleration.
// the boundary between controller and power distribution.
// also the quantities that drive the dynamics simulation.
struct quad_accel
{
	float linear;
	struct vec angular;
};

// the immutable physical parameters of the quadrotor.
struct quad_physical_params
{
	float mass; // kg
	float arm_length; // length of one arm. m
	struct vec inertia; // diagonal of axis-aligned ineria matrix. kg * m^2
	float thrust_to_torque; // property of the propeller. torque = this * thrust.
	float max_thrust; // max thrust per motor in newtons
	float drag; // constant c such that drag force is -c * |v|^2 * \hat{v}
	            // TODO: non-isotropic drag? angular drag?

	// motor_layout is one of 'x' or '+'.
	// if 'x', motor order starts at forward left (+x +y) and goes clockwise.
	// if '+', motor order starts at forward (+x) and goes clockwise.
	// motor_0_ccw is set true if the first motor in this order spins counterclockwise.
	// assumes diagonally opposed motors spin in same direction.
	char motor_layout;
	bool motor_0_ccw;
};

// Notes on indentifying quad_physical_params:
// - if thrust-to-weight ratio of quadrotor is known, then:
//   thrust_max = GRAV * t/w * mass / 4
