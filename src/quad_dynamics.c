#include "quad_dynamics.h"

// simulate the dynamics of the quadrotor.
void quad_dynamics(
	struct quad_physical_params const *params,
	struct quad_state const *now,
	struct quad_accel const *force,
	float const dt,
	struct quad_state *next)
{
	// rotational dynamics
	struct vec const omega_dot = veltdiv(
		vadd(
			vcross(
				vneg(now->omega),
				veltmul(params->inertia, now->omega)
			),
			force->angular
		),
		params->inertia
	);
	struct vec const omega = vadd(
		vscl((1.0f - 0.1f * dt), now->omega), // TODO: make this a parameter
		vscl(dt, omega_dot)
	);
	next->quat = quat_gyro_update(now->quat, omega, dt);
	next->omega = omega;
	next->ticks_since_qnormalize = now->ticks_since_qnormalize + 1;
	// TODO make this a parameter? do based on qdot? do always?
	if (now->ticks_since_qnormalize >= 10) {
		next->quat = qnormalize(next->quat);
		next->ticks_since_qnormalize = 0;
	}

	// translational dynamics
	// TODO: RK4? Verlet?
	float const drag = params->drag * vmag2(now->vel);
	struct vec const acc = vadd3(
		qvrot(next->quat, mkvec(0.0f, 0.0f, force->linear / params->mass)),
		mkvec(0.0f, 0.0f, -GRAV),
		vscl(-drag / params->mass, now->vel)
	);
	next->vel = vadd(now->vel, vscl(dt, acc));
	next->pos = vadd(now->pos, vscl(dt, next->vel));
	next->acc = acc;
}


// note: this computes the actual moment and thrust,
// not the acc + gravity and angular acceleration that the controllers output
void quad_motor_forces(
	struct quad_physical_params const *params,
	float const motors[4],
	struct quad_accel *force)
{
	force->linear = motors[0] + motors[1] + motors[2] + motors[3];

	if (params->motor_layout == 'x') {
		float const arm = 0.707106781f * params->arm_length;
		force->angular.x = arm * (motors[0] - motors[1] - motors[2] + motors[3]);
		force->angular.y = arm * (-motors[0] - motors[1] + motors[2] + motors[3]);
	}
	else {
		float const arm = params->arm_length;
		force->angular.x = arm * (-motors[1] + motors[3]);
		force->angular.y = arm * (-motors[0] + motors[2]);
	}

	float const yawscl =
		params->thrust_to_torque * (params->motor_0_ccw ? -1.0f : 1.0f);
	force->angular.z = yawscl * (motors[0] - motors[1] + motors[2] - motors[3]);
}
