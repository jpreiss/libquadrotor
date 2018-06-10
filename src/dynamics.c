#include "dynamics.h"

// simulate the dynamics of the quadrotor.
void quad_dynamics(
	struct quad_physical_params const *params,
	struct quad_state const *now,
	struct accel const *force,
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
		vscl((1.0f - 0.1f * dt), omega), // TODO: make this a parameter
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

