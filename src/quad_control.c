#include "quad_control.h"
#include <math.h>

//
// local helper functions.
//

static inline void integrate_clamp(struct vec *integral,
	struct vec const *bound, struct vec const *value, float dt)
{
	*integral = vclampabs(vadd(*integral, vscl(dt, *value)), *bound);
}

static inline struct vec xy_zmul(struct xy_z xy_z, struct vec v)
{
	return mkvec(xy_z.xy * v.x, xy_z.xy * v.y, xy_z.z * v.z);
}

static struct vec qzaxis(struct quat q)
{
	float x = q.x, y = q.y, z = q.z, w = q.w;
	return mkvec(
		2*x*z + 2*y*w,
		2*y*z - 2*x*w,
		1 - 2*x*x - 2*y*y
	);
}

static float sign(float x) {
	return (x >= 0.0f) ? 1.0f : -1.0f;
}

//
// library functions, see quad_control.h for descriptions.
//

void quad_ctrl_default_params(struct quad_ctrl_params *params)
{
	struct quad_ctrl_params p = {
		.linear = {
			.kp = {.xy = 12.5f, .z = 39.0f},
			.ki = {.xy = 1.55f, .z = 1.55f},
			.kd = {.xy = 6.25f, .z = 12.5f},
			.int_bound = {.x = 2.0f, .y = 2.0f, .z = 0.4f},
		},
		// TODO: does it make sense to integrate omega error,
		// when it's basically the same as the attitude error?
		
		// the CF implementation integrates the attitude error eR,
		// so that's like an actual integral term on attitude...
		// should we do that too?
		// but it's a little weird because attitude is already
		// the double integral of the controlled quantity,
		// and it's unusual to control based on a triple integral...
		// ... if we do it, then i_range_m_z = 1500 in the firmware is huge
		//     and probably a mistake.

		// TODO: gyroscope lowpass filter in EKF

		// TODO: should we really have a d term on omega?
		// some pilot-oriented flight controllers have it...
		.attitude = {
			.kp = {.xy = 130.3f, .z = 111.7f},
			.ki = {.xy = 0.00f, .z = 0.00f},
			.kd = {.xy = 18.6f, .z = 11.2f},
			.int_bound = {.x = 1.0f, .y = 1.0f, .z = 0.0f}, // TODO
		},
	};
	*params = p;
}


void physical_params_crazyflie2(struct quad_physical_params *params)
{
	// from "System Identification of the Crazyflie 2.0 Nano Quadrocopter".
	// J. Foerster, M. Hamer, R. D'Andrea. Bachelor Thesis, ETH Zurich, 2015.
	struct quad_physical_params p = {
		.mass = 0.027f,
		.arm_length = 0.046f,
		.inertia = { 1.66e-5f, 1.66e-5f, 2.92e-5f },
		.thrust_to_torque = 0.006f,
		.max_thrust = 0.15f,
		.drag = 0.0f, // TODO
		.motor_layout = 'x',
		.motor_0_ccw = true, // TODO
	};
	*params = p;
}


void quad_ctrl_init(struct quad_ctrl_state *state)
{
	state->int_linear_err = vzero();
	state->int_attitude_err = vzero();
}


struct quad_accel quad_ctrl_full(
	struct quad_ctrl_state *state,
	struct quad_ctrl_params const *param,
	struct quad_state const *s, struct quad_state const *set, float dt)
{
	// -------------------- Linear part --------------------
	struct vec const pos_error = vsub(set->pos, s->pos);
	struct vec const vel_error = vsub(set->vel, s->vel);
	integrate_clamp(&state->int_linear_err, &param->linear.int_bound, &pos_error, dt);
	struct vec const target_acc = vadd4(
		set->acc,
		xy_zmul(param->linear.kp, pos_error),
		xy_zmul(param->linear.ki, state->int_linear_err),
		xy_zmul(param->linear.kd, vel_error)
	);
	struct vec const target_thrust = vadd(target_acc, mkvec(0.0f, 0.0f, GRAV));

	struct vec const z_axis = qzaxis(s->quat);
	float const thrust_mag = vdot(target_thrust, z_axis);

	// -------------------- Angular part --------------------
	// construct the desired attitude, avoiding yaw singularity by using quats.
	// TODO: we could carry yaw angle as a setpoint directly,
	// so the user doesn't need to manually construct a quat setpoint
	float const target_yaw = quat2rpy(set->quat).z;
	struct quat q_yaw = qaxisangle(mkvec(0,0,1), target_yaw);
	struct quat q_thrust = qvectovec(mkvec(0,0,1), vnormalize(target_thrust));
	struct quat q_des = qqmul(q_thrust, q_yaw);
	struct quad_state set2;
	set2.quat = q_des;
	set2.omega = set->omega;

	struct quad_accel output = {
		.angular = quad_ctrl_attitude(state, param, s, &set2, dt),
		.linear = thrust_mag
	};
	return output;
}


struct vec quad_ctrl_attitude(
	struct quad_ctrl_state *state,
	struct quad_ctrl_params const *param,
	struct quad_state const *s, struct quad_state const *set, float dt)
{
	// this controller implements the paper:
	// "Nonlinear Quadrocopter Attitude Control Technical Report"
	// D. Brescianini, M. Hehn, R. D'Andrea., ETH Zurich, 2013.
	// notation is mostly same but we use standard quaternion multiplication order
	// and do more computations in body frame instead of inertial frame.

	// vector in inertial frame that gives desired thrust direction with no regard for yaw.
	struct vec const zbody_des_world = qzaxis(set->quat);

	// expression of zbody_des in body frame.
	struct vec const zbody_des = qvrot(qinv(s->quat), zbody_des_world);

	// rotation in body frame to reach thrust-only attitude.
	// TODO: can simplify by manual const-propagation of [0,0,1]
	//struct quat const q_reduced = qnormalize(qvectovec(mkvec(0,0,1), zbody_des));
	struct quat const q_reduced = qvectovec(mkvec(0,0,1), zbody_des);

	// rotation in body frame to full setpoint attitude.
	struct quat const q_full = qqmul(set->quat, qinv(s->quat));

	// rotation from q_reduced to q_full.
	struct quat q_mix = qnormalize(qqmul(q_full, qinv(q_reduced)));

	// fix double-covering by forcing q_mix to be a rotation about +z
	// so the angle expressed by q_mix can be compared to current z angular velocity
	if (q_mix.z < 0) {
		q_mix = qneg(q_mix);
	}
	// angle of yaw rotation
	float const alpha_mix = quat2angle(q_mix);

	float const yaw_priority = param->attitude.kp.z / param->attitude.kp.xy;

	// if body yaw rate is high, we should prefer to do a yaw correction
	// in the same direction even if it's more than 180 degrees.
	// copied verbatim from ETHZ paper. would like to simplify a little --
	// do we really need the trig and all the pi factors?
	// also, unless yaw priority is quite low, the threshold ends up being huge
	float const omega_z_min = param->attitude.kp.xy * sinf(M_PI_2_F * yaw_priority);
	float const threshold = ((M_PI_F - alpha_mix) * M_1_PI_F) * omega_z_min;
	float const omega_clash = sign(q_mix.w) * s->omega.z;

	float yaw_correction;
	if (omega_clash < -threshold) {
		// our yaw angular velocity is opposite of desired attitude correction,
		// but it's fast enough that we prefer continuing in this direction
		// even though we'll need to rotate more than 180 degrees.
		// this correction is only needed for yaw because control authority is low.
		yaw_correction = M_PI_F * sign(s->omega.z);
	}
	else {
		yaw_correction = alpha_mix;
	}
	yaw_correction *= yaw_priority;
	struct quat const q_yaw = qaxisangle(mkvec(0, 0, 1), yaw_correction);

	// final desired rotation
	struct quat const q_err = qqmul(q_yaw, q_reduced);

	// scale such that the desired angular acceleration for
	// 90 degrees error is exactly params.angular.kp
	float const scale = (2.0f / sqrtf(2.0f));
	struct vec const moment = vadd(
		vscl(param->attitude.kp.xy * scale, quatimagpart(qposreal(q_err))),
		xy_zmul(param->attitude.kd, vsub(set->omega, s->omega))
	);
	return moment;
}


struct vec quad_ctrl_attitude_rate(
	struct quad_ctrl_state *state,
	struct quad_ctrl_params const *param,
	struct vec s, struct vec set, float thrust, float dt)
{
	struct vec const omega_error = vsub(set, s);
	integrate_clamp(&state->int_attitude_err, &param->attitude.int_bound, &omega_error, dt);
	// TODO add feedforward term for angular acceleration in input?
	// TODO implement derivative term for angular velocity error
	struct vec moment = vadd(
		xy_zmul(param->attitude.kd, omega_error),
		xy_zmul(param->attitude.ki, state->int_attitude_err)
	);
	return moment;
}


void quad_power_distribute(
	struct quad_accel const *acc,
	struct quad_physical_params const *params,
	float prop_forces[4])
{
	struct vec const moment = veltmul(params->inertia, acc->angular);
	float const thrustpart = 0.25f * (params->mass * acc->linear);
	float const yawpart = 0.25f * (params->motor_0_ccw ? -1.0f : 1.0f) *
		moment.z / params->thrust_to_torque;

	if (params->motor_layout == 'x') {
		float const arm = 0.707106781f * params->arm_length;
		struct vec const moment_scl = vscl(0.25f / arm, moment);
		prop_forces[0] = thrustpart + moment_scl.x - moment_scl.y + yawpart;
		prop_forces[1] = thrustpart - moment_scl.x - moment_scl.y - yawpart;
		prop_forces[2] = thrustpart - moment_scl.x + moment_scl.y + yawpart;
		prop_forces[3] = thrustpart + moment_scl.x + moment_scl.y - yawpart;
	}
	else {
		struct vec const moment_scl = vscl(0.5f / params->arm_length, moment);
		prop_forces[0] = thrustpart - moment_scl.y + yawpart;
		prop_forces[1] = thrustpart - moment_scl.x - yawpart;
		prop_forces[2] = thrustpart + moment_scl.y + yawpart;
		prop_forces[3] = thrustpart + moment_scl.x - yawpart;
	}

	for (int i = 0; i < 4; ++i) {
		prop_forces[i] = clamp(prop_forces[i], 0.0f, params->max_thrust);
	}
}

