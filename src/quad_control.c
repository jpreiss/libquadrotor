#include "quad_control.h"
#include <math.h>
#include <assert.h>

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

struct vec vee_map_rot_error(struct mat33 const *R, struct mat33 const *Rdes)
{
	struct mat33 R_transpose = mtranspose(*R);
	struct mat33 Rdes_transpose = mtranspose(*Rdes);
	struct mat33 eRM = msub(mmul(Rdes_transpose, *R), mmul(R_transpose, *Rdes));
	struct vec eR = vscl(0.5, mkvec(eRM.m[2][1], eRM.m[0][2], eRM.m[1][0]));
	return eR;
}

struct vec vee_map_rot_error_quat2mat(struct quat const *quat, struct mat33 const *Rdes)
{
	// (generated using Mathematica, see paper for details)
	struct vec eR;
	float x = quat->x;
	float y = quat->y;
	float z = quat->z;
	float w = quat->w;
	struct vec xb_des = mcolumn(*Rdes, 0);
	struct vec yb_des = mcolumn(*Rdes, 1);
	struct vec zb_des = mcolumn(*Rdes, 2);
	eR.x = (-1 + 2*fsqr(x) + 2*fsqr(y))*yb_des.z + zb_des.y - 2*(x*yb_des.x*z + y*yb_des.y*z - x*y*zb_des.x + fsqr(x)*zb_des.y + fsqr(z)*zb_des.y - y*z*zb_des.z) +2*w*(-(y*yb_des.x) - z*zb_des.x + x*(yb_des.y + zb_des.z));
	eR.y = xb_des.z - zb_des.x - 2*(fsqr(x)*xb_des.z + y*(xb_des.z*y - xb_des.y*z) - (fsqr(y) + fsqr(z))*zb_des.x + x*(-(xb_des.x*z) + y*zb_des.y + z*zb_des.z) + w*(x*xb_des.y + z*zb_des.y - y*(xb_des.x + zb_des.z)));
	eR.z = yb_des.x - 2*(y*(x*xb_des.x + y*yb_des.x - x*yb_des.y) + w*(x*xb_des.z + y*yb_des.z)) + 2*(-(xb_des.z*y) + w*(xb_des.x + yb_des.y) + x*yb_des.z)*z - 2*yb_des.x*fsqr(z) + xb_des.y*(-1 + 2*fsqr(x) + 2*fsqr(z));
	return vscl(0.5, eR);
}

void quad_ctrl_SE3_default_params(struct quad_ctrl_SE3_params *params)
{
	struct quad_ctrl_SE3_params p = {
		.linear = {
			.kp = {.xy = 12.5f, .z = 39.0f},
			.ki = {.xy = 1.55f, .z = 1.55f},
			.kd = {.xy = 6.25f, .z = 12.5f},
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
		.omega = {
			.kp = {.xy = 18.6f, .z = 11.2f},
			.ki = {.xy = 0.00f, .z = 0.00f},
			//.kd = {.xy = 0.0f, .z = 0.0f}, // TODO
		},
		.attitude = {
			.kp = {.xy = 130.3f, .z = 111.7f},
		},

		.int_linear_bound = {.x = 2.0f, .y = 2.0f, .z = 0.4f},
		.int_omega_bound = {.x = 1.0f, .y = 1.0f, .z = 0.0f}, // TODO
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

void quad_ctrl_SE3_init(struct quad_ctrl_SE3_state *state)
{
	state->int_linear_err = vzero();
	state->int_omega_err = vzero();
}

struct quad_accel quad_ctrl_SE3(
	struct quad_ctrl_SE3_state *state,
	struct quad_ctrl_SE3_params const *param,
	struct quad_state const *s, struct quad_state const *set, float dt)
{
	struct mat33 const R = quat2rotmat(s->quat);
	struct quad_accel output;

	// -------------------- Linear part --------------------
	struct vec const pos_error = vsub(set->pos, s->pos);
	struct vec const vel_error = vsub(set->vel, s->vel);
	integrate_clamp(&state->int_linear_err, &param->int_linear_bound, &pos_error, dt);
	struct vec const target_acc = vadd4(
		set->acc,
		xy_zmul(param->linear.kp, pos_error),
		xy_zmul(param->linear.ki, state->int_linear_err),
		xy_zmul(param->linear.kd, vel_error)
	);
	struct vec const target_thrust = vadd(target_acc, mkvec(0.0f, 0.0f, GRAV));

	struct vec const z_axis = mcolumn(R, 2);
	output.linear = vdot(target_thrust, z_axis);

	// TODO: check for saturation and prioritize altitude hold, a la
	// Mike Hamer's (github.com/mikehamer/crazyflie-firmware/blob/K31-NL/src/modules/src/controller_new.c)

	// -------------------- Angular part --------------------
	// construct the desired attitude as a rotation matrix.
	// (sigularity when the thrust vector is near the horizontal plane -
	//  world yaw can only lie in one line on the xy plane.
	//  when near, becomes unstable.)
	float const target_yaw = quat2rpy(set->quat).z;
	struct vec const xc_des = mkvec(cosf(target_yaw), sinf(target_yaw), 0.0f);
	struct vec const zb_des = vnormalize(target_thrust);
	struct vec const yb_des = vnormalize(vcross(zb_des, xc_des));
	struct vec const xb_des = vcross(yb_des, zb_des);

	// compute the body frame angular velocity vector to correct the attitude
	struct mat33 Rdes = mcolumns(xb_des, yb_des, zb_des);
	struct quat const q_des = mat2quat(Rdes);
	struct quat const q_err = qposreal(qqmul(q_des, qinv(s->quat)));
	// scale to match output of vee map
	struct vec const eR = vscl(2.0f / sqrtf(2.0f), quatimagpart(q_err));

	struct vec const omega_error = vsub(set->omega, s->omega);
	integrate_clamp(&state->int_omega_err, &param->int_omega_bound, &omega_error, dt);
	// TODO add feedforward term for angular acceleration in input?
	// TODO implement derivative term for angular velocity error
	output.angular = vadd3(
		xy_zmul(param->attitude.kp, eR),
		xy_zmul(param->omega.kp, omega_error),
		xy_zmul(param->omega.ki, state->int_omega_err)
	);

	return output;
}

void quad_ctrl_attitude_default_params(struct quad_ctrl_attitude_params *params)
{
	struct quad_ctrl_attitude_params p = {
		// TODO all comments in in SE3 controller apply here too
		.omega = {
			.kp = {.xy = 18.6f, .z = 11.2f},
			.ki = {.xy = 0.00f, .z = 0.00f},
			//.kd = {.xy = 0.0f, .z = 0.0f}, // TODO
		},
		.attitude = {
			.kp = {.xy = 130.3f, .z = 111.7f},
		},

		.int_linear_bound = {.x = 2.0f, .y = 2.0f, .z = 0.4f},
		.int_omega_bound = {.x = 1.0f, .y = 1.0f, .z = 0.0f}, // TODO
	};
	*params = p;
}

void quad_ctrl_attitude_init(struct quad_ctrl_attitude_state *state)
{
	state->int_omega_err = vzero();
}

static float sign(float x) {
	return (x >= 0.0f) ? 1.0f : -1.0f;
}

struct quad_accel quad_ctrl_attitude(
	struct quad_ctrl_attitude_state *state,
	struct quad_ctrl_attitude_params const *param,
	struct quad_state const *s, struct quad_state const *set, float thrust_set, float dt)
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
	struct quat q_mix = qqmul(q_full, qinv(q_reduced));

	// fix double-covering by forcing q_mix to be a rotation about +z
	// so the angle expressed by q_mix can be compared to current z angular velocity
	assert(q_mix.x < 1e-4);
	assert(q_mix.y < 1e-4);
	if (q_mix.z < 0) {
		q_mix = qneg(q_mix);
	}
	// angle of yaw rotation
	float const alpha_mix = quat2angle(q_mix);

	// convert our PID-style params to time constant used in paper
	float const tau = 1.0 / param->attitude.kp.xy;
	float const yaw_priority = param->attitude.kp.z * tau;

	// if body yaw rate is high, we should prefer to do a yaw correction
	// in the same direction even if it's more than 180 degrees.
	// copied verbatim from ETHZ paper. would like to simplify a little --
	// do we really need the trig and all the pi factors?
	// also, unless yaw priority is quite low, the threshold ends up being huge
	float const omega_z_min = sinf(M_PI_2_F * yaw_priority) / tau;
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

	// mix feedforward with feedback
	struct vec const omega_feedback = vscl(2.0 / tau * sign(q_err.w),
		quatimagpart(q_err));
	struct quad_accel const output = {
		.angular = vadd(omega_feedback, set->omega),
		.linear = thrust_set
	};
	return output;
}


struct quad_accel quad_ctrl_attitude_rate(
	struct quad_ctrl_attitude_rate_state *state,
	struct quad_ctrl_attitude_rate_params const *param,
	struct vec s, struct vec set, float thrust, float dt)
{
	struct quad_accel output;

	struct vec const omega_error = vsub(set, s);
	integrate_clamp(&state->int_omega_err, &param->int_omega_bound, &omega_error, dt);
	// TODO add feedforward term for angular acceleration in input?
	// TODO implement derivative term for angular velocity error
	output.angular = vadd(
		xy_zmul(param->kp, omega_error),
		xy_zmul(param->ki, state->int_omega_err)
	);
	output.linear = thrust;
	return output;
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

