#include "libquadrotorcontrol.h"
#include <math.h>

static inline void integrate_clamp(vec_t *integral,
	vec_t const *bound, vec_t const *value, float dt)
{
	*integral = vclampabs(vadd(*integral, vscl(dt, *value)), *bound);
}

static inline vec_t xy_zmul(struct xy_z xy_z, struct vec v)
{
	return mkvec(xy_z.xy * v.x, xy_z.xy * v.y, xy_z.z * v.z);
}

struct vec vee_map_rot_error(struct mat33 const *R, struct mat33 const *Rdes)
{
	struct mat33 R_transpose = mtranspose(*R);
	struct mat33 Rdes_transpose = mtranspose(*Rdes);
	struct mat33 eRM = msub(mmult(Rdes_transpose, *R), mmult(R_transpose, *Rdes));
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



static float const GRAV = 9.81;

void ctrl_SE3_default_params(struct ctrl_SE3_params *params)
{
	struct ctrl_SE3_params p = {
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

void physical_params_crazyflie2(struct quadrotor_physical_params *params)
{
	// from "System Identification of the Crazyflie 2.0 Nano Quadrocopter".
	// J. Foerster, M. Hamer, R. D'Andrea. Bachelor Thesis, ETH Zurich, 2015.
	struct quadrotor_physical_params p = {
		.mass = 0.027f,
		.arm_length = 0.046f,
		.inertia = { 1.66e-5f, 1.66e-5f, 2.92e-5f },
		.thrust_to_torque = 0.006,
		.motor_layout = 'x',
		.motor_0_ccw = true, // TODO
	};
	*params = p;
}

void init_ctrl_SE3(struct ctrl_SE3_state *state)
{
	state->int_linear_err = vzero();
	state->int_omega_err = vzero();
}

struct accel ctrl_SE3(
	struct ctrl_SE3_state *state,
	struct ctrl_SE3_params const *param,
	struct state const *s, struct state const *set, float dt)
{
	struct mat33 const R = quat2rotmat(s->quat);
	struct accel output;

	// -------------------- Linear part --------------------
	vec_t const pos_error = vsub(set->pos, s->pos);
	vec_t const vel_error = vsub(set->vel, s->vel);
	integrate_clamp(&state->int_linear_err, &param->int_linear_bound, &pos_error, dt);
	vec_t const target_acc = vadd4(
		set->acc,
		xy_zmul(param->linear.kp, pos_error),
		xy_zmul(param->linear.ki, state->int_linear_err),
		xy_zmul(param->linear.kd, vel_error)
	);
	vec_t const target_thrust = vadd(target_acc, mkvec(0.0f, 0.0f, GRAV));

	struct vec const z_axis = mcolumn(R, 2);
	output.linear = vdot(target_thrust, z_axis);
	// TODO: check for saturation and prioritize altitude hold, a la
	// Mike Hamer's (github.com/mikehamer/crazyflie-firmware/blob/K31-NL/src/modules/src/controller_new.c)

	// -------------------- Angular part --------------------
	// construct the desired attitude as a rotation matrix.
	// (sigularity when the thrust vector is near the horizontal plane!)
	float const target_yaw = quat2rpy(set->quat).z;
	struct vec const xc_des = mkvec(cosf(target_yaw), sinf(target_yaw), 0.0f);
	struct vec const zb_des = vnormalize(target_thrust);
	struct vec const yb_des = vnormalize(vcross(zb_des, xc_des));
	struct vec const xb_des = vcross(yb_des, zb_des);

	// compute the body frame angular velocity vector to correct the attitude
	// (note: unlike other error vars, this one is (state - setpoint))
	struct mat33 Rdes = mcolumns(xb_des, yb_des, zb_des);
	struct vec eR = vee_map_rot_error_quat2mat(&s->quat, &Rdes);
	eR = vscl(-1.0, eR);

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

struct accel ctrl_attitude_rate(
	struct ctrl_attitude_rate_state *state,
	struct ctrl_attitude_rate_params const *param,
	struct vec s, struct vec set, float thrust, float dt)
{
	struct accel output;

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

void power_distribute_quad(
	struct accel const *acc,
	struct quadrotor_physical_params const *params,
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
}

#ifdef TEST
#include <assert.h>
#include <signal.h>
#include <stdio.h>
bool close(float a, float b)
{
	double const big = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
	double denom = big > 0.001 ? big : 1.0;
	return fabs((double)a - (double)b) / denom < 0.000001;
}
bool vclose(struct vec a, struct vec b)
{
	return close(a.x, b.x) && close(a.y, b.y) && close(a.z, b.z);
	/*
	if (!(close(a.x, b.x) && close(a.y, b.y) && close(a.z, b.z))) {
		puts("a =");
		printvec(a);
		puts("\nb =");
		printvec(b);
		puts("");
		return false;
	}
	return true;
	*/
}
void printvec(struct vec v)
{
	printf("(%f, %f, %f)", (double)v.x, (double)v.y, (double)v.z);
}

#define ANSI_RED     "\x1b[31m"
#define ANSI_GREEN   "\x1b[32m"
#define ANSI_RESET   "\x1b[0m"


static void sigabort(int unused)
{
	puts(ANSI_RED "libquadrotorcontrol: some tests failed." ANSI_RESET);
}

int main()
{
	signal(SIGABRT, sigabort);

	struct ctrl_SE3_params params;
	ctrl_SE3_default_params(&params);
	struct state zero;
	zero_state(&zero);
	float const dt = 0.01;

	// go straight up
	{
		struct ctrl_SE3_state ctrlstate;
		init_ctrl_SE3(&ctrlstate);
		struct state s = zero, set = zero;
		set.pos = mkvec(0.0f, 0.0f, 1.0f);
		struct accel acc = ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(vclose(acc.angular, vzero()));
		assert(acc.linear > params.linear.kp.z / 2.0f);
	}

	// go forward
	{
		struct ctrl_SE3_state ctrlstate;
		init_ctrl_SE3(&ctrlstate);
		struct state s = zero, set = zero;
		set.pos = mkvec(1.0f, 0.0f, 0.0f);
		struct accel acc = ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(close(acc.angular.x, 0.0f));
		assert(acc.angular.y > 0.0001f);
		assert(close(acc.angular.z, 0.0f));
		assert(close(acc.linear, GRAV)); // it should be dotted with the current up vec

		init_ctrl_SE3(&ctrlstate);
		s.quat = qaxisangle(mkvec(0,1,0), 0.1); // rotate a little bit towards goal
		acc = ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(close(acc.angular.x, 0.0f));
		assert(close(acc.angular.z, 0.0f));
		assert(acc.linear > GRAV + 1.0f); // now, we should want to hold alt and move
	}

	// go left
	{
		struct ctrl_SE3_state ctrlstate;
		init_ctrl_SE3(&ctrlstate);
		struct state s = zero, set = zero;
		set.pos = mkvec(0.0f, 1.0f, 0.0f);
		struct accel acc = ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(acc.angular.x < -0.0001f);
		assert(close(acc.angular.y, 0.0f));
		assert(close(acc.angular.z, 0.0f));
		assert(close(acc.linear, GRAV)); // it should be dotted with the current up vec
	}

	// yaw at hover
	{
		struct ctrl_SE3_state ctrlstate;
		init_ctrl_SE3(&ctrlstate);
		struct state s = zero, set = zero;
		set.quat = qaxisangle(mkvec(0,0,1), radians(90));
		struct accel acc = ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(close(acc.angular.z, params.attitude.kp.z));
		assert(close(acc.angular.x, 0));
		assert(close(acc.angular.y, 0));
		assert(close(acc.linear, GRAV));
	}

	// check vee map
	{
		struct mat33 R = meye();

		// roll
		struct mat33 Rdes = mcolumns(mkvec(1,0,0), mkvec(0,0,1), mkvec(0,-1,0));
		struct vec eR = vee_map_rot_error(&R, &Rdes);
		assert(vclose(eR, mkvec(-1, 0, 0)));

		// pitch
		Rdes = mcolumns(mkvec(0,0,-1), mkvec(0,1,0), mkvec(1,0,0));
		eR = vee_map_rot_error(&R, &Rdes);
		assert(vclose(eR, mkvec(0, -1, 0)));

		// yaw
		Rdes = mcolumns(mkvec(0,1,0), mkvec(-1,0,0), mkvec(0,0,1));
		eR = vee_map_rot_error(&R, &Rdes);
		assert(vclose(eR, mkvec(0, 0, -1)));
	}

	// check vee map fast
	{
		//struct mat33 R = meye();
		struct quat q = qeye();

		// roll
		struct mat33 Rdes = mcolumns(mkvec(1,0,0), mkvec(0,0,1), mkvec(0,-1,0));
		struct vec eR = vee_map_rot_error_quat2mat(&q, &Rdes);
		assert(vclose(eR, mkvec(-1, 0, 0)));

		// pitch
		Rdes = mcolumns(mkvec(0,0,-1), mkvec(0,1,0), mkvec(1,0,0));
		eR = vee_map_rot_error_quat2mat(&q, &Rdes);
		assert(vclose(eR, mkvec(0, -1, 0)));

		// yaw
		Rdes = mcolumns(mkvec(0,1,0), mkvec(-1,0,0), mkvec(0,0,1));
		eR = vee_map_rot_error_quat2mat(&q, &Rdes);
		assert(vclose(eR, mkvec(0, 0, -1)));
	}

	// check power distribution
	{
		struct quadrotor_physical_params params;
		physical_params_crazyflie2(&params);
		float motors[4];

		// free fall sanity check
		struct accel acc = { .linear = 0.0f, .angular = { 0.0f, 0.0f, 0.0f }};
		power_distribute_quad(&acc, &params, motors);
		for (int i = 1; i < 4; i++) {
			assert(close(motors[i], 0.0f));
		}

		// hover - total force should be exactly gravity compensation
		acc.linear = GRAV;
		power_distribute_quad(&acc, &params, motors);
		float const tot = motors[0] + motors[1] + motors[2] + motors[3];
		assert(close(tot, GRAV * params.mass));
		for (int i = 1; i < 4; i++) {
			assert(close(motors[i], motors[0]));
		}

		// TODO: actually check that they match the moment of inertia

		// yaw CCW
		acc.angular.z = 1.0f;
		power_distribute_quad(&acc, &params, motors);
		float const grav_comp_thrust = tot / 4.0f; // `tot` from previous test
		assert(params.motor_0_ccw);
		// motor exerts opposite reaction torque on body
		assert(motors[0] < grav_comp_thrust);
		assert(motors[1] > grav_comp_thrust);
		assert(motors[2] < grav_comp_thrust);
		assert(motors[3] > grav_comp_thrust);

		// pitch forward
		acc.angular = mkvec(0.0f, 1.0f, 0.0f);
		power_distribute_quad(&acc, &params, motors);
		assert(motors[0] < grav_comp_thrust);
		assert(motors[1] < grav_comp_thrust);
		assert(motors[2] > grav_comp_thrust);
		assert(motors[3] > grav_comp_thrust);

		// roll left
		acc.angular = mkvec(-1.0f, 0.0f, 0.0f);
		power_distribute_quad(&acc, &params, motors);
		assert(motors[0] < grav_comp_thrust);
		assert(motors[1] > grav_comp_thrust);
		assert(motors[2] > grav_comp_thrust);
		assert(motors[3] < grav_comp_thrust);
	}

	puts(ANSI_GREEN "libquadrotorcontrol: all tests passed." ANSI_RESET);

}
#endif // TEST
