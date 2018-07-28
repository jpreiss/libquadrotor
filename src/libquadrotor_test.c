// Putting together tests in simplest way possible:
// including all source files into test file.
// This makes it possible to test internal static functions
// and lets us avoid dealing with object files at all.

#include "quad_control.c"
#include "quad_dynamics.c"
#include "quad_ekf.c"

#include <assert.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


bool closeeps(float a, float b, float eps)
{
	double const big = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
	double denom = big > 0.001 ? big : 1.0;
	return fabs((double)a - (double)b) / denom < (double)eps;
}

bool vcloseeps(struct vec a, struct vec b, float eps)
{
	return closeeps(a.x, b.x, eps) && closeeps(a.y, b.y, eps) && closeeps(a.z, b.z, eps);
}

bool close(float a, float b) { return closeeps(a, b, 1e-6f); }
bool vclose(struct vec a, struct vec b) { return vcloseeps(a, b, 1e-6f); }

void printvec(struct vec v)
{
	printf("(%f, %f, %f)", (double)v.x, (double)v.y, (double)v.z);
}
float randu(float low, float high)
{
	return low + (rand() / (float)RAND_MAX) * (high - low);
}
struct vec randvecbox(float low, float high)
{
	return mkvec(randu(low, high), randu(low, high), randu(low, high));
}
struct quat randquat()
{
	return qnormalize(mkquat(randu(0,1), randu(0,1), randu(0,1), randu(0,1)));
}


#define RESET "\x1b[0m"
#define RED(str) "\x1b[31m" str RESET
#define GREEN(str) "\x1b[32m" str RESET

char const *hrule = "------------------------------------------------------------------";
static char const *testname = NULL;
void test(char const *name)
{
	if (testname != NULL) {
		printf("%36s: " GREEN("passed") "\n", testname);
	}
	testname = name;
}

void test_SE3_control()
{
	struct quad_ctrl_params params;
	quad_ctrl_default_params(&params);
	struct quad_state zero;
	quad_zero_state(&zero);
	float const dt = 0.01;

	test("SE3 ctrl straight up");
	{
		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.pos = mkvec(0.0f, 0.0f, 1.0f);
		struct quad_accel acc = quad_ctrl_full(&ctrlstate, &params, &s, &set, dt);
		assert(vclose(acc.angular, vzero()));
		assert(acc.linear > params.linear.kp.z / 2.0f);
	}

	test("SE3 ctrl forward");
	{
		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.pos = mkvec(1.0f, 0.0f, 0.0f);
		struct quad_accel acc = quad_ctrl_full(&ctrlstate, &params, &s, &set, dt);
		assert(close(acc.angular.x, 0.0f));
		assert(acc.angular.y > 0.0001f);
		assert(close(acc.angular.z, 0.0f));
		assert(close(acc.linear, GRAV)); // it should be dotted with the current up vec

		quad_ctrl_init(&ctrlstate);
		s.quat = qaxisangle(mkvec(0,1,0), 0.1); // rotate a little bit towards goal
		acc = quad_ctrl_full(&ctrlstate, &params, &s, &set, dt);
		assert(close(acc.angular.x, 0.0f));
		assert(close(acc.angular.z, 0.0f));
		assert(acc.linear > GRAV + 1.0f); // now, we should want to hold alt and move
	}

	test("SE3 ctrl left");
	{
		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.pos = mkvec(0.0f, 1.0f, 0.0f);
		struct quad_accel acc = quad_ctrl_full(&ctrlstate, &params, &s, &set, dt);
		assert(acc.angular.x < -0.0001f);
		assert(close(acc.angular.y, 0.0f));
		assert(close(acc.angular.z, 0.0f));
		assert(close(acc.linear, GRAV)); // it should be dotted with the current up vec
	}

	test("SE3 ctrl yaw at hover");
	{
		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.quat = qaxisangle(mkvec(0,0,1), radians(90));
		struct quad_accel acc = quad_ctrl_full(&ctrlstate, &params, &s, &set, dt);
		// TODO: FP error should not be so big, I think
		assert(closeeps(acc.angular.z, params.attitude.kp.z, 0.1));
		assert(close(acc.angular.x, 0));
		assert(close(acc.angular.y, 0));
		assert(close(acc.linear, GRAV));
	}

	test("SE3 ctrl yaw singularity");
	{
		// this fails using the desired rot mtx from Mellinger's paper
		// but succeeds with quaternions
		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.acc = mkvec(GRAV, 0, -GRAV);
		set.quat = qaxisangle(mkvec(0,1,0), radians(89.9999));
		struct quad_accel acc = quad_ctrl_full(&ctrlstate, &params, &s, &set, dt);
		assert(close(acc.angular.x, 0));
		assert(close(acc.angular.y, params.attitude.kp.xy));
		assert(close(acc.angular.z, 0));
	}
}

void test_quaternion_control()
{
	struct quad_state zero;
	quad_zero_state(&zero);
	float const dt = 0.01;

	test("quat ctrl hover");
	{
		struct quad_ctrl_params params;
		quad_ctrl_default_params(&params);

		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;

		struct vec moment = quad_ctrl_attitude(
			&ctrlstate, &params, &s, &set, dt);
		assert(vclose(moment, vzero()));
	}

	// TODO: should be possible to make the ratio bounds
	// in roll and pitch tests bigger than 1e3. is f32 error really that bad?
	test("quat ctrl roll");
	{
		struct quad_ctrl_params params;
		quad_ctrl_default_params(&params);

		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.quat = qaxisangle(mkvec(1,0,0), M_PI_2_F);

		struct vec moment = quad_ctrl_attitude(
			&ctrlstate, &params, &s, &set, dt);
		assert(moment.x > 1);
		assert(fabs(moment.x / moment.y) > 1e3);
		assert(fabs(moment.x / moment.z) > 1e3);
	}

	float pitch_90deg_omega;
	test("quat ctrl pitch");
	{
		struct quad_ctrl_params params;
		quad_ctrl_default_params(&params);

		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.quat = qaxisangle(mkvec(0,1,0), M_PI_2_F);

		struct vec moment = quad_ctrl_attitude(
			&ctrlstate, &params, &s, &set, dt);
		assert(fabs(moment.y / moment.x) > 1e3);
		assert(moment.y > 1);
		assert(fabs(moment.y / moment.z) > 1e3);
		pitch_90deg_omega = moment.y;
	}

	test("quat ctrl yaw");
	{
		struct quad_ctrl_params params;
		quad_ctrl_default_params(&params);

		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.quat = qaxisangle(mkvec(0,0,1), M_PI_2_F);

		struct vec moment = quad_ctrl_attitude(
			&ctrlstate, &params, &s, &set, dt);
		assert(close(moment.x, 0.0f));
		assert(close(moment.y, 0.0f));
		assert(moment.z > 1e-4);
	}

	test("quat ctrl pitch flip");
	{
		struct quad_ctrl_params params;
		quad_ctrl_default_params(&params);

		// this time we ask for a near-180 degree pitch flip
		// so we verify that the commanded angular accel is greater
		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.quat = qaxisangle(mkvec(0,1,0), M_PI_F - 0.1);

		struct vec moment = quad_ctrl_attitude(
			&ctrlstate, &params, &s, &set, dt);
		assert(close(moment.x, 0.0f));
		assert(moment.y > pitch_90deg_omega);
		assert(close(moment.z, 0.0f));
	}

	test("quat ctrl yaw priority");
	{
		struct quad_ctrl_params params;
		quad_ctrl_default_params(&params);

		float const kpz_orig = params.attitude.kp.z;

		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		// positive pitch after yaw ends up pointing towards +y
		set.quat = qqmul(
			qaxisangle(mkvec(0,1,0), M_PI_2_F),
			qaxisangle(mkvec(0,0,1), M_PI_2_F));
		struct vec moment;

		// typical case - yaw priority 0.4
		params.attitude.kp.z = 0.4 * params.attitude.kp.xy;
		moment = quad_ctrl_attitude(
			&ctrlstate, &params, &s, &set, dt);
		assert(fabs(moment.x) > fabs(moment.z));

		// extreme case, yaw priority almost 0
		params.attitude.kp.z = 1.0 * params.attitude.kp.xy;
		moment = quad_ctrl_attitude(
			&ctrlstate, &params, &s, &set, dt);
		assert(close(fabs(moment.x), fabs(moment.z)));

		// extreme case, equal yaw priority
		params.attitude.kp.z = 1e-6 * params.attitude.kp.xy;
		moment = quad_ctrl_attitude(
			&ctrlstate, &params, &s, &set, dt);
		assert(fabs(moment.z) < 1e-3);

		params.attitude.kp.z = kpz_orig;
	}

	test("quat ctrl random one step");
	{
		struct quad_ctrl_params params;
		quad_ctrl_default_params(&params);

		struct quad_ctrl_state ctrlstate;

		srand(1); // deterministic
		int const TRIALS = 1000;
		for (int i = 0; i < TRIALS; ++i) {

			quad_ctrl_init(&ctrlstate);
			struct quad_state s = zero, set = zero;

			// who knows what is the probability distribution of this quat???
			s.quat = randquat();
			set.quat = randquat();
			struct vec dummy = quat2axis(s.quat); // to get quat2axis compiled in binary
			#pragma unused(dummy)
			float const before_err = quat2angle(qqmul(set.quat, qinv(s.quat)));

			struct vec moment = quad_ctrl_attitude(
				&ctrlstate, &params, &s, &set, dt);
			struct quat const next = qnormalize(
				quat_gyro_update(s.quat, moment, 1e-4));
			float const next_err = quat2angle(qqmul(set.quat, qinv(next)));

			assert(fabs(next_err) < fabs(before_err));
		}
	}

	test("quat ctrl yaw vel heuristic");
	{
		struct quad_ctrl_params params;
		quad_ctrl_default_params(&params);
		// to avoid omega feedback interfering
		params.attitude.kd.z = 0.0f;

		struct quad_ctrl_state ctrlstate;
		quad_ctrl_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.quat = qaxisangle(mkvec(0,0,1), -M_PI_2_F);

		// heuristic should not happen for low velocity
		s.omega.z = 1.0f;
		struct vec moment = quad_ctrl_attitude(
			&ctrlstate, &params, &s, &set, dt);
		assert(moment.z < 0.0f);

		// heuristic should definitely happen for high velocity
		s.omega.z = 1000.0f;
		moment = quad_ctrl_attitude(
			&ctrlstate, &params, &s, &set, dt);
		assert(moment.z > 0.0f);
	}
}

void test_powerdist()
{
	test("power dist basic");
	{
		struct quad_physical_params params;
		physical_params_crazyflie2(&params);
		float motors[4];

		// free fall sanity check
		struct quad_accel acc = { .linear = 0.0f, .angular = { 0.0f, 0.0f, 0.0f }};
		quad_power_distribute(&acc, &params, motors);
		for (int i = 1; i < 4; i++) {
			assert(close(motors[i], 0.0f));
		}

		// hover - total force should be exactly gravity compensation
		acc.linear = GRAV;
		quad_power_distribute(&acc, &params, motors);
		float const tot = motors[0] + motors[1] + motors[2] + motors[3];
		assert(close(tot, GRAV * params.mass));
		for (int i = 1; i < 4; i++) {
			assert(close(motors[i], motors[0]));
		}

		// TODO: actually check that they match the moment of inertia

		// yaw CCW
		acc.angular.z = 1.0f;
		quad_power_distribute(&acc, &params, motors);
		float const grav_comp_thrust = tot / 4.0f; // `tot` from previous test
		assert(params.motor_0_ccw);
		// motor exerts opposite reaction torque on body
		assert(motors[0] < grav_comp_thrust);
		assert(motors[1] > grav_comp_thrust);
		assert(motors[2] < grav_comp_thrust);
		assert(motors[3] > grav_comp_thrust);

		// pitch forward
		acc.angular = mkvec(0.0f, 1.0f, 0.0f);
		quad_power_distribute(&acc, &params, motors);
		assert(motors[0] < grav_comp_thrust);
		assert(motors[1] < grav_comp_thrust);
		assert(motors[2] > grav_comp_thrust);
		assert(motors[3] > grav_comp_thrust);

		// roll left
		acc.angular = mkvec(-1.0f, 0.0f, 0.0f);
		quad_power_distribute(&acc, &params, motors);
		assert(motors[0] < grav_comp_thrust);
		assert(motors[1] > grav_comp_thrust);
		assert(motors[2] > grav_comp_thrust);
		assert(motors[3] < grav_comp_thrust);
	}

	test("power dist limits");
	{
		struct quad_physical_params params;
		physical_params_crazyflie2(&params);
		float motors[4];

		// high side clipping - request too high thrust
		struct quad_accel acc = { .linear = 10.0f * GRAV, .angular = vzero() };
		quad_power_distribute(&acc, &params, motors);
		for (int i = 1; i < 4; i++) {
			assert(motors[i] == params.max_thrust);
		}

		// low side clipping - request negative thrust
		// (TODO: should this even be considered valid input?)
		acc.linear = -GRAV;
		quad_power_distribute(&acc, &params, motors);
		for (int i = 1; i < 4; i++) {
			assert(motors[i] == 0.0f);
		}

		// roll hard
		acc.linear = GRAV;
		acc.angular = mkvec(1e5f, 0, 0);
		quad_power_distribute(&acc, &params, motors);
		assert(motors[0] == params.max_thrust);
		assert(motors[1] == 0.0f);
		assert(motors[2] == 0.0f);
		assert(motors[3] == params.max_thrust);

		// pitch hard
		acc.linear = GRAV;
		acc.angular = mkvec(0, 1e5f, 0);
		quad_power_distribute(&acc, &params, motors);
		assert(motors[0] == 0.0f);
		assert(motors[1] == 0.0f);
		assert(motors[2] == params.max_thrust);
		assert(motors[3] == params.max_thrust);
	}

	test("motors to moment/thrust");
	{
		struct quad_physical_params param;
		physical_params_crazyflie2(&param);
		struct quad_accel force;
		float motors[4];

		float const grav_comp = (GRAV * param.mass / 4.0f);

		// hovering
		motors[0] = motors[1] = motors[2] = motors[3] = grav_comp;
		quad_motor_forces(&param, motors, &force);
		assert(close(force.linear, GRAV * param.mass));
		assert(vclose(force.angular, vzero()));

		// rolling
		motors[0] = motors[3] = 1.1f * grav_comp;
		motors[1] = motors[2] = 0.9f * grav_comp;
		quad_motor_forces(&param, motors, &force);
		assert(close(force.linear, GRAV * param.mass));
		assert(force.angular.x > 0.0001f); // TODO exact amount
		assert(close(force.angular.y, 0.0f));
		assert(close(force.angular.z, 0.0f));

		// pitching
		motors[0] = motors[1] = 0.9f * grav_comp;
		motors[2] = motors[3] = 1.1f * grav_comp;
		quad_motor_forces(&param, motors, &force);
		assert(close(force.linear, GRAV * param.mass));
		assert(close(force.angular.x, 0.0f));
		assert(force.angular.y > 0.0001f); // TODO exact amount
		assert(close(force.angular.z, 0.0f));

		// yawing
		motors[0] = motors[2] = 0.9f * grav_comp;
		motors[1] = motors[3] = 1.1f * grav_comp;
		quad_motor_forces(&param, motors, &force);
		assert(close(force.linear, GRAV * param.mass));
		assert(close(force.angular.x, 0.0f));
		assert(close(force.angular.y, 0.0f));
		assert(force.angular.z > 0.0001f); // TODO exact amount
	}

	test("power dist roundtrip");
	{
		struct quad_physical_params params;
		physical_params_crazyflie2(&params);
		float motors[4];
		struct quad_accel force;

		// only test small forces/moments to avoid motor clipping
		float const linear_range = 1.0f;
		float const angular_range = 1.0f;

		srand(100);
		int const TRIALS = 100;
		for (int i = 0; i < TRIALS; ++i) {
			struct quad_accel acc = {
				.linear = GRAV + randu(-linear_range, linear_range),
				.angular = randvecbox(-angular_range, angular_range)
			};
			quad_power_distribute(&acc, &params, motors);
			quad_motor_forces(&params, motors, &force);

			force.linear /= params.mass;
			force.angular = veltdiv(force.angular, params.inertia);
			assert(close(force.linear, acc.linear));
			// TODO: figure out why a lot of precision is being lost here.
			// I would hope we can get better than 1e-2f.
			assert(vcloseeps(force.angular, acc.angular, 1e-2f));
		}
	}
}

void test_dynamics()
{
	int const HZ = 100;
	float const dt = 1.0f / HZ;
	struct quad_physical_params param;
	physical_params_crazyflie2(&param);
	float const grav_comp = param.mass * GRAV;

	test("freefall");
	{
		struct quad_state now, next;
		quad_zero_state(&now);
		struct quad_accel force = { .linear = 0.0f, .angular = vzero() };
		for (int i = 0; i < HZ; ++i) {
			quad_dynamics(&param, &now, &force, dt, &next);
			now = next;
		}
		assert(now.pos.z < -0.95f * 0.5f * GRAV);
		assert(now.pos.z > -1.05f * 0.5f * GRAV);
	}

	test("hover");
	{
		struct quad_state now, next;
		quad_zero_state(&now);
		struct quad_accel force = { .linear = grav_comp, .angular = vzero() };
		for (int i = 0; i < HZ; ++i) {
			quad_dynamics(&param, &now, &force, dt, &next);
			now = next;
		}
		assert(vclose(next.pos, vzero()));
	}

	test("freefall rotations");
	{
		struct quad_state now, next;
		struct quad_accel force = { .linear = 0.0f, .angular = vzero() };
		// TODO: union-ize struct vec to allow?
		float const *inertia_arr = (float const *)&param.inertia;
		float *moment_arr = (float *)&force.angular;

		for (int axis = 0; axis < 2; ++axis) {
			quad_zero_state(&now);
			force.angular = vzero();
			moment_arr[axis] = M_PI_F * inertia_arr[axis];
			for (int i = 0; i < HZ; ++i) {
				quad_dynamics(&param, &now, &force, dt, &next);
				now = next;
			}
			struct vec rotaxis = vzero();
			((float *)&rotaxis)[axis] = 1.0f;
			struct quat qgoal = qaxisangle(rotaxis, M_PI_F / 2.0f);
			// slop must be allowed due to damping
			assert(qanglebetween(qgoal, next.quat) < radians(1.0f));
		}
	}
}

// use the controller and the dynamics in a loop.
void test_closedloop()
{
	int const HZ = 100;
	float const dt = 1.0f / HZ;
	struct quad_physical_params param;
	physical_params_crazyflie2(&param);
	struct quad_ctrl_params ctrl_params;
	quad_ctrl_default_params(&ctrl_params);

	test("closed loop attitude correction");
	{
		srand(100);
		struct quad_ctrl_state ctrl_state;
		struct quad_state now, next, goal;
		quad_zero_state(&goal);

		int const TRIALS = 100;
		for (int i = 0; i < TRIALS; ++i) {
			quad_zero_state(&now);
			quad_ctrl_init(&ctrl_state);

			// no, this is not uniformly distributed on the sphere
			float const max_angle = radians(60.0);
			struct vec const axis = vnormalize(randvecbox(-1.0, 1.0));
			float const angle = randu(-max_angle, max_angle);

			now.quat = qaxisangle(axis, angle);
			now.omega = randvecbox(-0.5, 0.5);

			for (int t = 0; t < 300; ++t) {
				struct quad_accel acc = quad_ctrl_full(
					&ctrl_state, &ctrl_params, &now, &goal, dt);
				// assume the control is perfectly realized.
				// TODO: motor clipping & possibly motor inertia simulation.
				if (acc.linear < 0.0f) acc.linear = 0.0f;
				acc.linear *= param.mass;
				acc.angular = veltmul(acc.angular, param.inertia);
				quad_dynamics(&param, &now, &acc, dt, &next);
				now = next;
			}

			float const pos_err = vmag(vsub(now.pos, goal.pos));
			assert(pos_err < 0.05f);
			float const angle_err = degrees(qanglebetween(now.quat, goal.quat));
			assert(angle_err < 5.0f);
		}
	}

	test("closed loop position correction");
	{
		srand(100);
		struct quad_ctrl_state ctrl_state;
		struct quad_state now, next, goal;
		quad_zero_state(&goal);

		int const TRIALS = 100;
		for (int i = 0; i < TRIALS; ++i) {
			quad_zero_state(&now);
			quad_ctrl_init(&ctrl_state);

			struct vec const axis = vnormalize(randvecbox(-1.0, 1.0));
			float const angle = randu(-0.9, 0.9);
			now.quat = qaxisangle(axis, angle);
			now.omega = randvecbox(-0.5, 0.5);

			now.pos = randvecbox(-0.5, 0.5);
			now.vel = randvecbox(-0.5, 0.5);

			for (int t = 0; t < 1000; ++t) {
				struct quad_accel acc = quad_ctrl_full(
					&ctrl_state, &ctrl_params, &now, &goal, dt);
				// assume the control is perfectly realized.
				// TODO: motor clipping & possibly motor inertia simulation.
				if (acc.linear < 0.0f) acc.linear = 0.0f;
				acc.linear *= param.mass;
				acc.angular = veltmul(acc.angular, param.inertia);
				quad_dynamics(&param, &now, &acc, dt, &next);
				now = next;
			}

			float const pos_err = vmag(vsub(now.pos, goal.pos));
			assert(pos_err < 0.03f);
			assert(vmag(now.vel) < 0.01f);
			float const angle_err = degrees(qanglebetween(now.quat, goal.quat));
			assert(angle_err < 5.0f);
		}
	}
}

static struct quad_ekf ekf_repeat(struct quad_ekf const *init,
	struct vec gyro, struct vec acc, float dt, int steps)
{
	struct quad_ekf ekfs[2] = { *init, *init };
	for (int i = 0; i < steps; ++i) {
		int back = i & 0x1;
		quad_ekf_imu(&ekfs[back], &ekfs[!back], acc, gyro, dt);
	}
	return ekfs[steps & 0x1];
}

void test_ekf()
{
	struct quad_ekf ekf_zero;
	quad_ekf_init(&ekf_zero, vzero(), vzero(), qeye());

	test("EKF hover");
	{
		float dt = 0.01;
		struct vec gyro = vzero();
		struct vec accelerometer = mkvec(0, 0, GRAV);
		struct quad_ekf final = ekf_repeat(&ekf_zero, gyro, accelerometer, dt, 100);
		assert(vclose(final.state.pos, vzero()));
		assert(vclose(final.state.vel, vzero()));
		assert(vclose(final.state.acc, vzero()));
		assert(qanglebetween(final.state.quat, qeye()) < 1e-6f);
	}

	test("EKF integral translate");
	{
		float dt = 0.001;
		struct vec gyro = vzero();

		srand(0);
		int const TRIALS = 100;
		for (int i = 0; i < TRIALS; ++i) {
			struct vec dir = randvecbox(-1.0f, 1.0f);
			struct vec accel = vadd(dir, mkvec(0, 0, GRAV));

			int const steps = 100;
			struct quad_ekf moving = ekf_repeat(&ekf_zero, gyro, accel, dt, steps);
			float t = steps * dt;
			// tolerance must be fairly loose due to discretization error
			// (TODO: could we use Verlet integration in the EKF?)
			assert(vcloseeps(moving.state.vel, vscl(t, dir), 1e-2));
			assert(vcloseeps(moving.state.pos, vscl(0.5 * fsqr(t), dir), 1e-2));
		}
	}

	test("EKF integral axis aligned rotate");
	{
		float dt = 0.001;
		struct vec grav_world = mkvec(0, 0, GRAV);

		srand(0);
		int const TRIALS = 100;
		int const steps = 100;

		// the tolerances in this test need to be pretty loose
		// because the ekf gyro update sacrifices accuracy
		// to avoid trig functions.

		for (int trial = 0; trial < TRIALS; ++trial) {
			int const axis = rand() % 3;
			float const max_speed = 10.0f;
			assert(max_speed * steps * dt < M_PI_F);
			float const speed = randu(-max_speed, max_speed);
			if (fabs(speed) < 1e-2) continue;
			float g[3] = { 0.0f, 0.0f, 0.0f };
			g[axis] = speed;
			struct vec gyro = vloadf(g);

			struct quad_ekf ekfs[2] = { ekf_zero, ekf_zero };
			for (int i = 0; i < steps; ++i) {
				int back = i & 0x1;
				struct quat q_world2body = qinv(ekfs[back].state.quat);
				struct vec acc = qvrot(q_world2body, grav_world);
				quad_ekf_imu(&ekfs[back], &ekfs[!back], acc, gyro, dt);
			}
			struct quad_ekf ekf_end = ekfs[steps & 0x1];
			struct quat q = qposreal(ekf_end.state.quat);
			struct vec q_axis = quat2axis(q);
			float q_angle = quat2angle(q);
			assert(q_angle > 0.0f);
			assert(vclose(vnormalize(q_axis), vnormalize(gyro)));
			assert(closeeps(q_angle, fabs(steps * dt * speed), 1e-2));
		}
	}

	test("EKF variance grows");
	{
		float dt = 0.1;
		struct vec gyro = vzero();
		struct vec accelerometer = mkvec(0, 0, GRAV);
		struct quad_ekf final = ekf_repeat(&ekf_zero, gyro, accelerometer, dt, 1000);
		for (int i = 0; i < 6; ++i) {
			assert(final.P[i][i] > 1.1f * ekf_zero.P[i][i]);
		}
		// gyro variance is much smaller
		for (int i = 6; i < QUAD_EKF_N; ++i) {
			assert(final.P[i][i] > 1.001f * ekf_zero.P[i][i]);
		}
	}

	test("EKF update reduces variance");
	{
		float dt = 0.1;
		struct vec gyro = vzero();
		struct vec accelerometer = mkvec(0, 0, GRAV);
		struct quad_ekf propagated = ekf_repeat(
			&ekf_zero, gyro, accelerometer, dt, 100);
		struct quad_ekf updated;
		quad_ekf_fullstate(&propagated, &updated, vzero(), vzero(), qeye());
		for (int i = 0; i < QUAD_EKF_N; ++i) {
			assert(updated.P[i][i] < 0.1f * propagated.P[i][i]);
		}
	}
}

static void sigabort(int unused)
{
	puts(hrule);
	printf(RED("libquadrotor: test \"%s\" failed.\n"), testname);
	puts("");
}

int main()
{
	signal(SIGABRT, sigabort);

	puts("");
	puts("testing libquadrotor...");
	puts(hrule);

	test_SE3_control();
	test_quaternion_control();
	test_powerdist();
	test_dynamics();
	test_closedloop();
	test_ekf();

	test("dummy");
	puts(hrule);
	puts(GREEN("libquadrotor: all tests passed."));
	puts("");
	return 0;
}

