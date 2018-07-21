// Putting together tests in simplest way possible:
// including all source files into test file.
// This makes it possible to test internal static functions
// and lets us avoid dealing with object files at all.

#include "quad_control.c"
#include "quad_dynamics.c"

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


#define RESET "\x1b[0m"
#define RED(str) "\x1b[31m" str RESET
#define GREEN(str) "\x1b[32m" str RESET

char const *hrule = "------------------------------------------------------------";
static char const *testname = NULL;
void test(char const *name)
{
	if (testname != NULL) {
		printf("%30s: " GREEN("passed") "\n", testname);
	}
	testname = name;
}

void test_control()
{
	struct quad_ctrl_SE3_params params;
	quad_ctrl_SE3_default_params(&params);
	struct quad_state zero;
	quad_zero_state(&zero);
	float const dt = 0.01;

	test("go straight up");
	{
		struct quad_ctrl_SE3_state ctrlstate;
		quad_ctrl_SE3_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.pos = mkvec(0.0f, 0.0f, 1.0f);
		struct quad_accel acc = quad_ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(vclose(acc.angular, vzero()));
		assert(acc.linear > params.linear.kp.z / 2.0f);
	}

	test("go forward");
	{
		struct quad_ctrl_SE3_state ctrlstate;
		quad_ctrl_SE3_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.pos = mkvec(1.0f, 0.0f, 0.0f);
		struct quad_accel acc = quad_ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(close(acc.angular.x, 0.0f));
		assert(acc.angular.y > 0.0001f);
		assert(close(acc.angular.z, 0.0f));
		assert(close(acc.linear, GRAV)); // it should be dotted with the current up vec

		quad_ctrl_SE3_init(&ctrlstate);
		s.quat = qaxisangle(mkvec(0,1,0), 0.1); // rotate a little bit towards goal
		acc = quad_ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(close(acc.angular.x, 0.0f));
		assert(close(acc.angular.z, 0.0f));
		assert(acc.linear > GRAV + 1.0f); // now, we should want to hold alt and move
	}

	test("go left");
	{
		struct quad_ctrl_SE3_state ctrlstate;
		quad_ctrl_SE3_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.pos = mkvec(0.0f, 1.0f, 0.0f);
		struct quad_accel acc = quad_ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(acc.angular.x < -0.0001f);
		assert(close(acc.angular.y, 0.0f));
		assert(close(acc.angular.z, 0.0f));
		assert(close(acc.linear, GRAV)); // it should be dotted with the current up vec
	}

	test("yaw at hover");
	{
		struct quad_ctrl_SE3_state ctrlstate;
		quad_ctrl_SE3_init(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.quat = qaxisangle(mkvec(0,0,1), radians(90));
		struct quad_accel acc = quad_ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(close(acc.angular.z, params.attitude.kp.z));
		assert(close(acc.angular.x, 0));
		assert(close(acc.angular.y, 0));
		assert(close(acc.linear, GRAV));
	}

	test("check vee map");
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

	test("check vee map fast");
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
	struct quad_ctrl_SE3_params ctrl_params;
	quad_ctrl_SE3_default_params(&ctrl_params);

	test("hover attitude correction");
	{
		srand(100);
		struct quad_ctrl_SE3_state ctrl_state;
		struct quad_state now, next, goal;
		quad_zero_state(&goal);

		int const TRIALS = 100;
		for (int i = 0; i < TRIALS; ++i) {
			quad_zero_state(&now);
			quad_ctrl_SE3_init(&ctrl_state);

			// no, this is not uniformly distributed on the sphere
			float const max_angle = radians(60.0);
			struct vec const axis = vnormalize(randvecbox(-1.0, 1.0));
			float const angle = randu(-max_angle, max_angle);

			now.quat = qaxisangle(axis, angle);
			now.omega = randvecbox(-0.5, 0.5);

			for (int t = 0; t < 200; ++t) {
				struct quad_accel acc = quad_ctrl_SE3(
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

	test("hover position correction");
	{
		srand(100);
		struct quad_ctrl_SE3_state ctrl_state;
		struct quad_state now, next, goal;
		quad_zero_state(&goal);

		int const TRIALS = 100;
		for (int i = 0; i < TRIALS; ++i) {
			quad_zero_state(&now);
			quad_ctrl_SE3_init(&ctrl_state);

			struct vec const axis = vnormalize(randvecbox(-1.0, 1.0));
			float const angle = randu(-0.9, 0.9);
			now.quat = qaxisangle(axis, angle);
			now.omega = randvecbox(-0.5, 0.5);

			now.pos = randvecbox(-0.5, 0.5);
			now.vel = randvecbox(-0.5, 0.5);

			for (int t = 0; t < 500; ++t) {
				struct quad_accel acc = quad_ctrl_SE3(
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

	test_control();
	test_powerdist();
	test_dynamics();
	test_closedloop();

	test("dummy");
	puts(hrule);
	puts(GREEN("libquadrotorcontrol: all tests passed."));
	puts("");
	return 0;
}

