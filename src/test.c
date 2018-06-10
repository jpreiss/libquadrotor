// Putting together tests in simplest way possible:
// including all source files into test file.
// This makes it possible to test internal static functions
// and lets us avoid dealing with object files at all.

#include "control.c"
#include "dynamics.c"

#include <assert.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>


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

#define RESET "\x1b[0m"
#define RED(str) "\x1b[31m" str RESET
#define GREEN(str) "\x1b[32m" str RESET

char const *hrule = "------------------------------------------------------------";
static char const *testname = NULL;
void test(char const *name)
{
	if (testname != NULL) {
		printf("%20s: " GREEN("passed") "\n", testname);
	}
	testname = name;
}

void test_control()
{
	struct ctrl_SE3_params params;
	ctrl_SE3_default_params(&params);
	struct quad_state zero;
	zero_state(&zero);
	float const dt = 0.01;

	test("go straight up");
	{
		struct ctrl_SE3_state ctrlstate;
		init_ctrl_SE3(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.pos = mkvec(0.0f, 0.0f, 1.0f);
		struct accel acc = ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(vclose(acc.angular, vzero()));
		assert(acc.linear > params.linear.kp.z / 2.0f);
	}

	test("go forward");
	{
		struct ctrl_SE3_state ctrlstate;
		init_ctrl_SE3(&ctrlstate);
		struct quad_state s = zero, set = zero;
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

	test("go left");
	{
		struct ctrl_SE3_state ctrlstate;
		init_ctrl_SE3(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.pos = mkvec(0.0f, 1.0f, 0.0f);
		struct accel acc = ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
		assert(acc.angular.x < -0.0001f);
		assert(close(acc.angular.y, 0.0f));
		assert(close(acc.angular.z, 0.0f));
		assert(close(acc.linear, GRAV)); // it should be dotted with the current up vec
	}

	test("yaw at hover");
	{
		struct ctrl_SE3_state ctrlstate;
		init_ctrl_SE3(&ctrlstate);
		struct quad_state s = zero, set = zero;
		set.quat = qaxisangle(mkvec(0,0,1), radians(90));
		struct accel acc = ctrl_SE3(&ctrlstate, &params, &s, &set, dt);
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

	test("power distribution");
	{
		struct quad_physical_params params;
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
		struct quad_state now;
		struct quad_state next;
		zero_state(&now);
		struct accel force = { .linear = 0.0f, .angular = vzero() };
		for (int i = 0; i < HZ; ++i) {
			quad_dynamics(&param, &now, &force, dt, &next);
			now = next;
		}
		assert(now.pos.z < -0.95f * 0.5f * GRAV);
		assert(now.pos.z > -1.05f * 0.5f * GRAV);
	}

	test("hover");
	{
		struct quad_state now;
		struct quad_state next;
		zero_state(&now);
		struct accel force = { .linear = grav_comp, .angular = vzero() };
		for (int i = 0; i < HZ; ++i) {
			quad_dynamics(&param, &now, &force, dt, &next);
			now = next;
		}
		assert(vclose(next.pos, vzero()));
	}

	test("freefall rotations");
	{
		struct quad_state now;
		struct quad_state next;
		struct accel force = { .linear = 0.0f, .angular = vzero() };
		// TODO: union-ize struct vec to allow?
		float const *inertia_arr = (float const *)&param.inertia;
		float *moment_arr = (float *)&force.angular;

		for (int axis = 0; axis < 2; ++axis) {
			zero_state(&now);
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
	test_dynamics();

	test("dummy");
	puts(hrule);
	puts(GREEN("libquadrotorcontrol: all tests passed."));
	puts("");
	return 0;
}
