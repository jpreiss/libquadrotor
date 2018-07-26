#pragma once

#include "quad_common_types.h"

// Indirect (a.k.a. error-state) EKF.
//
// Rather than tracking the navigation states directly, this filter contains a "dumb"
// dead-reckoning IMU integrator and uses EKF to estimate the error of this integrator.
// When a measurement arrives, the filter uses it to correct the integrator state
// and updates the error estimate correspondingly.
//
// The advantage of this formulation is that it permits us to track the orientation error
// with RPY angles, which eliminates the problem of the quaternion having 4 dimensions
// but only 3 degrees of freedom. Also, it keeps the EKF quantities close to 0,
// which is a good thing [citation needed].

#define QUAD_EKF_N 9           // state dimension (pos, vel, quat err)
#define QUAD_EKF_M 9           // measurement dimension
#define QUAD_EKF_DISTURBANCE 6 // control (IMU) noise dimension

struct quad_ekf
{
	struct quad_state state;

	// TODO symmetric matrix storage optimization?
	float P[QUAD_EKF_N][QUAD_EKF_N]; // error state covariance
};

void quad_ekf_init(struct quad_ekf *ekf, struct vec pos, struct vec vel, struct quat quat);

void quad_ekf_imu(struct quad_ekf const *ekf_prev, struct quad_ekf *ekf,
	struct vec acc, struct vec gyro, float dt);

void quad_ekf_fullpose(struct quad_ekf const *ekf_prev, struct quad_ekf *ekf,
	struct vec pos, struct vec vel, struct quat quat);
