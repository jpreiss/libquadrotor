#pragma once

#include <math.h>
#include <stdbool.h> // bool

#include "cholsl.h"
#include "quad_ekf.h"


// --------------------------------------------- //
// in-place hacks to reduce math3d.h stack usage //
// --------------------------------------------- //


// sadly, GCC does not optimize the pass-by-value well in mmult
static inline void mmultp(struct mat33 const *a, struct mat33 const *b, struct mat33 *ab) {
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			float accum = 0;
			for (int k = 0; k < 3; ++k) {
				accum += a->m[i][k] * b->m[k][j];
			}
			ab->m[i][j] = accum;
		}
	}
}

static inline void maddridgei(struct mat33 *a, float d) {
	a->m[0][0] += d;
	a->m[1][1] += d;
	a->m[2][2] += d;
}

// this operation is used a lot in the EKF and gcc was not optimizing it well enough :(
static inline struct mat33 aXplusbI(float a, struct mat33 const *X, float b)
{
	struct mat33 m;

	m.m[0][0] = a * X->m[0][0] + b;
	m.m[0][1] = a * X->m[0][1];
	m.m[0][2] = a * X->m[0][2];
	
	m.m[1][0] = a * X->m[1][0];
	m.m[1][1] = a * X->m[1][1] + b;
	m.m[1][2] = a * X->m[1][2];

	m.m[2][0] = a * X->m[2][0];
	m.m[2][1] = a * X->m[2][1];
	m.m[2][2] = a * X->m[2][2] + b;

	return m;
}


// --------------------------------------------- //
//             big matrix operations             //
// --------------------------------------------- //


#define AS_1D(x) (&(x)[0][0])

#define ZEROARR(a) \
do { \
	unsigned char *ac = (unsigned char *)a; \
	unsigned char *end = ac + sizeof(a); \
	for (; ac != end; ++ac) { \
		*ac = 0; \
	} \
} while(false);

#define COPYMAT(dst, src) \
do { \
	float const *from = AS_1D(src); \
	float const *end = from + sizeof(src) / sizeof(float); \
	float *to = AS_1D(dst); \
	for (; from != end; ++from) { \
		*to = *from; \
		++to; \
	} \
} while(false);

#define IND(A, m, n, r, c) ((A)[(r) * (n) + (c)])
#define INDT(A, m, n, r, c) IND(A, n, m, c, r)

#define SGEMM_LOOP(IA, IB) \
for (int i = 0; i < m; ++i) { \
	for (int j = 0; j < n; ++j) { \
		float accum = 0; \
		for (int w = 0; w < k; ++w) { \
			accum += IA(a, m, k, i, w) * IB(b, k, n, w, j); \
		} \
		IND(c, m, n, i, j) = beta * IND(c, m, n, i, j) + alpha * accum; \
	} \
}

static inline void sgemm(char atrans, char btrans, int m, int n, int k, float alpha, float const *a, float const *b, float beta, float *c)
{
	if (atrans == 'n' && btrans == 'n') {
		SGEMM_LOOP(IND, IND);
	}
	if (atrans == 'n' && btrans == 't') {
		SGEMM_LOOP(IND, INDT);
	}
	if (atrans == 't' && btrans == 'n') {
		SGEMM_LOOP(INDT, IND);
	}
	if (atrans == 't' && btrans == 't') {
		SGEMM_LOOP(INDT, INDT);
	}
}

#define SGEMM2D(at, bt, m, n, k, alpha, a, b, beta, c) sgemm(at, bt, m, n, k, alpha, AS_1D(a), AS_1D(b), beta, AS_1D(c))

static inline void zeromat(float *a, int m, int n)
{
	for (int i = 0; i < m * n; ++i) {
		a[i] = 0;
	}
}

static inline void eyeN(float *a, int n)
{
	zeromat(a, n, n);
	for (int i = 0; i < n; ++i) {
		a[n * i + i] = 1.0f;
	}
}

// measured constants
static float ext_var_xy = 1.5e-7;
static float ext_var_vel = 2e-4;
static float ext_var_q = 4.5e-3;

static float gyro_var_xyz = 0.2e-4;

// the accelerometer variance in z was quite a bit higher
// but to keep the code simple for now we just average them
static float acc_var_xyz = 2.4e-3;

// ------ utility functions for manipulating blocks of the EKF matrices ------

static void set_K_block33(float m[EKF_N][EKF_N], int row, int col, struct mat33 const *block)
{
	float *blockptr = &m[row][col];
	set_block33_rowmaj(blockptr, EKF_N, block);
}

static void mult_K_block33(float m[EKF_N][EKF_N], int row, int col, struct mat33 const *a, struct mat33 const *b)
{
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			float accum = 0;
			for (int k = 0; k < 3; ++k) {
				accum += a->m[i][k] * b->m[k][j];
			}
			m[row + i][col + j] = accum;
		}
	}
}

// static void set_H_block33(float h[EKF_M][EKF_N], int row, int col, struct mat33 const *block)
// {
// 	float *blockptr = &h[row][col];
// 	set_block33_rowmaj(blockptr, EKF_N, block);
// }

static void set_G_block33(float G[EKF_N][EKF_DISTURBANCE], int row, int col, struct mat33 const *block)
{
	float *blockptr = &G[row][col];
	set_block33_rowmaj(blockptr, EKF_DISTURBANCE, block);
}

// -------------------------- EKF implementation -----------------------------

// initialize the EKF struct
void ekf_init(struct ekf *ekf, float const pos[3], float const vel[3], float const quat[4])
{
	ekf->pos = vloadf(pos);
	ekf->vel = vloadf(vel);
	ekf->quat = qloadf(quat);
	//memcpy(ekf->P, ekf_cov_init, sizeof(ekf_cov_init));
	eyeN(AS_1D(ekf->P), EKF_N);
	//initUsecTimer();
}

void dynamic_matrix(struct quat const q, struct vec const omega, struct vec const acc, float const dt, float F[EKF_N][EKF_N])
{
	float const dt_p2_2 = dt * dt * 0.5f;
	float const dt_p3_6 = dt_p2_2 * dt / 3.0f;
	//float const dt_p4_24 = dt_p3_6 * dt * 0.25;
	//float const dt_p5_120 = dt_p4_24 * dt * 0.2;

	struct mat33 C_eq = quat2rotmat(q);
	struct mat33 w_sk = mcrossmat(omega);
	struct mat33 a_sk = mcrossmat(acc);
	// TEMP DEBUG
	//struct vec acc_nograv = vsub(acc, qvrot(q, mkvec(0,0,GRAV)));
	//struct mat33 const a_sk = mcrossmat(acc_nograv);

	// TODO S.Weiss doesn't subtract gravity from the accelerations here, is that right?
	struct mat33 Ca3;
	mmultp(&C_eq, &a_sk, &Ca3);

	struct mat33 A = aXplusbI(dt_p3_6, &w_sk, -dt_p2_2); // position by quaternion
	struct mat33 E = aXplusbI(-dt, &w_sk, 1.0f); // quat by quat
	struct mat33 FF = aXplusbI(dt_p2_2, &w_sk, -dt); // quat by gyro bias

	eyeN(AS_1D(F), EKF_N);

	//set_K_block33(F, 0, 3, meyescl(dt));
	F[0][3] = F[1][4] = F[2][5] = dt;

	mult_K_block33(F, 0, 6, &Ca3, &A);

	mult_K_block33(F, 3, 6, &Ca3, &FF);

	set_K_block33(F, 6, 6, &E);
}

/* currently unused
static void symmetricize(float a[EKF_N][EKF_N])
{
	for (int i = 0; i < EKF_N; ++i) {
		for (int j = 0; j < EKF_N; ++j) {
			float mean = (a[i][j] + a[j][i]) / 2;
			a[i][j] = mean;
			a[j][i] = mean;
		}
	}
}
*/

void addQ(double dt, struct quat q, struct vec ew, struct vec ea, float Q[EKF_N][EKF_N])
{
	// optimize diagonal mmult or maybe A Diag A' ?
	static float Qdiag[EKF_DISTURBANCE][EKF_DISTURBANCE];
	ZEROARR(Qdiag);
	Qdiag[0][0] = Qdiag[1][1] = Qdiag[2][2] = gyro_var_xyz;
	Qdiag[3][3] = Qdiag[4][4] = Qdiag[5][5] = acc_var_xyz;

	static float G[EKF_N][EKF_DISTURBANCE];
	ZEROARR(G);
	struct mat33 quat_by_gyro = meyescl(-1);
	struct mat33 vel_by_acc = mneg(quat2rotmat(qinv(q)));
	set_G_block33(G, 6, 0, &quat_by_gyro);
	set_G_block33(G, 3, 3, &vel_by_acc);

	static float QGt[EKF_DISTURBANCE][EKF_N];
	ZEROARR(QGt);
	SGEMM2D('n', 't', EKF_DISTURBANCE, EKF_N, EKF_DISTURBANCE, 1.0, Qdiag, G, 0.0, QGt);

	SGEMM2D('n', 'n', EKF_N, EKF_N, EKF_DISTURBANCE, dt, G, QGt, 1.0, Q);

	// debug only
	//float Qd[EKF_N][EKF_N];
	//ZEROARR(Qd);
	//SGEMM2D('n', 'n', EKF_N, EKF_N, EKF_DISTURBANCE, 1.0, G, QGt, 0.0, Qd);
	//checksym("Q", AS_1D(Qd), EKF_N);
}

void ekf_imu(struct ekf const *ekf_prev, struct ekf *ekf, float const acc[3], float const gyro[3], float dt)
{
	//------------------------- integrate dynamics --------------------------//
	// propagate constant states
	//ekf->bias_acc = ekf_prev->bias_acc;
	//ekf->bias_gyro = ekf_prev->bias_gyro;

	// TODO TEMP avoid dumb bugs
	*ekf = *ekf_prev;


	// TODO follow S.Weiss and use means with the previous IMU measurements?
	// or stick with the ultra-simple propagation here?
	// (I think the means is not that helpful)

	// propagate rotation
	// TODO only normalize every N steps to save computation?
	//struct vec const omega = vsub(vloadf(gyro), ekf->bias_gyro);
	struct vec const omega = vloadf(gyro);
	ekf->quat = qnormalize(quat_gyro_update(ekf_prev->quat, omega, dt));

	// compute true acceleration
	//struct vec const acc_imu = vsub(float2vec(acc), ekf->bias_acc);
	struct vec const acc_imu = vloadf(acc);
  struct vec acc_world = qvrot(ekf->quat, acc_imu);
  acc_world.z -= GRAV;

	// propagate position + velocity
	ekf->vel = vadd(ekf_prev->vel, vscl(dt, acc_world));
	ekf->pos = vadd(ekf_prev->pos, vscl(dt, ekf->vel));

  // provide smoothed acceleration as a convenience for other system components
  ekf->acc = vadd(vscl(0.7, ekf_prev->acc), vscl(0.3, acc_world));

	//-------------------------- update covariance --------------------------//
	// TODO should use old quat??
	static float F[EKF_N][EKF_N];
	dynamic_matrix(ekf->quat, omega, acc_imu, dt, F);

	// Pnew = F P Ft + Q
	static float PFt[EKF_N][EKF_N];
	ZEROARR(PFt);
	SGEMM2D('n', 't', EKF_N, EKF_N, EKF_N, 1.0, ekf_prev->P, F, 0.0, PFt);
	SGEMM2D('n', 'n', EKF_N, EKF_N, EKF_N, 1.0, F, PFt, 0.0, ekf->P);
	addQ(dt, ekf->quat, omega, acc_imu, ekf->P);
	//symmetricize(ekf->P);
}

void ekf_vicon(struct ekf const *old, struct ekf *new, float const pos_vicon[3], float const vel_vicon[3], float const quat_vicon[4])//, double *debug)
{
	*new = *old;

	struct vec const p_vicon = vloadf(pos_vicon);
	struct vec const v_vicon = vloadf(vel_vicon);
	struct quat const q_vicon = qloadf(quat_vicon);

	struct quat const q_residual = qqmul(qinv(old->quat), q_vicon);
	struct vec const err_quat = vscl(2.0f / q_residual.w, quatimagpart(q_residual));
	struct vec const err_pos = vsub(p_vicon, old->pos);
	struct vec const err_vel = vsub(v_vicon, old->vel);

	float residual[EKF_M];
	vstoref(err_pos, residual);
	vstoref(err_vel, residual + 3);
	vstoref(err_quat, residual + 6);

	// S = (H P H' + R)  :  innovation
	// H = Identity, so S = (P + R)
	static float S[EKF_M][EKF_M];
	COPYMAT(S, old->P);
	float const Rdiag[EKF_M] =
		{ ext_var_xy, ext_var_xy, ext_var_xy, 
		  ext_var_vel, ext_var_vel, ext_var_vel,
		  ext_var_q, ext_var_q, ext_var_q };
	float R[EKF_M][EKF_M];
	ZEROARR(R);
	for (int i = 0; i < EKF_M; ++i) {
		S[i][i] += Rdiag[i];
		R[i][i] = Rdiag[i];
	}

	// K = P H' S^-1  :  gain
	// H = Identity, so K = P S^-1

	static float Sinv[EKF_M][EKF_M];
	static float scratch[EKF_M];
	cholsl(AS_1D(S), AS_1D(Sinv), scratch, EKF_M);

	static float K[EKF_N][EKF_M];
	ZEROARR(K);
	SGEMM2D('n', 'n', EKF_N, EKF_M, EKF_N, 1.0, old->P, Sinv, 0.0, K);


	// K residual : correction

	static float correction[EKF_N];
	ZEROARR(correction);
	sgemm('n', 'n', EKF_N, 1, EKF_M, 1.0, AS_1D(K), residual, 0.0, correction);

	new->pos = vadd(old->pos, vloadf(correction + 0));
	new->vel = vadd(old->vel, vloadf(correction + 3));
	struct quat error_quat = rpy2quat_small(vloadf(correction + 6));
	new->quat = qnormalize(qqmul(old->quat, error_quat));
	// TODO biases, if we use dem


	// Pnew = (I - KH) P (I - KH)^T + KRK^T  :  covariance update
	// TODO optimize KRK^T with R diagonal
	static float RKt[EKF_M][EKF_N];
	ZEROARR(RKt);
	SGEMM2D('n', 't', EKF_M, EKF_N, EKF_M, 1.0, R, K, 0.0, RKt);

	// KRKt - store in P so we can add to it in-place with SGEMM later
	SGEMM2D('n', 'n', EKF_N, EKF_N, EKF_M, 1.0, K, RKt, 0.0, new->P);

	// I - KH
	// H = Identity, so I - K
	static float IMKH[EKF_N][EKF_N];
	eyeN(AS_1D(IMKH), EKF_N);
	for (int i = 0; i < EKF_N; ++i) {
		for (int j = 0; j < EKF_N; ++j) {
			IMKH[i][j] -= K[i][j];
		}
	}

	static float PIMKHt[EKF_N][EKF_N];
	ZEROARR(PIMKHt);
	SGEMM2D('n', 't', EKF_N, EKF_N, EKF_N, 1.0, old->P, IMKH, 0.0, PIMKHt);

	// recall that new->P already contains KRK^T, and we use beta=1.0 to add in-place
	SGEMM2D('n', 'n', EKF_N, EKF_N, EKF_N, 1.0, IMKH, PIMKHt, 1.0, new->P);
}
