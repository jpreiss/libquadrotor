void solve_poly7_boundary(
	float const vals_init[4],
	float const vals_end[4],
	float const time_end,
	float coefs[8])
{
	float const T1 = time_end;

	float const T2 = T1 * T1;
	float const T3 = T2 * T1;
	float const T4 = T3 * T1;
	float const T5 = T4 * T1;
	float const T6 = T5 * T1;
	float const T7 = T6 * T1;

	float const d0_0 = vals_init[0];
	float const d0_1 = vals_init[1];
	float const d0_2 = vals_init[2];
	float const d0_3 = vals_init[3];

	float const d1_0 = vals_end[0];
	float const d1_1 = vals_end[1];
	float const d1_2 = vals_end[2];
	float const d1_3 = vals_end[3];

	coefs[0] = d0_0;
	coefs[1] = d0_1;
	coefs[2] = d0_2/2;
	coefs[3] = d0_3/6;
	coefs[4] = (-120*T1*d0_1 - 90*T1*d1_1 - 30*T2*d0_2 + 15*T2*d1_2 - 4*T3*d0_3 - T3*d1_3 - 210*d0_0 + 210*d1_0)/(6*T4);
	coefs[5] = (45*T1*d0_1 + 39*T1*d1_1 + 10*T2*d0_2 - 7*T2*d1_2 + T3*d0_3 + T3*d1_3/2 + 84*d0_0 - 84*d1_0)/T5;
	coefs[6] = (-216*T1*d0_1 - 204*T1*d1_1 - 45*T2*d0_2 + 39*T2*d1_2 - 4*T3*d0_3 - 3*T3*d1_3 - 420*d0_0 + 420*d1_0)/(6*T6);
	coefs[7] = (60*T1*d0_1 + 60*T1*d1_1 + 12*T2*d0_2 - 12*T2*d1_2 + T3*d0_3 + T3*d1_3 + 120*d0_0 - 120*d1_0)/(6*T7);

}