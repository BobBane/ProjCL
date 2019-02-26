#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define EPS 5e-14
#define EPS10   1.e-10
#define EPS8    1.e-8
#define EPS7    1.e-7
#define EPS6    1.e-6
#define ITOL 1.e-12
#define TOL7 1.e-7
#define TOL6 1.e-6
#define TOL5 1.e-5
#define I_ITER 20

#define AMERICAN_POLYCONIC_N_ITER 6
#define ALBERS_EQUAL_AREA_N_ITER  6
#define OBLIQUE_STEREOGRAPHIC_N_ITER 6
#define WINKEL_TRIPEL_N_ITER 4

#define M_PIF	         3.1415926535897932384626
#define M_PI_2F          1.570796326794896557999
#define M_PI_4F			 0.78539816339744833

/* Robinson */
#define FXC     0.8487
#define FYC     1.3523
#define C1      11.45915590261646417544
#define RC1     0.08726646259971647884
#define NODES   18
#define ONEEPS  1.000001

double8 pl_qsfn(double8 sinphi, double e, double one_es);
double8 pl_phi2(double8 log_ts, double e);
double8 pl_mod_pi(double8 phi);
double4 pl_interpolate_cubic4(double X, double4 A, double4 B, double4 C, double4 D);

double8 pl_qsfn(double8 sinphi, double e, double one_es) {
	double8 con = e * sinphi;
	return one_es * (sinphi / (1. - con * con) + atanh(con) / e);
}

double8 pl_phi2(double8 log_ts, double e) {
	double8 Phi, con, dphi;
	int i;
	
	Phi = -atan(sinh(log_ts));
	for (i = I_ITER; i; --i) {
		con = e * sin(Phi);
		dphi = -atan(sinh(log_ts - e * atanh(con))) - Phi;
		Phi += dphi;
		if (all(fabs(dphi) <= ITOL))
			break;
	}
	
	return Phi;
}

double8 pl_mod_pi(double8 phi) {
	return select(phi, phi - copysign(2.*M_PIF, phi), fabs(phi) > M_PIF);
}

double4 pl_interpolate_cubic4(double X, double4 A, double4 B, double4 C, double4 D) {
    return B + 0.5 * X * (C - A + X * (2. * A - 5. * B + 4. * C - D + X * (3. * (B - C) + D - A)));
}


