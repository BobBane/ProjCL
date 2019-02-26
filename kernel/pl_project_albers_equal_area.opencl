
double8 phi1_(double8 qs, double Te, double Tone_es);

double8 phi1_(double8 qs, double Te, double Tone_es) {
	int i;
	double8 Phi, sinphi, cosphi, con, com, dphi;

	Phi = asin(.5 * qs);
	if (Te < EPS7)
		return( Phi );
	
	i = ALBERS_EQUAL_AREA_N_ITER;
	
	do {
		sinphi = sincos(Phi, &cosphi);
		con = Te * sinphi;
		com = 1. - con * con;
		dphi = .5 * com * com / cosphi * (qs / Tone_es - sinphi / com - atanh(con) / Te);
		/* dphi = .5 * com * com / cosphi * (qs / Tone_es -
		   sinphi / com + .5 / Te * log ((1. - con) /
		   (1. + con))); */
		Phi += dphi;
	} while (any(fabs(dphi) > TOL7) && --i);
	
	return Phi;
}

__kernel void pl_project_albers_equal_area_s(
	__global double16 *xy_in,
	__global double16 *xy_out,
	const unsigned int count,
	
    double scale,
    double x0,
    double y0,

	double lambda0,
	double rho0,
	double c,
	double n
) {
	int i = get_global_id(0);
	
	double8 lambda = radians(xy_in[i].even) - lambda0;
	double8 phi    = radians(xy_in[i].odd);
	
	double8 x, y;
	
	double8 rho, sinLambda, cosLambda;

	rho = sqrt(c - 2 * n * sin(phi)) / n;

	sinLambda = sincos(lambda * n, &cosLambda);
	
	x = rho * sinLambda;
	y = rho0 - rho * cosLambda;
	
	xy_out[i].even = x0 + scale * x;
	xy_out[i].odd  = y0 + scale * y;
}

__kernel void pl_unproject_albers_equal_area_s(
	__global double16 *xy_in,
	__global double16 *xy_out,
	const unsigned int count,
	
    double scale,
    double x0,
    double y0,

	double lambda0,
	double rho0,
	double c,
	double n
) {
	int i = get_global_id(0);
	
	double8 x = (xy_in[i].even - x0) / scale;
	double8 y = (xy_in[i].odd - y0) / scale;
	
	double8 lambda, phi;
	
	y = rho0 - y;

    phi = 0.5 * (c / n - (x * x + y * y) * n);
	
	phi = select(copysign(M_PI_2F, phi), asin(phi), fabs(phi) <= 1.);
	lambda = atan2(x * copysign((double)1., n), y * copysign((double)1., n)) / n;
	
	xy_out[i].even = degrees(pl_mod_pi(lambda + lambda0));
	xy_out[i].odd = degrees(phi);
}

__kernel void pl_project_albers_equal_area_e(
	__global double16 *xy_in,
	__global double16 *xy_out,
	const unsigned int count,
	
	double ecc,
	double ecc2,
	double one_ecc2,
	
	double ec,
	
    double scale,
    double x0,
    double y0,
	double lambda0,
	double rho0,
	double c,
	double n
) {
	int i = get_global_id(0);
	
	double8 lambda = radians(xy_in[i].even) - lambda0;
	double8 phi    = radians(xy_in[i].odd);
	
	double8 x, y;
	
	double8 rho, sinLambda, cosLambda;

	rho = sqrt(c - n * pl_qsfn(sin(phi), ecc, one_ecc2)) / n;

	sinLambda = sincos(lambda * n, &cosLambda);
	
	x = rho * sinLambda;
	y = rho0 - rho * cosLambda;
	
	xy_out[i].even = x0 + scale * x;
	xy_out[i].odd  = y0 + scale * y;
}

__kernel void pl_unproject_albers_equal_area_e(
	__global double16 *xy_in,
	__global double16 *xy_out,
	const unsigned int count,
	
	double ecc,
	double ecc2,
	double one_ecc2,
	
	double ec,
	
    double scale,
    double x0,
    double y0,

	double lambda0,
	double rho0,
	double c,
	double n
) {
	int i = get_global_id(0);
	
	double8 x = (xy_in[i].even - x0) / scale;
	double8 y = (xy_in[i].odd - y0) / scale;
	
	double8 lambda, phi;
	
	y = rho0 - y;

	phi = (c / n - (x * x + y * y) * n);
	
	phi = select(copysign(M_PI_2F, phi), phi1_(phi, ecc, one_ecc2), fabs(ec - fabs(phi)) > TOL7);
	lambda = atan2(x * copysign((double)1., n), y * copysign((double)1., n)) / n;
	
	xy_out[i].even = degrees(pl_mod_pi(lambda + lambda0));
	xy_out[i].odd = degrees(phi);
}
