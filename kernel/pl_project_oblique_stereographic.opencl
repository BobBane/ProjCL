double8 srat(double8 esinp, double exp1);
double8 phi_sph2ell(double8 phi, double ecc, double log_k0, double c0);

double8 phi_sph2ell(double8 phi, double ecc, double log_k0, double c0) {
    int i;
    double8 log_num;
    double8 phi_ell;
    
    i = OBLIQUE_STEREOGRAPHIC_N_ITER;
    log_num = (asinh(tan(phi)) - log_k0)/c0;
    phi_ell = phi;

    do {
        phi = phi_ell;
        phi_ell = atan(sinh(log_num + ecc * atanh(ecc * sin(phi))));
    } while (any(fabs(phi_ell - phi) > TOL7) && --i);

    return phi_ell;
}

__kernel void pl_project_oblique_stereographic_e(
	__global double16 *xy_in,
	__global double16 *xy_out,
	const unsigned int count,

	double ecc,
	double ecc2,
	double one_ecc2,

    double scale,
    double x0,
    double y0,

    double c0,
    double log_k0,

    double lambda0,
    double sinPhiC0,
    double cosPhiC0
    ) {
	int i = get_global_id(0);

	double8 lambda_ell = radians(xy_in[i].even) - lambda0;
	double8 phi_ell    = radians(xy_in[i].odd);

    /* Project ellipsoid onto sphere */
    double8 lambda = c0 * lambda_ell;
    double8 phi = atan(sinh(log_k0 + c0 * (asinh(tan(phi_ell)) - ecc * atanh(ecc * sin(phi_ell)))));
    // Gudermannian Function gd(x) = 2 atan(exp(x)) - M_PI_2 = asin(tanh(x)) = atan(sinh(x))
    // also log(tan(0.5 * x + M_PI_4)) = asinh(tan(x))

    /* Project sphere onto plane */
    double8 sinPhi, cosPhi;
    double8 sinLambda, cosLambda;
    double8 x, y, k;

    sinPhi = sincos(phi, &cosPhi);
    sinLambda = sincos(lambda, &cosLambda);

    k = scale / (1. + sinPhiC0 * sinPhi + cosPhiC0 * cosPhi * cosLambda);

    x = cosPhi * sinLambda;
    y = cosPhiC0 * sinPhi - sinPhiC0 * cosPhi * cosLambda;

    xy_out[i].even = x0 + k * x;
    xy_out[i].odd = y0 + k * y;
}

__kernel void pl_unproject_oblique_stereographic_e(
	__global double16 *xy_in,
	__global double16 *xy_out,
	const unsigned int count,

	double ecc,
	double ecc2,
	double one_ecc2,

    double scale,
    double x0,
    double y0,

    double c0,
    double log_k0,

	double lambda0,
    double sinPhiC0,
    double cosPhiC0
    ) {
	int i = get_global_id(0);

	double8 x = (xy_in[i].even - x0) / scale;
	double8 y = (xy_in[i].odd - y0) / scale;
	
	double8 lambda, phi;
    double8 rho;
    double8 sinPhiC, cosPhiC;

    rho = hypot(x, y);
    sinPhiC = sincos(2. * atan(rho), &cosPhiC);

    /* Project plane onto sphere */
    phi = asin(select(cosPhiC * sinPhiC0 + y * sinPhiC * cosPhiC0 / rho,
                sinPhiC0, rho == 0.));
    lambda = atan2(x * sinPhiC, rho * cosPhiC0 * cosPhiC - y * sinPhiC0 * sinPhiC);

    /* Project sphere onto ellipsoid */
    lambda = lambda / c0;
    phi = phi_sph2ell(phi, ecc, log_k0, c0);

    xy_out[i].even = degrees(pl_mod_pi(lambda + lambda0));
    xy_out[i].odd = degrees(phi);
}
