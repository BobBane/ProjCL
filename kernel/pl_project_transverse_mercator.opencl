
/* See Karney, "Transverse Mercator with an accuracy of a few nanometers"
 * J. Geodesy 85(8), 475-485 (2011)
 * https://doi.org/10.1007/s00190-011-0445-3
 * https://arxiv.org/pdf/1002.1417.pdf
 */

__kernel void pl_project_transverse_mercator_s (
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double scale, double x0, double y0,
    double lambda0
) {
    int i = get_global_id(0);
    
    double8 lambda = radians(xy_in[i].even) - lambda0;
    double8 phi    = radians(xy_in[i].odd);
    
    double8 x, y, tau, sinLambda, cosLambda;
    
    sinLambda = sincos(lambda, &cosLambda);
    
    tau = tan(phi);
    y = atan2(tau, cosLambda);
    x = asinh(sinLambda / hypot(tau, cosLambda));
    
    xy_out[i].even = x0 + scale * x;
    xy_out[i].odd = y0 + scale * y;
}

__kernel void pl_unproject_transverse_mercator_s (
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double scale, double x0, double y0,
    double lambda0
) {
    int i = get_global_id(0);

    double8 x = (xy_in[i].even - x0) / scale;
    double8 y = (xy_in[i].odd - y0) / scale;
    
    double8 phi, lambda, sinhX, sinY, cosY;
    
    sinhX = sinh(x);
    sinY = sincos(y, &cosY);

    lambda = atan2(sinhX, cosY);
    phi = atan2(sinY, hypot(sinhX, cosY));
    
    xy_out[i].even = degrees(pl_mod_pi(lambda + lambda0));
    xy_out[i].odd  = degrees(phi);
}

__kernel void pl_project_transverse_mercator_e (
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double ecc,
    double ecc2,
    double one_ecc2,

    double scale, double x0, double y0,
    double lambda0,
    double8 krueger_alpha,
    double8 krueger_beta
) {
    int i = get_global_id(0);

    double8 lambda = radians(xy_in[i].even) - lambda0;
    double8 phi    = radians(xy_in[i].odd);

    double8 x, y, xi, eta, tau, tau1, sigma, sinLambda, cosLambda;
    double f, n;

    double8 sin2, cos2, sinh2, cosh2;
    double8 sin4, cos4, sinh4, cosh4;
    double8 sin6, cos6, sinh6, cosh6;
    double8 sin8, cos8, sinh8, cosh8;

    f = 1. - sqrt(one_ecc2);
    n = f / (2. - f);

    sinLambda = sincos(lambda, &cosLambda);
    
    tau = tan(phi);
    sigma = sinh(ecc * atanh(ecc * tau / hypot(1., tau)));

    tau1 = tau * hypot(1., sigma) - sigma * hypot(1., tau);

    xi = atan2(tau1, cosLambda);
    eta = asinh(sinLambda / hypot(tau1, cosLambda));

    sin2 = sincos(2. * xi, &cos2);

    sin4 = 2. * sin2 * cos2;
    cos4 = 2. * cos2 * cos2 - 1.;

    sin6 = sin4 * cos2 + cos4 * sin2;
    cos6 = cos4 * cos2 - sin4 * sin2;

    sin8 = 2. * sin4 * cos4;
    cos8 = 2. * cos4 * cos4 - 1.;

    sinh2 = sinh(2. * eta);
    cosh2 = cosh(2. * eta);

    sinh4 = 2. * sinh2 * cosh2;
    cosh4 = 2. * cosh2 * cosh2 - 1.;

    sinh6 = sinh4 * cosh2 + cosh4 * sinh2;
    cosh6 = cosh4 * cosh2 + sinh4 * sinh2;

    sinh8 = 2. * sinh4 * cosh4;
    cosh8 = 2. * cosh4 * cosh4 - 1.;

    y = xi;
    y += krueger_alpha.s0 * sin2 * cosh2;
    y += krueger_alpha.s1 * sin4 * cosh4;
    y += krueger_alpha.s2 * sin6 * cosh6;
    y += krueger_alpha.s3 * sin8 * cosh8;

    x = eta;
    x += krueger_alpha.s0 * cos2 * sinh2;
    x += krueger_alpha.s1 * cos4 * sinh4;
    x += krueger_alpha.s2 * cos6 * sinh6;
    x += krueger_alpha.s3 * cos8 * sinh8;

    xy_out[i].even = x0 + scale * x;
    xy_out[i].odd  = y0 + scale * y;
}

__kernel void pl_unproject_transverse_mercator_e (
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double ecc,
    double ecc2,
    double one_ecc2,

    double scale, double x0, double y0,
    double lambda0,
    double8 krueger_alpha,
    double8 krueger_beta
) {
    int i = get_global_id(0);

    double8 x = (xy_in[i].even - x0) / scale;
    double8 y = (xy_in[i].odd - y0) / scale;

    double8 phi, lambda, sinhX, sinY, cosY;
    double8 xi, eta;
    double8 tau0, tau, sigma, tauP, dtau;

    double8 sin2, cos2, sinh2, cosh2;
    double8 sin4, cos4, sinh4, cosh4;
    double8 sin6, cos6, sinh6, cosh6;
    double8 sin8, cos8, sinh8, cosh8;

    sin2 = sincos(2. * y, &cos2);

    sin4 = 2. * sin2 * cos2;
    cos4 = 2. * cos2 * cos2 - 1.;

    sin6 = sin4 * cos2 + cos4 * sin2;
    cos6 = cos4 * cos2 - sin4 * sin2;

    sin8 = 2. * sin4 * cos4;
    cos8 = 2. * cos4 * cos4 - 1.;

    sinh2 = sinh(2. * x);
    cosh2 = cosh(2. * x);

    sinh4 = 2. * sinh2 * cosh2;
    cosh4 = 2. * cosh2 * cosh2 - 1.;

    sinh6 = sinh4 * cosh2 + cosh4 * sinh2;
    cosh6 = cosh4 * cosh2 + sinh4 * sinh2;

    sinh8 = 2. * sinh4 * cosh4;
    cosh8 = 2. * cosh4 * cosh4 - 1.;

    xi = y;
    xi -= krueger_beta.s0 * sin2 * cosh2;
    xi -= krueger_beta.s1 * sin4 * cosh4;
    xi -= krueger_beta.s2 * sin6 * cosh6;
    xi -= krueger_beta.s3 * sin8 * cosh8;

    eta = x;
    eta -= krueger_beta.s0 * cos2 * sinh2;
    eta -= krueger_beta.s1 * cos4 * sinh4;
    eta -= krueger_beta.s2 * cos6 * sinh6;
    eta -= krueger_beta.s3 * cos8 * sinh8;

    sinhX = sinh(eta);
    sinY = sincos(xi, &cosY);

    tau = tau0 = sinY / hypot(sinhX, cosY);

    /* Newton's method (1 iteration) */
    sigma = sinh(ecc * atanh(ecc * tau / hypot(1., tau)));
    tauP = tau * hypot(1., sigma) - sigma * hypot(1., tau);
    dtau = (tau0 - tauP) / hypot(1., tauP) * (1. + one_ecc2 * tau * tau) / (one_ecc2 * hypot(1., tau));
    tau += dtau;

    lambda = atan2(sinhX, cosY);
    phi = atan(tau);

    xy_out[i].even = degrees(pl_mod_pi(lambda + lambda0));
    xy_out[i].odd  = degrees(phi);
}
