
__kernel void pl_project_lambert_azimuthal_equal_area_s(
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double scale,
    double x0,
    double y0,

    double phi0,
    double lambda0,

    double sinPhi0,
    double cosPhi0
) {
    int i = get_global_id(0);

    double8 lambda = radians(xy_in[i].even) - lambda0;
    double8 phi    = radians(xy_in[i].odd);

    double8 x, y;
    double8 b, sinLambda, cosLambda, sinPhi, cosPhi;

    sinLambda = sincos(lambda, &cosLambda);
    sinPhi = sincos(phi, &cosPhi);

    b = sqrt(2. / (1. + sinPhi0 * sinPhi + cosPhi0 * cosPhi * cosLambda));
    x = b * cosPhi * sinLambda;
    y = b * (cosPhi0 * sinPhi - sinPhi0 * cosPhi * cosLambda);

    xy_out[i].even = x0 + scale * x;
    xy_out[i].odd  = y0 + scale * y;
}

__kernel void pl_unproject_lambert_azimuthal_equal_area_s(
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double scale,
    double x0,
    double y0,

    double phi0,
    double lambda0,

    double sinPhi0,
    double cosPhi0
) {
    int i = get_global_id(0);

    double8 x = (xy_in[i].even - x0) / scale;
    double8 y = (xy_in[i].odd - y0) / scale;

    double8 lambda, phi;

    double8 rho2;
    double8 sinC, cosC;

    rho2 = x*x + y*y;

    cosC = 1. - 0.5 * rho2;
    sinC = sqrt(1. - 0.25 * rho2); // actually, sin(c) / rho

    phi = asin(cosC * sinPhi0 + y * sinC * cosPhi0);
    lambda = atan2(x * sinC, cosPhi0 * cosC - y * sinPhi0 * sinC);

    xy_out[i].even = degrees(pl_mod_pi(lambda + lambda0));
    xy_out[i].odd = degrees(phi);
}
 
__kernel void pl_project_lambert_azimuthal_equal_area_e(
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double ecc,
    double ecc2,
    double one_ecc2,

    double scale,
    double x0,
    double y0,

    double phi0,
    double lambda0,
    double qp,
    double sinB1,
    double cosB1,

    double rq,
    double4 apa,
    double dd,
    double xmf,
    double ymf
) {
    int i = get_global_id(0);

    double8 lambda = radians(xy_in[i].even) - lambda0;
    double8 phi    = radians(xy_in[i].odd);

    double8 x, y;
    double8 b, sinLambda, cosLambda, sinPhi, sinB, cosB;

    sinLambda = sincos(lambda, &cosLambda);
    sinPhi = sin(phi);

    sinB = pl_qsfn(sinPhi, ecc, one_ecc2) / qp;
    cosB = sqrt(1. - sinB * sinB);

    b = sqrt(2. / (1. + sinB1 * sinB + cosB1 * cosB * cosLambda));

    x = xmf * b * cosB * sinLambda;
    y = ymf * b * (cosB1 * sinB - sinB1 * cosB * cosLambda);

    xy_out[i].even = x0 + scale * x;
    xy_out[i].odd  = y0 + scale * y;
}

__kernel void pl_unproject_lambert_azimuthal_equal_area_e(
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double ecc,
    double ecc2,
    double one_ecc2,

    double scale,
    double x0,
    double y0,

    double phi0,
    double lambda0,
                                                          
    double qp,

    double sinB1,
    double cosB1,

    double rq,
    double4 apa,

    double dd,
    double xmf,
    double ymf
) {
    int i = get_global_id(0);

    double8 x = (xy_in[i].even - x0) / scale;
    double8 y = (xy_in[i].odd - y0) / scale;

    double8 lambda, phi;

    double8 cosCe, sinCe, rho2, beta;
        
    x /= dd;
    y *= dd;

    rho2 = (x*x + y*y) / rq / rq;

    cosCe = 1. - 0.5 * rho2;
    sinCe = sqrt(1. - 0.25 * rho2) / rq; // rather, sin(Ce) / rho

    beta = asin(cosCe * sinB1 + y * sinCe * cosB1);

    lambda = atan2(x * sinCe, cosB1 * cosCe - y * sinB1 * sinCe);
    
    phi = (beta + apa.s0 * sin(2. * beta) + apa.s1 * sin(4. * beta) + apa.s2 * sin(6. * beta));
     
    xy_out[i].even = degrees(pl_mod_pi(lambda + lambda0));
    xy_out[i].odd = degrees(phi);
}

