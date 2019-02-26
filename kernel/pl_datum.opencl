
#define AD_C       1.0026000 /* Toms region 1 constant */


__kernel void pl_cartesian_apply_affine_transform(
    __global double8 *x_rw,
    __global double8 *y_rw,
    __global double8 *z_rw,
    double16 matrix
) {
    int i = get_global_id(0);

    double8 x = x_rw[i];
    double8 y = y_rw[i];
    double8 z = z_rw[i];
    
    x_rw[i] = matrix.lo.lo.x * x + matrix.lo.lo.y * y + matrix.lo.lo.z * z + matrix.lo.lo.w;
    y_rw[i] = matrix.lo.hi.x * x + matrix.lo.hi.y * y + matrix.lo.hi.z * z + matrix.lo.hi.w;
    z_rw[i] = matrix.hi.lo.x * x + matrix.hi.lo.y * y + matrix.hi.lo.z * z + matrix.hi.lo.w;
}

__kernel void pl_geodesic_to_cartesian(
    __global double16 *lp_in,
    __global double8 *x_out,
    __global double8 *y_out,
    __global double8 *z_out,
                                    
    double ecc,
    double ecc2,
    double one_ecc2,
                                    
    double major_axis,
    double minor_axis
) {
    int i = get_global_id(0);
    
    double8 lam = radians(lp_in[i].even);
    double8 phi = radians(lp_in[i].odd);
    
    double8 x, y, z, r;
    
    double8 sinPhi, cosPhi, sinLambda, cosLambda;
    
    sinPhi = sincos(phi, &cosPhi);
    sinLambda = sincos(lam, &cosLambda);
    
    r = major_axis / sqrt(1. - ecc2 * sinPhi * sinPhi);
    x = r * cosPhi * cosLambda;
    y = r * cosPhi * sinLambda;
    z = r * one_ecc2 * sinPhi;
    
    x_out[i] = x;
    y_out[i] = y;
    z_out[i] = z;
}

__kernel void pl_cartesian_to_geodesic(
    __global double8 *x_in,
    __global double8 *y_in,
    __global double8 *z_in,
    __global double16 *lp_out,
                                    
    double ecc,
    double ecc2,
    double one_ecc2,

    double major_axis,
    double minor_axis                                    
) {
    int i = get_global_id(0);

    double8 X = x_in[i];
    double8 Y = y_in[i];
    double8 Z = z_in[i];
    /*
     * The method used here is derived from 'An Improved Algorithm for
     * Geocentric to Geodetic Coordinate Conversion', by Ralph Toms, Feb 1996
     */
    
    /* Note: Variable names follow the notation used in Toms, Feb 1996 */
    
    double8 W;        /* distance from Z axis */
    double8 T0;       /* initial estimate of vertical component */
    double8 T1;       /* corrected estimate of vertical component */
    double8 S0;       /* initial estimate of horizontal component */
    double8 Sin_B0;   /* sin(B0), B0 is estimate of Bowring aux variable */
    double8 Sin3_B0;  /* cube of sin(B0) */
    double8 Cos_B0;   /* cos(B0) */
    double8 Sum;      /* numerator of cos(phi1) */
    
    double8 lambda, phi;
   
    lambda = select(select((double8)M_PI_2F, (double8)-M_PI_2F, Y <= 0.), atan2(Y, X), X != 0.);

    W = hypot(X, Y);
    T0 = Z * AD_C;
    S0 = hypot(T0, W);
    Sin_B0 = T0 / S0;
    Cos_B0 = W / S0;
    Sin3_B0 = Sin_B0 * Sin_B0 * Sin_B0;
    T1 = Z + minor_axis * ecc2 / one_ecc2 * Sin3_B0;
    Sum = W - major_axis * ecc2 * Cos_B0 * Cos_B0 * Cos_B0;

    phi = atan2(T1, Sum);
    
    lp_out[i].even = degrees(lambda);
    lp_out[i].odd = degrees(phi);
}
