
/* Adapted from FORTRAN code in Cengizhan Ipbuker and I. Ozturg Bildirici, 
 * "Computer Program for the Inverse Transformation of the Winkel Projection"
 * Journal of Survey Engineering, Vol. 131 No. 4 (2005) pp. 125-129
 *
 * Modified to use a better initial guess for longitude.
 * And to fix a stupid bug in df1lm in their program.
 */

__kernel void pl_project_winkel_tripel_s(
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double scale, double x0, double y0,
    double lambda0,                                                   
    double cosphi1)
{
    int i = get_global_id(0);
    
    double8 lambda = radians(xy_in[i].even) - lambda0;
    double8 phi    = radians(xy_in[i].odd);
    
    double8 x, y, lambda2, d, cosD, dOverSinD;
    
    double8 cosPhi, sinPhi, cosLambda2, sinLambda2;

    lambda2 = .5 * lambda;
    sinPhi = sincos(phi, &cosPhi);
    sinLambda2 = sincos(lambda2, &cosLambda2);
    
    cosD = cosPhi * cosLambda2;
    d = acos(cosD);
    dOverSinD = select(d / sqrt(1. - cosD * cosD), 1., d == 0.);

    x = lambda2 * cosphi1 + dOverSinD * cosPhi * sinLambda2;
    y =    0.5 * (phi + dOverSinD * sinPhi);

    xy_out[i].even = x0 + scale * x;
    xy_out[i].odd = y0 + scale * y;
}

__kernel void pl_unproject_winkel_tripel_s(
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double scale, double x0, double y0,
    double lambda0,                                                   
    double cosphi1)
{
    int i = get_global_id(0);
    
    double8 x = (xy_in[i].even - x0) / scale;
    double8 y = (xy_in[i].odd - y0) / scale;
    
    double8 phi = y;
    double8 cosPhi;
    double8 sinPhi = sincos(phi, &cosPhi);
    double8 lambda = 2. * x / (cosPhi + cosphi1);
    
    double8 f1, f2;
    double8 c, d;
    double8 invC12, invC;
    double8 sinLambda;
    double8 sinLambda2, cosLambda2; 
    double8 sin2Phi;
    double8 df1phi, df1lam, df2phi, df2lam;
    double8 dPhi, dLam;
    double8 invDet;
    double8 dInvC32;
    
    int iter = WINKEL_TRIPEL_N_ITER;
    do {
        sin2Phi = 2. * sinPhi * cosPhi;
        sinLambda2 = sincos(.5 * lambda, &cosLambda2);
        sinLambda = 2. * sinLambda2 * cosLambda2;

        d = acos(cosPhi*cosLambda2);
        c = sin(d);
        invC = 1. / c / c;
        invC12 = 1. / c;
        dInvC32 = d * invC * invC12;
        f1 = d * cosPhi * sinLambda2 * invC12 + .5 * lambda * cosphi1 - x;
        f2 = .5 * d * sinPhi * invC12 + .5 * phi - y;
        
        df1phi = .25 * sinLambda * sin2Phi * invC 
            - dInvC32 * sinPhi * sinLambda2;
        df1lam = .5 * (cosPhi * cosPhi * sinLambda2 * sinLambda2 * invC
                        + dInvC32 * cosPhi * cosLambda2 * sinPhi * sinPhi
                        + cosphi1);
        
        df2phi = .5 * (sinPhi * sinPhi * cosLambda2 * invC
                        + dInvC32 * sinLambda2 * sinLambda2 * cosPhi
                        + 1.);
        df2lam = .125 * (sin2Phi * sinLambda2 * invC
                          - dInvC32 * sinPhi * cosPhi * cosPhi * sinLambda);

        invDet = 1. / (df1phi * df2lam - df2phi * df1lam);
        
        dPhi = -(f1 * df2lam - f2 * df1lam) * invDet;
        dLam = -(f2 * df1phi - f1 * df2phi) * invDet;
        
        phi += dPhi;
        lambda += dLam; 
        
        sinPhi = sincos(phi, &cosPhi);
    } while (--iter);
    
    xy_out[i].even = degrees(lambda + lambda0);
    xy_out[i].odd = degrees(phi);
}
