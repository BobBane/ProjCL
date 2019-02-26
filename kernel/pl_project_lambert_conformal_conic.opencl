__kernel void pl_project_lambert_conformal_conic_s(
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double scale,
    double x0,
    double y0,

    double lambda0,                                                   
    double rho0,
    double c,
    double n)
{
    int i = get_global_id(0);
    
    double8 lambda = radians(xy_in[i].even) - lambda0;
    double8 phi    = radians(xy_in[i].odd);
    
    double8 x, y;
    
    double8 rho, sinLambda, cosLambda;

    rho = c * powr((1.+tan(.5 * phi))/(1.-tan(.5 * phi)), -n);
    sinLambda = sincos(lambda * n, &cosLambda);
    
    x = rho * sinLambda;
    y = rho0 - rho * cosLambda;
    
    xy_out[i].even = x0 + scale * x;
	xy_out[i].odd  = y0 + scale * y;
}

__kernel void pl_unproject_lambert_conformal_conic_s(
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double scale,
    double x0,
    double y0,

    double lambda0,                                                     
    double rho0,
    double c,
    double n)
{
	int i = get_global_id(0);
	
	double8 x = (xy_in[i].even - x0) / scale;
	double8 y = (xy_in[i].odd - y0) / scale;

	double8 lambda, phi;
    
    double8 rho;
    
    y = rho0 - y;
    
    rho = copysign(hypot(x, y), n);
    
    phi = select(copysign((double)M_PI_2F, n), -atan(sinh(log(rho/c)/n)), rho != 0.);
    lambda = atan2(x * copysign((double)1., n), y * copysign((double)1., n)) / n;
    
	xy_out[i].even = degrees(pl_mod_pi(lambda + lambda0));
	xy_out[i].odd = degrees(phi);
}

__kernel void pl_project_lambert_conformal_conic_e(
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double ecc,
    double ecc2,
    double one_ecc2,

    double scale,
    double x0,
    double y0,

    double lambda0,                                                   
    double rho0,
    double c,
    double n)
{
    int i = get_global_id(0);

    double8 lambda = radians(xy_in[i].even) - lambda0;
    double8 phi    = radians(xy_in[i].odd);
    
    double8 x, y;
    
    double8 rho, sinLambda, cosLambda;
    double8 esinphi = ecc * sin(phi);
    
    rho = c * powr((1.-tan(0.5 * phi))/(1.+tan(0.5 * phi)), n) * powr((1.+esinphi)/(1.-esinphi), .5 * ecc * n);

    sinLambda = sincos(lambda * n, &cosLambda);
    
    x = rho * sinLambda;
    y = rho0 - rho * cosLambda;
    
    xy_out[i].even = x0 + scale * x;
	xy_out[i].odd  = y0 + scale * y;
}

__kernel void pl_unproject_lambert_conformal_conic_e(
    __global double16 *xy_in,
    __global double16 *xy_out,
    const unsigned int count,

    double ecc,
    double ecc2,
    double one_ecc2,

    double scale,
    double x0,
    double y0,

    double lambda0,                                                     
    double rho0,
    double c,
    double n)
{
	int i = get_global_id(0);
	
	double8 x = (xy_in[i].even - x0) / scale;
	double8 y = (xy_in[i].odd - y0) / scale;
    
	double8 lambda, phi;
    
    double8 rho;
        
    y = rho0 - y;
    
    rho = copysign(hypot(x, y), n);
    
    phi = select(copysign((double)M_PI_2F, n), pl_phi2(log(rho/c)/n, ecc), rho != 0.);
    lambda = atan2(x * copysign((double)1., n), y * copysign((double)1., n)) / n;
    
	xy_out[i].even = degrees(pl_mod_pi(lambda + lambda0));
	xy_out[i].odd = degrees(phi);
}

