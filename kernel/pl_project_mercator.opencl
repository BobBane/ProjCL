__kernel void pl_project_mercator_s(
	__global double16 *xy_in,
	__global double16 *xy_out,
	const unsigned int count,
		
	double scale, double x0, double y0)
{
	int i = get_global_id(0);
	
	double8 lambda = radians(xy_in[i].even);
	double8 phi    = radians(xy_in[i].odd);
	
	double8 x, y;
	
	x = lambda;
	y = asinh(tan(phi));

	xy_out[i].even = x0 + scale * x;
	xy_out[i].odd  = y0 + scale * y;
}

__kernel void pl_unproject_mercator_s(
	__global double16 *xy_in,
	__global double16 *xy_out,
	const unsigned int count,
		
	double scale, double x0, double y0)
{
	int i = get_global_id(0);
	
	double8 x = (xy_in[i].even - x0) / scale;
	double8 y = (xy_in[i].odd - y0) / scale;
	
	double8 lambda, phi;
	
	phi = atan(sinh(y));
	lambda = x;

	xy_out[i].even = degrees(lambda);
	xy_out[i].odd = degrees(phi);
}

__kernel void pl_project_mercator_e(
	__global double16 *xy_in,
	__global double16 *xy_out,
	const unsigned int count,
	
	double ecc,
	double ecc2,
	double one_ecc2,
	
	double scale, double x0, double y0)
{
	int i = get_global_id(0);
	
	double8 lambda = radians(xy_in[i].even);
	double8 phi    = radians(xy_in[i].odd);
	
	double8 x, y;
	
	x = lambda;
	y = asinh(tan(phi)) - ecc * atanh(ecc * sin(phi));
	
	xy_out[i].even = x0 + scale * x;
	xy_out[i].odd  = y0 + scale * y;
}

__kernel void pl_unproject_mercator_e(
	__global double16 *xy_in,
	__global double16 *xy_out,
	const unsigned int count,
	
	double ecc,
	double ecc2,
	double one_ecc2,
	
	double scale, double x0, double y0)
{
	int i = get_global_id(0);
	
	double8 x = (xy_in[i].even - x0) / scale;
	double8 y = (xy_in[i].odd - y0) / scale;
	
	double8 lambda, phi;
	
	phi = pl_phi2(-y, ecc);
	lambda = x;
	
	xy_out[i].even = degrees(lambda);
	xy_out[i].odd = degrees(phi);
}
