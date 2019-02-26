__kernel void pl_sincos(
	__global double8 *phi_in,
	__global double16 *sincos_out
) {
	int i = get_global_id(0);
	
	double8 phi = radians(phi_in[i]);
	double8 sinPhi, cosPhi;
	
	sinPhi = sincos(phi, &cosPhi);
	
	sincos_out[i].even = sinPhi;
	sincos_out[i].odd = cosPhi;
}

__kernel void pl_sincos1(
	__global double16 *phi_in,
	__global double16 *sincos_out
) {
	int i = get_global_id(0);
	
	double8 phi = radians(phi_in[i].odd);
	double8 sinPhi, cosPhi;
	
	sinPhi = sincos(phi, &cosPhi);
	
	sincos_out[i].even = sinPhi;
	sincos_out[i].odd = cosPhi;
}

__kernel void pl_inverse_geodesic_s(
	__global double2 *lp1_in,
	__global double16 *lp2_in,
	__global double8 *dist_out,
	
	double radius
) {
	int i = get_global_id(0);
	int j = get_global_id(1);
	
	int j_size = get_global_size(1);
	
	double lam1 = radians(lp1_in[i].even);
	double phi1 = radians(lp1_in[i].odd);
	
	double8 lam2 = radians(lp2_in[j].even);
	double8 phi2 = radians(lp2_in[j].odd);
	
	double cosPhi1 = cos(phi1);
	double8 cosPhi2 = cos(phi2);
	
    double8 dlam = lam2 - lam1;
    double8 dphi = phi2 - phi1;

    double8 shp2 = sin(0.5 * dphi);
    double8 shl2 = sin(0.5 * dlam);

    dist_out[i*j_size+j] = 2. * radius * asin(sqrt(shp2 * shp2 + cosPhi1 * cosPhi2 * shl2 * shl2));
}

__kernel void pl_forward_geodesic_fixed_distance_s(
	__global double2 *lp_in,
	__global double2 *phi_sincos,
	__global double16 *az_sincos,
	__global double16 *lp_out,
	
	double distance,
	double sinDistance,
	double cosDistance
) {
	int i = get_global_id(0);
	int j = get_global_id(1);
	
	int j_size = get_global_size(1);
	
	double lam1 = radians(lp_in[i].s0);
	
	double sinPhi = phi_sincos[i].s0;
	double cosPhi = phi_sincos[i].s1;
	
	double8 sinAz = az_sincos[j].even;
	double8 cosAz = az_sincos[j].odd;
	
	double8 lam2, phi2;
	
	phi2 = asin(sinPhi * cosDistance + cosPhi * sinDistance * cosAz);
	lam2 = lam1 + atan2(sinDistance * sinAz, 
		cosPhi * cosDistance - sinPhi * sinDistance * cosAz);

	lp_out[i*j_size+j].even = degrees(pl_mod_pi(lam2));
	lp_out[i*j_size+j].odd = degrees(phi2);
}

__kernel void pl_forward_geodesic_fixed_angle_s(
	double2 lp_in,
	__global double8 *dist,
	__global double16 *lp_out,
	
	double azimuth,
	double sinAz,
	double cosAz
) {
	int i = get_global_id(0);
	
	double lam1 = radians(lp_in.s0);
	double phi1 = radians(lp_in.s1);
	
	double sinPhi, cosPhi;
	double8 sinDistance, cosDistance;
	double8 lam2, phi2;
	
    sinPhi = sincos(phi1, &cosPhi);
	
    sinDistance = sincos(dist[i], &cosDistance);
	
	phi2 = asin(sinPhi * cosDistance + cosPhi * sinDistance * cosAz);
	lam2 = lam1 + atan2(sinDistance * sinAz, 
		cosPhi * cosDistance - sinPhi * sinDistance * cosAz);
	
	lp_out[i].even = degrees(pl_mod_pi(lam2));
	lp_out[i].odd = degrees(phi2);
}

/* *** SOLUTION OF THE GEODETIC DIRECT PROBLEM AFTER T.VINCENTY */
/* *** MODIFIED RAINSFORD'S METHOD WITH HELMERT'S ELLIPTICAL TERMS */
/* *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL */

/* *** A IS THE SEMI-MAJOR AXIS OF THE REFERENCE ELLIPSOID */
/* *** F IS THE FLATTENING OF THE REFERENCE ELLIPSOID */
/* *** LATITUDES AND LONGITUDES IN RADIANS POSITIVE NORTH AND EAST */
/* *** AZIMUTHS IN RADIANS CLOCKWISE FROM NORTH */
/* *** GEODESIC DISTANCE S ASSUMED IN UNITS OF SEMI-MAJOR AXIS A */

/* *** PROGRAMMED FOR CDC-6600 BY LCDR L.PFEIFER NGS ROCKVILLE MD 20FEB75 */
/* *** MODIFIED FOR SYSTEM 360 BY JOHN G GERGEN NGS ROCKVILLE MD 750608 */

/* *** ...modified for OpenCL by evan miller (CHICAGO IL) */

/* kernel void pl_forward_geodesic_e(
	__global double2 *lp_in,
	__global double2 *phi_sincos,
	__global double16 *az_sincos,
	__global double16 *lp_out,
	
	double distance,
	double sinDistance,
	double cosDistance,
	
	double flattening
) {
	int i = get_global_id(0);
	int j = get_global_id(1);
	
	int j_size = get_global_size(1);
	
	double lam1 = lp_in[i].s0;
	double phi1 = lp_in[i].s1;
	
	double sinPhi = phi_sincos[i].s0;
	double cosPhi = phi_sincos[i].s1;
	
	double8 sinAz = az_sincos[j].even;
	double8 cosAz = az_sincos[j].odd;
	
	double8 az21, lam2, phi2;
	
	double r;

    double8 c, d, e, x, y, cu, cy, cz, sa, su, tu, sy, c2a;



    r = 1. - flattening;
    tu = r * sinPhi / cosPhi;
	az21 = select(atan2(tu, cosAz) * 2., 0., cosAz == 0.);
    cu = 1. / sqrt(tu * tu + 1.);
    su = tu * cu;
    sa = cu * sinAz;
    c2a = -sa * sa + 1.;
    x = sqrt((1. / r / r - 1.) * c2a + 1.) + 1.;
    x = (x - 2.) / x;
    c = 1. - x;
    c = (x * x / 4. + 1.) / c;
    d = (x * .375 * x - 1.) * x;
    tu = distance / r / c;
    y = tu;
    do {
        sy = sincos(y, &cy);
        cz = cos(az21 + y);
        e = cz * cz * 2. - 1.;
        c = y;
        x = e * cy;
        y = e + e - 1.;
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) * d / 4. - cz) *
            sy * d + tu;
    } while (any(fabs(y - c) > EPS7));
    az21 = cu * cy * cosAz - su * sy;
    phi2 = atan2(su * cy + cu * sy * cosAz, r * hypot(sa, az21));
    x = atan2(sy * sinAz, cu * cy - su * sy * cosAz);
    c = ((c2a * -3. + 4.) * flattening + 4.) * c2a * flattening / 16.;
    d = ((e * cy * c + cz) * sy * c + y) * sa;
    lam2 = lam1 + x - (1. - c) * d * flattening;
	
//    az21 = atan2(sa, az21) + M_PIF;
	
	lp_out[i*j_size+j].even = lam2;
	lp_out[i*j_size+j].odd = phi2;
}
*/

