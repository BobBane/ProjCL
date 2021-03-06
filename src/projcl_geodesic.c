//
//  projcl_geodesic.c
//  Magic Maps
//
//  Created by Evan Miller on 3/31/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdlib.h>
#include <projcl/projcl.h>
#include "projcl_run.h"
#include "projcl_util.h"
#include "projcl_kernel.h"
#include "projcl_spheroid.h"

PLForwardGeodesicFixedDistanceBuffer *pl_load_forward_geodesic_fixed_distance_data(PLContext *pl_ctx,
    const double *xy_in, size_t xy_count, const double *az_in, size_t az_count, cl_int *outError)
{
	cl_int error;
	
	PLForwardGeodesicFixedDistanceBuffer *pl_buf;
	if ((pl_buf = malloc(sizeof(PLForwardGeodesicFixedDistanceBuffer))) == NULL) {
		return NULL;
	}
	pl_buf->xy_in = NULL;
	pl_buf->az_in = NULL;
	pl_buf->xy_count = xy_count;
	pl_buf->az_count = az_count;
	pl_buf->phi_sincos = NULL;
	pl_buf->az_sincos = NULL;
	
	
	double *az_pad_in = NULL, *xy_pad_in = NULL;
	
	int az_pad_count = ck_padding(az_count, PL_DOUBLE_VECTOR_SIZE);
	int xy_pad_count = ck_padding(xy_count, PL_DOUBLE_VECTOR_SIZE);
	
	if ((az_pad_in = malloc(az_pad_count * sizeof(double))) == NULL) {
		error = CL_OUT_OF_HOST_MEMORY;
		goto cleanup;
	}
	if ((xy_pad_in = malloc(xy_pad_count * sizeof(double) * 2)) == NULL) {
		error = CL_OUT_OF_HOST_MEMORY;
		goto cleanup;
	}
	
	_pl_copy_pad(az_pad_in, az_pad_count, az_in, az_count);
	_pl_copy_pad(xy_pad_in, xy_pad_count * 2, xy_in, xy_count * 2);
	
	pl_buf->xy_in = clCreateBuffer(pl_ctx->ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
                                   sizeof(cl_double) * xy_pad_count * 2, xy_pad_in, &error);
	if (error != CL_SUCCESS)
		goto cleanup;
	
	pl_buf->az_in = clCreateBuffer(pl_ctx->ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
                                   sizeof(cl_double) * az_pad_count, az_pad_in, &error);
	if (error != CL_SUCCESS)
		goto cleanup;
	
	pl_buf->phi_sincos = clCreateBuffer(pl_ctx->ctx, CL_MEM_READ_WRITE, 
                                        sizeof(cl_double) * xy_pad_count * 2, NULL, &error);
	if (error != CL_SUCCESS)
		goto cleanup;
	
	pl_buf->az_sincos = clCreateBuffer(pl_ctx->ctx, CL_MEM_READ_WRITE, 
                                       sizeof(cl_double) * az_pad_count * 2, NULL, &error);
	if (error != CL_SUCCESS)
		goto cleanup;
	
	pl_buf->xy_out = clCreateBuffer(pl_ctx->ctx, CL_MEM_WRITE_ONLY, 
									sizeof(cl_double) * xy_count * az_pad_count * 2, NULL, &error);
	if (error != CL_SUCCESS)
		goto cleanup;
	
cleanup:
	if (az_pad_in)
		free(az_pad_in);
	if (xy_pad_in)
		free(xy_pad_in);
	if (error != CL_SUCCESS) {
		if (pl_buf->xy_in)
			clReleaseMemObject(pl_buf->xy_in);
		if (pl_buf->az_in)
			clReleaseMemObject(pl_buf->az_in);
		if (pl_buf->phi_sincos)
			clReleaseMemObject(pl_buf->phi_sincos);
		if (pl_buf->az_sincos)
			clReleaseMemObject(pl_buf->az_sincos);
		if (pl_buf->xy_out)
			clReleaseMemObject(pl_buf->xy_out);
		free(pl_buf);
		if (outError != NULL)
			*outError = error;
		return NULL;
	}
	
	return pl_buf;
}

void pl_unload_forward_geodesic_fixed_distance_data(PLForwardGeodesicFixedDistanceBuffer *pl_buf) {
	clReleaseMemObject(pl_buf->xy_in);
	clReleaseMemObject(pl_buf->az_in);
	clReleaseMemObject(pl_buf->phi_sincos);
	clReleaseMemObject(pl_buf->az_sincos);
	clReleaseMemObject(pl_buf->xy_out);
	free(pl_buf);
}

PLForwardGeodesicFixedAngleBuffer *pl_load_forward_geodesic_fixed_angle_data(PLContext *pl_ctx,
    const double *dist_in, size_t dist_count, cl_int *outError) {
    cl_int error;
    
    PLForwardGeodesicFixedAngleBuffer *pl_buf;
	if ((pl_buf = malloc(sizeof(PLForwardGeodesicFixedAngleBuffer))) == NULL) {
		return NULL;
	}
    
    pl_buf->dist_in = NULL;
    pl_buf->dist_count = dist_count;
    pl_buf->xy_out = NULL;
    
    double *dist_pad_in;
    
    int dist_pad_count = ck_padding(dist_count, PL_DOUBLE_VECTOR_SIZE);
    
    if ((dist_pad_in = malloc(dist_pad_count * sizeof(double))) == NULL) {
        error = CL_OUT_OF_HOST_MEMORY;
        goto cleanup;
    }
    
    _pl_copy_pad(dist_pad_in, dist_pad_count, dist_in, dist_count);
    
    pl_buf->dist_in = clCreateBuffer(pl_ctx->ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
                                     sizeof(cl_double) * dist_pad_count, dist_pad_in, &error);
    if (error != CL_SUCCESS)
        goto cleanup;
    
    pl_buf->xy_out = clCreateBuffer(pl_ctx->ctx, CL_MEM_WRITE_ONLY, 
                                    sizeof(cl_double) * dist_pad_count * 2, NULL, &error);
    if (error != CL_SUCCESS)
        goto cleanup;
    
cleanup:
    if (dist_pad_in)
        free(dist_pad_in);
    if (error != CL_SUCCESS) {
        if (pl_buf->dist_in)
            clReleaseMemObject(pl_buf->dist_in);
        if (pl_buf->xy_out)
            clReleaseMemObject(pl_buf->xy_out);
        free(pl_buf);
        if (outError != NULL)
            *outError = error;
        return NULL;
    }
    return pl_buf;
}

void pl_unload_forward_geodesic_fixed_angle_data(PLForwardGeodesicFixedAngleBuffer *pl_buf) {
    clReleaseMemObject(pl_buf->dist_in);
    clReleaseMemObject(pl_buf->xy_out);
    free(pl_buf);
}

cl_int pl_forward_geodesic_fixed_distance(PLContext *pl_ctx, PLForwardGeodesicFixedDistanceBuffer *pl_buf,
        double *xy_out, PLSpheroid pl_ell, double distance) {
	cl_kernel sincos_kernel = NULL, sincos1_kernel = NULL, fwd_kernel = NULL;
	
	if (_pl_spheroid_is_spherical(pl_ell)) {
		fwd_kernel = pl_find_kernel(pl_ctx, "pl_forward_geodesic_fixed_distance_s");
	} else {
		fwd_kernel = pl_find_kernel(pl_ctx, "pl_forward_geodesic_fixed_distance_e");
	}
	if (fwd_kernel == NULL) {
		return CL_INVALID_KERNEL_NAME;
	}
	
	sincos_kernel = pl_find_kernel(pl_ctx, "pl_sincos");
	if (sincos_kernel == NULL) {
		return CL_INVALID_KERNEL_NAME;
	}
	
	sincos1_kernel = pl_find_kernel(pl_ctx, "pl_sincos1");
	if (sincos_kernel == NULL) {
		return CL_INVALID_KERNEL_NAME;
	}
	
	size_t phiVecCount = ck_padding(pl_buf->xy_count, PL_DOUBLE_VECTOR_SIZE) / PL_DOUBLE_VECTOR_SIZE;
	size_t azVecCount = ck_padding(pl_buf->az_count, PL_DOUBLE_VECTOR_SIZE) / PL_DOUBLE_VECTOR_SIZE;
	
	cl_int error = CL_SUCCESS;
	
	error |= pl_set_kernel_arg_mem(pl_ctx, sincos1_kernel, 0, pl_buf->xy_in);
	error |= pl_set_kernel_arg_mem(pl_ctx, sincos1_kernel, 1, pl_buf->phi_sincos);	
	if (error != CL_SUCCESS) {
		return error;
	}
	
	error = clEnqueueNDRangeKernel(pl_ctx->queue, sincos1_kernel, 1, NULL, 
								   &phiVecCount, NULL, 0, NULL, NULL);
	if (error != CL_SUCCESS) {
		return error;
	}
	
	error |= pl_set_kernel_arg_mem(pl_ctx, sincos_kernel, 0, pl_buf->az_in);
	error |= pl_set_kernel_arg_mem(pl_ctx, sincos_kernel, 1, pl_buf->az_sincos);
	if (error != CL_SUCCESS) {
		return error;
	}
	
	error = clEnqueueNDRangeKernel(pl_ctx->queue, sincos_kernel, 1, NULL, 
								   &azVecCount, NULL, 0, NULL, NULL);
	if (error != CL_SUCCESS) {
		return error;
	}
	
	return pl_run_kernel_forward_geodesic_fixed_distance(
            pl_ctx, fwd_kernel, pl_buf, xy_out, pl_ell, distance);
}

int pl_forward_geodesic_fixed_angle(PLContext *pl_ctx, PLForwardGeodesicFixedAngleBuffer *pl_buf,
        double *xy_in, double *xy_out, PLSpheroid pl_ell, double angle) {
    double start[2] = { xy_in[0], xy_in[1] };
    cl_kernel fwd_kernel = NULL;
    if (_pl_spheroid_is_spherical(pl_ell)) {
        fwd_kernel = pl_find_kernel(pl_ctx, "pl_forward_geodesic_fixed_angle_s");
    } else {
        fwd_kernel = pl_find_kernel(pl_ctx, "pl_forward_geodesic_fixed_angle_e");
    }
    if (fwd_kernel == NULL) {
        return CL_INVALID_KERNEL_NAME;
    }
    
    return pl_run_kernel_forward_geodesic_fixed_angle(
            pl_ctx, fwd_kernel, pl_buf, start, xy_out, pl_ell, angle);
}

PLInverseGeodesicBuffer *pl_load_inverse_geodesic_data(PLContext *pl_ctx,
        const double *xy1_in, size_t xy1_count, cl_bool xy1_copy,
        const double *xy2_in, size_t xy2_count, 
        cl_int *outError)
{
	cl_int error = CL_SUCCESS;
	
	PLInverseGeodesicBuffer *pl_buf;
	if ((pl_buf = malloc(sizeof(PLInverseGeodesicBuffer))) == NULL) {
		return NULL;
	}
	pl_buf->xy1_in = NULL;
	pl_buf->xy2_in = NULL;
	pl_buf->xy1_count = xy1_count;
	pl_buf->xy2_count = xy2_count;
	pl_buf->dist_out = NULL;
	
	double *xy2_pad_in = NULL;
	int xy2_pad_count = ck_padding(xy2_count, PL_DOUBLE_VECTOR_SIZE);
	
	if ((xy2_pad_in = malloc(xy2_pad_count * sizeof(double) * 2)) == NULL) {
		error = CL_OUT_OF_HOST_MEMORY;
		goto cleanup;
	}
	
	_pl_copy_pad(xy2_pad_in, xy2_pad_count * 2, xy2_in, xy2_count * 2);
	
	pl_buf->xy1_in = clCreateBuffer(pl_ctx->ctx, CL_MEM_READ_ONLY | (xy1_copy ? CL_MEM_COPY_HOST_PTR : CL_MEM_USE_HOST_PTR), 
									sizeof(cl_double) * xy1_count * 2, (void *)xy1_in, &error);
	if (error != CL_SUCCESS) {
		goto cleanup;
	}
	
	pl_buf->xy2_in = clCreateBuffer(pl_ctx->ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
									sizeof(cl_double) * xy2_pad_count * 2, xy2_pad_in, &error);
	if (error != CL_SUCCESS) {
		goto cleanup;
	}
	
	pl_buf->dist_out = clCreateBuffer(pl_ctx->ctx, CL_MEM_READ_WRITE, 
                                      sizeof(cl_double) * xy1_count * xy2_pad_count * 2, NULL, &error);
	if (error != CL_SUCCESS) {
		goto cleanup;
	}
	
cleanup:
	if (xy2_pad_in)
		free(xy2_pad_in);
	if (error != CL_SUCCESS) {
		if (pl_buf->xy1_in)
			clReleaseMemObject(pl_buf->xy1_in);
		if (pl_buf->xy2_in)
			clReleaseMemObject(pl_buf->xy2_in);
		if (pl_buf->dist_out)
			clReleaseMemObject(pl_buf->dist_out);
		free(pl_buf);
		if (outError != NULL)
			*outError = error;
		return NULL;
	}
	
	return pl_buf;	
}

void pl_unload_inverse_geodesic_data(PLInverseGeodesicBuffer *pl_buf) {
	clReleaseMemObject(pl_buf->xy1_in);
	clReleaseMemObject(pl_buf->xy2_in);
	clReleaseMemObject(pl_buf->dist_out);
	free(pl_buf);
}

cl_int pl_inverse_geodesic(PLContext *pl_ctx, PLInverseGeodesicBuffer *pl_buf, double *dist_out,
        PLSpheroid pl_ell, double scale) {
	cl_kernel inv_kernel = NULL;
	
	if (_pl_spheroid_is_spherical(pl_ell)) {
		inv_kernel = pl_find_kernel(pl_ctx, "pl_inverse_geodesic_s");
	} else {
		inv_kernel = pl_find_kernel(pl_ctx, "pl_inverse_geodesic_e");
	}
	if (inv_kernel == NULL) {
		return CL_INVALID_KERNEL_NAME;
	}
	
	return pl_run_kernel_inverse_geodesic(pl_ctx, inv_kernel, pl_buf, dist_out, pl_ell, scale);
}
