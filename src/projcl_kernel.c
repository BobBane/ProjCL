
#include <projcl/projcl.h>
#include "projcl_kernel.h"

#include <string.h>
#include <stdio.h>

cl_kernel pl_find_kernel(PLContext *pl_ctx, const char *requested_name) {
	char buf[128];
	size_t len;
	cl_int error;
	cl_kernel kernel;
	int i;
	for (i=0; i<pl_ctx->kernel_count; i++) {
		kernel = pl_ctx->kernels[i];
		error = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 
								sizeof(buf), buf, &len);
		if (strncmp(requested_name, buf, len) == 0) {
			return kernel;
		}
	}
	fprintf(stderr, "No kernel named %s\n", requested_name);
	fprintf(stderr, "KNOWN KERNELS: %d\n", pl_ctx->kernel_count);
	for (i=0; i<pl_ctx->kernel_count; i++) {
		kernel = pl_ctx->kernels[i];
		error = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 
								sizeof(buf), buf, &len);
		fprintf(stderr, "%3d: %*s\n", i, len, buf);
	}
	return NULL;
}

cl_int pl_set_kernel_arg_mem(PLContext *ctx, cl_kernel kernel, int argc, cl_mem buffer) {
    return clSetKernelArg(kernel, argc, sizeof(cl_mem), &buffer);
}

cl_int pl_set_kernel_arg_double16(PLContext *ctx, cl_kernel kernel, int argc, double value[16]) {
    double value_f[16] = { 
        value[0], value[1], value[2], value[3],
        value[4], value[5], value[6], value[7],
        value[8], value[9], value[10], value[11],
        value[12], value[13], value[14], value[15] };
    return clSetKernelArg(kernel, argc, sizeof(cl_double16), value_f);
}

cl_int pl_set_kernel_arg_double8(PLContext *ctx, cl_kernel kernel, int argc, double value[8]) {
    double value_f[8] = { value[0], value[1], value[2], value[3], value[4], value[5], value[6], value[7] };
    return clSetKernelArg(kernel, argc, sizeof(cl_double8), value_f);
}

cl_int pl_set_kernel_arg_double4(PLContext *ctx, cl_kernel kernel, int argc, double value[4]) {
    double value_f[4] = { value[0], value[1], value[2], value[3] };
    return clSetKernelArg(kernel, argc, sizeof(cl_double4), value_f);
}

cl_int pl_set_kernel_arg_double2(PLContext *ctx, cl_kernel kernel, int argc, double value[2]) {
    double value_f[2] = { value[0], value[1] };
    return clSetKernelArg(kernel, argc, sizeof(cl_double2), value_f);
}

cl_int pl_set_kernel_arg_double(PLContext *ctx, cl_kernel kernel, int argc, double value) {
    double value_f = value;
    return clSetKernelArg(kernel, argc, sizeof(cl_double), &value_f);
}

cl_int pl_set_kernel_arg_uint(PLContext *ctx, cl_kernel kernel, int argc, cl_uint value) {
    return clSetKernelArg(kernel, argc, sizeof(cl_uint), &value);
}
