/*
 *  projcl_types.h
 *  Magic Maps
 *
 *  Created by Evan Miller on 2/12/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

// This number is actually the number of entries processed by one
// kernel run, which is the same for both single and double
#define PL_FLOAT_VECTOR_SIZE 8
#define PL_DOUBLE_VECTOR_SIZE 8
#define ck_padding(n, size) (((n)+(size)-1)/(size)*(size))

typedef enum PLDatum {
    PL_DATUM_NONE = -1,
    PL_DATUM_WGS_84,
    PL_DATUM_WGS_72,
    PL_DATUM_ED_50,
    PL_DATUM_ED_79,
    PL_DATUM_ED_87,
    PL_DATUM_AUSTRIA_NS,
    PL_DATUM_BELGIUM_50,
    PL_DATUM_BERNE_1873,
    PL_DATUM_CH_1903,
    PL_DATUM_DANISH_GI_1934,
    PL_DATUM_NOUV_TRIG_DE_FRANCE_GREENWICH,
    PL_DATUM_NOUV_TRIG_DE_FRANCE_PARIS,
    PL_DATUM_POTSDAM,
    PL_DATUM_GGRS_87,
    PL_DATUM_HJORSEY_55,
    PL_DATUM_IRELAND_65,
    PL_DATUM_ITALY_1940,
    PL_DATUM_NOUV_TRIG_DE_LUX,
    PL_DATUM_NETHERLANDS_1921,
    PL_DATUM_OSGB_36,
    PL_DATUM_PORTUGAL_DLX,
    PL_DATUM_PORTUGAL_1973,
    PL_DATUM_RNB_72,
    PL_DATUM_RT_90,
    PL_DATUM_NAD_27,
    PL_DATUM_NAD_83,
    PL_DATUM_ETRS_89
} PLDatum;

typedef enum PLSpheroid {
	PL_SPHEROID_SPHERE,
	PL_SPHEROID_WGS_84,
    PL_SPHEROID_GRS_80,
    PL_SPHEROID_AIRY_1830,
    PL_SPHEROID_AIRY_1848,
    PL_SPHEROID_MODIFIED_AIRY,
    PL_SPHEROID_BESSEL_1841,
    PL_SPHEROID_CLARKE_1866,
    PL_SPHEROID_CLARKE_1880_RGS,
    PL_SPHEROID_GRS_1967_TRUNCATED,
    PL_SPHEROID_WGS_84_MAJOR_AUXILIARY_SPHERE,
    PL_SPHEROID_INTERNATIONAL_1924
} PLSpheroid;

typedef enum PLProjection_e {
    PL_PROJECT_ALBERS_EQUAL_AREA,
    PL_PROJECT_AMERICAN_POLYCONIC,
    PL_PROJECT_LAMBERT_CONFORMAL_CONIC,
    PL_PROJECT_LAMBERT_AZIMUTHAL_EQUAL_AREA,
    PL_PROJECT_MERCATOR,
    PL_PROJECT_OBLIQUE_STEREOGRAPHIC,
    PL_PROJECT_ROBINSON,
    PL_PROJECT_TRANSVERSE_MERCATOR,
    PL_PROJECT_WINKEL_TRIPEL,
    PL_PROJECT_STEREOGRAPHIC
} PLProjection;

struct pl_spheroid_info_s {
    double major_axis;
    double minor_axis;
};

#define PL_MODULE_DATUM                         (1 << 1)
#define PL_MODULE_GEODESIC                      (1 << 2)
#define PL_MODULE_WARP                          (1 << 3)
#define PL_MODULE_ALBERS_EQUAL_AREA             (1 << 4)
#define PL_MODULE_AMERICAN_POLYCONIC            (1 << 5)
#define PL_MODULE_LAMBERT_AZIMUTHAL_EQUAL_AREA  (1 << 6)
#define PL_MODULE_LAMBERT_CONFORMAL_CONIC       (1 << 7)
#define PL_MODULE_MERCATOR                      (1 << 8)
#define PL_MODULE_OBLIQUE_STEREOGRAPHIC         (1 << 9)
#define PL_MODULE_ROBINSON                      (1 << 10)
#define PL_MODULE_TRANSVERSE_MERCATOR           (1 << 11)
#define PL_MODULE_WINKEL_TRIPEL                 (1 << 12)
#define PL_MODULE_STEREOGRAPHIC                 (1 << 13)

#define PL_MODULE_PROJECTION                    (0xFFFF - 7)

typedef struct PLCode_s {
    cl_program   program;
	cl_uint      kernel_count;
} PLCode;

typedef struct PLSpheroidInfo_s {
	PLSpheroid   tag;
	double		major_axis;
	double		minor_axis;
	double		ecc;
	double		ecc2;
	double		one_ecc2;
	double		ec;
	double		inverse_flattening;
	double		en[8];
    double       apa[4];
    double       krueger_A;
	double		krueger_alpha[8];
	double		krueger_beta[8];
} PLSpheroidInfo;

typedef struct PLContext_s {
	cl_context		 ctx;
	cl_command_queue queue;
	cl_device_id     device_id;
	cl_uint		     kernel_count;
	cl_kernel	    *kernels;
    double           last_time;
    double           split_time;
} PLContext;

typedef struct PLProjectionBuffer_s {
	cl_mem      xy_in;
	cl_mem      xy_out;
	cl_uint	    count;
} PLProjectionBuffer;

typedef struct PLProjectionParams_s {
    double scale;
    double x0;
    double y0;
    double lon0;
    double lat0;
    double rlat1;
    double rlat2;

    PLSpheroid spheroid;
} PLProjectionParams;

typedef struct PLPointGridBuffer_s {
    cl_mem      grid;
    size_t      width;
    size_t      height;
} PLPointGridBuffer;

typedef struct PLDatumShiftBuffer_s {
    cl_mem      x_rw;
    cl_mem      y_rw;
    cl_mem      z_rw;
    cl_mem      xy_in;
    cl_mem      xy_out;
    cl_uint     count;
} PLDatumShiftBuffer;

typedef struct PLForwardGeodesicFixedDistanceBuffer_s {
	cl_mem	    xy_in;
	cl_uint     xy_count;
	cl_mem	    az_in;
	cl_uint     az_count;
	cl_mem	    phi_sincos;
	cl_mem      az_sincos;
	cl_mem      xy_out;
} PLForwardGeodesicFixedDistanceBuffer;

typedef struct PLForwardGeodesicFixedAngleBuffer_s {
	cl_mem	    dist_in;
	cl_uint     dist_count;
	cl_mem      xy_out;
} PLForwardGeodesicFixedAngleBuffer;

typedef struct PLInverseGeodesicBuffer_s {
	cl_mem      xy1_in;
	cl_uint     xy1_count;
	cl_mem      xy2_in;
	cl_uint     xy2_count;
	cl_mem      dist_out;
} PLInverseGeodesicBuffer;
