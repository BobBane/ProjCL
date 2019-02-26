//
//  projcl_jni.c
//
// Created by Bob Bane on 2019-01-16
//


#include <projcl/projcl.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#define _GNU_SOURCE     /* To get pthread_getattr_np() declaration */
#include <pthread.h>

#include "projcl_jni.h"

#define JNI_TIMERS
#ifdef JNI_TIMERS
struct pltimer {
  char *label;
  double times;
  int calls;
};

struct pltimer pltimers[] = {
  {"load_projection_data", 0.0, 0},
  {"project and unload", 0.0, 0}
};

#define RECORD_TIME(n) (pltimers[n].calls++, pltimers[n].times += ctx->last_time)
#else
#define RECORD_TIME(n)
#endif

// Global error buffer
char errbuf[255];

// Little anaphoric helper for returning Java strings
#define JS(str) ((*env)->NewStringUTF(env, str))

// And common error handler for pthread stuff below
#define handle_error_en(en, msg) \
  do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)

/*
void print_stack_stuff()
{
  // Check the stack - how much do we really have left?
  // First, which one are we?
  pthread_t ptid = pthread_self();
  pthread_attr_t pat;
  int result = pthread_getattr_np (ptid, &pat);
  if (result != 0)
    handle_error(result, "pthread_getattr_np");

  size_t stack_size, guard_size;
  void *stack_addr;
  result = pthread_attr_getguardsize(&pat, &guard_size);
  if (result != 0)
    handle_error(result, "pthread_attr_getguardsize");
  result = pthread_attr_getstack(&pat, &stack_addr, &stack_size);
  if (result != 0)
    handle_error(result, "pthread_attr_getstack");
  printf("Stack size (bytes) = 0x%x\n", stack_size);
  printf("Guard size (bytes) = 0x%x\n", guard_size);

  result = pthread_attr_destroy(&pat);
  if (result != 0)
    handle_error(result, "pthread_attr_destroy");
}
*/

// Our OpenCL context handle
PLContext *ctx = NULL;

// The table of known proj.4 projections and a map to ProjCL's types.
// You would think this would be a service in ProjCL...
// Yes, ProjCL does some stuff outside of proj.4...
typedef struct proj4ToProjCL_s {
  char *p4name;
  PLProjection projection;
  int module;
} proj4ToProjCL;

proj4ToProjCL projectionMap[] = {
  {"stere",
   PL_PROJECT_OBLIQUE_STEREOGRAPHIC,
   PL_MODULE_OBLIQUE_STEREOGRAPHIC},
  // We map both proj.4 stereographic projections to
  // ProjCL oblique stereographic.  Close enough?
  {"sterea",
   PL_PROJECT_OBLIQUE_STEREOGRAPHIC,
   PL_MODULE_OBLIQUE_STEREOGRAPHIC},
  {"lcc",
   PL_PROJECT_LAMBERT_CONFORMAL_CONIC,
   PL_MODULE_LAMBERT_CONFORMAL_CONIC},
  {"laea",
   PL_PROJECT_LAMBERT_AZIMUTHAL_EQUAL_AREA,
   PL_MODULE_LAMBERT_AZIMUTHAL_EQUAL_AREA},
  {"tmerc",
   PL_PROJECT_TRANSVERSE_MERCATOR,
   PL_MODULE_TRANSVERSE_MERCATOR},
  {NULL, 0, 0}
};

// And code to search it

 proj4ToProjCL *findPLProjection (const char *p4name)
{
  proj4ToProjCL *pm;
  for(pm = projectionMap;
      pm->p4name != NULL;
      pm++) {
    if (!strcmp(p4name, pm->p4name)) {
      return pm;
    }
  }
  return NULL;
}

// And our memory of cmopiled modules
int COMPILED_MODULES = 0;
#define isCompiled(module) (COMPILED_MODULES & (module))
#define markCompiled(module) (COMPILED_MODULES |= (module))

// Helper to ensure a module is compiled, loaded, and ready to go
int compile_module_code(PLContext *ctx, unsigned int module, char *name) {
  PLCode *code = NULL;
  cl_int error = CL_SUCCESS;
  code = pl_compile_code(ctx, module, &error);
  if (code == NULL) {
    sprintf(errbuf, "projcl_jni - compile_module(%s) failed with code %d",
	    name, error);
    return 1;
  }
  error = pl_load_code(ctx, code);
  pl_release_code(code);
  if (error != CL_SUCCESS)
    printf("Failed to load code: %d\n", error);
  return error;
}

JNIEXPORT jstring JNICALL Java_gov_nasa_gsfc_drl_opencl_projcl_ProjCL_init
(JNIEnv *env, jclass jc)
{
  if (ctx != NULL) {
    return JS("projcl_jni - already initialized");
  }
  cl_int error = CL_SUCCESS;

  ctx = pl_context_init(CL_DEVICE_TYPE_GPU, &error);
  if (ctx == NULL) {
    sprintf(errbuf, "projcl_jni - failed to initialize context: %d", error);
    return JS(errbuf);
  }
  // Check the stack limits before diving into the compiler
  // Would you believe at this point, we have 0x101000 stack?
  // print_stack_stuff();

  // for now, ensure that our fave projections are compiled
  // This will turn into something lazy eventually...
  
  //if(compile_module_code(ctx, PL_MODULE_OBLIQUE_STEREOGRAPHIC, "Oblique Stereographic"))
  // return JS(errbuf);
  
  return NULL;
}


JNIEXPORT void JNICALL Java_gov_nasa_gsfc_drl_opencl_projcl_ProjCL_quit
(JNIEnv *env, jclass jc)
{
  if (ctx != NULL) {
    pl_unload_code(ctx);
    pl_context_free(ctx);
    ctx = NULL;
  }
#ifdef JNI_TIMERS
  int i;
  for(i=0; i< sizeof(pltimers)/sizeof(struct pltimer); i++)
    fprintf(stderr, "%s: CALLS: %d MILLIS: %f\n",
	    pltimers[i].label,
	    pltimers[i].calls,
	    pltimers[i].times * 1000.0);
#endif
}


// Exported projection lookup function    
JNIEXPORT jint JNICALL
Java_gov_nasa_gsfc_drl_opencl_projcl_ProjCL_getProjection
(JNIEnv *env, jclass jc,
 jstring p4name
 )
{
  jint result = 0;

  const char *cstring = (*env)->GetStringUTFChars(env, p4name, 0);

  proj4ToProjCL *plp = findPLProjection(cstring);

  if (plp != NULL)
    result = plp->projection;

  (*env)->ReleaseStringUTFChars(env, p4name, cstring);
  return result;
}

// Helper function to check array lengths
int check_array_length(JNIEnv *env, jdoubleArray arr, jsize count)
{
  jsize c2 = (*env)->GetArrayLength(env, arr);
  if(c2 != count) {
    sprintf(errbuf, "projcl_jni: array sizes do not match - %d versus %d",
	    count, c2);
    return 1;
  }
  return 0;
}

JNIEXPORT jstring JNICALL
Java_gov_nasa_gsfc_drl_opencl_projcl_ProjCL_projectForward
(JNIEnv *env, jclass jc,
 jstring p4name,
 jdoubleArray xyin,
 jdoubleArray xyout,
 jdouble latOrigin, jdouble lonOrigin,
 double standardParallel1,
 double standardParallel2
 )
{
  
  const char *cstring = (*env)->GetStringUTFChars(env, p4name, 0);

  proj4ToProjCL *plp =  findPLProjection(cstring);

  sprintf(errbuf, "projcl_jni: no ProjCL projection for %s", cstring);

  (*env)->ReleaseStringUTFChars(env, p4name, cstring);
  
  if (plp == 0)
    return JS(errbuf);
  
  // Found projection - make sure it's compiled
  if (!isCompiled(plp->module)) {
    if (compile_module_code(ctx, plp->module, plp->p4name))
      return JS(errbuf);
    markCompiled(plp->module);
  }
  
  cl_int error = CL_SUCCESS;

  // Grab the incoming data and push it into the GPU
  // Geet the length of the xin array
  jsize count = (*env)->GetArrayLength(env, xyin);
  // Check the rest of the arrays - they need to match
  if (check_array_length(env, xyout, count))
    return JS(errbuf);

  jdouble *XYin = (*env)->GetDoubleArrayElements(env, xyin, 0);
  
  PLProjectionBuffer *inbuf = NULL;
  inbuf  = pl_load_projection_data(ctx, XYin, count/2, JNI_TRUE, &error);
  RECORD_TIME(0);

  // can now release the xyin pointers
  (*env)->ReleaseDoubleArrayElements(env, xyin, XYin, 0);

  // And check for errors
  if (error != CL_SUCCESS) {
    sprintf(errbuf, "projcl_jni: pl_load_projection_data returned %d",
	    error);
    if(inbuf != NULL)
      pl_unload_projection_data(inbuf);
    return JS(errbuf);
  }

  // Set up the transformation
  PLProjectionParams *params = pl_params_init();

  // Default values striaght out of ProjCL test suite
  pl_params_set_false_northing(params, 0.0);
  pl_params_set_false_easting(params, 0.0);
  pl_params_set_scale(params, 1.0);
  pl_params_set_spheroid(params, PL_SPHEROID_WGS_84);
  pl_params_set_latitude_of_origin(params, latOrigin);
  pl_params_set_longitude_of_origin(params, lonOrigin);
  pl_params_set_standard_parallels(params,standardParallel1, standardParallel2);

  // Grab the output array and run the GPU
  jdouble *XYout = (*env)->GetDoubleArrayElements(env, xyout, 0);
  error = pl_project_points_forward(ctx, plp->projection,
				      params, inbuf, XYout);
  RECORD_TIME(1);

  // can now release the xyout pointers
  (*env)->ReleaseDoubleArrayElements(env, xyout, XYout, 0);
  // And check for errors
  if (error != CL_SUCCESS) {
    sprintf(errbuf, "projcl_jni: pl_project_points_forward_2 returned %d",
	    error);
    pl_unload_projection_data(inbuf);
    return JS(errbuf);
  }

  // We are done - clean up and return
  pl_unload_projection_data(inbuf);

  return NULL;
}

// This one takes two arrays and writes results into two arrays
// (it bottoms out in two calls to clEnqueueWriteBufferRect)
JNIEXPORT jstring JNICALL
Java_gov_nasa_gsfc_drl_opencl_projcl_ProjCL_projectForward2
(JNIEnv *env, jclass jc,
 jstring p4name,
 jdoubleArray xin, jdoubleArray yin,
 jdoubleArray xout, jdoubleArray yout,
 jdouble latOrigin, jdouble lonOrigin,
 double standardParallel1,
 double standardParallel2
 )
{
  
  const char *cstring = (*env)->GetStringUTFChars(env, p4name, 0);

  proj4ToProjCL *plp =  findPLProjection(cstring);

  sprintf(errbuf, "projcl_jni: no ProjCL projection for %s", cstring);

  (*env)->ReleaseStringUTFChars(env, p4name, cstring);
  
  if (plp == 0)
    return JS(errbuf);
  
  // Found projection - make sure it's compiled
  if (!isCompiled(plp->module)) {
    if (compile_module_code(ctx, plp->module, plp->p4name))
      return JS(errbuf);
    markCompiled(plp->module);
  }
  
  cl_int error = CL_SUCCESS;

  // Grab the incoming data and push it into the GPU
  // Geet the length of the xin array
  jsize count = (*env)->GetArrayLength(env, xin);
  // Check the rest of the arrays - they need to match
  if (check_array_length(env, yin, count))
    return JS(errbuf);
  if (check_array_length(env, xout, count))
    return JS(errbuf);
  if (check_array_length(env, yout, count))
    return JS(errbuf);

  jdouble *Xin = (*env)->GetDoubleArrayElements(env, xin, 0);
  jdouble *Yin = (*env)->GetDoubleArrayElements(env, yin, 0);
  
  PLProjectionBuffer *inbuf = NULL;
  inbuf  = pl_load_projection_data_2(ctx, Xin, Yin, count, &error);
  RECORD_TIME(0);

  // can now release the xin/yin pointers
  (*env)->ReleaseDoubleArrayElements(env, xin, Xin, 0);
  (*env)->ReleaseDoubleArrayElements(env, yin, Yin, 0);

  // And check for errors
  if (error != CL_SUCCESS) {
    sprintf(errbuf, "projcl_jni: pl_load_projection_data returned %d",
	    error);
    if(inbuf != NULL)
      pl_unload_projection_data(inbuf);
    return JS(errbuf);
  }

  // Set up the transformation
  PLProjectionParams *params = pl_params_init();

  // Default values striaght out of ProjCL test suite
  pl_params_set_false_northing(params, 0.0);
  pl_params_set_false_easting(params, 0.0);
  pl_params_set_scale(params, 1.0);
  pl_params_set_spheroid(params, PL_SPHEROID_WGS_84);
  pl_params_set_latitude_of_origin(params, latOrigin);
  pl_params_set_longitude_of_origin(params, lonOrigin);
  pl_params_set_standard_parallels(params,standardParallel1, standardParallel2);

  // Grab the output arrays and run the GPU
  jdouble *Xout = (*env)->GetDoubleArrayElements(env, xout, 0);
  jdouble *Yout = (*env)->GetDoubleArrayElements(env, yout, 0);
  error = pl_project_points_forward_2(ctx, plp->projection,
				      params, inbuf, Xout, Yout);
  RECORD_TIME(1);

  // can now release the xout/youy pointers
  (*env)->ReleaseDoubleArrayElements(env, xout, Xout, 0);
  (*env)->ReleaseDoubleArrayElements(env, yout, Yout, 0);
  // And check for errors
  if (error != CL_SUCCESS) {
    sprintf(errbuf, "projcl_jni: pl_project_points_forward_2 returned %d",
	    error);
    pl_unload_projection_data(inbuf);
    return JS(errbuf);
  }

  // We are done - clean up and return
  pl_unload_projection_data(inbuf);

  return NULL;
}
