#
# Trying to go with "modern" cmake structure for JNI
# (the GENERATE_NATIVE_HEADERS below)
#
cmake_minimum_required(VERSION 3.11.4)


find_package(Java REQUIRED COMPONENTS Development )
find_package(JNI REQUIRED)
include(UseJava)

project(projcl_jni)
set(CMAKE_JNI_TARGET TRUE)

set(PROJCL_JAVA_STUBS gov/nasa/gsfc/drl/opencl/projcl/ProjCL.java)
add_jar ( projcl_jni
	  SOURCES ${PROJCL_JAVA_STUBS}
	  GENERATE_NATIVE_HEADERS projcl_jni_native
)
add_library(projcl_jni_c SHARED projcl_jni.c)
target_link_libraries(projcl_jni_c PRIVATE projcl_jni_native ${PROJCL_LIBRARY})
include_directories(../include ${JNI_INCLUDE_DIRS})
