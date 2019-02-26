set(PROJCL_HEADERS_DIR ${PROJECT_SOURCE_DIR}/include/projcl)

set(PROJCL_H
  ${PROJCL_HEADERS_DIR}/projcl.h
  ${PROJCL_HEADERS_DIR}/projcl_warp.h 
  ${PROJCL_HEADERS_DIR}/projcl_types.h 

)

set(PROJCL_C
  projcl_datum.c
  projcl_geodesic.c
  projcl_kernel.c
  projcl_kernel.h
  projcl_kernel_sources.c
  projcl_kernel_sources.h
  projcl_kernel_sources_text.h
  projcl_params.c
  projcl_project.c
  projcl_run.c
  projcl_run.h
  projcl_spheroid.c
  projcl_spheroid.h
  projcl_util.c
  projcl_util.h
  projcl_warp.c

  projcl.c
)

# Standard include directory of PROJCL library
include_directories(../include)

set(PROJCL_SOURCES
  ${PROJCL_H}
  ${PROJCL_C})

add_library(${PROJCL_LIBRARY} SHARED ${PROJCL_SOURCES})
target_link_libraries(${PROJCL_LIBRARY} ${PROJCL_DEPENDENT_LIBRARIES})

# projcl_kernel_sources_text.h is special - it is created from the
# kernel/*.opencl files next door, like this:

FILE(GLOB KernelSources ../kernel/*.opencl)
# message(STATUS "Kernel sources are ${KernelSources}")
add_custom_command( OUTPUT projcl_kernel_sources_text.h
		    DEPENDS ${KernelSources}
		    COMMAND ../kernel/make_projcl_kernel_sources_text.sh
		    )
