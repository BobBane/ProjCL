//
//  projcl_kernel_sources.c
//  Magic Maps
//
//  Created by Bob Bane, NASA Direct Readout Laboratory
//  Copyright 2019
//

#include "projcl_kernel_sources.h"
#include <stdio.h>

typedef struct PLKernelSources_s {
    char *filename;
    char *data;
} PLKernelSources;

#include "projcl_kernel_sources_text.h"

const char *find_kernel_source(const char *filename) {
  PLKernelSources *pks = _pl_kernelsources;
  for (pks = _pl_kernelsources;
       pks->filename != NULL;
       pks++) {
    if (!strcmp(filename, pks->filename)) {
      // fprintf(stderr, "Found kernel %s\n", pks->filename);
      return pks->data;
    }
  }
  return NULL;
}

void dump_kernel_sources() {
  PLKernelSources *pks = _pl_kernelsources;
  for (pks = _pl_kernelsources;
       pks->filename != NULL;
       pks++) {
    fprintf(stderr, "Kernel %s\n", pks->filename);
    fprintf(stderr, "%s\n\n", pks->data);
  }
}
