
#ifndef BLAS_INCLUDE
#define BLAS_INCLUDE

/* Data types specific to BLAS implementation */
typedef struct { float r, i; } fcomplex;
typedef struct { double r, i; } dcomplex;
typedef int blasbool;

#include "blasp.h"    /* Prototypes for all BLAS functions */

#define FALSE 0
#define TRUE  1

/* Macro functions */
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))

#endif
