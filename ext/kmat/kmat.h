#ifndef KMAT_H
#define KMAT_H 1

#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <ruby.h>
#include "lapack_headers/blas.h"
#include "lapack_headers/lapacke.h"
#include "lapack_headers/lapacke_utils.h"

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif
#define M_2PI 6.28318530717958647692
#define COMPLEX double _Complex
#ifdef CMPLX /* C11 macro is available */
#  define cpack(x,y) CMPLX(x,y)
#elif (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)) && !defined(__INTEL_COMPILER)
#  define cpack(x,y) __builtin_complex((double) (x), (double) (y))
#else
#  define cpack(x,y) ((double)(x)+I*(double)(y))
#endif

#ifndef __GNUC__
#  define  __attribute__(x)
#endif

typedef enum {
	VT_DOUBLE,
	VT_COMPLEX,
	VT_INT,
	VT_BOOL,
	VT_VALUE,
	VT_END
} VTYPE;
typedef enum {
	ST_FULL,
	ST_SSUB,
	ST_RSUB,
	ST_END
} STYPE;
typedef struct {
	union {
		void *body;
		double *dbody;
		COMPLEX *zbody;
		int *ibody;
		bool *bbody;
		VALUE *vbody;
		void **pbody;
		double **dpbody;
		COMPLEX **zpbody;
		int **ipbody;
		bool **bpbody;
		VALUE **vpbody;
	};
	int ld, m, n;
	VTYPE vtype;
	STYPE stype;
	bool trans, may_have_sub;
	VALUE parent;
} SMAT;
typedef union {
	double d;
	COMPLEX z;
	int i;
	bool b;
	VALUE v;
} ELEMENT;
typedef struct {
	union {
		void *body;
		double *d;
		COMPLEX *z;
		int *i;
		bool *b;
		VALUE *v;
	};
	int ld;
	bool need_to_free;
} LAWORK;

#include "id_sym.h"
#include "global_variables.h"
extern const rb_data_type_t km_mat_data_type;
extern VALUE kmgv_mat_random;

#include "km_util.h"

#include "auto_collected.h"
#include "elementwise_function.h"

#endif /* KMAT_H */
