#include "../kmat.h"

// sort rows or columns by ascending order
static int
km_int_comp(const void *a, const void *b)
{
	return *(int *)a-*(int *)b;
}
static int
km_bool_comp(const void *a, const void *b)
{
	if ( *(bool *)a ) {
		if ( *(bool *)b ) {
			return 0;
		} else {
			return 1;
		}
	} else {
		if ( *(bool *)b ) {
			return -1;
		} else {
			return 0;
		}
	}
}
static int
km_value_comp(const void *a, const void *b)
{
	return NUM2INT(rb_funcall(*(VALUE *)a, id_op_comp, 1, *(VALUE *)b));
}
static void
km_sort_seg(size_t size_seg, size_t ld, size_t num_seg, void *body, VTYPE vt)
{
	if ( vt == VT_DOUBLE ) {
		for ( size_t i=0; i<num_seg; i++ ) {
			int info, sseg=s2i(size_seg); char id[] = "I";
			dlasrt_(id, &sseg, ((double *)body)+(i*ld), &info);
		}
	} else if ( vt == VT_INT ) {
		for ( size_t i=0; i<num_seg; i++ ) {
			qsort(((int *)body)+(i*ld), size_seg, sizeof(int), km_int_comp);
		}
	} else if ( vt == VT_BOOL ) {
		for ( size_t i=0; i<num_seg; i++ ) {
			qsort(((bool *)body)+(i*ld), size_seg, sizeof(bool), km_bool_comp);
		}
	} else {
		for ( size_t i=0; i<num_seg; i++ ) {
			qsort(((VALUE *)body)+(i*ld), size_seg, sizeof(VALUE), km_value_comp);
		}
	}
}
#define SKEWER(type) type *b = (type *)body; type *b2 = ZALLOC_N(type, size_sk*num_sk); \
	for ( size_t i=0; i<size_sk; i++ ) { \
		for ( size_t j=0; j<num_sk; j++ ) { \
			b2[i+j*size_sk] = b[j+i*ld]; \
		} \
	} \
	km_sort_seg(size_sk, size_sk, num_sk, b2, vt); \
	for ( size_t i=0; i<size_sk; i++ ) { \
		for ( size_t j=0; j<num_sk; j++ ) { \
			b[j+i*ld] = b2[i+j*size_sk]; \
		} \
	} \
	ruby_xfree(b2)

static void
km_sort_skewer(size_t size_sk, size_t ld, size_t num_sk, void *body, VTYPE vt)
{
	if ( vt == VT_DOUBLE ) {
		SKEWER(double);
	} else if ( vt == VT_INT ) {
		SKEWER(int);
	} else {
		SKEWER(VALUE);
	}
}
static void
km_smat_sort(SMAT *smat, VALUE asym)
{
	if ( rb_respond_to(asym, id_to_sym) ) {
		asym = rb_funcall(asym, id_to_sym, 0);
	}
	if ( asym == sym_col || asym == sym_c || asym == sym_column ) {
		if ( smat->trans ) {
			km_sort_skewer(smat->m, smat->ld, smat->n, smat->body, smat->vtype);
		} else {
			km_sort_seg(smat->m, smat->ld, smat->n, smat->body, smat->vtype);
		}
	} else if ( asym == sym_row || asym == sym_r ) {
		if ( smat->trans ) {
			km_sort_seg(smat->n, smat->ld, smat->m, smat->body, smat->vtype);
		} else {
			km_sort_skewer(smat->n, smat->ld, smat->m, smat->body, smat->vtype);
		}
	} else if ( asym == sym_none || asym == sym_n ) {
		return;
	} else {
		rb_raise(rb_eArgError, "unknown axis name");
	}
}
VALUE
kmm_mat_sort_destl(VALUE self, VALUE asym)
{
	km_check_frozen(self);
	VALUE foo=self;
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype != VT_COMPLEX ) {
		if ( smat->stype == ST_RSUB ) {
			foo = rb_obj_dup(self);
			smat = km_mat2smat(foo);
		}
	} else {
		foo = kmm_mat_to_omat(self);
		smat = km_mat2smat(foo);
	}
	km_smat_sort(smat, asym);
	if ( foo != self ) {
		kmm_mat_copy_from(self, foo);
	}
	return self;
}
VALUE
kmm_mat_sort(VALUE self, VALUE asym)
{
	SMAT *smat = km_mat2smat(self);
	VALUE ret;
	if ( smat->vtype != VT_COMPLEX ) {
		ret = rb_obj_dup(self);
	} else {
		ret = kmm_mat_to_omat(self);
	}
	km_smat_sort(km_mat2smat(ret), asym);
	if ( smat->vtype == VT_COMPLEX ) {
		return kmm_mat_to_cmat(ret);
	} else {
		return ret;
	}
}

// sort rows or columns by descending order
static int
km_int_rcomp(const void *a, const void *b)
{
	return *(int *)b-*(int *)a;
}
static int
km_bool_rcomp(const void *a, const void *b)
{
	if ( *(bool *)a ) {
		if ( *(bool *)b ) {
			return 0;
		} else {
			return -1;
		}
	} else {
		if ( *(bool *)b ) {
			return 1;
		} else {
			return 0;
		}
	}
}
static int
km_value_rcomp(const void *a, const void *b)
{
	return NUM2INT(rb_funcall(*(VALUE *)b, id_op_comp, 1, *(VALUE *)a));
}
static void
km_rsort_seg(size_t size_seg, size_t ld, size_t num_seg, void *body, VTYPE vt)
{
	if ( vt == VT_DOUBLE ) {
		for ( size_t i=0; i<num_seg; i++ ) {
			int info, sseg=s2i(size_seg); char id[] = "D";
			dlasrt_(id, &sseg, ((double *)body)+(i*ld), &info);
		}
	} else if ( vt == VT_INT ) {
		for ( size_t i=0; i<num_seg; i++ ) {
			qsort(((int *)body)+(i*ld), size_seg, sizeof(int), km_int_rcomp);
		}
	} else if ( vt == VT_BOOL ) {
		for ( size_t i=0; i<num_seg; i++ ) {
			qsort(((bool *)body)+(i*ld), size_seg, sizeof(bool), km_bool_rcomp);
		}
	} else {
		for ( size_t i=0; i<num_seg; i++ ) {
			qsort(((VALUE *)body)+(i*ld), size_seg, sizeof(VALUE), km_value_rcomp);
		}
	}
}
#define rSKEWER(type) type *b = (type *)body; type *b2 = ZALLOC_N(type, size_sk*num_sk); \
	for ( size_t i=0; i<size_sk; i++ ) { \
		for ( size_t j=0; j<num_sk; j++ ) { \
			b2[i+j*size_sk] = b[j+i*ld]; \
		} \
	} \
	km_rsort_seg(size_sk, size_sk, num_sk, b2, vt); \
	for ( size_t i=0; i<size_sk; i++ ) { \
		for ( size_t j=0; j<num_sk; j++ ) { \
			b[j+i*ld] = b2[i+j*size_sk]; \
		} \
	} \
	ruby_xfree(b2)

static void
km_rsort_skewer(size_t size_sk, size_t ld, size_t num_sk, void *body, VTYPE vt)
{
	if ( vt == VT_DOUBLE ) {
		rSKEWER(double);
	} else if ( vt == VT_INT ) {
		rSKEWER(int);
	} else {
		rSKEWER(VALUE);
	}
}
static void
km_smat_rsort(SMAT *smat, VALUE asym)
{
	if ( rb_respond_to(asym, id_to_sym) ) {
		asym = rb_funcall(asym, id_to_sym, 0);
	}
	if ( asym == sym_col || asym == sym_c || asym == sym_column ) {
		if ( smat->trans ) {
			km_rsort_skewer(smat->m, smat->ld, smat->n, smat->body, smat->vtype);
		} else {
			km_rsort_seg(smat->m, smat->ld, smat->n, smat->body, smat->vtype);
		}
	} else if ( asym == sym_row || asym == sym_r ) {
		if ( smat->trans ) {
			km_rsort_seg(smat->n, smat->ld, smat->m, smat->body, smat->vtype);
		} else {
			km_rsort_skewer(smat->n, smat->ld, smat->m, smat->body, smat->vtype);
		}
	} else if ( asym == sym_none || asym == sym_n ) {
		return;
	} else {
		rb_raise(rb_eArgError, "unknown axis name");
	}
}
VALUE
kmm_mat_rsort_destl(VALUE self, VALUE asym)
{
	km_check_frozen(self);
	VALUE foo=self;
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype != VT_COMPLEX ) {
		if ( smat->stype == ST_RSUB ) {
			foo = rb_obj_dup(self);
			smat = km_mat2smat(foo);
		}
	} else {
		foo = kmm_mat_to_omat(self);
		smat = km_mat2smat(foo);
	}
	km_smat_rsort(smat, asym);
	if ( foo != self ) {
		kmm_mat_copy_from(self, foo);
	}
	return self;
}
VALUE
kmm_mat_rsort(VALUE self, VALUE asym)
{
	SMAT *smat = km_mat2smat(self);
	VALUE ret;
	if ( smat->vtype != VT_COMPLEX ) {
		ret = rb_obj_dup(self);
	} else {
		ret = kmm_mat_to_omat(self);
	}
	km_smat_rsort(km_mat2smat(ret), asym);
	if ( smat->vtype == VT_COMPLEX ) {
		return kmm_mat_to_cmat(ret);
	} else if ( smat->vtype == VT_BOOL ) {
		return kmm_mat_to_bmat(ret);
	} else {
		return ret;
	}
}

// if `self' is a row vector, `self' is replaced by a vector which satisfy `A'[`p', `self'] is sorted row vector
// if `self' is a column vector, `self' is replaced by a vector which satisfy `A'[`self', `p'] is sorted column vector
static union {
	double *d;
	int *i;
	bool *b;
	VALUE *v;
} as_body;
static size_t as_step;
static int
as_dcomp(const void *a, const void *b)
{
	return (int)copysign(1.0, as_body.d[(*(size_t *)a)*as_step]-as_body.d[(*(size_t *)b)*as_step]);
}
static int
as_icomp(const void *a, const void *b)
{
	return as_body.i[(*(size_t *)a)*as_step] - as_body.i[(*(size_t *)b)*as_step];
}
static int
as_bcomp(const void *a, const void *b)
{
	if ( as_body.b[(*(size_t *)a)*as_step] ) {
		if ( as_body.b[(*(size_t *)b)*as_step] ) {
			return 0;
		} else {
			return 1;
		}
	} else {
		if ( as_body.b[(*(size_t *)b)*as_step] ) {
			return -1;
		} else {
			return 0;
		}
	}
}
static int
as_vcomp(const void *a, const void *b)
{
	return NUM2INT(rb_funcall(as_body.v[(*(size_t *)a)*as_step], id_op_comp, 1, as_body.v[(*(size_t *)b)*as_step]));
}
VALUE
kmm_mat_argsort_destl(VALUE self, VALUE va, VALUE vp)
{
	km_check_frozen(self);
	SMAT *si = km_mat2smat(self), *sa = km_mat2smat(va);
	if ( si->vtype != VT_INT || si->stype != ST_FULL || !VECTOR_P(si) ) {
		rb_raise(rb_eRuntimeError, "self must be an int full vector");
	}
	si->trans = false; si->ld = si->m; // this is available because si is a vector
	int *idx = si->ibody;
	size_t p=NUM2ZU(vp), offset;
	if ( si->m == sa->m && si->n == 1 ) {
		if ( sa->n <= p ) { rb_raise(rb_eIndexError, "p is out of range"); }
		for ( size_t i=0; i<si->m; i++ ) { idx[i] = s2i(i); }
		if ( sa->trans ) {
			offset = p; as_step = sa->ld;
		} else {
			offset = p*si->ld; as_step = 1;
		}
	} else if ( si->m == 1 && si->n == sa->n ) {
		if ( sa->m <= p ) { rb_raise(rb_eIndexError, "p is out of range"); }
		for ( size_t i=0; i<si->n; i++ ) { idx[i] = s2i(i); }
		if ( sa->trans ) {
			offset = p*si->ld; as_step = 1;
		} else {
			offset = p; as_step = sa->ld;
		}
	} else {
		rb_raise(km_eDim, "dimension mismatched");
	}
	if ( sa->vtype == VT_DOUBLE ) {
		as_body.d = sa->dbody; as_body.d += offset;
		qsort(idx, LENGTH(si), sizeof(int), as_dcomp);
	} else if ( sa->vtype == VT_INT ) {
		as_body.i = sa->ibody; as_body.i += offset;
		qsort(idx, LENGTH(si), sizeof(int), as_icomp);
	} else if ( sa->vtype == VT_BOOL ) {
		as_body.b = sa->bbody; as_body.b += offset;
		qsort(idx, LENGTH(si), sizeof(int), as_bcomp);
	} else if ( sa->vtype == VT_VALUE ) {
		as_body.v = sa->vbody; as_body.v += offset;
		qsort(idx, LENGTH(si), sizeof(int), as_vcomp);
	} else {
		SMAT *sa2 = km_mat2smat(kmm_mat_to_omat(va));
		as_body.v = sa2->vbody; as_body.v += offset;
		qsort(idx, LENGTH(si), sizeof(int), as_vcomp);
	}
	return self;
}
// select column-sort or row-sort by the argument Symbol
// the argument can be omitted only if `self' is a vector
VALUE
kmm_mat_argsort(int argc, VALUE *argv, VALUE va)
{
	rb_check_arity(argc, 0, 2);
	if ( argc == 0 ) {
		SMAT *sa = km_mat2smat(va);
		if ( sa->n == 1 ) {
			return kmm_mat_argsort_destl(km_Mat(sa->m, 1, VT_INT), va, INT2NUM(0));
		} else if ( sa->m == 1 ) {
			return kmm_mat_argsort_destl(km_Mat(1, sa->n, VT_INT), va, INT2NUM(0));
		} else {
			rb_raise(rb_eArgError, "argsort without arguments is available only for vectors");
		}
	} else if ( argc == 2 ) {
		VALUE asym = argv[0];
		if ( rb_respond_to(asym, id_to_sym) ) {
			asym = rb_funcall(asym, id_to_sym, 0);
		}
		SMAT *sa = km_mat2smat(va);
		if ( asym == sym_col || asym == sym_c || asym == sym_column ) {
			return kmm_mat_argsort_destl(km_Mat(sa->m, 1, VT_INT), va, argv[1]);
		} else if ( asym == sym_row || asym == sym_c ) {
			return kmm_mat_argsort_destl(km_Mat(1, sa->n, VT_INT), va, argv[1]);
		} else {
			rb_raise(rb_eArgError, "unknown axis name");
		}
	} else {
		rb_raise(rb_eArgError, "argsort with 1 argument is not available");
	}
}

// reversed argsort
static int
ras_dcomp(const void *a, const void *b)
{
	return (int)copysign(1.0, as_body.d[(*(size_t *)a)*as_step]-as_body.d[(*(size_t *)b)*as_step]);
}
static int
ras_icomp(const void *a, const void *b)
{
	return as_body.i[(*(size_t *)a)*as_step] - as_body.i[(*(size_t *)b)*as_step];
}
static int
ras_bcomp(const void *a, const void *b)
{
	if ( as_body.b[(*(size_t *)a)*as_step] ) {
		if ( as_body.b[(*(size_t *)b)*as_step] ) {
			return 0;
		} else {
			return -1;
		}
	} else {
		if ( as_body.b[(*(size_t *)b)*as_step] ) {
			return 1;
		} else {
			return 0;
		}
	}
}
static int
ras_vcomp(const void *a, const void *b)
{
	return NUM2INT(rb_funcall(as_body.v[(*(size_t *)a)*as_step], id_op_comp, 1, as_body.v[(*(size_t *)b)*as_step]));
}
VALUE
kmm_mat_rargsort_destl(VALUE self, VALUE va, VALUE vp)
{
	km_check_frozen(self);
	SMAT *si = km_mat2smat(self), *sa = km_mat2smat(va);
	if ( si->vtype != VT_INT || si->stype != ST_FULL || !VECTOR_P(si) ) {
		rb_raise(rb_eRuntimeError, "self must be an int full vector");
	}
	si->trans = false; si->ld = si->m; // this is available because si is a vector
	int *idx = si->ibody;
	size_t p=NUM2ZU(vp), offset;
	if ( si->m == sa->m && si->n == 1 ) {
		if ( sa->n <= p ) { rb_raise(rb_eIndexError, "p is out of range"); }
		for ( size_t i=0; i<si->m; i++ ) { idx[i] = s2i(i); }
		if ( sa->trans ) {
			offset = p; as_step = sa->ld;
		} else {
			offset = p*si->ld; as_step = 1;
		}
	} else if ( si->m == 1 && si->n == sa->n ) {
		if ( sa->m <= p ) { rb_raise(rb_eIndexError, "p is out of range"); }
		for ( size_t i=0; i<si->n; i++ ) { idx[i] = s2i(i); }
		if ( sa->trans ) {
			offset = p*si->ld; as_step = 1;
		} else {
			offset = p; as_step = sa->ld;
		}
	} else {
		rb_raise(km_eDim, "dimension mismatched");
	}
	if ( sa->vtype == VT_DOUBLE ) {
		as_body.d = sa->dbody; as_body.d += offset;
		qsort(idx, LENGTH(si), sizeof(int), ras_dcomp);
	} else if ( sa->vtype == VT_INT ) {
		as_body.i = sa->ibody; as_body.i += offset;
		qsort(idx, LENGTH(si), sizeof(int), ras_icomp);
	} else if ( sa->vtype == VT_BOOL ) {
		as_body.b = sa->bbody; as_body.b += offset;
		qsort(idx, LENGTH(si), sizeof(int), ras_bcomp);
	} else if ( sa->vtype == VT_VALUE ) {
		as_body.v = sa->vbody; as_body.v += offset;
		qsort(idx, LENGTH(si), sizeof(int), ras_vcomp);
	} else {
		SMAT *sa2 = km_mat2smat(kmm_mat_to_omat(va));
		as_body.v = sa2->vbody; as_body.v += offset;
		qsort(idx, LENGTH(si), sizeof(int), ras_vcomp);
	}
	return self;
}
VALUE
kmm_mat_rargsort(int argc, VALUE *argv, VALUE va)
{
	rb_check_arity(argc, 0, 2);
	if ( argc == 0 ) {
		SMAT *sa = km_mat2smat(va);
		if ( sa->n == 1 ) {
			return kmm_mat_rargsort_destl(km_Mat(sa->m, 1, VT_INT), va, INT2NUM(0));
		} else if ( sa->m == 1 ) {
			return kmm_mat_rargsort_destl(km_Mat(1, sa->n, VT_INT), va, INT2NUM(0));
		} else {
			rb_raise(rb_eArgError, "argsort without arguments is available only for vectors");
		}
	} else if ( argc == 2 ) {
		VALUE asym = argv[0];
		if ( rb_respond_to(asym, id_to_sym) ) {
			asym = rb_funcall(asym, id_to_sym, 0);
		}
		SMAT *sa = km_mat2smat(va);
		if ( asym == sym_col || asym == sym_c || asym == sym_column ) {
			return kmm_mat_rargsort_destl(km_Mat(sa->m, 1, VT_INT), va, argv[1]);
		} else if ( asym == sym_row || asym == sym_c ) {
			return kmm_mat_rargsort_destl(km_Mat(1, sa->n, VT_INT), va, argv[1]);
		} else {
			rb_raise(rb_eArgError, "unknown axis name");
		}
	} else {
		rb_raise(rb_eArgError, "argsort with 1 argument is not available");
	}
}

// store the maximum/minimum elements of `self' to the argument
// the axis is selcted by the shape of the argument
// if the argument is (1, 1)-matrix, return argument[0, 0]
// otherwise return the argument
static void
km_max_func_d(double *ea, void *data)
{
	double *r = (double *)data;
	if ( *r < *ea ) { *r = *ea; }
}
static void
km_max_func_i(int *ea, void *data)
{
	int *r = (int *)data;
	if ( *r < *ea ) { *r = *ea; }
}
static void
km_max_func_b(bool *ea, void *data)
{
	bool *r = (bool *)data;
	if ( *ea ) { *r = true; }
}
static void
km_max_func_v(VALUE *ea, void *data)
{
	VALUE *r = (VALUE *)data;
	if ( RTEST(rb_funcall(*r, id_op_lt, 1, *ea)) ) {
		*r = *ea;
	}
}
static void
km_min_func_d(double *ea, void *data)
{
	double *r = (double *)data;
	if ( *r > *ea ) { *r = *ea; }
}
static void
km_min_func_i(int *ea, void *data)
{
	int *r = (int *)data;
	if ( *r > *ea ) { *r = *ea; }
}
static void
km_min_func_b(bool *ea, void *data)
{
	bool *r = (bool *)data;
	if ( !*ea ) { *r = false; }
}
static void
km_min_func_v(VALUE *ea, void *data)
{
	VALUE *r = (VALUE *)data;
	if ( RTEST(rb_funcall(*r, id_op_gt, 1, *ea)) ) {
		*r = *ea;
	}
}
#define DEFINE_MM(mm, op, op_id, opb) \
static VALUE \
km_mat_##mm##_all(SMAT *sr, SMAT *sa, VALUE va) \
{ \
	if ( sr->vtype == VT_DOUBLE ) { \
		sr->dbody[0] = sa->dbody[0]; \
		km_smat_each_d(sa, km_##mm##_func_d, sr->body); \
		return rb_float_new(sr->dbody[0]); \
	} else if ( sr->vtype == VT_INT ) { \
		sr->ibody[0] = sa->ibody[0]; \
		km_smat_each_i(sa, km_##mm##_func_i, sr->body); \
		return INT2NUM(sr->ibody[0]); \
	} else if ( sr->vtype == VT_BOOL ) { \
		sr->bbody[0] = sa->bbody[0]; \
		km_smat_each_b(sa, km_##mm##_func_b, sr->body); \
		return TF2V(sr->bbody[0]); \
	} else if ( sr->vtype == VT_VALUE ) { \
		sr->vbody[0] = sa->vbody[0]; \
		km_smat_each_v(sa, km_##mm##_func_v, sr->body); \
		return sr->vbody[0]; \
	} else { \
		sa = km_mat2smat(kmm_mat_to_omat(va)); \
		VALUE ret=sa->vbody[0]; \
		km_smat_each_v(sa, km_##mm##_func_v, &ret); \
		if ( sr->vtype == VT_COMPLEX ) { \
			sr->zbody[0] = km_v2c(ret); \
		} else { \
			rb_raise(km_eInternal, "unexpected reach"); \
		} \
		return ret; \
	} \
} \
static VALUE \
km_mat_##mm##_col(VALUE self, SMAT *sr, SMAT *sa, VALUE va) \
{ \
	if ( sr->vtype == VT_DOUBLE ) { \
		for ( size_t j=0; j<sr->n; j++ ) { \
			sr->dbody[j] = ENTITY(sa, d, 0, j); \
			for ( size_t i=1; i<sa->m; i++ ) { \
				double tmp = ENTITY(sa, d, i, j); \
				if ( sr->dbody[j] op tmp ) { \
					sr->dbody[j] = tmp; \
				} \
			} \
		} \
	} else if ( sr->vtype == VT_INT ) { \
		for ( size_t j=0; j<sr->n; j++ ) { \
			sr->ibody[j] = ENTITY(sa, i, 0, j); \
			for ( size_t i=1; i<sa->m; i++ ) { \
				int tmp = ENTITY(sa, i, i, j); \
				if ( sr->ibody[j] op tmp ) { \
					sr->ibody[j] = tmp; \
				} \
			} \
		} \
	} else if ( sr->vtype == VT_BOOL ) { \
		for ( size_t j=0; j<sr->n; j++ ) { \
			sr->bbody[j] = ENTITY(sa, b, 0, j); \
			if ( !(opb sr->bbody[j]) ) { \
				for ( size_t i=1; i<sa->m; i++ ) { \
					bool tmp = ENTITY(sa, b, i, j); \
					if ( opb tmp ) { \
						sr->bbody[j] = tmp; \
					} \
				} \
			} \
		} \
	} else if ( sr->vtype == VT_VALUE ) { \
		for ( size_t j=0; j<sr->n; j++ ) { \
			sr->vbody[j] = ENTITY(sa, v, 0, j); \
			for ( size_t i=1; i<sa->m; i++ ) { \
				VALUE tmp = ENTITY(sa, v, i, j); \
				if ( RTEST(rb_funcall(sr->vbody[j], op_id, 1, tmp)) ) { \
					sr->vbody[j] = tmp; \
				} \
			} \
		} \
	} else { \
		VALUE tmp = kmm_mat_##mm##_destl(km_Mat(1, sr->n, VT_VALUE), kmm_mat_to_omat(va)); \
		if ( sr->vtype == VT_COMPLEX ) { \
			return kmm_mat_copy_from(self, kmm_mat_to_cmat(tmp)); \
		} else { \
			rb_raise(km_eInternal, "unexpected reach"); \
		} \
	} \
	return self; \
} \
static VALUE \
km_mat_##mm##_row(VALUE self, SMAT *sr, SMAT *sa, VALUE va) \
{ \
	if ( sr->vtype == VT_DOUBLE ) { \
		for ( size_t i=0; i<sr->m; i++ ) { \
			sr->dbody[i] = ENTITY(sa, d, i, 0); \
			for ( size_t j=1; j<sa->n; j++ ) { \
				double tmp = ENTITY(sa, d, i, j); \
				if ( sr->dbody[i] op tmp ) { \
					sr->dbody[i] = tmp; \
				} \
			} \
		} \
	} else if ( sr->vtype == VT_INT ) { \
		for ( size_t i=0; i<sr->m; i++ ) { \
			sr->ibody[i] = ENTITY(sa, i, i, 0); \
			for ( size_t j=1; j<sa->n; j++ ) { \
				int tmp = ENTITY(sa, i, i, j); \
				if ( sr->ibody[i] op tmp ) { \
					sr->ibody[i] = tmp; \
				} \
			} \
		} \
	} else if ( sr->vtype == VT_BOOL ) { \
		for ( size_t i=0; i<sr->m; i++ ) { \
			sr->bbody[i] = ENTITY(sa, b, i, 0); \
			if ( !(opb sr->bbody[i]) ) { \
				for ( size_t j=0; j<sa->n; j++ ) { \
					bool tmp = ENTITY(sa, b, i, j); \
					if ( opb tmp ) { \
						sr->bbody[i] = tmp; \
					} \
				} \
			} \
		} \
	} else if ( sr->vtype == VT_VALUE ) { \
		for ( size_t i=0; i<sr->m; i++ ) { \
			sr->vbody[i] = ENTITY(sa, v, i, 0); \
			for ( size_t j=1; j<sa->n; j++ ) { \
				VALUE tmp = ENTITY(sa, v, i, j); \
				if ( RTEST(rb_funcall(sr->vbody[i], op_id, 1, tmp)) ) { \
					sr->vbody[i] = tmp; \
				} \
			} \
		} \
	} else { \
		VALUE tmp = kmm_mat_##mm##_destl(km_Mat(sr->m, 1, VT_VALUE), kmm_mat_to_omat(va)); \
		if ( sr->vtype == VT_COMPLEX ) { \
			return kmm_mat_copy_from(self, kmm_mat_to_cmat(tmp)); \
		} else { \
			rb_raise(km_eInternal, "unexpected reach"); \
		} \
	} \
	return self; \
}
DEFINE_MM(max, <, id_op_lt, )
DEFINE_MM(min, <, id_op_gt, !)
VALUE
kmm_mat_max_destl(VALUE self, VALUE vr)
{
	km_check_frozen(self);
	SMAT *sr = km_mat2smat(vr), *sa = km_mat2smat(self);
	if ( sr->vtype != sa->vtype ) { rb_raise(km_eVT, "value types must be the same"); }
	if ( sr->stype != ST_FULL ) { rb_raise(km_eShare, "self must be full matrix"); }
	if ( sr->trans ) {
		if ( kmm_mat_have_submatrix_p(vr) ) {
			VALUE mm = kmm_mat_max_destl(self, rb_obj_dup(vr));
			if ( sr->m == 1 && sr->n == 1 ) {
				kmm_mat_set_value(vr, INT2NUM(0), INT2NUM(0), mm);
				return mm;
			} else {
				kmm_mat_copy_from(vr, mm);
				return vr;
			}
		} else {
			sr->trans = false; sr->ld = sr->m; // this is available because the contents of sr are ignored
		}
	}
	if ( sr->m == 1 && sr->n == 1 ) {
		return km_mat_max_all(sr, sa, self);
	} else if ( sr->m==1 && sr->n==sa->n ) {
		return km_mat_max_col(vr, sr, sa, self);
	} else if ( sr->m==sa->m && sr->n==1 ) {
		return km_mat_max_row(vr, sr, sa, self);
	} else if ( SAME_SIZE(sr, sa) ) {
		km_smat_copy(sr, sa);
		return self;
	} else {
		rb_raise(km_eDim, "max/min from size (%zu, %zu) to size (%zu, %zu) is not defined", sa->m, sa->n, sr->m, sr->n);
	}
}
// the argument is a Symbol which specify the axis
// if the argument is :col take a max in column axis and return row vector
// if the argument is :row take a max in row axis and return column vector
// if the argument is :all or omitted, take a max of all the elements
VALUE
kmm_mat_max(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 0, 1);
	if ( argc == 0 ) {
		return kmm_mat_max_destl(self, km_Mat(1, 1, km_mat2smat(self)->vtype));
	} else {
		VALUE asym = argv[0];
		if ( rb_respond_to(asym, id_to_sym) ) {
			asym = rb_funcall(asym, id_to_sym, 0);
		}
		SMAT *sa = km_mat2smat(self);
		if ( asym == sym_col || asym == sym_c ) {
			return kmm_mat_max_destl(self, km_Mat(1, sa->n, sa->vtype));
		} else if ( asym == sym_row || asym == sym_r ) {
			return kmm_mat_max_destl(self, km_Mat(sa->m, 1, sa->vtype));
		} else if ( asym == sym_all || asym == sym_a ) {
			return kmm_mat_max_destl(self, km_Mat(1, 1, sa->vtype));
		} else if ( asym == sym_none || asym == sym_n ) {
			return kmm_mat_max_destl(self, km_Mat(sa->m, sa->n, sa->vtype));
		} else {
			rb_raise(rb_eArgError, "unknown axis name");
		}
	}
}

VALUE
kmm_mat_min_destl(VALUE self, VALUE vr)
{
	km_check_frozen(self);
	SMAT *sr = km_mat2smat(vr), *sa = km_mat2smat(self);
	if ( sr->vtype != sa->vtype ) { rb_raise(km_eVT, "value types must be the same"); }
	if ( sr->stype != ST_FULL ) { rb_raise(km_eShare, "self must be full matrix"); }
	if ( sr->trans ) {
		if ( kmm_mat_have_submatrix_p(vr) ) {
			VALUE mm = kmm_mat_min_destl(self, rb_obj_dup(vr));
			if ( sr->m == 1 && sr->n == 1 ) {
				kmm_mat_set_value(vr, INT2NUM(0), INT2NUM(0), mm);
				return mm;
			} else {
				kmm_mat_copy_from(vr, mm);
				return vr;
			}
		} else {
			sr->trans = false; sr->ld = sr->m; // this is available because the contents of sr are ignored
		}
	}
	if ( sr->m == 1 && sr->n == 1 ) {
		return km_mat_min_all(sr, sa, self);
	} else if ( sr->m==1 && sr->n==sa->n ) {
		return km_mat_min_col(vr, sr, sa, self);
	} else if ( sr->m==sa->m && sr->n==1 ) {
		return km_mat_min_row(vr, sr, sa, self);
	} else if ( SAME_SIZE(sr, sa) ) {
		km_smat_copy(sr, sa);
		return self;
	} else {
		rb_raise(km_eDim, "max/min from size (%zu, %zu) to size (%zu, %zu) is not defined", sa->m, sa->n, sr->m, sr->n);
	}
}
VALUE
kmm_mat_min(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 0, 1);
	if ( argc == 0 ) {
		return kmm_mat_min_destl(self, km_Mat(1, 1, km_mat2smat(self)->vtype));
	} else {
		VALUE asym = argv[0];
		if ( rb_respond_to(asym, id_to_sym) ) {
			asym = rb_funcall(asym, id_to_sym, 0);
		}
		SMAT *sa = km_mat2smat(self);
		if ( asym == sym_col || asym == sym_c ) {
			return kmm_mat_min_destl(self, km_Mat(1, sa->n, sa->vtype));
		} else if ( asym == sym_row || asym == sym_r ) {
			return kmm_mat_min_destl(self, km_Mat(sa->m, 1, sa->vtype));
		} else if ( asym == sym_all || asym == sym_a ) {
			return kmm_mat_min_destl(self, km_Mat(1, 1, sa->vtype));
		} else if ( asym == sym_none || asym == sym_n ) {
			return kmm_mat_min_destl(self, km_Mat(sa->m, sa->n, sa->vtype));
		} else {
			rb_raise(rb_eArgError, "unknown axis name");
		}
	}
}

struct km_amm_arg {
	size_t i;
	size_t j;
	union {
		double dmm;
		int imm;
		bool bmm;
		VALUE vmm;
	};
};
static void
km_argmax_func_d(double *ea, size_t i, size_t j, void *data)
{
	struct km_amm_arg *arg = (struct km_amm_arg *)data;
	if ( arg->dmm < *ea ) {
		arg->i = i; arg->j = j; arg->dmm = *ea;
	}
}
static void
km_argmax_func_i(int *ea, size_t i, size_t j, void *data)
{
	struct km_amm_arg *arg = (struct km_amm_arg *)data;
	if ( arg->imm < *ea ) {
		arg->i = i; arg->j = j; arg->imm = *ea;
	}
}
static void
km_argmax_func_b(bool *ea, size_t i, size_t j, void *data)
{
	struct km_amm_arg *arg = (struct km_amm_arg *)data;
	if ( (!(arg->bmm)) && *ea ) {
		arg->i = i; arg->j = j; arg->bmm = true;
	}
}
static void
km_argmax_func_v(VALUE *ea, size_t i, size_t j, void *data)
{
	struct km_amm_arg *arg = (struct km_amm_arg *)data;
	if ( RTEST(rb_funcall(arg->vmm, id_op_lt, 1, *ea)) ) {
		arg->i = i; arg->j = j; arg->vmm = *ea;
	}
}
#define DEFINE_AMM(mm, op, op_id, opb) \
static VALUE \
km_mat_arg##mm##_all(SMAT *sr, SMAT *sa, VALUE va) \
{ \
	if ( sa->vtype == VT_DOUBLE ) { \
		struct km_amm_arg data = {0, 0, {.dmm = sa->dbody[0]}}; \
		km_smat_each_with_index_d(sa, km_arg##mm##_func_d, &data); \
		sr->ibody[0] = INDEXi(sa, data.i, data.j); \
		if ( VECTOR_P(sa) ) { \
			return ZU2NUM(MAX(data.i, data.j)); \
		} else { \
			return rb_ary_new3(2, ZU2NUM(data.i), ZU2NUM(data.j)); \
		} \
	} else if ( sa->vtype == VT_INT ) { \
		struct km_amm_arg data = {0, 0, {.imm = sa->ibody[0]}}; \
		km_smat_each_with_index_i(sa, km_arg##mm##_func_i, &data); \
		sr->ibody[0] = INDEXi(sa, data.i, data.j); \
		if ( VECTOR_P(sa) ) { \
			return ZU2NUM(MAX(data.i, data.j)); \
		} else { \
			return rb_ary_new3(2, ZU2NUM(data.i), ZU2NUM(data.j)); \
		} \
	} else if ( sa->vtype == VT_BOOL ) { \
		struct km_amm_arg data = {0, 0, {.bmm = sa->bbody[0]}}; \
		km_smat_each_with_index_b(sa, km_arg##mm##_func_b, &data); \
		sr->ibody[0] = INDEXi(sa, data.i, data.j); \
		if ( VECTOR_P(sa) ) { \
			return ZU2NUM(MAX(data.i, data.j)); \
		} else { \
			return rb_ary_new3(2, ZU2NUM(data.i), ZU2NUM(data.j)); \
		} \
	} else if ( sa->vtype == VT_VALUE ) { \
		struct km_amm_arg data = {0, 0, {.vmm = sa->vbody[0]}}; \
		km_smat_each_with_index_v(sa, km_arg##mm##_func_v, &data); \
		sr->ibody[0] = INDEXi(sa, data.i, data.j); \
		if ( VECTOR_P(sa) ) { \
			return ZU2NUM(MAX(data.i, data.j)); \
		} else { \
			return rb_ary_new3(2, ZU2NUM(data.i), ZU2NUM(data.j)); \
		} \
	} else { \
		VALUE omat = kmm_mat_to_omat(va); \
		return km_mat_arg##mm##_all(sr, km_mat2smat(omat), omat); \
	} \
} \
static VALUE \
km_mat_arg##mm##_col(VALUE self, SMAT *sr, SMAT *sa, VALUE va) \
{ \
	if ( sa->vtype == VT_DOUBLE ) { \
		for ( size_t j=0; j<sr->n; j++ ) { \
			sr->ibody[j] = 0; double mm = ENTITY(sa, d, 0, j); \
			for ( size_t i=0; i<sa->m; i++ ) { \
				double tmp = ENTITY(sa, d, i, j); \
				if ( mm op tmp ) { \
					sr->ibody[j] = s2i(i); mm = tmp; \
				} \
			} \
		} \
	} else if ( sa->vtype == VT_INT ) { \
		for ( size_t j=0; j<sr->n; j++ ) { \
			sr->ibody[j] = 0; int mm = ENTITY(sa, i, 0, j); \
			for ( size_t i=0; i<sa->m; i++ ) { \
				int tmp = ENTITY(sa, i, i, j); \
				if ( mm op tmp ) { \
					sr->ibody[j] = s2i(i); mm = tmp; \
				} \
			} \
		} \
	} else if ( sa->vtype == VT_BOOL ) { \
		for ( size_t j=0; j<sr->n; j++ ) { \
			sr->ibody[j] = 0; bool mm = ENTITY(sa, b, 0, j); \
			if ( !( opb mm ) ) { \
				for ( size_t i=0; i<sa->m; i++ ) { \
					int tmp = ENTITY(sa, b, i, j); \
					if ( opb tmp ) { \
						sr->bbody[j] = s2i(i); break; \
					} \
				} \
			} \
		} \
	} else if ( sa->vtype == VT_VALUE ) { \
		for ( size_t j=0; j<sr->n; j++ ) { \
			sr->ibody[j] = 0; VALUE mm = ENTITY(sa, v, 0, j); \
			for ( size_t i=0; i<sa->m; i++ ) { \
				VALUE tmp = ENTITY(sa, v, i, j); \
				if ( rb_funcall(mm, op_id, 1, tmp) ) { \
					sr->ibody[j] = s2i(i); mm = tmp; \
				} \
			} \
		} \
	} else { \
		VALUE omat = kmm_mat_to_omat(va); \
		return km_mat_arg##mm##_col(self, sr, km_mat2smat(omat), omat); \
	} \
	return self; \
} \
static VALUE \
km_mat_arg##mm##_row(VALUE self, SMAT *sr, SMAT *sa, VALUE va) \
{ \
	if ( sa->vtype == VT_DOUBLE ) { \
		for ( size_t i=0; i<sr->m; i++ ) { \
			sr->ibody[i] = 0; double mm = ENTITY(sa, d, i, 0); \
			for ( size_t j=1; j<sa->n; j++ ) { \
				double tmp = ENTITY(sa, d, i, j); \
				if ( mm op tmp ) { \
					sr->ibody[i] = s2i(j); mm = tmp; \
				} \
			} \
		} \
	} else if ( sa->vtype == VT_INT ) { \
		for ( size_t i=0; i<sr->m; i++ ) { \
			sr->ibody[i] = 0; int mm = ENTITY(sa, i, i, 0); \
			for ( size_t j=1; j<sa->n; j++ ) { \
				int tmp = ENTITY(sa, i, i, j); \
				if ( mm op tmp ) { \
					sr->ibody[i] = s2i(j); mm = tmp; \
				} \
			} \
		} \
	} else if ( sa->vtype == VT_BOOL ) { \
		for ( size_t i=0; i<sr->m; i++ ) { \
			sr->ibody[i] = 0; int mm = ENTITY(sa, b, i, 0); \
			if ( ! (opb mm) ) { \
				for ( size_t j=1; j<sa->n; j++ ) { \
					int tmp = ENTITY(sa, i, i, j); \
					if ( opb tmp ) { \
						sr->ibody[i] = s2i(j); break; \
					} \
				} \
			} \
		} \
	} else if ( sa->vtype == VT_VALUE ) { \
		for ( size_t i=0; i<sr->m; i++ ) { \
			sr->ibody[i] = 0; VALUE mm = ENTITY(sa, v, i, 0); \
			for ( size_t j=1; j<sa->n; j++ ) { \
				VALUE tmp = ENTITY(sa, v, i, j); \
				if ( rb_funcall(mm, op_id, 1, tmp) ) { \
					sr->ibody[i] = s2i(j); mm = tmp; \
				} \
			} \
		} \
	} else { \
		VALUE omat = kmm_mat_to_omat(va); \
		return km_mat_arg##mm##_row(self, sr, km_mat2smat(omat), omat); \
	} \
	return self; \
}
DEFINE_AMM(max, <, id_op_lt, )
VALUE
kmm_mat_argmax_destl(VALUE self, VALUE va)
{
	km_check_frozen(self); SMAT *sr = km_mat2smat(self), *sa = km_mat2smat(va);
	if ( sr->vtype != VT_INT || sr->stype != ST_FULL ) {
		rb_raise(rb_eRuntimeError, "self must be an int full matrix");
	}
	if ( sr->trans ) {
		if ( kmm_mat_have_submatrix_p(self) ) {
			return kmm_mat_copy_from(self, kmm_mat_argmax_destl(rb_obj_dup(self), va));
		} else {
			sr->trans = false; sr->ld = sr->m;
		}
	}
	if ( sr->m == 1 && sr->n == 1 ) {
		return km_mat_argmax_all(sr, sa, va);
	} else if ( sr->m == 1 && sr->n == sa->n ) {
		return km_mat_argmax_col(self, sr, sa, va);
	} else if ( sr->m == sa->m && sr->n == 1 ) {
		return km_mat_argmax_row(self, sr, sa, va);
	} else {
		rb_raise(km_eDim, "argmax/min from size (%zu, %zu) to size (%zu, %zu) is not defined", sa->m, sa->n, sr->m, sr->n);
	}
}
VALUE
kmm_mat_argmax(int argc, VALUE *argv, VALUE va)
{
	rb_check_arity(argc, 0, 1);
	if ( argc == 0 ) {
		return kmm_mat_argmax_destl(km_Mat(1, 1, VT_INT), va);
	} else {
		VALUE asym = argv[0];
		if ( rb_respond_to(asym, id_to_sym) ) {
			asym = rb_funcall(asym, id_to_sym, 0);
		}
		SMAT *sa = km_mat2smat(va);
		if ( asym == sym_col || asym == sym_c ) {
			return kmm_mat_argmax_destl(km_Mat(1, sa->n, VT_INT), va);
		} else if ( asym == sym_row || asym == sym_r ) {
			return kmm_mat_argmax_destl(km_Mat(sa->m, 1, VT_INT), va);
		} else if ( asym == sym_all || asym == sym_a ) {
			return kmm_mat_argmax_destl(km_Mat(1, 1, VT_INT), va);
		} else {
			rb_raise(rb_eArgError, "unknown axis name");
		}
	}
}

static void
km_argmin_func_d(double *ea, size_t i, size_t j, void *data)
{
	struct km_amm_arg *arg = (struct km_amm_arg *)data;
	if ( arg->dmm > *ea ) {
		arg->i = i; arg->j = j; arg->dmm = *ea;
	}
}
static void
km_argmin_func_i(int *ea, size_t i, size_t j, void *data)
{
	struct km_amm_arg *arg = (struct km_amm_arg *)data;
	if ( arg->imm > *ea ) {
		arg->i = i; arg->j = j; arg->imm = *ea;
	}
}
static void
km_argmin_func_b(bool *ea, size_t i, size_t j, void *data)
{
	struct km_amm_arg *arg = (struct km_amm_arg *)data;
	if ( (arg->bmm) && (!(*ea)) ) {
		arg->i = i; arg->j = j; arg->bmm = false;
	}
}
static void
km_argmin_func_v(VALUE *ea, size_t i, size_t j, void *data)
{
	struct km_amm_arg *arg = (struct km_amm_arg *)data;
	if ( RTEST(rb_funcall(arg->vmm, id_op_gt, 1, *ea)) ) {
		arg->i = i; arg->j = j; arg->vmm = *ea;
	}
}
DEFINE_AMM(min, >, id_op_gt, !)
VALUE
kmm_mat_argmin_destl(VALUE self, VALUE va)
{
	km_check_frozen(self); SMAT *sr = km_mat2smat(self), *sa = km_mat2smat(va);
	if ( sr->vtype != VT_INT || sr->stype != ST_FULL ) {
		rb_raise(rb_eRuntimeError, "self must be an int full matrix");
	}
	if ( sr->trans ) {
		if ( kmm_mat_have_submatrix_p(self) ) {
			return kmm_mat_copy_from(self, kmm_mat_argmin_destl(rb_obj_dup(self), va));
		} else {
			sr->trans = false; sr->ld = sr->m;
		}
	}
	if ( sr->m == 1 && sr->n == 1 ) {
		return km_mat_argmin_all(sr, sa, va);
	} else if ( sr->m == 1 && sr->n == sa->n ) {
		return km_mat_argmin_col(self, sr, sa, va);
	} else if ( sr->m == sa->m && sr->n == 1 ) {
		return km_mat_argmin_row(self, sr, sa, va);
	} else {
		rb_raise(km_eDim, "argmax/min from size (%zu, %zu) to size (%zu, %zu) is not defined", sa->m, sa->n, sr->m, sr->n);
	}
}
VALUE
kmm_mat_argmin(int argc, VALUE *argv, VALUE va)
{
	rb_check_arity(argc, 0, 1);
	if ( argc == 0 ) {
		return kmm_mat_argmin_destl(km_Mat(1, 1, VT_INT), va);
	} else {
		VALUE asym = argv[0];
		if ( rb_respond_to(asym, id_to_sym) ) {
			asym = rb_funcall(asym, id_to_sym, 0);
		}
		SMAT *sa = km_mat2smat(va);
		if ( asym == sym_col || asym == sym_c ) {
			return kmm_mat_argmin_destl(km_Mat(1, sa->n, VT_INT), va);
		} else if ( asym == sym_row || asym == sym_r ) {
			return kmm_mat_argmin_destl(km_Mat(sa->m, 1, VT_INT), va);
		} else if ( asym == sym_all || asym == sym_a ) {
			return kmm_mat_argmin_destl(km_Mat(1, 1, VT_INT), va);
		} else {
			rb_raise(rb_eArgError, "unknown axis name");
		}
	}
}
