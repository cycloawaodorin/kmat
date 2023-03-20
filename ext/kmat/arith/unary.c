#include "../kmat.h"

static void
km_abs_func_d(double *ent, void *data)
{
	*ent = fabs(*ent);
}
static void
km_abs_func_z(COMPLEX *ent, void *data)
{
	*ent = cpack(cabs(*ent), 0.0);
}
static void
km_abs_func_i(int *ent, void *data)
{
	*ent = ABS(*ent);
}
static void
km_abs_func_v(VALUE *ent, void *data)
{
	*ent = rb_funcall(*ent, id_abs, 0);
}
VALUE
kmm_mat_abs_destl(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	VT_SWITCH( smat->vtype,
		km_smat_each_d(smat, km_abs_func_d, NULL);,
		km_smat_each_z(smat, km_abs_func_z, NULL);,
		km_smat_each_i(smat, km_abs_func_i, NULL);,
		/* nothing to do */,
		km_smat_each_v(smat, km_abs_func_v, NULL);
	);
	return self;
}
static void
km_cabs_func(double *dest, const COMPLEX *src,  void *data)
{
	*dest = cabs(*src);
}
VALUE
kmm_mat_abs(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype == VT_COMPLEX ) {
		VALUE ret = km_Mat(smat->m, smat->n, VT_DOUBLE);
		km_smat_each2_dcz(km_mat2smat(ret), smat, km_cabs_func, NULL);
		return ret;
	} else {
		return kmm_mat_abs_destl(rb_obj_dup(self));
	}
}

static void
km_uminus_func_d(double *ent, void *data)
{
	*ent = -(*ent);
}
static void
km_uminus_func_z(COMPLEX *ent, void *data)
{
	*ent = -(*ent);
}
static void
km_uminus_func_i(int *ent, void *data)
{
	*ent = -(*ent);
}
static void
km_uminus_func_v(VALUE *ent, void *data)
{
	*ent = rb_funcall(*ent, id_op_uminus, 0);
}
VALUE
kmm_mat_uminus_destl(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	VT_SWITCH( smat->vtype,
		km_smat_each_d(smat, km_uminus_func_d, NULL);,
		km_smat_each_z(smat, km_uminus_func_z, NULL);,
		km_smat_each_i(smat, km_uminus_func_i, NULL);,
		/* nothing to do */,
		km_smat_each_v(smat, km_uminus_func_v, NULL);
	);
	return self;
}
// alias -@
VALUE
kmm_mat_uminus(VALUE self)
{
	return kmm_mat_uminus_destl(rb_obj_dup(self));
}

static void
km_conj_func_z(COMPLEX *ent, void *data)
{
	*ent = conj(*ent);
}
static void
km_conj_func_v(VALUE *ent, void *data)
{
	*ent = rb_funcall(*ent, id_conj, 0);
}
// alias conjugate
VALUE
kmm_mat_conj_dest(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype == VT_COMPLEX ) {
		km_smat_each_z(smat, km_conj_func_z, NULL);
	} else if ( smat->vtype == VT_VALUE ) {
		km_smat_each_v(smat, km_conj_func_v, NULL);
	}
	return self;
}

// trace
// alias tr
VALUE
kmm_mat_trace(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	const size_t len = MIN(smat->m, smat->n);
	if ( smat->vtype == VT_DOUBLE ) {
		double ret = 0.0;
		for ( size_t i=0; i<len; i++ ) {
			ret += ENTITY(smat, d, i, i);
		}
		return rb_float_new(ret);
	} else if ( smat->vtype == VT_COMPLEX ) {
		COMPLEX ret = cpack(0.0, 0.0);
		for ( size_t i=0; i<len; i++ ) {
			ret += ENTITY(smat, z, i, i);
		}
		return km_c2v(ret);
	} else if ( smat->vtype == VT_INT ) {
		int ret = 0;
		for ( size_t i=0; i<len; i++ ) {
			ret += ENTITY(smat, i, i, i);
		}
		return INT2NUM(ret);
	} else if ( smat->vtype == VT_BOOL ) {
		bool ret = false;
		for ( size_t i=0; i<len; i++ ) {
			bool ent = ENTITY(smat, b, i, i);
			ret = XOR(ret, ent);
		}
		return TF2V(ret);
	} else if ( smat->vtype == VT_VALUE ) {
		if ( len == 0 ) {
			return INT2NUM(0);
		} else {
			VALUE ret = ENTITY(smat, v, 0, 0);
			for ( size_t i=1; i<len; i++ ) {
				ret = rb_funcall(ret, id_op_plus, 1, ENTITY(smat, v, i, i));
			}
			return ret;
		}
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
}
