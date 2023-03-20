#include "../kmat.h"

// returns sizeof values of the elements with value type `vt'
size_t
km_sizeof_vt(VTYPE vt)
{
	VT_SWITCH( vt,
		return sizeof(double);,
		return sizeof(COMPLEX);,
		return sizeof(int);,
		return sizeof(bool);,
		return sizeof(VALUE);
	);
}

// replace the body of `dest' by complex values with real parts `real' and imaginary parts `imag'
static void
km_smat_complex_func(COMPLEX *ed, const double *er, const double *ei, void *null)
{
	*ed = cpack(*er, *ei);
}

void
km_smat_complex(SMAT *dest, const SMAT *real, const SMAT *imag)
{
	if ( dest->vtype != VT_COMPLEX || real->vtype != VT_DOUBLE || imag->vtype != VT_DOUBLE ) {
		rb_raise(km_eVT, "unmatched value type");
	}
	CHECK_SAME_SIZE(dest, real);
	CHECK_SAME_SIZE(dest, imag);
	dest->trans = false;
	km_smat_each3_zcdcd(dest, real, imag, km_smat_complex_func, NULL);
}

// returns value type from a Symbol `vsym'
VTYPE
km_sym2vt(VALUE vsym)
{
	if ( rb_respond_to(vsym, id_to_sym) ) {
		vsym = rb_funcall(vsym, id_to_sym, 0);
	}
	if ( vsym == sym_float || vsym == sym_f || vsym == sym_double || vsym == sym_d ) {
		return VT_DOUBLE;
	} else if ( vsym == sym_complex || vsym == sym_c || vsym == sym_z ) {
		return VT_COMPLEX;
	} else if ( vsym == sym_int || vsym == sym_integer || vsym == sym_i ) {
		return VT_INT;
	} else if ( vsym == sym_bool || vsym == sym_boolean || vsym == sym_b ) {
		return VT_BOOL;
	} else if ( vsym == sym_object || vsym == sym_o || vsym == sym_value || vsym == sym_v ) {
		return VT_VALUE;
	} else {
		rb_raise(rb_eArgError, "unknown value type");
	}
}

// convert ruby-object `obj' to COMPLEX value in C
COMPLEX
km_v2c(VALUE obj)
{
	VALUE c_obj = rb_funcall(obj, id_to_c, 0);
	const double x = NUM2DBL(rb_funcall(c_obj, id_real, 0));
	const double y = NUM2DBL(rb_funcall(c_obj, id_imag, 0));
	return cpack(x, y);
}

// convert COMPLEX value `c' in C to an instance of Complex in Ruby
VALUE
km_c2v(COMPLEX c)
{
	return rb_complex_new(rb_float_new(creal(c)), rb_float_new(cimag(c)));
}

VALUE
kmm_mat_row_size(VALUE self)
{
	return ZU2NUM(km_mat2smat(self)->m);
}

// alias column_size
VALUE
kmm_mat_col_size(VALUE self)
{
	return ZU2NUM(km_mat2smat(self)->n);
}

VALUE
kmm_mat_length(VALUE self)
{
	return ZU2NUM(LENGTH(km_mat2smat(self)));
}

// alias size
VALUE
kmm_mat_shape(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	return rb_ary_new3(2, ZU2NUM(smat->m), ZU2NUM(smat->n));
}

VALUE
kmm_mat_vtype(VALUE self)
{
	return km_vt2sym(km_mat2smat(self)->vtype);
}

VALUE
km_vt2sym(VTYPE vt)
{
	VT_SWITCH( vt,
		return sym_float;,
		return sym_complex;,
		return sym_int;,
		return sym_bool;,
		return sym_object;
	);
}

VALUE
kmm_Mat_vt_refine(VALUE klass, VALUE vsym)
{
	VTYPE vt = km_sym2vt(vsym);
	VT_SWITCH( vt,
		return sym_float;,
		return sym_complex;,
		return sym_int;,
		return sym_bool;,
		return sym_object;
	);
}

// transform the value type of `self' by a Symbol vsym
// alias vtype=
VALUE
kmm_mat_set_vtype(VALUE self, VALUE vsym)
{
	VT_SWITCH( km_sym2vt(vsym),
		return kmm_mat_to_fmat_destl(self);,
		return kmm_mat_to_cmat_destl(self);,
		return kmm_mat_to_imat_destl(self);,
		return kmm_mat_to_bmat_destl(self);,
		return kmm_mat_to_omat_destl(self);
	);
}

VALUE
kmm_mat_stype(VALUE self)
{
	const STYPE st = km_mat2smat(self)->stype;
	if ( st == ST_FULL ) {
		return sym_full;
	} else if ( st == ST_SSUB ) {
		return sym_serial;
	} else if ( st == ST_RSUB ) {
		return sym_random;
	} else {
		rb_raise(km_eInternal, "unknown storage type");
	}
}

VALUE
kmm_mat_vector_p(VALUE self)
{
	return TF2V(VECTOR_P(km_mat2smat(self)));
}
VALUE
kmm_mat_square_p(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	return TF2V(smat->m==smat->n);
}
