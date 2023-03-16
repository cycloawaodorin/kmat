#include "../kmat.h"

static VALUE
km_block(RB_BLOCK_CALL_FUNC_ARGLIST(idx, ary))
{
	VALUE row = rb_ary_entry(ary, NUM2LONG(rb_ary_entry(idx, 0)));
	return rb_ary_entry(row, NUM2LONG(rb_ary_entry(idx, 1)));
}
static VALUE
km_block_vector(RB_BLOCK_CALL_FUNC_ARGLIST(idx, ary))
{
	return rb_ary_entry(ary, NUM2LONG(rb_ary_entry(idx, 0)));
}
static long
km_check_col_size(VALUE ary)
{
	long ret = RARRAY_LEN(rb_ary_entry(ary, 0));
	for ( long i=1; i<RARRAY_LEN(ary); i++ ) {
		if ( RARRAY_LEN(rb_ary_entry(ary, i)) != ret ) {
			rb_raise(km_eDim, "the length of rows must be the same");
		}
	}
	return ret;
}

VALUE
kmm_ary_to_fmat(VALUE self)
{
	if ( RARRAY_LEN(self) == 0 ) {
		return km_Mat(0, 0, VT_DOUBLE);
	} else if ( TYPE(rb_ary_entry(self, 0)) == T_ARRAY ) {
		VALUE new_arg[3] = {LONG2NUM(RARRAY_LEN(self)), LONG2NUM(km_check_col_size(self)), sym_float};
		return rb_block_call(km_cMat, id_new, 3, new_arg, km_block, self);
	} else {
		VALUE new_arg[3] = {LONG2NUM(RARRAY_LEN(self)), INT2NUM(1), sym_float};
		return rb_block_call(km_cMat, id_new, 3, new_arg, km_block_vector, self);
	}
}
VALUE
kmm_ary_to_cmat(VALUE self)
{
	if ( RARRAY_LEN(self) == 0 ) {
		return km_Mat(0, 0, VT_COMPLEX);
	} else if ( TYPE(rb_ary_entry(self, 0)) == T_ARRAY ) {
		VALUE new_arg[3] = {LONG2NUM(RARRAY_LEN(self)), LONG2NUM(km_check_col_size(self)), sym_complex};
		return rb_block_call(km_cMat, id_new, 3, new_arg, km_block, self);
	} else {
		VALUE new_arg[3] = {LONG2NUM(RARRAY_LEN(self)), INT2NUM(1), sym_complex};
		return rb_block_call(km_cMat, id_new, 3, new_arg, km_block_vector, self);
	}
}
VALUE
kmm_ary_to_imat(VALUE self)
{
	if ( RARRAY_LEN(self) == 0 ) {
		return km_Mat(0, 0, VT_INT);
	} else if ( TYPE(rb_ary_entry(self, 0)) == T_ARRAY ) {
		VALUE new_arg[3] = {LONG2NUM(RARRAY_LEN(self)), LONG2NUM(km_check_col_size(self)), sym_int};
		return rb_block_call(km_cMat, id_new, 3, new_arg, km_block, self);
	} else {
		VALUE new_arg[3] = {LONG2NUM(RARRAY_LEN(self)), INT2NUM(1), sym_int};
		return rb_block_call(km_cMat, id_new, 3, new_arg, km_block_vector, self);
	}
}
VALUE
kmm_ary_to_bmat(VALUE self)
{
	if ( RARRAY_LEN(self) == 0 ) {
		return km_Mat(0, 0, VT_BOOL);
	} else if ( TYPE(rb_ary_entry(self, 0)) == T_ARRAY ) {
		VALUE new_arg[3] = {LONG2NUM(RARRAY_LEN(self)), LONG2NUM(km_check_col_size(self)), sym_bool};
		return rb_block_call(km_cMat, id_new, 3, new_arg, km_block, self);
	} else {
		VALUE new_arg[3] = {LONG2NUM(RARRAY_LEN(self)), INT2NUM(1), sym_bool};
		return rb_block_call(km_cMat, id_new, 3, new_arg, km_block_vector, self);
	}
}
VALUE
kmm_ary_to_omat(VALUE self)
{
	if ( RARRAY_LEN(self) == 0 ) {
		return km_Mat(0, 0, VT_VALUE);
	} else if ( TYPE(rb_ary_entry(self, 0)) == T_ARRAY ) {
		VALUE new_arg[3] = {LONG2NUM(RARRAY_LEN(self)), LONG2NUM(km_check_col_size(self)), sym_object};
		return rb_block_call(km_cMat, id_new, 3, new_arg, km_block, self);
	} else {
		VALUE new_arg[3] = {LONG2NUM(RARRAY_LEN(self)), INT2NUM(1), sym_object};
		return rb_block_call(km_cMat, id_new, 3, new_arg, km_block_vector, self);
	}
}

// alias to_m
VALUE
kmm_ary_to_mat(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 0, 1);
	if ( argc == 0 ) {
		return kmm_ary_to_fmat(self);
	} else {
		VT_SWITCH( km_sym2vt(argv[0]),
			return kmm_ary_to_fmat(self);,
			return kmm_ary_to_cmat(self);,
			return kmm_ary_to_imat(self);,
			return kmm_ary_to_bmat(self);,
			return kmm_ary_to_omat(self);
		);
	}
}
