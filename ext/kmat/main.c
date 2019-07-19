#include "kmat.h"

// km_global_variables_begin

VALUE km_cMat;
VALUE km_sMat;
VALUE km_eDim;
VALUE km_eVT;
VALUE km_eUncomp;
VALUE km_eInternal;
VALUE km_eShare;
VALUE km_eNotImp;

VALUE rb_mMarshal;
VALUE rb_mObjSpace;
VALUE rb_sMath;

// km_global_variables_end

#include "id_sym.c"

// return VALUE as Integer
// this is useful for debugging
static VALUE
kmm_obj_value(VALUE self)
{
	return INT2NUM((int)self);
}

#include "method_definitions.c"
#include "elementwise_function.c"
#include "elementwise_function_definitions.c"

void
Init_kmat(void)
{
	rb_mMarshal = rb_const_get(rb_cObject, rb_intern("Marshal"));
	rb_mObjSpace = rb_const_get(rb_cObject, rb_intern("ObjectSpace"));
	rb_sMath = rb_singleton_class(rb_mMath);

	// km_cMat = rb_define_class("Mat", rb_cObject); # Mat is defined in /lib/kmat.rb
	km_cMat = rb_const_get(rb_cObject, rb_intern("Mat"));
	km_sMat = rb_singleton_class(km_cMat);
	rb_undef_alloc_func(km_cMat);
	rb_define_alloc_func(km_cMat, km_Mat_alloc);
	
	km_eDim = rb_define_class_under(km_cMat, "MismatchedDimensionError", rb_eStandardError);
	km_eVT = rb_define_class_under(km_cMat, "ValueTypeError", rb_eStandardError);
	km_eUncomp = rb_define_class_under(km_cMat, "UncomputableMatrixError", rb_eStandardError);
	km_eInternal = rb_define_class_under(km_cMat, "InternalError", rb_eException);
	km_eShare = rb_define_class_under(km_cMat, "SharingError", rb_eStandardError);
	km_eNotImp = rb_define_class_under(km_cMat, "NotImplementedYetError", rb_eStandardError);

	km_init_id(); km_init_sym();

	rb_define_alias(km_cMat, "+@", "itself");
	km_define_methods();
	km_define_efs();
	km_Mat_rand_init();
}

// invoke `func'(`data') disabling GC
// the state of GC will be restored before return
// return the return value of `func'(`data') unless an exception is raised
VALUE
km_gc_escape(VALUE (*func)(), VALUE data)
{
	int status;
	VALUE old = rb_gc_disable();
	VALUE ret = rb_protect(func, data, &status);
	if ( old == Qfalse ) {
		rb_gc_enable();
	}
	if ( status != 0 ) {
		rb_jump_tag(status);
	} else {
		return ret;
	}
}

// invoke `b_proc'(`data1') and ensure invoking `e_proc'(`data2')
// if an exception is raised in `b_proc', the exception is re-raised after invoking `e_proc'
// otherwise, return the return value of `b_proc'
VALUE
km_ensure(VALUE (* b_proc)(ANYARGS), VALUE data1, VALUE (* e_proc)(ANYARGS), VALUE data2)
{
	int status;
	VALUE ret = rb_protect(b_proc, data1, &status);
	(*e_proc)(data2);
	if ( status != 0 ) {
		rb_jump_tag(status);
	} else {
		return ret;
	}
}
