#include "../kmat.h"

VALUE
kmm_MATH_exp2(VALUE self, VALUE x)
{
	return rb_float_new(exp2(NUM2DBL(x)));
}

VALUE
kmm_MATH_expm1(VALUE self, VALUE x)
{
	return rb_float_new(expm1(NUM2DBL(x)));
}

VALUE
kmm_MATH_log1p(VALUE self, VALUE x)
{
	return rb_float_new(log1p(NUM2DBL(x)));
}

VALUE
kmm_float_sign(VALUE self)
{
	const double x = NUM2DBL(self);
	if ( x == 0.0 ) {
		return INT2NUM(0);
	} else if ( x > 0.0 ) {
		return INT2NUM(1);
	} else if ( x < 0.0 ) {
		return INT2NUM(-1);
	} else {
		return self;
	}
}
