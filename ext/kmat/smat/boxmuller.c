#include "../kmat.h"

VALUE kmgv_mat_random;

// return a random number which follows N(0, 1)
// `random' is an instance of Random
double
km_rand_normal(VALUE random)
{
	if (RTEST(rb_ivar_get(random, id_iv_kmat_stored))) {
		rb_ivar_set(random, id_iv_kmat_stored, Qfalse);
		return NUM2DBL(rb_ivar_get(random, id_iv_kmat_stored_value));
	} else {
		rb_ivar_set(random, id_iv_kmat_stored, Qtrue);
		double sqrt_m2_log_x = sqrt(-2.0*log1p(-NUM2DBL(rb_funcall(random, id_rand, 0))));
		double two_pi_y = M_2PI*NUM2DBL(rb_funcall(random, id_rand, 0));
		rb_ivar_set(random, id_iv_kmat_stored_value, rb_float_new(sqrt_m2_log_x*sin(two_pi_y)));
		return sqrt_m2_log_x*cos(two_pi_y);
	}
}

// replace ary[0..(len-1)] by random numbers which follow N(0, 1)
// `random' is an instance of Random
void
km_fill_normal(int len, double *ary, VALUE random)
{
	int i;
	double sqrt_m2_log_x, two_pi_y;
	len--;
	if (RTEST(rb_ivar_get(random, id_iv_kmat_stored))) {
		ary[0] = NUM2DBL(rb_ivar_get(random, id_iv_kmat_stored_value));
		i = 1;
	} else {
		i = 0;
	}
	for ( ; i<len; i+= 2 ) {
		sqrt_m2_log_x = sqrt(-2.0*log1p(-NUM2DBL(rb_funcall(random, id_rand, 0))));
		two_pi_y = M_2PI*NUM2DBL(rb_funcall(random, id_rand, 0));
		ary[i] = sqrt_m2_log_x*cos(two_pi_y);
		ary[i+1] = sqrt_m2_log_x*sin(two_pi_y);
	}
	if ( i==len ) {
		sqrt_m2_log_x = sqrt(-2.0*log1p(-NUM2DBL(rb_funcall(random, id_rand, 0))));
		two_pi_y = M_2PI*NUM2DBL(rb_funcall(random, id_rand, 0));
		ary[i] = sqrt_m2_log_x*cos(two_pi_y);
		rb_ivar_set(random, id_iv_kmat_stored, Qtrue);
		rb_ivar_set(random, id_iv_kmat_stored_value, rb_float_new(sqrt_m2_log_x*sin(two_pi_y)));
	} else {
		rb_ivar_set(random, id_iv_kmat_stored, Qfalse);
	}
}

static VALUE
km_random_randn(int argc, VALUE *argv, VALUE random)
{
	rb_check_arity(argc, 0, 2);
	if ( argc == 0 ) {
		return rb_float_new(km_rand_normal(random));
	} else if ( argc == 2 ) {
		return rb_float_new(NUM2DBL(argv[0]) + km_rand_normal(random) * sqrt(NUM2DBL(argv[1])));
	} else {
		rb_raise(rb_eArgError, "when mean or standard deviation is specified, the other must be specified");
	}
}

void
km_Mat_rand_init(void)
{
	kmgv_mat_random = rb_funcall(rb_cRandom, id_new, 0);
	rb_define_variable("$MatRandom", &kmgv_mat_random);
	rb_define_method(rb_cRandom, "randn", km_random_randn, -1);
}
