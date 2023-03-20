#include "../kmat.h"

// alias []
VALUE
kmm_Mat_array_m2(VALUE self, VALUE args)
{
	return kmm_ary_to_fmat(args);
}

static VALUE
km_Mat_uninitialized(int argc, VALUE *argv)
{
	rb_check_arity(argc, 0, 2);
	if ( argc == 2 ) {
		return km_Mat(NUM2ZU(argv[0]), NUM2ZU(argv[1]), VT_DOUBLE);
	} else if ( argc == 1 ) {
		return km_Mat(NUM2ZU(argv[0]), 1, VT_DOUBLE);
	} else {
		return km_Mat(0, 0, VT_DOUBLE);
	}
}

// alias O
VALUE
kmm_Mat_zeros(int argc, VALUE *argv, VALUE self)
{
	return km_Mat_uninitialized(argc, argv);
}

static void
km_ones_func(double *ent, void *null)
{
	*ent = 1.0;
}
// alias J
VALUE
kmm_Mat_ones(int argc, VALUE *argv, VALUE self)
{
	VALUE ret = km_Mat_uninitialized(argc, argv);
	SMAT *smat = km_mat2smat(ret);
	km_smat_each_d(smat, km_ones_func, NULL);
	return ret;
}

// make a Mat which filled the last argument
// the value type of the return is determined by the class of the last argument
VALUE
kmm_Mat_fill(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 1, 3);
	VTYPE vt; VALUE ret;
	const VALUE v = argv[argc-1];
	const int type = TYPE(v);
	if ( type == T_FLOAT ) {
		vt = VT_DOUBLE;
	} else if ( type == T_COMPLEX ) {
		const VALUE x = rb_funcall(v, id_real, 0);
		const VALUE y = rb_funcall(v, id_imag, 0);
		if ( ( TYPE(x) == T_FLOAT || TYPE(x) == T_FIXNUM ) && ( TYPE(y) == T_FLOAT || TYPE(y) == T_FIXNUM ) ) {
			vt = VT_COMPLEX;
		} else {
			vt = VT_VALUE;
		}
	} else if ( type == T_FIXNUM ) {
		vt = VT_INT;
	} else if ( type == T_TRUE || type == T_FALSE ) {
		vt = VT_BOOL;
	} else {
		vt = VT_VALUE;
	}
	if ( argc == 3 ) {
		ret = km_Mat(NUM2ZU(argv[0]), NUM2ZU(argv[1]), vt);
	} else if ( argc == 2 ) {
		ret = km_Mat(NUM2ZU(argv[0]), 1, vt);
	} else {
		ret = km_Mat(1, 1, vt);
	}
	return kmm_mat_fill(ret, v);
}

// return a Mat which diagonals are 1.0 and the other elements are 0.0
VALUE
kmm_Mat_eye(int argc, VALUE *argv, VALUE self)
{
	VALUE ret = kmm_Mat_zeros(argc, argv, self);
	SMAT *smat = km_mat2smat(ret);
	const size_t n = MIN(smat->m, smat->n);
	for ( size_t i=0; i<n; i++ ) {
		smat->dbody[i+i*(smat->ld)] = 1.0;
	}
	return ret;
}

// almost the same as Mat.eye but Mat.identity returns only square matrix
// alias I
VALUE
kmm_Mat_identity(VALUE self, VALUE vn)
{
	const size_t n = NUM2ZU(vn);
	VALUE ret = km_Mat(n, n, VT_DOUBLE);
	SMAT *smat = km_mat2smat(ret);
	for ( size_t i=0; i<n; i++ ) {
		smat->dbody[i+i*n] = 1.0;
	}
	return ret;
}

// make a diagonal matrix, the diagonal elements are specified by the arguments
// if the first argument is a vector or an Array, the diagonal elements are the elements of the first argument
// if the first argument is a matrix, the diagonal elements are the diagonal elements of the first argument
// if the second argument `k' given, the elements are set to k-th diagonal
VALUE
kmm_Mat_diag(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 1, 2);
	VALUE vsrc = argv[0];
	if ( argc == 1 ) {
		if ( rb_obj_is_kind_of(vsrc, km_cMat) ) {
			VALUE ret;
			SMAT *src = km_mat2smat(vsrc);
			if ( src->vtype != VT_DOUBLE ) {
				vsrc = kmm_mat_to_fmat(argv[0]);
				src = km_mat2smat(vsrc);
			}
			if ( VECTOR_P(src) ) { // the argument is a vector
				const size_t len = LENGTH(src);
				ret = km_Mat(len, len, VT_DOUBLE);
				SMAT *dest = km_mat2smat(ret);
				if ( src->stype == ST_RSUB ) {
					for ( size_t i=0; i<len; i++ ) {
						dest->dbody[i+i*len] = *((src->dpbody)[i]);
					}
				} else {
					for ( size_t i=0; i<len; i++ ) {
						dest->dbody[i+i*len] = src->dbody[i];
					}
				}
			} else { // the argument is a matrix
				ret = km_Mat(src->m, src->n, VT_DOUBLE);
				SMAT *dest = km_mat2smat(ret);
				const size_t n = MIN(src->m, src->n);
				if ( src->stype == ST_RSUB ) {
					for ( size_t i=0; i<n; i++ ) {
						dest->dbody[i+i*(dest->ld)] = *((src->dpbody)[i+i*(src->ld)]);
					}
				} else {
					for ( size_t i=0; i<n; i++ ) {
						dest->dbody[i+i*(dest->ld)] = src->dbody[i+i*(src->ld)];
					}
				}
			}
			return ret;
		} else if ( TYPE(vsrc) == T_ARRAY ) { // the argument is an Array
			const long len = RARRAY_LEN(argv[0]);
			VALUE ret = km_Mat((size_t)len, (size_t)len, VT_DOUBLE);
			double *body = km_mat2smat(ret)->dbody;
			for ( long i=0; i<len; i++ ) {
				body[l2s(i+i*len)] = NUM2DBL(rb_ary_entry(vsrc, i));
			}
			return ret;
		} else { // (1, 1)-matrix. the (0, 0)-th element is the argument
			VALUE ret = km_Mat(1, 1, VT_DOUBLE);
			km_mat2smat(ret)->dbody[0] = NUM2DBL(vsrc);
			return ret;
		}
	} else { // k is given
		const long k = NUM2LONG(argv[1]);
		if ( rb_obj_is_kind_of(vsrc, km_cMat) ) {
			VALUE ret;
			SMAT *src = km_mat2smat(vsrc);
			if ( src->vtype != VT_DOUBLE ) {
				vsrc = kmm_mat_to_fmat(argv[0]);
				src = km_mat2smat(vsrc);
			}
			if ( VECTOR_P(src) ) { // the argument is a vector
				const size_t len_s = LENGTH(src);
				const size_t len = len_s + l2s(ABS(k));
				ret = km_Mat(len, len, VT_DOUBLE);
				SMAT *dest = km_mat2smat(ret);
				if ( 0 < k ) {
					if ( src->stype == ST_RSUB ) {
						for ( size_t i=0; i<len_s; i++ ) {
							dest->dbody[i+(i+l2s(k))*len] = *((src->dpbody)[i]);
						}
					} else {
						for ( size_t i=0; i<len_s; i++ ) {
							dest->dbody[i+(i+l2s(k))*len] = src->dbody[i];
						}
					}
				} else {
					if ( src->stype == ST_RSUB ) {
						for ( size_t i=0; i<len_s; i++ ) {
							dest->dbody[(i+l2s(-k))+i*len] = *((src->dpbody)[i]);
						}
					} else {
						for ( size_t i=0; i<len_s; i++ ) {
							dest->dbody[(i+l2s(-k))+i*len] = src->dbody[i];
						}
					}
				}
			} else { // the argument is a matrix
				ret = km_Mat(src->m, src->n, VT_DOUBLE);
				SMAT *dest = km_mat2smat(ret);
				size_t len, i_s, j_s;
				if ( 0 < k ) {
					if ( k < s2l(src->n)-s2l(src->m) ) {
						len = src->m;
					} else {
						len = src->n+l2s(-k);
					}
					i_s = 0; j_s = l2s(k);
				} else {
					if ( k < s2l(src->n)-s2l(src->m) ) {
						len = src->m+l2s(k);
					} else {
						len = src->n;
					}
					i_s = l2s(-k); j_s = 0;
				}
				for ( size_t i=0; i<len; i++ ) {
					dest->dbody[i_s+i+(j_s+i)*(dest->ld)] = ENTITY(src, d, i_s+i, j_s+i);
				}
			}
			return ret;
		} else if ( TYPE(vsrc) == T_ARRAY ) { // the argument is an Array
			const long len_s = RARRAY_LEN(argv[0]);
			const size_t len = l2s(len_s+ABS(k));
			VALUE ret = km_Mat(len, len, VT_DOUBLE);
			double *body = km_mat2smat(ret)->dbody;
			if ( 0 < k ) {
				for ( long i=0; i<len_s; i++ ) {
					body[l2s(i)+l2s(i+k)*len] = NUM2DBL(rb_ary_entry(vsrc, i));
				}
			} else {
				for ( long i=0; i<len_s; i++ ) {
					body[l2s(i-k)+l2s(i)*len] = NUM2DBL(rb_ary_entry(vsrc, i));
				}
			}
			return ret;
		} else { // (1+abs(k), 1+abs(k))-matrix. the (0, k)-th element (0<k) or (-k, 0)-th element (k<0) is the argument
			const size_t len = l2s(1+ABS(k));
			VALUE ret = km_Mat(len, len, VT_DOUBLE);
			if ( 0 < k ) {
				km_mat2smat(ret)->dbody[l2s(k)*len] = NUM2DBL(vsrc);
			} else {
				km_mat2smat(ret)->dbody[l2s(-k)] = NUM2DBL(vsrc);
			}
			return ret;
		}
	}
}

// make a vector which contains 0..(vn-1)
VALUE
kmm_Mat_range(VALUE self, VALUE vn)
{
	const size_t n = NUM2ZU(vn);
	VALUE ret = km_Mat(n, 1, VT_DOUBLE);
	double *body = km_mat2smat(ret)->dbody;
	for ( size_t i=0; i<n; i++ ) {
		body[i] = (double)i;
	}
	return ret;
}
VALUE
kmm_Mat_irange(VALUE self, VALUE vn)
{
	const size_t n = NUM2ZU(vn);
	VALUE ret = km_Mat(n, 1, VT_INT);
	int *body = km_mat2smat(ret)->ibody;
	for ( size_t i=0; i<n; i++ ) {
		body[i] = s2i(i);
	}
	return ret;
}

// make a matrix which elements follows U(0, 1)
VALUE
kmm_Mat__rand0(VALUE self, VALUE vm, VALUE vn, VALUE random)
{
	VALUE ret = km_Mat(NUM2ZU(vm), NUM2ZU(vn), VT_DOUBLE);
	return kmm_mat__rand0(ret, random, Qnil);
}

// make a matrix which elements follows N(0, 1)
VALUE
kmm_Mat__randn0(VALUE self, VALUE vm, VALUE vn, VALUE random)
{
	VALUE ret = km_Mat(NUM2ZU(vm), NUM2ZU(vn), VT_DOUBLE);
	return kmm_mat__randn0(ret, random);
}

// make a complex matrix. the real parts are `vr' and the imaginary parts are `vi'
VALUE
kmm_Mat_complex(VALUE self, VALUE vr, VALUE vi)
{
	SMAT *real = km_mat2smat(kmm_mat_to_fmat(vr));
	SMAT *imag = km_mat2smat(kmm_mat_to_fmat(vi));
	VALUE ret = km_Mat(real->m, real->n, VT_COMPLEX);
	km_smat_complex(km_mat2smat(ret), real, imag);
	return ret;
}
