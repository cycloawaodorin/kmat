#include "../kmat.h"

// pairc に，続けて渡す int の組の数を指定(引数の数は1+parc*2)
// いずれかの組の値が異なれば，MismatchedDimensionError

// the number of pairs of SMATs is specified by `argc'
// check whether elements of each pair are the same or not
// raise Mat::MismatchedDimensionError if at least one of the pair has different elements
void
km_check_size(int pairc, ...)
{
	va_list argp;
	va_start(argp, pairc);
	for (int i=0; i<pairc; i++ ) {
		int a = va_arg(argp, int); int b = va_arg(argp, int);
		if ( a != b ) {
			va_end(argp);
			rb_raise(km_eDim, "dimension mismatched (%d != %d)", a, b);
		}
	}
	va_end(argp);
}

// the number of SMATs is specified by `argc'
// check whether the value types of the arguments are VT_DOUBLE or not
// raise Mat::ValueTypeError if at least one of the arguments is not a float matrix
void
km_check_double(int argc, ...)
{
	va_list argp;
	va_start(argp, argc);
	for ( int i=0; i<argc; i++ ) {
		SMAT *smat = va_arg(argp, SMAT *);
		if ( smat->vtype != VT_DOUBLE ) {
			va_end(argp);
			rb_raise(km_eVT, "float matrix is expected");
		}
	}
	va_end(argp);
}

// the number of SMATs is specified by `argc'
// check whether the value types of the arguments are VT_COMPLEX or not
// raise Mat::ValueTypeError if at least one of the arguments is not a complex matrix
void
km_check_complex(int argc, ...)
{
	va_list argp;
	va_start(argp, argc);
	for ( int i=0; i<argc; i++ ) {
		SMAT *smat = va_arg(argp, SMAT *);
		if ( smat->vtype != VT_COMPLEX ) {
			va_end(argp);
			rb_raise(km_eVT, "complex matrix is expected");
		}
	}
	va_end(argp);
}

// the number of SMATs is specified by `argc'
// check whether the value types of the arguments are VT_VALUE or not
// raise Mat::ValueTypeError if at least one of the arguments is not a ruby-object matrix
void
km_check_value(int argc, ...)
{
	va_list argp;
	va_start(argp, argc);
	for ( int i=0; i<argc; i++ ) {
		SMAT *smat = va_arg(argp, SMAT *);
		if ( smat->vtype != VT_VALUE ) {
			va_end(argp);
			rb_raise(km_eVT, "Ruby-object matrix is expected");
		}
	}
	va_end(argp);
}

// check whether LAPACK routine's output `info' is 0 or not
// raise Mat::InternalError if `info' < 0
// raise `klass' with message `mess' if `info' > 0
// if `klass' is Qnil, no exception will raised when `info' > 0
// on entry, `func' is the name of LAPACK routine that outputs `info'
void
km_check_info(int info, VALUE klass, const char *mess, const char *func)
{
	if ( info == 0 ) {
		return;
	} else if ( info > 0 ) {
		if ( klass != Qnil ) {
			rb_raise(klass, "%s", mess);
		} else {
			return;
		}
	} else {
		rb_raise(km_eInternal, "%d-th argument of %s_() is illegal", -info, func);
	}
}

// compute X which satisfy AX=B. `self' is output
VALUE
kmm_mat_solve_destl(VALUE self, VALUE va, VALUE vb)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	VT_SWITCH( smat->vtype,
		return km_dmat_solve(self, va, vb);,
		rb_raise(km_eNotImp, "solve for complex matricies is not available yet");,
		rb_raise(km_eVT, "can't solve int matricies");,
		rb_raise(km_eVT, "can't solve boolean matricies");,
		return km_vmat_solve(self, va, vb);
	);
}
VALUE
kmm_mat_solve(VALUE va, VALUE vb)
{
	SMAT *sb = km_mat2smat(vb);
	return kmm_mat_solve_destl(km_Mat(sb->m, sb->n, sb->vtype), va, vb);
}

// compute X which satisfy A'X=B'. `self' is output
VALUE
km_recover_trans(VALUE vab)
{
	SMAT *sa = km_mat2smat(rb_ary_entry(vab, 0));
	SMAT *sb = km_mat2smat(rb_ary_entry(vab, 1));
	sa->trans = !sa->trans; SWAP(int, sa->m, sa->n);
	sb->trans = !sb->trans; SWAP(int, sb->m, sb->n);
	return vab;
}
static VALUE
km_mat_solve_wrap(VALUE data)
{
	return kmm_mat_solve_destl(rb_ary_entry(data, 0), rb_ary_entry(data, 1), rb_ary_entry(data, 2));
}
VALUE
kmm_mat_tsolve_destl(VALUE self, VALUE va, VALUE vb)
{
	VALUE vab = rb_ary_new3(2, va, vb);
	km_recover_trans(vab);
	km_ensure(km_mat_solve_wrap, rb_ary_new3(3, self, va, vb), km_recover_trans, vab);
	return self;
}
VALUE
kmm_mat_tsolve(VALUE va, VALUE vb)
{
	SMAT *sb = km_mat2smat(vb);
	return kmm_mat_tsolve_destl(km_Mat(sb->n, sb->m, sb->vtype), va, vb);
}

// compute matrix inverse of `self'. `self' is replaced by the result
// alias inv!
VALUE
kmm_mat_inverse_destl(VALUE self, VALUE va)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	VT_SWITCH( smat->vtype,
		return km_dmat_inverse(self, va);,
		rb_raise(km_eNotImp, "inverse for complex matrix is not available yet");,
		rb_raise(km_eVT, "can't inverse int matricies");,
		rb_raise(km_eVT, "can't inverse boolean matricies");,
		return km_vmat_inverse(self, va);
	);
}
// alias inv
VALUE
kmm_mat_inverse(VALUE va)
{
	SMAT *sa = km_mat2smat(va);
	return kmm_mat_inverse_destl(km_Mat(sa->m, sa->m, sa->vtype), va);
}

// matrix product
void
km_dmprod(int m, int n, int k, double *a, int lda, double *b, int ldb, double *r, int ldr)
{
	static char trans[] = "N";
	static double alpha=1.0, beta=0.0;
	dgemm_(trans, trans, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, r, &ldr);
}
struct km_mprod_arg {
	SMAT *sr, *sa, *sb;
	LAWORK r, a, b;
};
static VALUE
km_mprod_body(VALUE data)
{
	struct km_mprod_arg *a = (struct km_mprod_arg *)data;
	
	int m = a->sr->m, n = a->sr->n, k = a->sa->n;
	km_check_size(3, a->sb->m,k, a->sa->m,m, a->sb->n,n);
	
	KALLOCz(a->r, a->sr);
	KALLOCn(a->a, a->sa);
	KALLOCn(a->b, a->sb);
	
	VTYPE vt = a->sr->vtype;
	if ( vt == VT_DOUBLE ) {
		char trans[] = "N";
		double alpha=1.0, beta=0.0;
		dgemm_(trans, trans, &m, &n, &k, &alpha, a->a.d, &(a->a.ld), a->b.d, &(a->b.ld), &beta, a->r.d, &(a->r.ld));
	} else if ( vt == VT_COMPLEX ) {
		char trans[] = "N";
		COMPLEX alpha=cpack(1.0, 0.0), beta=cpack(0.0, 0.0);
		zgemm_(trans, trans, &m, &n, &k, &alpha, a->a.z, &(a->a.ld), a->b.z, &(a->b.ld), &beta, a->r.z, &(a->r.ld));
	} else if ( vt == VT_INT ) {
		for ( int i=0; i<m; i++ ) { for ( int j=0; j<n; j++ ) {
			int *foo=(a->r.i)+(i+j*(a->r.ld));
			int jldb = j*(a->b.ld);
			for ( int l=0; l<k; l++ ) {
				*foo += (a->a.i)[i+l*(a->a.ld)] * (a->b.i)[l+jldb];
			}
		} }
	} else if ( vt == VT_BOOL ) {
		for ( int i=0; i<m; i++ ) { for ( int j=0; j<n; j++ ) {
			bool *foo=(a->r.b)+(i+j*(a->r.ld));
			int jldb = j*(a->b.ld);
			for ( int l=0; l<k; l++ ) {
				bool bar = ( (a->a.b)[i+l*(a->a.ld)] && (a->b.b)[l+jldb] );
				*foo = XOR(*foo, bar);
			}
		} }
	} else if ( vt == VT_VALUE ) {
		for ( int i=0; i<m; i++ ) { for ( int j=0; j<n; j++ ) {
			VALUE *foo=(a->r.v)+(i+j*(a->r.ld));
			*foo = INT2NUM(0);
			int jldb = j*(a->b.ld);
			for ( int l=0; l<k; l++ ) {
				VALUE bar = rb_funcall( (a->a.v)[i+l*(a->a.ld)], id_op_mul, 1, (a->b.v)[l+jldb] );
				*foo = rb_funcall(*foo, id_op_plus, 1, bar);
			}
		} }
	}
	
	return Qnil;
}
static VALUE
km_mprod_ensure(VALUE data)
{
	struct km_mprod_arg *a = (struct km_mprod_arg *)data;
	
	km_free_if_needed(&(a->b));
	km_free_if_needed(&(a->a));
	km_copy_and_free_if_needed(a->sr, &(a->r));
	
	return Qnil;
}
VALUE
kmm_mat_mprod_destl(VALUE self, VALUE va, VALUE vb)
{
	km_check_frozen(self);
	struct km_mprod_arg a;
	a.sr = km_mat2smat(self); a.sa = km_mat2smat(va); a.sb = km_mat2smat(vb);
	if ( a.sr->vtype != a.sa->vtype || a.sr->vtype != a.sb->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	
	km_ensure(km_mprod_body, (VALUE)&a, km_mprod_ensure, (VALUE)&a);
	
	return self;
}
VALUE
kmm_mat_mprod(VALUE va, VALUE vb)
{
	SMAT *sa = km_mat2smat(va), *sb = km_mat2smat(vb);
	return kmm_mat_mprod_destl(km_Mat(sa->m, sb->n, sa->vtype), va, vb);
}
