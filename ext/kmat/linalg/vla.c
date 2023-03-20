#include "../kmat.h"

struct km_vmat_solve_arg {
	SMAT *sx, *sa, *sb;
	VALUE *a;
	LAWORK x;
};
static VALUE
km_vmat_solve_body(VALUE data)
{
	struct km_vmat_solve_arg *a = (struct km_vmat_solve_arg *)data;
	const size_t n=a->sx->m, nrhs=a->sx->n;
	km_check_size_s(4, a->sa->m,n, a->sa->n,n, a->sb->m,n, a->sb->n,nrhs);
	
	KALLOCc(a->a, a->sa);
	km_alloc_if_needed(a->sx, &(a->x));
	km_copy2work(a->x.v, a->x.ld, a->sb);
	const VALUE zero=INT2NUM(0);
	
	// forward elimination
	for ( size_t k=0; k<n; k++ ) {
		if ( rb_funcall((a->a)[k+k*n], id_op_eq, 1, zero) ) {
			for ( size_t i=k+1; i<n; i++ ) {
				if ( !rb_funcall((a->a)[i+k*n], id_op_eq, 1, zero) ) {
					for ( size_t j=k; j<n; j++ ) {
						SWAP(VALUE, (a->a)[k+j*n], (a->a)[i+j*n]);
					}
					for ( size_t j=0; j<nrhs; j++ ) {
						SWAP(VALUE, (a->x.v)[k+j*(a->x.ld)], (a->x.v)[i+j*(a->x.ld)]);
					}
					goto nonsingular;
				}
			}
			rb_raise(km_eUncomp, "matrix is singular");
			nonsingular: ;
		}
		const VALUE akk = (a->a)[k+k*n];
		for ( size_t j=k+1; j<n; j++ ) {
			(a->a)[k+j*n] = rb_funcall((a->a)[k+j*n], id_quo, 1, akk);
		}
		for ( size_t j=0; j<nrhs; j++ ) {
			(a->x.v)[k+j*(a->x.ld)] = rb_funcall((a->x.v)[k+j*(a->x.ld)], id_quo, 1, akk);
		}
		for( size_t i=k+1; i<n; i++ ) {
			const VALUE aik = (a->a)[i+k*n];
			for ( size_t j=k+1; j<n; j++ ) {
				const VALUE tmp = rb_funcall(aik, id_op_mul, 1, (a->a)[k+j*n]);
				(a->a)[i+j*n] = rb_funcall((a->a)[i+j*n], id_op_minus, 1, tmp);
			}
			for ( size_t j=0; j<nrhs; j++ ) {
				const VALUE tmp = rb_funcall(aik, id_op_mul, 1, (a->x.v)[k+j*(a->x.ld)]);
				(a->x.v)[i+j*(a->x.ld)] = rb_funcall((a->x.v)[i+j*(a->x.ld)], id_op_minus, 1, tmp);
			}
		}
	}
	
	// back substitution
	for ( size_t k=n-1; k>0; k-- ) {
		for ( size_t i=0; i<k; i++ ) {
			const VALUE aik = (a->a)[i+k*n];
			for ( size_t j=0; j<nrhs; j++ ) {
				const VALUE tmp = rb_funcall(aik, id_op_mul, 1, (a->x.v)[k+j*(a->x.ld)]);
				(a->x.v)[i+j*(a->x.ld)] = rb_funcall((a->x.v)[i+j*(a->x.ld)], id_op_minus, 1, tmp);
			}
		}
	}
	
	return Qnil;
}
static VALUE
km_vmat_solve_ensure(VALUE data)
{
	struct km_vmat_solve_arg *a = (struct km_vmat_solve_arg *)data;
	
	ruby_xfree(a->a);
	km_copy_and_free_if_needed(a->sx, &(a->x));
	
	return Qnil;
}
VALUE
km_vmat_solve(VALUE self, VALUE va, VALUE vb)
{
	struct km_vmat_solve_arg a; memset(&a, 0, sizeof(a));
	a.sx=km_mat2smat(self); a.sa=km_mat2smat(va); a.sb=km_mat2smat(vb);
	km_check_value(3, a.sx, a.sa, a.sb);
	
	km_ensure(km_vmat_solve_body, (VALUE)(&a), km_vmat_solve_ensure, (VALUE)(&a));
	
	return self;
}

VALUE
km_vmat_inverse(VALUE self, VALUE va)
{
	SMAT *sx = km_mat2smat(self), *sa = km_mat2smat(va);
	km_check_value(2, sx, sa);
	const size_t n = sa->m;
	km_check_size_s(3, sx->m,n, sx->n,n, sa->n,n);
	VALUE ident = km_Mat(n, n, VT_VALUE);
	kmm_mat_eye(ident);
	return km_vmat_solve(self, va, ident);
}
