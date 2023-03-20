#include "../kmat.h"

static inline void
km_check_info_opt(int info, const char *funcname)
{
	km_check_info(info, rb_eRuntimeError, "error occured while computing optimal lwork", funcname);
}

// move elements of the first column of `body' to the diagonal position
// set 0 to the original positions
static void
km_vec2diag(int m, int n, double *body, int ld)
{
	const int lm1 = MIN(m, n)-1, mp1=ld+1; int i=1;
	dcopy_(&lm1, body+1, &i, body+mp1, &mp1);
	for ( ; i<=lm1; i++ ) { body[i] = 0.0; }
}

// restore the modification caused by `dgebal'
// it means, smat->body = D * AA / D
// on entry, AA is `smat' and D can be construct by scale, ilo and ihi
struct km_unbal_arg {
	SMAT *smat;
	double *scale;
	int ilo, ihi;
	int n;
	LAWORK body;
	double *at;
	SMAT *sat;
};
static VALUE
km_unbal_body(VALUE data)
{
	struct km_unbal_arg *a = (struct km_unbal_arg *)data;
	
	int info;
	char job[]="B", side[]="R";
	KALLOC(a->at, (a->n)*(a->n));
	KALLOCn(a->body, a->smat);
	int ld = s2i(a->body.ld);
	dgebak_(job, side, &(a->n), &(a->ilo), &(a->ihi), a->scale, &(a->n), a->body.d, &ld, &info); // multiply D from the left
	km_check_info(info, rb_eRuntimeError, "unexpected info value", "dgebak");
	a->sat = km_smat_alloc_with(i2s(a->n), i2s(a->n), VT_DOUBLE, a->body.body);
	a->sat->trans = true;
	km_copy2work(a->at, i2s(a->n), a->sat); // take transpose
	side[0] = 'L';
	dgebak_(job, side, &(a->n), &(a->ilo), &(a->ihi), a->scale, &(a->n), a->at, &(a->n), &info); // multiply D^(-T) from the left
	km_check_info(info, rb_eRuntimeError, "unexpected info value", "dgebak");
	a->sat->body = a->at; km_copy2work(a->body.body, a->body.ld, a->sat); // take transpose again
	
	return Qnil;
}
static VALUE
km_unbal_ensure(VALUE data)
{
	struct km_unbal_arg *a = (struct km_unbal_arg *)data;
	
	km_copy_and_free_if_needed(a->smat, &(a->body));
	ruby_xfree(a->sat);
	ruby_xfree(a->at);
	
	return Qnil;
}
static void
km_unbal(SMAT *smat, double *scale, int ilo, int ihi)
{
	struct km_unbal_arg a = {smat, scale, ilo, ihi, s2i(smat->m), {{NULL}, 0, false}, NULL, NULL};
	
	km_ensure(km_unbal_body, (VALUE)(&a), km_unbal_ensure, (VALUE)(&a));
}

static void
km_check_finite_func_d(double *ent, size_t i, size_t j, void *null)
{
	if ( !isfinite(*ent) ) {
		rb_raise(km_eUncomp, "the matrix has an illegal (infinite or nan) element at (%zu, %zu)", i, j);
	}
}
static inline void
km_check_finite(SMAT *smat)
{
	km_smat_each_with_index_d(smat, km_check_finite_func_d, NULL);
}

// compute X which satisfies AX=B. `self' is output
struct km_dmat_solve_arg {
	SMAT *sx;
	LAWORK x;
	SMAT *sa, *sb;
	double *a, *af;
	int *ipiv;
	double *r, *c, *b;
	double *ferr, *berr;
	double *work;
	int *iwork;
};
VALUE
km_dmat_solve_body(VALUE data)
{
	struct km_dmat_solve_arg *a = (struct km_dmat_solve_arg *)data;
	
	int n=s2i(a->sx->m), nrhs=s2i(a->sx->n);
	km_check_size(4, s2i(a->sa->m),n, s2i(a->sa->n),n, s2i(a->sb->m),n, s2i(a->sb->n),nrhs);
	
	int info;
	char equed[2];
	char fact[]="E";
	char trans[]="N";
	KALLOCc(a->a, a->sa);
	KALLOC(a->af, n*n);
	KALLOC(a->ipiv, n);
	KALLOC(a->r, n);
	KALLOC(a->c, n);
	KALLOCc(a->b, a->sb);
	KALLOC(a->ferr, nrhs);
	KALLOC(a->berr, nrhs);
	KALLOC(a->work, 4*n);
	KALLOC(a->iwork, n);
	KALLOCn(a->x, a->sx);
	double rcond;
	int ldx=s2i(a->x.ld);
	
	dgesvx_(fact, trans, &n, &nrhs, a->a, &n, a->af, &n, a->ipiv, equed, a->r, a->c, a->b, &n,
		a->x.d, &ldx, &rcond, a->ferr, a->berr, a->work, a->iwork, &info);
	km_check_info(info, km_eUncomp, "A is singular or near singular", "dgesvx");
	
	return Qnil;
}
VALUE
km_dmat_solve_ensure(VALUE data)
{
	struct km_dmat_solve_arg *a = (struct km_dmat_solve_arg *)data;
	
	km_copy_and_free_if_needed(a->sx, &(a->x));	
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	ruby_xfree(a->berr);
	ruby_xfree(a->ferr);
	ruby_xfree(a->b);
	ruby_xfree(a->c);
	ruby_xfree(a->r);
	ruby_xfree(a->ipiv);
	ruby_xfree(a->af);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
km_dmat_solve(VALUE self, VALUE va, VALUE vb)
{
	struct km_dmat_solve_arg a; memset(&a, 0, sizeof(a));
	a.sx = km_mat2smat(self);
	a.sa = km_mat2smat(va);
	a.sb = km_mat2smat(vb);
	km_check_double(3, a.sx, a.sa, a.sb);
	km_check_finite(a.sa); km_check_finite(a.sb);
	
	km_ensure(km_dmat_solve_body, (VALUE)&a, km_dmat_solve_ensure, (VALUE)&a);
	
	return self;
}

// compute the inverse matrix using Mat#solve
VALUE
km_dmat_inverse(VALUE self, VALUE va)
{
	SMAT *sx = km_mat2smat(self), *sa = km_mat2smat(va);
	km_check_double(2, sx, sa);
	km_check_finite(sa);
	const size_t n = sa->m;
	km_check_size_s(3, sx->m,n, sx->n,n, sa->n,n);
	return km_dmat_solve(self, va, kmm_Mat_identity(km_cMat, ZU2NUM(n)));
}

// compute X with the smallest norm which minimize ||AX-B||. `self' is output
struct km_ls_arg {
	SMAT *sx, *sa, *sb;
	double *a, *b, *s, *work;
	int *iwork;
};
static VALUE
km_mat_ls_body(VALUE data)
{
	struct km_ls_arg *a = (struct km_ls_arg *)data;
	int m=s2i(a->sa->m), n=s2i(a->sa->n), nrhs=s2i(a->sx->n);
	int ldb = MAX(m, n);
	km_check_size(3, s2i(a->sx->m),n, s2i(a->sb->m),m, s2i(a->sb->n),nrhs);
	
	int rank, lwork=-1, liwork, info;
	double rcond=-1.0, opt;
	dgelsd_(&m, &n, &nrhs, NULL, &m, NULL, &ldb, NULL, NULL, &rank, &opt, &lwork, &liwork, &info);
	km_check_info_opt(info, "dgelsd");
	lwork = (int)opt;
	KALLOCc(a->a, a->sa);
	KALLOC(a->b, ldb*nrhs); // don't use KALLOCc becaus ldb can be greater than sb->m
	km_copy2work(a->b, i2s(ldb), a->sb);
	KALLOC(a->s, MIN(m, n));
	KALLOC(a->work, lwork);
	KALLOC(a->iwork, liwork);
	dgelsd_(&m, &n, &nrhs, a->a, &m, a->b, &ldb, a->s, &rcond, &rank, a->work, &lwork, a->iwork, &info);
	km_check_info(info, rb_eRuntimeError, "compution the SVD faild to converge", "dgelsd");
	km_copy_from_work(a->sx, a->b, i2s(ldb));
	
	return Qnil;
}
static VALUE
km_mat_ls_ensure(VALUE data)
{
	struct km_ls_arg *a = (struct km_ls_arg *)data;
	
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	ruby_xfree(a->s);
	ruby_xfree(a->b);
	ruby_xfree(a->a);
	
	return Qnil;
}

VALUE
kmm_mat__ls(VALUE self, VALUE va, VALUE vb)
{
	km_check_frozen(self);
	struct km_ls_arg a; memset(&a, 0, sizeof(a));
	a.sx = km_mat2smat(self); a.sa = km_mat2smat(va); a.sb = km_mat2smat(vb);
	km_check_double(3, a.sx, a.sa, a.sb);
	km_check_finite(a.sa); km_check_finite(a.sb);
	
	km_ensure(km_mat_ls_body, (VALUE)&a, km_mat_ls_ensure, (VALUE)&a);
	
	return self;
}

// compute X with the smallest norm which minimize ||A'X-B'||. `self' is output
static VALUE
km_mat_ls_wrap(VALUE data)
{
	return kmm_mat__ls(rb_ary_entry(data, 0), rb_ary_entry(data, 1), rb_ary_entry(data, 2));
}
VALUE
kmm_mat_tls_destl(VALUE self, VALUE va, VALUE vb)
{
	SMAT *sa = km_mat2smat(va), *sb = km_mat2smat(vb);
	sa->trans = !(sa->trans); SWAP(size_t, sa->m, sa->n);
	sb->trans = !(sb->trans); SWAP(size_t, sb->m, sb->n);
	VALUE vab = rb_ary_new3(2, va, vb);
	km_ensure(km_mat_ls_wrap, rb_ary_new3(3, self, va, vb), km_recover_trans, vab);
	return self;
}
VALUE
kmm_mat_tls(VALUE va, VALUE vb)
{
	SMAT *sa = km_mat2smat(va), *sb = km_mat2smat(vb);
	return kmm_mat_tls_destl(km_Mat(sb->n, sa->n, VT_DOUBLE), va, vb);
}

// compute x which minimize ||A'(b-Ax)|| using the conjugate gradient method. `self' is output
struct km_ls_conj_arg {
	SMAT *sx, *sa, *sb;
	double tor;
	double *r, *p, *ap, *aa;
	LAWORK a, b, x;
};
static VALUE
km_mat_ls_conj_body(VALUE data)
{
	struct km_ls_conj_arg *a = (struct km_ls_conj_arg *)data;
	
	const int m = s2i(a->sa->m), n = s2i(a->sa->n);
	if ( m < n ) {
		rb_raise(km_eDim, "A must be row-full-rank");
	}
	km_check_size(4, a->sx->m,n, a->sx->n,1, a->sb->m,m, a->sb->n,1);
	
	const int ione=1; const double dzero=0.0, done=1.0;
	KALLOC(a->r, n);
	KALLOC(a->p, n);
	KALLOC(a->ap, n);
	KALLOC(a->aa, n*n);
	KALLOCn(a->a, a->sa);
	KALLOCn(a->b, a->sb);
	KALLOCz(a->x, a->sx);
	dgemm_("T", "N", &n, &ione, &m, &done, a->a.d, &m, a->b.d, &m, &dzero, a->r, &n); // r = A'b
	dgemm_("T", "N", &n, &n, &m, &done, a->a.d, &m, a->a.d, &m, &dzero, a->aa, &n); // aa = A'A
	
	dcopy_(&n, a->r, &ione, a->p, &ione); // p = r
	double rr = ddot_(&n, a->r, &ione, a->r, &ione); // rr = ||r||^2
	double th = rr*(a->tor)*(a->tor);
	if ( th < DBL_MIN ) {
		th = DBL_MIN;
	}
	double rrp;
	int i=0;
	for (;;) {
		dgemm_("N", "N", &n, &ione, &n,  &done, a->aa, &n, a->p, &n, &dzero, a->ap, &n); // ap = Ap
		double alp = rr/ddot_(&n, a->p, &ione, a->ap, &ione); // alp = ||r||^2/(p'Ap)
		daxpy_(&n, &alp, a->p, &ione, a->x.d, &ione); // x += alp*p
		alp = -alp;
		daxpy_(&n, &alp, a->ap, &ione, a->r, &ione); // r -= alp*Ap
		rrp = ddot_(&n, a->r, &ione, a->r, &ione); // rrp = ||r||^2
		if ( rrp < th ) {
			break;
		} else if ( rr < rrp ) {
			i++;
			if ( n < i ) {
				rb_raise(km_eUncomp, "A may not be positive definite");
			}
		}
		const double bet = rrp/rr;
		for ( int j=0; j<n; j++ ) {
			a->p[j] = bet*a->p[j]+a->r[j];
		}
		// the above is faster than the below
		// dscal_(&n, &bet, a->p, &ione);
		// daxpy_(&n, &done, a->r, &ione, a->p, &ione);
		rr = rrp;
	}
	
	return Qnil;
}
static VALUE
km_mat_ls_conj_ensure(VALUE data)
{
	struct km_ls_conj_arg *a = (struct km_ls_conj_arg *)data;

	km_copy_and_free_if_needed(a->sx, &(a->x));
	km_free_if_needed(&(a->b));
	km_free_if_needed(&(a->a));
	ruby_xfree(a->aa);
	ruby_xfree(a->ap);
	ruby_xfree(a->p);
	ruby_xfree(a->r);

	return Qnil;
}
VALUE
kmm_mat__ls_conj(VALUE self, VALUE va, VALUE vb, VALUE vtor)
{
	km_check_frozen(self);
	struct km_ls_conj_arg a; memset(&a, 0, sizeof(a));
	a.sx = km_mat2smat(self); a.sa = km_mat2smat(va); a.sb = km_mat2smat(vb);
	km_check_double(3, a.sx, a.sa, a.sb);
	km_check_finite(a.sa); km_check_finite(a.sb);
	
	a.tor = NUM2DBL(vtor);
	km_ensure(km_mat_ls_conj_body, (VALUE)&a, km_mat_ls_conj_ensure, (VALUE)&a);
	
	return self;
}


// compute the solution x of minimize_x ||y|| s.t. d=Ax+By. `self' is output
// self: m-vector
// A: (n, m)-matrix
// B: (n, p)-matrix
// d: n-vector
// if B is regular square matrix, it is equivalent to minimize_x ||B\(d-Ax)||
struct km_glm_arg {
	SMAT *sx, *sa, *sb, *sd;
	double *a, *b, *d, *y, *work;
	LAWORK x;
};
static VALUE
km_mat_glm_body(VALUE data)
{
	struct km_glm_arg *a = (struct km_glm_arg *)data;
	
	int n=s2i(a->sa->m), m=s2i(a->sa->n), p=s2i(a->sb->n); // m, n are swapped from those of sa
	if ( n < m || m+p < n ) {
		rb_raise(km_eDim, "m <= n <= m+p must be satisfied for glm, given are (m, n, p) = (%d, %d, %d)", m, n, p);
	}
	km_check_size(5, s2i(a->sb->m),n, s2i(MIN(a->sx->m,a->sx->n)),1, LENGTHi(a->sx),m, s2i(MIN(a->sd->m, a->sd->n)),1, LENGTHi(a->sd),n);
	
	double opt; int lwork=-1, info;
	dggglm_(&n, &m, &p, NULL, &n, NULL, &n, NULL, NULL, NULL, &opt, &lwork, &info);
	km_check_info_opt(info, "dggglm");
	lwork = (int)opt;
	
	KALLOCc(a->a, a->sa);
	KALLOCc(a->b, a->sb);
	KALLOCc(a->d, a->sd);
	KALLOC(a->y, p);
	KALLOC(a->work, lwork);
	KALLOCn(a->x, a->sx);
	
	dggglm_(&n, &m, &p, a->a, &n, a->b, &n, a->d, a->x.d, a->y, a->work, &lwork, &info);
	km_check_info(info, km_eUncomp, "pair (A, B) is not full-rank", "dgglm");
	
	return Qnil;
}
static VALUE
km_mat_glm_ensure(VALUE data)
{
	struct km_glm_arg *a = (struct km_glm_arg *)data;

	km_copy_and_free_if_needed(a->sx, &(a->x));
	ruby_xfree(a->work);
	ruby_xfree(a->y);
	ruby_xfree(a->d);
	ruby_xfree(a->b);
	ruby_xfree(a->a);

	return Qnil;
}
VALUE
kmm_mat_glm_destl(VALUE self, VALUE va, VALUE vb, VALUE vd)
{
	km_check_frozen(self);
	struct km_glm_arg a; memset(&a, 0, sizeof(a));
	a.sx = km_mat2smat(self); a.sa = km_mat2smat(va); a.sb = km_mat2smat(vb); a.sd = km_mat2smat(vd);
	km_check_double(4, a.sx, a.sa, a.sb, a.sd);
	km_check_finite(a.sa); km_check_finite(a.sb); km_check_finite(a.sd);
	
	km_ensure(km_mat_glm_body, (VALUE)&a, km_mat_glm_ensure, (VALUE)&a);
	
	return self;
}
VALUE
kmm_mat_glm(VALUE va, VALUE vb, VALUE vd)
{
	return kmm_mat_glm_destl(km_Mat(km_mat2smat(va)->n, 1, VT_DOUBLE), va, vb, vd);
}

// compute the eigenvalues of a symmetric matrix A. `self' is output
struct km_sym_ev_arg {
	SMAT *sd, *sa;
	double *a, *work;
	int *iwork;
	LAWORK w;
};
static VALUE
km_sym_ev_body(VALUE data)
{
	struct km_sym_ev_arg *a = (struct km_sym_ev_arg *)data;
	
	int n = LENGTHi(a->sd);
	km_check_size(3, s2i(MIN(a->sd->m, a->sd->n)),1, s2i(a->sa->m),n, s2i(a->sa->n),n);
	
	char cmach[] = "S";
	double abstol = dlamch_(cmach);
	int lwork=-1, liwork=-1;
	double dopt; int m, iopt, info;
	char jobz[]="N", range[] = "A", upto[] = "U";
	dsyevr_(jobz, range, upto, &n, NULL, &n, NULL, NULL, NULL, NULL, &abstol, &m, NULL, NULL, &n, NULL, &dopt, &lwork, &iopt, &liwork, &info);
	km_check_info_opt(info, "dsyevr");
	lwork = (int)dopt; liwork = iopt;
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->work, lwork);
	KALLOC(a->iwork, liwork);
	KALLOCn(a->w, a->sd);
	dsyevr_(jobz, range, upto, &n, a->a, &n, NULL, NULL, NULL, NULL, &abstol, &m, a->w.d, NULL, &n, NULL, a->work, &lwork, a->iwork, &liwork, &info);
	km_check_info(info, rb_eRuntimeError, "internal error occured while invoking dsyevr", "dsyevr");
	
	return Qnil;
}
static VALUE
km_sym_ev_ensure(VALUE data)
{
	struct km_sym_ev_arg *a = (struct km_sym_ev_arg *)data;
	
	km_copy_and_free_if_needed(a->sd, &(a->w));
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	ruby_xfree(a->a);

	return Qnil;
}
VALUE
kmm_mat_sym_eigen_values_destl(VALUE self, VALUE va)
{
	km_check_frozen(self);
	struct km_sym_ev_arg a; memset(&a, 0, sizeof(a));
	a.sd = km_mat2smat(self); a.sa = km_mat2smat(va);
	km_check_double(2, a.sd, a.sa);
	km_check_finite(a.sa);
	
	km_ensure(km_sym_ev_body, (VALUE)&a, km_sym_ev_ensure, (VALUE)&a);
	
	return self;
}
VALUE
kmm_mat_sym_eigen_values(VALUE va)
{
	return kmm_mat_sym_eigen_values_destl(km_Mat(km_mat2smat(va)->m, 1, VT_DOUBLE), va);
}

// invoke eigen-decomposition A=VDV' of a symmetric matrix A
// `self' is A and the arguments are outputs
struct km_sym_evd_arg {
	SMAT *sa, *sv, *sd;
	double *a, *work;
	int *isuppz, *iwork;
	LAWORK w, z;
};
static VALUE
km_sym_evd_body(VALUE data)
{
	struct km_sym_evd_arg *a = (struct km_sym_evd_arg *)data;
	
	int n = s2i(a->sa->m);
	km_check_size(5, s2i(a->sa->n),n, s2i(a->sv->m),n, s2i(a->sv->n),n, s2i(a->sd->m),n, s2i(a->sd->n),n);
	
	char cmach[] = "S";
	double dopt, abstol = dlamch_(cmach);
	int m, lwork=-1, liwork=-1, iopt, info;
	char jobz[] = "V", range[] = "A", upto[] = "U";
	dsyevr_(jobz, range, upto, &n, NULL, &n, NULL, NULL, NULL, NULL, &abstol, &m, NULL, NULL, &n, NULL, &dopt, &lwork, &iopt, &liwork, &info);
	km_check_info_opt(info, "dsyevr");
	lwork = (int)dopt; liwork = iopt;
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->isuppz, 2*n);
	KALLOC(a->work, lwork);
	KALLOC(a->iwork, liwork);
	KALLOCz(a->w, a->sd);
	KALLOCn(a->z, a->sv);
	int ldz=s2i(a->z.ld);
	
	dsyevr_(jobz, range, upto, &n, a->a, &n, NULL, NULL, NULL, NULL, &abstol, &m,
		a->w.d, a->z.d, &ldz, a->isuppz, a->work, &lwork, a->iwork, &liwork, &info);
	km_check_info(info, rb_eRuntimeError, "internal error occured while invoking dsyevr", "dsyevr");
	km_vec2diag(n, n, a->w.d, s2i(a->w.ld));
	
	return Qnil;
}
static VALUE
km_sym_evd_ensure(VALUE data)
{
	struct km_sym_evd_arg *a = (struct km_sym_evd_arg *)data;
	
	km_copy_and_free_if_needed(a->sv, &(a->z));
	km_copy_and_free_if_needed(a->sd, &(a->w));
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	ruby_xfree(a->isuppz);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat_sym_evd_destl(VALUE self, VALUE vv, VALUE vd)
{
	km_check_frozen(vv); km_check_frozen(vd);
	struct km_sym_evd_arg a; memset(&a, 0, sizeof(a));
	a.sa = km_mat2smat(self); a.sv = km_mat2smat(vv); a.sd = km_mat2smat(vd);
	km_check_double(3, a.sa, a.sv, a.sd);
	km_check_finite(a.sa);
	
	km_ensure(km_sym_evd_body, (VALUE)&a, km_sym_evd_ensure, (VALUE)&a);
	
	return Qnil;
}
VALUE
kmm_mat_sym_evd(VALUE self)
{
	const size_t n = km_mat2smat(self)->n;
	VALUE vv = km_Mat(n, n, VT_DOUBLE);
	VALUE vd = km_Mat(n, n, VT_DOUBLE);
	kmm_mat_sym_evd_destl(self, vv, vd);
	return rb_ary_new3(2, vv, vd);
}

// compute eigenvalues of a non-symmetric matrix A. the output `self' is a complex matrix
struct km_ge_eigen_values_arg {
	SMAT *sd, *sa;
	double *a, *wr, *wi, *scale, *work;
};
static void
km_ge_ev_cpack(COMPLEX *ent, size_t i, size_t j, void *data)
{
	struct km_ge_eigen_values_arg *a = (struct km_ge_eigen_values_arg *)data;
	*ent = cpack(a->wr[i+j], a->wr[i+j]); // i==0 for column-vector or j==0 for row-vector
}
static VALUE
km_ge_eigen_values_body(VALUE data)
{
	struct km_ge_eigen_values_arg *a = (struct km_ge_eigen_values_arg *)data;
	
	int n = LENGTHi(a->sd);
	km_check_size(3, s2i(MIN(a->sd->m, a->sd->n)),1, s2i(a->sa->m),n, s2i(a->sa->n),n);
	
	double opt;
	int lwork=-1, ilo, ihi, info;
	char balanc[]="B", jobvl[]="N", jobvr[]="N", sense[]="N";
	dgeevx_(balanc, jobvl, jobvr, sense, &n, NULL, &n, NULL, NULL, NULL, &n, NULL, &n, &ilo, &ihi, NULL, NULL, NULL, NULL, &opt, &lwork, NULL, &info);
	km_check_info_opt(info, "dgeevx");
	lwork = (int)opt;
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->wr, n);
	KALLOC(a->wi, n);
	KALLOC(a->scale, n);
	KALLOC(a->work, lwork);
	
	double abnrm;
	dgeevx_(balanc, jobvl, jobvr, sense, &n, a->a, &n, a->wr, a->wi, NULL, &n, NULL, &n, &ilo, &ihi, a->scale, &abnrm, NULL, NULL, a->work, &lwork, NULL, &info);
	km_check_info(info, rb_eRuntimeError, "the QR algorithm failed to compute all the eigenvalues", "dgeevx");
	km_smat_each_with_index_z(a->sd, km_ge_ev_cpack, a);
	
	return Qnil;
}
static VALUE
km_ge_eigen_values_ensure(VALUE data)
{
	struct km_ge_eigen_values_arg *a = (struct km_ge_eigen_values_arg *)data;
	
	ruby_xfree(a->work);
	ruby_xfree(a->scale);
	ruby_xfree(a->wi);
	ruby_xfree(a->wr);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat_ge_eigen_values_destl(VALUE self, VALUE va)
{
	km_check_frozen(self);
	struct km_ge_eigen_values_arg a; memset(&a, 0, sizeof(a));
	a.sd = km_mat2smat(self); a.sa = km_mat2smat(va);
	km_check_complex(1, a.sd); km_check_double(1, a.sa);
	km_check_finite(a.sa);
	
	km_ensure(km_ge_eigen_values_body, (VALUE)&a, km_ge_eigen_values_ensure, (VALUE)&a);

	return self;
}
VALUE
kmm_mat_ge_eigen_values(VALUE va)
{
	return kmm_mat_ge_eigen_values_destl(km_Mat(km_mat2smat(va)->m, 1, VT_COMPLEX), va);
}

// compute a matrix consists of right-eigenvectors V and a diagonal matrix consists of right-eigenvalues D, AV=VD of a non-symmetric matrix A.
// the arguments are outputs
struct km_ge_evd_arg {
	SMAT *sa, *sv, *sd;
	double *a, *wr, *wi, *scale, *work, *vr;
};
static VALUE
km_ge_evd_body(VALUE data)
{
	struct km_ge_evd_arg *a = (struct km_ge_evd_arg *)data;
	
	const size_t n_s = a->sa->m;
	int n = s2i(n_s);
	km_check_size(5, s2i(a->sa->n),n, s2i(a->sv->m),n, s2i(a->sv->n),n, s2i(a->sd->m),n, s2i(a->sd->n),n);
	
	double opt; int lwork=-1, ilo, ihi, info;
	char balanc[]="B", jobvl[]="N", jobvr[]="V", sense[]="N";
	dgeevx_(balanc, jobvl, jobvr, sense, &n, NULL, &n, NULL, NULL, NULL, &n, NULL, &n,
		&ilo, &ihi, NULL, NULL, NULL, NULL, &opt, &lwork, NULL, &info);
	km_check_info_opt(info, "dgeevx");
	lwork = (int)opt;
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->wr, n);
	KALLOC(a->wi, n);
	KALLOC(a->scale, n);
	KALLOC(a->work, lwork);
	KALLOC(a->vr, n*n);
	
	double abnrm;
	dgeevx_(balanc, jobvl, jobvr, sense, &n, a->a, &n, a->wr, a->wi, NULL, &n, a->vr, &n,
		&ilo, &ihi, a->scale, &abnrm, NULL, NULL, a->work, &lwork, NULL, &info);
	km_check_info(info, rb_eRuntimeError, "the QR algorithm failed to compute all the eigenvalues", "dgeevx");
	if ( a->sd->stype == ST_RSUB ) {
		for ( size_t j=0; j<n_s; j++ ) {
			if (a->wi[j] == 0.0) {
				for ( size_t i=0; i<n_s; i++ ) {
					ENTITYr0(a->sv, z, INDEX(a->sv, i, j)) = cpack(a->vr[i+j*n_s], 0.0);
				}
			} else if (j == 0) {
				for ( size_t i=0; i<n_s; i++ ) {
					ENTITYr0(a->sv, z, INDEX(a->sv, i, j)) = cpack(a->vr[i], a->vr[i+n_s]);
				}
			} else if (a->wr[j] == a->wr[j+1] && a->wi[j] == -a->wi[j+1]) {
				for ( size_t i=0; i<n_s; i++ ) {
					ENTITYr0(a->sv, z, INDEX(a->sv, i, j)) = cpack(a->vr[i+j*n_s], a->vr[i+(j+1)*n_s]);
				}
			} else {
				for ( size_t i=0; i<n_s; i++ ) {
					ENTITYr0(a->sv, z, INDEX(a->sv, i, j)) = cpack(a->vr[i+(j-1)*n_s], -a->vr[i+j*n_s]);
				}
			}
		}
	} else {
		for ( size_t j=0; j<n_s; j++ ) {
			if (a->wi[j] == 0.0) {
				for ( size_t i=0; i<n_s; i++ ) {
					ENTITYd0(a->sv, z, INDEX(a->sv, i, j)) = cpack(a->vr[i+j*n_s], 0.0);
				}
			} else if (j == 0) {
				for ( size_t i=0; i<n_s; i++ ) {
					ENTITYd0(a->sv, z, INDEX(a->sv, i, j)) = cpack(a->vr[i], a->vr[i+n_s]);
				}
			} else if (a->wr[j] == a->wr[j+1] && a->wi[j] == -a->wi[j+1]) {
				for ( size_t i=0; i<n_s; i++ ) {
					ENTITYd0(a->sv, z, INDEX(a->sv, i, j)) = cpack(a->vr[i+j*n_s], a->vr[i+(j+1)*n_s]);
				}
			} else {
				for ( size_t i=0; i<n_s; i++ ) {
					ENTITYd0(a->sv, z, INDEX(a->sv, i, j)) = cpack(a->vr[i+(j-1)*n_s], -a->vr[i+j*n_s]);
				}
			}
		}
	}
	if ( a->sd->stype == ST_RSUB ) {
		for ( size_t i=0; i<n_s; i++ ) {
			ENTITYr0(a->sd, z, i+i*(a->sd->ld)) = cpack(a->wr[i], a->wi[i]);
		}
	} else {
		for ( size_t i=0; i<n_s; i++ ) {
			ENTITYd0(a->sd, z, i+i*(a->sd->ld)) = cpack(a->wr[i], a->wi[i]);
		}
	}
	
	return Qnil;
}
static VALUE
km_ge_evd_ensure(VALUE data)
{
	struct km_ge_evd_arg *a = (struct km_ge_evd_arg *)data;
	
	ruby_xfree(a->vr);
	ruby_xfree(a->work);
	ruby_xfree(a->scale);
	ruby_xfree(a->wi);
	ruby_xfree(a->wr);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat_ge_evd_destl(VALUE self, VALUE vv, VALUE vd)
{
	km_check_frozen(vv); km_check_frozen(vd);
	struct km_ge_evd_arg a; memset(&a, 0, sizeof(a));
	a.sa = km_mat2smat(self); a.sv = km_mat2smat(vv); a.sd = km_mat2smat(vd);
	km_check_double(1, a.sa, a.sv); km_check_complex(2, a.sv, a.sd);
	km_check_finite(a.sa);
	
	kmm_mat_zero(vd);
	km_ensure(km_ge_evd_body, (VALUE)&a, km_ge_evd_ensure, (VALUE)&a);
	
	return self;
}
VALUE
kmm_mat_ge_evd(VALUE self)
{
	const size_t n = km_mat2smat(self)->n;
	VALUE vv = km_Mat(n, n, VT_COMPLEX);
	VALUE vd = km_Mat(n, n, VT_COMPLEX);
	kmm_mat_ge_evd_destl(self, vv, vd);
	return rb_ary_new3(2, vv, vd);
}

// A の特異値を計算し，self に格納する
// compute singular values of A. `self' is output
struct km_singular_values_arg {
	SMAT *ss, *sa;
	double *a, *work;
	int *iwork;
	LAWORK s;
};
static VALUE
km_singular_values_body(VALUE data)
{
	struct km_singular_values_arg *a = (struct km_singular_values_arg *)data;
	
	int m = s2i(a->sa->m), n = s2i(a->sa->n);
	km_check_size(2, LENGTHi(a->ss),MIN(m, n), s2i(MIN(a->ss->m,a->ss->n)),1);
	
	double opt; int lwork=-1, info;
	char jobz[] = "N";
	dgesdd_(jobz, &m, &n, NULL, &m, NULL, NULL, &m, NULL, &n, &opt, &lwork, NULL, &info);
	km_check_info_opt(info, "dgesdd");
	lwork = (int)opt;
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->work, lwork);
	KALLOC(a->iwork, 8*MIN(m, n));
	KALLOCn(a->s, a->ss);
	
	dgesdd_(jobz, &m, &n, a->a, &m, a->s.d, NULL, &m, NULL, &n, a->work, &lwork, a->iwork, &info);
	km_check_info(info, rb_eRuntimeError, "DBDSDC did not converge", "dgesdd");
	
	return Qnil;
}
static VALUE
km_singular_values_ensure(VALUE data)
{
	struct km_singular_values_arg *a = (struct km_singular_values_arg *)data;
	
	km_copy_and_free_if_needed(a->ss, &(a->s));
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat_singular_values_destl(VALUE self, VALUE va)
{
	km_check_frozen(self);
	struct km_singular_values_arg a; memset(&a, 0, sizeof(a));
	a.ss = km_mat2smat(self); a.sa = km_mat2smat(va);
	km_check_double(2, a.ss, a.sa);
	km_check_finite(a.sa);
	
	km_ensure(km_singular_values_body, (VALUE)&a, km_singular_values_ensure, (VALUE)&a);
	
	return self;
}
VALUE
kmm_mat_singular_values(VALUE va)
{
	SMAT *sa = km_mat2smat(va);
	return kmm_mat_singular_values_destl(km_Mat(MIN(sa->m, sa->n), 1, VT_DOUBLE), va);
}

// invoke singular value decomposition A=USV'
// `self' is A and the arguments are the outputs U, S and V
struct km_svd_arg {
	SMAT *sa, *ss, *su, *sv;
	double *a, *work;
	int *iwork;
	LAWORK s, u, vt;
};
static VALUE
km_svd_body(VALUE data)
{
	struct km_svd_arg *a = (struct km_svd_arg *)data;
	
	int m = s2i(a->sa->m), n = s2i(a->sa->n);
	km_check_size(6, s2i(a->su->m),m, s2i(a->su->n),m, s2i(a->ss->m),m, s2i(a->ss->n),n, s2i(a->sv->m),n, s2i(a->sv->n),n);
	
	double opt; int lwork=-1, info;
	char jobz[] = "A";
	dgesdd_(jobz, &m, &n, NULL, &m, NULL, NULL, &m, NULL, &n, &opt, &lwork, NULL, &info);
	km_check_info_opt(info, "dgesdd");
	lwork = (int)opt;
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->work, lwork);
	KALLOC(a->iwork, 8*MIN(m, n));
	KALLOCz(a->s, a->ss);
	KALLOCn(a->u, a->su);
	KALLOCn(a->vt, a->sv);
	int ldu=s2i(a->u.ld), ldvt=s2i(a->vt.ld);
	
	dgesdd_(jobz, &m, &n, a->a, &m, a->s.d, a->u.d, &ldu, a->vt.d, &ldvt, a->work, &lwork, a->iwork, &info);
	km_check_info(info, rb_eRuntimeError, "DBDSDC did not converge", "dgesvd");
	km_vec2diag(m, n, a->s.d, s2i(a->s.ld));
	
	return Qnil;
}
static VALUE
km_svd_ensure(VALUE data)
{
	struct km_svd_arg *a = (struct km_svd_arg *)data;
	
	km_copy_and_free_if_needed(a->sv, &(a->vt));
	a->sv->trans = !(a->sv->trans);
	km_copy_and_free_if_needed(a->su, &(a->u));
	km_copy_and_free_if_needed(a->ss, &(a->s));
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat_svd_destl(VALUE self, VALUE vu, VALUE vs, VALUE vv)
{
	km_check_frozen(vu); km_check_frozen(vs); km_check_frozen(vv);
	struct km_svd_arg a; memset(&a, 0, sizeof(a));
	a.sa = km_mat2smat(self); a.su = km_mat2smat(vu); a.ss = km_mat2smat(vs); a.sv = km_mat2smat(vv);
	km_check_double(4, a.sa, a.su, a.ss, a.sv);
	km_check_finite(a.sa);
	
	km_ensure(km_svd_body, (VALUE)&a, km_svd_ensure, (VALUE)&a);
	
	return self;
}
VALUE
kmm_mat_svd(VALUE self)
{
	SMAT *sa = km_mat2smat(self);
	const size_t m = sa->m, n = sa->n;
	VALUE vu = km_Mat(m, m, VT_DOUBLE);
	VALUE vs = km_Mat(m, n, VT_DOUBLE);
	VALUE vv = km_Mat(n, n, VT_DOUBLE);
	kmm_mat_svd_destl(self, vu, vs, vv);
	return rb_ary_new3(3, vu, vs, vv);
}

// symmetrize keeping AA' identical, where `self' is A.
// this is using singular vale decompsition A=USV' and returns USU'.
struct km_svd_symmetrize_arg {
	SMAT *sa;
	double *work, *u, *s, *vt;
	int *iwork;
	LAWORK a;
};
VALUE km_svd_symmetrize_body(VALUE data)
{
	struct km_svd_symmetrize_arg *a = (struct km_svd_symmetrize_arg *)data;
	
	int n = s2i(a->sa->m);
	km_check_size(1, s2i(a->sa->n),n);
	
	double opt; int lwork=-1, info;
	char jobz[] = "A";
	dgesdd_(jobz, &n, &n, NULL, &n, NULL, NULL, &n, NULL, &n, &opt, &lwork, NULL, &info);
	km_check_info_opt(info, "dgesdd");
	lwork = (int)opt;
	
	KALLOCn(a->a, a->sa);
	KALLOC(a->work, lwork);
	KALLOC(a->iwork, 8*n);
	KALLOC(a->s, n);
	KALLOC(a->u, n*n);
	KALLOC(a->vt, n*n);
	int lda = s2i(a->a.ld);
	
	dgesdd_(jobz, &n, &n, a->a.d, &lda, a->s, a->u, &n, a->vt, &n, a->work, &lwork, a->iwork, &info);
	km_check_info(info, rb_eRuntimeError, "DBDSDC did not converge", "dgesvd");
	
	int one=1;
	memset(a->vt, 0, sizeof(double)*(size_t)(n*n));
	for (int i=0; i<n; i++) {
		daxpy_(&n, a->s+i, a->u+i*n, &one, a->vt+i*n, &one);
	}
	char ta[] = "N";
	char tb[] = "T";
	double alpha = 1.0, beta=0.0;
	dgemm_(ta, tb, &n, &n, &n, &alpha, a->vt, &n, a->u, &n, &beta, a->a.d, &lda);
	
	return Qnil;
}
VALUE km_svd_symmetrize_ensure(VALUE data)
{
	struct km_svd_symmetrize_arg *a = (struct km_svd_symmetrize_arg *)data;
	
	ruby_xfree(a->vt);
	ruby_xfree(a->u);
	ruby_xfree(a->s);
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	km_copy_and_free_if_needed(a->sa, &(a->a));
	
	return Qnil;
}
VALUE
kmm_mat_svd_symmetrize_dest(VALUE self)
{
	km_check_frozen(self);
	struct km_svd_symmetrize_arg a;
	a.sa = km_mat2smat(self);
	km_check_double(1, a.sa);
	km_check_finite(a.sa);
	
	km_ensure(km_svd_symmetrize_body, (VALUE)&a, km_svd_symmetrize_ensure, (VALUE)&a);
	
	return kmm_mat_symmetrize_dest(self);
}

// invoke a LU decomposition A=LU. `self' is A
// L is a permutated lower triangular matrix and U is a upper triangular matrix
struct km_lu_arg {
	SMAT *sa, *spl, *su;
	double *a;
	int *ipiv, *perm;
	LAWORK pl, u;
};
static VALUE
km_lu_body(VALUE data)
{
	struct km_lu_arg *a = (struct km_lu_arg *)data;
	
	int m = s2i(a->sa->m), n = s2i(a->sa->n);
	const int k = MIN(m, n);
	km_check_size(4, s2i(a->spl->m),m, s2i(a->spl->n),k, s2i(a->su->m),k, s2i(a->su->n),n);
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->ipiv, k);
	KALLOC(a->perm, m);
	KALLOCn(a->pl, a->spl);
	KALLOCn(a->u, a->su);
	for (int i=0; i<m; i++ ) {
		(a->perm)[i] = i;
	}
	
	int info;
	dgetrf_(&m, &n, a->a, &m, a->ipiv, &info);
	km_check_info(info, Qnil, NULL, "dgetrf");
	for ( int i=0; i<k; i++ ) {
		int s = (a->ipiv)[i]-1;
		if ( s != i ) {
			SWAP(int, (a->perm)[i], (a->perm)[s]);
		}
	}
	const int ldpl=s2i(a->pl.ld), ldu = s2i(a->u.ld);
	for ( int i=0; i<k; i++ ) {
		int s = (a->perm)[i];
		dcopy_(&i, (a->a)+i, &m, (a->pl.d)+s, &ldpl);
		(a->pl.d)[s+i*ldpl] = 1.0;
		int izero = 0, len = k-i-1; double dzero = 0.0;
		dcopy_(&len, &dzero, &izero, (a->pl.d)+(s+(i+1)*s2i(a->pl.ld)), &ldpl);
		dcopy_(&i, &dzero, &izero, (a->u.d)+i, &ldu);
		len = n-i;
		dcopy_(&len, (a->a)+(i+i*m), &m, (a->u.d)+(i+i*ldu), &ldu);
	}
	for ( int i=m; i<k; i++ ) {
		dcopy_(&n, a->a+i, &m, (a->pl.d)+(a->perm)[i], &ldpl);
	}
	
	return Qnil;
}
static VALUE
km_lu_ensure(VALUE data)
{
	struct km_lu_arg *a = (struct km_lu_arg *)data;
	
	km_copy_and_free_if_needed(a->su, &(a->u));
	km_copy_and_free_if_needed(a->spl, &(a->pl));
	ruby_xfree(a->perm);
	ruby_xfree(a->ipiv);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat_lu_destl(VALUE self, VALUE vpl, VALUE vu)
{
	km_check_frozen(vpl); km_check_frozen(vu);
	struct km_lu_arg a; memset(&a, 0, sizeof(a));
	a.sa = km_mat2smat(self); a.spl = km_mat2smat(vpl); a.su = km_mat2smat(vu);
	km_check_double(3, a.sa, a.spl, a.su);
	km_check_finite(a.sa);
	
	km_ensure(km_lu_body, (VALUE)&a, km_lu_ensure, (VALUE)&a);
	
	return Qnil;
}
VALUE
kmm_mat_lu(VALUE self)
{
	SMAT *sa = km_mat2smat(self);
	const size_t m = sa->m, n = sa->n;
	const size_t k = MIN(m, n);
	VALUE vpl = km_Mat(m, k, VT_DOUBLE);
	VALUE vu = km_Mat(k, n, VT_DOUBLE);
	kmm_mat_lu_destl(self, vpl, vu);
	return rb_ary_new3(2, vpl, vu);
}

// invoke a LUP decomposition PA=LU. `self' is A
// L is a lower triangular matrix, U is a upper triangular matrix and P is a permutation matrix
struct km_lup_arg {
	SMAT *sa, *sl, *su, *sp;
	double *a;
	int *ipiv, *perm;
	LAWORK p, l, u;
};
static VALUE
km_lup_body(VALUE data)
{
	struct km_lup_arg *a = (struct km_lup_arg *)data;
	
	int m = s2i(a->sa->m), n = s2i(a->sa->n);
	const int k = MIN(m, n);
	km_check_size(6, s2i(a->sl->m),m, s2i(a->sl->n),k, s2i(a->su->m),k, s2i(a->su->n),n, s2i(a->sp->m),m, s2i(a->sp->n),m);
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->ipiv, k);
	KALLOC(a->perm, m);
	KALLOCz(a->p, a->sp);
	KALLOCn(a->l, a->sl);
	KALLOCn(a->u, a->su);
	for ( int i=0; i<m; i++ ) {
		(a->perm)[i] = i;
	}
	
	int info;
	dgetrf_(&m, &n, a->a, &m, a->ipiv, &info);
	km_check_info(info, Qnil, NULL, "dgetrf");
	for ( int i=0; i<k; i++ ) {
		const int s = (a->ipiv)[i]-1;
		if ( s != i ) {
			SWAP(int, (a->perm)[i], (a->perm)[s]);
		}
	}
	const int ldl = s2i(a->l.ld), ldu = s2i(a->u.ld);
	for ( int i=0; i<k; i++ ) {
		(a->p.d)[i+(a->perm)[i]*s2i(a->p.ld)] = 1.0;
		dcopy_(&i, (a->a)+i, &m, (a->l.d)+i, &ldl);
		(a->l.d)[i+i*ldl] = 1.0;
		int len = k-i-1, izero = 0; const double dzero = 0.0;
		dcopy_(&len, &dzero, &izero, (a->l.d)+(i+(i+1)*ldl), &ldl);
		dcopy_(&i, &dzero, &izero, (a->u.d)+i, &ldu);
		len = n-i;
		dcopy_(&len, (a->a)+(i+i*m), &m, (a->u.d)+(i+i*ldu), &ldu);
	}
	for ( int i=m; i<k; i++ ) {
		(a->p.d)[i+(a->perm)[i]*s2i(a->p.ld)] = 1.0;
		dcopy_(&n, (a->a)+i, &m, (a->l.d)+i, &ldl);
	}
	
	return Qnil;
}
static VALUE
km_lup_ensure(VALUE data)
{
	struct km_lup_arg *a = (struct km_lup_arg *)data;
	
	km_copy_and_free_if_needed(a->su, &(a->u));
	km_copy_and_free_if_needed(a->sl, &(a->l));
	km_copy_and_free_if_needed(a->sp, &(a->p));
	ruby_xfree(a->perm);
	ruby_xfree(a->ipiv);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat_lup_destl(VALUE self, VALUE vl, VALUE vu, VALUE vp)
{
	km_check_frozen(vl); km_check_frozen(vu); km_check_frozen(vp);
	struct km_lup_arg a; memset(&a, 0, sizeof(a));
	a.sa = km_mat2smat(self); a.sl = km_mat2smat(vl);
	a.su = km_mat2smat(vu); a.sp = km_mat2smat(vp);
	km_check_double(4, a.sa, a.sl, a.su, a.sp);
	km_check_finite(a.sa);
	
	km_ensure(km_lup_body, (VALUE)&a, km_lup_ensure, (VALUE)&a);
	
	return Qnil;
}
VALUE
kmm_mat_lup(VALUE self)
{
	SMAT *sa = km_mat2smat(self);
	const size_t m = sa->m, n = sa->n; const size_t k = MIN(m, n);
	VALUE vl = km_Mat(m, k, VT_DOUBLE);
	VALUE vu = km_Mat(k, n, VT_DOUBLE);
	VALUE vp = km_Mat(m, m, VT_DOUBLE);
	kmm_mat_lup_destl(self, vl, vu, vp);
	return rb_ary_new3(3, vl, vu, vp);
}

// compute the determinant via LU decomposition
struct km_det_arg {
	SMAT *sa;
	double *a;
	int *ipiv;
};
static VALUE
km_det_body(VALUE data)
{
	struct km_det_arg *a = (struct km_det_arg *)data;
	
	int n = s2i(a->sa->m);
	km_check_size(1, s2i(a->sa->n),n);
	if ( n == 0 ) { return rb_float_new(1.0); }
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->ipiv, n);
	int info;
	dgetrf_(&n, &n, a->a, &n, a->ipiv, &info);
	km_check_info(info, Qnil, NULL, "dgetrf");
	bool neg = ((a->ipiv)[0] != 1);
	for ( int i=1; i<n; i++ ) {
		(a->a)[0] *= (a->a)[i+i*n];
		if ( (a->ipiv)[i]-1 != i ) {
			neg = !neg;
		}
	}
	
	return rb_float_new( neg ? -(a->a)[0] : (a->a)[0] );
}
static VALUE
km_det_ensure(VALUE data)
{
	struct km_det_arg *a = (struct km_det_arg *)data;
	
	ruby_xfree(a->ipiv);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat_det(VALUE self)
{
	struct km_det_arg a; memset(&a, 0, sizeof(a));
	a.sa = km_mat2smat(self);
	km_check_double(1, a.sa);
	km_check_finite(a.sa);
	
	return km_ensure(km_det_body, (VALUE)&a, km_det_ensure, (VALUE)&a);
}

// invoke a QR decomposition A=QR. the arguments are outputs
struct km_qr_arg {
	SMAT *sa, *sq, *sr;
	double *a, *tau, *work;
	LAWORK q, r;
};
static VALUE
km_qr_body(VALUE data)
{
	struct km_qr_arg *a = (struct km_qr_arg *)data;
	
	int m = s2i(a->sa->m), n = s2i(a->sa->n);
	int k = MIN(m, n);
	km_check_size(4, s2i(a->sq->m),m, s2i(a->sq->n),m, s2i(a->sr->m),m, s2i(a->sr->n),n);
	
	double opt; int lwork=-1, info;
	dgeqrf_(&m, &n, NULL, &m, NULL, &opt, &lwork, &info);
	km_check_info_opt(info, "dgeqrf");
	lwork = (int)opt;
	KALLOCc(a->a, a->sa);
	KALLOC(a->tau, k);
	KALLOC(a->work, lwork);
	dgeqrf_(&m, &n, a->a, &m, a->tau, a->work, &lwork, &info);
	km_check_info(info, rb_eRuntimeError, "unexpected info value", "dgeqrf");
	ruby_xfree(a->work);
	
	lwork = -1;
	dorgqr_(&m, &m, &k, NULL, &m, NULL, &opt, &lwork, &info); // the second argument N is &m
	km_check_info_opt(info, "dorgqr");
	lwork = (int)opt;
	KALLOC(a->work, lwork);
	KALLOCn(a->q, a->sq);
	KALLOCn(a->r, a->sr);
	int ldr=s2i(a->r.ld), ldq=s2i(a->q.ld);
	for ( int i=m; i<k; i++ ) { // 0 clear (i>=k)-th rows of r and (i>=k)-th columns of q
		const int izero=0; const double dzero=0.0;
		dcopy_(&n, &dzero, &izero, (a->r.d)+i, &ldr);
		const int ione=1;
		dcopy_(&m, &dzero, &izero, (a->q.d)+i*ldq, &(ione));
	}
	for ( int i=0; i<k; i++ ) { // copy the results of `dgeqrf' to q.d, rd, construct R and prepare to call `dorgqr'
		const int izero=0; const double dzero=0.0;
		dcopy_(&i, &dzero, &izero, (a->r.d)+i, &ldr);
		int l = n-i;
		dcopy_(&l, (a->a)+(i+i*m), &m, (a->r.d)+(i+i*ldr), &ldr);
		const int ione = 1;
		dcopy_(&i, &dzero, &izero, (a->q.d)+(i*ldq), &ione);
		(a->q.d)[i+i*ldq] = 1.0;
		l = m-i-1;
		dcopy_(&l, (a->a)+(i+i*m+1), &ione, (a->q.d)+(i+i*ldq+1), &ione);
	}
	dorgqr_(&m, &m, &k, a->q.d, &ldq, a->tau, a->work, &lwork, &info);
	km_check_info(info, rb_eRuntimeError, "unexpected info value", "dorgqr");
	
	return Qnil;
}
static VALUE
km_qr_ensure(VALUE data)
{
	struct km_qr_arg *a = (struct km_qr_arg *)data;
	
	km_copy_and_free_if_needed(a->sr, &(a->r));
	km_copy_and_free_if_needed(a->sq, &(a->q));
	ruby_xfree(a->work);
	ruby_xfree(a->tau);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat_qr_destl(VALUE self, VALUE vq, VALUE vr)
{
	km_check_frozen(vq); km_check_frozen(vr);
	struct km_qr_arg a; memset(&a, 0, sizeof(a));
	a.sa = km_mat2smat(self); a.sq = km_mat2smat(vq), a.sr = km_mat2smat(vr);
	km_check_double(3, a.sa, a.sq, a.sr);
	km_check_finite(a.sa);
	
	km_ensure(km_qr_body, (VALUE)&a, km_qr_ensure, (VALUE)&a);
	
	return Qnil;
}
VALUE
kmm_mat_qr(VALUE self)
{
	SMAT *sa = km_mat2smat(self);
	VALUE vq = km_Mat(sa->m, sa->m, VT_DOUBLE);
	VALUE vr = km_Mat(sa->m, sa->n, VT_DOUBLE);
	kmm_mat_qr_destl(self, vq, vr);
	return rb_ary_new3(2, vq, vr);
}

// construct random orthogonal matrix via QR decomposition
struct km_rand_orth_arg {
	SMAT *smat;
	VALUE random;
	double *tau, *work;
	int lwork;
	LAWORK a;
};
static VALUE
km_rand_orth_body(VALUE data)
{
	struct km_rand_orth_arg *a = (struct km_rand_orth_arg *)data;
	
	int n = s2i(a->smat->m);
	km_check_size(1, s2i(a->smat->n),n);
	
	double opt; int lwork=-1, info;
	dgeqrf_(&n, &n, NULL, &n, NULL, &opt, &lwork, &info);
	km_check_info_opt(info, "dgeqrf");
	lwork = (int)opt;
	
	KALLOC(a->tau, n);
	KALLOC(a->work, lwork);
	KALLOCn(a->a, a->smat);
	
	km_fill_normal(i2s(n*n), a->a.d, a->random);
	dgeqrf_(&n, &n, a->a.d, &n, a->tau, a->work, &lwork, &info);
	km_check_info(info, rb_eRuntimeError, "unexpected info value", "dgeqrf");
	ruby_xfree(a->work);
	
	lwork = -1;
	dorgqr_(&n, &n, &n, NULL, &n, NULL, &opt, &lwork, &info);
	km_check_info_opt(info, "dorgqr");
	lwork = (int)opt;
	
	KALLOC(a->work, lwork);
	int lda = s2i(a->a.ld);
	
	for (int i=0; i<n; i++ ) { // clear R to prepare to call `dorgqr'
		const int l = n-i, izero = 0; const double dzero = 0.0; const int x=i+i*lda;
		dcopy_(&l, &dzero, &izero, (a->a.d)+x, &lda);
		(a->a.d)[x] = 1.0;
	}
	dorgqr_(&n, &n, &n, a->a.d, &lda, a->tau, a->work, &lwork, &info);
	km_check_info(info, rb_eRuntimeError, "unexpected info value", "dorgqr");
	// multiply by scalar -1 with probability 1/2
	if ( rb_funcall(a->random, id_rand, 1, INT2NUM(2)) == INT2NUM(0) ) {
		const int l = n*n, ione=1; const double m1=-1.0;
		dscal_(&l, &m1, a->a.d, &ione);
	}
	
	return Qnil;
}
static VALUE
km_rand_orth_ensure(VALUE data)
{
	struct km_rand_orth_arg *a = (struct km_rand_orth_arg *)data;
	
	km_copy_and_free_if_needed(a->smat, &(a->a));
	ruby_xfree(a->work);
	ruby_xfree(a->tau);
	
	return Qnil;
}
VALUE
kmm_mat__rand_orth(VALUE self, VALUE random)
{
	km_check_frozen(self);
	struct km_rand_orth_arg a; memset(&a, 0, sizeof(a));
	a.smat = km_mat2smat(self);
	km_check_double(1, a.smat);
	
	a.random = random;
	km_ensure(km_rand_orth_body, (VALUE)&a, km_rand_orth_ensure, (VALUE)&a);
	
	return self;
}

// ノルムのバランスを取った AA = DD \ A * DD となる AA および DD を計算し，それぞれ引数に格納する
// compute norm balanced AA and factor DD which satisfy AA = DD \ A * DD
// `self' is the input A and arguments are the outputs
struct km_balance_arg {
	SMAT *sa, *sd, *saa;
	double *scale;
	int *perm;
	LAWORK a, d;
};
static VALUE
km_balance_body(VALUE data)
{
	struct km_balance_arg *a = (struct km_balance_arg *)data;
	
	int n = s2i(a->sa->m);
	km_check_size(5, s2i(a->sa->n),n, s2i(a->sd->m),n, s2i(a->sd->n),n, s2i(a->saa->m),n, s2i(a->saa->n),n);
	
	KALLOC(a->scale, n);
	KALLOC(a->perm, n);
	KALLOCn(a->a, a->saa);
	KALLOCz(a->d, a->sd);
	for ( int i=0; i<n; i++ ) {
		(a->perm)[i] = i;
	}
	
	int ilo, ihi, info;
	char job[] = "B";
	int lda=s2i(a->a.ld);
	dgebal_(job, &n, a->a.d, &lda, &ilo, &ihi, a->scale, &info);
	km_check_info(info, rb_eRuntimeError, "unexpected info value", "dgebal");
	for ( int i=ilo-2; 0 <= i; i-- ) {
		SWAP(int, (a->perm)[i], (a->perm)[(int)((a->scale)[i])-1]);
	}
	for ( int i=n-1; ihi <= i; i-- ) {
		SWAP(int, (a->perm)[i], (a->perm)[(int)((a->scale)[i])-1]);
	}
	for ( int i=0; i<n; i++ ) {
		if ( i < ilo-1 || ihi <= i ) { // permutation
			(a->d.d)[(a->perm)[i]+i*s2i(a->d.ld)] = 1.0;
		} else {
			(a->d.d)[(a->perm)[i]+i*s2i(a->d.ld)] = (a->scale)[i];
		}
	}
	
	return Qnil;
}
static VALUE
km_balance_ensure(VALUE data)
{
	struct km_balance_arg *a = (struct km_balance_arg *)data;
	
	km_copy_and_free_if_needed(a->sd, &(a->d));
	km_copy_and_free_if_needed(a->saa, &(a->a));
	ruby_xfree(a->perm);
	ruby_xfree(a->scale);
	
	return Qnil;
}
VALUE
kmm_mat_balance_destl(VALUE self, VALUE vd, VALUE vaa)
{
	km_check_frozen(vd); km_check_frozen(vaa);
	struct km_balance_arg a; memset(&a, 0, sizeof(a));
	a.sa = km_mat2smat(self); a.sd = km_mat2smat(vd); a.saa = km_mat2smat(vaa);
	km_check_double(3, a.sa, a.sd, a.saa);
	km_check_finite(a.sa);
	
	km_ensure(km_balance_body, (VALUE)&a, km_balance_ensure, (VALUE)&a);
	
	return Qnil;
}
VALUE
kmm_mat_balance(VALUE self)
{
	const size_t n = km_mat2smat(self)->m;
	VALUE vd = km_Mat(n, n, VT_DOUBLE);
	VALUE vaa = km_Mat(n, n, VT_DOUBLE);
	kmm_mat_balance_destl(self, vd, vaa);
	return rb_ary_new3(2, vd, vaa);
}

// compute matrix exponential of `self' and replace `self' by the result
// this method is a manual translation of GNU Octave's expm.m
struct km_expm_arg {
	SMAT *sa;
	LAWORK a;
	double *scale;
	double *a2, *x, *y, *foo;
	double *r, *c, *ferr, *berr, *work;
	int *ipiv, *iwork;
	VALUE self;
};
static VALUE
km_expm_body(VALUE data)
{
	struct km_expm_arg *a = (struct km_expm_arg *)data;
	
	int n = s2i(a->sa->m); const int n2 = n*n;
	km_check_size(1, s2i(a->sa->n),n);
	
	// trace reduction
	double neg_max = -DBL_MAX;
	KALLOCn(a->a, a->sa);
	int lda=s2i(a->a.ld);
	for ( int i=0; i<n; i++ ) { for ( int j=0; j<n; j++ ) {
		if ( (a->a.d)[i+j*lda] < neg_max ) { (a->a.d)[i+j*lda] = neg_max; }
	} }
	double trshift = 0.0;
	for ( int i=0; i<n; i++ ) {
		trshift += (a->a.d)[i+i*lda];
	}
	if ( 0 < trshift ) {
		trshift /= n;
		for ( int i=0; i<n; i++ ) {
			(a->a.d)[i+i*lda] -= trshift;
		}
	}
	
	// balancing
	int ilo, ihi, info;
	char job[] = "B";
	KALLOC(a->scale, n);
	dgebal_(job, &n, a->a.d, &lda, &ilo, &ihi, a->scale, &info);
	km_check_info(info, rb_eRuntimeError, "unexpected info value", "dgebal");
	
	// scaling
	const int ione=1;
	double s = 0.0;
	for ( int i=0; i<n; i++ ) {
		double foo = dasum_(&n, (a->a.d)+i, &lda);
		if ( s < foo ) { s = foo; }
	}
	s = logb(s);
	if ( s < 0.0 ) {
		s = 0.0;
	} else if ( 1023.0 < s ) {
		s = 1023.0;
	}
	const double ps = exp2(-s);
	for ( int i=0; i<n; i++ ) {
		dscal_(&n, &ps, (a->a.d)+(i*lda), &ione);
	}
	
	// Pade approximation
	static const double c[] = { 5.0000000000000000e-1, 1.1666666666666667e-1, 1.6666666666666667e-2, 1.6025641025641026e-3,
		1.0683760683760684e-4, 4.8562548562548563e-6, 1.3875013875013875e-7, 1.9270852604185938e-9 };
	const int np1 = n+1;
	KALLOC(a->a2, n2);
	KALLOC(a->x, n2);
	KALLOC(a->y, n2);
	KALLOC(a->foo, n2);
#define MPROD(_r, _a, _b) km_dmprod(n, n, n, _a, _b, _r)
#define A a->a.d, lda
#define A2 a->a2, n
#define X a->x, n
#define Y a->y, n
#define FOO a->foo, n
	MPROD(A2, A, A); // a2 = a*a
	memcpy( a->foo, a->a2, sizeof(double)*i2s(n2) ); // foo = a2
	dscal_( &n2, c+7, a->foo, &ione ); // foo *= c[7]
	for ( int i=0; i<n; i++ ) { a->foo[i*np1] += c[5]; } // foo += c[5]*I
	MPROD(X, FOO, A2); // x = foo*a2 = a^4*c[7]+a^2*c[5]
	for ( int i=0; i<n; i++ ) { a->x[i*np1] += c[3]; } // x += c[3]*I
	MPROD(FOO, X, A2); // foo = x*a2 = a^6*c[7]+a^4*c[5]+a^2*c[3]
	for ( int i=0; i<n; i++ ) { a->foo[i*np1] += c[1]; } // foo += c[1]*I
	MPROD(X, FOO, A2); // x = foo*a2 = a^8*c[7]+a^6*c[5]+a^4*c[3]+a^2*c[1]
	for ( int i=0; i<n; i++ ) { a->x[i*np1] += 1.0; } // x += I
	memcpy( a->foo, a->a2, sizeof(double)*((size_t)n2) ); // foo = a2
	dscal_(&n2, c+6, a->foo, &ione); // foo *= c[6]
	for ( int i=0; i<n; i++ ) { a->foo[i*np1] += c[4]; } // foo += c[4]*I
	MPROD(Y, FOO, A2); // y = foo*a2 = a^4*c[6]+a^2*c[4]
	for ( int i=0; i<n; i++ ) { a->y[i*np1] += c[2]; } // y += c[2]*I
	MPROD(FOO, Y, A2); // foo = y*a2 = a^6*c[6]+a^4*c[4]+a^2*c[2]
	for ( int i=0; i<n; i++ ) { a->foo[i*np1] += c[0]; } // foo += c[0]*I
	MPROD(Y, FOO, A); // y = foo*a = a^7*c[6]+a^5*c[4]+a^3*c[2]+a*c[0]
	double alp=1.0;
	daxpy_(&n2, &alp, a->y, &ione, a->x, &ione); // x = X+Y
	alp = -2.0;
	dscal_(&n2, &alp, a->y, &ione);
	alp = 1.0;
	daxpy_(&n2, &alp, a->x, &ione, a->y, &ione); // y = -2y+x i.e. -2Y+(X+Y) = X-Y
	
	// r = y\x i.e. (X-Y)\(X+Y)
	KALLOC(a->ipiv, n);
	KALLOC(a->r, n);
	KALLOC(a->c, n);
	KALLOC(a->ferr, n);
	KALLOC(a->berr, n);
	KALLOC(a->work, 4*n);
	KALLOC(a->iwork, n);
	char equed[2];
	char fact[]="E";
	char trans[]="N";
	dgesvx_(fact, trans, &n, &n, a->y, &n, a->a2, &n, a->ipiv, equed, a->r, a->c, a->x, &n,
		a->a.d, &lda, &alp, a->ferr, a->berr, a->work, a->iwork, &info);
	km_check_info(info, km_eUncomp, "an internal matrix is singular or near singular", "dgesvx");
	
	// undo scaling by repeated squaring
	const int is = (int)s;
	if ( is & 1 ) { // if is is odd, then r = r^2
		MPROD(FOO, A, A);
		dlacpy_(job, &n, &n, a->foo, &n, a->a.d, &lda);
	}
	for ( int i=0; i<is/2; i++ ) {
		MPROD(FOO, A, A);
		MPROD(A, FOO, FOO);
	}
	
	// inverse balancing
	km_copy_and_free_if_needed(a->sa, &(a->a));
	km_unbal(a->sa, a->scale, ilo, ihi);
	
	// inverse trace reduction
	if ( 0 < trshift ) {
		kmm_mat_s_mul_destl(a->self, rb_float_new(exp(trshift)));
	}
	
	return Qnil;
}
static VALUE
km_expm_ensure(VALUE data)
{
	struct km_expm_arg *a = (struct km_expm_arg *)data;
	
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	ruby_xfree(a->berr);
	ruby_xfree(a->ferr);
	ruby_xfree(a->c);
	ruby_xfree(a->r);
	ruby_xfree(a->ipiv);
	
	ruby_xfree(a->foo);
	ruby_xfree(a->y);
	ruby_xfree(a->x);
	ruby_xfree(a->a2);
	km_free_if_needed(&(a->a));
	ruby_xfree(a->scale);
	
	return Qnil;
}
VALUE
kmm_mat_expm_dest(VALUE self)
{
	km_check_frozen(self);
	struct km_expm_arg a; memset(&a, 0, sizeof(a));
	a.sa = km_mat2smat(self);
	km_check_double(1, a.sa);
	km_check_finite(a.sa);
	
	a.self = self;
	km_ensure(km_expm_body, (VALUE)&a, km_expm_ensure, (VALUE)&a);
	
	return self;
}

// compute a Cholesky decomposition A = U' * U
// on entry, `self' is A. on exit, `self' is U
struct km_chol_arg {
	SMAT *sa;
	LAWORK a;
};
static VALUE
km_chol_body(VALUE data)
{
	struct km_chol_arg *a = (struct km_chol_arg *)data;
	
	int n = s2i(a->sa->m), info;
	
	KALLOCn(a->a, a->sa);
	char uplo[] = "U";
	int lda=s2i(a->a.ld);
	dpotrf_(uplo, &n, a->a.d, &lda, &info);
	km_check_info(info, km_eUncomp, "self is not positive definite", "dpotrf");
	int ione=1, izero=0; double dzero=0.0;
	for ( int i=0; i<n-1; i++ ) {
		const int len=n-i-1;
		dcopy_(&len, &dzero, &izero, (a->a.d)+(i+i*lda+1), &ione);
	}
	
	return Qnil;
}
static VALUE
km_chol_ensure(VALUE data)
{
	struct km_chol_arg *a = (struct km_chol_arg *)data;
	
	km_copy_and_free_if_needed(a->sa, &(a->a));
	
	return Qnil;
}
VALUE
kmm_mat_chol_dest(VALUE self)
{
	km_check_frozen(self);
	if ( !kmm_mat_symmetry_p(0, NULL, self) ) {
		rb_raise(km_eUncomp, "self must be a symmetry matrix");
	}
	struct km_chol_arg a; memset(&a, 0, sizeof(a));
	a.sa = km_mat2smat(self);
	km_check_double(1, a.sa);
	km_check_finite(a.sa);
	
	km_ensure(km_chol_body, (VALUE)&a, km_chol_ensure, (VALUE)&a);
	
	return self;
}
