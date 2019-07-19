#include "../kmat.h"

#define SET_LD(ldx, smat) int ldx; \
	if ( smat->n == 1 ) { \
		if ( smat->trans ) { \
			ldx = smat->ld; \
		} else { \
			ldx = 1; \
		} \
	} else { \
		if ( smat->trans ) { \
			ldx = 1; \
		} else { \
			ldx = smat->ld; \
		} \
	}
#define DEF_IPROD(id, type, mul, add, zero) static type \
km_##id##mat_iprod(SMAT *sa, SMAT *sb) { \
	type ret; \
	SET_LD(lda, sa) \
	SET_LD(ldb, sb) \
	int len = LENGTH(sa); \
	if ( len == 0 ) { return zero; } \
	if ( sa->stype == ST_RSUB ) { \
		if ( sb->stype == ST_RSUB ) { \
			ret = mul( *(sa->id##pbody[0]), *(sb->id##pbody[0]) ); \
			for ( int i=1; i<len; i++ ) { \
				type temp = mul( *(sa->id##pbody[i*lda]), *(sb->id##pbody[i*ldb]) ); \
				ret = add(ret, temp); \
			} \
		} else { \
			ret = mul( *(sa->id##pbody[0]), (sb->id##body[0]) ); \
			for ( int i=1; i<len; i++ ) { \
				type temp = mul( *(sa->id##pbody[i*lda]), (sb->id##body[i*ldb]) ); \
				ret = add(ret, temp); \
			} \
		} \
	} else { \
		if ( sb->stype == ST_RSUB ) { \
			ret = mul( (sa->id##body[0]), *(sb->id##pbody[0]) ); \
			for ( int i=1; i<len; i++ ) { \
				type temp = mul( (sa->id##body[i*lda]), *(sb->id##pbody[i*ldb]) ); \
				ret = add(ret, temp); \
			} \
		} else { \
			ret = mul( (sa->id##body[0]), (sb->id##body[0]) ); \
			for ( int i=1; i<len; i++ ) { \
				type temp = mul( (sa->id##body[i*lda]), (sb->id##body[i*ldb]) ); \
				ret = add(ret, temp); \
			} \
		} \
	} \
	return ret; \
}
#define NMUL(a, b) ((a)*(b))
#define NADD(a, b) ((a)+(b))
#define BMUL(a, b) ((a)&&(b))
#define VMUL(a, b) rb_funcall((a), id_op_mul, 1, (b))
#define VADD(a, b) rb_funcall((a), id_op_plus, 1, (b))
DEF_IPROD(d, double, NMUL, NADD, 0.0)
DEF_IPROD(z, COMPLEX, NMUL, NADD, cpack(0.0, 0.0))
DEF_IPROD(i, int, NMUL, NADD, 0)
DEF_IPROD(b, bool, BMUL, XOR, false)
DEF_IPROD(v, VALUE, VMUL, VADD, INT2NUM(0))
static double
km_dmat_iprodb(SMAT *sa, SMAT *sb)
{
	if ( sa->stype == ST_RSUB || sb->stype == ST_RSUB ) {
		return km_dmat_iprod(sa, sb);
	}
	SET_LD(lda, sa)
	SET_LD(ldb, sb)
	int n = LENGTH(sa);
	return ddot_(&n, sa->dbody, &lda, sb->dbody, &ldb);
}
static COMPLEX
km_zmat_iprodb(SMAT *sa, SMAT *sb)
{
	if ( sa->stype == ST_RSUB || sb->stype == ST_RSUB ) {
		return km_zmat_iprod(sa, sb);
	}
	SET_LD(lda, sa)
	SET_LD(ldb, sb)
	int n = LENGTH(sa);
	return zdotu_(&n, sa->zbody, &lda, sb->zbody, &ldb);
}
VALUE
kmm_mat__iprod(VALUE self, VALUE vb)
{
	SMAT *sa = km_mat2smat(self), *sb = km_mat2smat(vb);
	if ( !VECTOR_P(sa) || !VECTOR_P(sb) ) {
		rb_raise(km_eDim, "both self and the argument must be vectors");
	} else if ( LENGTH(sa) != LENGTH(sb) ) {
		rb_raise(km_eDim, "dimensions must be the same (%d != %d)", LENGTH(sa), LENGTH(sb));
	} else if ( sa->vtype != sb->vtype ) {
		rb_raise(km_eVT, "value types must be the same");
	}
	VT_SWITCH( sa->vtype,
		return rb_float_new(km_dmat_iprodb(sa, sb));,
		return km_c2v(km_zmat_iprodb(sa, sb));,
		return INT2NUM(km_imat_iprod(sa, sb));,
		return TF2V(km_bmat_iprod(sa, sb));,
		return km_vmat_iprod(sa, sb);
	);
}

// Frobenius norm
static void
km_normf_func_d(double *ent, void *data)
{
	double *ret = (double *)data;
	*ret = hypot(*ret, *ent);
}
static void
km_normf_func_z(COMPLEX *ent, void *data)
{
	double *ret = (double *)data;
	*ret = hypot(*ret, cabs(*ent));
}
static void
km_normf_func_i(int *ent, void *data)
{
	double *ret = (double *)data;
	*ret = hypot(*ret, (double)(*ent));
}
static void
km_normf_func_b(bool *ent, void *data)
{
	bool *ret = (bool *)data;
	*ret = XOR(*ret, *ent);
}
static void
km_normf_func_v(VALUE *ent, void *data)
{
	VALUE *ret = (VALUE *)data;
	*ret = rb_funcall(rb_mMath, id_hypot, 2, *ret, *ent);
}
VALUE
kmm_mat_normf(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype == VT_DOUBLE ) {
		if ( smat->stype == ST_FULL ) {
			int ione=1, n = LENGTH(smat);
			return rb_float_new(dnrm2_(&n, smat->dbody, &ione));
		} else if ( smat->stype == ST_SSUB ) {
			int ione=1, lps, size; double ret=0.0;
			if ( smat->trans ) {
				lps = smat->m; size = smat->n;
			} else {
				lps = smat->n; size = smat->m;
			}
			for ( int i=0; i<lps; i++ ) {
				ret = hypot(ret, dnrm2_(&size, smat->dbody+i*(smat->ld), &ione));
			}
			return rb_float_new(ret);
		} else {
			double ret = 0.0;
			km_smat_each_d(smat, km_normf_func_d, (void *)(&ret));
			return rb_float_new(ret);
		}
	} else if ( smat->vtype == VT_COMPLEX ) {
		if ( smat->stype == ST_FULL ) {
			int ione=1, n = LENGTH(smat);
			return rb_float_new(dznrm2_(&n, smat->zbody, &ione));
		} else if ( smat->stype == ST_SSUB ) {
			int ione=1, lps, size; double ret=0.0;
			if ( smat->trans ) {
				lps = smat->m; size = smat->n;
			} else {
				lps = smat->n; size = smat->m;
			}
			for ( int i=0; i<lps; i++ ) {
				ret = hypot(ret, dznrm2_(&size, smat->zbody+i*(smat->ld), &ione));
			}
			return rb_float_new(ret);
		} else {
			double ret = 0.0;
			km_smat_each_z(smat, km_normf_func_z, (void *)(&ret));
			return rb_float_new(ret);
		}
	} else if ( smat->vtype == VT_INT ) {
		double ret = 0.0;
		km_smat_each_i(smat, km_normf_func_i, (void *)(&ret));
		return rb_float_new(ret);
	} else if ( smat->vtype == VT_BOOL ) {
		bool ret = false;
		km_smat_each_b(smat, km_normf_func_b, (void *)(&ret));
		return TF2V(ret);
	} else if ( smat->vtype == VT_VALUE ) {
		VALUE ret = INT2NUM(0);
		km_smat_each_v(smat, km_normf_func_v, (void *)(&ret));
		return ret;
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
}

// L-1, L-infinity norm
#define NORM1(type, id, m_ent, m, n, m_abs, m_add, m_cmp, x2v) \
	type ret = m_abs(m_ent(smat, id, 0)); \
	for ( int i=1; i<m; i++ ) { \
		type foo = m_abs(m_ent(smat, id, i)); \
		ret = m_add(ret, foo); \
	} \
	for ( int j=1; j<n; j++ ) { \
		type bar = m_abs(m_ent(smat, id, j*(smat->ld))); \
		for ( int i=1; i<m; i++ ) { \
			type foo = m_abs(m_ent(smat, id, i+j*(smat->ld))); \
			bar = m_add(bar, foo); \
		} \
		if ( m_cmp(ret, bar) ) { ret = bar; } \
	} \
	return x2v(ret)
#define NORM1_DZ(id, bid, m,  n, m_abs) do { \
	if ( smat->stype == ST_RSUB ) { \
		NORM1(double, id, ENTITYr0, m, n, m_abs, NUM_ADD, NUM_CMP, rb_float_new); \
	} else { \
		const int ione=1; \
		double ret = bid##asum_(&(m), smat->id##body, &ione); \
		for ( int j=1; j<n; j++ ) { \
			double temp = bid##asum_(&(m), smat->id##body+(j*(smat->ld)), &ione); \
			if ( ret < temp ) { ret = temp; } \
		} \
		return rb_float_new(ret); \
	} \
} while (0)
#define NORMi(type, id, m_ent, m, n, m_abs, m_add, m_cmp, x2v) \
	type ret = m_abs(m_ent(smat, id, 0)); \
	for ( int j=1; j<n; j++ ) { \
		type foo = m_abs(m_ent(smat, id, j*(smat->ld))); \
		ret = m_add(ret, foo); \
	} \
	for ( int i=1; i<m; i++ ) { \
		type bar = m_abs(m_ent(smat, id, i)); \
		for ( int j=1; j<n; j++ ) { \
			type foo = m_abs(m_ent(smat, id, i+j*(smat->ld))); \
			bar = m_add(bar, foo); \
		} \
		if ( m_cmp(ret, bar) ) { ret = bar; } \
	} \
	return x2v(ret)
#define NORMi_DZ(id, bid, m,  n, m_abs) do { \
	if ( smat->stype == ST_RSUB ) { \
		NORMi(double, id, ENTITYr0, m, n, m_abs, NUM_ADD, NUM_CMP, rb_float_new); \
	} else { \
		double ret = bid##asum_(&(n), smat->id##body, &(smat->ld)); \
		for ( int i=1; i<m; i++ ) { \
			double temp = bid##asum_(&(n), smat->id##body+i, &(smat->ld)); \
			if ( ret < temp ) { ret = temp; } \
		} \
		return rb_float_new(ret); \
	} \
} while (0)
#define NORM1i_OTHER(oi, type, id, m, n, m_abs, m_add, m_cmp, x2v) do { \
	if ( smat->stype == ST_RSUB ) { \
		NORM##oi(type, id, ENTITYr0, m, n, m_abs, m_add, m_cmp, x2v); \
	} else { \
		NORM##oi(type, id, ENTITYd0, m, n, m_abs, m_add, m_cmp, x2v); \
	} \
} while (0)
#define NUM_ADD(a, b) a+b
#define NUM_CMP(a, b) a<b
#define VAL_ABS(a) rb_funcall(a, id_abs, 0)
#define VAL_ADD(a, b) rb_funcall(a, id_op_plus, 1, b)
#define VAL_CMP(a, b) RTEST(rb_funcall(a, id_op_lt, 1, b))

// induced L-1 norm (maximum of absolute summation of each column)
VALUE
kmm_mat_norm1(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->trans ) {
		VT_SWITCH( smat->vtype,
			NORMi_DZ(d, d, smat->n, smat->m, fabs);,
			NORMi_DZ(z, dz, smat->n, smat->m, cabs);,
			NORM1i_OTHER(i, int, i, smat->n, smat->m, ABS, NUM_ADD, NUM_CMP, INT2NUM);,
			NORM1i_OTHER(i, bool, b, smat->n, smat->m, ITSELF, XOR, NUM_CMP, TF2V);,
			NORM1i_OTHER(i, VALUE, v, smat->n, smat->m, VAL_ABS, VAL_ADD, VAL_CMP, ITSELF);
		);
	} else {
		VT_SWITCH( smat->vtype,
			NORM1_DZ(d, d, smat->m, smat->n, fabs);,
			NORM1_DZ(z, dz, smat->m, smat->n, cabs);,
			NORM1i_OTHER(1, int, i, smat->m, smat->n, ABS, NUM_ADD, NUM_CMP, INT2NUM);,
			NORM1i_OTHER(1, bool, b, smat->m, smat->n, ITSELF, XOR, NUM_CMP, TF2V);,
			NORM1i_OTHER(1, VALUE, v, smat->m, smat->n, VAL_ABS, VAL_ADD, VAL_CMP, ITSELF);
		);
	}
}
// induced L-infinity norm (maximum of absolute summation of each row)
VALUE
kmm_mat_normi(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->trans ) {
		VT_SWITCH( smat->vtype,
			NORM1_DZ(d, d, smat->n, smat->m, fabs);,
			NORM1_DZ(z, dz, smat->n, smat->m, cabs);,
			NORM1i_OTHER(1, int, i, smat->n, smat->m, ABS, NUM_ADD, NUM_CMP, INT2NUM);,
			NORM1i_OTHER(1, bool, b, smat->n, smat->m, ITSELF, XOR, NUM_CMP, TF2V);,
			NORM1i_OTHER(1, VALUE, v, smat->n, smat->m, VAL_ABS, VAL_ADD, VAL_CMP, ITSELF);
		);
	} else {
		VT_SWITCH( smat->vtype,
			NORMi_DZ(d, d, smat->m, smat->n, fabs);,
			NORMi_DZ(z, dz, smat->m, smat->n, cabs);,
			NORM1i_OTHER(i, int, i, smat->m, smat->n, ABS, NUM_ADD, NUM_CMP, INT2NUM);,
			NORM1i_OTHER(i, bool, b, smat->m, smat->n, ITSELF, XOR, NUM_CMP, TF2V);,
			NORM1i_OTHER(i, VALUE, v, smat->m, smat->n, VAL_ABS, VAL_ADD, VAL_CMP, ITSELF);
		);
	}
}

// element-wise L-1 norm (absolute summation of all elements)
static void
km_asum_d(double *ent, void *data)
{
	double *ret = (double *)data;
	*ret += fabs(*ent);
}
static void
km_asum_z(COMPLEX *ent, void *data)
{
	double *ret = (double *)data;
	*ret += cabs(*ent);
}
static void
km_asum_i(int *ent, void *data)
{
	int *ret = (int *)data;
	*ret += ABS(*ent);
}
static void
km_asum_b(bool *ent, void *data)
{
	bool *ret = (bool *)data;
	*ret = XOR(*ret, *ent);
}
static void
km_asum_v(VALUE *ent, void *data)
{
	VALUE *ret = (VALUE *)data;
	VALUE temp = rb_funcall(*ent, id_abs, 0);
	*ret = rb_funcall(*ret, id_op_plus, 1, temp);
}
VALUE
kmm_mat_norm_e1(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype == VT_DOUBLE ) {
		if ( smat->stype == ST_FULL ) {
			int ione = 1, n = LENGTH(smat);
			return rb_float_new(dasum_(&n, smat->dbody, &ione));
		} else if ( smat->stype == ST_SSUB ) {
			int ione = 1, lps, size; double ret = 0.0;
			if ( smat->trans ) {
				lps = smat->m; size = smat->n;
			} else {
				lps = smat->n; size = smat->m;
			}
			for ( int i=0; i<lps; i++ ) {
				ret += dasum_(&size, smat->dbody+(i*(smat->ld)), &ione);
			}
			return rb_float_new(ret);
		} else {
			double ret = 0.0;
			km_smat_each_d(smat, km_asum_d, (void *)&ret);
			return rb_float_new(ret);
		}
	} else if ( smat->vtype == VT_COMPLEX ) {
		if ( smat->stype == ST_FULL ) {
			int ione = 1, n = LENGTH(smat);
			return rb_float_new(dzasum_(&n, smat->zbody, &ione));
		} else if ( smat->stype == ST_SSUB ) {
			int ione = 1, lps, size; double ret = 0.0;
			if ( smat->trans ) {
				lps = smat->m; size = smat->n;
			} else {
				lps = smat->n; size = smat->m;
			}
			for ( int i=0; i<lps; i++ ) {
				ret += dzasum_(&size, smat->zbody+(i*(smat->ld)), &ione);
			}
			return rb_float_new(ret);
		} else {
			double ret = 0.0;
			km_smat_each_z(smat, km_asum_z, (void *)&ret);
			return rb_float_new(ret);
		}
	} else if ( smat->vtype == VT_INT ) {
		int ret = 0;
		km_smat_each_i(smat, km_asum_i, (void *)&ret);
		return INT2NUM(ret);
	} else if ( smat->vtype == VT_BOOL ) {
		bool ret = false;
		km_smat_each_b(smat, km_asum_b, (void *)&ret);
		return TF2V(ret);
	} else if ( smat->vtype == VT_VALUE ) {
		VALUE ret = INT2NUM(0);
		km_smat_each_v(smat, km_asum_v, (void *)&ret);
		return ret;
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
}

// element-wise L-infinity norm (maximum of absolute value of each element)
static void
km_amax_d(double *ent, void *data)
{
	double *ret = (double *)data;
	double temp = fabs(*ent);
	if ( *ret < temp ) {
		*ret = temp;
	}
}
static void
km_amax_z(COMPLEX *ent, void *data)
{
	double *ret = (double *)data;
	double temp = cabs(*ent);
	if ( *ret < temp ) {
		*ret = temp;
	}
}
static void
km_amax_i(int *ent, void *data)
{
	int *ret = (int *)data;
	int temp = ABS(*ent);
	if ( *ret < temp ) {
		*ret = temp;
	}
}
static void
km_amax_b(bool *ent, void *data)
{
	bool *ret = (bool *)data;
	*ret = ( *ret || *ent );
}
static void
km_amax_v(VALUE *ent, void *data)
{
	VALUE *ret = (VALUE *)data;
	VALUE temp = rb_funcall(*ent, id_abs, 0);
	if ( rb_funcall(*ret, id_op_lt, 1, temp) ) {
		*ret = temp;
	}
}
VALUE
kmm_mat_norm_einf(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype == VT_DOUBLE ) {
		if ( smat->stype == ST_FULL ) {
			int ione=1, n=LENGTH(smat);
			return rb_float_new(fabs(smat->dbody[idamax_(&n, smat->dbody, &ione)-1]));
		} else if ( smat->stype == ST_SSUB ) {
			int ione=1, lps, size; double ret = 0.0;
			if ( smat->trans ) {
				lps = smat->m; size = smat->n;
			} else {
				lps = smat->n; size = smat->m;
			}
			for ( int i=0; i<lps; i++ ) {
				double temp = fabs(smat->dbody[idamax_(&size, smat->dbody+(i*(smat->ld)), &ione)-1]);
				if ( ret < temp ) { ret = temp; }
			}
			return rb_float_new(ret);
		} else {
			double ret = 0.0;
			km_smat_each_d(smat, km_amax_d, (void *)&ret);
			return rb_float_new(ret);
		}
	} else if ( smat->vtype == VT_COMPLEX ) {
		if ( smat->stype == ST_FULL ) {
			int ione=1, n=LENGTH(smat);
			return rb_float_new(cabs(smat->zbody[izamax_(&n, smat->zbody, &ione)-1]));
		} else if ( smat->stype == ST_SSUB ) {
			int ione=1, lps, size; double ret = 0.0;
			if ( smat->trans ) {
				lps = smat->m; size = smat->n;
			} else {
				lps = smat->n; size = smat->m;
			}
			for ( int i=0; i<lps; i++ ) {
				double temp = cabs(smat->zbody[izamax_(&size, smat->zbody+(i*(smat->ld)), &ione)-1]);
				if ( ret < temp ) { ret = temp; }
			}
			return rb_float_new(ret);
		} else {
			double ret = 0.0;
			km_smat_each_z(smat, km_amax_z, (void *)&ret);
			return rb_float_new(ret);
		}
	} else if ( smat->vtype == VT_INT ) {
		int ret = 0;
		km_smat_each_i(smat, km_amax_i, (void *)&ret);
		return INT2NUM(ret);
	} else if ( smat->vtype == VT_BOOL ) {
		bool ret = false;
		km_smat_each_b(smat, km_amax_b, (void *)&ret);
		return TF2V(ret);
	} else if ( smat->vtype == VT_VALUE ) {
		VALUE ret = INT2NUM(0);
		km_smat_each_v(smat, km_amax_v, (void *)&ret);
		return ret;
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
}

// induced L-2 norm (maximum of singular values)
struct km_norm2_arg {
	SMAT *sa;
	double *a, *s, *work;
	int *iwork;
};
static VALUE
km_norm2_body(VALUE data)
{
	struct km_norm2_arg *a = (struct km_norm2_arg *)data;
	
	int m = a->sa->m, n = a->sa->n;
	
	double opt; int lwork=-1, info;
	char jobz[] = "N";
	dgesdd_(jobz, &m, &n, NULL, &m, NULL, NULL, &m, NULL, &n, &opt, &lwork, NULL, &info);
	km_check_info(info, rb_eRuntimeError, "error occured while computing optimal lwork", "dgesdd");
	lwork = (int)opt;
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->s, MIN(m, n));
	KALLOC(a->work, lwork);
	KALLOC(a->iwork, 8*MIN(m, n));
	
	dgesdd_(jobz, &m, &n, a->a, &m, a->s, NULL, &m, NULL, &n, a->work, &lwork, a->iwork, &info);
	km_check_info(info, rb_eRuntimeError, "DBDSDC did not converge", "dgesdd");
	
	return rb_float_new((a->s)[0]);
}
static VALUE
km_norm2_ensure(VALUE data)
{
	struct km_norm2_arg *a = (struct km_norm2_arg *)data;
	
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	ruby_xfree(a->s);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat_norm2(VALUE self)
{
	struct km_norm2_arg a;
	a.sa = km_mat2smat(self);
	km_check_double(1, a.sa);
	
	return km_ensure(km_norm2_body, (VALUE)&a, km_norm2_ensure, (VALUE)&a);
}

// estimate for the reciprocal condition of `self'
struct km__rcond2_arg {
	SMAT *sa;
	double *a, *s, *work;
	int *iwork;
};
static VALUE
km__rcond2_body(VALUE data)
{
	struct km__rcond2_arg *a = (struct km__rcond2_arg *)data;
	
	int n = a->sa->m;
	km_check_size(1, a->sa->n,n);
	
	double opt; int lwork=-1, info;
	char jobz[]="N";
	dgesdd_(jobz, &n, &n, NULL, &n, NULL, NULL, &n, NULL, &n, &opt, &lwork, NULL, &info);
	km_check_info(info, rb_eRuntimeError, "error occured while computing optimal lwork", "dgesdd");
	lwork = (int)opt;
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->s, n);
	KALLOC(a->work, lwork);
	KALLOC(a->iwork, 8*n);
	
	dgesdd_(jobz, &n, &n, a->a, &n, a->s, NULL, &n, NULL, &n, a->work, &lwork, a->iwork, &info);
	km_check_info(info, rb_eRuntimeError, "DBDSDC did not converge", "dgesdd");
	
	return rb_float_new((a->s)[n-1]/(a->s)[0]);
}
static VALUE
km__rcond2_ensure(VALUE data)
{
	struct km__rcond2_arg *a = (struct km__rcond2_arg *)data;
	
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	ruby_xfree(a->s);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat__rcond2(VALUE self)
{
	struct km__rcond2_arg a;
	a.sa = km_mat2smat(self);
	km_check_double(1, a.sa);
	
	return km_ensure(km__rcond2_body, (VALUE)&a, km__rcond2_ensure, (VALUE)&a);
}

struct km__rcondf_arg {
	SMAT *sa;
	double *a, *work;
	int *ipiv;
	double normf;
};
static VALUE
km__rcondf_body(VALUE data)
{
	struct km__rcondf_arg *a = (struct km__rcondf_arg *)data;
	
	int n = a->sa->m;
	km_check_size(1, a->sa->n,n);
	
	double opt; int lwork=-1, info;
	dgetri_(&n, NULL, &n, NULL, &opt, &lwork, &info);
	km_check_info(info, rb_eRuntimeError, "error occured while computing optimal lwork", "dgetri");
	lwork = (int)opt;
	KALLOCc(a->a, a->sa);
	KALLOC(a->ipiv, n);
	KALLOC(a->work, lwork);
	
	dgetrf_(&n, &n, a->a, &n, a->ipiv, &info);
	km_check_info(info, km_eUncomp, "the matrix is exactry singular", "dgetrf");
	dgetri_(&n, a->a, &n, a->ipiv, a->work, &lwork, &info);
	km_check_info(info, rb_eRuntimeError, "unexpected info value", "dgetri");
	int len = n*n, ione = 1;
	return rb_float_new(1.0/((a->normf)*dnrm2_(&len, a->a, &ione)));
}
static VALUE
km__rcondf_ensure(VALUE data)
{
	struct km__rcondf_arg *a = (struct km__rcondf_arg *)data;
	
	ruby_xfree(a->work);
	ruby_xfree(a->ipiv);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat__rcondf(VALUE self)
{
	struct km__rcondf_arg a;
	a.sa = km_mat2smat(self);
	km_check_double(1, a.sa);
	
	a.normf = NUM2DBL(kmm_mat_normf(self));
	
	return km_ensure(km__rcondf_body, (VALUE)&a, km__rcondf_ensure, (VALUE)&a);
}

struct km__rcondoi_arg {
	SMAT *sa;
	char oi[2];
	double *a, *work;
	int *iwork;
	double anorm;
};
static VALUE
km__rcondoi_body(VALUE data)
{
	struct km__rcondoi_arg *a = (struct km__rcondoi_arg *)data;
	
	int n = a->sa->m;
	km_check_size(1, a->sa->n,n);
	
	KALLOCc(a->a, a->sa);
	KALLOC(a->iwork, n);
	double rcond;
	int info;
	dgetrf_(&n, &n, a->a, &n, a->iwork, &info);
	km_check_info(info, km_eUncomp, "the matrix is exactry singular", "dgetrf");
	KALLOC(a->work, 4*n);
	dgecon_(a->oi, &n, a->a, &n, &(a->anorm), &rcond, a->work, a->iwork, &info);
	km_check_info(info, rb_eRuntimeError, "unexpected info value", "dgecon");
	
	return rb_float_new(rcond);
}
static VALUE
km__rcondoi_ensure(VALUE data)
{
	struct km__rcondoi_arg *a = (struct km__rcondoi_arg *)data;
	
	ruby_xfree(a->iwork);
	ruby_xfree(a->work);
	ruby_xfree(a->a);
	
	return Qnil;
}
VALUE
kmm_mat__rcondoi(VALUE self, VALUE sym)
{
	struct km__rcondoi_arg a;
	a.sa = km_mat2smat(self);
	km_check_double(1, a.sa);
	if ( sym == sym_one ) {
		strcpy(a.oi, "O");
		a.anorm = NUM2DBL(kmm_mat_norm1(self));
	} else if ( sym == sym_infinity ) {
		strcpy(a.oi, "I");
		a.anorm = NUM2DBL(kmm_mat_normi(self));
	} else {
		rb_raise(rb_eArgError, "unknown norm type");
	}
	
	km_ensure(km__rcondoi_body, (VALUE)&a, km__rcondoi_ensure, (VALUE)&a);
	
	return self;
}

