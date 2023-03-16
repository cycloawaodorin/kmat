#include "../kmat.h"

// product of all elements
void
km_prod_d(double *ent, void *data)
{
	double *ret = (double *)data;
	*ret *= *ent;
}
void
km_prod_z(COMPLEX *ent, void *data)
{
	COMPLEX *ret = (COMPLEX *)data;
	*ret *= *ent;
}
void
km_prod_i(int *ent, void *data)
{
	int *ret = (int *)data;
	*ret *= *ent;
}
void
km_prod_v(VALUE *ent, void *data)
{
	VALUE *ret = (VALUE *)data;
	*ret = rb_funcall(*ret, id_op_mul, 1, *ent);
}
VALUE
kmm_mat_prod(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	VT_SWITCH( smat->vtype,
		double ret=1.0; km_smat_each_d(smat, km_prod_d, (void *)(&ret)); return rb_float_new(ret);,
		COMPLEX ret=cpack(1.0, 0.0); km_smat_each_z(smat, km_prod_z, (void *)(&ret)); return km_c2v(ret);,
		int ret=1; km_smat_each_i(smat, km_prod_i, (void *)(&ret)); return INT2NUM(ret);,
		return rb_funcall(self, id_all_p, 0);,
		VALUE ret=INT2NUM(1); km_smat_each_v(smat, km_prod_v, (void *)(&ret)); return ret;
	);
}

// summation
// the shape of the output-argument determines the axis of summations
void
km_sum_all_d(double *ent, void *data)
{
	double *ret = (double *)data;
	*ret += *ent;
}
void
km_sum_all_z(COMPLEX *ent, void *data)
{
	COMPLEX *ret = (COMPLEX *)data;
	*ret += *ent;
}
void
km_sum_all_i(int *ent, void *data)
{
	int *ret = (int *)data;
	*ret += *ent;
}
void
km_sum_all_b(bool *ent, void *data)
{
	bool *ret = (bool *)data;
	*ret = XOR(*ret, *ent);
}
void
km_sum_all_v(VALUE *ent, void *data)
{
	VALUE *ret = (VALUE *)data;
	*ret = rb_funcall(*ret, id_op_plus, 1, *ent);
}
void
km_sum_col_d(double *ent, size_t i, size_t j, void *data)
{
	double *r = (double *)data;
	r[j] += *ent;
}
void
km_sum_col_z(COMPLEX *ent, size_t i, size_t j, void *data)
{
	COMPLEX *r = (COMPLEX *)data;
	r[j] += *ent;
}
void
km_sum_col_i(int *ent, size_t i, size_t j, void *data)
{
	int *r = (int *)data;
	r[j] += *ent;
}
void
km_sum_col_b(bool *ent, size_t i, size_t j, void *data)
{
	bool *r = (bool *)data;
	r[j] = XOR(r[j], *ent);
}
void
km_sum_col_v(VALUE *ent, size_t i, size_t j, void *data)
{
	VALUE *r = (VALUE *)data;
	r[j] = rb_funcall(r[j], id_op_plus, 1, *ent);
}
void
km_sum_row_d(double *ent, size_t i, size_t j, void *data)
{
	double *r = (double *)data;
	r[i] += *ent;
}
void
km_sum_row_z(COMPLEX *ent, size_t i, size_t j, void *data)
{
	COMPLEX *r = (COMPLEX *)data;
	r[i] += *ent;
}
void
km_sum_row_i(int *ent, size_t i, size_t j, void *data)
{
	int *r = (int *)data;
	r[i] += *ent;
}
void
km_sum_row_b(bool *ent, size_t i, size_t j, void *data)
{
	bool *r = (bool *)data;
	r[i] = XOR(r[i], *ent);
}
void
km_sum_row_v(VALUE *ent, size_t i, size_t j, void *data)
{
	VALUE *r = (VALUE *)data;
	r[i] = rb_funcall(r[i], id_op_plus, 1, *ent);
}
VALUE
kmm_mat__sum(VALUE self, VALUE dest)
{
	km_check_frozen(dest);
	SMAT *sr = km_mat2smat(dest), *sa = km_mat2smat(self);
	if ( sr->vtype != sa->vtype ) { rb_raise(km_eVT, "value types must be same"); }
	if ( sr->stype != ST_FULL ) { rb_raise(km_eShare, "the argument must be a full matrix"); }
	sr->trans = false; sr->ld = sr->m;
	if ( sr->m == 1 && sr->n == 1 ) {
		VT_SWITCH( sr->vtype,
			sr->dbody[0]=0.0; km_smat_each_d(sa, km_sum_all_d, sr->body); return rb_float_new(sr->dbody[0]);,
			sr->zbody[0]=cpack(0.0, 0.0); km_smat_each_z(sa, km_sum_all_z, sr->body); return km_c2v(sr->zbody[0]);,
			sr->ibody[0]=0; km_smat_each_i(sa, km_sum_all_i, sr->body); return INT2NUM(sr->ibody[0]);,
			sr->bbody[0]=false; km_smat_each_b(sa, km_sum_all_b, sr->body); return TF2V(sr->bbody[0]);,
			sr->vbody[0]=INT2NUM(0); km_smat_each_v(sa, km_sum_all_v, sr->body); return sr->vbody[0];
		);
	} else if ( sr->m == 1 && sr->n == sa->n ) {
		kmm_mat_zero(dest);
		VT_SWITCH( sr->vtype,
			km_smat_each_with_index_d(sa, km_sum_col_d, sr->body);,
			km_smat_each_with_index_z(sa, km_sum_col_z, sr->body);,
			km_smat_each_with_index_i(sa, km_sum_col_i, sr->body);,
			km_smat_each_with_index_b(sa, km_sum_col_b, sr->body);,
			km_smat_each_with_index_v(sa, km_sum_col_v, sr->body);
		);
	} else if ( sr->m == sa->m && sr->n == 1 ) {
		kmm_mat_zero(dest);
		VT_SWITCH( sr->vtype,
			km_smat_each_with_index_d(sa, km_sum_row_d, sr->body);,
			km_smat_each_with_index_z(sa, km_sum_row_z, sr->body);,
			km_smat_each_with_index_i(sa, km_sum_row_i, sr->body);,
			km_smat_each_with_index_b(sa, km_sum_row_b, sr->body);,
			km_smat_each_with_index_v(sa, km_sum_row_v, sr->body);
		);
	} else if ( SAME_SIZE(sr, sa) ) {
		km_smat_copy(sr, sa);
	} else {
		rb_raise(km_eDim, "sum from size (%zu, %zu) to size (%zu, %zu) is not defined", sa->m, sa->n, sr->m, sr->n);
	}
	return dest;
}
