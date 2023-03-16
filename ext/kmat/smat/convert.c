#include "../kmat.h"

// alias to_m!
VALUE
kmm_mat_to_mat_destl(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 0, 1);
	if ( argc == 0 ) {
		return self;
	} else {
		VT_SWITCH( km_sym2vt(argv[0]),
			return kmm_mat_to_fmat_destl(self);,
			return kmm_mat_to_cmat_destl(self);,
			return kmm_mat_to_imat_destl(self);,
			return kmm_mat_to_bmat_destl(self);,
			return kmm_mat_to_omat_destl(self);
		);
	}
}

// alias to_m
VALUE
kmm_mat_to_mat(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 0, 1);
	if ( argc == 0 ) {
		return rb_funcall(self, id_dup, 0);
	} else {
		VT_SWITCH( km_sym2vt(argv[0]),
			return kmm_mat_to_fmat(self);,
			return kmm_mat_to_cmat(self);,
			return kmm_mat_to_imat(self);,
			return kmm_mat_to_bmat(self);,
			return kmm_mat_to_omat(self);
		);
	}
}


struct km_conv_arg {
	union {
		void *body;
		double *dbody;
		COMPLEX *zbody;
		int *ibody;
		bool *bbody;
		VALUE *vbody;
	};
	size_t ldb;
};

static void
z2d_func(COMPLEX *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->dbody)[i+j*(data->ldb)] = creal(*ent);
}
static void
i2d_func(int *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->dbody)[i+j*(data->ldb)] = (double)*ent;
}
static void
b2d_func(bool *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->dbody)[i+j*(data->ldb)] = (*ent ? 1.0 : 0.0);
}
static void
v2d_func(VALUE *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->dbody)[i+j*(data->ldb)] = NUM2DBL(*ent);
}
VALUE
kmm_mat_to_fmat_destl(VALUE self)
{
	km_check_frozen(self);
	SMAT *src = km_mat2smat(self);
	if ( src->vtype == VT_DOUBLE ) {
		return self;
	}
	SMAT *dest = km_smat_alloc(src->m, src->n, VT_DOUBLE);
	struct km_conv_arg data = {{dest->body}, dest->ld};
	VT_SWITCH( src->vtype,
		rb_raise(km_eInternal, "unexpected reach");,
		km_smat_each_with_index_z(src, z2d_func, &data);,
		km_smat_each_with_index_i(src, i2d_func, &data);,
		km_smat_each_with_index_b(src, b2d_func, &data);,
		km_smat_each_with_index_v(src, v2d_func, &data);
	);
	km_smat_copy(src, dest);
	return self;
}
VALUE
kmm_mat_to_fmat(VALUE self)
{
	km_check_frozen(self);
	SMAT *src = km_mat2smat(self);
	if ( src->vtype == VT_DOUBLE ) {
		return rb_funcall(self, id_dup, 0);
	}
	VALUE ret = km_Mat(src->m, src->n, VT_DOUBLE);
	SMAT *dest = km_mat2smat(ret);
	struct km_conv_arg data = {{dest->body}, dest->ld};
	VT_SWITCH( src->vtype,
		rb_raise(km_eInternal, "unexpected reach");,
		km_smat_each_with_index_z(src, z2d_func, &data);,
		km_smat_each_with_index_i(src, i2d_func, &data);,
		km_smat_each_with_index_b(src, b2d_func, &data);,
		km_smat_each_with_index_v(src, v2d_func, &data);
	);
	return ret;
}

static void
d2z_func(double *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->zbody)[i+j*(data->ldb)] = cpack(*ent, 0.0);
}
static void
i2z_func(int *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->zbody)[i+j*(data->ldb)] = cpack((double)*ent, 0.0);
}
static void
b2z_func(bool *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->zbody)[i+j*(data->ldb)] = (*ent ? cpack(1.0, 0.0) : cpack(0.0, 0.0));
}
static void
v2z_func(VALUE *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->zbody)[i+j*(data->ldb)] = km_v2c(*ent);
}
VALUE
kmm_mat_to_cmat_destl(VALUE self)
{
	km_check_frozen(self);
	SMAT *src = km_mat2smat(self);
	if ( src->vtype == VT_COMPLEX ) {
		return self;
	}
	SMAT *dest = km_smat_alloc(src->m, src->n, VT_COMPLEX);
	struct km_conv_arg data = {{dest->body}, dest->ld};
	VT_SWITCH( src->vtype,
		km_smat_each_with_index_d(src, d2z_func, &data);,
		rb_raise(km_eInternal, "unexpected reach");,
		km_smat_each_with_index_i(src, i2z_func, &data);,
		km_smat_each_with_index_b(src, b2z_func, &data);,
		km_smat_each_with_index_v(src, v2z_func, &data);
	);
	km_smat_copy(src, dest);
	return self;
}
VALUE
kmm_mat_to_cmat(VALUE self)
{
	km_check_frozen(self);
	SMAT *src = km_mat2smat(self);
	if ( src->vtype == VT_COMPLEX ) {
		return rb_funcall(self, id_dup, 0);
	}
	VALUE ret = km_Mat(src->m, src->n, VT_COMPLEX);
	SMAT *dest = km_mat2smat(ret);
	struct km_conv_arg data = {{dest->body}, dest->ld};
	VT_SWITCH( src->vtype,
		km_smat_each_with_index_d(src, d2z_func, &data);,
		rb_raise(km_eInternal, "unexpected reach");,
		km_smat_each_with_index_i(src, i2z_func, &data);,
		km_smat_each_with_index_b(src, b2z_func, &data);,
		km_smat_each_with_index_v(src, v2z_func, &data);
	);
	return ret;
}

static void
d2i_func(double *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->ibody)[i+j*(data->ldb)] = (int)(*ent);
}
static void
z2i_func(COMPLEX *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->ibody)[i+j*(data->ldb)] = (int)creal(*ent);
}
static void
b2i_func(bool *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->ibody)[i+j*(data->ldb)] = (*ent ? 1 : 0);
}
static void
v2i_func(VALUE *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->ibody)[i+j*(data->ldb)] = NUM2INT(*ent);
}
VALUE
kmm_mat_to_imat_destl(VALUE self)
{
	km_check_frozen(self);
	SMAT *src = km_mat2smat(self);
	if ( src->vtype == VT_INT ) {
		return self;
	}
	SMAT *dest = km_smat_alloc(src->m, src->n, VT_INT);
	struct km_conv_arg data = {{dest->body}, dest->ld};
	VT_SWITCH( src->vtype,
		km_smat_each_with_index_d(src, d2i_func, &data);,
		km_smat_each_with_index_z(src, z2i_func, &data);,
		rb_raise(km_eInternal, "unexpected reach");,
		km_smat_each_with_index_b(src, b2i_func, &data);,
		km_smat_each_with_index_v(src, v2i_func, &data);
	);
	km_smat_copy(src, dest);
	return self;
}
VALUE
kmm_mat_to_imat(VALUE self)
{
	km_check_frozen(self);
	SMAT *src = km_mat2smat(self);
	if ( src->vtype == VT_INT ) {
		return rb_funcall(self, id_dup, 0);
	}
	VALUE ret = km_Mat(src->m, src->n, VT_INT);
	SMAT *dest = km_mat2smat(ret);
	struct km_conv_arg data = {{dest->body}, dest->ld};
	VT_SWITCH( src->vtype,
		km_smat_each_with_index_d(src, d2i_func, &data);,
		km_smat_each_with_index_z(src, z2i_func, &data);,
		rb_raise(km_eInternal, "unexpected reach");,
		km_smat_each_with_index_b(src, b2i_func, &data);,
		km_smat_each_with_index_v(src, v2i_func, &data);
	);
	return ret;
}

static void
d2b_func(double *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->bbody)[i+j*(data->ldb)] = (*ent != 0.0);
}
static void
z2b_func(COMPLEX *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->bbody)[i+j*(data->ldb)] = (*ent != cpack(0.0, 0.0));
}
static void
i2b_func(int *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->bbody)[i+j*(data->ldb)] = (*ent != 0);
}
static void
v2b_func(VALUE *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->bbody)[i+j*(data->ldb)] = RTEST(*ent);
}
VALUE
kmm_mat_to_bmat_destl(VALUE self)
{
	km_check_frozen(self);
	SMAT *src = km_mat2smat(self);
	if ( src->vtype == VT_BOOL ) {
		return self;
	}
	SMAT *dest = km_smat_alloc(src->m, src->n, VT_BOOL);
	struct km_conv_arg data = {{dest->body}, dest->ld};
	VT_SWITCH( src->vtype,
		km_smat_each_with_index_d(src, d2b_func, &data);,
		km_smat_each_with_index_z(src, z2b_func, &data);,
		km_smat_each_with_index_i(src, i2b_func, &data);,
		rb_raise(km_eInternal, "unexpected reach");,
		km_smat_each_with_index_v(src, v2b_func, &data);
	);
	km_smat_copy(src, dest);
	return self;
}
VALUE
kmm_mat_to_bmat(VALUE self)
{
	km_check_frozen(self);
	SMAT *src = km_mat2smat(self);
	if ( src->vtype == VT_BOOL ) {
		return rb_funcall(self, id_dup, 0);
	}
	VALUE ret = km_Mat(src->m, src->n, VT_BOOL);
	SMAT *dest = km_mat2smat(ret);
	struct km_conv_arg data = {{dest->body}, dest->ld};
	VT_SWITCH( src->vtype,
		km_smat_each_with_index_d(src, d2b_func, &data);,
		km_smat_each_with_index_z(src, z2b_func, &data);,
		km_smat_each_with_index_i(src, i2b_func, &data);,
		rb_raise(km_eInternal, "unexpected reach");,
		km_smat_each_with_index_v(src, v2b_func, &data);
	);
	return ret;
}

static void
d2v_func(double *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->vbody)[i+j*(data->ldb)] = rb_float_new(*ent);
}
static void
z2v_func(COMPLEX *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->vbody)[i+j*(data->ldb)] = km_c2v(*ent);
}
static void
i2v_func(int *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->vbody)[i+j*(data->ldb)] = INT2NUM(*ent);
}
static void
b2v_func(bool *ent, size_t i, size_t j, void *data_)
{
	struct km_conv_arg *data = (struct km_conv_arg *)data_;
	(data->vbody)[i+j*(data->ldb)] = TF2V(*ent);
}
VALUE
kmm_mat_to_omat_destl(VALUE self)
{
	km_check_frozen(self);
	SMAT *src = km_mat2smat(self);
	if ( src->vtype == VT_VALUE ) {
		return self;
	}
	SMAT *dest = km_smat_alloc(src->m, src->n, VT_VALUE);
	struct km_conv_arg data = {{dest->body}, dest->ld};
	VT_SWITCH( src->vtype,
		km_smat_each_with_index_d(src, d2v_func, &data);,
		km_smat_each_with_index_z(src, z2v_func, &data);,
		km_smat_each_with_index_i(src, i2v_func, &data);,
		km_smat_each_with_index_b(src, b2v_func, &data);,
		rb_raise(km_eInternal, "unexpected reach");
	);
	km_smat_copy(src, dest);
	return self;
}
VALUE
kmm_mat_to_omat(VALUE self)
{
	km_check_frozen(self);
	SMAT *src = km_mat2smat(self);
	if ( src->vtype == VT_VALUE ) {
		return rb_funcall(self, id_dup, 0);
	}
	VALUE ret = km_Mat(src->m, src->n, VT_VALUE);
	SMAT *dest = km_mat2smat(ret);
	struct km_conv_arg data = {{dest->body}, dest->ld};
	VT_SWITCH( src->vtype,
		km_smat_each_with_index_d(src, d2v_func, &data);,
		km_smat_each_with_index_z(src, z2v_func, &data);,
		km_smat_each_with_index_i(src, i2v_func, &data);,
		km_smat_each_with_index_b(src, b2v_func, &data);,
		rb_raise(km_eInternal, "unexpected reach");
	);
	return ret;
}
