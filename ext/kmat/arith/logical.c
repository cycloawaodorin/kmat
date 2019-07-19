#include "../kmat.h"

// the number of SMATs is specified by `argc'
// check whether the value types of the arguments are VT_BOOL or not
// raise Mat::ValueTypeError if at least one of the arguments is not a boolean matrix
void
km_check_bool(int argc, ...)
{
	va_list argp;
	va_start(argp, argc);
	for ( int i=0; i<argc; i++ ) {
		SMAT *smat = va_arg(argp, SMAT *);
		if ( smat->vtype != VT_BOOL ) {
			va_end(argp);
			rb_raise(km_eVT, "boolean matrix is expected");
		}
	}
	va_end(argp);
}


// element-wise logical operations
static void
km_not_func(bool *ent, void *null)
{
	*ent = !(*ent);
}
VALUE
kmm_mat_not_dest(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	km_check_bool(1, smat);
	km_smat_each_b(smat, km_not_func, NULL);
	return self;
}

static void
km_and_func(bool *dest, bool *other, void *null)
{
	*dest = ( *dest && *other );
}
VALUE
kmm_mat_and_dest(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *ss = km_mat2smat(self), *so = km_mat2smat(other);
	km_check_bool(2, ss, so);
	CHECK_SAME_SIZE(ss, so);
	km_smat_each2_b(ss, so, km_and_func, NULL);
	return self;
}

static void
km_or_func(bool *dest, bool *other, void *null)
{
	*dest = ( *dest || *other );
}
VALUE
kmm_mat_or_dest(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *ss = km_mat2smat(self), *so = km_mat2smat(other);
	km_check_bool(2, ss, so);
	CHECK_SAME_SIZE(ss, so);
	km_smat_each2_b(ss, so, km_or_func, NULL);
	return self;
}

static void
km_xor_func(bool *dest, bool *other, void *null)
{
	*dest = XOR( *dest, *other );
}
VALUE
kmm_mat_xor_dest(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *ss = km_mat2smat(self), *so = km_mat2smat(other);
	km_check_bool(2, ss, so);
	CHECK_SAME_SIZE(ss, so);
	km_smat_each2_b(ss, so, km_xor_func, NULL);
	return self;
}

// compare operations
static void
km_eq_d(bool *r, const double *a, const double *b, void *null)
{
	*r = (*a == *b);
}
static void
km_eq_z(bool *r, const COMPLEX *a, const COMPLEX *b, void *null)
{
	*r = (*a == *b);
}
static void
km_eq_i(bool *r, const int *a, const int *b, void *null)
{
	*r = (*a == *b);
}
static void
km_eq_b(bool *r, const bool *a, const bool *b, void *null)
{
	*r = (*a == *b);
}
static void
km_eq_v(bool *r, const VALUE *a, const VALUE *b, void *null)
{
	*r = rb_funcall(*a, id_op_eq, 1, *b);
}
VALUE
kmm_mat_eq_destl(VALUE self, VALUE va, VALUE vb)
{
	km_check_frozen(self);
	SMAT *sr=km_mat2smat(self), *sa=km_mat2smat(va), *sb=km_mat2smat(vb);
	km_check_bool(1, sr);
	CHECK_SAME_SIZE(sr, sa);
	CHECK_SAME_SIZE(sr, sb);
	if ( sa->vtype != sb->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	VT_SWITCH( sa->vtype,
		km_smat_each3_bcdcd(sr, sa, sb, km_eq_d, NULL);,
		km_smat_each3_bczcz(sr, sa, sb, km_eq_z, NULL);,
		km_smat_each3_bcici(sr, sa, sb, km_eq_i, NULL);,
		km_smat_each3_bcbcb(sr, sa, sb, km_eq_b, NULL);,
		km_smat_each3_bcvcv(sr, sa, sb, km_eq_v, NULL);
	);
	return self;
}
static void
km_ne_d(bool *r, const double *a, const double *b, void *null)
{
	*r = (*a != *b);
}
static void
km_ne_z(bool *r, const COMPLEX *a, const COMPLEX *b, void *null)
{
	*r = (*a != *b);
}
static void
km_ne_i(bool *r, const int *a, const int *b, void *null)
{
	*r = (*a != *b);
}
static void
km_ne_b(bool *r, const bool *a, const bool *b, void *null)
{
	*r = (*a != *b);
}
static void
km_ne_v(bool *r, const VALUE *a, const VALUE *b, void *null)
{
	*r = rb_funcall(*a, id_op_ne, 1, *b);
}
VALUE
kmm_mat_ne_destl(VALUE self, VALUE va, VALUE vb)
{
	km_check_frozen(self);
	SMAT *sr=km_mat2smat(self), *sa=km_mat2smat(va), *sb=km_mat2smat(vb);
	km_check_bool(1, sr);
	CHECK_SAME_SIZE(sr, sa);
	CHECK_SAME_SIZE(sr, sb);
	if ( sa->vtype != sb->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	VT_SWITCH( sa->vtype,
		km_smat_each3_bcdcd(sr, sa, sb, km_ne_d, NULL);,
		km_smat_each3_bczcz(sr, sa, sb, km_ne_z, NULL);,
		km_smat_each3_bcici(sr, sa, sb, km_ne_i, NULL);,
		km_smat_each3_bcbcb(sr, sa, sb, km_ne_b, NULL);,
		km_smat_each3_bcvcv(sr, sa, sb, km_ne_v, NULL);
	);
	return self;
}
static void
km_lt_d(bool *r, const double *a, const double *b, void *null)
{
	*r = (*a < *b);
}
static void
km_lt_i(bool *r, const int *a, const int *b, void *null)
{
	*r = (*a < *b);
}
static void
km_lt_v(bool *r, const VALUE *a, const VALUE *b, void *null)
{
	*r = rb_funcall(*a, id_op_lt, 1, *b);
}
VALUE
kmm_mat_lt_destl(VALUE self, VALUE va, VALUE vb)
{
	km_check_frozen(self);
	SMAT *sr=km_mat2smat(self), *sa=km_mat2smat(va), *sb=km_mat2smat(vb);
	km_check_bool(1, sr);
	CHECK_SAME_SIZE(sr, sa);
	CHECK_SAME_SIZE(sr, sb);
	if ( sa->vtype != sb->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	if ( sa->vtype == VT_DOUBLE ) {
		km_smat_each3_bcdcd(sr, sa, sb, km_lt_d, NULL);
	} else if ( sa->vtype == VT_INT ) {
		km_smat_each3_bcici(sr, sa, sb, km_lt_i, NULL);
	} else {
		if ( sa->vtype != VT_VALUE ) {
			sa = km_mat2smat(kmm_mat_to_omat(va));
			sb = km_mat2smat(kmm_mat_to_omat(vb));
		}
		km_smat_each3_bcvcv(sr, sa, sb, km_lt_v, NULL);
	}
	return self;
}
static void
km_le_d(bool *r, const double *a, const double *b, void *null)
{
	*r = (*a <= *b);
}
static void
km_le_i(bool *r, const int *a, const int *b, void *null)
{
	*r = (*a <= *b);
}
static void
km_le_v(bool *r, const VALUE *a, const VALUE *b, void *null)
{
	*r = rb_funcall(*a, id_op_le, 1, *b);
}
VALUE
kmm_mat_le_destl(VALUE self, VALUE va, VALUE vb)
{
	km_check_frozen(self);
	SMAT *sr=km_mat2smat(self), *sa=km_mat2smat(va), *sb=km_mat2smat(vb);
	km_check_bool(1, sr);
	CHECK_SAME_SIZE(sr, sa);
	CHECK_SAME_SIZE(sr, sb);
	if ( sa->vtype != sb->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	if ( sa->vtype == VT_DOUBLE ) {
		km_smat_each3_bcdcd(sr, sa, sb, km_le_d, NULL);
	} else if ( sa->vtype == VT_INT ) {
		km_smat_each3_bcici(sr, sa, sb, km_le_i, NULL);
	} else {
		if ( sa->vtype != VT_VALUE ) {
			sa = km_mat2smat(kmm_mat_to_omat(va));
			sb = km_mat2smat(kmm_mat_to_omat(vb));
		}
		km_smat_each3_bcvcv(sr, sa, sb, km_le_v, NULL);
	}
	return self;
}
static void
km_gt_d(bool *r, const double *a, const double *b, void *null)
{
	*r = (*a > *b);
}
static void
km_gt_i(bool *r, const int *a, const int *b, void *null)
{
	*r = (*a > *b);
}
static void
km_gt_v(bool *r, const VALUE *a, const VALUE *b, void *null)
{
	*r = rb_funcall(*a, id_op_gt, 1, *b);
}
VALUE
kmm_mat_gt_destl(VALUE self, VALUE va, VALUE vb)
{
	km_check_frozen(self);
	SMAT *sr=km_mat2smat(self), *sa=km_mat2smat(va), *sb=km_mat2smat(vb);
	km_check_bool(1, sr);
	CHECK_SAME_SIZE(sr, sa);
	CHECK_SAME_SIZE(sr, sb);
	if ( sa->vtype != sb->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	if ( sa->vtype == VT_DOUBLE ) {
		km_smat_each3_bcdcd(sr, sa, sb, km_gt_d, NULL);
	} else if ( sa->vtype == VT_INT ) {
		km_smat_each3_bcici(sr, sa, sb, km_gt_i, NULL);
	} else {
		if ( sa->vtype != VT_VALUE ) {
			sa = km_mat2smat(kmm_mat_to_omat(va));
			sb = km_mat2smat(kmm_mat_to_omat(vb));
		}
		km_smat_each3_bcvcv(sr, sa, sb, km_gt_v, NULL);
	}
	return self;
}
static void
km_ge_d(bool *r, const double *a, const double *b, void *null)
{
	*r = (*a >= *b);
}
static void
km_ge_i(bool *r, const int *a, const int *b, void *null)
{
	*r = (*a >= *b);
}
static void
km_ge_v(bool *r, const VALUE *a, const VALUE *b, void *null)
{
	*r = rb_funcall(*a, id_op_ge, 1, *b);
}
VALUE
kmm_mat_ge_destl(VALUE self, VALUE va, VALUE vb)
{
	km_check_frozen(self);
	SMAT *sr=km_mat2smat(self), *sa=km_mat2smat(va), *sb=km_mat2smat(vb);
	km_check_bool(1, sr);
	CHECK_SAME_SIZE(sr, sa);
	CHECK_SAME_SIZE(sr, sb);
	if ( sa->vtype != sb->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	if ( sa->vtype == VT_DOUBLE ) {
		km_smat_each3_bcdcd(sr, sa, sb, km_ge_d, NULL);
	} else if ( sa->vtype == VT_INT ) {
		km_smat_each3_bcici(sr, sa, sb, km_ge_i, NULL);
	} else {
		if ( sa->vtype != VT_VALUE ) {
			sa = km_mat2smat(kmm_mat_to_omat(va));
			sb = km_mat2smat(kmm_mat_to_omat(vb));
		}
		km_smat_each3_bcvcv(sr, sa, sb, km_ge_v, NULL);
	}
	return self;
}
