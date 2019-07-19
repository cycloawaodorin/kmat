#include "../kmat.h"

// element-wise arithmetic operations
static void
km_add_d(double *y, double *x, void *null)
{
	*y += *x;
}
static void
km_add_z(COMPLEX *y, COMPLEX *x, void *null)
{
	*y += *x;
}
static void
km_add_i(int *y, int *x, void *null)
{
	*y += *x;
}
static void
km_add_b(bool *y, bool *x, void *null)
{
	*y = XOR(*y, *x);
}
static void
km_add_v(VALUE *y, VALUE *x, void *null)
{
	*y = rb_funcall(*y, id_op_plus, 1, *x);
}
VALUE
kmm_mat_add_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sy=km_mat2smat(self), *sx=km_mat2smat(other);
	if ( sy->vtype != sx->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	CHECK_SAME_SIZE(sy, sx);
	if ( sy->vtype == VT_DOUBLE ) {
		if ( sy->stype == ST_FULL && sx->stype == ST_FULL && sy->trans == sx->trans ) {
			double alpha = 1.0; int len = LENGTH(sy), ione=1;
			daxpy_(&len, &alpha, sx->dbody, &ione, sy->dbody, &ione);
		} else if ( sy->stype != ST_RSUB && sx->stype != ST_RSUB ) {
			double alpha = 1.0; int incy, stepy, incx, stepx;
			if ( sy->trans ) {
				incy = sy->ld; stepy = 1;
			} else {
				incy = 1; stepy = sy->ld;
			}
			if ( sx->trans ) {
				incx = sx->ld; stepx = 1;
			} else {
				incx = 1; stepx = sx->ld;
			}
			for ( int i=0; i<(sy->n); i++ ) {
				daxpy_(&(sy->m), &alpha, (sx->dbody)+(i*stepx), &incx, (sy->dbody)+(i*stepy), &incy);
			}
		} else {
			km_smat_each2_d(sy, sx, km_add_d, NULL);
		}
	} else if ( sy->vtype == VT_COMPLEX ) {
		if ( sy->stype == ST_FULL && sx->stype == ST_FULL && sy->trans == sx->trans ) {
			COMPLEX alpha = cpack(1.0, 0.0); int len = LENGTH(sy), ione=1;
			zaxpy_(&len, &alpha, sx->zbody, &ione, sy->zbody, &ione);
		} else if ( sy->stype != ST_RSUB && sx->stype != ST_RSUB ) {
			COMPLEX alpha = cpack(1.0, 0.0); int incy, stepy, incx, stepx;
			if ( sy->trans ) {
				incy = sy->ld; stepy = 1;
			} else {
				incy = 1; stepy = sy->ld;
			}
			if ( sx->trans ) {
				incx = sx->ld; stepx = 1;
			} else {
				incx = 1; stepx = sx->ld;
			}
			for ( int i=0; i<(sy->n); i++ ) {
				zaxpy_(&(sy->m), &alpha, (sx->zbody)+(i*stepx), &incx, (sy->zbody)+(i*stepy), &incy);
			}
		} else {
			km_smat_each2_z(sy, sx, km_add_z, NULL);
		}
	} else if ( sy->vtype == VT_INT ) {
		km_smat_each2_i(sy, sx, km_add_i, NULL);
	} else if ( sy->vtype == VT_BOOL ) {
		km_smat_each2_b(sy, sx, km_add_b, NULL);
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each2_v(sy, sx, km_add_v, NULL);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

static void
km_sub_d(double *y, double *x, void *null)
{
	*y -= *x;
}
static void
km_sub_z(COMPLEX *y, COMPLEX *x, void *null)
{
	*y -= *x;
}
static void
km_sub_i(int *y, int *x, void *null)
{
	*y -= *x;
}
static void
km_sub_v(VALUE *y, VALUE *x, void *null)
{
	*y = rb_funcall(*y, id_op_minus, 1, *x);
}
VALUE
kmm_mat_sub_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sy=km_mat2smat(self), *sx=km_mat2smat(other);
	if ( sy->vtype != sx->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	CHECK_SAME_SIZE(sy, sx);
	if ( sy->vtype == VT_DOUBLE ) {
		if ( sy->stype == ST_FULL && sx->stype == ST_FULL && sy->trans == sx->trans ) {
			double alpha = -1.0; int len = LENGTH(sy), ione=1;
			daxpy_(&len, &alpha, sx->dbody, &ione, sy->dbody, &ione);
		} else if ( sy->stype != ST_RSUB && sx->stype != ST_RSUB ) {
			double alpha = -1.0; int incy, stepy, incx, stepx;
			if ( sy->trans ) {
				incy = sy->ld; stepy = 1;
			} else {
				incy = 1; stepy = sy->ld;
			}
			if ( sx->trans ) {
				incx = sx->ld; stepx = 1;
			} else {
				incx = 1; stepx = sx->ld;
			}
			for ( int i=0; i<(sy->n); i++ ) {
				daxpy_(&(sy->m), &alpha, (sx->dbody)+(i*stepx), &incx, (sy->dbody)+(i*stepy), &incy);
			}
		} else {
			km_smat_each2_d(sy, sx, km_sub_d, NULL);
		}
	} else if ( sy->vtype == VT_COMPLEX ) {
		if ( sy->stype == ST_FULL && sx->stype == ST_FULL && sy->trans == sx->trans ) {
			COMPLEX alpha = cpack(-1.0, 0.0); int len = LENGTH(sy), ione=1;
			zaxpy_(&len, &alpha, sx->zbody, &ione, sy->zbody, &ione);
		} else if ( sy->stype != ST_RSUB && sx->stype != ST_RSUB ) {
			COMPLEX alpha = cpack(-1.0, 0.0); int incy, stepy, incx, stepx;
			if ( sy->trans ) {
				incy = sy->ld; stepy = 1;
			} else {
				incy = 1; stepy = sy->ld;
			}
			if ( sx->trans ) {
				incx = sx->ld; stepx = 1;
			} else {
				incx = 1; stepx = sx->ld;
			}
			for ( int i=0; i<(sy->n); i++ ) {
				zaxpy_(&(sy->m), &alpha, (sx->zbody)+(i*stepx), &incx, (sy->zbody)+(i*stepy), &incy);
			}
		} else {
			km_smat_each2_z(sy, sx, km_sub_z, NULL);
		}
	} else if ( sy->vtype == VT_INT ) {
		km_smat_each2_i(sy, sx, km_sub_i, NULL);
	} else if ( sy->vtype == VT_BOOL ) {
		km_smat_each2_b(sy, sx, km_add_b, NULL); // subtract is the same as add for boolean
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each2_v(sy, sx, km_sub_v, NULL);
	} else {
		rb_raise(km_eInternal, "unkown value type");
	}
	return self;
}

static void
km_mul_d(double *y, double *x, void *null)
{
	*y *= *x;
}
static void
km_mul_z(COMPLEX *y, COMPLEX *x, void *null)
{
	*y *= *x;
}
static void
km_mul_i(int *y, int *x, void *null)
{
	*y *= *x;
}
static void
km_mul_b(bool *y, bool *x, void *null)
{
	*y = ( *y && *x );
}
static void
km_mul_v(VALUE *y, VALUE *x, void *null)
{
	*y = rb_funcall(*y, id_op_mul, 1, *x);
}
VALUE
kmm_mat_e_mul_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sy=km_mat2smat(self), *sx=km_mat2smat(other);
	if ( sy->vtype != sx->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	CHECK_SAME_SIZE(sy, sx);
	if ( sy->vtype == VT_DOUBLE ) {
		km_smat_each2_d(sy, sx, km_mul_d, NULL);
	} else if ( sy->vtype == VT_COMPLEX ) {
		km_smat_each2_z(sy, sx, km_mul_z, NULL);
	} else if ( sy->vtype == VT_INT ) {
		km_smat_each2_i(sy, sx, km_mul_i, NULL);
	} else if ( sy->vtype == VT_BOOL ) {
		km_smat_each2_b(sy, sx, km_mul_b, NULL);
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each2_v(sy, sx, km_mul_v, NULL);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

static void
km_div_d(double *y, double *x, void *null)
{
	if ( *x == 0.0 ) {
		rb_raise(rb_eZeroDivError, "divided by 0");
	}
	*y /= *x;
}
static void
km_div_z(COMPLEX *y, COMPLEX *x, void *null)
{
	if ( *x == cpack(0.0, 0.0) ) {
		rb_raise(rb_eZeroDivError, "divided by 0");
	}
	*y /= *x;
}
static void
km_div_i(int *y, int *x, void *null)
{
	if ( *x == 0 ) {
		rb_raise(rb_eZeroDivError, "divided by 0");
	}
	*y /= *x;
}
static void
km_div_v(VALUE *y, VALUE *x, void *null)
{
	*y = rb_funcall(*y, id_op_div, 1, *x);
}
VALUE
kmm_mat_e_div_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sy=km_mat2smat(self), *sx=km_mat2smat(other);
	if ( sy->vtype != sx->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	CHECK_SAME_SIZE(sy, sx);
	if ( sy->vtype == VT_DOUBLE ) {
		km_smat_each2_d(sy, sx, km_div_d, NULL);
	} else if ( sy->vtype == VT_COMPLEX ) {
		km_smat_each2_z(sy, sx, km_div_z, NULL);
	} else if ( sy->vtype == VT_INT ) {
		km_smat_each2_i(sy, sx, km_div_i, NULL);
	} else if ( sy->vtype == VT_BOOL ) {
		rb_raise(km_eVT, "division for boolean is meaningless");
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each2_v(sy, sx, km_div_v, NULL);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

static void
km_s_add_d(double *ent, void *data)
{
	double *val = (double *)data;
	*ent += *val;
}
static void
km_s_add_z(COMPLEX *ent, void *data)
{
	COMPLEX *val = (COMPLEX *)data;
	*ent += *val;
}
static void
km_s_add_i(int *ent, void *data)
{
	int *val = (int *)data;
	*ent += *val;
}
static void
km_s_add_b(bool *ent, void *data)
{
	bool *val = (bool *)data;
	*ent = XOR(*ent, *val);
}
static void
km_s_add_v(VALUE *ent, void *data)
{
	VALUE val = (VALUE)data;
	*ent = rb_funcall(*ent, id_op_plus, 1, val);
}
VALUE
kmm_mat_s_add_destl(VALUE self, VALUE vval)
{
	km_check_frozen(self);
	SMAT *sy = km_mat2smat(self);
	if ( sy->vtype == VT_DOUBLE ) {
		double val = NUM2DBL(vval);
		if ( sy->stype == ST_FULL ) {
			double alpha=1.0;
			int len=LENGTH(sy), ione=1, izero=0;
			daxpy_(&len, &alpha, &val, &izero, sy->dbody, &ione);
		} else if ( sy->stype == ST_SSUB ) {
			double alpha=1.0;
			int ione=1, izero=0;
			if ( sy->trans ) {
				for ( int i=0; i<(sy->m); i++ ) {
					daxpy_(&(sy->n), &alpha, &val, &izero, (sy->dbody)+(i*(sy->ld)), &ione);
				}
			} else {
				for ( int i=0; i<(sy->n); i++ ) {
					daxpy_(&(sy->m), &alpha, &val, &izero, (sy->dbody)+(i*(sy->ld)), &ione);
				}
			}
		} else {
			km_smat_each_d(sy, km_s_add_d, (void *)&val);
		}
	} else if ( sy->vtype == VT_COMPLEX ) {
		COMPLEX val = km_v2c(vval);
		if ( sy->stype == ST_FULL ) {
			COMPLEX alpha=cpack(1.0, 0.0);
			int len=LENGTH(sy), ione=1, izero=0;
			zaxpy_(&len, &alpha, &val, &izero, sy->zbody, &ione);
		} else if ( sy->stype == ST_SSUB ) {
			COMPLEX alpha=cpack(1.0, 0.0);
			int ione=1, izero=0;
			if ( sy->trans ) {
				for ( int i=0; i<(sy->m); i++ ) {
					zaxpy_(&(sy->n), &alpha, &val, &izero, (sy->zbody)+(i*(sy->ld)), &ione);
				}
			} else {
				for ( int i=0; i<(sy->n); i++ ) {
					zaxpy_(&(sy->m), &alpha, &val, &izero, (sy->zbody)+(i*(sy->ld)), &ione);
				}
			}
		} else {
			km_smat_each_z(sy, km_s_add_z, (void *)&val);
		}
	} else if ( sy->vtype == VT_INT ) {
		int val = NUM2INT(vval);
		km_smat_each_i(sy, km_s_add_i, (void *)&val);
	} else if ( sy->vtype == VT_BOOL ) {
		bool val = RTEST(vval);
		km_smat_each_b(sy, km_s_add_b, (void *)&val);
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each_v(sy, km_s_add_v, (void *)vval);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

static void
km_s_sub_d(double *ent, void *data)
{
	double *val = (double *)data;
	*ent -= *val;
}
static void
km_s_sub_z(COMPLEX *ent, void *data)
{
	COMPLEX *val = (COMPLEX *)data;
	*ent -= *val;
}
static void
km_s_sub_i(int *ent, void *data)
{
	int *val = (int *)data;
	*ent -= *val;
}
static void
km_s_sub_v(VALUE *ent, void *data)
{
	VALUE val = (VALUE)data;
	*ent = rb_funcall(*ent, id_op_minus, 1, val);
}
VALUE
kmm_mat_s_sub_destl(VALUE self, VALUE vval)
{
	km_check_frozen(self);
	SMAT *sy = km_mat2smat(self);
	if ( sy->vtype == VT_DOUBLE ) {
		double val = NUM2DBL(vval);
		if ( sy->stype == ST_FULL ) {
			double alpha=-1.0;
			int len=LENGTH(sy), ione=1, izero=0;
			daxpy_(&len, &alpha, &val, &izero, sy->dbody, &ione);
		} else if ( sy->stype == ST_SSUB ) {
			double alpha=-1.0;
			int ione=1, izero=0;
			if ( sy->trans ) {
				for ( int i=0; i<(sy->m); i++ ) {
					daxpy_(&(sy->n), &alpha, &val, &izero, (sy->dbody)+(i*(sy->ld)), &ione);
				}
			} else {
				for ( int i=0; i<(sy->n); i++ ) {
					daxpy_(&(sy->m), &alpha, &val, &izero, (sy->dbody)+(i*(sy->ld)), &ione);
				}
			}
		} else {
			km_smat_each_d(sy, km_s_sub_d, (void *)&val);
		}
	} else if ( sy->vtype == VT_COMPLEX ) {
		COMPLEX val = km_v2c(vval);
		if ( sy->stype == ST_FULL ) {
			COMPLEX alpha=cpack(-1.0, 0.0);
			int len=LENGTH(sy), ione=1, izero=0;
			zaxpy_(&len, &alpha, &val, &izero, sy->zbody, &ione);
		} else if ( sy->stype == ST_SSUB ) {
			COMPLEX alpha=cpack(-1.0, 0.0);
			int ione=1, izero=0;
			if ( sy->trans ) {
				for ( int i=0; i<(sy->m); i++ ) {
					zaxpy_(&(sy->n), &alpha, &val, &izero, (sy->zbody)+(i*(sy->ld)), &ione);
				}
			} else {
				for ( int i=0; i<(sy->n); i++ ) {
					zaxpy_(&(sy->m), &alpha, &val, &izero, (sy->zbody)+(i*(sy->ld)), &ione);
				}
			}
		} else {
			km_smat_each_z(sy, km_s_sub_z, (void *)&val);
		}
	} else if ( sy->vtype == VT_INT ) {
		int val = NUM2INT(vval);
		km_smat_each_i(sy, km_s_sub_i, (void *)&val);
	} else if ( sy->vtype == VT_BOOL ) {
		bool val = RTEST(vval);
		km_smat_each_b(sy, km_s_add_b, (void *)&val); // subtract is the same as add for boolean
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each_v(sy, km_s_sub_v, (void *)vval);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}


static void
km_s_mul_d(double *ent, void *data)
{
	double *val = (double *)data;
	*ent *= *val;
}
static void
km_s_mul_z(COMPLEX *ent, void *data)
{
	COMPLEX *val = (COMPLEX *)data;
	*ent *= *val;
}
static void
km_s_mul_i(int *ent, void *data)
{
	int *val = (int *)data;
	*ent *= *val;
}
static void
km_s_mul_v(VALUE *ent, void *data)
{
	VALUE val = (VALUE)data;
	*ent = rb_funcall(*ent, id_op_mul, 1, val);
}
VALUE
kmm_mat_s_mul_destl(VALUE self, VALUE vval)
{
	km_check_frozen(self);
	SMAT *sy = km_mat2smat(self);
	if ( sy->vtype == VT_DOUBLE ) {
		double val = NUM2DBL(vval);
		if ( sy->stype == ST_FULL ) {
			int len=LENGTH(sy), ione=1;
			dscal_(&len, &val, sy->dbody, &ione);
		} else if ( sy->stype == ST_SSUB ) {
			int ione=1;
			if ( sy->trans ) {
				for ( int i=0; i<(sy->m); i++ ) {
					dscal_(&(sy->n), &val, (sy->dbody)+(i*(sy->ld)), &ione);
				}
			} else {
				for ( int i=0; i<(sy->n); i++ ) {
					dscal_(&(sy->m), &val, (sy->dbody)+(i*(sy->ld)), &ione);
				}
			}
		} else {
			km_smat_each_d(sy, km_s_mul_d, (void *)&val);
		}
	} else if ( sy->vtype == VT_COMPLEX ) {
		COMPLEX val = km_v2c(vval);
		if ( sy->stype == ST_FULL ) {
			int len=LENGTH(sy), ione=1;
			zscal_(&len, &val, sy->zbody, &ione);
		} else if ( sy->stype == ST_SSUB ) {
			int ione=1;
			if ( sy->trans ) {
				for ( int i=0; i<(sy->m); i++ ) {
					zscal_(&(sy->n), &val, (sy->zbody)+(i*(sy->ld)), &ione);
				}
			} else {
				for ( int i=0; i<(sy->n); i++ ) {
					zscal_(&(sy->m), &val, (sy->zbody)+(i*(sy->ld)), &ione);
				}
			}
		} else {
			km_smat_each_z(sy, km_s_mul_z, (void *)&val);
		}
	} else if ( sy->vtype == VT_INT ) {
		int val = NUM2INT(vval);
		km_smat_each_i(sy, km_s_mul_i, (void *)&val);
	} else if ( sy->vtype == VT_BOOL ) {
		if ( RTEST(vval) ) {
			// nothing to do
		} else {
			kmm_mat_fill(self, Qfalse);
		}
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each_v(sy, km_s_mul_v, (void *)vval);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

static void
km_s_div_i(int *ent, void *data)
{
	int *val = (int *)data;
	*ent /= *val;
}
static void
km_s_div_v(VALUE *ent, void *data)
{
	VALUE val = (VALUE)data;
	*ent = rb_funcall(*ent, id_op_div, 1, val);
}
VALUE
kmm_mat_s_div_destl(VALUE self, VALUE vval)
{
	km_check_frozen(self);
	SMAT *sy = km_mat2smat(self);
	if ( sy->vtype == VT_DOUBLE ) {
		double val = NUM2DBL(vval);
		if ( val == 0.0 ) { rb_raise(rb_eZeroDivError, "divided by 0"); }
		val = 1.0/val;
		if ( sy->stype == ST_FULL ) {
			int len=LENGTH(sy), ione=1;
			dscal_(&len, &val, sy->dbody, &ione);
		} else if ( sy->stype == ST_SSUB ) {
			int ione=1;
			if ( sy->trans ) {
				for ( int i=0; i<(sy->m); i++ ) {
					dscal_(&(sy->n), &val, (sy->dbody)+(i*(sy->ld)), &ione);
				}
			} else {
				for ( int i=0; i<(sy->n); i++ ) {
					dscal_(&(sy->m), &val, (sy->dbody)+(i*(sy->ld)), &ione);
				}
			}
		} else {
			km_smat_each_d(sy, km_s_mul_d, (void *)&val);
		}
	} else if ( sy->vtype == VT_COMPLEX ) {
		COMPLEX val = km_v2c(vval);
		if ( val == cpack(0.0, 0.0) ) { rb_raise(rb_eZeroDivError, "divided by 0"); }
		val = cpack(1.0, 0.0)/val;
		if ( sy->stype == ST_FULL ) {
			int len=LENGTH(sy), ione=1;
			zscal_(&len, &val, sy->zbody, &ione);
		} else if ( sy->stype == ST_SSUB ) {
			int ione=1;
			if ( sy->trans ) {
				for ( int i=0; i<(sy->m); i++ ) {
					zscal_(&(sy->n), &val, (sy->zbody)+(i*(sy->ld)), &ione);
				}
			} else {
				for ( int i=0; i<(sy->n); i++ ) {
					zscal_(&(sy->m), &val, (sy->zbody)+(i*(sy->ld)), &ione);
				}
			}
		} else {
			km_smat_each_z(sy, km_s_mul_z, (void *)&val);
		}
	} else if ( sy->vtype == VT_INT ) {
		int val = NUM2INT(vval);
		if ( val == 0 ) { rb_raise(rb_eZeroDivError, "divided by 0"); }
		km_smat_each_i(sy, km_s_div_i, (void *)&val);
	} else if ( sy->vtype == VT_BOOL ) {
		rb_raise(km_eVT, "division for boolean is meaningless");
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each_v(sy, km_s_div_v, (void *)vval);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

static void
km_add_times_d(double *y, double *x, void *data)
{
	double *a = (double *)data;
	*y += ( (*x) * (*a) );
}
static void
km_add_times_z(COMPLEX *y, COMPLEX *x, void *data)
{
	COMPLEX *a = (COMPLEX *)data;
	*y += ( (*x) * (*a) );
}
static void
km_add_times_i(int *y, int *x, void *data)
{
	int *a = (int *)data;
	*y += ( (*x) * (*a) );
}
static void
km_add_times_v(VALUE *y, VALUE *x, void *data)
{
	VALUE a = (VALUE)data;
	a = rb_funcall(*x, id_op_mul, 1, a);
	*y = rb_funcall(*y, id_op_plus, 1, a);
}
VALUE
kmm_mat_add_times_destl(VALUE self, VALUE other, VALUE valpha)
{
	km_check_frozen(self);
	SMAT *sy=km_mat2smat(self), *sx=km_mat2smat(other);
	if ( sy->vtype != sx->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	CHECK_SAME_SIZE(sx, sy);
	if ( sy->vtype == VT_DOUBLE ) {
		double alpha = NUM2DBL(valpha);
		if ( sy->stype == ST_FULL && sx->stype == ST_FULL && sy->trans == sx->trans ) {
			int len = LENGTH(sy), ione=1;
			daxpy_(&len, &alpha, sx->dbody, &ione, sy->dbody, &ione);
		} else if ( sy->stype != ST_RSUB &&  sx->stype != ST_RSUB ) {
			int incy, stepy, incx, stepx;
			if ( sy->trans ) {
				incy = sy->ld; stepy = 1;
			} else {
				incy = 1; stepy = sy->ld;
			}
			if ( sx->trans ) {
				incx = sx->ld; stepx = 1;
			} else {
				incx = 1; stepx = sx->ld;
			}
			for ( int i=0; i<(sy->n); i++ ) {
				daxpy_(&(sy->m), &alpha, (sx->dbody)+(i*stepx), &incx, (sy->dbody)+(i*stepy), &incy);
			}
		} else {
			km_smat_each2_d(sy, sx, km_add_times_d, (void *)&alpha);
		}
	} else if ( sy->vtype == VT_COMPLEX ) {
		COMPLEX alpha = km_v2c(valpha);
		if ( sy->stype == ST_FULL && sx->stype == ST_FULL && sy->trans == sx->trans ) {
			int len = LENGTH(sy), ione=1;
			zaxpy_(&len, &alpha, sx->zbody, &ione, sy->zbody, &ione);
		} else if ( sy->stype != ST_RSUB &&  sx->stype != ST_RSUB ) {
			int incy, stepy, incx, stepx;
			if ( sy->trans ) {
				incy = sy->ld; stepy = 1;
			} else {
				incy = 1; stepy = sy->ld;
			}
			if ( sx->trans ) {
				incx = sx->ld; stepx = 1;
			} else {
				incx = 1; stepx = sx->ld;
			}
			for ( int i=0; i<(sy->n); i++ ) {
				zaxpy_(&(sy->m), &alpha, (sx->zbody)+(i*stepx), &incx, (sy->zbody)+(i*stepy), &incy);
			}
		} else {
			km_smat_each2_z(sy, sx, km_add_times_z, (void *)&alpha);
		}
	} else if ( sy->vtype == VT_INT ) {
		int alpha = NUM2INT(valpha);
		km_smat_each2_i(sy, sx, km_add_times_i, (void *)&alpha);
	} else if ( sy->vtype == VT_BOOL ) {
		if ( RTEST(valpha) ) {
			kmm_mat_add_destl(self, other);
		} else {
			// nothing to do
		}
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each2_v(sy, sx, km_add_times_v, (void *)valpha);
	} else {
		rb_raise(km_eVT, "unknown value type");
	}
	return self;
}

// element-wise max/min operations of two matricies
static void
km_maximum_d(double *ey, double *ex, void *null)
{
	if ( *ey < *ex ) {
		*ey = *ex;
	}
}
static void
km_maximum_i(int *ey, int *ex, void *null)
{
	if ( *ey < *ex ) {
		*ey = *ex;
	}
}
static void
km_maximum_b(bool *ey, bool *ex, void *null)
{
	if ( *ex ) {
		*ey = true;
	}
}
static void
km_maximum_v(VALUE *ey, VALUE *ex, void *null)
{
	if ( rb_funcall(*ey, id_op_lt, 1, *ex) ) {
		*ey = *ex;
	}
}
VALUE
kmm_mat_maximum_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sy = km_mat2smat(self), *sx = km_mat2smat(other);
	if ( sy->vtype != sx->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	CHECK_SAME_SIZE(sx, sy);
	if ( sy->vtype == VT_DOUBLE ) {
		km_smat_each2_d(sy, sx, km_maximum_d, NULL);
	} else if ( sy->vtype == VT_COMPLEX ) {
		VALUE oself = kmm_mat_to_omat(self);
		sy = km_mat2smat(oself);
		sx = km_mat2smat(kmm_mat_to_omat(other));
		km_smat_each2_v(sy, sx, km_maximum_v, NULL);
		kmm_mat_copy_from(self, oself);
	} else if ( sy->vtype == VT_INT ) {
		km_smat_each2_i(sy, sx, km_maximum_i, NULL);
	} else if ( sy->vtype == VT_BOOL ) {
		km_smat_each2_b(sy, sx, km_maximum_b, NULL);
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each2_v(sy, sx, km_maximum_v, NULL);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

static void
km_s_maximum_d(double *ey, void *data)
{
	double *val = (double *)data;
	if ( *ey < *val ) {
		*ey = *val;
	}
}
static void
km_s_maximum_i(int *ey, void *data)
{
	int *val = (int *)data;
	if ( *ey < *val ) {
		*ey = *val;
	}
}
static void
km_s_maximum_b(bool *ey, void *data)
{
	bool *val = (bool *)data;
	if ( *val ) {
		*ey = true;
	}
}
static void
km_s_maximum_v(VALUE *ey, void *data)
{
	VALUE val = (VALUE)data;
	if ( rb_funcall(*ey, id_op_lt, 1, val) ) {
		*ey = val;
	}
}
VALUE
kmm_mat_s_maximum_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sy = km_mat2smat(self);
	if ( sy->vtype == VT_DOUBLE ) {
		double data = NUM2DBL(other);
		km_smat_each_d(sy, km_s_maximum_d, (void *)&data);
	} else if ( sy->vtype == VT_COMPLEX ) {
		VALUE oself = kmm_mat_to_omat(self);
		sy = km_mat2smat(oself);
		km_smat_each_v(sy, km_s_maximum_v, (void *)other);
		kmm_mat_copy_from(self, oself);
	} else if ( sy->vtype == VT_INT ) {
		int data = NUM2INT(other);
		km_smat_each_i(sy, km_s_maximum_i, (void *)&data);
	} else if ( sy->vtype == VT_BOOL ) {
		bool data = RTEST(other);
		km_smat_each_b(sy, km_s_maximum_b, (void *)&data);
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each_v(sy, km_s_maximum_v, (void *)other);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

static void
km_minimum_d(double *ey, double *ex, void *null)
{
	if ( *ex < *ey ) {
		*ey = *ex;
	}
}
static void
km_minimum_i(int *ey, int *ex, void *null)
{
	if ( *ex < *ey ) {
		*ey = *ex;
	}
}
static void
km_minimum_b(bool *ey, bool *ex, void *null)
{
	if ( !*ex ) {
		*ey = false;
	}
}
static void
km_minimum_v(VALUE *ey, VALUE *ex, void *null)
{
	if ( rb_funcall(*ex, id_op_lt, 1, *ey) ) {
		*ey = *ex;
	}
}
VALUE
kmm_mat_minimum_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sy = km_mat2smat(self), *sx = km_mat2smat(other);
	if ( sy->vtype != sx->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	CHECK_SAME_SIZE(sx, sy);
	if ( sy->vtype == VT_DOUBLE ) {
		km_smat_each2_d(sy, sx, km_minimum_d, NULL);
	} else if ( sy->vtype == VT_COMPLEX ) {
		VALUE oself = kmm_mat_to_omat(self);
		sy = km_mat2smat(oself);
		sx = km_mat2smat(kmm_mat_to_omat(other));
		km_smat_each2_v(sy, sx, km_minimum_v, NULL);
		kmm_mat_copy_from(self, oself);
	} else if ( sy->vtype == VT_INT ) {
		km_smat_each2_i(sy, sx, km_minimum_i, NULL);
	} else if ( sy->vtype == VT_BOOL ) {
		km_smat_each2_b(sy, sx, km_minimum_b, NULL);
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each2_v(sy, sx, km_minimum_v, NULL);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

static void
km_s_minimum_d(double *ey, void *data)
{
	double *val = (double *)data;
	if ( *val < *ey ) {
		*ey = *val;
	}
}
static void
km_s_minimum_i(int *ey, void *data)
{
	int *val = (int *)data;
	if ( *val < *ey ) {
		*ey = *val;
	}
}
static void
km_s_minimum_b(bool *ey, void *data)
{
	bool *val = (bool *)data;
	if ( !*val ) {
		*ey = false;
	}
}
static void
km_s_minimum_v(VALUE *ey, void *data)
{
	VALUE val = (VALUE)data;
	if ( rb_funcall(val, id_op_lt, 1, *ey) ) {
		*ey = val;
	}
}
VALUE
kmm_mat_s_minimum_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sy = km_mat2smat(self);
	if ( sy->vtype == VT_DOUBLE ) {
		double data = NUM2DBL(other);
		km_smat_each_d(sy, km_s_minimum_d, (void *)&data);
	} else if ( sy->vtype == VT_COMPLEX ) {
		VALUE oself = kmm_mat_to_omat(self);
		sy = km_mat2smat(oself);
		km_smat_each_v(sy, km_s_minimum_v, (void *)other);
		kmm_mat_copy_from(self, oself);
	} else if ( sy->vtype == VT_INT ) {
		int data = NUM2INT(other);
		km_smat_each_i(sy, km_s_minimum_i, (void *)&data);
	} else if ( sy->vtype == VT_BOOL ) {
		bool data = RTEST(other);
		km_smat_each_b(sy, km_s_minimum_b, (void *)&data);
	} else if ( sy->vtype == VT_VALUE ) {
		km_smat_each_v(sy, km_s_minimum_v, (void *)other);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

// pow, hypot
static void
km_pow_d(double *b, double *e, void *null)
{
	*b = pow(*b, *e);
}
static void
km_pow_z(COMPLEX *b, COMPLEX *e, void *null)
{
	*b = cpow(*b, *e);
}
static void
km_pow_i(int *b, int *e, void *null)
{
	if ( *e < 0 ) {
		rb_raise(rb_const_get(rb_mMath, id_DomainError), "the exponent must be non-negative");
	}
	int a = *b;
	*b = 1;
	for (int i=*e; i>0; i>>=1) {
		if ( i & 1 ) {
			*b *= a;
		}
		a *= a;
	}
}
static void
km_pow_b(bool *b, bool *e, void *null)
{
	*b = ( *b || !*e );
}
static void
km_pow_v(VALUE *b, VALUE *e, void *null)
{
	*b = rb_funcall(*b, id_op_pow, 1, *e);
}
VALUE
kmm_mat_pow_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sb = km_mat2smat(self), *se = km_mat2smat(other);
	if ( sb->vtype != se->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	CHECK_SAME_SIZE(sb, se);
	VT_SWITCH( sb->vtype,
		km_smat_each2_d(sb, se, km_pow_d, NULL);,
		km_smat_each2_z(sb, se, km_pow_z, NULL);,
		km_smat_each2_i(sb, se, km_pow_i, NULL);,
		km_smat_each2_b(sb, se, km_pow_b, NULL);,
		km_smat_each2_v(sb, se, km_pow_v, NULL);
	);
	return self;
}

static void
km_s_pow_d(double *b, void *data)
{
	double *e = (double *)data;
	*b = pow(*b, *e);
}
static void
km_s_pow_z(COMPLEX *b, void *data)
{
	COMPLEX *e = (COMPLEX *)data;
	*b = cpow(*b, *e);
}
static void
km_s_pow_i(int *b, void *data)
{
	int *e = (int *)data;
	if ( *e < 0 ) {
		rb_raise(rb_const_get(rb_mMath, id_DomainError), "the exponent must be non-negative");
	}
	int a = *b;
	*b = 1;
	for (int i=*e; i>0; i>>=1) {
		if ( i & 1 ) {
			*b *= a;
		}
		a *= a;
	}
}
static void
km_s_pow_b(bool *b, void *data)
{
	bool *e = (bool *)data;
	*b = ( *b || !*e );
}
static void
km_s_pow_v(VALUE *b, void *data)
{
	VALUE e = (VALUE)data;
	*b = rb_funcall(*b, id_op_pow, 1, e);
}
VALUE
kmm_mat_s_pow_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sb = km_mat2smat(self);
	VT_SWITCH( sb->vtype,
		double data=NUM2DBL(other); km_smat_each_d(sb, km_s_pow_d, (void *)&data);,
		COMPLEX data=km_v2c(other); km_smat_each_z(sb, km_s_pow_z, (void *)&data);,
		int data=NUM2INT(other); km_smat_each_i(sb, km_s_pow_i, (void *)&data);,
		bool data=RTEST(other); km_smat_each_b(sb, km_s_pow_b, (void *)&data);,
		km_smat_each_v(sb, km_s_pow_v, (void *)other);
	);
	return self;
}

static void
km_hypot_d(double *y, double *x, void *null)
{
	*y = hypot(*y, *x);
}
static void
km_hypot_z(COMPLEX *y, COMPLEX *x, void *null)
{
	*y = csqrt((*y)*(*y)+(*x)*(*x));
}
static void
km_hypot_v(VALUE *y, VALUE *x, void *null)
{
	*y = rb_funcall(rb_mMath, id_hypot, 2, *y, *x);
}
VALUE
kmm_mat_hypot_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sy = km_mat2smat(self), *sx = km_mat2smat(other);
	if ( sy->vtype != sx->vtype ) {
		rb_raise(km_eVT, "value types must be same");
	}
	CHECK_SAME_SIZE(sy, sx);
	VT_SWITCH( sy->vtype,
		km_smat_each2_d(sy, sx, km_hypot_d, NULL);,
		km_smat_each2_z(sy, sx, km_hypot_z, NULL);,
		rb_raise(km_eVT, "Mat#hypot is not available for int matricies");,
		km_smat_each2_b(sy, sx, km_add_b, NULL);,
		km_smat_each2_v(sy, sx, km_hypot_v, NULL);
	);
	return self;
}

static void
km_s_hypot_d(double *y, void *data)
{
	double *x = (double *)data;
	*y = hypot(*y, *x);
}
static void
km_s_hypot_z(COMPLEX *y, void *data)
{
	COMPLEX *x = (COMPLEX *)data;
	*y = csqrt((*y)*(*y)+(*x)*(*x));
}
static void
km_s_hypot_v(VALUE *y, void *data)
{
	*y = rb_funcall(rb_mMath, id_hypot, 2, *y, (VALUE)data);
}
VALUE
kmm_mat_s_hypot_destl(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *sy = km_mat2smat(self);
	VT_SWITCH( sy->vtype,
		double data=NUM2DBL(other); km_smat_each_d(sy, km_s_hypot_d, (void *)&data);,
		COMPLEX data=km_v2c(other); km_smat_each_z(sy, km_s_hypot_z, (void *)&data);,
		rb_raise(km_eVT, "Mat#hypot is not available for int matricies");,
		bool data=RTEST(other); km_smat_each_b(sy, km_s_add_b, (void *)&data);,
		km_smat_each_v(sy, km_s_hypot_v, (void *)other);
	);
	return self;
}
