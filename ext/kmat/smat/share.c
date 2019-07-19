#include "../kmat.h"

// freeze `self' and all submatricies of `self'
// `self' must be a full-matrix (can't deep_freeze submatrix)
static void
km_mat_deep_freeze_func(VALUE *ent, void *null)
{
	if ( rb_respond_to(*ent, id_deep_freeze) ) {
		rb_funcall(*ent, id_deep_freeze, 0);
	} else {
		rb_obj_freeze(*ent);
	}
}
static void
km_mat_deep_freeze0(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype == VT_VALUE ) {
		km_smat_each_v(smat, km_mat_deep_freeze_func, NULL);
	}
}
VALUE
kmm_mat_deep_freeze(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->stype != ST_FULL ) {
		rb_raise(km_eShare, "can't deep_freeze submatricies directory. deep_freeze supermatrix instead");
	}
	if ( smat->may_have_sub ) {
		VALUE subs = kmm_mat_submatricies(self);
		for ( long i=0; i<RARRAY_LEN(subs); i++ ) {
			km_mat_deep_freeze0(rb_ary_entry(subs, i));
		}
	}
	km_mat_deep_freeze0(self);
	return self;
}

typedef int each_obj_callback(void *, void *, size_t, void *);
void rb_objspace_each_objects(each_obj_callback *callback, void *data);
int rb_objspace_internal_object_p(VALUE obj);

// return true if `self' has at least one submatrix, otherwise, return false
struct km_have_sub {
	bool found;
	VALUE target;
};
static int
km_have_sub_p(void *vstart, void *vend, size_t stride, void *data)
{
	struct km_have_sub *have_sub = (struct km_have_sub *)data;
	for (VALUE v = (VALUE)vstart; v != (VALUE)vend; v += stride) {
		if (!rb_objspace_internal_object_p(v)) {
			if (rb_obj_is_kind_of(v, km_cMat)) {
				SMAT *smat = km_mat2smat(v);
				if ( have_sub->target == smat->parent ) {
					have_sub->found = true;
					return 1;
				}
			}
		}
	}
	return 0;
}
VALUE
kmm_mat_have_submatrix_p(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->may_have_sub ) {
		struct km_have_sub data = {false, self};
		rb_objspace_each_objects(km_have_sub_p, (void *)&data);
		if ( data.found ) {
			return Qtrue;
		} else {
			smat->may_have_sub = false;
			return Qfalse;
		}
	} else {
		return Qfalse;
	}
}

// the SMAT version of above
bool
km_smat_have_submatrix_p(SMAT *smat)
{
	if ( smat->may_have_sub ) {
		VALUE v = km_smat_find_value(smat);
		if ( v == Qnil ) {
			rb_raise(km_eInternal, "km_smat_have_submatrix_p has been called with out-of-ruby-controlled Mat struct");
		} else {
			return (kmm_mat_have_submatrix_p(v) != Qfalse);
		}
	} else {
		return false;
	}
}

// find VALUE of `smat' from ObjectSpace.each_object
// return Qnil if not found
struct km_find_value {
	SMAT *target;
	VALUE out;
};
static int
km_find_value(void *vstart, void *vend, size_t stride, void *data)
{
	struct km_find_value *fv = (struct km_find_value *)data;
	for (VALUE v = (VALUE)vstart; v != (VALUE)vend; v += stride) {
		if (!rb_objspace_internal_object_p(v)) {
			if ( rb_obj_is_kind_of(v, km_cMat) && km_mat2smat(v) == fv->target ) {
				fv->out = v;
				return 1;
			}
		}
	}
	return 0;
}
VALUE
km_smat_find_value(SMAT *smat)
{
	struct km_find_value data = {smat, Qnil};
	rb_objspace_each_objects(km_find_value, (void *)&data);
	return data.out;
}

// return an Array submatricies of `self'
// the return will be nil if `self' is a submatrix (submatrix cannot have a submatrix of itself)
// the return will be an empty Array ([]) if `self' is full-matrix but does not have any submatricies
struct km_submx {
	VALUE parent;
	VALUE list;
};
static int
km_mat_submx(void *vstart, void *vend, size_t stride, void *data)
{
	struct km_submx *sm = (struct km_submx *)data;
	for ( VALUE v = (VALUE)vstart; v != (VALUE)vend; v += stride ) {
		if (!rb_objspace_internal_object_p(v)) {
			if ( rb_obj_is_kind_of(v, km_cMat) && km_mat2smat(v)->parent == sm->parent ) {
				rb_ary_push(sm->list, v);
			}
		}
	}
	return 0;
}
VALUE
kmm_mat_submatricies(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->stype != ST_FULL ) {
		return Qnil;
	} else if ( !smat->may_have_sub ) {
		return rb_ary_new();
	}
	struct km_submx data = {self, rb_ary_new()};
	rb_objspace_each_objects(km_mat_submx, (void *)&data);
	if ( RARRAY_LEN(data.list) == 0 ) {
		smat->may_have_sub = false;
	}
	return data.list;
}

// return supermatrix of `self'
// the return will be nil if `self' is full-matrix
VALUE
kmm_mat_supermatrix(VALUE self)
{
	SMAT *data = km_mat2smat(self);
	return data->parent;
}

// if `self' is a supermatrix of some submatricies, detach all the submatricies by replacing them by copied full-matrix
// if `self' is a submatrix of a supermatrix, replace `self' by a copied full-matrix of 'self'
// return self if at least a matrix is replaced, otherwise, return nil
VALUE
kmm_mat_detach(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 0, 1);
	bool check_other = ( ( argc == 1 ) && RTEST(argv[0]) );
	SMAT *smat = km_mat2smat(self);
	if ( smat->stype == ST_FULL ) {
		VALUE list = kmm_mat_submatricies(self);
		if ( !smat->may_have_sub ) {
			return Qnil;
		}
		for (long i=0; i<RARRAY_LEN(list); i++) {
			VALUE elm = rb_ary_entry(list, i);
			kmm_mat_replace(elm, rb_obj_dup(elm));
		}
		smat->may_have_sub = false;
		return self;
	} else {
		VALUE old_parent = smat->parent;
		kmm_mat_replace(self, rb_obj_dup(self));
		if ( check_other ) {
			kmm_mat_submatricies(old_parent);
		}
		return self;
	}
}

// if `self' is a supermatrix of some submatricies, detach all the submatricies by replacing them by (0, 0)-matrix
// if `self' is a submatrix of a supermatrix, detach from the supermatrix
// in any case, `self' is replaced by (0, 0)-matrix, and return `self'
static VALUE
km_kill_func(VALUE elm, VALUE nil, VALUE self)
{
	SMAT *ssub = km_mat2smat(elm);
	if ( ssub->stype != ST_SSUB ) {
		ruby_xfree(ssub->body);
	}
	ssub->stype = ST_FULL;
	ssub->body = NULL;
	ssub->m = 0; ssub->n = 0; ssub->ld = 0;
	ssub->parent = Qnil;
	return elm;
}
VALUE
kmm_mat__kill(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( smat->stype == ST_FULL ) {
		VALUE list = kmm_mat_submatricies(self);
		rb_iterate(rb_ary_each, list, km_kill_func, Qnil);
	}
	if ( smat->stype != ST_SSUB ) {
		ruby_xfree(smat->body);
	}
	smat->stype = ST_FULL;
	smat->body = NULL;
	smat->m = 0; smat->n = 0; smat->ld = 0;
	smat->may_have_sub = false; smat->parent = Qnil;
	return self;
}

// replace the substance of `self' by `val'
// if `val' is a submatrix, `self' is a submatrix of the same supermatrix
// if `val' is a full-matrix, `self' is a copy of it (submatricies are not succeeded)
VALUE
kmm_mat_replace(VALUE self, VALUE val)
{
	km_check_frozen(self);
	SMAT *dest = km_mat2smat(self);
	SMAT *src = km_mat2smat(val);
	if ( kmm_mat_have_submatrix_p(self) ) {
		rb_raise(km_eShare, "can't replace supermatrix. try detach before replacing");
	}
	if ( src->stype == ST_FULL ) {
		switch ( dest->stype ) {
		case ST_RSUB:
			ruby_xfree(dest->body);
			// FALL THROUGH
		case ST_SSUB:
			dest->vtype = VT_END;
			dest->stype = ST_FULL;
			dest->parent = Qnil;
			dest->may_have_sub = false;
			dest->body = NULL;
			break;
		default:
			break;
		}
		km_smat_copy(dest, src);
	} else if ( src->stype == ST_SSUB ) {
		if ( dest->stype != ST_SSUB ) {
			ruby_xfree(dest->body);
		}
		dest->body = src->body;
		dest->ld = src->ld; dest->m = src->m; dest->n = src->n;
		dest->vtype = src->vtype; dest->stype = ST_SSUB;
		dest->trans = src->trans;
		dest->parent = src->parent;
	} else if ( src->stype == ST_RSUB ) {
		if ( dest->stype != ST_RSUB || LENGTH(dest) != LENGTH(src) ) {
			if ( dest->stype != ST_SSUB ) {
				ruby_xfree(dest->body);
			}
			dest->body = ruby_xcalloc(LENGTHs(src), sizeof(void*));
		}
		memcpy(dest->body, src->body, sizeof(void*)*LENGTHs(src));
		dest->ld = src->ld; dest->m = src->m; dest->n = src->n;
		dest->vtype = src->vtype; dest->stype = ST_RSUB;
		dest->trans = src->trans;
		dest->parent = src->parent;
	} else {
		rb_raise(km_eInternal, "unknown storage type");
	}
	return self;
}

VALUE
km_Mat_ssub(int i, int j, int m, int n, VALUE super)
{
	SMAT *ssup = km_mat2smat(super);
	if ( ssup->m < i+m || ssup->n < j+n ) {
		rb_raise(rb_eIndexError, "given index+size (%d+%d, %d+%d) is out of range (%d, %d)", i, m, j, n, ssup->m, ssup->n);
	}
	SMAT *sret;
	if ( ssup->stype == ST_FULL ) {
		VT_SWITCH( ssup->vtype,
			sret = km_smat_alloc_with(m, n, VT_DOUBLE, ssup->dbody+INDEX(ssup, i, j));,
			sret = km_smat_alloc_with(m, n, VT_COMPLEX, ssup->zbody+INDEX(ssup, i, j));,
			sret = km_smat_alloc_with(m, n, VT_INT, ssup->ibody+INDEX(ssup, i, j));,
			sret = km_smat_alloc_with(m, n, VT_BOOL, ssup->bbody+INDEX(ssup, i, j));,
			sret = km_smat_alloc_with(m, n, VT_VALUE, ssup->vbody+INDEX(ssup, i, j));
		);
		sret->parent = super;
		ssup->may_have_sub = true;
	} else if ( ssup->stype == ST_SSUB ) {
		VT_SWITCH( ssup->vtype,
			sret = km_smat_alloc_with(m, n, VT_DOUBLE, ssup->dbody+INDEX(ssup, i, j));,
			sret = km_smat_alloc_with(m, n, VT_COMPLEX, ssup->zbody+INDEX(ssup, i, j));,
			sret = km_smat_alloc_with(m, n, VT_INT, ssup->ibody+INDEX(ssup, i, j));,
			sret = km_smat_alloc_with(m, n, VT_BOOL, ssup->bbody+INDEX(ssup, i, j));,
			sret = km_smat_alloc_with(m, n, VT_VALUE, ssup->vbody+INDEX(ssup, i, j));
		);
		sret->parent = ssup->parent;
	} else if ( ssup->stype == ST_RSUB ) {
		int is[m]; int js[n];
		for ( int k=0; k<m; k++ ) { is[k] = i+k; }
		for ( int k=0; k<n; k++ ) { js[k] = j+k; }
		return km_Mat_rsub1(m, n, is, js, super);
	} else {
		rb_raise(km_eInternal, "unknown storage type");
	}
	sret->trans=ssup->trans; sret->ld = ssup->ld; sret->stype = ST_SSUB;
	VALUE ret = TypedData_Wrap_Struct(km_cMat, &km_mat_data_type, sret);
	km_infect_frozen(super, ret);
	return ret;
}

static void
km_rsub_check_range(SMAT *ssup, int ii, int jj)
{
	if ( ii < 0 || ssup->m <= ii || jj < 0 || ssup->n <= jj ) {
		rb_raise(rb_eIndexError, "given index (%d, %d) is out of range (%d, %d)", ii, jj, ssup->m, ssup->n);
	}
}
#define RSUB1_LOOPr(id) for ( int i=0; i<m; i++ ) { for ( int j=0; j<n; j++ ) { \
	km_rsub_check_range(ssup, is[i], js[j]); \
	sret->id##pbody[i+j*m] = ssup->id##pbody[INDEX(ssup, is[i], js[j])]; \
} }
#define RSUB1_LOOP(id) for ( int i=0; i<m; i++ ) { for ( int j=0; j<n; j++ ) { \
	km_rsub_check_range(ssup, is[i], js[j]); \
	sret->id##pbody[i+j*m] = ssup->id##body + INDEX(ssup, is[i], js[j]); \
} }
VALUE
km_Mat_rsub1(int m, int n, int *is, int *js, VALUE super)
{
	km_check_positive(m, n);
	VALUE ret = km_Mat_alloc(km_cMat);
	SMAT *ssup = km_mat2smat(super);
	SMAT *sret = km_mat2smat(ret);
	km_smat_alloc_pbody(sret, m, n, ssup->vtype);
	if ( ssup->stype == ST_RSUB ) {
		VT_SWITCH( ssup->vtype,
			RSUB1_LOOPr(d);,
			RSUB1_LOOPr(z);,
			RSUB1_LOOPr(i);,
			RSUB1_LOOPr(b);,
			RSUB1_LOOPr(v);
		);
		sret->parent = ssup->parent;
	} else {
		VT_SWITCH( ssup->vtype,
			RSUB1_LOOP(d);,
			RSUB1_LOOP(z);,
			RSUB1_LOOP(i);,
			RSUB1_LOOP(b);,
			RSUB1_LOOP(v);
		);
		if ( ssup->stype == ST_FULL ) {
			sret->parent = super;
			ssup->may_have_sub = true;
		} else {
			sret->parent = ssup->parent;
		}
	}
	km_infect_frozen(super, ret);
	return ret;
}

#define RSUB2_LOOPr(id) for( int i=0; i<m; i++ ) { for ( int j=0; j<n; j++ ) { \
	km_rsub_check_range(ssup, is[i+j*m], js[i+j*m]); \
	sret->id##pbody[i+j*m] = ssup->id##pbody[INDEX(ssup, is[i+j*m], js[i+j*m])]; \
} }
#define RSUB2_LOOP(id) for( int i=0; i<m; i++ ) { for ( int j=0; j<n; j++ ) { \
	km_rsub_check_range(ssup, is[i+j*m], js[i+j*m]); \
	sret->id##pbody[i+j*m] = ssup->id##body+INDEX(ssup, is[i+j*m], js[i+j*m]); \
} }
VALUE
km_Mat_rsub2(int m, int n, int *is, int *js, VALUE super)
{
	km_check_positive(m, n);
	VALUE ret = km_Mat_alloc(km_cMat);
	SMAT *ssup = km_mat2smat(super);
	SMAT *sret = km_mat2smat(ret);
	km_smat_alloc_pbody(sret, m, n, ssup->vtype);
	if ( ssup->stype == ST_RSUB ) {
		VT_SWITCH( ssup->vtype,
			RSUB2_LOOPr(d);,
			RSUB2_LOOPr(z);,
			RSUB2_LOOPr(i);,
			RSUB2_LOOPr(b);,
			RSUB2_LOOPr(v);
		);
		sret->parent = ssup->parent;
	} else {
		VT_SWITCH( ssup->vtype,
			RSUB2_LOOP(d);,
			RSUB2_LOOP(z);,
			RSUB2_LOOP(i);,
			RSUB2_LOOP(b);,
			RSUB2_LOOP(v);
		);
		if ( ssup->stype == ST_FULL ) {
			sret->parent = super;
			ssup->may_have_sub = true;
		} else {
			sret->parent = ssup->parent;
		}
	}
	km_infect_frozen(super, ret);
	return ret;
}
