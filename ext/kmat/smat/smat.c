#include "../kmat.h"

static void
km_rb_gc_mark_wrap(VALUE *obj, void *null)
{
	rb_gc_mark(*obj);
}
static void
km_smat_mark(void *_data)
{
	SMAT *data = (SMAT *)_data;
	rb_gc_mark(data->parent);
	if ( data->vtype == VT_VALUE && data->body != NULL ) {
		km_smat_each_v(data, km_rb_gc_mark_wrap, NULL);
	}
}

static void
km_smat_free(void *_data)
{
	SMAT *data = (SMAT *)_data;
	if ( data != NULL ) {
		if ( data->stype != ST_SSUB && data->body != NULL ) {
			ruby_xfree(data->body);
		}
		ruby_xfree(data);
	}
}

static size_t
km_smat_size(const void *_data)
{
	const SMAT *data = (const SMAT *)_data;
	if ( data != NULL ) {
		size_t ret = sizeof(SMAT);
		if ( data->body != NULL ) {
			if ( data->stype == ST_FULL ) {
				ret += km_sizeof_vt(data->vtype)*int2size_t(data->m*data->n);
			} else if ( data->stype == ST_RSUB ) {
				ret += sizeof(void *)*int2size_t(data->m*data->n);
			}
		}
		return ret;
	} else {
		return (size_t)0;
	}
}

// functions for allocation
const rb_data_type_t km_mat_data_type = {
	"kmat-Mat",
	{
		km_smat_mark,
		km_smat_free,
		km_smat_size,
		{ (void *)0, (void *)0 }
	},
	(rb_data_type_t *)0,
	(void *)0,
	(VALUE)0
};
VALUE
km_Mat_alloc(VALUE klass)
{
	return TypedData_Wrap_Struct(klass, &km_mat_data_type, km_smat_alloc_with(0, 0, VT_END, NULL));
}

// calloc a side of SMAT
// return->body is the argument `body'
SMAT *
km_smat_alloc_with(int m, int n, VTYPE vt, void *body)
{
	km_check_positive(m, n);
	SMAT *data = ZALLOC(SMAT);
	data->body = body; data->ld = data->m = m; data->n = n;
	data->vtype = vt; data->stype = ST_FULL; data->trans = false;
	data->may_have_sub = false; data->parent = Qnil;
	return data;
}

SMAT *
km_smat_alloc(int m, int n, VTYPE vt)
{
	km_check_positive(m, n);
	void *body = ruby_xcalloc(int2size_t(n*m), km_sizeof_vt(vt));
	return km_smat_alloc_with(m, n, vt, body);
}

// calloc a SMAT of (m, n)-matrix
// return->body will be calloc-ed
void
km_smat_alloc_body(SMAT *data, int m, int n, VTYPE vt)
{
	km_check_positive(m, n);
	if ( data->stype != ST_FULL ) {
		rb_raise(km_eShare, "can't re-alloc submatrix body");
	} else if ( km_smat_have_submatrix_p(data) ) {
		rb_raise(km_eShare, "can't re-alloc supermatrix body");
	}
	data->ld = data->m = m; data->n = n; data->vtype = vt; data->stype = ST_FULL;
	ruby_xfree(data->body);
	data->body = ruby_xcalloc(int2size_t(n*m), km_sizeof_vt(vt));
}
void
km_smat_alloc_pbody(SMAT *data, int m, int n, VTYPE vt)
{
	km_check_positive(m, n);
	if ( data->stype != ST_FULL ) {
		rb_raise(km_eShare, "can't re-alloc submatrix body");
	} else if ( km_smat_have_submatrix_p(data) ) {
		rb_raise(km_eShare, "can't re-alloc supermatrix body");
	}
	data->ld = data->m = m; data->n = n; data->vtype = vt; data->stype = ST_RSUB;
	ruby_xfree(data->body);
	data->body = ruby_xcalloc(int2size_t(n*m), sizeof(void*));
}


// extract and return SMAT from mat_obj
SMAT *
km_mat2smat(VALUE mat_obj)
{
	SMAT *data;
	if ( !rb_obj_is_kind_of(mat_obj, km_cMat) ) {
		rb_raise(rb_eTypeError, "Mat is expected, not %s", rb_obj_classname(mat_obj));
	}
	TypedData_Get_Struct(mat_obj, SMAT, &km_mat_data_type, data);
	return data;
}

// copy the content of SMAT
// if the size are not same, calloc before copying
#define DEFINE_SMAT_COPY(id, type) static void \
km_smat_copy_##id(type *ed, const type *es, void *null) \
{ \
	*ed = *es; \
}
DEFINE_SMAT_COPY(d, double)
DEFINE_SMAT_COPY(z, COMPLEX)
DEFINE_SMAT_COPY(i, int)
DEFINE_SMAT_COPY(b, bool)
DEFINE_SMAT_COPY(v, VALUE)
void
km_smat_copy(SMAT *dest, const SMAT *src)
{
	if ( src->body == NULL ) {
		ruby_xfree(dest->body);
		dest->body = NULL;
	} else {
		if ( dest->vtype != src->vtype || LENGTH(dest) != LENGTH(src) ) { // re-calloc
			if ( dest->stype != ST_FULL ) { // the destination must not be a submatrix
				rb_raise(km_eShare, "can't copy to value-type mismatched or dimension mismatched submatrix");
			}
			ruby_xfree(dest->body);
			dest->body = ruby_xcalloc(LENGTHs(src), km_sizeof_vt(src->vtype));
			dest->ld = src->m; dest->vtype = src->vtype;
		} else if ( dest->m != src->m ) { // need not to resize but reshape is needed
			if ( dest->stype != ST_FULL ) { // the destination must not be a submatrix
				rb_raise(km_eShare, "can't reshape submatrix");
			}
			dest->ld = src->m;
		}
		dest->m = src->m; dest->n = src->n;
		if ( dest->stype==ST_FULL && src->stype==ST_FULL ) {
			dest->ld = src->ld; dest->trans = src->trans;
			memcpy(dest->body, src->body, km_sizeof_vt(src->vtype)*LENGTHs(src));
		} else {
			VT_SWITCH( dest->vtype,
				km_smat_each2_dcd(dest, src, km_smat_copy_d, NULL);,
				km_smat_each2_zcz(dest, src, km_smat_copy_z, NULL);,
				km_smat_each2_ici(dest, src, km_smat_copy_i, NULL);,
				km_smat_each2_bcb(dest, src, km_smat_copy_b, NULL);,
				km_smat_each2_vcv(dest, src, km_smat_copy_v, NULL);
			);
		}
	}
}

// Mat#marhsal_dump/load for Marshal.dump/load
// return an Array with 2 elements
// the first element is a bytes (String) which contains shape and value type
// the second element is `self'.to_ary if the value type is ruby-object
// otherwise, the second element is a bytes (String) which is the same as `self'->body
// the supermatrix-submatrix information will not be saved
VALUE
kmm_mat_marshal_dump(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	if ( smat->stype != ST_FULL ) {
		smat = km_mat2smat(rb_obj_dup(self));
	}
	VALUE headder = rb_str_new((char *)&(smat->vtype), sizeof(VTYPE));
	VALUE body;
	if ( smat->vtype == VT_VALUE ) {
		body = kmm_mat_to_ary(self);
	} else {
		rb_str_cat(headder, (char *)&(smat->m), sizeof(int));
		rb_str_cat(headder, (char *)&(smat->n), sizeof(int));
		rb_str_cat(headder, (char *)&(smat->trans), sizeof(bool));
		body = rb_str_new((char *)(smat->body), (long)km_sizeof_vt(smat->vtype)*LENGTH(smat));
	}
	return rb_ary_new3(2, headder, body);
}
VALUE
kmm_mat_marshal_load(VALUE self, VALUE dump)
{
	SMAT *smat = km_mat2smat(self);
	char *hptr = RSTRING_PTR(rb_ary_entry(dump, 0));
	VALUE body = rb_ary_entry(dump, 1);
	memcpy(&(smat->vtype), hptr, sizeof(VTYPE));
	if ( smat->vtype == VT_VALUE ) {
		km_smat_copy(smat, km_mat2smat(kmm_ary_to_omat(body)));
	} else {
		hptr += sizeof(VTYPE);
		int m, n; bool t;
		memcpy(&m, hptr, sizeof(int)); hptr += sizeof(int);
		memcpy(&n, hptr, sizeof(int)); hptr += sizeof(int);
		memcpy(&t, hptr, sizeof(bool)); hptr += sizeof(bool);
		km_smat_alloc_body(smat, m, n, smat->vtype);
		smat->trans = t;
		size_t sb = km_sizeof_vt(smat->vtype)*LENGTHs(smat);
		if ( RSTRING_LEN(body) != (long)sb ) {
			rb_raise(rb_eArgError, "wrong object given");
		}
		memcpy(smat->body, RSTRING_PTR(body), sb);
	}
	return self;
}

// call `func'(&element, `data') for each elements of `smat'
#define DEFINE_KM_SMAT_EACH_ID(id, type) void \
km_smat_each_##id(SMAT *smat, void (*func)(type *, void *), void *data) \
{ \
	if ( smat->stype == ST_RSUB ) { \
		if ( smat->trans ) { \
			for ( int i=0; i<smat->m; i++ ) { for ( int j=0; j<smat->n; j++ ) { \
				func(smat->id##pbody[j+i*(smat->ld)], data); \
			} } \
		} else { \
			for ( int i=0; i<smat->m; i++ ) { for ( int j=0; j<smat->n; j++ ) { \
				func(smat->id##pbody[i+j*(smat->ld)], data); \
			} } \
		} \
	} else { \
		if ( smat->trans ) { \
			for ( int i=0; i<smat->m; i++ ) { for ( int j=0; j<smat->n; j++ ) { \
				func(&(smat->id##body[j+i*(smat->ld)]), data); \
			} } \
		} else { \
			for ( int i=0; i<smat->m; i++ ) { for ( int j=0; j<smat->n; j++ ) { \
				func(&(smat->id##body[i+j*(smat->ld)]), data); \
			} } \
		} \
	} \
}
DEFINE_KM_SMAT_EACH_ID(d, double)
DEFINE_KM_SMAT_EACH_ID(z, COMPLEX)
DEFINE_KM_SMAT_EACH_ID(i, int)
DEFINE_KM_SMAT_EACH_ID(b, bool)
DEFINE_KM_SMAT_EACH_ID(v, VALUE)

// call `func'(&element, `data', i, j) for each elements of `smat'
#define DEFINE_KM_SMAT_EACH_WI_ID(id, type) void \
km_smat_each_with_index_##id(SMAT *smat, void (*func)(type *, int, int, void *), void *data) \
{ \
	if ( smat->stype == ST_RSUB ) { \
		if ( smat->trans ) { \
			for ( int i=0; i<smat->m; i++ ) { for ( int j=0; j<smat->n; j++ ) { \
				func(smat->id##pbody[j+i*(smat->ld)], i, j, data); \
			} } \
		} else { \
			for ( int i=0; i<smat->m; i++ ) { for ( int j=0; j<smat->n; j++ ) { \
				func(smat->id##pbody[i+j*(smat->ld)], i, j, data); \
			} } \
		} \
	} else { \
		if ( smat->trans ) { \
			for ( int i=0; i<smat->m; i++ ) { for ( int j=0; j<smat->n; j++ ) { \
				func(&(smat->id##body[j+i*(smat->ld)]), i, j, data); \
			} } \
		} else { \
			for ( int i=0; i<smat->m; i++ ) { for ( int j=0; j<smat->n; j++ ) { \
				func(&(smat->id##body[i+j*(smat->ld)]), i, j, data); \
			} } \
		} \
	} \
}
DEFINE_KM_SMAT_EACH_WI_ID(d, double)
DEFINE_KM_SMAT_EACH_WI_ID(z, COMPLEX)
DEFINE_KM_SMAT_EACH_WI_ID(i, int)
DEFINE_KM_SMAT_EACH_WI_ID(b, bool)
DEFINE_KM_SMAT_EACH_WI_ID(v, VALUE)

// call `func'(&(element of sa), &(element of sb), `data') for each elements of `sa'
// you must check the shape of `sb' is the same as that of `sa'
// SEGV will occur if `sb' is smaller than `sa'
#define KM_SMAT_EACH2_BLOOP(elma, id) if ( sb->stype == ST_RSUB ) { \
	if ( sb->trans ) { \
		for ( int i=0; i<sa->m; i++ ) { for ( int j=0; j<sa->n; j++ ) { \
			func(elma, sb->id##pbody[j+i*(sb->ld)], data); \
		} } \
	} else { \
		for ( int i=0; i<sa->m; i++ ) { for ( int j=0; j<sa->n; j++ ) { \
			func(elma, sb->id##pbody[i+j*(sb->ld)], data); \
		} } \
	} \
} else { \
	if ( sb->trans ) { \
		for ( int i=0; i<sa->m; i++ ) { for ( int j=0; j<sa->n; j++ ) { \
			func(elma, &(sb->id##body[j+i*(sb->ld)]), data); \
		} } \
	} else { \
		for ( int i=0; i<sa->m; i++ ) { for ( int j=0; j<sa->n; j++ ) { \
			func(elma, &(sb->id##body[i+j*(sb->ld)]), data); \
		} } \
	} \
}
#define DEFINE_KM_SMAT_EACH2_ID(id, type) void \
km_smat_each2_##id(SMAT *sa, SMAT *sb, void (*func)(type *, type *, void *), void *data) \
{ \
	if ( sa->stype == ST_RSUB ) { \
		if ( sa->trans ) { \
			KM_SMAT_EACH2_BLOOP(sa->id##pbody[j+i*(sa->ld)], id) \
		} else { \
			KM_SMAT_EACH2_BLOOP(sa->id##pbody[i+j*(sa->ld)], id) \
		} \
	} else { \
		if ( sa->trans ) { \
			KM_SMAT_EACH2_BLOOP(&(sa->id##body[j+i*(sa->ld)]), id) \
		} else { \
			KM_SMAT_EACH2_BLOOP(&(sa->id##body[i+j*(sa->ld)]), id) \
		} \
	} \
}
#define DEFINE_KM_SMAT_EACH2_ID_CONST_ID(id, type, id2, type2) void \
km_smat_each2_##id##c##id2(SMAT *sa, const SMAT *sb, void (*func)(type *, const type2 *, void *), void *data) \
{ \
	if ( sa->stype == ST_RSUB ) { \
		if ( sa->trans ) { \
			KM_SMAT_EACH2_BLOOP(sa->id##pbody[j+i*(sa->ld)], id2) \
		} else { \
			KM_SMAT_EACH2_BLOOP(sa->id##pbody[i+j*(sa->ld)], id2) \
		} \
	} else { \
		if ( sa->trans ) { \
			KM_SMAT_EACH2_BLOOP(&(sa->id##body[j+i*(sa->ld)]), id2) \
		} else { \
			KM_SMAT_EACH2_BLOOP(&(sa->id##body[i+j*(sa->ld)]), id2) \
		} \
	} \
}
DEFINE_KM_SMAT_EACH2_ID(d, double)
DEFINE_KM_SMAT_EACH2_ID(z, COMPLEX)
DEFINE_KM_SMAT_EACH2_ID(i, int)
DEFINE_KM_SMAT_EACH2_ID(b, bool)
DEFINE_KM_SMAT_EACH2_ID(v, VALUE)
DEFINE_KM_SMAT_EACH2_ID_CONST_ID(d, double, d, double)
DEFINE_KM_SMAT_EACH2_ID_CONST_ID(z, COMPLEX, z, COMPLEX)
DEFINE_KM_SMAT_EACH2_ID_CONST_ID(i, int, i, int)
DEFINE_KM_SMAT_EACH2_ID_CONST_ID(b, bool, b, bool)
DEFINE_KM_SMAT_EACH2_ID_CONST_ID(v, VALUE, v, VALUE)
DEFINE_KM_SMAT_EACH2_ID_CONST_ID(d, double, b, bool)
DEFINE_KM_SMAT_EACH2_ID_CONST_ID(z, COMPLEX, b, bool)
DEFINE_KM_SMAT_EACH2_ID_CONST_ID(i, int, b, bool)
DEFINE_KM_SMAT_EACH2_ID_CONST_ID(v, VALUE, b, bool)
DEFINE_KM_SMAT_EACH2_ID_CONST_ID(d, double, z, COMPLEX)

// call `func'(&(element of sa), &(element of sb), &(element of sc), `data') for each elements of `sa'
// you must check the shape of `sb' and that of `sc' are the same as that of `sa'
// SEGV will occur if `sb' or `sc' is smaller than `sa'
#define KM_SMAT_EACH3_CLOOP(elma, elmb, id) if ( sc->stype == ST_RSUB ) { \
	if ( sc->trans ) { \
		for ( int i=0; i<sa->m; i++ ) { for ( int j=0; j<sa->n; j++ ) { \
			func(elma, elmb, sc->id##pbody[j+i*(sc->ld)], data); \
		} } \
	} else { \
		for ( int i=0; i<sa->m; i++ ) { for ( int j=0; j<sa->n; j++ ) { \
			func(elma, elmb, sc->id##pbody[i+j*(sc->ld)], data); \
		} } \
	} \
} else { \
	if ( sc->trans ) { \
		for ( int i=0; i<sa->m; i++ ) { for ( int j=0; j<sa->n; j++ ) { \
			func(elma, elmb, &(sc->id##body[j+i*(sc->ld)]), data); \
		} } \
	} else { \
		for ( int i=0; i<sa->m; i++ ) { for ( int j=0; j<sa->n; j++ ) { \
			func(elma, elmb, &(sc->id##body[i+j*(sc->ld)]), data); \
		} } \
	} \
}
#define KM_SMAT_EACH3_BLOOP(elma, id) if ( sb->stype == ST_RSUB ) { \
	if ( sb->trans ) { \
		KM_SMAT_EACH3_CLOOP(elma, sb->id##pbody[j+i*(sb->ld)], id) \
	} else { \
		KM_SMAT_EACH3_CLOOP(elma, sb->id##pbody[i+j*(sb->ld)], id) \
	} \
} else { \
	if ( sb->trans ) { \
		KM_SMAT_EACH3_CLOOP(elma, &(sb->id##body[j+i*(sb->ld)]), id) \
	} else { \
		KM_SMAT_EACH3_CLOOP(elma, &(sb->id##body[i+j*(sb->ld)]), id) \
	} \
}
#define DEFINE_KM_SMAT_EACH3_ID(id, type) void \
km_smat_each3_##id(SMAT *sa, SMAT *sb, SMAT *sc, void (*func)(type *, type *, type *, void *), void *data) \
{ \
	if ( sa->stype == ST_RSUB ) { \
		if ( sa->trans ) { \
			KM_SMAT_EACH3_BLOOP(sa->id##pbody[j+i*(sa->ld)], id) \
		} else { \
			KM_SMAT_EACH3_BLOOP(sa->id##pbody[i+j*(sa->ld)], id) \
		} \
	} else { \
		if ( sa->trans ) { \
			KM_SMAT_EACH3_BLOOP(&(sa->id##body[j+i*(sa->ld)]), id) \
		} else { \
			KM_SMAT_EACH3_BLOOP(&(sa->id##body[i+j*(sa->ld)]), id) \
		} \
	} \
}
// DEFINE_KM_SMAT_EACH3_ID(d, double)
// DEFINE_KM_SMAT_EACH3_ID(z, COMPLEX)
// DEFINE_KM_SMAT_EACH3_ID(i, int)
// DEFINE_KM_SMAT_EACH3_ID(b, bool)
// DEFINE_KM_SMAT_EACH3_ID(v, VALUE)
void
km_smat_each3_zcdcd(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(COMPLEX *, const double *, const double *, void *), void *data)
{
	if ( sa->stype == ST_RSUB ) {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(sa->zpbody[j+i*(sa->ld)], d)
		} else {
			KM_SMAT_EACH3_BLOOP(sa->zpbody[i+j*(sa->ld)], d)
		}
	} else {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(&(sa->zbody[j+i*(sa->ld)]), d)
		} else {
			KM_SMAT_EACH3_BLOOP(&(sa->zbody[i+j*(sa->ld)]), d)
		}
	}
}

void
km_smat_each3_bcdcd(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(bool *, const double *, const double *, void *), void *data)
{
	if ( sa->stype == ST_RSUB ) {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(sa->bpbody[j+i*(sa->ld)], d)
		} else {
			KM_SMAT_EACH3_BLOOP(sa->bpbody[i+j*(sa->ld)], d)
		}
	} else {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(&(sa->bbody[j+i*(sa->ld)]), d)
		} else {
			KM_SMAT_EACH3_BLOOP(&(sa->bbody[i+j*(sa->ld)]), d)
		}
	}
}
void
km_smat_each3_bczcz(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(bool *, const COMPLEX *, const COMPLEX *, void *), void *data)
{
	if ( sa->stype == ST_RSUB ) {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(sa->bpbody[j+i*(sa->ld)], z)
		} else {
			KM_SMAT_EACH3_BLOOP(sa->bpbody[i+j*(sa->ld)], z)
		}
	} else {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(&(sa->bbody[j+i*(sa->ld)]), z)
		} else {
			KM_SMAT_EACH3_BLOOP(&(sa->bbody[i+j*(sa->ld)]), z)
		}
	}
}
void
km_smat_each3_bcici(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(bool *, const int *, const int *, void *), void *data)
{
	if ( sa->stype == ST_RSUB ) {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(sa->bpbody[j+i*(sa->ld)], i)
		} else {
			KM_SMAT_EACH3_BLOOP(sa->bpbody[i+j*(sa->ld)], i)
		}
	} else {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(&(sa->bbody[j+i*(sa->ld)]), i)
		} else {
			KM_SMAT_EACH3_BLOOP(&(sa->bbody[i+j*(sa->ld)]), i)
		}
	}
}
void
km_smat_each3_bcbcb(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(bool *, const bool *, const bool *, void *), void *data)
{
	if ( sa->stype == ST_RSUB ) {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(sa->bpbody[j+i*(sa->ld)], b)
		} else {
			KM_SMAT_EACH3_BLOOP(sa->bpbody[i+j*(sa->ld)], b)
		}
	} else {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(&(sa->bbody[j+i*(sa->ld)]), b)
		} else {
			KM_SMAT_EACH3_BLOOP(&(sa->bbody[i+j*(sa->ld)]), b)
		}
	}
}
void
km_smat_each3_bcvcv(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(bool *, const VALUE *, const VALUE *, void *), void *data)
{
	if ( sa->stype == ST_RSUB ) {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(sa->bpbody[j+i*(sa->ld)], v)
		} else {
			KM_SMAT_EACH3_BLOOP(sa->bpbody[i+j*(sa->ld)], v)
		}
	} else {
		if ( sa->trans ) {
			KM_SMAT_EACH3_BLOOP(&(sa->bbody[j+i*(sa->ld)]), v)
		} else {
			KM_SMAT_EACH3_BLOOP(&(sa->bbody[i+j*(sa->ld)]), v)
		}
	}
}

