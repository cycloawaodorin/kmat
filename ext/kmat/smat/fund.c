#include "../kmat.h"

// make a (m, n)-matrix with value type vt
VALUE
km_Mat(int m, int n, VTYPE vt)
{
	VALUE ret = km_Mat_alloc(km_cMat);
	SMAT *smat = km_mat2smat(ret);
	km_smat_alloc_body(smat, m, n, vt);
	return ret;
}

// Mat#initialize
// the arguments specifies the shape and the value type
// if the value type is ruby-object, it is filled by nil
// otherwise, the return is filled by 0
// if a block given, the (i, j)-th element is a return of yield(i, j)
#define DEFINE_INIT_LOOP_FUNC(id, type, func) static void \
km_initialize_loop_func_##id(type *elm, int i, int j, void *null) \
{\
	*elm = func( rb_yield( rb_ary_new3(2, INT2NUM(i), INT2NUM(j)) ) ); \
}
DEFINE_INIT_LOOP_FUNC(d, double, NUM2DBL)
DEFINE_INIT_LOOP_FUNC(z, COMPLEX, km_v2c)
DEFINE_INIT_LOOP_FUNC(i, int, NUM2INT)
DEFINE_INIT_LOOP_FUNC(b, bool, RTEST)
DEFINE_INIT_LOOP_FUNC(v, VALUE, ITSELF)
struct km_initialize_loop_arg {
	SMAT *smat;
	VALUE *argv;
	VTYPE vt;
};
static VALUE
km_initialize_loop(VALUE varg)
{
	struct km_initialize_loop_arg *arg = (struct km_initialize_loop_arg *)varg;
	km_smat_alloc_body(arg->smat, NUM2INT(arg->argv[0]), NUM2INT(arg->argv[1]), arg->vt);
	if ( rb_block_given_p() ) {
		VT_SWITCH( arg->vt,
			km_smat_each_with_index_d(arg->smat, km_initialize_loop_func_d, NULL);,
			km_smat_each_with_index_z(arg->smat, km_initialize_loop_func_z, NULL);,
			km_smat_each_with_index_i(arg->smat, km_initialize_loop_func_i, NULL);,
			km_smat_each_with_index_b(arg->smat, km_initialize_loop_func_b, NULL);,
			km_smat_each_with_index_v(arg->smat, km_initialize_loop_func_v, NULL);
		);
	} else {
		if ( arg->vt == VT_VALUE ) {
			for ( int i=0; i<arg->smat->m; i++ ) {
				for ( int j=0; j<arg->smat->n; j++ ) {
					arg->smat->vbody[i+j*arg->smat->ld] = Qnil;
				}
			}
		}
	}
	return Qnil;
}
VALUE
kmm_mat_initialize(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 2, 3);
	SMAT *smat = km_mat2smat(self); VTYPE vt;
	if ( smat->stype != ST_FULL ) {
		rb_raise(km_eShare, "can't re-initialize submatrix. try detach before re-initializing");
	} else if ( kmm_mat_have_submatrix_p(self) ) {
		rb_raise(km_eShare, "can't re-initialize supermatrix. try detach before re-initializing");
	}
	if ( argc == 2 ) {
		vt = VT_DOUBLE;
	} else {
		vt = km_sym2vt(argv[2]);
	}
	struct km_initialize_loop_arg arg;
	arg.smat = smat; arg.argv = argv; arg.vt = vt;
	km_gc_escape(km_initialize_loop, (VALUE)&arg);
	return self;
}

// Mat#initialize_copy. this will used by Mat#dup or Mat#clone
// copy shape, value type and the elements from `obj'
// `self' will be a full-matrix even if `obj' is a submatrix
VALUE
kmm_mat_initialize_copy(VALUE self, VALUE obj)
{
	rb_obj_init_copy(self, obj);
	SMAT *ss = km_mat2smat(obj);
	SMAT *sd = km_mat2smat(self);
	if ( sd->stype != ST_FULL ) {
		rb_raise(km_eShare, "can't re-initialize submatrix. try detach before re-initializing");
	} else if ( kmm_mat_have_submatrix_p(self) ) {
		rb_raise(km_eShare, "can't re-initialize supermatrix. try detach before re-initializing");
	}
	sd->vtype = VT_END;
	ruby_xfree(sd->body);
	km_smat_copy(sd, ss);
	return self;
}

// reshpae to (`vm', `vn')-matrix `vm'*`vn' must be the same as `self'.row_size*'self'.col_size
VALUE
kmm_mat_reshape(VALUE self, VALUE vm, VALUE vn)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( smat->stype == ST_SSUB ) {
		rb_raise(km_eShare, "can't reshape submatrix. try detach before reshaping");
	}
	int m = NUM2INT(vm), n = NUM2INT(vn);
	km_check_positive(m, n);
	if ( m*n == LENGTH(smat) ) {
		smat->trans = false;
		smat->m = m; smat->n = n;
		smat->ld = m;
		return self;
	} else {
		rb_raise(km_eDim, "the length of the matrix must no be changed, (%d, %d) -> (%d, %d) is unavailable",
			smat->m, smat->n, m, n);
	}
}

// change shape to (`vm', `vn')
// calloc() will by called if needed
VALUE
kmm_mat_resize(VALUE self, VALUE vm, VALUE vn)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( smat->stype != ST_FULL ) {
		rb_raise(km_eShare, "can't resize submatrix. try detach before reshaping");
	} else if ( kmm_mat_have_submatrix_p(self) ) {
		rb_raise(km_eShare, "can't resize supermatrix. try detach before reshaping");
	}
	int m = NUM2INT(vm), n = NUM2INT(vn);
	if ( m*n == LENGTH(smat) ) {
		return kmm_mat_reshape(self, vm, vn);
	} else {
		VALUE argv[3] = {vm, vn, kmm_mat_vtype(self)};
		return kmm_mat_initialize(3, argv, self);
	}
}

// transpose `self'
// alias t
VALUE
kmm_mat_transpose_dest(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	smat->trans = !(smat->trans);
	SWAP(int, smat->m, smat->n);
	return self;
}

// fill elements of `self' by `val'
#define DEFINE_FILL_FUNC(id, type) static void \
km_fill_func_##id(type *elm, void *data) \
{\
	*elm = *(type *)data; \
}
DEFINE_FILL_FUNC(d, double)
DEFINE_FILL_FUNC(z, COMPLEX)
DEFINE_FILL_FUNC(i, int)
DEFINE_FILL_FUNC(b, bool)
DEFINE_FILL_FUNC(v, VALUE)
VALUE
kmm_mat_fill(VALUE self, VALUE val)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	VT_SWITCH( smat->vtype,
		double data=NUM2DBL(val); km_smat_each_d(smat, km_fill_func_d, (void *)&data);,
		COMPLEX data=km_v2c(val); km_smat_each_z(smat, km_fill_func_z, (void *)&data);,
		int data=NUM2INT(val); km_smat_each_i(smat, km_fill_func_i, (void *)&data);,
		bool data=RTEST(val); km_smat_each_b(smat, km_fill_func_b, (void *)&data);,
		VALUE data=val; km_smat_each_v(smat, km_fill_func_v, (void *)&data);
	);
	return self;
}

// fill elements of `self' by 0
// if `self' is a boolean matrix, the filler is false
// if `self' is a ruby-object matrix, the filler is 0, an instance of Integer
VALUE
kmm_mat_zero(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	VT_SWITCH( smat->vtype,
		double data=0.0; km_smat_each_d(smat, km_fill_func_d, (void *)&data);,
		COMPLEX data=cpack(0.0, 0.0); km_smat_each_z(smat, km_fill_func_z, (void *)&data);,
		int data=0; km_smat_each_i(smat, km_fill_func_i, (void *)&data);,
		bool data=false; km_smat_each_b(smat, km_fill_func_b, (void *)&data);,
		VALUE data=INT2NUM(0); km_smat_each_v(smat, km_fill_func_v, (void *)&data);
	);
	return self;
}

// fill elements of `self' by eye()
#define DEFINE_EYE_FUNC(id, type) struct km_eye_##id { \
	type zero; \
	type one; \
}; \
static void \
km_eye_func_##id(type *elm, int i, int j, void *data_) \
{\
	struct km_eye_##id *data = (struct km_eye_##id *)data_; \
	if ( i == j ) { \
		*elm = data->one; \
	} else { \
		*elm = data->zero; \
	} \
}
DEFINE_EYE_FUNC(d, double)
DEFINE_EYE_FUNC(z, COMPLEX)
DEFINE_EYE_FUNC(i, int)
DEFINE_EYE_FUNC(b, bool)
DEFINE_EYE_FUNC(v, VALUE)
VALUE
kmm_mat_eye(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype == VT_DOUBLE ) {
		struct km_eye_d data = {0.0, 1.0};
		km_smat_each_with_index_d(smat, km_eye_func_d, (void *)&data);
	} else if ( smat->vtype == VT_COMPLEX ) {
		struct km_eye_z data = {cpack(0.0, 0.0), cpack(1.0, 0.0)};
		km_smat_each_with_index_z(smat, km_eye_func_z, (void *)&data);
	} else if ( smat->vtype == VT_INT ) {
		struct km_eye_i data = {0, 1};
		km_smat_each_with_index_i(smat, km_eye_func_i, (void *)&data);
	} else if ( smat->vtype == VT_BOOL ) {
		struct km_eye_b data = {false, true};
		km_smat_each_with_index_b(smat, km_eye_func_b, (void *)&data);
	} else if ( smat->vtype == VT_VALUE ) {
		struct km_eye_v data = {INT2NUM(0), INT2NUM(1)};
		km_smat_each_with_index_v(smat, km_eye_func_v, (void *)&data);
	} else {
		rb_raise(km_eInternal, "unknown value type");
	}
	return self;
}

// fill the elements of `self' by kmat_vt_rand()
// float, object: U(0, 1), U(0...max), U(range)
// int: U({0, 1}), U({0...max}), U({range}), "1 with probability arg and 0 with probability 1-arg"
// bool: U({false, true}), "true with probability arg and false with probability 1-arg"
struct km_rand_arg {
	VALUE random;
	VALUE arg;
};
#define DEFINE_RAND_FUNC(id, type, cast_func) static void \
km_rand_func_##id(type *elm, void *data_) \
{\
	struct km_rand_arg *data = (struct km_rand_arg *)data_; \
	*elm = cast_func(rb_funcall(data->random, id__kmat_rand_##id, 1, data->arg)); \
}
DEFINE_RAND_FUNC(d, double, NUM2DBL)
DEFINE_RAND_FUNC(z, COMPLEX, km_v2c)
DEFINE_RAND_FUNC(i, int, NUM2INT)
DEFINE_RAND_FUNC(b, bool, RTEST)
DEFINE_RAND_FUNC(v, VALUE, ITSELF)
VALUE
kmm_mat__rand0(VALUE self, VALUE random, VALUE arg)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	struct km_rand_arg data = {random, arg};
	VT_SWITCH( smat->vtype,
		km_smat_each_d(smat, km_rand_func_d, (void *)&data);,
		km_smat_each_z(smat, km_rand_func_z, (void *)&data);,
		km_smat_each_i(smat, km_rand_func_i, (void *)&data);,
		km_smat_each_b(smat, km_rand_func_b, (void *)&data);,
		km_smat_each_v(smat, km_rand_func_v, (void *)&data);
	);
	return self;
}

// fill the elments of `self' by random values with follows N(0, 1)
static void
km_randn_func(double *elm, void *data)
{
	*elm = km_rand_normal(*(VALUE *)data);
}
VALUE
kmm_mat__randn0(VALUE self, VALUE random)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype != VT_DOUBLE ) {
		rb_raise(km_eVT, "Mat#randn is available for float matricies, change value-type before call this");
	}
	if ( smat->stype == ST_FULL ) {
		km_fill_normal(LENGTH(smat), smat->dbody, random);
	} else if ( smat->stype == ST_SSUB ) {
		if ( smat->trans ) {
			for ( int i=0; i<smat->m; i++ ) {
				km_fill_normal(smat->n, smat->dbody+i*smat->ld, random);
			}
		} else {
			for ( int i=0; i<smat->n; i++ ) {
				km_fill_normal(smat->m, smat->dbody+i*smat->ld, random);
			}
		}
	} else if ( smat->stype == ST_RSUB ) {
		VALUE data = random;
		km_smat_each_d(smat, km_randn_func, (void *)&data);
	} else {
		rb_raise(km_eInternal, "unknown storage type");
	}
	return self;
}

// invoke yield(element) for all the elements
#define DEFINE_EACH_FUNC(id, type, func) static void \
km_each_func_##id(type *elm, void *null) \
{ \
	rb_yield(func(*elm)); \
}
DEFINE_EACH_FUNC(d, double, rb_float_new)
DEFINE_EACH_FUNC(z, COMPLEX, km_c2v)
DEFINE_EACH_FUNC(i, int, INT2NUM)
DEFINE_EACH_FUNC(b, bool, TF2V)
DEFINE_EACH_FUNC(v, VALUE, ITSELF)
VALUE
kmm_mat_each(VALUE self)
{
	if ( rb_block_given_p() ) {
		SMAT *smat = km_mat2smat(self);
		VT_SWITCH( smat->vtype,
			km_smat_each_d(smat, km_each_func_d, NULL);,
			km_smat_each_z(smat, km_each_func_z, NULL);,
			km_smat_each_i(smat, km_each_func_i, NULL);,
			km_smat_each_b(smat, km_each_func_b, NULL);,
			km_smat_each_v(smat, km_each_func_v, NULL);
		);
		return self;
	} else {
		return rb_enumeratorize(self, sym_each, 0, NULL);
	}
}

// invoke yield(element, i, j) for all the elements
#define DEFINE_EACH_WI2_FUNC(id, type, func) static void \
km_each_with_index2_func_##id(type *elm, int i, int j, void *null) \
{ \
	rb_yield(rb_ary_new3(3, func(*elm), INT2NUM(i), INT2NUM(j))); \
}
DEFINE_EACH_WI2_FUNC(d, double, rb_float_new)
DEFINE_EACH_WI2_FUNC(z, COMPLEX, km_c2v)
DEFINE_EACH_WI2_FUNC(i, int, INT2NUM)
DEFINE_EACH_WI2_FUNC(b, bool, TF2V)
DEFINE_EACH_WI2_FUNC(v, VALUE, ITSELF)
// alias each_with_index
VALUE
kmm_mat_each_with_index2(VALUE self)
{
	if ( rb_block_given_p() ) {
		SMAT *smat = km_mat2smat(self);
		VT_SWITCH( smat->vtype,
			km_smat_each_with_index_d(smat, km_each_with_index2_func_d, NULL);,
			km_smat_each_with_index_z(smat, km_each_with_index2_func_z, NULL);,
			km_smat_each_with_index_i(smat, km_each_with_index2_func_i, NULL);,
			km_smat_each_with_index_b(smat, km_each_with_index2_func_b, NULL);,
			km_smat_each_with_index_v(smat, km_each_with_index2_func_v, NULL);
		);
		return self;
	} else {
		return rb_enumeratorize(self, sym_each_with_index2, 0, NULL);
	}
}

// invoke yield(element) and replace the element by the return for all the elements
#define DEFINE_MMAP_FUNC(id, type, func, func2) static void \
km_mmap_func_##id(type *elm, void *null) \
{ \
	*elm = func2(rb_yield(func(*elm))); \
}
DEFINE_MMAP_FUNC(d, double, rb_float_new, NUM2DBL)
DEFINE_MMAP_FUNC(z, COMPLEX, km_c2v, km_v2c)
DEFINE_MMAP_FUNC(i, int, INT2NUM, NUM2INT)
DEFINE_MMAP_FUNC(b, bool, TF2V, RTEST)
DEFINE_MMAP_FUNC(v, VALUE, ITSELF, ITSELF)
// alias map
VALUE
kmm_mat_mmap_dest(VALUE self)
{
	if ( rb_block_given_p() ) {
		SMAT *smat = km_mat2smat(self);
		VT_SWITCH( smat->vtype,
			km_smat_each_d(smat, km_mmap_func_d, NULL);,
			km_smat_each_z(smat, km_mmap_func_z, NULL);,
			km_smat_each_i(smat, km_mmap_func_i, NULL);,
			km_smat_each_b(smat, km_mmap_func_b, NULL);,
			km_smat_each_v(smat, km_mmap_func_v, NULL);
		);
		return self;
	} else {
		return rb_enumeratorize(self, sym_mmap, 0, NULL);
	}
}

// invoke yield(element, i, j) and replace the element by the return for all the elements
#define DEFINE_MMAP_WI2_FUNC(id, type, func, func2) static void \
km_mmap_with_index2_func_##id(type *elm, int i, int j, void *null) \
{ \
	*elm = func2(rb_yield(rb_ary_new3(3, func(*elm), INT2NUM(i), INT2NUM(j)))); \
}
DEFINE_MMAP_WI2_FUNC(d, double, rb_float_new, NUM2DBL)
DEFINE_MMAP_WI2_FUNC(z, COMPLEX, km_c2v, km_v2c)
DEFINE_MMAP_WI2_FUNC(i, int, INT2NUM, NUM2INT)
DEFINE_MMAP_WI2_FUNC(b, bool, TF2V, RTEST)
DEFINE_MMAP_WI2_FUNC(v, VALUE, ITSELF, ITSELF)
// alias map_with_index
VALUE
kmm_mat_mmap_with_index2_dest(VALUE self)
{
	if ( rb_block_given_p() ) {
		SMAT *smat = km_mat2smat(self);
		VT_SWITCH( smat->vtype,
			km_smat_each_with_index_d(smat, km_mmap_with_index2_func_d, NULL);,
			km_smat_each_with_index_z(smat, km_mmap_with_index2_func_z, NULL);,
			km_smat_each_with_index_i(smat, km_mmap_with_index2_func_i, NULL);,
			km_smat_each_with_index_b(smat, km_mmap_with_index2_func_b, NULL);,
			km_smat_each_with_index_v(smat, km_mmap_with_index2_func_v, NULL);
		);
		return Qnil;
	} else {
		return rb_enumeratorize(self, sym_mmap_with_index2, 0, NULL);
	}
}

// make a row-major Array of Arrays (Mat#to_a is Enumerable#to_a and it is not the same as this)
#define TA_LOOP(id, func) for ( int i=0; i<smat->m; i++ ) { \
	VALUE row = rb_ary_new2(smat->n); \
	for ( int j=0; j<smat->n; j++ ) { \
		rb_ary_push(row, func(smat->id##body[INDEX(smat, i, j)])); \
	} \
	rb_ary_push(ret, row); \
}
VALUE
kmm_mat_to_ary(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	VALUE ret = rb_ary_new2(smat->m);
	VT_SWITCH( smat->vtype,
		TA_LOOP(d, rb_float_new),
		TA_LOOP(z, km_c2v),
		TA_LOOP(i, INT2NUM),
		TA_LOOP(b, TF2V),
		TA_LOOP(v, ITSELF)
	);
	return ret;
}

// make a String
#define TS_LOOP(tid, func) for ( int i=0; i<smat->m; i++ ) { \
	VALUE row = rb_ary_new2(smat->n); \
	for ( int j=0; j<smat->n; j++ ) { \
		rb_ary_push(row, rb_funcall(func(ENTITY(smat, tid, i, j)), id, 0)); \
	} \
	rb_str_cat2(ret, " ["); \
	rb_str_concat(ret, rb_funcall(row, id_join, 1, sep)); \
	rb_str_cat2(ret, "],\n"); \
}
static VALUE
km_smat_to_s(const SMAT *smat, ID id, VALUE vtsym, VALUE stsym)
{
	VALUE ret = rb_utf8_str_new_cstr("Mat[\n");
	VALUE sep = rb_utf8_str_new_cstr(", ");
	VT_SWITCH( smat->vtype,
		TS_LOOP(d, rb_float_new),
		TS_LOOP(z, km_c2v),
		TS_LOOP(i, INT2NUM),
		TS_LOOP(b, TF2V),
		TS_LOOP(v, ITSELF)
	);
	rb_str_cat2(ret, "] (");
	rb_str_concat(ret, rb_funcall(vtsym, id_inspect, 0));
	rb_str_cat2(ret, ", ");
	rb_str_concat(ret, rb_funcall(stsym, id_inspect, 0));
	rb_str_cat2(ret, ")");
	return ret;
}
#define TSP_LOOP(tid, func) for ( int i=0; i<smat->m; i++ ) { \
	VALUE row = rb_ary_new2(smat->n); \
	for ( int j=0; j<smat->n; j++ ) { \
		rb_ary_push(row, rb_funcall(format, id_op_percent, 1, func(ENTITY(smat, tid, i, j)))); \
	} \
	rb_str_cat2(ret, " ["); \
	rb_str_concat(ret, rb_funcall(row, id_join, 1, sep)); \
	rb_str_cat2(ret, "],\n"); \
}
static VALUE
km_smat_to_s_percent(const SMAT *smat, VALUE vtsym, VALUE stsym, VALUE format)
{
	VALUE ret = rb_utf8_str_new_cstr("Mat[\n");
	VALUE sep = rb_utf8_str_new_cstr(", ");
	VT_SWITCH( smat->vtype,
		TSP_LOOP(d, rb_float_new),
		TSP_LOOP(z, km_c2v),
		TSP_LOOP(i, INT2NUM),
		TSP_LOOP(b, TF2V),
		TSP_LOOP(v, ITSELF)
	);
	rb_str_cat2(ret, "] (");
	rb_str_concat(ret, rb_funcall(vtsym, id_inspect, 0));
	rb_str_cat2(ret, ", ");
	rb_str_concat(ret, rb_funcall(stsym, id_inspect, 0));
	rb_str_cat2(ret, ")");
	return ret;
}
// alias to_str
VALUE
kmm_mat_to_s(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 0, 1);
	if ( argc == 0 ) {
		return km_smat_to_s(km_mat2smat(self), id_to_s, kmm_mat_vtype(self), kmm_mat_stype(self));
	} else {
		return km_smat_to_s_percent(km_mat2smat(self), kmm_mat_vtype(self), kmm_mat_stype(self), argv[0]);
	}
}
VALUE
kmm_mat_inspect(VALUE self)
{
	return km_smat_to_s(km_mat2smat(self), id_inspect, kmm_mat_vtype(self), kmm_mat_stype(self));
}

// symmetrize `self' by (A+A')/2 where A is `self' on entry
VALUE
kmm_mat_symmetrize_dest(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( smat->m != smat->n ) {
		rb_raise(km_eDim, "non-square matrix cannot be symmetrized");
	}
	if ( smat->vtype == VT_DOUBLE ) {
		for ( int i=0; i<smat->m-1; i++ ) {
			for ( int j=i+1; j<smat->m; j++ ) {
				smat->dbody[i+j*smat->ld] = ( smat->dbody[i+j*smat->ld] + smat->dbody[j+i*smat->ld] ) * 0.5;
				smat->dbody[j+i*smat->ld] = smat->dbody[i+j*smat->ld];
			}
		}
	} else if ( smat->vtype == VT_COMPLEX ) {
		for ( int i=0; i<smat->m-1; i++ ) {
			for ( int j=i+1; j<smat->m; j++ ) {
				smat->zbody[i+j*smat->ld] = ( smat->zbody[i+j*smat->ld] + smat->zbody[j+i*smat->ld] ) * 0.5;
				smat->zbody[j+i*smat->ld] = smat->zbody[i+j*smat->ld];
			}
		}
	} else {
		rb_raise(km_eVT, "symmetrize is available only for float or complex matricies");
	}
	return self;
}

// judge whether `self' is symmetric or not
// if the argument is true `self' will be completely symmetric by Mat#symmetrize!
// float, complex: consider that a is sufficiently near to b iff |a-b|/(|a|+|b|) < eps*8
// int, bool: a is equal to b iff a == b
// object: consider that a is sufficiently near to b iff a.respond_to?(:near?) ? a.near?(b) : a==b
static bool
km_near_rel_d(double a, double b)
{
	if ( a == b ) {
		return true;
	} else {
		return ( fabs(a-b)/(fabs(a)+fabs(b)) < DBL_EPSILON*8.0 );
	}
}
static bool
km_near_rel_z(COMPLEX a, COMPLEX b)
{
	if ( a == b ) {
		return true;
	} else {
		return ( cabs(a-b)/(cabs(a)+cabs(b)) < DBL_EPSILON*8.0 );
	}
}
static bool
km_near_v(VALUE a, VALUE b)
{
	if ( rb_respond_to(a, id_near_p) ) {
		return RTEST(rb_funcall(a, id_near_p, 1, b));
	} else {
		return RTEST(rb_funcall(a, id_op_eq, 1, b));
	}
}
#define SYM_BODY(id, func) for ( int i=0; i<(smat->n)-1; i++ ) { \
	for ( int j=i+1; j<(smat->n); j++ ) { \
		if ( !func(smat->id##body[i+j*(smat->ld)], smat->id##body[j+i*(smat->ld)]) ) { \
			return Qfalse; \
		} \
	} \
}
VALUE
kmm_mat_symmetry_p(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 0, 1);
	bool symmetrize;
	if ( argc == 0 ) {
		symmetrize = false;
	} else {
		symmetrize = RTEST(argv[0]);
	}
	SMAT *smat = km_mat2smat(self);
	if ( smat->m != smat->n ) {
		return Qfalse;
	}
	VT_SWITCH( smat->vtype,
		SYM_BODY(d, km_near_rel_d)
		if ( symmetrize ) { kmm_mat_symmetrize_dest(self); },
		SYM_BODY(z, km_near_rel_z)
		if ( symmetrize ) { kmm_mat_symmetrize_dest(self); },
		SYM_BODY(i, SAME),
		SYM_BODY(b, SAME),
		SYM_BODY(v, km_near_v)
	);
	return Qtrue;
}

// judge whether `self' is near to `other'
// consider a is near to b iff |a-b| < tol
// the default value of tol is (maximub absolute value of each element) * max(m, n) * eps
static bool
km_near_abs_d(double a, double b, double tol)
{
	if ( a == b ) {
		return true;
	} else {
		return ( fabs(a-b) < tol );
	}
}
static bool
km_near_abs_z(COMPLEX a, COMPLEX b, double tol)
{
	if ( a == b ) {
		return true;
	} else {
		return ( cabs(a-b) < tol );
	}
}
struct km_near_arg {
	bool ret;
	double tol;
};
static void
km_near_func_d(double *a, double *b, void *data)
{
	struct km_near_arg *arg = (struct km_near_arg *)data;
	arg->ret = (arg->ret && km_near_abs_d(*a, *b, arg->tol));
}
static void
km_near_func_z(COMPLEX *a, COMPLEX *b, void *data)
{
	struct km_near_arg *arg = (struct km_near_arg *)data;
	arg->ret = (arg->ret && km_near_abs_z(*a, *b, arg->tol));
}
static void
km_near_func_i(int *a, int *b, void *data)
{
	struct km_near_arg *arg = (struct km_near_arg *)data;
	arg->ret = (arg->ret && *a==*b);
}
static void
km_near_func_b(bool *a, bool *b, void *data)
{
	struct km_near_arg *arg = (struct km_near_arg *)data;
	arg->ret = (arg->ret && *a==*b);
}
static void
km_near_func_v(VALUE *a, VALUE *b, void *data)
{
	struct km_near_arg *arg = (struct km_near_arg *)data;
	arg->ret = (arg->ret && km_near_v(*a, *b));
}
VALUE
kmm_mat_near_p(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 1, 2);
	SMAT *ss = km_mat2smat(self);
	SMAT *so = km_mat2smat(argv[0]);
	if ( ss->vtype != so->vtype || ss->m != so->m || ss->n != so->n ) {
		return Qfalse;
	}
	struct km_near_arg data = {true, 0.0};
	if ( argc == 1 ) {
		double ns = NUM2DBL(kmm_mat_norm_einf(self)), no = NUM2DBL(kmm_mat_norm_einf(argv[0]));
		data.tol = MAX(ns, no) * MAX(ss->m, ss->n) * DBL_EPSILON;
	} else {
		data.tol = NUM2DBL(argv[1]);
	}
	VT_SWITCH( ss->vtype,
		km_smat_each2_d(ss, so, km_near_func_d, (void *)&data);,
		km_smat_each2_z(ss, so, km_near_func_z, (void *)&data);,
		km_smat_each2_i(ss, so, km_near_func_i, (void *)&data);,
		km_smat_each2_b(ss, so, km_near_func_b, (void *)&data);,
		km_smat_each2_v(ss, so, km_near_func_v, (void *)&data);
	);
	return data.ret;
}
