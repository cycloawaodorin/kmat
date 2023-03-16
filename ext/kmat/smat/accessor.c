#include "../kmat.h"

static int
mod(int a, int b)
{
	int r = a % b;
	if ( r < 0 ) {
		return r+b;
	} else {
		return r;
	}
}

// returns (i, j)-th element
VALUE
kmm_mat_get_value(VALUE self, VALUE vi, VALUE vj)
{
	SMAT *smat = km_mat2smat(self);
	size_t i = NUM2ZU(vi), j = NUM2ZU(vj);
	if ( smat->m <= i || smat->n <= j ) {
		rb_raise(rb_eIndexError, "index (%zu, %zu) is out of range (%zu, %zu)", i, j, smat->m, smat->n);
	}
	VT_SWITCH( smat->vtype,
		return rb_float_new(ENTITY(smat, d, i, j));,
		return km_c2v(ENTITY(smat, z, i, j));,
		return INT2NUM(ENTITY(smat, i, i, j));,
		return TF2V(ENTITY(smat, b, i, j));,
		return ENTITY(smat, v, i, j);
	);
}

// make a serial-submatrix with shape (`vm', `vn')
// (0, 0)-th element of the return is (`vi', `vj')-th element of `self'
VALUE
kmm_mat_get_ssub(VALUE self, VALUE vi, VALUE vj, VALUE vm, VALUE vn)
{
	return km_Mat_ssub(NUM2ZU(vi), NUM2ZU(vj), NUM2ZU(vm), NUM2ZU(vn), self);
}

// make a random-submatrix by :ivec indecies `vi', `vj'
// (i, j)-th element of the return is (`vi'[i], `vj'[j])-th element of `self'
// if both `vi' and `vj' are serial, serial-submatrix is made
VALUE
kmm_mat_get_rsub(VALUE self, VALUE vi, VALUE vj)
{
	SMAT *si = km_mat2smat(vi), *sj = km_mat2smat(vj);
	if ( si->vtype != VT_INT || sj->vtype != VT_INT ) {
		rb_raise(km_eVT, "index arguments must be int vectors");
	}
	if ( !VECTOR_P(si) || !VECTOR_P(sj) ) {
		rb_raise(km_eDim, "index arguments must be int vectors");
	}
	if ( si->stype != ST_FULL ) {
		si = km_mat2smat(rb_obj_dup(vi));
	}
	if ( sj->stype != ST_FULL ) {
		sj = km_mat2smat(rb_obj_dup(vj));
	}
	int *is = si->ibody, *js = sj->ibody;
	bool flg = true;
	size_t m = LENGTH(si), n = LENGTH(sj);
	for (size_t i=0; i<m-1; i++) {
		if ( is[i+1]-is[i] != 1 ) {
			flg = false;
			break;
		}
	}
	if ( flg ) {
		for (size_t i=0; i<n-1; i++) {
			if ( js[i+1]-js[i] != 1 ) {
				flg = false;
				break;
			}
		}
	}
	if ( flg ) {
		return km_Mat_ssub(i2s(is[0]), i2s(js[0]), m, n, self);
	} else {
		return km_Mat_rsub1(m, n, is, js, self);
	}
}

// make a random-submatrix by :vmat index `vi'
// (i, j)-th element of the return is (`vi'[i, j][0], `vi'[i, j][1])-th element of `self'
VALUE
kmm_mat_get_rsub2(VALUE self, VALUE vi)
{
	SMAT *smat = km_mat2smat(self);
	SMAT *vmat = km_mat2smat(vi);
	if ( vmat->vtype != VT_VALUE ) {
		rb_raise(km_eVT, "index argument must be an object matrix");
	}
	size_t m = vmat->m, n = vmat->n, idx;
	size_t m2 = smat->m, n2 = smat->n;
	size_t *is, *js;
	is = ALLOCA_N(size_t, n*m);
	js = ALLOCA_N(size_t, n*m);
	VALUE ij;
	for ( size_t i=0; i<m; i++ ) { for ( size_t j=0; j<n; j++ ) {
		ij = ENTITY(vmat, v, i, j);
		if ( (TYPE(ij) != T_ARRAY) || (RARRAY_LEN(ij) != 2) ) {
			rb_raise(rb_eIndexError, "value of vmat-index must be an Array with 2 elements");
		}
		idx = i+j*m;
		is[idx] = i2s(mod( NUM2INT(rb_ary_entry(ij, 0)), s2i(m2) ));
		js[idx] = i2s(mod( NUM2INT(rb_ary_entry(ij, 1)), s2i(n2) ));
	} }
	return km_Mat_rsub2(m, n, is, js, self);
}

// make a random-submatrix by :bmat index `vi'
// if `vi'[i, j] is true, the return contains `self'[i, j]
// if `vi'[i, j] is false, the return does not contain `self'[i, j]
// the return is a vector which length is the number of trues in `vi'
static void
km_gbb_len(bool *eb, void *data)
{
	if ( *eb ) { *((size_t *)data) += 1; };
}
struct km_gbb_arg {
	size_t l;
	union {
		double **dpbody;
		COMPLEX **zpbody;
		int **ipbody;
		bool **bpbody;
		VALUE **vpbody;
	};
};
static void
km_gbb_func_d(double *ea, const bool *eb, void *data_)
{
	struct km_gbb_arg *data = (struct km_gbb_arg *)data_;
	if ( *eb ) {
		data->dpbody[data->l] = ea;
		data->l += 1;
	}
}
static void
km_gbb_func_z(COMPLEX *ea, const bool *eb, void *data_)
{
	struct km_gbb_arg *data = (struct km_gbb_arg *)data_;
	if ( *eb ) {
		data->zpbody[data->l] = ea;
		data->l += 1;
	}
}
static void
km_gbb_func_i(int *ea, const bool *eb, void *data_)
{
	struct km_gbb_arg *data = (struct km_gbb_arg *)data_;
	if ( *eb ) {
		data->ipbody[data->l] = ea;
		data->l += 1;
	}
}
static void
km_gbb_func_b(bool *ea, const bool *eb, void *data_)
{
	struct km_gbb_arg *data = (struct km_gbb_arg *)data_;
	if ( *eb ) {
		data->bpbody[data->l] = ea;
		data->l += 1;
	}
}
static void
km_gbb_func_v(VALUE *ea, const bool *eb, void *data_)
{
	struct km_gbb_arg *data = (struct km_gbb_arg *)data_;
	if ( *eb ) {
		data->vpbody[data->l] = ea;
		data->l += 1;
	}
}
VALUE
kmm_mat_get_by_bmat(VALUE self, VALUE vi)
{
	SMAT *smat = km_mat2smat(self);
	SMAT *bmat = km_mat2smat(vi);
	struct km_gbb_arg data = {0, {NULL}};
	km_smat_each_b(bmat, km_gbb_len, (void *)&(data.l));
	VALUE ret = km_Mat_alloc(km_cMat);
	SMAT *sr = km_mat2smat(ret);
	km_smat_alloc_pbody(sr, data.l, 1, smat->vtype);
	data.l = 0;
	VT_SWITCH( smat->vtype,
		data.dpbody = sr->dpbody; km_smat_each2_dcb(smat, bmat, km_gbb_func_d, &data);,
		data.zpbody = sr->zpbody; km_smat_each2_zcb(smat, bmat, km_gbb_func_z, &data);,
		data.ipbody = sr->ipbody; km_smat_each2_icb(smat, bmat, km_gbb_func_i, &data);,
		data.bpbody = sr->bpbody; km_smat_each2_bcb(smat, bmat, km_gbb_func_b, &data);,
		data.vpbody = sr->vpbody; km_smat_each2_vcb(smat, bmat, km_gbb_func_v, &data);
	);
	if ( smat->stype == ST_FULL ) {
		smat->may_have_sub = true;
		sr->parent = self;
	} else {
		sr->parent = smat->parent;
	}
	km_infect_frozen(self, ret);
	return ret;
}

// set val to (i, j)-th element
VALUE
kmm_mat_set_value(VALUE self, VALUE vi, VALUE vj, VALUE val)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	size_t i = NUM2ZU(vi), j = NUM2ZU(vj);
	if ( smat->m<i || smat->n<j ) {
		rb_raise(rb_eIndexError, "index (%zu, %zu) is out of range (%zu, %zu)", i, j, smat->m, smat->n);
	}
	if ( smat->stype == ST_RSUB ) {
		VT_SWITCH( smat->vtype,
			ENTITYr0(smat, d, INDEX(smat, i, j)) = NUM2DBL(val);,
			ENTITYr0(smat, z, INDEX(smat, i, j)) = km_v2c(val);,
			ENTITYr0(smat, i, INDEX(smat, i, j)) = NUM2INT(val);,
			ENTITYr0(smat, b, INDEX(smat, i, j)) = RTEST(val);,
			ENTITYr0(smat, v, INDEX(smat, i, j)) = val;
		);
	} else {
		VT_SWITCH( smat->vtype,
			ENTITYd0(smat, d, INDEX(smat, i, j)) = NUM2DBL(val);,
			ENTITYd0(smat, z, INDEX(smat, i, j)) = km_v2c(val);,
			ENTITYd0(smat, i, INDEX(smat, i, j)) = NUM2INT(val);,
			ENTITYd0(smat, b, INDEX(smat, i, j)) = RTEST(val);,
			ENTITYd0(smat, v, INDEX(smat, i, j)) = val;
		);
	}
	return self;
}

// copy the values of `other' to `self'
VALUE
kmm_mat_copy_from(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *dest = km_mat2smat(self), *src = km_mat2smat(other);
	CHECK_SAME_SIZE(dest, src);
	if ( dest->vtype == src->vtype ) {
		if ( ( dest->stype != ST_FULL && ( dest->parent ==  src->parent || dest->parent == other ) ) || src->parent == self ) {
			km_smat_copy(dest, km_mat2smat(rb_obj_dup(other)));
		} else {
			km_smat_copy(dest, src);
		}
	} else {
		VT_SWITCH( dest->vtype,
			km_smat_copy(dest, km_mat2smat(kmm_mat_to_fmat(other)));,
			km_smat_copy(dest, km_mat2smat(kmm_mat_to_cmat(other)));,
			km_smat_copy(dest, km_mat2smat(kmm_mat_to_imat(other)));,
			km_smat_copy(dest, km_mat2smat(kmm_mat_to_bmat(other)));,
			km_smat_copy(dest, km_mat2smat(kmm_mat_to_omat(other)));
		);
	}
	return self;
}
VALUE
kmm_mat_tcopy_from(VALUE self, VALUE other)
{
	km_check_frozen(self);
	SMAT *dest = km_mat2smat(self), *src = km_mat2smat(other);
	if ( dest->m!=src->n || dest->n!=src->m ) {
		rb_raise(km_eDim, "transposed sizes must be the same, (%zu, %zu) != (%zu, %zu)", (dest)->m, (dest)->n, (src)->n, (src)->m);
	}
	dest->trans = !dest->trans;
	SWAP(size_t, dest->m, dest->n);
	if ( dest->vtype == src->vtype ) {
		if ( ( dest->stype != ST_FULL && ( dest->parent ==  src->parent || dest->parent == other ) ) || src->parent == self ) {
			km_smat_copy(dest, km_mat2smat(rb_obj_dup(other)));
		} else {
			km_smat_copy(dest, src);
		}
	} else {
		VT_SWITCH( dest->vtype,
			km_smat_copy(dest, km_mat2smat(kmm_mat_to_fmat(other)));,
			km_smat_copy(dest, km_mat2smat(kmm_mat_to_cmat(other)));,
			km_smat_copy(dest, km_mat2smat(kmm_mat_to_imat(other)));,
			km_smat_copy(dest, km_mat2smat(kmm_mat_to_bmat(other)));,
			km_smat_copy(dest, km_mat2smat(kmm_mat_to_omat(other)));
		);
	}
	dest->trans = !dest->trans;
	SWAP(size_t, dest->m, dest->n);
	return self;
}

// make a random-submatrix vector which contains diagonal elements of `self'
VALUE
kmm_mat__diag(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	VALUE ret = km_Mat_alloc(km_cMat);
	SMAT *sr = km_mat2smat(ret);
	km_smat_alloc_pbody(sr, MIN(smat->m, smat->n), 1, smat->vtype);
	if ( smat->stype == ST_RSUB ) {
		for ( size_t i=0; i<sr->m; i++ ) {
			sr->pbody[i] = smat->pbody[i+i*(smat->ld)];
		}
	} else {
		VT_SWITCH( smat->vtype,
			for ( size_t i=0; i<sr->m; i++ ) {
				sr->pbody[i] = smat->dbody+(i+i*smat->ld);
			},
			for ( size_t i=0; i<sr->m; i++ ) {
				sr->pbody[i] = smat->zbody+(i+i*smat->ld);
			},
			for ( size_t i=0; i<sr->m; i++ ) {
				sr->pbody[i] = smat->ibody+(i+i*smat->ld);
			},
			for ( size_t i=0; i<sr->m; i++ ) {
				sr->pbody[i] = smat->bbody+(i+i*smat->ld);
			},
			for ( size_t i=0; i<sr->m; i++ ) {
				sr->pbody[i] = smat->vbody+(i+i*smat->ld);
			}
		);
	}
	if ( smat->stype == ST_FULL ) {
		sr->parent = self;
		smat->may_have_sub = true;
	} else {
		sr->parent = smat->parent;
	}
	km_infect_frozen(self, ret);
	return ret;
}

VALUE
kmm_mat__diag_ul(VALUE self, VALUE vk)
{
	SMAT *smat = km_mat2smat(self);
	VALUE ret = km_Mat_alloc(km_cMat);
	SMAT *sr = km_mat2smat(ret);
	long k = NUM2LONG(vk);
	if ( k == 0 ) { return kmm_mat__diag(self); }
	long len;
	size_t i_s, j_s;
	if ( 0 < k ) {
		if ( k < s2l(smat->n)-s2l(smat->m) ) {
			len = s2l(smat->m);
		} else {
			len = s2l(smat->n)-k;
		}
		if ( smat->trans ) {
			i_s = l2s(k); j_s = 0;
		} else {
			i_s = 0; j_s = l2s(k);
		}
	} else {
		if ( k < s2l(smat->n)-s2l(smat->m) ) {
			len = s2l(smat->m)+k;
		} else {
			len = s2l(smat->n);
		}
		if ( smat->trans ) {
			i_s = 0; j_s = l2s(-k);
		} else {
			i_s = l2s(-k); j_s = 0;
		}
	}
	if ( len <= 0 ) {
		rb_raise(rb_eArgError, "given offset %ld exceeds range [-%ld, %ld]", k, s2l(smat->m)-1, s2l(smat->n)-1);
	}
	size_t len_s = l2s(len);
	km_smat_alloc_pbody(sr, len_s, 1, smat->vtype);
	if ( smat->stype == ST_RSUB ) {
		for ( size_t i=0; i<len_s; i++ ) {
			sr->pbody[i] = smat->pbody[i_s+i+(j_s+i)*smat->ld];
		}
	} else {
		VT_SWITCH( smat->vtype,
			for ( size_t i=0; i<len_s; i++ ) {
				sr->pbody[i] = smat->dbody+(i_s+i+(j_s+i)*smat->ld);
			},
			for ( size_t i=0; i<len_s; i++ ) {
				sr->pbody[i] = smat->zbody+(i_s+i+(j_s+i)*smat->ld);
			},
			for ( size_t i=0; i<len_s; i++ ) {
				sr->pbody[i] = smat->ibody+(i_s+i+(j_s+i)*smat->ld);
			},
			for ( size_t i=0; i<len_s; i++ ) {
				sr->pbody[i] = smat->bbody+(i_s+i+(j_s+i)*smat->ld);
			},
			for ( size_t i=0; i<len_s; i++ ) {
				sr->pbody[i] = smat->vbody+(i_s+i+(j_s+i)*smat->ld);
			}
		);
	}
	if ( smat->stype == ST_FULL ) {
		sr->parent = self;
		smat->may_have_sub = true;
	} else {
		sr->parent = smat->parent;
	}
	km_infect_frozen(self, ret);
	return ret;
}

// Mat#[]
// get an element or make a submatrix
typedef enum {
	IT_SERIAL,
	IT_IVEC,
	IT_BMAT,
	IT_INT,
	IT_ZERO,
	IT_VMAT,
	IT_END
} ITSYM;
static void
km_index_treatment(VALUE *oidx, ITSYM *itsym, VALUE iidx, size_t size, bool convert)
{
	if ( iidx == Qnil ) {
		*oidx = rb_ary_new3(2, INT2NUM(0), ZU2NUM(size));
		*itsym = IT_SERIAL;
	} else if ( rb_obj_is_kind_of(iidx, km_cMat) ) {
		SMAT *si = km_mat2smat(iidx);
		if ( si->vtype == VT_INT ) {
			if ( VECTOR_P(si) ) {
				if ( LENGTH(si) == 0 ) {
					*oidx = Qnil;
					*itsym = IT_ZERO;
				} else {
					*oidx = iidx;
					*itsym = IT_IVEC;
				}
			} else {
				rb_raise(rb_eArgError, "int-Mat-index must be a vector");
			}
		} else if ( si->vtype == VT_BOOL ) {
			if ( convert && VECTOR_P(si) ) {
				if ( LENGTH(si) == size ) {
					size_t oi_size = 0;
					if ( si->stype == ST_RSUB ) {
						for ( size_t i=0; i<size; i++ ) {
							if ( *((si->bpbody)[i]) ) {
								oi_size++;
							}
						}
						*oidx = km_Mat(oi_size,1, VT_INT);
						SMAT *oi = km_mat2smat(*oidx);
						oi_size = 0;
						for ( size_t i=0; i<size; i++ ) {
							if ( *((si->bpbody)[i]) ) {
								oi->ibody[oi_size] = s2i(i);
								oi_size++;
							}
						}
					} else {
						for ( size_t i=0; i<size; i++ ) {
							if ( (si->bbody)[i] ) {
								oi_size++;
							}
						}
						*oidx = km_Mat(oi_size, 1, VT_INT);
						SMAT *oi = km_mat2smat(*oidx);
						oi_size = 0;
						for ( size_t i=0; i<size; i++ ) {
							if ( (si->bbody)[i] ) {
								oi->ibody[oi_size] = s2i(i);
								oi_size++;
							}
						}
					}
					*itsym = IT_IVEC;
				} else {
					rb_raise(km_eDim, "boolean-Mat-index length must match with object size");
				}
			} else {
				*oidx = iidx;
				*itsym = IT_BMAT;
			}
		} else if ( si->vtype == VT_VALUE ) {
			*oidx = iidx;
			*itsym = IT_VMAT;
		} else {
			rb_raise(km_eVT, "float or complex matrix cannot be an index");
		}
	} else if ( rb_obj_is_kind_of(iidx, rb_cArray) ) {
		long ii_size = RARRAY_LEN(iidx), size_l = s2l(size);
		if ( ii_size == 0 ) {
			*oidx = rb_ary_new3(2, INT2NUM(0), ZU2NUM(size));
			*itsym = IT_SERIAL;
		} else if ( rb_ary_entry(iidx, 0) == Qtrue || rb_ary_entry(iidx, 0) == Qfalse ) {
			if ( ii_size == size_l ) {
				size_t oi_size = 0;
				for ( long i=0; i<size_l; i++ ) {
					if ( RTEST(rb_ary_entry(iidx, i)) ) {
						oi_size++;
					}
				}
				*oidx = km_Mat(oi_size, 1, VT_INT);
				SMAT *oi = km_mat2smat(*oidx);
				oi_size = 0;
				for ( long i=0; i<size_l; i++ ) {
					if ( RTEST(rb_ary_entry(iidx, i)) ) {
						oi->ibody[oi_size] = (int)i;
						oi_size++;
					}
				}
				*itsym = IT_IVEC;
			} else {
				rb_raise(km_eDim, "boolean-Array-index length must match with object size");
			}
		} else {
			*oidx = km_Mat(l2s(ii_size), 1, VT_INT);
			SMAT *oi = km_mat2smat(*oidx);
			for ( long i=0; i<ii_size; i++ ) {
				(oi->ibody)[i] = mod(NUM2INT(rb_ary_entry(iidx, i)), s2i(size));
			}
			*itsym = IT_IVEC;
		}
	} else if ( rb_obj_is_kind_of(iidx, rb_cInteger) ) {
		*oidx = INT2NUM(mod(NUM2INT(iidx), s2i(size)));
		*itsym = IT_INT;
	} else if ( rb_obj_is_kind_of(iidx, rb_cRange) ) {
		int f = mod(NUM2INT(rb_funcall(iidx, id_first, 0)), s2i(size));
		int l = NUM2INT(rb_funcall(iidx, id_last, 0));
		if ( RTEST(rb_funcall(iidx, id_exclude_end_p, 0)) ) {
			l--;
		}
		l = mod(l, s2i(size));
		*oidx = rb_ary_new3(2, INT2NUM(f), INT2NUM(l-f+1));
		*itsym = IT_SERIAL;
	} else {
		rb_raise(rb_eArgError, "unknown index type");
	}
}

VALUE
kmm_mat_bracket(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 0, 2);
	SMAT *smat = km_mat2smat(self);
	if ( LENGTH(smat) == 0 ) {
		rb_raise(km_eDim, "Mat#[] is not available for 0-size matricies");
	}
	if ( argc == 2 ) {
		VALUE ri, ci;
		ITSYM rit, cit;
		km_index_treatment(&ri, &rit, argv[0], smat->m, true);
		km_index_treatment(&ci, &cit, argv[1], smat->n, true);
		if ( rit == IT_INT ) {
			if ( cit == IT_INT ) {
				return kmm_mat_get_value(self, ri, ci);
			} else if ( cit == IT_SERIAL ) {
				return km_Mat_ssub(NUM2ZU(ri), NUM2ZU(rb_ary_entry(ci, 0)), 1, NUM2ZU(rb_ary_entry(ci, 1)), self);
			} else if (  cit == IT_IVEC ) {
				VALUE riv = km_Mat(1, 1, VT_INT);
				(km_mat2smat(riv)->ibody)[0] = NUM2INT(ri);
				return kmm_mat_get_rsub(self, riv, ci);
			}
		} else if ( rit == IT_SERIAL ) {
			if ( cit == IT_INT ) {
				return km_Mat_ssub(NUM2ZU(rb_ary_entry(ri, 0)), NUM2ZU(ci), NUM2ZU(rb_ary_entry(ri, 1)), 1, self);
			} else if ( cit == IT_SERIAL ) {
				return kmm_mat_get_ssub(self, rb_ary_entry(ri, 0), rb_ary_entry(ci, 0), rb_ary_entry(ri, 1), rb_ary_entry(ci, 1));
			} else if ( cit == IT_IVEC ) {
				size_t ri1 = NUM2ZU(rb_ary_entry(ri, 1));
				VALUE riv = km_Mat(ri1, 1, VT_INT);
				size_t ri0 = NUM2ZU(rb_ary_entry(ri, 0));
				SMAT *sri = km_mat2smat(riv);
				for ( size_t i=0; i<ri1; i++ ) {
					(sri->ibody)[i] = s2i(i+ri0);
				}
				return kmm_mat_get_rsub(self, riv, ci);
			}
		} else if ( rit == IT_IVEC ) {
			if ( cit == IT_INT ) {
				VALUE civ = km_Mat(1, 1, VT_INT);
				(km_mat2smat(civ)->ibody)[0] = NUM2INT(ci);
				return kmm_mat_get_rsub(self, ri, civ);
			} else if ( cit == IT_SERIAL ) {
				size_t ci1 = NUM2ZU(rb_ary_entry(ci, 1));
				VALUE civ = km_Mat(ci1, 1, VT_INT);
				size_t ci0 = NUM2ZU(rb_ary_entry(ci, 0));
				SMAT *sci = km_mat2smat(civ);
				for ( size_t i=0; i<ci1; i++ ) {
					(sci->ibody)[i] = s2i(i+ci0);
				}
				return kmm_mat_get_rsub(self, ri, civ);
			} else if ( cit == IT_IVEC ) {
				return kmm_mat_get_rsub(self, ri, ci);
			}
		} else if ( rit == IT_ZERO ) {
			size_t cs;
			if ( cit == IT_INT ) {
				cs = 1;
			} else if ( cit == IT_SERIAL ) {
				cs = NUM2ZU(rb_ary_entry(ci, 1));
			} else if ( cit == IT_IVEC ) {
				cs = LENGTH(km_mat2smat(ci));
			} else if ( cit == IT_ZERO ) {
				cs = 0;
			} else {
				rb_raise(rb_eArgError, "illegal index type combination [%d, %d]", rit, cit);
			}
			return km_Mat(0, cs, smat->vtype);
		}
		rb_raise(rb_eArgError, "illegal index type combination [%d, %d]", rit, cit);
	} else if ( argc == 1 ) {
		VALUE i;
		ITSYM it;
		km_index_treatment(&i, &it, argv[0], LENGTH(smat), false);
		if ( it == IT_BMAT ) {
			SMAT *si  = km_mat2smat(i);
			if ( SAME_SIZE(smat, si) ) {
				return kmm_mat_get_by_bmat(self, i);
			} else {
				rb_raise(km_eDim, "boolean index size (%zu, %zu) must be the same as matrix size (%zu, %zu)", si->m, si->n, smat->m, smat->n);
			}
		} else if ( it == IT_INT ) {
			if ( smat->n == 1 ) {
				return kmm_mat_get_value(self, i, INT2NUM(0));
			} else if ( smat->m == 1 ) {
				return kmm_mat_get_value(self, INT2NUM(0), i);
			} else {
				rb_raise(rb_eArgError, "single integer index is available only for vectors, not for (%zu, %zu) matricies", smat->m, smat->n);
			}
		} else if ( it == IT_SERIAL ) {
			if ( smat->n == 1 ) {
				return km_Mat_ssub(NUM2ZU(rb_ary_entry(i, 0)), 0, NUM2ZU(rb_ary_entry(i, 1)), 1, self);
			} else if ( smat->m == 1 ) {
				return km_Mat_ssub(0, NUM2ZU(rb_ary_entry(i, 0)), 1, NUM2ZU(rb_ary_entry(i, 1)), self);
			} else {
				rb_raise(rb_eArgError, "single serial index is available only for vectors, not for (%zu, %zu) matricies", smat->m, smat->n);
			}
		} else if ( it == IT_IVEC ) {
			if ( smat->n == 1 ) {
				VALUE i2 = km_Mat(1, 1, VT_INT);
				(km_mat2smat(i2)->ibody)[0] = 0;
				return kmm_mat_get_rsub(self, i, i2);
			} else if ( smat->m == 1 ) {
				VALUE i2 = km_Mat(1, 1, VT_INT);
				(km_mat2smat(i2)->ibody)[0] = 0;
				return kmm_mat_get_rsub(self, i2, i);
			} else {
				rb_raise(rb_eArgError, "single ivec index is available only for vectors, not for (%zu, %zu) matricies", smat->m, smat->n);
			}
		} else if ( it == IT_ZERO ) {
			if ( smat->n == 1 ) {
				return km_Mat(1, 0, smat->vtype);
			} else if ( smat->m == 1 ) {
				return km_Mat(0, 1, smat->vtype);
			} else {
				rb_raise(rb_eArgError, "single serial index is available only for vectors, not for (%zu, %zu) matricies", smat->m, smat->n);
			}
		} else if ( it == IT_VMAT ) {
			return kmm_mat_get_rsub2(self, i);
		} else {
			rb_raise(km_eInternal, "unknown index type %d", it);
		}
	} else {
		return km_Mat_ssub(0, 0, smat->m, smat->n, self);
	}
}

// alias []=
VALUE
kmm_mat_bracket_set(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 1, 3);
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( LENGTH(smat) == 0 ) {
		rb_raise(km_eDim, "Mat#[]= is not avialable for 0-size matricies");
	}
	if ( argc == 3 ) {
		if ( rb_obj_is_kind_of(argv[0], rb_cInteger) && rb_obj_is_kind_of(argv[1], rb_cInteger) ) {
			kmm_mat_set_value(self, INT2NUM(mod(NUM2INT(argv[0]), s2i(smat->m))), INT2NUM(mod(NUM2INT(argv[1]), s2i(smat->n))), argv[2]);
		} else {
			bool flg = smat->may_have_sub;
			VALUE target = kmm_mat_bracket(2, argv, self);
			if ( rb_obj_is_kind_of(argv[2], km_cMat) ) {
				kmm_mat_copy_from(target, argv[2]);
			} else {
				kmm_mat_fill(target, argv[2]);
			}
			kmm_mat__kill(target);
			smat->may_have_sub = flg;
		}
	} else if ( argc == 2 ) {
		if ( rb_obj_is_kind_of(argv[0], rb_cInteger) ) {
			if ( smat->n == 1 ) {
				kmm_mat_set_value(self, INT2NUM(mod(NUM2INT(argv[0]), s2i(smat->m))), INT2NUM(0), argv[1]);
			} else if ( smat->m == 1 ) {
				kmm_mat_set_value(self, INT2NUM(0), INT2NUM(mod(NUM2INT(argv[0]), s2i(smat->n))), argv[1]);
			} else {
				rb_raise(rb_eArgError, "setting value with single integer index is available only for vectors, not (%zu, %zu) matricies", smat->m, smat->n);
			}
		} else {
			bool flg = smat->may_have_sub;
			VALUE target = kmm_mat_bracket(1, argv, self);
			if ( rb_obj_is_kind_of(argv[1], km_cMat) ) {
				if ( rb_obj_is_kind_of(argv[0], km_cMat) && km_mat2smat(argv[0])->vtype == VT_BOOL ) {
					SMAT *smat2 = km_mat2smat(argv[1]);
					bool flg2 = smat2->may_have_sub;
					VALUE vsrc = kmm_mat_bracket(1, argv, argv[1]);
					kmm_mat_copy_from(target, vsrc);
					kmm_mat__kill(vsrc);
					smat2->may_have_sub = flg2;
				} else {
					kmm_mat_copy_from(target, argv[1]);
				}
			} else {
				kmm_mat_fill(target, argv[1]);
			}
			kmm_mat__kill(target);
			smat->may_have_sub = flg;
		}
	} else {
		if ( rb_obj_is_kind_of(argv[0], km_cMat) ) {
			kmm_mat_copy_from(self, argv[0]);
		} else {
			kmm_mat_fill(self, argv[0]);
		}
	}
	return argv[argc-1];
}
