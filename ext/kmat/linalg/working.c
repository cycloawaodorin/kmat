#include "../kmat.h"

void *
km_alloc_and_copy(SMAT *smat)
{
	void *ret = ruby_xcalloc(LENGTH(smat), km_sizeof_vt(smat->vtype));
	km_copy2work(ret, smat->m, smat);
	return ret;
}

// use `smat'->body if it can be used
// otherwise, calloc()
void
km_alloc_if_needed(SMAT *smat, LAWORK *lawork)
{
	if ( smat->stype != ST_RSUB && !(smat->trans) ) {
		lawork->body = smat->body;
		lawork->ld = smat->ld;
		lawork->need_to_free = false;
	} else {
		lawork->body = ruby_xcalloc(smat->m*smat->n, km_sizeof_vt(smat->vtype));
		lawork->ld = smat->m;
		lawork->need_to_free = true;
	}
}

// km_alloc_if_needed() and copy from `smat'->body if calloc() is called
void
km_alloc_and_copy_if_needed(SMAT *smat, LAWORK *lawork)
{
	km_alloc_if_needed(smat, lawork);
	if ( lawork->need_to_free ) {
		km_copy2work(lawork->body, lawork->ld, smat);
	}
}

// km_alloc_if_needed() and 0 clear in any case
void
km_alloc_if_needed_and_0clear(SMAT *smat, LAWORK *lawork)
{
	km_alloc_if_needed(smat, lawork);
	if ( lawork->need_to_free ) { return; }
	if ( smat->vtype == VT_DOUBLE ) {
		for ( size_t i=0; i<smat->m; i++ ) { for ( size_t j=0; j<smat->n; j++ ) {
			lawork->d[i+j*lawork->ld] = 0.0;
		} }
	} else if ( smat->vtype == VT_COMPLEX ) {
		for ( size_t i=0; i<smat->m; i++ ) { for ( size_t j=0; j<smat->n; j++ ) {
			lawork->z[i+j*lawork->ld] = cpack(0.0, 0.0);
		} }
	} else if ( smat->vtype == VT_INT ) {
		for ( size_t i=0; i<smat->m; i++ ) { for ( size_t j=0; j<smat->n; j++ ) {
			lawork->i[i+j*lawork->ld] = 0;
		} }
	} else if ( smat->vtype == VT_BOOL ) {
		for ( size_t i=0; i<smat->m; i++ ) { for ( size_t j=0; j<smat->n; j++ ) {
			lawork->b[i+j*lawork->ld] = false;
		} }
	} else {
		for ( size_t i=0; i<smat->m; i++ ) { for ( size_t j=0; j<smat->n; j++ ) {
			lawork->v[i+j*lawork->ld] = Qfalse;
		} }
	}
}

struct km_c2w_arg {
	union {
		void *work;
		double *dwork;
		COMPLEX *zwork;
		int *iwork;
		bool *bwork;
		VALUE *vwork;
	};
	size_t ld;
};
void
km_c2w_func_d(double *ent, size_t i, size_t j, void *data)
{
	struct km_c2w_arg *arg = (struct km_c2w_arg *)data;
	arg->dwork[i+j*arg->ld] = *ent;
}
void
km_c2w_func_z(COMPLEX *ent, size_t i, size_t j, void *data)
{
	struct km_c2w_arg *arg = (struct km_c2w_arg *)data;
	arg->zwork[i+j*arg->ld] = *ent;
}
void
km_c2w_func_i(int *ent, size_t i, size_t j, void *data)
{
	struct km_c2w_arg *arg = (struct km_c2w_arg *)data;
	arg->iwork[i+j*arg->ld] = *ent;
}
void
km_c2w_func_b(bool *ent, size_t i, size_t j, void *data)
{
	struct km_c2w_arg *arg = (struct km_c2w_arg *)data;
	arg->bwork[i+j*arg->ld] = *ent;
}
void
km_c2w_func_v(VALUE *ent, size_t i, size_t j, void *data)
{
	struct km_c2w_arg *arg = (struct km_c2w_arg *)data;
	arg->vwork[i+j*arg->ld] = *ent;
}
void
km_copy2work(void *work, size_t ldw, SMAT *smat)
{
	if ( smat->stype == ST_FULL && !(smat->trans) && ldw==smat->ld ) {
		memcpy(work, smat->body, km_sizeof_vt(smat->vtype)*LENGTH(smat));
	} else {
		struct km_c2w_arg data = {{work}, ldw};
		if ( smat->vtype == VT_DOUBLE ) {
			if ( smat->stype == ST_RSUB ) {
				km_smat_each_with_index_d(smat, km_c2w_func_d, &data);
			} else if ( smat->trans ) {
				const int one=1, n=s2i(smat->n), ldw_i=s2i(ldw);
				for ( size_t i=0; i<smat->m; i++ ) {
					dcopy_(&n, smat->dbody+(i*smat->ld), &one, data.dwork+i, &ldw_i);
				}
			} else {
				char str_a[] = "A";
				int m=s2i(smat->m), n=s2i(smat->n), ld=s2i(smat->ld), ldw_i=s2i(ldw);
				dlacpy_(str_a, &m, &n, smat->dbody, &ld, work, &ldw_i);
			}
		} else if ( smat->vtype == VT_COMPLEX ) {
			if ( smat->stype == ST_RSUB ) {
				km_smat_each_with_index_z(smat, km_c2w_func_z, &data);
			} else if ( smat->trans ) {
				const int one=1, n=s2i(smat->n), ldw_i=s2i(ldw);
				for ( size_t i=0; i<smat->m; i++ ) {
					zcopy_(&n, smat->zbody+(i*smat->ld), &one, data.zwork+i, &ldw_i);
				}
			} else {
				char str_a[] = "A";
				int m=s2i(smat->m), n=s2i(smat->n), ld=s2i(smat->ld), ldw_i=s2i(ldw);
				zlacpy_(str_a, &m, &n, smat->zbody, &ld, work, &ldw_i);
			}
		} else if ( smat->vtype == VT_INT ) {
			km_smat_each_with_index_i(smat, km_c2w_func_i, &data);
		} else if ( smat->vtype == VT_BOOL ) {
			km_smat_each_with_index_b(smat, km_c2w_func_b, &data);
		} else if ( smat->vtype == VT_VALUE ) {
			km_smat_each_with_index_v(smat, km_c2w_func_v, &data);
		} else {
			rb_raise(km_eInternal, "unknown value type");
		}
	}
}
void
km_cfw_func_d(double *ent, size_t i, size_t j, void *data)
{
	struct km_c2w_arg *arg = (struct km_c2w_arg *)data;
	*ent = arg->dwork[i+j*arg->ld];
}
void
km_cfw_func_z(COMPLEX *ent, size_t i, size_t j, void *data)
{
	struct km_c2w_arg *arg = (struct km_c2w_arg *)data;
	*ent = arg->zwork[i+j*arg->ld];
}
void
km_cfw_func_i(int *ent, size_t i, size_t j, void *data)
{
	struct km_c2w_arg *arg = (struct km_c2w_arg *)data;
	*ent = arg->iwork[i+j*arg->ld];
}
void
km_cfw_func_b(bool *ent, size_t i, size_t j, void *data)
{
	struct km_c2w_arg *arg = (struct km_c2w_arg *)data;
	*ent = arg->bwork[i+j*arg->ld];
}
void
km_cfw_func_v(VALUE *ent, size_t i, size_t j, void *data)
{
	struct km_c2w_arg *arg = (struct km_c2w_arg *)data;
	*ent = arg->vwork[i+j*arg->ld];
}
void
km_copy_from_work(SMAT *smat, void *work, size_t ldw)
{
	if ( smat->stype == ST_FULL && !(smat->trans) && ldw == smat->ld ) {
		memcpy(smat->body, work, km_sizeof_vt(smat->vtype)*LENGTH(smat));
	} else {
		struct km_c2w_arg data = {{work}, ldw};
		if ( smat->vtype == VT_DOUBLE ) {
			if ( smat->stype == ST_RSUB ) {
				km_smat_each_with_index_d(smat, km_cfw_func_d, &data);
			} else if ( smat->trans ) {
				const int one=1, n=s2i(smat->n), ldw_i=s2i(ldw);
				for ( size_t i=0; i<smat->m; i++ ) {
					dcopy_(&n, data.dwork+i, &ldw_i, smat->dbody+(i*smat->ld), &one);
				}
			} else {
				char str_a[] = "A";
				int m=s2i(smat->m), n=s2i(smat->n), ld=s2i(smat->ld), ldw_i=s2i(ldw);
				dlacpy_(str_a, &m, &n, work, &ldw_i, smat->dbody, &ld);
			}
		} else if ( smat->vtype == VT_COMPLEX ) {
			if ( smat->stype == ST_RSUB ) {
				km_smat_each_with_index_z(smat, km_cfw_func_z, &data);
			} else if ( smat->trans ) {
				const int one=1, n=s2i(smat->n), ldw_i=s2i(ldw);
				for ( size_t i=0; i<smat->m; i++ ) {
					zcopy_(&n, data.zwork+i, &ldw_i, smat->zbody+(i*smat->ld), &one);
				}
			} else {
				char str_a[] = "A";
				int m=s2i(smat->m), n=s2i(smat->n), ld=s2i(smat->ld), ldw_i=s2i(ldw);
				zlacpy_(str_a, &m, &n, work, &ldw_i, smat->zbody, &ld);
			}
		} else if ( smat->vtype == VT_INT ) {
			km_smat_each_with_index_i(smat, km_cfw_func_i, &data);
		} else if ( smat->vtype == VT_BOOL ) {
			km_smat_each_with_index_b(smat, km_cfw_func_b, &data);
		} else if ( smat->vtype == VT_VALUE ) {
			km_smat_each_with_index_v(smat, km_cfw_func_v, &data);
		} else {
			rb_raise(km_eInternal, "unknown value type");
		}
	}
}

void
km_free_if_needed(LAWORK *work)
{
	if ( work->need_to_free ) {
		ruby_xfree(work->body);
		work->body = NULL;
		work->need_to_free = false;
	}
}
void
km_copy_and_free_if_needed(SMAT *smat, LAWORK *work)
{
	if ( work->need_to_free ) {
		km_copy_from_work(smat, work->body, work->ld);
		ruby_xfree(work->body);
		work->body = NULL;
		work->need_to_free = false;
	}
}
