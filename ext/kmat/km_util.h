static inline size_t
int2size_t(int i)
{
	return (size_t)i;
}

static inline void
km_check_frozen(VALUE obj)
{
	if ( OBJ_FROZEN(obj) ) {
		rb_raise(rb_eFrozenError, "can't modify frozen %s", rb_obj_classname(obj));
	}
}

static inline void
km_infect_frozen(VALUE src, VALUE dest)
{
	if ( OBJ_FROZEN(src) ) {
		rb_obj_freeze(dest);
	}
}

static inline void
km_check_positive(int m, int n)
{
	if ( m < 0 || n < 0 ) {
		rb_raise(km_eDim, "matrix size must not be negative");
	} else if ( ((long long)m)*((long long)n) != (long long)(m*n) ) {
		rb_raise(km_eDim, "matrix length must be within int range");
	}
}

#define TODO rb_raise(km_eNotImp, "comming soon")
#define TF2V(p) ( (p) ? Qtrue : Qfalse )
#define ITSELF(x) (x)
#define SWAP(type, a, b) do { type __temp__swap__=a; a=b; b=__temp__swap__; } while (0)
#define SAME(a, b) ( (a) == (b) )
#define XOR(p, q) ( ( (p)&&(!(q)) ) || ( (!(p))&&(q) ) )

#define KALLOC(val, n) val = ruby_xcalloc(int2size_t(n), sizeof(*(val)))
#define KALLOCc(work, smat) work = km_alloc_and_copy(smat)
#define KALLOCn(work, smat) km_alloc_and_copy_if_needed(smat, &(work))
#define KALLOCz(work, smat) km_alloc_if_needed_and_0clear(smat, &(work))

#define INDEX(smat, i, j) ( (smat)->trans ? ((j)+(i)*(smat)->ld) : ((i)+(j)*(smat)->ld) )
#define ENTITYd0(smat, id, idx) ( ((smat)->id##body)[idx] )
#define ENTITYr0(smat, id, idx) ( *( ((smat)->id##pbody)[idx] ) )
#define ENTITY(smat, id, i, j) ( ( (smat)->stype==ST_RSUB ) ? ENTITYr0(smat, id, INDEX(smat, i, j)) : ENTITYd0(smat, id, INDEX(smat, i, j)) )

#define LENGTH(smat) ( (smat)->m*(smat)->n )
#define LENGTHs(smat) int2size_t( (smat)->m*(smat)->n )
#define VECTOR_P(smat) ( (smat)->m==1 || (smat)->n==1 )
#define SAME_SIZE(sa, sb) ( ( (sa)->m==(sb)->m ) && ( (sa)->n==(sb)->n ) )
#define CHECK_SAME_SIZE(sa, sb) do { if ( !SAME_SIZE(sa, sb) ) { rb_raise(km_eDim, "sizes must be the same, (%d, %d) != (%d, %d)", (sa)->m, (sa)->n, (sb)->m, (sb)->n); } } while (0)

#define VT_SWITCH(vt_, dstate, zstate, istate, bstate, vstate) do { \
	VTYPE __vt__ = (vt_); \
	if ( __vt__ == VT_DOUBLE ) { dstate } \
	else if ( __vt__ == VT_COMPLEX ) { zstate } \
	else if ( __vt__ == VT_INT ) { istate } \
	else if ( __vt__ == VT_BOOL ) { bstate } \
	else if ( __vt__ == VT_VALUE ) { vstate } \
	else { rb_raise(km_eInternal, "unknown value type"); } \
} while (0)

// definition is macro within smat/smat.c
void km_smat_each_d(SMAT *smat, void (*func)(double *, void *), void *data);
void km_smat_each_z(SMAT *smat, void (*func)(COMPLEX *, void *), void *data);
void km_smat_each_i(SMAT *smat, void (*func)(int *, void *), void *data);
void km_smat_each_b(SMAT *smat, void (*func)(bool *, void *), void *data);
void km_smat_each_v(SMAT *smat, void (*func)(VALUE *, void *), void *data);
void km_smat_each_with_index_d(SMAT *smat, void (*func)(double *, int, int, void *), void *data);
void km_smat_each_with_index_z(SMAT *smat, void (*func)(COMPLEX *, int, int, void *), void *data);
void km_smat_each_with_index_i(SMAT *smat, void (*func)(int *, int, int, void *), void *data);
void km_smat_each_with_index_b(SMAT *smat, void (*func)(bool *, int, int, void *), void *data);
void km_smat_each_with_index_v(SMAT *smat, void (*func)(VALUE *, int, int, void *), void *data);
void km_smat_each2_d(SMAT *sa, SMAT *sb, void (*func)(double *, double *, void *), void *data);
void km_smat_each2_z(SMAT *sa, SMAT *sb, void (*func)(COMPLEX *, COMPLEX *, void *), void *data);
void km_smat_each2_i(SMAT *sa, SMAT *sb, void (*func)(int *, int *, void *), void *data);
void km_smat_each2_b(SMAT *sa, SMAT *sb, void (*func)(bool *, bool *, void *), void *data);
void km_smat_each2_v(SMAT *sa, SMAT *sb, void (*func)(VALUE *, VALUE *, void *), void *data);
void km_smat_each2_dcd(SMAT *sa, const SMAT *sb, void (*func)(double *, const double *, void *), void *data);
void km_smat_each2_zcz(SMAT *sa, const SMAT *sb, void (*func)(COMPLEX *, const COMPLEX *, void *), void *data);
void km_smat_each2_ici(SMAT *sa, const SMAT *sb, void (*func)(int *, const int *, void *), void *data);
void km_smat_each2_bcb(SMAT *sa, const SMAT *sb, void (*func)(bool *, const bool *, void *), void *data);
void km_smat_each2_vcv(SMAT *sa, const SMAT *sb, void (*func)(VALUE *, const VALUE *, void *), void *data);
void km_smat_each2_dcb(SMAT *sa, const SMAT *sb, void (*func)(double *, const bool *, void *), void *data);
void km_smat_each2_zcb(SMAT *sa, const SMAT *sb, void (*func)(COMPLEX *, const bool *, void *), void *data);
void km_smat_each2_icb(SMAT *sa, const SMAT *sb, void (*func)(int *, const bool *, void *), void *data);
void km_smat_each2_vcb(SMAT *sa, const SMAT *sb, void (*func)(VALUE *, const bool *, void *), void *data);
void km_smat_each2_dcz(SMAT *sa, const SMAT *sb, void (*func)(double *, const COMPLEX *, void *), void *data);
void km_smat_each3_zcdcd(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(COMPLEX *, const double *, const double *, void *), void *data);
void km_smat_each3_bcdcd(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(bool *, const double *, const double *, void *), void *data);
void km_smat_each3_bczcz(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(bool *, const COMPLEX *, const COMPLEX *, void *), void *data);
void km_smat_each3_bcici(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(bool *, const int *, const int *, void *), void *data);
void km_smat_each3_bcbcb(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(bool *, const bool *, const bool *, void *), void *data);
void km_smat_each3_bcvcv(SMAT *sa, const SMAT *sb, const SMAT *sc, void (*func)(bool *, const VALUE *, const VALUE *, void *), void *data);
