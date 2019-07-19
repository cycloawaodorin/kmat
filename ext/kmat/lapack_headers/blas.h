// The original version of this header file has been written by Azalea Clive.
// The download link is available at <http://azalea.s35.xrea.com/blas/sample.html>.

// The following modification has been made by cycloawaodorin.
// * Use complex.h as structs for complex.
// * Add const modifer for appropriate arguments.

#include <complex.h>

#ifndef BLAS_INTERFACE_HEADER
#define BLAS_INTERFACE_HEADER

#ifdef __cplusplus 
extern "C"{
#endif

//Structs
#undef complex
#define complex float _Complex
#define doublecomplex double _Complex

//int xrebla_(char *srname, int *info);

//Level1

//AXPY
void saxpy_(const int *n, const float         *alpha, const float         *x, const int *incx, float         *y, const int *incy);
void daxpy_(const int *n, const double        *alpha, const double        *x, const int *incx, double        *y, const int *incy);
void caxpy_(const int *n, const complex       *alpha, const complex       *x, const int *incx, complex       *y, const int *incy);
void zaxpy_(const int *n, const doublecomplex *alpha, const doublecomplex *x, const int *incx, doublecomplex *y, const int *incy);

//SUM
float   sasum_(const int *n, const float         *x, const int *incx);
float  scasum_(const int *n, const complex       *x, const int *incx);
double  dasum_(const int *n, const double        *x, const int *incx);
double dzasum_(const int *n, const doublecomplex *x, const int *incx);

//COPY
void scopy_(const int *n, const float         *x, const int *incx, float         *y, const int *incy);
void dcopy_(const int *n, const double        *x, const int *incx, double        *y, const int *incy);
void ccopy_(const int *n, const complex       *x, const int *incx, complex       *y, const int *incy);
void zcopy_(const int *n, const doublecomplex *x, const int *incx, doublecomplex *y, const int *incy);

//DOT
float  sdot_(const int *n, const float  *x, const int *incx, const float  *y, const int *incy);
double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);

//DOTC
complex       cdotc_(const int *n, const complex       *x, const int *incx, const complex       *y, const int *incy);
doublecomplex zdotc_(const int *n, const doublecomplex *x, const int *incx, const doublecomplex *y, const int *incy);

//DOTU
complex       cdotu_(const int *n, const complex       *x, const int *incx, const complex       *y, const int *incy);
doublecomplex zdotu_(const int *n, const doublecomplex *x, const int *incx, const doublecomplex *y, const int *incy);

//NRM2
float   snrm2_(const int *n, const float         *x, const int *incx);
double  dnrm2_(const int *n, const double        *x, const int *incx);
float  scnrm2_(const int *n, const complex       *x, const int *incx);
double dznrm2_(const int *n, const doublecomplex *x, const int *incx);

//ROT
void  srot_(const int *n, float         *x, const int *incx, float         *y, const int *incy, const float  *c, const float  *s);
void  drot_(const int *n, double        *x, const int *incx, double        *y, const int *incy, const double *c, const double *s);
void csrot_(const int *n, complex       *x, const int *incx, complex       *y, const int *incy, const float  *c, const float  *s);
void zdrot_(const int *n, doublecomplex *x, const int *incx, doublecomplex *y, const int *incy, const double *c, const double *s);

//ROTG
void srotg_(float         *a,float         *b, float  *c, float  *s);
void drotg_(double        *a,double        *b, double *c, double *s);
void crotg_(complex       *a,complex       *b, float  *c, float  *s);
void zrotg_(doublecomplex *a,doublecomplex *b, double *c, double *s);

//Stub
//ROTMG
//ROTM


//SCAL
void  sscal_(const int *n,  const float         *a, float          *x, const int *incx);
void  dscal_(const int *n,  const double        *a, double         *x, const int *incx);
void  cscal_(const int *n,  const complex       *a, complex        *x, const int *incx);
void  zscal_(const int *n,  const doublecomplex *a, doublecomplex  *x, const int *incx);
void csscal_(const int *n,  const float         *a, complex        *x, const int *incx);
void zdscal_(const int *n,  const double        *a, doublecomplex  *x, const int *incx);

//SWAP
void sswap_(const int *n, float         *x, const int *incx, float         *y, const int *incy);
void dswap_(const int *n, double        *x, const int *incx, double        *y, const int *incy);
void cswap_(const int *n, complex       *x, const int *incx, complex       *y, const int *incy);
void zswap_(const int *n, doublecomplex *x, const int *incx, doublecomplex *y, const int *incy);

//IAMAX
int isamax_(const int *n, const float         *x, const int *incx);
int idamax_(const int *n, const double        *x, const int *incx);
int icamax_(const int *n, const complex       *x, const int *incx);
int izamax_(const int *n, const doublecomplex *x, const int *incx);

//IAMIN
int isamin_(const int *n, const float         *x, const int *incx);
int idamin_(const int *n, const double        *x, const int *incx);
int icamin_(const int *n, const complex       *x, const int *incx);
int izamin_(const int *n, const doublecomplex *x, const int *incx);

//IMAX
int ismax_(const int *n, const float  *x, const int *incx);
int idmax_(const int *n, const double *x, const int *incx);

//IMIN
int ismin_(const int *n, const float  *x, const int *incx);
int idmin_(const int *n, const double *x, const int *incx);

//Level2

//GBMV
void sgbmv_(char *trans, int *m, int *n, int *kl, int *ku, 
            float         *alpha, float         *A, int *ldA, float         *x, int *incx,
            float         *beta , float         *y, int *incy);
void dgbmv_(char *trans, int *m, int *n, int *kl, int *ku, 
            double        *alpha, double        *A, int *ldA, double        *x, int *incx, 
            double        *beta , double        *y, int *incy);
void cgbmv_(char *trans, int *m, int *n, int *kl, int *ku, 
            complex       *alpha, complex       *A, int *ldA, complex       *x, int *incx, 
            complex       *beta , complex       *y, int *incy);
void zgbmv_(char *trans, int *m, int *n, int *kl, int *ku, 
            doublecomplex *alpha, doublecomplex *A, int *ldA, doublecomplex *x, int *incx, 
            doublecomplex *beta , doublecomplex *y, int *incy);

//GEMV
void sgemv_(char *trans, int *m, int *n, 
            float         *alpha, float         *A, int *ldA, float         *x, int *incx,
            float         *beta , float         *y, int *incy);
void dgemv_(char *trans, int *m, int *n, 
            double        *alpha, double        *A, int *ldA, double        *x, int *incx, 
            double        *beta , double        *y, int *incy);
void cgemv_(char *trans, int *m, int *n, 
            complex       *alpha, complex       *A, int *ldA, complex       *x, int *incx, 
            complex       *beta , complex       *y, int *incy);
void zgemv_(char *trans, int *m, int *n, 
            doublecomplex *alpha, doublecomplex *A, int *ldA, doublecomplex *x, int *incx, 
            doublecomplex *beta , doublecomplex *y, int *incy);

//GER
void sger_(int *m, int *n, float  *alpha,  float *x, int *incx,  float *y, int *incy,  float *A, int *ldA);
void dger_(int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *A, int *ldA);

//GERC
void cgerc_(int *m, int *n, complex       *alpha, complex       *x, int *incx,
            complex       *y, int *incy, complex       *A, int *ldA);
void zgerc_(int *m, int *n, doublecomplex *alpha, doublecomplex *x, int *incx, 
            doublecomplex *y, int *incy, doublecomplex *A, int *ldA);

//GREU
void cgeru_(int *m, int *n, complex       *alpha, complex       *x, int *incx,
            complex       *y, int *incy, complex       *A, int *ldA);
void zgeru_(int *m, int *n, doublecomplex *alpha, doublecomplex *x, int *incx, 
            doublecomplex *y, int *incy, doublecomplex *A, int *ldA);

//HBMV
void chbmv_(char *uplo, int *n, int *k, complex       *alpha, complex       *A, int *ldA,
            complex       *x, int *incx, complex       *beta, complex       *y, int *incy);
void zhbmv_(char *uplo, int *n, int *k, doublecomplex *alpha, doublecomplex *A, int *ldA,
            doublecomplex *x, int *incx, doublecomplex *beta, doublecomplex *y, int *incy);

//HEMV
void chemv_(char *uplo, int *n, complex       *alpha, complex       *A, int *ldA,
            complex       *x, int *incx, complex       *beta, complex       *y, int *incy);
void zhemv_(char *uplo, int *n, doublecomplex *alpha, doublecomplex *A, int *ldA,
            doublecomplex *x, int *incx, doublecomplex *beta, doublecomplex *y, int *incy);

//HER
void cher_(char *uplo, int *n, float  *alpha, complex       *x, int *incx, complex       *A, int *ldA);
void zher_(char *uplo, int *n, double *alpha, doublecomplex *x, int *incx, doublecomplex *A, int *ldA);

//Stub
//HER2

//HPMV
void chpmv_(char *uplo, int *n, complex       *alpha, complex       *A,
            complex       *x, int *incx, complex       *beta, complex       *y, int *incy);
void zhpmv_(char *uplo, int *n, doublecomplex *alpha, doublecomplex *A,
            doublecomplex *x, int *incx, doublecomplex *beta, doublecomplex *y, int *incy);

//HPR
void chpr_ (char *uplo, int *n, float   *alpha, complex       *x, int *incx, complex       *A);
void zhpr_ (char *uplo, int *n, double  *alpha, doublecomplex *x, int *incx, doublecomplex *A);

//Stub
//HPR2

//SBMV
void ssbmv_(char *uplo, int *n, int *k, float  *alpha, float  *A, int *ldA,
            float  *x, int *incx, float  *beta, float  *y, int *incy);
void dsbmv_(char *uplo, int *n, int *k, double *alpha, double *A, int *ldA,
            double *x, int *incx, double *beta, double *y, int *incy);

//SPMV
void sspmv_(char *uplo, int *n, float  *alpha, float  *A, float  *x, int *incx, float  *beta, float  *y, int *incy);
void dspmv_(char *uplo, int *n, double *alpha, double *A, double *x, int *incx, double *beta, double *y, int *incy);

//SPR
void sspr_(char *uplo, int *n, float  *alpha, float  *x, int *incx, float  *A);
void dspr_(char *uplo, int *n, double *alpha, double *x, int *incx, double *A);

//Stub
//SPR2

//SYMV
void ssymv_(char *uplo, int *n, float  *alpha, float  *A, int *ldA,
            float  *x, int *incx, float  *beta, float  *y, int *incy);
void dsymv_(char *uplo, int *n, double *alpha, double *A, int *ldA,
            double *x, int *incx, double *beta, double *y, int *incy);

//SYR
void ssyr_(char *uplo, int *n, float  *alpha, float  *x, int *incx, float  *A, int *ldA);
void dsyr_(char *uplo, int *n, double *alpha, double *x, int *incx, double *A, int *ldA);

//Stub
//SYR2

//TBMV
void stbmv_(char *uplo, char *trans, char *diag, int *n, int *k, float         *A, int *ldA, float         *x, int *incx);
void dtbmv_(char *uplo, char *trans, char *diag, int *n, int *k, double        *A, int *ldA, double        *x, int *incx);
void ctbmv_(char *uplo, char *trans, char *diag, int *n, int *k, complex       *A, int *ldA, complex       *x, int *incx);
void ztbmv_(char *uplo, char *trans, char *diag, int *n, int *k, doublecomplex *A, int *ldA, doublecomplex *x, int *incx);

//TBSV
void stbsv_(char *uplo, char *trans, char *diag, int *n, int *k, float         *A, int *ldA, float         *x, int *incx);
void dtbsv_(char *uplo, char *trans, char *diag, int *n, int *k, double        *A, int *ldA, double        *x, int *incx);
void ctbsv_(char *uplo, char *trans, char *diag, int *n, int *k, complex       *A, int *ldA, complex       *x, int *incx);
void ztbsv_(char *uplo, char *trans, char *diag, int *n, int *k, doublecomplex *A, int *ldA, doublecomplex *x, int *incx);

//TPMV
void stpmv_(char *uplo, char *trans, char *diag, int *n, float         *A, float         *x, int *incx);
void dtpmv_(char *uplo, char *trans, char *diag, int *n, double        *A, double        *x, int *incx);
void ctpmv_(char *uplo, char *trans, char *diag, int *n, complex       *A, complex       *x, int *incx);
void ztpmv_(char *uplo, char *trans, char *diag, int *n, doublecomplex *A, doublecomplex *x, int *incx);

//TPSV
void stpsv_(char *uplo, char *trans, char *diag, int *n, float         *A, float         *x, int *incx);
void dtpsv_(char *uplo, char *trans, char *diag, int *n, double        *A, double        *x, int *incx);
void ctpsv_(char *uplo, char *trans, char *diag, int *n, complex       *A, complex       *x, int *incx);
void ztpsv_(char *uplo, char *trans, char *diag, int *n, doublecomplex *A, doublecomplex *x, int *incx);

//TRSV
void strsv_(char *uplo, char *trans, char *diag, int *n, float         *A, int *ldA, float         *x, int *incx);
void dtrsv_(char *uplo, char *trans, char *diag, int *n, double        *A, int *ldA, double        *x, int *incx);
void ctrsv_(char *uplo, char *trans, char *diag, int *n, complex       *A, int *ldA, complex       *x, int *incx);
void ztrsv_(char *uplo, char *trans, char *diag, int *n, doublecomplex *A, int *ldA, doublecomplex *x, int *incx);

//TRMV
void strmv_(char *uplo, char *trans, char *diag, int *n, float         *A, int *ldA, float         *x, int *incx);
void dtrmv_(char *uplo, char *trans, char *diag, int *n, double        *A, int *ldA, double        *x, int *incx);
void ctrmv_(char *uplo, char *trans, char *diag, int *n, complex       *A, int *ldA, complex       *x, int *incx);
void ztrmv_(char *uplo, char *trans, char *diag, int *n, doublecomplex *A, int *ldA, doublecomplex *x, int *incx);

//Level3

//GEMM
void sgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
            const float         *alpha, const float         *A, const int *ldA, const float         *B, const int *ldB, 
            const float         *beta , float         *C, const int *ldC);
void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, 
            const double        *alpha, const double        *A, const int *ldA, const double        *B, const int *ldB, 
            const double        *beta , double        *C, const int *ldC);
void cgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, 
            const complex       *alpha, const complex       *A, const int *ldA, const complex       *B, const int *ldB, 
            const complex       *beta , complex       *C, const int *ldC);
void zgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
            const doublecomplex *alpha, const doublecomplex *A, const int *ldA, const doublecomplex *B, const int *ldB, 
            const doublecomplex *beta , doublecomplex *C, const int *ldC);

//HEMM
void chemm_(char *side, char *uplo, int *m, int *n, complex       *alpha, complex       *A, int *ldA,
            complex       *B, int *ldB, complex       *beta, complex       *C, int *ldC);
void zhemm_(char *side, char *uplo, int *m, int *n, doublecomplex *alpha, doublecomplex *A, int *ldA,
            doublecomplex *B, int *ldB, doublecomplex *beta, doublecomplex *C, int *ldC);

//HERK
void cherk_(char *uplo, char *trans, int *n, int *k, float  *alpha, complex       *A, int *ldA,
            float  *beta , complex       *C, int *ldC);
void zherk_(char *uplo, char *trans, int *n, int *k, double *aplha, doublecomplex *A, int *ldA,
            double *beta , doublecomplex *C, int *ldC);

//HERK2
void cher2k_(char *uplo, char *trans, int *n, int *k, complex       *alpha, complex       *A, int *ldA,
             complex       *B, int *ldB,float  *beta, complex       *C, int *ldC);
void zher2k_(char *uplo, char *trans, int *n, int *k, doublecomplex *alpha, doublecomplex *A, int *ldA,
             doublecomplex *B, int *ldB,double *beta, doublecomplex *C, int *ldC);

//SYMM
void ssymm_(char *side, char *uplo, int *m, int *n, 
            float         *alpha, float         *A, int *ldA, float         *B, int *ldB, 
            float         *beta , float         *C, int *ldC);
void dsymm_(char *side, char *uplo, int *m, int *n, 
            double        *alpha, double        *A, int *ldA, double        *B, int *ldB, 
            double        *beta , double        *C, int *ldC);
void csymm_(char *side, char *uplo, int *m, int *n, 
            complex       *alpha, complex       *A, int *ldA, complex       *B, int *ldB, 
            complex       *beta , complex       *C, int *ldC);
void zsymm_(char *side, char *uplo, int *m, int *n, 
            doublecomplex *alpha, doublecomplex *A, int *ldA, doublecomplex *B, int *ldB, 
            doublecomplex *beta , doublecomplex *C, int *ldC);

//SYRK
void ssyrk_(char *uplo, char *trans, int *n, int *k, float  *alpha, float         *A, int *ldA,
            float  *beta , float         *C, int *ldC);
void dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double        *A, int *ldA,
            double *beta , double        *C, int *ldC);
void csyrk_(char *uplo, char *trans, int *n, int *k, float  *alpha, complex       *A, int *ldA,
            float  *beta , complex       *C, int *ldC);
void zsyrk_(char *uplo, char *trans, int *n, int *k, double *aplha, doublecomplex *A, int *ldA,
            double *beta , doublecomplex *C, int *ldC);

//SYR2K
void ssyr2k_(char *uplo, char *trans, int *n, int *k, float  *alpha, float         *A, int *ldA, float         *B, int *ldB,
            float  *beta , float         *C, int *ldC);
void dsyr2k_(char *uplo, char *trans, int *n, int *k, double *alpha, double        *A, int *ldA, double        *B, int *ldB,
            double *beta , double        *C, int *ldC);
void csyr2k_(char *uplo, char *trans, int *n, int *k, float  *alpha, complex       *A, int *ldA, complex       *B, int *ldB,
            float  *beta , complex       *C, int *ldC);
void zsyr2k_(char *uplo, char *trans, int *n, int *k, double *aplha, doublecomplex *A, int *ldA, doublecomplex *B, int *ldB,
            double *beta , doublecomplex *C, int *ldC);

//TRMM
void strmm_(char *side, char *uplo, char *trans, char *diag, int *m, int *n, 
            float         *alpha, float         *A, int *ldA, float         *B, int *ldB);
void dtrmm_(char *side, char *uplo, char *trans, char *diag, int *m, int *n, 
            double        *alpha, double        *A, int *ldA, double        *B, int *ldB);
void ctrmm_(char *side, char *uplo, char *trans, char *diag, int *m, int *n, 
            complex       *alpha, complex       *A, int *ldA, complex       *B, int *ldB);
void ztrmm_(char *side, char *uplo, char *trans, char *diag, int *m, int *n, 
            doublecomplex *alpha, doublecomplex *A, int *ldA, doublecomplex *B, int *ldB);

//TRSM
void strsm_(char *side, char *uplo, char *trans, char *diag, int *m, int *n, 
            float         *alpha, float         *A, int *ldA, float         *B, int *ldB);
void dtrsm_(char *side, char *uplo, char *trans, char *diag, int *m, int *n, 
            double        *alpha, double        *A, int *ldA, double        *B, int *ldB);
void ctrsm_(char *side, char *uplo, char *trans, char *diag, int *m, int *n, 
            complex       *alpha, complex       *A, int *ldA, complex       *B, int *ldB);
void ztrsm_(char *side, char *uplo, char *trans, char *diag, int *m, int *n, 
            doublecomplex *alpha, doublecomplex *A, int *ldA, doublecomplex *B, int *ldB);

#undef complex
#undef doublecomplex
#define complex _Complex

#ifdef __cplusplus 
}
#endif

#endif

