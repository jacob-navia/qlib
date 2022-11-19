#ifndef __qhead_h__
#define __qhead_h__
#include <string.h>
/* Type of the array elements in a Q-type number */
typedef unsigned long long QELT;
/* Number of bits in a word of that type */
#define WORDSIZE 64
#define MANTISSA_LENGTH	7
#define ACCUM_LENGTH 9
/* Most significant bit of that type.  */
#define SIGNBIT (0x8000000000000000ULL)
/* Largest exponent value */
#define MAXEXP 1048576 // 0x100000
/* The exponent of 1.0 */
#define EXPONE 0x80001
/* Number of bits of precision */
#define NBITS (MANTISSA_LENGTH*64)
/* Bits in accumulators */
#define NBITS_ACC	(ACCUM_LENGTH*64)
/* Maximum number of decimal digits in conversion */
#define NDEC (NBITS*8/27)
#define STDCALL
#define QLIB_OVERFLOW 1
#ifdef x86_64
#define EXTRAWORDS 6
#else
#define EXTRAWORDS 4
#endif

/* Constant definitions for math error conditions */

#define DOMAIN		1	/* argument domain error */
#define SING		2	/* argument singularity */
#define OVERFLOW	3	/* overflow range error */
#define UNDERFLOW	4	/* underflow range error */
#define TLOSS		5	/* total loss of precision */
#define PLOSS		6	/* partial loss of precision */

#define EDOM		33
#define ERANGE		34

/* Complex numeral.  */
typedef struct {
	double r;
	double i;
}cmplx;

#ifdef arm64
typedef struct _float128_t{
	unsigned long long low, high;
} float128_t;
#else 
#ifndef __LCC__
typedef long double float128_t;
#endif
#endif
typedef struct {
	unsigned long long mantissalow:64;
	unsigned long mantissahigh:48;
	unsigned int exponent:15;
	unsigned int sign:1;
}  ld113;
/* Long double complex numeral.  */
/* Comment out if your compiler does not have long double.  */
typedef struct
	{
	float128_t r;
	float128_t i;
	}cmplxl;
#define EXPORT
#define IMPORT

/* If you define UNK, then be sure to set BIGENDIAN properly. */
#define BIGENDIAN 0

/* Define to ask for infinity support, else undefine. */
#define INFINITIES 1 

/* Define to ask for support of numbers that are Not-a-Number,
   else undefine.  This may automatically define INFINITY in some files. */
/* #define NANS 1 */

/* Define to support tiny denormal numbers, else undefine. */
#define DENORMAL 1

/* Define to distinguish between -0.0 and +0.0.  */
#define MINUSZERO 1

/* Define 1 for ANSI C atan2() function
   See atan.c and clog.c. */
#define ANSIC 1

/* Variable for error reporting.  See mtherr.c.  */
extern int merror;

#define ACC_LEN 20
enum { ADDITION,SUBTRACTION};
#if 0
#define __shdn1 shdn1
#define __shup1 shup1
#define __addm addm
#define __subm subm
#define __mulin mulin
#define __mulm mulm
#define __divm divm
#define __qmovz	qmovz
#define __shup1 shup1
#define __shdn1 shdn1
#define __shiftupn shiftupn
#define __shiftdownn shiftdownn
#define _mdnorm mdnorm
#endif

#define signof(a) (a->sign)
#define setpositive(a) (a->sign=0)
#define setnegative(a) (a->sign=-1)
#define setsign(a,b) (a->sign = b)
#define exponent(a) (a->exponent)
#define setexponent(a,e) (a->exponent=e)
#define decreaseExponent(y) (y->exponent -=1)
#define iszero(x) (0==x->exponent)
#pragma pack(push,1)

typedef struct {
        unsigned int mantissa1;
        unsigned int mantissa0:20;
        unsigned int exponent:11;
        unsigned int sign:1;
} _Double;
union ud {
	_Double d;
	double dd;
} ;
typedef struct tagQfloat {
	int sign;
	unsigned int exponent;
	unsigned long long mantissa[MANTISSA_LENGTH];
} Qfloat;

typedef struct tagQfloatAccum {
	int sign;
	unsigned int exponent;
	unsigned long long mantissa[MANTISSA_LENGTH+EXTRAWORDS];
} QfloatAccum;
typedef struct {
        Qfloat n[1]; /* numerator */
        Qfloat d[1]; /* denominator */
        }qfract;
typedef Qfloat *Qfloatp;
typedef QfloatAccum *QfloatAccump;
typedef struct {
        Qfloat re[1];
        Qfloat im[1];
} qcmplx;

typedef union __int128Struct {
	struct  {
		unsigned long long low;
		long long high;
	} U ;
	__int128 int128;
} Int128S;
#pragma pack(pop)
#define NG 67

#include <math.h>
#undef INFINITY
double lgam(double);
double polevl(double x,double coef[],int N );
double p1evl(double x,double coef[],int N );
double p1evll (double x, void *P, int n );
double incbet(double aa,double bb,double xx );
double ndtri(double y0);
int cmpm ( QfloatAccump a, QfloatAccump b);
void mulm(Qfloatp,Qfloatp,QfloatAccump);
int dtoq ( unsigned short *d, Qfloatp y );
void e113toq ( float128_t e, Qfloatp y );
void e24toq ( unsigned short *pe, Qfloatp y );
void e64toq ( unsigned short *e, Qfloatp y );
int etoq (double e, Qfloatp y );
long ltoq ( int *lp, Qfloatp y );
void itoq(long long,Qfloatp);
int mtherr ( char *name, int code );
void qabs ( Qfloatp  x );
void qacos ( Qfloatp const x, Qfloatp y );
void qacosh ( Qfloatp const x, Qfloatp y );
int  qairy ( Qfloatp const x, Qfloatp ai, Qfloatp aip, Qfloatp bi, Qfloatp bip );
void qasin ( Qfloatp const x, Qfloatp y );
void qasinh ( Qfloatp const x, Qfloatp y );
void qatanh ( Qfloatp const x, Qfloatp y );
void qatn ( Qfloatp const x, Qfloatp y );
void qatn2 ( Qfloatp const x, Qfloatp y, Qfloatp z );
void qbdtr (Qfloatp const k, Qfloatp const n, Qfloatp const p, Qfloatp y );
void qbdtrc (Qfloatp const  k, Qfloatp const n, Qfloatp const p, Qfloatp y );
void qbdtri ( int k, int n, Qfloatp const y, Qfloatp p );
void qbeta ( Qfloatp const a, Qfloatp const b, Qfloatp y );
void qcabs (qcmplx *a, Qfloatp y );
void qcacos (qcmplx *z,qcmplx *w );
void qcadd (qcmplx *a,qcmplx *b,qcmplx *c );
void qcasin ( qcmplx *a, qcmplx *w );
void qcatan ( qcmplx *z, qcmplx *w );
void qcbrt ( Qfloatp const xx, Qfloatp y );
void qccos ( qcmplx *a, qcmplx *c );
void qccot ( qcmplx *z, qcmplx *w );
void qcdiv ( qcmplx *a, qcmplx *b, qcmplx *c );
void qcexp ( qcmplx *a, qcmplx *c );
void qcgamma ( qcmplx *x, qcmplx *y );
void qchdtc ( Qfloatp const df, Qfloatp constx, Qfloatp y );
void qChiSquareComp(Qfloatp const df, Qfloatp const y, Qfloatp x );
void qChiSquare( Qfloatp df, Qfloatp x, Qfloatp y );
void qclgam ( qcmplx *const x, qcmplx *y );
void qclog ( qcmplx *a, qcmplx *c );
void qcmov ( qcmplx *a, qcmplx *b );
int  qcmp(Qfloatp const p,const  Qfloatp q );
int qcmpto1(Qfloatp const x);
void qcmul ( qcmplx *a, qcmplx *b, qcmplx *c );
void qcneg ( qcmplx *a );
void qfcos ( Qfloatp const x, Qfloatp y );
void qcosdg ( Qfloatp const x, Qfloatp y );
void qcosh ( Qfloatp const x, Qfloatp y );
void qcosm1 ( Qfloatp x, Qfloatp y );
void qcot ( Qfloatp const x, Qfloatp y );
void qcpow ( qcmplx *x, qcmplx *y, qcmplx *z );
void qcsin ( qcmplx *a, qcmplx *c );
void qcsqrt ( qcmplx *z, qcmplx *w );
void qcsub ( qcmplx *a, qcmplx *b, qcmplx *c );
void qctan ( qcmplx *z, qcmplx *w );
int qcpolylog ( int n, qcmplx *x, qcmplx *y );
void qpolylog ( int n, Qfloatp x, Qfloatp y );
int qdawsn ( Qfloatp const xx, Qfloatp y );
void qdiv (const Qfloatp a, const Qfloatp b, Qfloatp c );
void qinv(Qfloatp const a,Qfloatp b);
void qdivi(long long a, const Qfloatp b, Qfloatp c );
void divi(unsigned long long a,Qfloatp b,QfloatAccump r);
int qellie ( Qfloatp const phi, Qfloatp const m, Qfloatp y );
void qellik ( Qfloatp const phi, const Qfloatp m, Qfloatp y );
void qellpe ( Qfloatp const x, Qfloatp y );
void qellpj (Qfloatp const u,Qfloatp const m,Qfloatp sn,Qfloatp cn,Qfloatp dn,Qfloatp ph);
int qellpk ( Qfloatp const x, Qfloatp y );
void qerf ( Qfloatp const x, Qfloatp y );
void qerfc ( Qfloatp const x, Qfloatp y );
void qeuclid ( Qfloatp const num, Qfloatp const den, Qfloatp gcda );
void qfexp ( Qfloatp const x, Qfloatp y );
void qexp10 ( Qfloatp const x, Qfloatp y );
void qexp2 ( Qfloatp const x, Qfloatp y );
void qexpn ( int n, Qfloatp const x, Qfloatp yy );
int qfac ( Qfloatp const x, Qfloatp y );
void qfdtr ( int ia, int ib, Qfloatp const x, Qfloatp y );
void qfdtrc ( int ia, int ib, Qfloatp const x, Qfloatp y );
void qfdtri ( int ia, int ib, Qfloatp const y, Qfloatp x );
void qfloor (const Qfloatp x, Qfloatp y );
void qfrexp ( Qfloatp const x, long *e, Qfloatp y );
int qfresnl ( Qfloatp const x, Qfloatp const ss, Qfloatp cc );
void qgamma ( Qfloatp const xx, Qfloatp y );
void qgdtr ( Qfloatp const a, Qfloatp const b, Qfloatp const x, Qfloatp y );
void qgdtrc ( Qfloatp const a, Qfloatp const b, Qfloatp const x, Qfloatp y );
void qhy2f1 ( Qfloatp const a, Qfloatp const b, Qfloatp const c, Qfloatp const x, Qfloatp y );
int hypergeomq ( Qfloatp const a, Qfloatp const b, Qfloatp const x, Qfloatp const y );
int qi1 ( Qfloatp x, Qfloatp y );
void qifrac ( Qfloatp x,long long *i, Qfloatp frac );
void qigam ( Qfloatp a, Qfloatp x, Qfloatp y );
void qigamc ( Qfloatp a, Qfloatp x, Qfloatp y );
int qigami ( Qfloatp a, Qfloatp y0, Qfloatp ans );
int bessel_I( Qfloatp n, Qfloatp x, Qfloatp y );
void qincb ( Qfloatp aa, Qfloatp bb, Qfloatp xx, Qfloatp y );
int qincbi ( Qfloatp a, Qfloatp b, Qfloatp yy, Qfloatp ans );
int qincg ( Qfloatp a, Qfloatp x, Qfloatp y );
int qine ( Qfloatp n, Qfloatp x, Qfloatp y );
void qinfin ( Qfloatp x );
int qisneg ( Qfloatp x );
void bessel_J( Qfloatp nn, Qfloatp xx, Qfloatp y );
void qk0 ( Qfloatp x, Qfloatp y );
void bessel_K( Qfloatp const nn, Qfloatp const x, Qfloatp y );
int qkne (Qfloatp nn, Qfloatp const x, Qfloatp y );
void qldexp ( Qfloatp const x, long n, Qfloatp y );
void qlgam ( Qfloatp const x, Qfloatp y );
void qflog ( Qfloatp const x, Qfloatp y );
void qlog1 ( Qfloatp const x, Qfloatp y );
void qlog10 ( Qfloatp const x, Qfloatp y );
void qlogtwo ( Qfloatp const x, Qfloatp y );
int __qmov ( Qfloatp const a, Qfloatp b);
void __swapm(QfloatAccump,QfloatAccump);
void qmul ( Qfloatp const a, Qfloatp const b, Qfloatp c );
int qfma(const Qfloatp a,const Qfloatp b,const Qfloatp addend, Qfloatp result);
void qmuli ( Qfloatp a, Qfloatp b, Qfloatp c );
void qnbdtc ( int k, int n, Qfloatp p, Qfloatp y );
void qnbdtr ( int k, int n, Qfloatp p, Qfloatp y );
void qndtr ( Qfloatp x, Qfloatp y );
void qndtri ( Qfloatp qy0, Qfloatp qx0 );
void qneg ( Qfloatp x );
void qnthroot( Qfloatp x, Qfloatp n, Qfloatp r );
int qPoissonDistribution( int k, Qfloatp m, Qfloatp y );
int qPoissonDistributionInv( int k, Qfloatp m, Qfloatp y );
int qPoissonDistributionComp( int k, Qfloatp y, Qfloatp m );
void qfpow ( Qfloatp x, Qfloatp y, Qfloatp z );
void qpowi ( Qfloatp x, Qfloatp y, Qfloatp z );
void qpsi ( Qfloatp x, Qfloatp y );
void qradd ( qfract *ff1, qfract *ff2, qfract *ff3 );
int qrand ( Qfloatp q );
void qrdiv ( qfract *ff1, qfract *ff2, qfract *ff3 );
void qredpi ( Qfloatp x, Qfloatp y );
void qremain ( Qfloatp a, Qfloatp b, Qfloatp c );
void qrmul ( qfract *ff1, qfract *ff2, qfract *ff3 );
void qround ( Qfloatp x, Qfloatp y );
void qrsub ( qfract *ff1, qfract *ff2, qfract *ff3 );
void qshici ( Qfloatp x, Qfloatp si, Qfloatp ci );
int qsici ( Qfloatp x, Qfloatp si, Qfloatp ci );
void qfsin ( Qfloatp x, Qfloatp y );
void qsindg ( Qfloatp x, Qfloatp y );
void qsinh ( Qfloatp x, Qfloatp y );
int qspenc ( Qfloatp x, Qfloatp y );
void qfsqrt (const Qfloatp x, Qfloatp y );
void qstdtri ( int k, Qfloatp p, Qfloatp t );
void qstudt ( int k, Qfloatp t, Qfloatp y );
void qftan ( Qfloatp x, Qfloatp y );
void qtandg ( Qfloatp x, Qfloatp y );
void qtanh(const Qfloatp x, Qfloatp y );
int qtoasc ( Qfloatp q, char *string, int width, int ndigs, int flags);
int qtod ( Qfloatp x, unsigned short *d );
double qtoe( Qfloatp x, int doround );
double qtoefast(Qfloatp x);
enum {NOROUNDING,DOROUNDING};
float128_t qtoe113 ( Qfloatp x);
int qtoe24 ( Qfloatp x, unsigned short *e );
void qtoe64 ( Qfloatp x, unsigned short *e );
void neumann_N( Qfloatp qn, Qfloatp x, Qfloatp y );
void qzetac ( Qfloatp x, Qfloatp y );
int simq(Qfloat qA[1],Qfloat qB[1],Qfloat qX[],int n,int flag,int IPS[]);
int isinfq(Qfloatp a);
int bsr64(long long);
void shiftupn(QfloatAccump ,int);
void shiftdownn(QfloatAccump ,int);
int asctoq(char * s,Qfloatp  y,char ** );
int qadd_subtract(const Qfloatp ,const Qfloatp , Qfloatp ,int );
void qincr(Qfloatp, Qfloatp);
void pack(QfloatAccump,Qfloatp);
void qhank( Qfloatp, Qfloatp, Qfloatp );
void qcatalan(Qfloatp n,Qfloatp result);
void qrecur( Qfloatp, Qfloatp, Qfloatp, Qfloatp );
void qmovz(Qfloatp, QfloatAccump);
void qclear(Qfloatp x);
void shup1(QfloatAccump);
void shup16(QfloatAccump );
void shdn1(QfloatAccump);
void subm(QfloatAccump ,QfloatAccump );
void addm(QfloatAccump,QfloatAccump);
void qmul(Qfloatp,Qfloatp,Qfloatp);
void qsquare(Qfloatp x,Qfloatp y);
void divm(Qfloatp,Qfloatp, QfloatAccump);
void mulin(QELT y, Qfloatp b,QfloatAccump ac3);
void lltoq(long long lp,Qfloatp y);
void qfact(Qfloatp,Qfloatp);
void e64toq(unsigned short *e,Qfloatp y );
void addm(QfloatAccump,QfloatAccump);
void qclear(Qfloat *);
void qmov(const Qfloat *,Qfloat *);
long long qtoll(Qfloatp q);
__int128 qtoi128(Qfloatp x);
void i128toq(__int128 src,Qfloatp result);
void qagm(Qfloatp,Qfloatp,Qfloatp);
void etoqfast(double,Qfloatp);
void qlstir(Qfloatp,Qfloatp);
int roundAccum(QfloatAccump plarger,int lost,int subflg);
void qdivadd(const QfloatAccump x,Qfloatp y,int doround,int iteration,Qfloatp c);
void qnrmlz(Qfloatp x);
void qmulAccum(const Qfloatp a,const QfloatAccump ac2,Qfloatp c);
void qshi(Qfloatp const x,Qfloatp y);
void qchi(Qfloatp const x,Qfloatp y);
void qci(Qfloatp const x,Qfloatp y);
void qsi(Qfloatp const x,Qfloatp y);
void qhypot(Qfloatp x, Qfloatp y, Qfloatp z);
void qchdti(Qfloatp df, Qfloatp y,Qfloatp x);
void qcatalanConstant(Qfloatp y);
int qkolmogorov(Qfloatp const x,Qfloatp y);
void qPochhammerUp(Qfloatp const x,Qfloatp const n,Qfloatp y);
void qPochhammerDown(Qfloatp const x,Qfloatp const n,Qfloatp y);

#define qadd(a,b,c) qadd_subtract(a,b,c,0)
#define qsub(a,b,c) qadd_subtract(a,b,c,1)
void addbit(QfloatAccump);
#ifndef qequal
int qequal(Qfloatp ,Qfloatp );
#endif
void beta_distribution_invQ(Qfloatp a,Qfloatp b,Qfloatp yy,Qfloatp ans);
int bsr(int);
int isinteger(Qfloatp );

extern Qfloat qone[1];
extern Qfloat qminus_one[1];
extern Qfloat qtwo[1];
extern Qfloat qlog2[1];
extern Qfloat qhalf[1];
extern Qfloat qPi_Div_2[1];
extern Qfloat qpi[1];
extern QfloatAccum qpiAccum[1];
extern Qfloat qthree[1];
extern Qfloat oneThird[1];
extern Qfloat qzero[1];
extern Qfloat qnine[1];
extern Qfloat qexp1[1];
extern Qfloat qlog10c[1];
extern Qfloat oneopi[1];
extern Qfloat hankpp[1];
extern Qfloat hanks[1];
extern Qfloat hankqq[1];
extern Qfloat hankc[1];
extern Qfloat hankcc[1];
extern Qfloat hankzz[1];
extern Qfloat qeul[1];
extern Qfloat qsqrt2[1];
extern Qfloat qinv_sqrt2[1];
extern Qfloat ql10e[1];
extern Qfloat qinv_log10[1];
extern QfloatAccum qinv_log2Accum[1];
extern Qfloat qminusone[1];
extern Qfloat invSqrt2pi[1];
extern Qfloat qepsilon[1];
extern Qfloat qinv_log2[1];
extern Qfloat qinv_pi[1];
extern Qfloat qten[1];
extern Qfloat qfifteen[1];
extern Qfloat qeight[1];
extern Qfloat qminus_three[1];
extern Qfloat C1[1];
extern Qfloat C2[1];
extern int qprint(Qfloatp);
#define INVERSE_FACTORIALS     160
#endif /* __QHEAD_H */
