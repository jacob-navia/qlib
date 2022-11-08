#ifndef __QHEAD_H
#define __QHEAD_H
#include <string.h>
#define NO_ASM
/* Type of the array elements in a Q-type number */
typedef unsigned int QELT;
/* Number of bits in a word of that type */
#define WORDSIZE 32
/* Most significant bit of that type.  */
#define SIGNBIT (0x80000000)
/* Largest exponent value */
#define MAXEXP 1048576
//#define MAXEXP 65536
/* The exponent of 1.0 */
#define EXPONE 0x80001
/* Number of WORDSIZE-bit words in a q type number (12 or 24) */
#define NQ 16
/* Number of bits of precision */
#define NBITS (7*64)
#define NBITS_ACC	(9*64)
/* Maximum number of decimal digits in conversion */
#define NDEC (NBITS*8/27)
#ifdef _WIN64
#define STDCALL
#else
#define STDCALL _stdcall
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
typedef struct
	{
	double r;
	double i;
	}cmplx;

/* Long double complex numeral.  */
/* Comment out if your compiler does not have long double.  */
#if 1
typedef struct
	{
	long double r;
	long double i;
	}cmplxl;
#endif
#ifdef DLL
#define EXPORT __declspec(dllexport)
#define IMPORT __declspec(dllimport)
#else
#define EXPORT
#define IMPORT
#endif
#define IBMPC 1

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

#include "defines.h"
#define signof(a) (a->sign)
#define setpositive(a) (a->sign=0)
#define setnegative(a) (a->sign=-1)
#define setsign(a,b) (a->sign = b)
#define exponent(a) (a->exponent)
#define setexponent(a,e) (a->exponent=e)
#define decreaseExponent(y) (y->exponent -=1)
#define iszero(x) (0==x->exponent)
#pragma pack(push,1)
#define MANTISSA_LENGTH	7
typedef struct tagQfloat {
	int sign;
	unsigned int exponent;
	unsigned long long mantissa[MANTISSA_LENGTH];
} Qfloat;

typedef struct tagQfloatAccum {
	int sign;
	unsigned int exponent;
	unsigned long long mantissa[MANTISSA_LENGTH+2];
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
#pragma pack(pop)
#define NG 67

int cmpm ( QfloatAccump a, QfloatAccump b );
int __divm ( QfloatAccump a, QfloatAccump b );
void __mulin(Qfloatp,QfloatAccump);
void __mulm(QfloatAccump,QfloatAccump);
int dtoq ( unsigned short *d, Qfloatp y );
void e113toq ( unsigned short *e, Qfloatp y );
void e24toq ( unsigned short *pe, Qfloatp y );
void __e64toq ( unsigned short *e, Qfloatp y );
int etoq ( unsigned short *e, Qfloatp y );
long ltoq ( int *lp, Qfloatp y );
void itoq(int,Qfloatp);
int mtherr ( char *name, int code );
void qabs ( Qfloatp x );
void qacos ( Qfloatp x, Qfloatp y );
void qacosh ( Qfloatp x, Qfloatp y );
void qincr(Qfloatp q,Qfloatp c,int subflg);
int qairy ( Qfloatp x, Qfloatp ai, Qfloatp aip, Qfloatp bi, Qfloatp bip );
void qasin ( Qfloatp x, Qfloatp y );
void qasinh ( Qfloatp x, Qfloatp y );
void qatanh ( Qfloatp x, Qfloatp y );
void qatn ( Qfloatp x, Qfloatp y );
void qatn2 ( Qfloatp x, Qfloatp y, Qfloatp z );
void qbdtr ( int k, int n, Qfloatp p, Qfloatp y );
void qbdtrc ( int k, int n, Qfloatp p, Qfloatp y );
void qbdtri ( int k, int n, Qfloatp y, Qfloatp p );
void qbeta ( Qfloatp a, Qfloatp b, Qfloatp y );
void qcabs ( qcmplx *a, Qfloatp y );
void qcacos ( qcmplx *z, qcmplx *w );
void qcadd ( qcmplx *a, qcmplx *b, qcmplx *c );
void qcasin ( qcmplx *a, qcmplx *w );
void qcatan ( qcmplx *z, qcmplx *w );
void qcbrt ( Qfloatp xx, Qfloatp y );
void qccos ( qcmplx *a, qcmplx *c );
void qccot ( qcmplx *z, qcmplx *w );
void qcdiv ( qcmplx *a, qcmplx *b, qcmplx *c );
void qcexp ( qcmplx *a, qcmplx *c );
void qcgamma ( qcmplx *x, qcmplx *y );
void qchdtc ( Qfloatp df, Qfloatp x, Qfloatp y );
void qchdti ( Qfloatp df, Qfloatp y, Qfloatp x );
void qchdtr ( Qfloatp df, Qfloatp x, Qfloatp y );
void qclgam ( qcmplx *x, qcmplx *y );
void qclog ( qcmplx *a, qcmplx *c );
void qcmov ( qcmplx *a, qcmplx *b );
int qcmp ( Qfloatp p, Qfloatp q );
void qcmul ( qcmplx *a, qcmplx *b, qcmplx *c );
void qcneg ( qcmplx *a );
void qfcos ( Qfloatp x, Qfloatp y );
void qcosdg ( Qfloatp x, Qfloatp y );
void qcosh ( Qfloatp x, Qfloatp y );
int qcosm1 ( Qfloatp x, Qfloatp y );
void qcot ( Qfloatp x, Qfloatp y );
void qcpow ( qcmplx *x, qcmplx *y, qcmplx *z );
void qcsin ( qcmplx *a, qcmplx *c );
void qcsqrt ( qcmplx *z, qcmplx *w );
void qcsub ( qcmplx *a, qcmplx *b, qcmplx *c );
void qctan ( qcmplx *z, qcmplx *w );
int qcpolylog ( int n, qcmplx *x, qcmplx *y );
void qpolylog ( int n, Qfloatp x, Qfloatp y );
int qdawsn ( Qfloatp xx, Qfloatp y );
void qdiv ( Qfloatp a, Qfloatp b, Qfloatp c );
int qellie ( Qfloatp phi, Qfloatp m, Qfloatp y );
void qellik ( Qfloatp phi, Qfloatp m, Qfloatp y );
void qellpe ( Qfloatp x, Qfloatp y );
void qellpj ( Qfloatp u, Qfloatp m, Qfloatp sn, Qfloatp cn, Qfloatp dn, Qfloatp ph );
int qellpk ( Qfloatp x, Qfloatp y );
void qerf ( Qfloatp x, Qfloatp y );
void qerfc ( Qfloatp x, Qfloatp y );
void qeuclid ( Qfloatp num, Qfloatp den, Qfloatp gcda );
void qfexp ( Qfloatp x, Qfloatp y );
void qexp10 ( Qfloatp x, Qfloatp y );
void qexp2 ( Qfloatp x, Qfloatp y );
void qexpn ( int n, Qfloatp x, Qfloatp yy );
int qfac ( Qfloatp x, Qfloatp y );
void qfdtr ( int ia, int ib, Qfloatp x, Qfloatp y );
void qfdtrc ( int ia, int ib, Qfloatp x, Qfloatp y );
void qfdtri ( int ia, int ib, Qfloatp y, Qfloatp x );
void qfloor ( Qfloatp x, Qfloatp y );
void qfrexp ( Qfloatp x, long *e, Qfloatp y );
int qfresnl ( Qfloatp x, Qfloatp ss, Qfloatp cc );
void qgamma ( Qfloatp xx, Qfloatp y );
void qgdtr ( Qfloatp a, Qfloatp b, Qfloatp x, Qfloatp y );
void qgdtrc ( Qfloatp a, Qfloatp b, Qfloatp x, Qfloatp y );
void qhy2f1 ( Qfloatp a, Qfloatp b, Qfloatp c, Qfloatp x, Qfloatp y );
int hypergeomq ( Qfloatp a, Qfloatp b, Qfloatp x, Qfloatp y );
int qi1 ( Qfloatp x, Qfloatp y );
void qifrac ( Qfloatp x, long long *i, Qfloatp frac );
void qigam ( Qfloatp a, Qfloatp x, Qfloatp y );
void qigamc ( Qfloatp a, Qfloatp x, Qfloatp y );
int qigami ( Qfloatp a, Qfloatp y0, Qfloatp ans );
int qin ( Qfloatp n, Qfloatp x, Qfloatp y );
void qincb ( Qfloatp aa, Qfloatp bb, Qfloatp xx, Qfloatp y );
int qincbi ( Qfloatp a, Qfloatp b, Qfloatp yy, Qfloatp ans );
int qincg ( Qfloatp a, Qfloatp x, Qfloatp y );
int qine ( Qfloatp n, Qfloatp x, Qfloatp y );
void qinfin ( Qfloatp x );
int qisneg ( Qfloatp x );
void qjn ( Qfloatp nn, Qfloatp xx, Qfloatp y );
void qk0 ( Qfloatp x, Qfloatp y );
void qkn ( int nn, Qfloatp x, Qfloatp y );
int qkne ( int nn, Qfloatp x, Qfloatp y );
void qldexp ( Qfloatp x, long n, Qfloatp y );
void qlgam ( Qfloatp x, Qfloatp y );
void qflog ( Qfloatp x, Qfloatp y );
void qlog1 ( Qfloatp x, Qfloatp y );
void qlog10 ( Qfloatp x, Qfloatp y );
void qlogtwo ( Qfloatp x, Qfloatp y );
int __qmov ( Qfloatp a, Qfloatp b);
void qmul ( Qfloatp a, Qfloatp b, Qfloatp c );
void qmuli ( Qfloatp a, Qfloatp b, Qfloatp c );
void qnbdtc ( int k, int n, Qfloatp p, Qfloatp y );
void qnbdtr ( int k, int n, Qfloatp p, Qfloatp y );
void qndtr ( Qfloatp x, Qfloatp y );
void qndtri ( Qfloatp qy0, Qfloatp qx0 );
void qneg ( Qfloatp x );
int qpdtr ( int k, Qfloatp m, Qfloatp y );
int qpdtrc ( int k, Qfloatp m, Qfloatp y );
int qpdtri ( int k, Qfloatp y, Qfloatp m );
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
void qfsqrt ( Qfloatp x, Qfloatp y );
void qstdtri ( int k, Qfloatp p, Qfloatp t );
void qstudt ( int k, Qfloatp t, Qfloatp y );
void qftan ( Qfloatp x, Qfloatp y );
void qtandg ( Qfloatp x, Qfloatp y );
void qtanh ( Qfloatp x, Qfloatp y );
int qtoasc ( Qfloatp q, char *string, int ndigs );
int qtod ( Qfloatp x, unsigned short *d );
int qtoe ( Qfloatp x, unsigned short *e );
void qtoe113 ( Qfloatp x, unsigned short *e );
int qtoe24 ( Qfloatp x, unsigned short *e );
void qtoe64 ( Qfloatp x, unsigned short *e );
void qyn ( Qfloatp qn, Qfloatp x, Qfloatp y );
void qzetac ( Qfloatp x, Qfloatp y );
int simq ( Qfloatp A, Qfloatp B, Qfloatp X, int n, int flag, int IPS );
int isinfq(Qfloatp a);
int bsr64(long long);
void shiftupn(QfloatAccump ,int);
void shiftdownn(QfloatAccump ,int);
int asctoq(char * s,Qfloatp  y,char ** );
void qadd_subtract( Qfloatp , Qfloatp , Qfloatp ,int );
void __pack(QfloatAccump,Qfloatp);
void qhank( Qfloatp, Qfloatp, Qfloatp );
void qrecur( Qfloatp, Qfloatp, Qfloatp, Qfloatp );
void qmovz(Qfloatp, QfloatAccump);
void __qclear(Qfloatp x);
void shup1(QfloatAccump);
void shdn1(QfloatAccump);
void subm(QfloatAccump,QfloatAccump);
void addm(QfloatAccump,QfloatAccump);
void __divi(Qfloatp,QfloatAccump);
void __qmul(Qfloatp,Qfloatp,Qfloatp);
void lltoq(long long *lp,Qfloatp y);
int normlz(QfloatAccump,int *);
void qfact(Qfloatp,Qfloatp);
//int shift(QfloatAccump,int *);
void e64toq(unsigned short *e,Qfloatp y );
void __addm(QfloatAccump,QfloatAccump);
#define qadd(a,b,c) __qadd_subtract(a,b,c,0)
#define qsub(a,b,c) __qadd_subtract(a,b,c,1)
void __addbit(QfloatAccump);
#ifndef qequal
int qequal(Qfloatp ,Qfloatp );
#endif
void beta_distribution_invQ(Qfloatp a,Qfloatp b,Qfloatp yy,Qfloatp ans);
int bsr(int);
int isinteger(Qfloatp );
long long qtoll(Qfloatp x);
extern Qfloat qone[1],qtwo[1],qlog2[1],qhalf[1],qPi_Div_2[1],qpi[1],qthree[1],
              qzero[1],qnine[1],qexp1[1],qlog10c[1],oneopi[1],hankpp[1],
              hanks[1],hankqq[1],hankc[1],hankcc[1],hankzz[1],qeul[1],
              qsqrt2[1],ql10e[1],qminusone[1],invSqrt2pi[1],qfive[1];
#endif /* __QHEAD_H */
