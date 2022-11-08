//							qfloat.h
//
#include "qfloat.h"
#include <string.h>
//
extern void qfrand(qfloat *);
extern void qflog(qfloat &,qfloat &);
extern void qlog10(qfloat &, qfloat &);
extern void qcbrt(qfloat &,qfloat&);
extern void qlogb(qfloat &,qfloat &);
extern long long qtoll(qfloat &);
extern void qacosh(qfloat &,qfloat&);
extern void qasinh(qfloat &,qfloat &);
extern void qatanh(qfloat &,qfloat &);
extern void qsinh(qfloat &,qfloat &);
extern void qtanh(qfloat &,qfloat &);
extern void qcosh(qfloat &,qfloat &);
extern void qsinh(qfloat &,qfloat &);
extern void qfrexp(qfloat &,int *,qfloat &);
extern void qldexp(qfloat &,long n,qfloat &);
extern void qexpm1(qfloat &,qfloat &);
extern void qlog1(qfloat &,qfloat &);
extern void qatn(qfloat &,qfloat &);
extern void qatn2(qfloat &,qfloat &,qfloat &);
extern void qerf(qfloat &,qfloat &);
extern void qigam(qfloat &,qfloat &,qfloat&);
extern void qerfc(qfloat &,qfloat &);
extern void qlgam(qfloat &,qfloat &);
extern void qgamma(qfloat &,qfloat &);
extern void qround(qfloat &,qfloat &);
extern void qasin(qfloat &,qfloat &);
extern void qacos(qfloat &,qfloat &);
extern void qlogtwo(qfloat &,qfloat &);
extern void qmov(qfloat &,qfloat &);
extern void qfsqrt(qfloat &,qfloat &);
extern void qexp10(qfloat &,qfloat &);
extern int qremquo(qfloat &,qfloat &,qfloat &);
extern void qhyp(qfloat &,qfloat &,qfloat &x,qfloat &);
extern void qpsi(qfloat &,qfloat &);
extern void qellpk(qfloat &,qfloat &);
extern void qhy2f1(qfloat &a,qfloat &b,qfloat &c,qfloat & x,qfloat & y);
extern void qstudt(int,qfloat &,qfloat &);
extern void qstdtri(int,qfloat &,qfloat &);
extern qfloat qminusone,ql10e;
extern void qincb(qfloat &a,qfloat &b,qfloat &x,qfloat &result);
extern void qbeta(qfloat &a,qfloat &b,qfloat &result);
extern void qspenc(qfloat &,qfloat &);
extern void qexpn(int n,qfloat &x,qfloat &r);
extern void beta_distribution_invQ(qfloat &a,qfloat &b,qfloat &p,qfloat &ans);
extern void qpolylog(int n,qfloat &x,qfloat &y);
extern void qzetac(qfloat &x,qfloat &y);
#ifndef _NQ_
#define _NQ_	12
#endif
#define EXPONENT_BIAS 0x8001
#define WORDSIZE 32
#define NBITS ((_NQ_-1)*WORDSIZE)
int EXPORT qtoi (qfloat &x)
{
	qfloat y;
	long l;
	qifrac(&x, &l, &y);
	return (int) l;
}

int EXPORT qtol (qfloat &x) {
	qfloat y;
	long l;
	qifrac(&x, &l, &y);
	return l;
}

EXPORT float qtof (qfloat &x) {
	float f;
	qtoe24(&x, &f);
	return f;
}

double qtod (qfloat &x) {
	double d;
	qtoe(&x, &d);
	return d;
}

long double qtold (qfloat &x) {
	long double d;
	qtoe64(&x, &d);
	return d;
}


// Arithmetic operators.

qfloat & EXPORT operator += (qfloat & x, const qfloat & y)
{
        qadd (x, y, x);
        return x;
}
qfloat & EXPORT operator +=(qfloat &x,const long double y)
{
        qfloat z;
        e64toq(&y,&z);
        qadd(x,z,x);
        return x;
}
long double & EXPORT operator +=(long double &x,qfloat &y)
{
	x += qtold(y);
	return x;
}

qfloat & EXPORT operator +=(qfloat &x,double y)
{
	qfloat z;
	etoq(&y,&z);
	qadd(x,z,x);
	return x;
}

qfloat EXPORT __declspec(naked) operator + (const qfloat & x, const qfloat & y)
{
}
void __qadd(qfloat &z,qfloat &x,qfloat &y)
{
	qadd (x, y, z);
}

qfloat EXPORT __declspec(naked) operator - (const qfloat & x, const qfloat & y)
{
}

void EXPORT __qsub(qfloat &z,const qfloat & x, const qfloat & y)
{
	qsub (y, x, z);
}

qfloat EXPORT __declspec(naked) operator++(qfloat &x)
{
}
void __qincr(qfloat &z,qfloat &x)
{
	qmov(x,z);
	qadd(x,qone,x);
}

qfloat EXPORT __declspec(naked) operator--(qfloat &x)
{
}
void __qdecr(qfloat &z,qfloat &x)
{
	qmov(x,z);
	qsub(x,qone,x);
}

// Unary negation.
qfloat EXPORT __declspec(naked) operator - (const qfloat & x)
{
}
void __qunaryminus(qfloat &z,qfloat &x)
{
	qmov(x,z);
	qneg (z);
}

qfloat EXPORT __declspec(naked) operator -= (qfloat & x, const qfloat & y)
{
}
void __minusasgn(qfloat &r,qfloat & x, const qfloat & y)
{
	qsub (y, x, x);
	qmov(x,r);
}

qfloat EXPORT __declspec(naked) operator * (const qfloat & x, const qfloat & y)
{
}
void ___qmul(qfloat &z,const qfloat &x,const qfloat &y)
{
      qmul(x, y, z);
}

qfloat EXPORT __declspec(naked) operator *= (qfloat & x, const qfloat & y)
{
}
void __masgn(qfloat &r,qfloat & x, const qfloat & y)
{
	qmul (y, x, x);
	qmov(x,r);
}

qfloat EXPORT __declspec(naked) operator / (const qfloat & x, const qfloat & y)
{
}
void ___qdiv(qfloat & z,qfloat &x,qfloat & y)
{
	qdiv (y, x, z);
}

qfloat EXPORT operator /= (qfloat & x, const qfloat & y)
{
	qdiv (y, x, x);
	return x;
}

// Comparisons.

int EXPORT operator == (const qfloat & x, const qfloat & y)
{
	return (qcmp (x, y) == 0);
}
int EXPORT operator == (const qfloat & x, double y)
{
	qfloat z;
	etoq(&y,&z);
	return (qcmp (x, z) == 0);
}
int EXPORT operator ==(double y,const qfloat &x)
{
	qfloat z;
	etoq(&y,&z);
	return qcmp(x,z) == 0;
}
int EXPORT operator == (const qfloat & x, int y)
{
	qfloat z;
	itoq(y,&z);
	return (qcmp (x, z) == 0);
}

int EXPORT operator!= (const qfloat & x, const qfloat & y)
{
	return (qcmp (x, y) != 0);
}
int EXPORT operator != (const qfloat & x, long double y)
{
	qfloat z;
	if (y == 0.0)
		return qcmp(x,qzero);
	else if (y == 1.0)
		return qcmp(x,qone);
	e64toq(&y,&z);
	return (qcmp (x, z));
}
int EXPORT operator != (const qfloat & x, double y)
{
	qfloat z;
	if (y == 0.0)
		return qcmp(x,qzero);
	else if (x == 1.0)
		return qcmp(x,qone);
	etoq(&y,&z);
	return (qcmp (x, z));
}
int EXPORT operator!=(double y,const qfloat &x)
{
	qfloat z;
	if (y == 0.0)
		return qcmp(x,qzero);
	etoq(&y,&z);
	return qcmp(x,z);
}
int EXPORT operator != (const qfloat & x, int y)
{
	qfloat z;

	if (y == 0)
		return qcmp(x,qzero);
	itoq(y,&z);
	return (qcmp (x, z));
}

int EXPORT operator < (const qfloat & x, const qfloat & y)
{
	return (qcmp (x, y) == -1);
}
int EXPORT operator < (const qfloat &a,const double b)
{
	qfloat x;

	if (b == 0.0)
		return qcmp(a,qzero) == -1;
	etoq(&b, &x);
	return qcmp(a,x) == -1;
}

int EXPORT operator <(const int x, const qfloat &y)
{
	qfloat z;

	if (x == 0)
		return qcmp(qzero,y) == -1;
	itoq(x,&z);
	return qcmp(z,y) == -1;
}

int EXPORT operator <(const qfloat &x,const int y)
{
	qfloat z;

	if (y == 0)
		return qcmp(x,qzero) == -1;
	itoq(y,&z);
	return qcmp(x,z) == -1;
}

int EXPORT operator <(const double x,const qfloat &y)
{
	qfloat z;

	if (x == 0.0)
		return qcmp(qzero,y) == -1;
	etoq(&x,&z);
	return qcmp(z,y) == -1;
}

int EXPORT operator > (const qfloat & x, const qfloat & y)
{
	return (qcmp (x, y) >0);
}

int EXPORT operator >(const qfloat &x,int y)
{
	qfloat z;

	itoq(y,&z);
	return qcmp(x,z) == 1;
}
int EXPORT operator >(const qfloat &x,double y)
{
	qfloat z;

	etoq(&y,&z);
	return qcmp(x,z) == 1;
}

int EXPORT operator >(const qfloat &x,long double y)
{
	qfloat z;
	e64toq(&y,&z);
	return qcmp(x,z) == 1;
}

int EXPORT operator >(const double x,const qfloat &y)
{
	qfloat z;

	etoq(&x,&z);
	return qcmp(z,y) == 1;
}
int EXPORT operator >(const int x,const qfloat &y)
{
	qfloat z;

	itoq(x,&z);
	return qcmp(z,y) == 1;
}

int EXPORT operator >= (const qfloat & x, const qfloat & y)
{
	return (qcmp (x, y) >= 0);
}
int EXPORT operator >=(const qfloat &x,int y)
{
	qfloat z;
	itoq(y,&z);
	return qcmp(x,z) >= 0;
}

int EXPORT operator <= (const qfloat & x, const qfloat & y)
{
	return (qcmp (x, y) <= 0);
}

qfloat EXPORT operator=(qfloat &a,int b)
{
	if (b == 0)
		qclear(a);
	else
		itoq(b,&a);
	return a;
}

qfloat EXPORT operator=(qfloat &a,double b)
{
	if (b == 0.0)
		qclear(a);
	else
		etoq(&b,&a);
	return a;
}
qfloat EXPORT operator=(qfloat &a,float b)
{
	double c = b;
	if (c == 0.0)
		qclear(a);
	else
		etoq(&c,&a);
	return a;
}
qfloat EXPORT operator=(qfloat &a,long double b)
{
	if (b == 0.0)
		qclear(a);
	else
		e64toq(&b,&a);
	return a;
}

qfloat EXPORT operator=(qfloat &a,char *b)
{
	asctoq(b,&a);
	return a;
}

double EXPORT operator=(double &a,qfloat &b)
{
	a = qtod(b);
	return a;
}

int EXPORT operator=(int &a,qfloat &b)
{
	qfloat z;
	a = qtod(b);
	return a;
}

long double EXPORT operator()(qfloat &x)
{
        long double d;
        qtoe64(&x, &d);
        return d;
}

qfloat EXPORT __declspec(naked) operator ()(int n)
{
}
void __itoq(qfloat &z,int n)
{
	itoq(n,&z);
}

qfloat EXPORT __declspec(naked) operator ()(unsigned long n)
{
}
void __longtoq(qfloat &z,unsigned long n)
{
	long long d=n;

	lltoq(&d,z);
}

qfloat EXPORT __declspec(naked) operator ()(double n)
{
}
void __doubletoq(qfloat *z,double n)
{
	etoq(&n,z);
}

qfloat EXPORT __declspec(naked) operator ()(long double n)
{
}
void __longdoubletoq(qfloat *z,long double n)
{
	e64toq(&n,z);
}

qfloat EXPORT __declspec(naked) operator ()(long long n)
{
}
void __longlongtoq(qfloat &z,long long n)
{
	lltoq(&n,z);
}

qfloat EXPORT __declspec(naked) operator()(unsigned int i)
{
}
void __uinttoq(qfloat &z,unsigned int i)
{
	long long ll = i;
	lltoq(&ll,z);
}

int EXPORT operator()(qfloat &x)
{
        qfloat y;
        long l;
        qifrac(&x, &l, &y);
        return l;
}

unsigned int EXPORT operator()(qfloat &x)
{
	long long ll;

	ll = qtoll(x);
	return ll;
}


double EXPORT operator()(qfloat &x)
{
        double d;
        qtoe(&x, &d);
        return d;
}

long long EXPORT operator()(qfloat &x)
{
	return qtoll(x);
}


// Exact remainder of x/y.
qfloat EXPORT __declspec(naked) fmodq(const qfloat & x,const qfloat &y)
{
}
void __modfq(qfloat &z,qfloat &x,qfloat &y)
{
	qremain (&y, &x, &z);
}

qfloat EXPORT overloaded __fabs(qfloat & x)
{
	qfloat y;
	y = x;
	y.sign = 0;
	return y;
}
qfloat EXPORT __declspec(naked) fabsq(qfloat &x)
{
}
void __fabsq(qfloat &y,qfloat &x)
{
	y = x;
	y.sign = 0;
}

qfloat EXPORT __declspec(naked) fdimq(qfloat &x,qfloat &y)
{
}
void __fdimq(qfloat &r,qfloat &x,qfloat &y)
{
	if (x > y)
		r = x - y;
	else
		r = 0;
}

// Print decimal value of X with N digits.  Precede value with string s1,
// follow with string s2.
#include <stdio.h>
void EXPORT qprint ( const int n, const char *s1, const qfloat &x,
	const char *s2)
{
	char str[200];
	qtoasc (&x, str, n);
	printf ("%s %s %s", s1, str, s2);
}

qfloat EXPORT overloaded __declspec(naked) __sqrt(qfloat &x)
{
}
qfloat EXPORT __declspec(naked) sqrtq(qfloat &x)
{
}
void __sqrtq(qfloat &z,qfloat &x)
{
	qfsqrt(x,z);
}

qfloat EXPORT __declspec(naked) expq(qfloat &x)
{
}
void __expq(qfloat &z,qfloat &x)
{
	qfexp(&x,&z);
}

qfloat EXPORT __declspec(naked) exp10q(qfloat &x)
{
}
void __exp10q(qfloat &r,qfloat &x)
{
	qexp10(x,r);
}

qfloat EXPORT exp2q(qfloat &e)
{
	qfloat result;
	long l;

	qifrac(&e,&l,&result);
	if (0 == qcmp(result,qzero)) {
		memset(&result,0,sizeof(qfloat));
		result.sign = (e<0)?1:0;
		result.exponent = e.exponent+0x8001;
		result.mantissa[1] = 0x80000000;
		return result;
	}
	qfpow(&qtwo,&e,&result);
	return result;
}

qfloat EXPORT __declspec(naked) logq(qfloat &x)
{
}
void __logq(qfloat &z,qfloat &x)
{
	qflog(x,z);
}

qfloat EXPORT overloaded __declspec(naked) __log(qfloat &x)
{
}
void __log_real(qfloat &z,qfloat &x)
{
	qflog(x,z);
}


qfloat EXPORT __declspec(naked) log10q(qfloat &x)
{
}
void __log10q(qfloat &z,qfloat &x)
{

        qflog( x, z );
        qmul( ql10e, z, z );
}

qfloat EXPORT __declspec(naked) log2q(qfloat &x)
{
}
void __log2q(qfloat &r,qfloat &x)
{
	qlogtwo(x,r);
}

qfloat EXPORT __declspec(naked) sinq(qfloat &x)
{
}
void __sinq(qfloat &z,qfloat &x)
{
	qfsin(&x,&z);
}

qfloat EXPORT __declspec(naked) tanq(qfloat &x)
{
}
void __tanq(qfloat &z,qfloat &x)
{
	qftan(&x,&z);
}

qfloat EXPORT __declspec(naked) atanq(qfloat &x)
{
}
void __atanq(qfloat &r,qfloat &x)
{
	qatn(x,r);
}
qfloat EXPORT __declspec(naked) atan2q(qfloat &x,qfloat &y)
{
}
void __atan2q(qfloat &r,qfloat &x,qfloat &y)
{
	qatn2(x,y,r);
}

qfloat EXPORT __declspec(naked) cosq(qfloat &x)
{
}
void __cosq(qfloat &z,qfloat &x)
{
	qfcos(&x,&z);
}

qfloat EXPORT __declspec(naked) asinq(qfloat &x)
{
}
void __asinq(qfloat &r,qfloat &x)
{
	qasin(x,r);
}

qfloat EXPORT __declspec(naked) acosq(qfloat &x)
{
}
void __acosq(qfloat &r,qfloat &x)
{
	qacos(x,r);
}

qfloat EXPORT __declspec(naked) qratio(int n,int d)
{
}
void __qratio(qfloat &z,int n,int d)
{
	qfloat nn,dd;

	itoq(n,&nn);
	itoq(d,&dd);
	qdiv(dd,nn,z);
}

qfloat EXPORT __declspec(naked) qratiol(long long n,long long d)
{
}
void __qratiol(qfloat &z,long long n,long long d)
{
	qfloat 	nn,dd;

	lltoq(&n,nn);
	lltoq(&d,dd);
	qdiv(dd,nn,z);
}

qfloat EXPORT __declspec(naked) qrand(void)
{
}
void __qrand(qfloat &r)
{
	qfrand(&r);
}

void __qfloor(qfloat &r,qfloat &x);
qfloat EXPORT overloaded __floor(qfloat &x)
{
	qfloat r;
	qfloor(r,x);
	return r;
}

static unsigned int bmask[] = {
        0xffffffff,
        0xfffffffe,
        0xfffffffc,
        0xfffffff8,
        0xfffffff0,
        0xffffffe0,
        0xffffffc0,
        0xffffff80,
        0xffffff00,
        0xfffffe00,
        0xfffffc00,
        0xfffff800,
        0xfffff000,
        0xffffe000,
        0xffffc000,
        0xffff8000,
        0xffff0000,
        0xfffe0000,
        0xfffc0000,
        0xfff80000,
        0xfff00000,
        0xffe00000,
        0xffc00000,
        0xff800000,
        0xff000000,
        0xfe000000,
        0xfc000000,
        0xf8000000,
        0xf0000000,
        0xe0000000,
        0xc0000000,
        0x80000000,
        0x00000000
};
qfloat EXPORT __declspec(naked) floorq(qfloat &x)
{
}
void __floorq(qfloat &y,qfloat &x)
{
	int e;
	unsigned int *p;

	if (x.exponent == 0) {
		memset(&y,0,sizeof(qfloat));
	}
	else {
		e = x.exponent - (EXPONENT_BIAS-1);
		if (e <= 0) {
			if (x.sign) {
				qmov(qone,y);
				y.sign ^= -1;
			}
			else
				memset(&y,0,sizeof(qfloat));
		}
		else {
			e = NBITS - e;
			qmov(x,y);
			p = &y.mantissa[_NQ_-1];
			while (e >= WORDSIZE) {
				*p-- = 0;
				e -= WORDSIZE;
			}
			*p &= bmask[e];
			if (x.sign) {
				if (qcmp(x,y) != 0) {
					qadd(qminusone,y,y);
				}
			}
		}
	}
}

qfloat EXPORT __declspec(naked) ceilq(qfloat &x)
{
}

void __ceilq(qfloat &r,qfloat &x)
{
	r = floorq(x);
	if (qcmp(r,x) < 0)
		qadd(r,qone,r);
}

qfloat EXPORT overloaded __pow(qfloat &x,qfloat &y)
{
	qfloat r;

	qfpow(&x,&y,&r);
	return r;
}

qfloat EXPORT powqul(qfloat &x,unsigned long &y)
{
	qfloat r;
	qfloat yy;
	itoq(y,&yy);
	qfpow(&x,&yy,&r);
	return r;
}
qfloat EXPORT powqd(qfloat &x,double y)
{
	qfloat r,yy;
	etoq(&y,&yy);
	qfpow(&x,&yy,&r);
	return r;
}

qfloat EXPORT __declspec(naked) powq(qfloat &x,qfloat &y)
{
}
void __powq(qfloat &r,qfloat &x,qfloat &y)
{
	qfpow(&x,&y,&r);
}

qfloat EXPORT __declspec(naked) cbrtq(qfloat &x)
{
}
qfloat EXPORT overloaded __declspec(naked) __cbrt(qfloat &x)
{
}
void __cbrtq(qfloat &r,qfloat &x)
{
	qcbrt(x,r);
}

qfloat EXPORT overloaded __declspec(naked) __copysign(qfloat &x,qfloat &y)
{
}
qfloat EXPORT __declspec(naked) copysignq(qfloat &x,qfloat &y)
{
}
void __copysignq(qfloat &r,qfloat &x,qfloat &y)
{
	r = x;
	r.sign = y.sign;
}

qfloat EXPORT __declspec(naked) logbq(qfloat &x)
{
}
void __logbq(qfloat &r,qfloat &x)
{
	qlogb(x,r);
}

qfloat EXPORT __declspec(naked) acoshq(qfloat &x)
{
}
void __acoshq(qfloat &r,qfloat &x)
{
	qacosh(x,r);
}

qfloat EXPORT __declspec(naked) asinhq(qfloat &x)
{
}
void __asinhq(qfloat &r,qfloat &x)
{
	qasinh(x,r);
}

qfloat EXPORT __declspec(naked) atanhq(qfloat &x)
{
}
void __atanhq(qfloat &r,qfloat &x)
{
	qatanh(x,r);
}

qfloat EXPORT __declspec(naked) sinhq(qfloat &x)
{
}
void __sinhq(qfloat &r,qfloat &x)
{
	qsinh(x,r);
}

qfloat EXPORT __declspec(naked) coshq(qfloat &x)
{
}
void __coshq(qfloat &r,qfloat &x)
{
	qcosh(x,r);
}


qfloat EXPORT __declspec(naked) tanhq(qfloat &x)
{
}
void __tanhq(qfloat &r,qfloat &x)
{
	qtanh(x,r);
}

extern void qpowi(qfloat &,qfloat&,qfloat&);

qfloat EXPORT __declspec(naked) qipow(qfloat &x, qfloat &y)
{
}
void __qipow(qfloat &result,qfloat &x,qfloat &y)
{
	qpowi(x,y,result);
}

qfloat EXPORT __declspec(naked) frexpq(qfloat &x,int *exp)
{
}
void __frexpq(qfloat &r,qfloat &x,int *exp)
{
	qfrexp(x,exp,r);
}

int EXPORT ilogbq(qfloat &x)
{
	if (x.exponent == 0)
		return 0x80000000;
	return x.exponent - 0x8001;
}
qfloat __declspec(naked) remquoq(qfloat &x,qfloat &y,int *pint)
{
}

void __remquoq(qfloat &r,qfloat &x,qfloat &y,int *pint)
{
	*pint = qremquo(x,y,r);
}


qfloat EXPORT __declspec(naked) ldexpq(qfloat &x,long n)
{
}
void __ldexpq(qfloat &r,qfloat &x,long n)
{
	qldexp(x,n,r);
}

qfloat EXPORT __declspec(naked) expm1q(qfloat &x)
{
}
void __expm1q(qfloat &r,qfloat &x)
{
	qexpm1(x,r);
}

qfloat EXPORT __declspec(naked) log1pq(qfloat &x)
{
}
void __log1pq(qfloat &r,qfloat &x)
{
	qlog1(x,r);
}

qfloat EXPORT __declspec(naked) erfq(qfloat &x)
{
}
void __erfq(qfloat &r,qfloat &x)
{
	qerf(x,r);
}

qfloat EXPORT __declspec(naked) erfcq(qfloat &x)
{
}
void __erfcq(qfloat &r,qfloat &x)
{
	qerfc(x,r);
}

qfloat EXPORT __declspec(naked) lgammaq(qfloat &x)
{
}
void __lgammaq(qfloat &r,qfloat &x)
{
	qlgam(x,r);
}

qfloat EXPORT __declspec(naked) tgammaq(qfloat &x)
{
}
void __tgammaq(qfloat &r,qfloat &x)
{
	qgamma(x,r);
}

qfloat EXPORT __declspec(naked) igamq(qfloat &a,qfloat &x)
{
}
void __igammaq(qfloat &r,qfloat &a,qfloat &x)
{
	qigam(a,x,r);
}

qfloat EXPORT __declspec(naked) roundq(qfloat &x)
{
}
void __roundq(qfloat &r,qfloat &x)
{
	qround(x,r);
}

qfloat EXPORT __declspec(naked) scalbnq(qfloat &x,int n)
{
}
void __scalbnq(qfloat &r,qfloat &x,int n)
{
	qmov(x,r);
	r.exponent += n;
}
qfloat EXPORT __declspec(naked) fmaxq(qfloat &x,qfloat &y)
{
}
void __fmaxq(qfloat &r,qfloat &x,qfloat &y)
{
	int d = qcmp(x,y);
	if (d <= 0)
		r = y;
	else
		r = x;
}
qfloat EXPORT __declspec(naked) fminq(qfloat &x,qfloat &y)
{
}
void __fminq(qfloat &r,qfloat &x,qfloat &y)
{
        int d = qcmp(x,y);
        if (d >= 0)
                r = y;
        else
                r = x;
}

long double EXPORT fmodpiq(long double ld)
{
	qfloat c;

	e64toq(&ld,&c);
	qremain(&qpi,&c,&c);
	qtoe64(&c,&ld);
	return ld;
}

int EXPORT signbitq(qfloat &x)
{
	return x.sign;
}

qfloat EXPORT __declspec(naked) hypotq(qfloat &x,qfloat &y)
{
}
void __hypotq(qfloat &r,qfloat &x,qfloat  &y)
{
	qfloat q;
	qmul(x,x,q);
	qmul(y,y,r);
	qadd(r,q,r);
	qfsqrt(r,r);
}

//************************************************************************
// qfloat complex numbers section
//************************************************************************

extern void qcadd(qfloat _Complex&,qfloat _Complex &,qfloat _Complex &);
extern void qcsub(qfloat _Complex&,qfloat _Complex &,qfloat _Complex &);
extern void qcmul(qfloat _Complex &,qfloat _Complex &,qfloat _Complex &);
extern void qcdiv(qfloat _Complex &,qfloat _Complex &,qfloat _Complex &);
extern void qcacos(qfloat _Complex &,qfloat _Complex &);
extern void qcabs(qfloat _Complex &,qfloat &);
extern void qcsin(qfloat _Complex &,qfloat _Complex &);
extern void qccos(qfloat _Complex &,qfloat _Complex &);
extern void qctan(qfloat _Complex &,qfloat _Complex &);
extern void qccosh(qfloat _Complex &,qfloat _Complex &);
extern void qcsinh(qfloat _Complex &,qfloat _Complex &);
extern void qctanh(qfloat _Complex &,qfloat _Complex &);
extern void qcasin(qfloat _Complex &,qfloat _Complex &);
extern void qcasinh(qfloat _Complex &,qfloat _Complex &);
extern void qcatan(qfloat _Complex &,qfloat _Complex &);
extern void qcacosh(qfloat _Complex &,qfloat _Complex &);
extern void qcatanh(qfloat _Complex &,qfloat _Complex &);
extern void qcexp(qfloat _Complex &,qfloat _Complex &);
extern void qclog(qfloat _Complex &,qfloat _Complex &);
extern void qcpow(qfloat _Complex &,qfloat _Complex &,qfloat _Complex &);
extern void qcsqrt(qfloat _Complex &,qfloat _Complex &);
qfloat _Complex EXPORT __declspec(naked) operator+(qfloat _Complex &a,qfloat _Complex &b)
{
}
void __qcadd(qfloat _Complex &r,qfloat _Complex &a,qfloat _Complex &b)
{
	qcadd(a,b,r);
}
qfloat _Complex EXPORT __declspec(naked) operator-(qfloat _Complex &a,qfloat _Complex &b)
{
}
void __qcsub(qfloat _Complex &r,qfloat _Complex &a,qfloat _Complex &b)
{
	qcsub(a,b,r);
}
qfloat _Complex EXPORT __declspec(naked) operator*(qfloat _Complex &a,qfloat _Complex &b)
{
}
void __qcmul(qfloat _Complex &r,qfloat _Complex &a,qfloat _Complex &b)
{
	qcmul(a,b,r);
}
qfloat _Complex EXPORT __declspec(naked) operator/(qfloat _Complex &a,qfloat _Complex &b)
{
}
void __qcdiv(qfloat _Complex &r,qfloat _Complex &a,qfloat _Complex &b)
{
        qcdiv(a,b,r);
}

qfloat _Complex EXPORT __declspec(naked) operator()(long double _Complex a)
{
}
void __ldctoq(qfloat _Complex &r,long double _Complex ld)
{
	r.re = ld.re;
	r.im = ld.im;
}
qfloat _Complex EXPORT __declspec(naked) operator=(qfloat _Complex &a,long double _Complex b)
{
}
void __ldtocqa(qfloat _Complex &r,qfloat _Complex &a,long double _Complex b)
{
	a.re = b.re;
	a.im = b.im;
	r.re = b.re;
	r.im = b.im;
}

qfloat _Complex EXPORT __declspec(naked) cacosq(qfloat _Complex &x)
{
}
qfloat _Complex EXPORT overloaded __declspec(naked) __cacos(qfloat _Complex &x)
{
}
void __cacosq(qfloat _Complex &r,qfloat _Complex &x)
{
	qcacos(x,r);
}
qfloat _Complex EXPORT __declspec(naked) casinq(qfloat _Complex &x)
{
}
qfloat _Complex EXPORT overloaded __declspec(naked) __casin(qfloat _Complex &x)
{
}
void __casinq(qfloat _Complex &r,qfloat _Complex &x)
{
        qcasin(x,r);
}
qfloat _Complex EXPORT overloaded __declspec(naked) __casinh(qfloat _Complex &x)
{
}
qfloat _Complex EXPORT __declspec(naked) casinhq(qfloat _Complex &x)
{
}
void __casinhq(qfloat _Complex &r,qfloat _Complex &x)
{
        qcasinh(x,r);
}
qfloat _Complex EXPORT overloaded __declspec(naked) __cacosh(qfloat _Complex &x)
{
}
qfloat _Complex EXPORT __declspec(naked) cacoshq(qfloat _Complex &x)
{
}
void __cacoshq(qfloat _Complex &r,qfloat _Complex &x)
{
        qcacosh(x,r);
}
//
// Assignment operators
//
qfloat _Complex EXPORT operator=(qfloat _Complex &a,int v)
{
        itoq(v,&a.re);
	qclear(a.im);
        return a;
}
qfloat _Complex EXPORT operator=(qfloat _Complex &a,long double v)
{
        e64toq(&v,&a.re);
	qclear(a.im);
        return a;
}

qfloat _Complex EXPORT operator=(qfloat _Complex &a,long double _Complex &v)
{
#ifdef _WIN64
	etoq(&v.re,&a.re);
	etoq(&v.im,&a.im);
#else
	e64toq(&v.re,&a.re);
	e64toq(&v.im,&a.im);
#endif
	return a;
}

qfloat EXPORT  __declspec(naked) cabsq(qfloat _Complex &a)
{
}
void __qcabs(qfloat &r,qfloat _Complex &a)
{
	qcabs(a,r);
}

qfloat _Complex EXPORT __declspec(naked) csinq(qfloat _Complex &a)
{
}
void __csinq(qfloat _Complex &r,qfloat _Complex &a)
{
	qcsin(a,r);
}
qfloat _Complex EXPORT __declspec(naked) ccosq(qfloat _Complex &a)
{
}
void __ccosq(qfloat _Complex &r,qfloat _Complex &a)
{
        qccos(a,r);
}
qfloat _Complex EXPORT __declspec(naked) ctanq(qfloat _Complex &a)
{
}
void __ctanq(qfloat _Complex &r,qfloat _Complex &a)
{
        qctan(a,r);
}
qfloat _Complex EXPORT __declspec(naked) catanq(qfloat _Complex &a)
{
}
void __catanq(qfloat _Complex &r,qfloat _Complex &a)
{
        qcatan(a,r);
}
qfloat _Complex EXPORT __declspec(naked) ccoshq(qfloat _Complex &a)
{
}
void __ccoshq(qfloat _Complex &r,qfloat _Complex &a)
{
        qccosh(a,r);
}
qfloat _Complex EXPORT __declspec(naked) csinhq(qfloat _Complex &a)
{
}
void __ccsinhq(qfloat _Complex &r,qfloat _Complex &a)
{
        qcsinh(a,r);
}

qfloat _Complex EXPORT __declspec(naked) ctanhq(qfloat _Complex &a)
{
}
void __cctanhq(qfloat _Complex &r,qfloat _Complex &a)
{
        qctanh(a,r);
}
qfloat _Complex EXPORT __declspec(naked) catanhq(qfloat _Complex &a)
{
}
void __catanhq(qfloat _Complex &r,qfloat _Complex &a)
{
        qcatanh(a,r);
}
qfloat _Complex EXPORT __declspec(naked) cexpq(qfloat _Complex &a)
{
}
void __cexpq(qfloat _Complex &r,qfloat _Complex &a)
{
        qcexp(a,r);
}

qfloat _Complex EXPORT __declspec(naked) clogq(qfloat _Complex &a)
{
}
void __clogq(qfloat _Complex &r,qfloat _Complex &a)
{
        qclog(a,r);
}
qfloat _Complex EXPORT __declspec(naked) cpowq(qfloat _Complex &a,qfloat _Complex &b)
{
}
void __cpowq(qfloat _Complex &r,qfloat _Complex &a,qfloat _Complex &b)
{
        qcpow(a,b,r);
}
qfloat _Complex EXPORT __declspec(naked) csqrtq(qfloat _Complex &a)
{
}
void __csqrtq(qfloat _Complex &r,qfloat _Complex &a)
{
        qcsqrt(a,r);
}

extern void qyn(qfloat &,qfloat &,qfloat &);
extern void qjn(qfloat &,qfloat &,qfloat &);
qfloat EXPORT __declspec(naked) j0q(qfloat &x)
{
}
void __j0q(qfloat &r,qfloat &a)
{
	qjn(qzero,a,r);
}
qfloat EXPORT __declspec(naked) j1q(qfloat &x)
{
}
void __j1q(qfloat &r,qfloat &a)
{
        qjn(qone,a,r);
}
qfloat EXPORT __declspec(naked) y0q(qfloat &x)
{
}
void __y0q(qfloat &r,qfloat &a)
{
        qyn(qzero,a,r);
}

qfloat EXPORT __declspec(naked) y1q(qfloat &x)
{
}
void __y1q(qfloat &r,qfloat &a)
{
        qyn(qone,a,r);
}

qfloat EXPORT __declspec(naked) ynq(qfloat &n,qfloat &x)
{
}
void __ynq(qfloat &r,qfloat &n,qfloat &x)
{
	qyn(n,x,r);
}
qfloat EXPORT __declspec(naked) betaq(qfloat &a,qfloat &b)
{
}
void __betaq(qfloat &r,qfloat &a,qfloat &b)
{
	qfloat s;

	qadd( a, b, r );
	qgamma( r, s );

	qgamma( a, r );
	qdiv( s, r, s );

	qgamma( b, r );
	qmul( s, r, r );

}
qfloat EXPORT __declspec(naked) hypergeomq(qfloat &a,qfloat &b,qfloat &x)
{
}
void __hypergeomq(qfloat &result,qfloat &a,qfloat &b,qfloat &x)
{
	qhyp(a,b,x,result);
}
qfloat EXPORT __declspec(naked) psiq(qfloat x)
{
}

void __psiq(qfloat &result,qfloat &x)
{
	qpsi(result,x);
}
qfloat EXPORT __declspec(naked) students_tq(int k,qfloat &t)
{
}

void __qstud(qfloat &result,int k,qfloat &t)
{
	qstudt(k,t,result);
}
qfloat EXPORT __declspec(naked) students_t_invq(int k,qfloat &t)
{
}

void __qstdtri(qfloat &result,int k,qfloat &t)
{
        qstdtri(k,t,result);
}
qfloat EXPORT __declspec(naked) beta_incompleteq(qfloat &a,qfloat &b,qfloat &x)
{
}

void __betainc(qfloat &result,qfloat &a,qfloat&b,qfloat &x)
{
        qincb(a,b,x,result);
}
qfloat EXPORT __declspec(naked) ellipticKq(qfloat &x)
{
}

void __elliptick(qfloat &result,qfloat &x)
{
        qellpk(x,result);
}
qfloat EXPORT __declspec(naked) cyl_bessel_jq(qfloat &n,qfloat &x)
{
}

void __jnq(qfloat &result,qfloat &n,qfloat &x)
{
        qjn(n,x,result);
}
qfloat EXPORT __declspec(naked) expintq(qfloat x)
{
}

void __expintq(qfloat &result,qfloat x)
{
	qexpn(1,x,result);
}

qfloat EXPORT __declspec(naked) expintNq(unsigned n,qfloat x)
{
}
void __expintNq(qfloat &r,int n,qfloat x)
{
	qexpn(n,x,r);
}


qfloat EXPORT hypergeom2f1q(qfloat &a,qfloat &b,qfloat &c,qfloat &x)
{
	qfloat sum_pos = qone;
	qfloat sum_neg = qzero;
	qfloat del_pos = qone;
	qfloat del_neg = qzero;
	qfloat del = qone;
	qfloat k = qzero;
	int i = 0;

	do {
		if(++i > 500) {
			mtherr("hypergeom2f1q",PLOSS);
			break;
		}
		del *= (a+k)*(b+k) * x / ((c+k) * (k+1.0));  /* Gauss series */

		if (del > 0) {
			del_pos  =  del;
			sum_pos +=  del;
		}
		else if (del == 0) {
			/* Exact termination (a or b was a negative integer).  */
			break;
		}
		else {
			del_neg  = -del;
			sum_neg -=  del;
		}

		k += 1.0;
	} while(fabsq((del_pos + del_neg)/(sum_pos-sum_neg)) > 1e-108Q);

	return  sum_pos - sum_neg;
}
qfloat EXPORT hypergeom1f1q(const qfloat a, const qfloat b, const qfloat x)
{
  qfloat an  = a;
  qfloat bn  = b;
  qfloat n   = 1.0;
  qfloat del = 1.0;
  qfloat abs_del = 1.0;
  qfloat max_abs_del = 1.0;
  qfloat sum_val = 1.0;

  while(abs_del/fabsq(sum_val) > 1e-105) {
    qfloat u, abs_u;

    if(bn == 0.0) {
      mtherr("hypergeom1f1", DOMAIN);
     return 0;
    }
    if(an == 0.0 || n > 1000.0) {
      return sum_val;
    }

    u = x * (an/(bn*n));
    abs_u = fabsq(u);
    if(abs_u > 1.0 && max_abs_del > 1e5000Q/abs_u) {
ovfl:
      mtherr("hypergeom1f1", OVERFLOW);
      return sum_val;
    }
    del *= u;
    sum_val += del;
    if(fabsq(sum_val) > 1e5000Q) {
        goto ovfl;
    }

    abs_del = fabsq(del);
    max_abs_del =(abs_del> max_abs_del) ? abs_del : max_abs_del;

    an += 1.0;
    bn += 1.0;
    n  += 1.0;
  }

  return sum_val;
}

qfloat EXPORT __declspec(naked) beta_distribution_invq(qfloat &a,qfloat &b, qfloat &p)
{
}
void __beta_distribution_invq(qfloat &r,qfloat &a,qfloat &b,qfloat &p)
{
	beta_distribution_invQ(a,b,p,r);
}
qfloat EXPORT __declspec(naked) PolyLogq(int n,qfloat &x)
{
}
void  __PolyLogq(qfloat &r,int n,qfloat &x)
{
	qpolylog(n,x,r);
}

/*                                                      qbdtrc
 *
 *      Complemented binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int qbdtrc( k, n, p, y );
 * int k, n;
 * QELT *p, *y;
 *
 * y = qbdtrc( k, n, p, y );
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms k+1 through n of the Binomial
 * probability density:
 *
 *   n
 *   --  ( n )   j      n-j
 *   >   (   )  p  (1-p)
 *   --  ( j )
 *  j=k+1
 *
 * The terms are not summed directly; instead the incomplete
 * beta integral is employed, according to the formula
 *
 * y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).
 *
 * The arguments must be positive, with p ranging from 0 to 1.
 *
 */
qfloat binomial_cq(int k,int n,qfloat &p)
{
	qfloat y,dk,dn;
	long li;

	if (k < 0) {
		qmov(qone,y);
	}
	else if (k == n) {
		qclear(y);
	}
	else {
		li = k+1;
		itoq(li,&dk);
		li = n-k;
		itoq(li,&dn);
		qincb(dk,dn,p,y);
	}
	return y;
}

int __asctoq(char *,qfloat *,char **);
qfloat strtoq(const char *p,char **pend)
{
	qfloat q;
	__asctoq(p,&q,pend);
	return q;
}


#if 0
static qfloat dilogseries(qfloat &x)
{
	/* series around x = 1.0 */
	qfloat eps = x - qone;
	qfloat lne;
	qflog(eps,lne);
	qfloat c0 = qpi*qpi/6.0Q;
	qfloat c1 =   qone - lne;
	qfloat c2 = -(qone - 2.0Q*lne)/4.0Q;
	qfloat c3 =  (1.0Q - 3.0Q*lne)/9.0Q;
	qfloat c4 = -(1.0Q - 4.0Q*lne)/16.0Q;
	qfloat c5 =  (qone - 5.0Q*lne)/25.0Q;
	qfloat c6 = -(qone - 6.0Q*lne)/36.0Q;
	qfloat c7 =  (qone - 7.0Q*lne)/49.0Q;
	qfloat c8 = -(qone - 8.0Q*lne)/64.0Q;
	qfloat c9 =  (qone - 9.0Q*lne)/81.0Q;
	qfloat result = c0+eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*(c7+eps*c8+eps*(c9*eps))))))));
}

qfloat EXPORT dilog(qfloat &x)
{
	qfloat r,tmp = qone-x;
	if (tmp <1.1 && tmp > 1)
		r = dilogseries(tmp);
	else qspenc(tmp,r);
	return r;
}
#else
qfloat EXPORT __declspec(naked) dilogq(qfloat &x)
{

}
void __dilogq(qfloat &r,qfloat &x)
{
	qpolylog(2,x,r);
}
#endif
#include <qfloat.h>

/* Implementation of Lamberts W-function which is defined as
 * w(x)*e^(w(x))=x
 * Implementation by Gunter Kuhnle, gk@uni-leipzig.de
 * Algorithm originally developed by
 * KEITH BRIGGS, DEPARTMENT OF PLANT SCIENCES,
 * ANSI C code for W(x).  K M Briggs 98 Feb 12
 * http://keithbriggs.info/W-ology.html

   Based on Halley iteration.  Converges rapidly for all valid x.

  double W(const double x) {
  int i; double p,e,t,w,eps=4.0e-16;  eps=desired precision
  if (x<-0.36787944117144232159552377016146086) {
    fprintf(stderr,"x=%g is < -1/e, exiting.\n",x); exit(1); }
  if (0.0==x) return 0.0;
   get initial approximation for iteration...
  if (x<1.0) {  series near 0
    p=sqrt(2.0*(2.7182818284590452353602874713526625*x+1.0));
    w=-1.0+p-p*p/3.0+11.0/72.0*p*p*p;
  } else w=log(x);  asymptotic
  if (x>3.0) w-=log(w);
  for (i=0; i<20; i++) {  Halley loop
    e=exp(w); t=w*e-x;
    t/=e*(w+1.0)-0.5*(w+2.0)*t/(w+1.0); w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w;  rel-abs error
  }  never gets here
  fprintf(stderr,"No convergence at x=%g\n",x); exit(1);
}
*/

extern qfloat qmem1;
qfloat __declspec(naked) lambertwq(qfloat &x) {}
void EXPORT __lambertwq(qfloat &result,qfloat &x)
{
    qfloat p, e, t, w, fabs_x,fabs_t,fabs_w,t1,t2,t3;
    int i;


    if (qcmp(x,qmem1) < 0)      /* if (x < -exp(-1.0))    */
	goto errexit;              /* error, value undefined */

     __fabsq(fabs_x,x);

    if (qcmp(fabs_x,qepsilon) <= 0) { //if (fabs(x) <= eps)
	memcpy(&result,&x,sizeof(qfloat));
	return;
    }

    if (qcmp(x,qone) < 0) {    // if (x < 1)
        qmul(twiceqexp1,x,p); // p = sqrtq(2.0q * (exp(1.0q) * x + 1.0q));
	qadd(p,qone,p);
	qfsqrt(p,p);
	w = -1.0q + p - p * p / 3.0 + 11.0q / 72.0q * p * p * p;
    } else {
	//w = log(x);
	qflog(x,w);
    }

    if (x > 3) {
	//w = w - log(w);
	qflog(w,t1);
	qsub(t1,w,w);
    }
    for (i = 0; i < 20; i++) {
		//e = exp(w);
		qfexp(&w,&e);

		t = w * e - x;
		qmul(w,e,t);
		qsub(x,t,t);

		//t = t / (e * (w + 1.0) - 0.5q * (w + 2.0q) * t / (w + 1.0));
		qadd(qone,w,t1);
		qmul(e,t1,t3);

		qadd(qtwo,w,t2);
		qmul(qhalf,t2,t2);
		qdiv(t1,t,t1);
		qmul(t1,t2,t2);
		qsub(t2,t3,t3);
		qdiv(t3,t,t);
		
		//w = w - t;
		qsub(t,w,w);
		__fabsq(fabs_t,t);
		__fabsq(fabs_w,w);
		qadd(fabs_w,qone,fabs_w);
		qmul(qepsilon,fabs_w,fabs_w);
		if (qcmp(fabs_t,fabs_w) < 0) {
		//if (fabs_t < qepsilon * (1.0 + fabs(w))) {
	            memcpy(&result,&w,sizeof(qfloat));
		    return;
        }
    }
errexit:
    memcpy(&result,&qminusone,sizeof(qfloat));                 /* error: iteration didn't converge */
}

qfloat EXPORT __declspec(naked) zetacq(qfloat &x) {}
void __qzetacq(qfloat &r,qfloat &x)
{
	qzetac(x,r);
}

/*------------------------------------------------------------------------
Procedure:     LibMain ID:1
Purpose:       Dll entry point.Called when a dll is loaded or
unloaded by a process, and when new threads are
created or destroyed.
Input:         hDllInst: Instance handle of the dll
fdwReason: event: attach/detach
lpvReserved: not used
Output:        The return value is used only when the fdwReason is
DLL_PROCESS_ATTACH. True means that the dll has
sucesfully loaded, False means that the dll is unable
to initialize and should be unloaded immediately.
Errors:
------------------------------------------------------------------------*/
int _stdcall EXPORT LibMain(void * hDLLInst, unsigned long fdwReason,void * lpvReserved)
{
	return 1;
}
