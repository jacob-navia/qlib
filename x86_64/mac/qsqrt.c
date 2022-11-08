/*							qfsqrt.c
 *      Square root check routine
 *
 *
 *
 * SYNOPSIS:
 *
 * int qfsqrt( x, y );
 * QELT *x, *y;
 *
 * qfsqrt( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the square root of x.
 *
 * Range reduction involves isolating the power of two of the
 * argument and using a polynomial approximation to obtain
 * a rough value for the square root.  Then Heron's iteration
 * is used to converge to an accurate value.
 *
 */

/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988, 1996 by Stephen L. Moshier
 */

/* full precision for 9 word mantissa */
#include "qhead.h"

/* constants for first approximation polynomial */

static Qfloat qsq2[1] = { -1,EXPONE-3,0xd14fc42fe79ba800ULL};

static Qfloat qsq1[1] = { 0,EXPONE-1,0xe3e3c2ae4c146700ULL};

static Qfloat qsq0[1] = { 0,EXPONE-2,0xa08bdc7dd5ffe300ULL};

//#define DEBUG
#ifdef DEBUG
static void prt(Qfloatp x,char *title)
{
        char buf[512];
        __qtoasc(x,buf,135);
        printf("%s: %s\n",title,buf);
}
#else
#define prt(a,b)
#endif


long double sqrtl(long double);
void qfsqrt(Qfloatp x,Qfloatp y)
{
	int i, m;
	long double d=0;
	Qfloat a[1], *b, tmp[1];

	if( signof(x) != 0 ) goto negative;
	if( exponent(x) == 0 ) goto iszero;
	if (exponent(x) > (EXPONE - 13600) && exponent(x) < EXPONE +0x3500) {
		unsigned short *e = (unsigned short *)&d;
		unsigned long long *p = (unsigned long long *)e;
		e += 4;
		*e |= x->exponent - (EXPONE-0x3fff);
		*p = x->mantissa[0];

		d = sqrtl(d);

		tmp[0].sign = 0;
		tmp[0].mantissa[0] = *p;
		tmp[0].exponent = *e +(EXPONE-0x3fff);
		memset(&tmp->mantissa[1],0,(MANTISSA_LENGTH-1)*sizeof(long long));

		prt(tmp,"x=");

		qdiv(tmp,x,a);
		qadd(a,tmp,tmp);
		tmp[0].exponent -= 1;

		qdiv(tmp,x,a);
		qadd(a,tmp,tmp);
		tmp[0].exponent -= 1;

		qdiv(tmp,x,a);
		qadd(a,tmp,y);
		y[0].exponent -= 1;
		return;
	}

	qmov( x, a );
	m = exponent(a);	/* save the exponent		*/
	a[0].exponent = EXPONE - 1;	/* change range to 0.5, 1 */

	/* b = ( a * qsq2 + qsq1) * a + qsq0		*/
	// This must be in intermediate storage to cope with the case of x == y
	if (x == y)
		b = tmp;
	else
		b = y;
	qmul( a, qsq2, b );	/* b = a * qsq2		*/
	qadd( qsq1, b, b );	/* b += qsq1		*/
	qmul( a, b, b );	/* b *= a;		*/
	qadd( qsq0, b, b );	/* b += qsq0;		*/

	/* divide exponent by 2 */
	m -= EXPONE - 1;
	b[0].exponent = (m / 2) + (EXPONE - 1);


	/* multiply by sqrt(2) if exponent odd */
	if( (m & 1) != 0 )
		qmul( b, qsqrt2, b );


	/* Newton iterations */

	for( i=0; i<9; i++ ) {
		qdiv( b, x, a );
		qadd( a, b, b );
		b[0].exponent -= 1;
	}
	if (x == y)
		qmov( b, y );
	return;
negative:
	qclear(y);
	if (qequal(qzero,x)) qneg(y);
	return;
iszero:
	qclear(y);
	return;
}


