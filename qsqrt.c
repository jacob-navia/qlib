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
float128_t sqrtf128(float128_t x); 
float128_t inversef128(float128_t);
#ifdef aarch64 // ARM64 linux
#define DO_LONG_DOUBLE
#endif
#ifdef DO_LONG_DOUBLE
#define DOUBLE long double
#define TO_DOUBLE qtoe113
#define FROM_DOUBLE e113toq
#define SQRT sqrtl
#else
#define DOUBLE float128_t
#define SQRT sqrtf128
#define TO_DOUBLE qtoe113
#define FROM_DOUBLE e113toq
#endif

/* constants for first approximation polynomial */

static Qfloat qsq2[1] = { {-1,EXPONE-3,{0xd14fc42fe79ba800ULL}}};

static Qfloat qsq1[1] = { {0,EXPONE-1,{0xe3e3c2ae4c146700ULL}}};

static Qfloat qsq0[1] = { {0,EXPONE-2,{0xa08bdc7dd5ffe300ULL}}};

/*
    X     =  0.5 * (X   + a/X   )
     n+1             n       n
*/

//static Qfloat epsilon={{0,7fe41,{0x8000000000000000ULL,0,0,0,0,0,0}}};
#include <stdio.h>
void qfsqrt(const Qfloatp x,Qfloatp y)
{
	int i, m;
	DOUBLE d;
	Qfloat a[1], *b, tmp[1];

	if( signof(x) != 0 ) goto negative;
	if( exponent(x) == 0 ) goto iszero;
	if (exponent(x) > (EXPONE - 1023) && exponent(x) < EXPONE +1023) {
		/*
                y * (3 - xy^2)
		y    = ------------
		 n+1       2

		*/
        Qfloat xsq[1];
        d = TO_DOUBLE(x);
        d = SQRT(d);
        d = inversef128(d);
        FROM_DOUBLE(d,tmp);

        qsquare(tmp,xsq);       // y^2
        qmul(x,xsq,a);          // x * y^2
        qsub(a,qthree,a);       // 3 - x * y^2
        qmul(tmp,a,tmp);        // y * (3 - x * y^2)
        tmp[0].exponent -= 1;   // ( y * (3 - x * y^2) ) / 2

        qsquare(tmp,xsq);
        qmul(x,xsq,a);
        qsub(a,qthree,a);
        qmul(tmp,a,tmp);
        tmp[0].exponent -= 1;

		qmul(tmp,x,y);
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

#ifdef aarch64
float128_t inversef128(float128_t n) { return 1.0/n; }
#endif
