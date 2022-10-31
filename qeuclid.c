/* Rational arithmetic routines
*
 * radd( a, b, c )	c = b + a
* rsub( a, b, c )	c = b - a
* rmul( a, b, c )	c = b * a
* rdiv( a, b, c )	c = b / a
* euclid( n, d )	Reduce n/d to lowest terms, return g.c.d.
*
 * Note: arguments are assumed,
* without checking,
* to be integer valued.
*/

#include "qhead.h"
#include <stdlib.h> // For exit()

/* Integer overflow threshold */
#define BIG (EXPONE + NBITS)
/*
typedef struct { Qfloat n[1]; Qfloat d[1]; }qfract;
*/

void qeuclid(Qfloatp,Qfloatp,Qfloatp);


/* Add fractions. */

void qradd(qfract *ff1,qfract *ff2,qfract *ff3 )
{
	Qfloat gcd[1];
	Qfloat d1[1];
	Qfloat d2[1];
	Qfloat gcn[1];
	Qfloat n1[1];
	Qfloat n2[1];

	qmov( ff1->n, n1 );
	qmov( ff1->d, d1 );
	qmov( ff2->n, n2 );
	qmov( ff2->d, d2 );
	if( n1[0].exponent == 0 )
	{
		qmov( n2, ff3->n );
		qmov( d2, ff3->d );
		return;
	}
	if( n2[0].exponent == 0 )
	{
		qmov( n1, ff3->n );
		qmov( d1, ff3->d );
		return;
	}
	qeuclid( d1, d2, gcd );
	qeuclid( n1, n2, gcn );

	/* f3->n = (n2 * d1 + n1 * d2) * gcn; */
	qmul( n2, d1, n2 );
	qmul( n1, d2, n1 );
	qadd( n1, n2, n2 );
	qmul( gcn, n2, ff3->n );

	/* f3->d = d1 * d2 * gcd;*/
	qmul( d1, d2, d2 );
	qmul( d2, gcd, ff3->d );
	qeuclid( ff3->n, ff3->d, gcd );
}




/* Subtract fractions. */

void qrsub(qfract *ff1,qfract *ff2,qfract *ff3 )
{
        Qfloat gcd[1];
        Qfloat d1[1];
        Qfloat d2[1];
        Qfloat gcn[1];
        Qfloat n1[1];
        Qfloat n2[1];

	qmov( ff1->n, n1 );
	qmov( ff1->d, d1 );
	qmov( ff2->n, n2 );
	qmov( ff2->d, d2 );
	if( n1[0].exponent == 0 )
	{
		qmov( n2, ff3->n );
		qmov( d2, ff3->d );
		return;
	}
	if( n2[0].exponent == 0 )
	{
		qneg( n1 );
		qmov( n1, ff3->n );
		qmov( d1, ff3->d );
		return;
	}
	qeuclid( d1, d2, gcd );
	qeuclid( n1, n2, gcn );

	/* f3->n = (n2 * d1 - n1 * d2) * gcn; */
	qmul( n2, d1, n2 );
	qmul( n1, d2, n1 );
	qsub( n1, n2, n2 );
	qmul( gcn, n2, ff3->n );

	/* f3->d = d1 * d2 * gcd;*/
	qmul( d1, d2, d2 );
	qmul( d2, gcd, ff3->d );
	qeuclid( ff3->n, ff3->d, gcd );
}

/* Multiply fractions. */

void qrmul(qfract *ff1,qfract *ff2,qfract *ff3 )
{
        Qfloat gcd[1];
        Qfloat d1[1];
        Qfloat d2[1];
        Qfloat n1[1];
        Qfloat n2[1];
	qmov( ff1->n, n1 );
	qmov( ff1->d, d1 );
	qmov( ff2->n, n2 );
	qmov( ff2->d, d2 );
	if( (n1[0].exponent == 0) || (n2[0].exponent == 0) )
	{
		qclear( ff3->n );
		qmov( qone, ff3->d );
		return;
	}

	qeuclid( n1, d2, gcd ); /* cross cancel any common divisors */
	qeuclid( n2, d1, gcd );
	qmul( n1, n2, ff3->n );
	qmul( d1, d2, ff3->d );

	/* Check for overflow. */
	if( (ff3->n[0].exponent >= (QELT) BIG) || (ff3->d[0].exponent >= (QELT) BIG) )
	{
		mtherr( "qrmul", OVERFLOW );
		return;	/* terminate program */
	}
}




/* Divide fractions. */

void qrdiv(qfract *ff1,qfract *ff2, qfract *ff3 )
{
        Qfloat gcd[1];
        Qfloat d1[1];
        Qfloat d2[1];
        Qfloat n1[1];
        Qfloat n2[1];
	/* Invert ff1, then multiply */
	qmov( ff1->d, n1 );
	qmov( ff1->n, d1 );
	if( n1[0].sign != 0 )
	{ /* keep denominator positive */
		qneg(n1);
		qneg(d1);
	}
	qmov( ff2->n, n2 );
	qmov( ff2->d, d2 );
	if( (n1[0].exponent == 0) || (n2[0].exponent == 0) )
	{
		qclear( ff3->n );
		qmov( qone, ff3->d );
		return;
	}
	qeuclid( n1, d2, gcd ); /* cross cancel any common divisors */
	qeuclid( n2, d1, gcd );
	qmul( n1, n2, ff3->n );
	qmul( d1, d2, ff3->d );

	/* Check for overflow. */
	if( (ff3->n[0].exponent >= (QELT) BIG) || (ff3->d[0].exponent >= (QELT) BIG) )
	{
		mtherr( "qrdiv", OVERFLOW );
		return;	/* terminate program */
	}
}

/* Euclidean algorithm
*   reduces fraction to lowest terms,
*   returns greatest common divisor.
*/

void qeuclid(Qfloatp num,Qfloatp den,Qfloatp gcda )
{
	Qfloat nn[1], dd[1], q[1], r[1];


	/* Numerator. */
	qmov( num, nn );
	/* Denominator. */
	qmov( den, dd );

	/* Make numbers positive, locally. */
	nn[0].sign = 0;
	dd[0].sign = 0;

	/* Abort if numbers are too big for integer arithmetic. */
	if( (nn[0].exponent >= (QELT) BIG) || (dd[0].exponent >= (QELT) BIG) )
	{
		mtherr( "qeuclid", OVERFLOW );
		qmov( qone, gcda );
		return;
	}


	/* Divide by zero, gcd = 1. */
	if( dd[0].exponent <= (QELT) (EXPONE - 1) )
	{
		qmov( qone, gcda );
		return;
	}

	/* Zero. Return 0/1, gcd = denominator. */
	if( nn[0].exponent <= (QELT) (EXPONE - 1) )
	{
		qmov( qone, den );
		qmov( dd, gcda );
		return;
	}

	while( dd[0].exponent > (QELT) (EXPONE - 1) )
	{
		/* Find integer part of n divided by d. */
		qdiv( dd, nn, r );
		qfloor( r, q );
		/* Find remainder = n - d*q after dividing n by d. */
		qmul( dd, q, q );
		qsub( q, nn, r );
		/* The next fraction is d/r. */
		qmov( dd, nn );
		qmov( r, dd );
	}

	qdiv( nn, num, num );
	qdiv( nn, den, den );
	qmov( nn, gcda );
}

