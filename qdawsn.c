/*							qdawsn.c
 *
 *	Dawson's Integral
 *
 *
 *
 * SYNOPSIS:
 *
 * int qdawsn( x, y );
 * QELT *x, *y;
 *
 * qdawsn( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *                             x
 *                             -
 *                      2     | |        2
 *  dawsn(x)  =  exp( -x  )   |    exp( t  ) dt
 *                          | |
 *                           -
 *                           0
 *
 *
 *
 * ACCURACY:
 *
 * Series expansions are truncated at NBITS/2.
 *
 */

/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
*/


/* In 144 bit arithmetic, accuracy at crossover between
 * ascending series and asymptotic series is 72 bits.
 *
 *	x	ascending	asymptotic
 *		series		series
 *		loss, bits	accuracy, bits
 *
 *	 5	 35		 36
 *	 7	 70		 71
 *	 8	 91		 92
 *	 9	117		117
 *	10	130		144
 */

#include "qhead.h"

#define DEBUG 0

int qdawsn(Qfloat xx[],Qfloat y[] )
{
Qfloat a[1], x[1], q[1], r[1], s[1], t[1], sum[1];

double dx;
QELT it;


qmul( xx, xx, x );	/* 	x = -xx**2	*/
x[0].sign = -1;
qmov( qone, a );	/*	a = 0.5	*/
a[0].exponent -= 1;

dx = qtoe( xx, NOROUNDING );

#if !DEBUG
#if NBITS > 144
if( dx > 10.0 )
#else
if( dx > 7.1 )
#endif
{
#endif

/* asymptotic expansion */

qmov( qone, sum );	/*	sum = 1.0	*/
qmov( qone, r );	/*	r = 1.0	*/
qmov( qone, s );	/*	s = 1.0	*/
qmov( a, q );		/*	q = a	*/
it = 32767;

do
	{
	qsub( qone, q, q );	/*	q -= 1.0	*/
	qdiv( x, q, t );	/*	r *= q/x	*/
	qmul( r, t, r );

	if( r[0].exponent > it )
		goto asydon;
	it = r[0].exponent;

	qadd( sum, r, sum );	/*	sum += r	*/
	}
while( ((int) sum[0].exponent - (int) r[0].exponent) < NBITS/2 );

asydon:

#if DEBUG
printf( "qdawsn asymp %d\n", qone[1] - it );
#endif
qdiv( xx, sum, y );
y[0].exponent -= 1;

#if !DEBUG
return(0);
}
#endif

/* ascending power series */

qmov( qone, r );	/*	r = 1.0	*/
qmov( qone, sum );	/*	sum = 1.0	*/
qmov( a, q );		/*	q = a	*/
#if DEBUG
it = 0;
#endif

do
	{
	qadd( qone, q, q );	/*	q += 1.0	*/
	qdiv( q, x, s );	/*	r *= x/q	*/
	qmul( r, s, r );
#if DEBUG
	if( r[0].exponent > it )
		it = r[0].exponent;
#endif
	qadd( sum, r, sum );	/*	sum += r	*/
	}
while( ((int) sum[0].exponent -(int) r[0].exponent) < NBITS/2 );

#if DEBUG
it = it - sum[1];
printf( "qdawsn ascend %d bits\n", NBITS - it );
#endif

qmul( xx, sum, y );
return 0;
}
