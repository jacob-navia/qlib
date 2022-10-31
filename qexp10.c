/*							exp10.c
 *
 *	Base 10 exponential function
 *      (Common antilogarithm)
 *
 *
 *
 * SYNOPSIS:
 *
 * int qexp10( x, y );
 * QELT *x, *y;
 *
 * qexp10( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns 10 raised to the x power.
 *
 *   x      x ln 10
 * 10   =  e
 *
 * Cephes Math Library Release 2.2:  January, 1991
 * Copyright 1984, 1991 by Stephen L. Moshier
 */

#include "qhead.h"

void qexp10(Qfloatp x,Qfloatp y)
{
	Qfloat a[1];

	qmul( x, qlog10c, a );
	qfexp( a, y );
}

#if 0
/* qexp11.c
 *
 *  10**x - 1 for small x.
 *
 *             1 + tanh x/2
 *  exp(x)  =  ------------
 *             1 - tanh x/2
 *
 *
 *                    2 tanh x/2
 *  exp(x) - 1   =   ------------
 *                   1 - tanh x/2
 *
 *  exp10(x)  =  exp( x log 10 )
 *
 */

void qexp11(Qfloatp xx,Qfloatp y )
{
	Qfloat num[1], den[1], x[1];

	if( xx->exponent == 0 ) {
		qclear(y);
		return;
	}

	qmul( xx, qlog10c, x );	/*   x * log(10)   */

	x[0].exponent -= 1;		/* x/2				*/
	qtanh( x, num );	/* tanh( x/2 )			*/

	qsub( num, qone, den );	/* 1 - tanh			*/
	num[0].exponent += 1;
	qdiv( den, num, y );	/* (2 * tanh)/(1 - tanh)	*/
}
#endif
