/*							qsinh.c
 *
 *	Hyperbolic sine check routine
 *
 *
 *
 * SYNOPSIS:
 *
 * int qsinh( x, y );
 * QELT *x, *y;
 *
 * qsinh( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * The range is partitioned into two segments.  If |x| <= 1/4,
 *
 *                3    5    7
 *               x    x    x
 * sinh(x) = x + -- + -- + -- + ...
 *               3!   5!   7!
 *
 * Otherwise the calculation is sinh(x) = ( exp(x) - exp(-x) )/2.
 *
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
Revised by J.N. 2018
 */
#include "qhead.h"

void EXPORT qsinh(Qfloatp x,Qfloatp y)
{
	long long n;
	int r;
	Qfloat xx[1], qf[1], qz[1];

	if( exponent(x) < 2 ) {
		qclear( y );
	}
	else if( (long)exponent(x) < (QELT) (EXPONE - 1) ) {
		qsquare( x, xx );
		qmov( qone, qz );
		qmov( qone, qf );
		n = 1;
		do {
			n++;
			qdivi( n, qf, qf );
			n++;
			qdivi( n, qf, qf );
			qmul( xx, qf, qf );
			r = qadd( qz, qf, qz );
		}
		while( r != 0 );
		qmul( x, qz, y );
	}
	else {
		qfexp( x, qz );
		qinv( qz, y );
		qsub( y, qz, y );
		decreaseExponent(y);
	}
}
