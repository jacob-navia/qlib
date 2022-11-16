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
*/
#include "qhead.h"

void EXPORT qsinh(Qfloatp x,Qfloatp y)
{
	Qfloat xx[1], qn[1], qf[1], qz[1];

	if( exponent(x) < 2 ) {
		qclear( y );
	}
	else if( (long)exponent(x) < (QELT) (EXPONE - 1) ) {
		qmul( x, x, xx );
		qmov( qone, qz );
		qmov( qone, qf );
		qmov( qone, qn );
		do {
			qadd( qone, qn, qn );
			qdiv( qn, qf, qf );
			qadd( qone, qn, qn );
			qdiv( qn, qf, qf );
			qmul( xx, qf, qf );
			qadd( qz, qf, qz );
		}
		while( (int)(exponent(qz) - (int)exponent(qf)) < NBITS );
		qmul( x, qz, y );
	}
	else {
		qfexp( x, qz );
		qdiv( qz, qone, y );
		qsub( y, qz, y );
		decreaseExponent(y);
	}
}
