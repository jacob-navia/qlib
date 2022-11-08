/*							qatanh.c
*
 *	Inverse hyperbolic tangent
*
 *
 *
 * SYNOPSIS:
*
 * int qatanh( x, y );
* QELT x[], y[];
*
 * qatanh( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns inverse hyperbolic tangent of argument.
*
 *        atanh(x) = 0.5 * log( (1+x)/(1-x) ).
*
 * For very small x, the first few terms of the Taylor series
* are summed.
*
* Cephes Math Library Release 2.3:  March, 1995
* Copyright (C) 1987, 1995 by Stephen L. Moshier
* Adapted to lcc-win by jacob navia
*/

#include "qhead.h"

void EXPORT qatanh(Qfloatp x,Qfloatp y)
{
	Qfloat a[1], b[1];
	int sign;

	qmov( x, a );
	sign = signof(a);
	setpositive(a);
	if( qcmp( a, qone ) >= 0 ) {
		mtherr( "qatanh", DOMAIN );
		qinfin(y);
		setsign(y,sign);
	}
	else if( (EXPONE - exponent(a)) >= NBITS/4 ) {
		/* x + x^3 / 3 + x^5 / 5  */
		qmul( a, a, b );
		qmul( a, b, y ); /* x^3 */

		qmul( b, y, b ); /* x^5 */
		qdiv( qthree, y, y );

		qdiv( qfive, b, b );
		qadd( b, y, y );
		qadd( a, y, y );
		setsign(y ,sign);
	}
	else {/* 0.5 * log((1+x)/(1-x)) */
		qadd( x, qone, a );
		qsub( x, qone, b );
		qdiv( b, a, a );
		qflog( a, y );
		decreaseExponent(y);
	}
}
