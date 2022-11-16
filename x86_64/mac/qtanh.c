/*							qtanh.c
*      Hyperbolic tangent check routine
*
 *
 *
 *
 * SYNOPSIS:
*
 * int qtanh( x, y );
* QELT *x, *y;
*
 * qtanh( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * For x >= 1 the program uses the definition
*
 *             exp(x) - exp(-x)
* tanh(x)  =  ----------------
*             exp(x) + exp(-x)
*
 *
 * For x < 1 the method is a continued fraction
*
 *                   2   2   2
*              x   x   x   x
* tanh(x)  =  --- --- --- --- ...
*              1+  3+  5+  7+
*
* Cephes Math Library Release 2.3:  March, 1995
* Copyright 1984, 1995, 1996 by Stephen L. Moshier
* Adapted to lcc-win by Jacob Navia, 1998-2008
*/


#include "qhead.h"

void qtanh(Qfloatp x,Qfloatp y)
{
	Qfloat e[1], r[1], j[1], xx[1], m2[1];
	int i, n, sign;
	int lj;

	sign = signof(x);
	qmov( x, r );
	setpositive(r);
	if( qcmp(r, qone) >= 0 ) {
		/* This subroutine is used by the exponential function routine.
		* tanh(x) = (exp(x) - exp(-x)) / (exp(x) + exp(-x))
		* Note qfexp() calls qtanh, but with an argument less than (1 + log 2)/2.
		*/
		qfexp( r, e );
		qdiv( e, qone, r );
		qsub( r, e, xx );
		qadd( r, e, j );
		qdiv( j, xx, y );
		goto done;
	}

	qmov( qtwo, m2 );
	qneg( m2 );

	/* Adjust loop count for convergence to working precision.  */
	n = NBITS/7 + 1; /*10;*/
	lj = 2 * n + 1;
	itoq(lj, j );

	qmov( j, e );
	qmul( x, x, xx );

	/* continued fraction */

	for( i=0; i<n; i++)
	{
		qdiv( e, xx, r );
		qadd( m2, j, j );
		qadd( r, j, e );
	}

	qdiv( e, x, y );

done:
	if( sign != 0 )
		setnegative(y);
}
