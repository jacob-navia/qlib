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
 * Adapted to lcc-win by Jacob Navia, 1998-2018
*/


#include "qhead.h"

/* ispow2 - if u > 1 && u == 2^n, return n, otherwise return 0 */
int ispow2(unsigned u)
{
	int n;

	if (u > 1 && (u&(u-1)) == 0)
		for (n = 0; u; u >>= 1, n++)
			if (u&1)
				return n;
	return 0;
}

void qtanh(Qfloatp x,Qfloatp y)
{
	Qfloat e[1], r[1], j[1], xx[1],*p;
	int i, n, sign;
	long lj,toshift;

	sign = signof(x);
	if (sign) {
		qmov( x, r );
		setpositive(r);
		p = r;
	}
	else p = x;
	if( qcmpto1(p) >= 0 ) {
		/* This subroutine is used by the exponential function routine.
		* tanh(x) = (exp(x) - exp(-x)) / (exp(x) + exp(-x))
		* Note qfexp() calls qtanh, but with an argument less than (1 + log 2)/2.
		*/
		qfexp( p, e );
		qinv( e, p );
		qsub( p, e, xx );
		qadd( p, e, j );
		qdiv( j, xx, y );
		goto done;
	}

	/* Adjust loop count for convergence to working precision.  */
	n = NBITS/7 + 1; /* 65 */
	lj = 2 * n + 1;  /* 131 */

	/* 
	Inlining the first pass of the loop we spare a move and a 480 bit add 
	This saves 5 seconds (!) in qmtst.
	*/
	qsquare( x, xx );
	lj = lj-2;
	itoq(lj, j );
	qdivi(lj, xx, r);
	qadd(r, j, e);

	//qmov( j, e );

	/* continued fraction */

	for( i=1; i<n; i++) {
		qdiv( e, xx, r );
		lj = lj-2;
//		itoq(lj,j);
		toshift=bsr64(lj);
		
		j[0].mantissa[0] = lj << toshift;
		j[0].exponent = (EXPONE-1) + 64 - toshift;
		qadd( r, j, e );
	}

	qdiv( e, x, y );

done:
	if( sign != 0 )
		setnegative(y);
}
