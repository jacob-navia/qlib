/*							qellpk.c
 *
 *	Complete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * int qellpk(x, y);
 * QELT *x, *y;
 *
 * qellpk(x, y);
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *            pi/2
 *             -
 *            | |
 *            |           dt
 * K(m)  =    |    ------------------
 *            |                   2
 *          | |    sqrt( 1 - m sin t )
 *           -
 *            0
 *
 * where m = 1 - m1, using the arithmetic-geometric mean method.
 *
 * The argument m1 is used rather than m so that the logarithmic
 * singularity at m = 1 will be shifted to the origin; this
 * preserves maximum accuracy.
 *
 * K(0) = pi/2.
 *
 * ACCURACY:
 *
 * Truncated at NBITS
 *
 */

/*
  Cephes Math Library, Release 2.0:  April, 1987
  Copyright 1984, 1987 by Stephen L. Moshier
  Revised 2018 J.N.
 */


/*	qellpk.c		*/

#include "qhead.h"

int EXPORT qellpk(Qfloatp x,Qfloatp y)
{
	int count=0;
	Qfloat a[1], b[1], c[1], temp[1];


	if( qcmpto1( x ) > 0 || x->sign != 0 )	{
		mtherr( "qellpk", DOMAIN );
		return(0);
	}
	qsub( x, qone, temp );
	//qfsqrt( temp, c );		/* c = sqrt(x)		*/
	qfsqrt(x,c);
	//qsub( temp, qone, b );	/* b = sqrt( 1 - x )	*/
	// This was doing the function for 1-x, like qellpe. At least according to
	// Matlab and Mathematica. Fixed. J.N.
	//qfsqrt( x, b );
	qfsqrt(temp,b);
	qmov( qone, a );	/* a = 1	*/

	while( ((int) a[0].exponent - (int) c[0].exponent) < (NBITS) )	{
		qsub( b, a, c );	/* c = (a - b)/2.0	*/
		c[0].exponent -= 1;
		qmul( a, b, temp );	/* temp = sqrt( a * b )	*/
		qadd( a, b, a );	/* a = (a + b)/2.0	*/
		a[0].exponent -= 1;
		qfsqrt(temp,b);
		count++;
		if (count > 600)
			break;
	}

	qdiv( a, qPi_Div_2, y );
	return(count);
}
