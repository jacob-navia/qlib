/*							qexp2.c
 *
 *	Check routine for base 2 exponential function
 *
 *
 *
 * SYNOPSIS:
 *
 * int qexp2( x, y );
 * QELT *x, *y;
 *
 * qexp2( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns 2 raised to the x power.
 *
 *        x      ln 2  x      x ln 2
 * y  =  2  = ( e     )   =  e
 *
 * Cephes Math Library Release 2.3:  March, 1995
 * Copyright 1984, 1995 by Stephen L. Moshier
 *
 * Adapted to lcc-win by J.N.
*/


#include "qhead.h"

void qexp2(Qfloatp x,Qfloatp y)
{
	Qfloat a[1];

	qmul( x, qlog2, a );
	qfexp( a, y );
}

/*							qlogtwo
 *
 *	Base 2 logarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * int qlogtwo(Qfloatp x, Qfloatp y );
 *
 * qlogtwo( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base 2 logarithm of x.
 *
 */

void qlogtwo(Qfloatp x,Qfloatp y)
{

	qflog( x, y );
	qmul( y, qinv_log2, y );
}
