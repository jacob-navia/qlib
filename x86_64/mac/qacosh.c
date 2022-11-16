/*							qacosh.c
 *
 *	Inverse hyperbolic cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * int qacosh( x, y )
 * QELT *x, *y;
 *
 * qacosh( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * acosh(x)  =  log( x + sqrt( (x-1)(x+1) ).
 *
 * Cephes Math Library Release 2.3:  March, 1995
 * Copyright 1984, 1995 by Stephen L. Moshier
   Adapted to lcc-win by jacob navia 1998-2008
*/

#include "qhead.h"


void qacosh(Qfloatp x,Qfloatp y)
{
	Qfloat a[1];

	if( qcmp( x, qone ) < 0 )
	{
		qclear(y);
		return;
	}
	if( (int)exponent(x) > (int) (EXPONE + NBITS))
	{
		qflog( x, y );
		qadd( qlog2, y, y );
		return;
	}
	qmul( x, x, a );	/* sqrt( x**2 - 1 )	*/
	qsub( qone, a, a );
	qfsqrt( a, a );
	qadd( x, a, a );
	qflog( a, y );		/* log( x + sqrt(...)	*/
}
