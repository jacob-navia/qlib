/*							qlog10.c
*
 *	Common logarithm
*
 *
 *
 * SYNOPSIS:
*
 * int qlog10( x, y );
* QELT *x, *y;
*
 * qlog10( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns base 10, or common, logarithm of x.
*
 * log  (x) = log  (e) log (x)
*    10         10       e
*
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

/*	qlog10	*/

#include "qhead.h"


void qlog10(Qfloatp x,Qfloatp y)
{
	qflog( x, y );
	qmul( ql10e, y, y );
}
