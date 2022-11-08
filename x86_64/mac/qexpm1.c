/*							qexpm1.c
*
 *	Exponential function check routine
*
 *
 *
 * SYNOPSIS:
*
 * int qexpm1( x, y );
* QELT *x, *y;
*
 * qexpm1( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns e (2.71828...) raised to the x power, minus 1.
*
 */

/*
Cephes Math Library Release 2.3:  January, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

#include "qhead.h"

/* C1 + C2 = ln 2 */
static Qfloat C1[1] = { EXPONE-1,0xb17217f7};


void qexpm1(Qfloatp x,Qfloatp y )
{
	Qfloat num[1], den[1], x2[1];


	/* goto use_exp; */

	qmov (C1, num);
	num[0].exponent -= 1;
	if (qcmp(x, num) > 0)
		goto use_exp;
	qmov(x, den);
	qneg(den);
	if (qcmp(den, num) > 0)
		goto use_exp;

	qmov(x, x2);
	x2[0].exponent -= 1;		/* x/2				*/
	qtanh( x2, x2 );	/* tanh( x/2 )			*/
	/* 2 tanh / (1 - tanh) */
	qmov( x2, num );	/* 2 tanh			*/
	num[0].exponent += 1;
	qneg( x2 );
	qadd( x2, qone, den );	/* 1 - tanh			*/
	qdiv( den, num, y );	/* (2 tanh)/(1 - tanh)	*/
	return ;

use_exp:

	qfexp(x, num);
	qsub(qone, num, y);
}
