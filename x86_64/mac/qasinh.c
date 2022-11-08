/*							qasinh.c
*
 *	Inverse hyperbolic sine
*
 *
 *
 * SYNOPSIS:
*
 * int qasinh( x, y );
* QELT *x, *y;
*
 * qasinh( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns inverse hyperbolic sine of argument.
*
 *     asinh(x) = log( x + sqrt(1 + x*x) ).
*
 * For very large x, asinh(x) = log x  +  log 2.
*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/


#include "qhead.h"

void EXPORT qasinh(Qfloatp xxx,Qfloatp y)
{
	Qfloat x[1], a[1], b[1], s[1], xx[1], z[1];
	int j, n, i, sign;

	sign = signof(xxx);
	qmov( xxx, x );
	setpositive(x);
	if( (long)exponent(x) < (QELT) (EXPONE - 4) )
		goto cfrac;

	if( ((int)exponent(x) - (int)EXPONE) >= ((int)MAXEXP - (int)EXPONE)/2 )
	{
		qflog( x, y );
		qadd( qlog2, y, y );
		if( sign )
			qneg( y );
		return;
	}
	qmul( x, x, a );	/* sqrt( x**2 + 1 )	*/
	qadd( qone, a, a );
	qfsqrt( a, a );
	qadd( x, a, a );
	qflog( a, y );		/* log( x + sqrt(...)	*/
	if( sign )
		qneg( y );
	return;

cfrac:


	qmul( x, x, z );	/* z = x * x */
	qadd( qone, z, a );	/* a = sqrt( z + 1.0 ) */
	qfsqrt( a, a );

	i = NBITS/6;     /* 8  */
	n = 2*i+1;     /* 17 for x < .05, convergence to 1e-17 */
	qclear( s );

	do	{
		j = i * (i-1);
		itoq( j, b );	/* b = j * z */
		qmul( b, z, b );
		itoq( n, xx );     /* s = b/(s+n) */
		qadd( s, xx, xx );
		qdiv( xx, b, s );
		n -= 2;
		itoq( n, xx );     /* s = b/(s+n) */
		qadd( s, xx, xx );
		qdiv( xx, b, s );
		n -= 2;
		i -= 2;
	}
	while( n > 1 );

	/* a * x / (1+s) */
	qadd( qone, s, s );
	qdiv( s, a, a );
	qmul( a, x, y );
	if( sign )
		qneg( y );
}
