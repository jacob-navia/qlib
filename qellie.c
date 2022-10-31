/*							qellie.c
 *
 *	Incomplete elliptic integral of the second kind
 *
 *
 *
 * SYNOPSIS:
 *
 * int qellie( phi, m, y );
 * QELT *phi, *m, *y;
 *
 * qellie( phi, m, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *                phi
 *                 -
 *                | |
 *                |                   2
 * E(phi_\m)  =    |    sqrt( 1 - m sin t ) dt
 *                |
 *              | |
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 * ACCURACY:
 *
 * Sequence terminates at NBITS/2.
 *
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987, 1993 by Stephen L. Moshier
 * Adapted to lcc-win by Jacob Navia
 * Incomplete elliptic integral of second kind
 * Arguments are phi and m
 */

#include "qhead.h"
static Qfloat a[1], b[1], c[1], e[1], t[1], temp[1];
static Qfloat lphi[1], temp2[1];

int qellie(Qfloatp phi,Qfloatp m,Qfloatp y )
{
	int sign;
	int d;
	long long mod;
	double dmod;

	if( m->exponent < (QELT) (EXPONE - 129) )	{
		qmov( phi, y );
		return(0);
	}

	if( qcmp( m, qone ) > 0 || m->sign != 0 )	{
		mtherr( "qellie", DOMAIN );
		return(0);
	}
	qsub( m, qone, b );	/* b = sqrt( 1 - m )	*/
	if( b[0].exponent < (QELT) (EXPONE - 129) )	{
		qfsin( phi, y );
		return(0);
	}
	qmov( phi, lphi );
	if( lphi[0].sign )  sign = -1;
	else
	    sign = 0;
	lphi[0].sign = 0;		/* make positive */
	qmov( qone, a );	/* a = 1		*/
	qfsqrt( b, b );
	qfsqrt( m, c );		/* c = sqrt( m )	*/
	d = 1;			/* d = 1		*/
	qmov( qone, e );	/* e = 0		*/
	e[0].exponent = 0;
	e[0].mantissa[0] = 0;
	qftan( lphi, t );			/* t = tan( phi )	*/
	qmov( qpi, temp );	/* temp = pi/2 + phi	*/
	temp[0].exponent -= 1;
	qadd( lphi, temp, temp );
	qdiv( qpi, temp, temp );	/* mod = (phi + pi/2)/pi	*/
	dmod = qtoe( temp, NOROUNDING);
	mod = dmod;

	while( ((int) a[0].exponent - (int) c[0].exponent) < (NBITS/2) )	{
			qdiv( a, b, temp );	/* temp = b/a		*/
			qmul( t, temp, temp2 );	/* phi += arctan(t*temp) + mod*pi	*/
			qatn( temp2, temp2 );
			qadd( lphi, temp2, lphi );
			lltoq( mod, temp2 );
			qmul( temp2, qpi, temp2 );
			qadd( lphi, temp2, lphi );
			qmov( qpi, temp2 );	/* mod = (phi + pi/2)/pi	*/
			temp2[0].exponent -= 1;
			qadd( lphi, temp2, temp2 );
			qdiv( qpi, temp2, temp2 );
			dmod = qtoe( temp2, NOROUNDING);
			mod = dmod;
			qmul( t, t, temp2 );	/* t *= (1+temp)/(1-temp*t*t)	*/
			qmul( temp, temp2, temp2 );
			qsub( temp2, qone, temp2 );
			qadd( qone, temp, temp );
			qmul( t, temp, temp );
			qdiv( temp2, temp, t );
			qsub( b, a, c );	/* c = (a - b)/2.0	*/
			c[0].exponent -= 1;
			qmul( a, b, temp );	/* temp = sqrt( a * b )	*/
			qfsqrt( temp, temp );
			qadd( a, b, a );	/* a = (a + b)/2.0	*/
			a[0].exponent -= 1;
			qmov( temp, b );	/* b = temp		*/
			d += d;			/* d += d 		*/
			qfsin( lphi, temp );	/* e += c * sin(lphi)	*/
			qmul( temp, c, temp );
			qadd( e, temp, e );
		}

	qsub( m, qone, b );	/* b = 1 - m			*/
	qellpe( b, temp );	/* ellpe(b)/ellpk(b)		*/
	qellpk( b, temp2 );
	qdiv( temp2, temp, c );
	lltoq( mod, temp);	/* (arctan(t) + mod*pi)/(d*a)	*/
	qmul( temp, qpi, temp );
	qatn( t, t );
	qadd( t, temp, t );
	itoq( d, temp2 );
	qmul( temp2, a, temp );
	qdiv( temp, t, t );
	qmul( c, t, c );
	qadd( c, e, y );
	y[0].sign = sign;
	return(0);
}
