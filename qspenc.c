/*							qspenc.c
*
*	Dilogarithm
*
*
*
* SYNOPSIS:
*
* int qspenc( x, y );
* QELT *x, *y;
*
* qspenc( x, y );
*
*
*
* DESCRIPTION:
*
* Computes the integral
*
*                    x
*                    -
*                   | | log t
* spence(x)  =  -   |   ----- dt
*                 | |   t - 1
*                  -
*                  1
*
 * for x >= 0.  A power series gives the integral in
* the interval (0.5, 1.5).  Transformation formulas for 1/x
* and 1-x are employed outside the basic expansion range.
*
 *
 *
 */

/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1985, 1987, 1989, 1999 by Stephen L. Moshier
*/

#include "qhead.h"

int qspenss( Qfloat *, Qfloat * );

int qspenc(Qfloat x[],Qfloat y[])
{
	Qfloat xxx[1], a[1], w[1], t[1];
	double xx;
	int flag;
	long ll;

	xx = qtoe( x, DOROUNDING );
	if( xx == 1.0 )
	{
		qclear( y );
		return( 0 );
	}
	if( xx == 0.0 )
	{
		qmul( qpi, qpi, y );	/* y = pi*pi/6 */
		qmov( qone, a );
		ll = 6;
		itoq(ll,a);
		/*
		a[3] = 0140000;
		a[1] = EXPONE + 2;
		*/
		qdiv( a, y, y );
	}

	qmov( x, xxx );
	flag = 0;
	if( xx > 2.0 )
	{
		xx = 1.0/xx;
		qdiv( xxx, qone, xxx );	/* x = 1.0/x */
		flag |= 2;
	}
	if( xx > 1.5 )
	{
		qdiv( xxx, qone, w );	/* w = (1.0/x) - 1.0 */
		qsub( qone, w, w );
		flag |= 2;
	}
	else if( xx < 0.5 )
	{
		qmov( xxx, w );		/* w = -x */
		w[0].sign = ~w[0].sign;
		flag |= 1;
	}
	else
		qsub( qone, xxx, w );	/* w = x - 1.0 */

	qspenss( w, y );

	if( flag & 1 )
	{
		/*	y = (PI * PI)/6.0  - log(x) * log(1.0-x) - y */
		qmul( qpi, qpi, t );	/* t = pi*pi/6 */
		qmov( qone, a );
		ll = 6;
		itoq (ll,a);
		/*
		a[3] = 0140000;
		a[1] = EXPONE + 2;
		*/
		qdiv( a, t, t );

		qflog( xxx, a );
		qsub( xxx, qone, w );
		qflog( w, w );
		qmul( w, a, a );
		qsub( a, t, t );
		qsub( y, t, y );
	}
	if( flag & 2 )
	{
		qflog( xxx, t );		/* t = log(x)	*/
		/* y = -0.5 * t * t  -  y */
		qmul( t, t, t );
		t[0].exponent -= 1;
		qadd( y, t, y );
		y[0].sign = ~y[0].sign;
	}

	return( 0 );
}


/* Power series for dilogarithm */

int qspenss(Qfloat w[],Qfloat y[])
{
	Qfloat a[1], k[1], t[1];

	qmov( qone, a );		/* a = 1.0 */
	qclear( k );			/* k = 0.0 */
	qclear( y );			/* y = 0.0 */

	do {
		qmul( a, w, a );	/*	a *= w		*/
		qadd( k, qone, k );	/*	k += 1.0	*/
		qmul( k, k, t );	/*	y -= a/(k*k)	*/
		qdiv( t, a, t );
		qsub( t, y, y );
		qmul( a, w, a );	/*	a *= w		*/
		qadd( k, qone, k );	/*	k += 1.0	*/
		qmul( k, k, t );	/*	y += a/(k*k)	*/
		qdiv( t, a, t );
		qadd( y, t, y );
	}
	while( ( y[0].exponent - a[0].exponent) < NBITS );

	/*printf( "%.3E ", k );*/
	return 0;
}

