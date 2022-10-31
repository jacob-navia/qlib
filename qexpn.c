/*							qexpn.c
*
 *		Exponential integral En
*
 *
 *
 * SYNOPSIS:
*
 * int qexpn( n, x, y );
* int n;
* QELT *x, *y;
*
 * qexpn( n, x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Evaluates the exponential integral
*
 *                 inf.
*                   -
*                  | |   -xt
*                  |    e
*      E (x)  =    |    ----  dt.
*       n          |      n
*                | |     t
*                 -
*                  1
*
 *
 * Both n and x must be nonnegative.
*
 *
 * ACCURACY:
*
 * 105 decimal digits
*

 * Cephes Math Library Release 1.1:  March, 1985
* Copyright 1985 by Stephen L. Moshier */

/* Exponential integral	*/
#include "qhead.h"
//int qtol (qfloat *x);

void Qexpn(Qfloatp n,Qfloatp x,Qfloatp y)
{
	long long N = qtoll(n);
	qexpn(N,x,y);
}

void qexpn(int n,Qfloatp x,Qfloatp yy)
{
	int i, k;
	int ln;
	Qfloat ans[1], r[1], t[1], yk[1], xk[1], qn[1];
	Qfloat pk[1], pkm1[1], pkm2[1], qk[1], qkm1[1], qkm2[1];
	Qfloat psi[1], z[1];

	if( n < 0 )	{
		mtherr("qexpn", DOMAIN );
		goto overf;
	}

	if( signof(x) != 0 )	{
		mtherr("qexpn", DOMAIN );
		goto overf;
	}

	if( exponent(x) < 3 )	{
		if( n < 2 )		goto overf;
		else {
			ln = n - 1;
			itoq( ln, ans );
			qinv( ans, yy );
			return;
		}
	}

	if( n == 0 )	{
		qfexp( x, ans );		/* exp(-x)/x */
		qmul( ans, x, ans );
		qinv( ans, yy );
		return;
	}

	/*							expn.c	*/
	/*		Expansion for large n		*/
	/*
	if( n > 5000 )
	{
	xk = x + n;
	yk = 1.0 / (xk * xk);
	t = n;
	ans = yk * t * (6.0 * x * x  -  8.0 * t * x  +  t * t);
	ans = yk * (ans + t * (t  -  2.0 * x));
	ans = yk * (ans + t);
	ans = (ans + 1.0) * exp( -x ) / xk;
	goto done;
	}
	*/

	if( exponent(x) > (QELT) (EXPONE+1) )	goto cfrac;

	/*							expn.c	*/

	/*		Power series expansion		*/


	qflog( x, psi );		/* psi = -EUL - log(x) */
	qneg( psi );
	qsub( qeul, psi, psi );

	for( i=1; i<n; i++ )
	{
		ln = i;
		itoq( ln, qn );
		qinv( qn, qn );	/* psi = psi + 1.0/i */
		qadd( qn, psi, psi );
	}

	qmov( x, z );		/* z = -x */
	qneg( z );
	qclear( xk );		/* xk = 0.0 */
	qmov( qone, yk );	/* yk = 1.0 */
	ln = n;
	itoq( ln, qn );
	qsub( qn, qone, pk );	/* pk = 1.0 - n */
	if( n == 1 )	qclear( ans );	/* ans = 0.0 */
	else
	    qinv( pk, ans );	/* ans = 1.0/pk */

	do
	    {
		qincr( xk, xk );	/* xk += 1.0 */
		qdiv( xk, z, qn );	/* yk *= z/xk */
		qmul( qn, yk, yk );
		qincr( pk, pk );	/* pk += 1.0 */
		if( exponent(pk) > 10 )		{
				qdiv( pk, yk, t );	/* ans += yk/pk */
				qadd( t, ans, ans );
			}
		else
			qmov( qone, t );	/* t = 1.0 */
	}
	while( ( exponent(t) - exponent(ans)) > -70 );

#ifdef NOTUSED
	temp = qtoe( xk, NOROUNDING );
	k = temp;	/* k = xk */
#endif
	ln = n;
	itoq( ln, t );	/* t = n */
	/* ans = (powi(z, n-1) * psi / gamma(t)) - ans */
	qgamma( t, t );
	qflog( x, qn );
	ln = n - 1;
	itoq( ln, yk );
	qmul( yk, qn, qn );
	qfexp( qn, qn );
	if( ((n-1) & 1) != 0 )	qneg( qn );
	qmul( psi, qn, qn );
	qdiv( t, qn, qn );
	qsub( ans, qn, ans );
	goto done;


	/*							expn.c	*/
	/*		continued fraction		*/
cfrac:
	k = 1;
	qmov( qone, pkm2 );	/* pkm2 = 1.0 */
	qmov( x, qkm2 );	/* qkm2 = x   */
	qmov( qone, pkm1 );	/* pkm1 = 1.0 */
	ln = n;
	itoq( ln, qn );
	qadd( qn, x, qkm1 );	/* qkm1 = x + n */
	qdiv( qkm1, pkm1, ans );	/* ans = pkm1/qkm1 */

	do
	    {
		k += 1;
		if( k & 1 )		{
				qmov( qone, yk );	/* yk = 1.0 */

				ln = n + (k-1)/2;	/* xk = n + (k-1)/2 */
				itoq( ln, xk );
			}
		else
		{
			qmov( x, yk );		/* yk = x */
			ln = k/2;		/* xk = k/2 */
			itoq( ln, xk );
		}
		qmul( yk, pkm1, qn );	/* pk = pkm1 * yk  +  pkm2 * xk */
		qmul( xk, pkm2, pk );
		qadd( qn, pk, pk );
		qmul( yk, qkm1, qn );	/* qk = qkm1 * yk  +  qkm2 * xk */
		qmul( xk, qkm2, qk );
		qadd( qn, qk, qk );
		if( qk[0].exponent > 2 )		{
				qdiv( qk, pk, r );	/* r = pk/qk */
				qsub( r, ans, t );	/* t = abs( (ans - r)/r ) */
				qmov( r, ans );		/* ans = r */
			}
		else
			qmov( qone, t );	/* t = 1.0 */
		qmov( pkm1, pkm2 );		/* pkm2 = pkm1 */
		qmov( pk, pkm1 );		/* pkm1 = pk   */
		qmov( qkm1, qkm2 );		/* qkm2 = qkm1 */
		qmov( qk, qkm1 );		/* qkm1 = qk   */
		if( pk[0].exponent > (QELT) (EXPONE + 64) )	{
				pkm2[0].exponent -= 64;
				pkm1[0].exponent -= 64;
				qkm2[0].exponent -= 64;
				qkm1[0].exponent -= 64;
		}
	}
	while( ( exponent(t) -  exponent(r)) > -NBITS );

	qfexp( x, qn );
	qdiv( qn, ans, ans );	/* ans *= exp( -x ) */

done:
	qmov( ans, yy );
	return;

overf:
	mtherr( "qexpn", OVERFLOW );
}
