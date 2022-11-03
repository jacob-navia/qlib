/*							qk0.c
*
 * Modified Bessel function K of order 0
*
 * K (z)  =  - (ln(z/2) + eul) I (z)
*  0                           0
*              2                    2   2                     2   3
*             z / 4               (z /4)                    (z /4)
*           + -----  + (1 + 1/2) -------- + (1 + 1/2 + 1/3) ------- + ...
*                 2                    2                          2
*             (1!)                 (2!)                       (3!)
*
 * (AMS55 #9.6.13)
*
 * Series expansions are set to terminate at less than full
* working precision.
*/

#include "qhead.h"

int qinzero(Qfloatp x,Qfloatp y )
{
	Qfloat t[1], u[1], z[1], k[1], ans[1];
	long i, ln;
	double dx;
	Qfloat nkf[1];
	Qfloat nk1f[1];
	Qfloat pk[1];
	Qfloat fn[1];
	Qfloat s[1];
	Qfloat pn[1];
	Qfloat z0[1];
	Qfloat t1[1];

	dx = qtoe( x, NOROUNDING );
	if( dx > 40.0 )
		goto asymp;

	qmul( x, x, z );	/*z = -x * x / 4.0;		*/

	if( z[0].exponent < 3 )		/* x = 0, n = 0 is special case	*/
	{
		qmov( qone, y );
		return(0);
	}

	if( z[0].exponent > 1 )
		z[0].exponent -= 2;

	qmov( qone, u );
	qmov( u, ans );		/*ans = u;*/
	qmov( qone, k );	/*k = 1.0;*/

	while( ((int) ans[0].exponent - (int) u[0].exponent) < NBITS/2 )	/* 70 */
	{
		qmov( k, t );
		qmul( t, k, t );	/*u *= z / (k * (n+k));*/
		qdiv( t, z, t );
		qmul( t, u, u );
		if( (int) u[0].exponent <= 0 ) {
			//printf( "qin overflow\n" );
			break;
		}
		qadd( u, ans, ans );	/*ans += u;*/
		qadd( qone, k, k );	/*k += 1.0;*/
	}

	/* ans *= exp( n * log( x/2.0 ) );*/

	qmov( ans, y );
	return 0;


	/* Asymptotic expansion for In(x) */

asymp:

	dx = 0.0;		/* ln = n */
	ln = dx;
	dx = 4.0 * dx * dx;		/* pn = 4.0 * n * n; */
	etoq( dx, pn );
	qmov( qone, pk );		/* pk = 1.0; */
	qmov( x, z0 );			/* z0 = -8.0 * x; */
	z0[0].exponent += 3;
	qneg( z0 );
	qmov( qone, fn );		/* fn = 1.0; */
	qmov( qone, t );		/* t = 1.0; */
	qmov( t, s );			/* s = t; */
	qmov( qone, nkf );		/* nkf = MAXNUM; */
	nkf[0].exponent += 16000;
	i = 0;
	do
	    {
		qmul( pk, pk, t1 );	/* z = pn - pk * pk; */
		qsub( t1, pn, z );
		qmul( fn, z0, t1 );	/* t = t * z /(fn * z0); */
		qdiv( t1, z, t1 );
		qmul( t, t1, t );
		qmov( t, nk1f );	/* nk1f = fabs(t); */
		nk1f[0].sign = 0;
		qsub( nkf, nk1f, t1 );
		if( (i >= ln) && (t1[0].sign == 0) ) /* nk1f > nkf */
		{
			//printf( "qin: i=%ld, %d\n", i, t[0].exponent-s[0].exponent );
			goto adone;
		}
		qmov( nk1f, nkf );	/* nkf = nk1f; */
		qadd( s, t, s );	/* s += t; */
		qadd( qone, fn, fn );	/* fn += 1.0; */
		qadd( qtwo, pk, pk );	/* pk += 2.0; */
		i += 1;
	}
	while( ((int) s[0].exponent - (int) t[0].exponent) < NBITS );	/* fabs(t/s) > MACHEP ); */

adone:

	qmul( invSqrt2pi, s, y );	/* ans = s exp(x) / sqrt( 2 pi x ) */

	qfexp( x, t );
	qfsqrt( x, z );
	qdiv( z, t, t );
	qmul( t, y, y );
	return 0;
}
void qk0(Qfloatp x,Qfloatp y)
{
	static Qfloat t[1], u[1], z[1], k[1], p[1];

	qmul( x, x, z );	/*z = x * x / 4.0;		*/
	z[0].exponent -= 2;

	qmov( qone, p );	/* psi function */
	qmov( qtwo, k );	/*k = 1.0;*/

	qmov( z, y );		/* y = u;*/
	qmov( z, u );
	qmov( z, t );

	while( ((int) y[0].exponent - (int) t[0].exponent) < NBITS )
	{
		qdiv( k, qone, t );
		qadd( t, p, p );	/* psi function */
		qmul( k, k, t );
		qdiv( t, u, u );	/* (k!)**2	*/
		qmul( z, u, u );	/* z**k		*/
		qmul( p, u, t );	/* psi * u	*/
		qadd( t, y, y );	/*ans += u;*/
		qadd( qone, k, k );	/*k += 1.0;*/
	}


	qmov( x, t );		/* log(x/2) + eul	*/
	t[0].exponent -= 1;
	qflog( t, u );
	qadd( qeul, u, u );

	qinzero(  x, t );
	qmul( t, u, u );
	qsub( u, y, y );
}
