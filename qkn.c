/*							kn.c
 *
 *	Modified Bessel function, third kind, integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * int qkn( n, x, y );
 * int n;
 * QELT *x, *y;
 *
 * qkn( n, x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of the third kind
 * of order n of the argument.
 *
 * The range is partitioned into the two intervals [0,9.55] and
 * (9.55, infinity).  An ascending power series is used in the
 * low range, and an asymptotic expansion in the high range.
 *
 * ACCURACY:
 *
 * Series expansions are set to terminate at less than full
 * working precision.
*
 */

/*
Cephes Math Library Release 2.1:  November, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
*/

/*   qkn.c  */
/* Relative accuracy about 22 decimals at crossover point
   that was set for 144-bit arithmetic.  */

/*
Algorithm for Kn.
                       n-1
                   -n   -  (n-k-1)!    2   k
K (x)  =  0.5 (x/2)     >  -------- (-x /4)
 n                      -     k!
                       k=0

                    inf.                                   2   k
       n         n   -                                   (x /4)
 + (-1)  0.5(x/2)    >  {p(k+1) + p(n+k+1) - 2log(x/2)} ---------
                     -                                  k! (n+k)!
                    k=0

where  p(m) is the psi function: p(1) = -EUL and

                      m-1
                       -
      p(m)  =  -EUL +  >  1/k
                       -
                      k=1

For large x,
                                         2        2     2
                                      u-1     (u-1 )(u-3 )
K (z)  =  sqrt(pi/2z) exp(-z) { 1 + ------- + ------------ + ...}
 v                                        1            2
                                    1! (8z)     2! (8z)
asymptotically, where

           2
    u = 4 v .

*/


#include "qhead.h"

#define MAXFAC 150

void qkn(Qfloatp const integer, Qfloatp const x,Qfloatp y)
{
	long long i, lk, lj;
	long long n,nn;
	double dx;
	Qfloat k[1];
	Qfloat kf[1];
	Qfloat nk1f[1];
	Qfloat nkf[1];
	Qfloat zn[1];
	Qfloat t[1];
	Qfloat s[1];
	Qfloat z0[1];
	Qfloat z[1];
	Qfloat ans[1];
	Qfloat fn[1];
	Qfloat pn[1];
	Qfloat pk[1];
	Qfloat zmn[1];
	Qfloat t1[1];
	Qfloat t2[1];
	Qfloat tlg[1];

	qifrac(integer,&nn,k);
	if( nn < 0 )
		n = -nn;
	else
		n = nn;

	if( (x[0].sign != 0) || (x[0].exponent < 3) || (n > MAXFAC) )
	{
		mtherr( "qkn", DOMAIN );
		qclear(y);
		return;
	}


	dx = qtoe( x, DOROUNDING );
	if( dx > 40.0 ) // was 24
		goto asymp;

	qclear(ans);		 /* ans = 0.0;*/
	qmul( x, x, z0 );	/* z0 = 0.25 * x * x; */
	z0[0].exponent -= 2;
	qmov( qone, fn );	/* fn = 1.0; */
	qclear(pn);		/* pn = 0.0; */
	qmov( qone, zmn );	/* zmn = 1.0; */

	if( n > 0 )
	{
		/* compute factorial of n and psi(n) */
		qmov( qeul, pn );	/* pn = -EUL; */
		qneg(pn);
		qmov( qone, k );	/* k = 1.0; */
		for( i=1; i<n; i++ )
		{
			qdiv( k, qone, t );	/* pn += 1.0/k; */
			qadd( pn, t, pn );
			qadd( qone, k, k );	/* k += 1.0; */
			qmul( fn, k, fn );	/* fn *= k; */
		}

		qdiv( x, qtwo, zmn );		/* zmn = 2.0/x; */

		if( n == 1 )
		{
			qdiv( x, qone, ans );	/* ans = 1.0/x; */
		}
		else
			{
			lltoq( n, t );
			qdiv( t, fn, nk1f );	/* nk1f = fn/n; */
			qmov( qone, kf );	/* kf = 1.0; */
			qmov( nk1f, s );	/* s = nk1f; */
			qmov( z0, z );		/* z = -z0; */
			qneg( z );
			qmov( qone, zn );	/* zn = 1.0; */
			for( i=1; i<n; i++ ) {
				lk = (int)(n - i);		/* nk1f = nk1f/(n-i); */
				itoq( lk, t );
				qdiv( t, nk1f, nk1f );
				itoq( i, t );		/* kf = kf * i; */
				qmul( kf, t, kf );
				qmul( zn, z, zn );	/* zn *= z; */
				qmul( nk1f, zn, t );	/* t = nk1f * zn / kf; */
				qdiv( kf, t, t );
				qadd( s, t, s );	/* s += t; */
				qdiv( x, zmn, zmn );	/* zmn *= 2.0/x; */
				zmn[0].exponent += 1;
			}
			qmul( s, zmn, ans );		/* ans = s * zmn * 0.5; */
			ans[0].exponent -= 1;
		}
	}


	qmov( x, s );	/* 2 log(x/2) */
	s[0].exponent -= 1;
	qflog( s, tlg );
	tlg[0].exponent += 1;

	qmov( qeul, pk );		/* pk = -EUL; */
	qneg( pk );
	if( n == 0 ) {
		qmov( pk, pn );		/* pn = pk; */
		qmov( qone, t );	/* t = 1.0; */
	}
	else {
		lltoq( n, t );		/* pn = pn + 1.0/n; */
		qdiv( t, qone, t );
		qadd( pn, t, pn );
		qdiv( fn, qone, t );	/* t = 1.0/fn; */
	}
	qadd( pk, pn, s );		/* s = (pk+pn)*t; */
	qsub( tlg, s, s );		/* pk + pn - 2log(x/2) */
	qmul( t, s, s );
	lk = 1;		/* k = 1.0; */
	do {
		lj = lk + n;		/* t *= z0 / (k * (k+n)); */
		itoq( lj, t1 );
		itoq( lk, t2 );
		qmul( t2, t1, z );
		qdiv( z, z0, z );
		qmul( t, z, t );
		qdiv( t2, qone, z );	/* pk += 1.0/k; */
		qadd( pk, z, pk );
		qdiv( t1, qone, z );	/* pn += 1.0/(k+n); */
		qadd( pn, z, pn );
		qadd( pk, pn, z );	/* s += (pk+pn)*t; */
		qsub( tlg, z, z );	/* pk + pn - 2log(x/2) */
		qmul( z, t, z );
		qadd( s, z, s );
		lk += 1.0;
	}
	while( ((int) s[0].exponent - (int) t[0].exponent) < NBITS ); /* fabs(t/s) > MACHEP ); */

	if( n > 0 )
		qdiv( zmn, s, s );		/* s = 0.5 * s / zmn; */
	s[0].exponent -= 1;
	if( n & 1 )
		qneg( s );		/* s = -s; */
	qadd( ans, s, y );		/* ans += s; */

	return;

	/* Asymptotic expansion for Kn(x) */
	/* Converges to 1.4e-17 for x > 18.4 */

asymp:

	lk = 4 * n * n;			/* pn = 4.0 * n * n; */
	itoq( lk, pn );
	qmov( qone, pk );		/* pk = 1.0; */
	qmov( x, z0 );			/* z0 = 8.0 * x; */
	z0[0].exponent += 3;
	qmov( qone, fn );		/* fn = 1.0; */
	qmov( qone, t );		/* t = 1.0; */
	qmov( t, s );			/* s = t; */
	qmov( qone, nkf );		/* nkf = MAXNUM; */
	nkf[0].exponent += 16000;
	i = 0;
	do {
		qmul( pk, pk, t1 );	/* z = pn - pk * pk; */
		qsub( t1, pn, z );
		qmul( fn, z0, t1 );	/* t = t * z /(fn * z0); */
		qdiv( t1, z, t1 );
		qmul( t, t1, t );
		qmov( t, nk1f );	/* nk1f = fabs(t); */
		nk1f[0].sign = 0;
		qsub( nkf, nk1f, t1 );
		if( (i >= n) && (t1[0].sign == 0) ) /* nk1f > nkf */ {
			/*		printf( "qkn: i=%D, %d\n", i, t[1]-s[1] );*/
			goto adone;
		}
		qmov( nk1f, nkf );	/* nkf = nk1f; */
		qadd( s, t, s );	/* s += t; */
		qadd( qone, fn, fn );	/* fn += 1.0; */
		qadd( qtwo, pk, pk );	/* pk += 2.0; */
		i += 1;
	}
	while( ((int) s[0].exponent - (int) t[0].exponent) < NBITS ); /* fabs(t/s) > MACHEP ); */

adone:
	qdiv( x, qpi, z );	/* ans = exp(-x) * sqrt( PI/(2.0*x) ) * s; */
	z[0].exponent -= 1;
	qfsqrt( z, z );
	qfexp( x, t );
	qdiv( t, z, ans );
	qmul( s, ans, y );
}
