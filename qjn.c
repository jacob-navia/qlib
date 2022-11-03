/*							qjn.c
*
 *	Bessel function of noninteger order
*
 *
 *
 * SYNOPSIS:
*
 * int qjn( v, x, y );
* QELT *v, *x, *y;
*
 * qjn( v, x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns Bessel function of order v of the argument,
* where v is real.  Negative x is allowed if v is an integer.
*
 * Two expansions are used: the ascending power series and the
* Hankel expansion for large v.  If v is not too large, it
* is reduced by recurrence to a region of better accuracy.
*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
 * jn.c		1 Dec 83
* Bessel function of order n
*/
#undef DEBUG
#define DEBUG 0
#define ERRCK 1
#include <stdio.h>
#include <math.h>
#include "qhead.h"
double fabs(double);
double floor(double);

Qfloat t[1],u[1],hankzz[1],k[1],ans[1],hankc[1],hanks[1],j[1],m[1];
Qfloat hankpp[1],hankqq[1],rans[1],ru[1];
/* Cephes former qjn */
void bessel_J(Qfloatp const nn,Qfloatp const xx,Qfloatp y)
{
	double dx, dn, an, du;
	int i, sign;
	int bt;
	Qfloat n[1], x[1];

	bt = 0;
	qmov( nn, n );
	qmov( xx, x );
	dn = qtoe( n, DOROUNDING );
	dx = qtoe( x, DOROUNDING );
	sign = 1;
	an = fabs( dn );
	if( an == floor(an) )
	{
		i = an - 16384.0 * floor( an/16384.0 );
		if( n[0].sign != 0 )
		{
			if( i & 1 )
				sign = -sign;
			n[0].sign = 0;
		}
		if( x[0].sign != 0 )
		{
			if( i & 1 )
				sign = -sign;
			x[0].sign = 0;
		}
	}
	else
	{
		if( x[0].sign != 0 )
		{
			mtherr( "qjv", DOMAIN );
			qclear( y );
			goto done;
		}
	}

	dx = fabs(dx);


	if( dx > 81.0 ) {
				/*	if( dx > 0.95*an )*/
				if( dx > 1.4*an )
				{
					qhank( n, x, y );
					goto done;
				}
				if( dx > 0.7 * an )
				{
					du = 3.6 * sqrt(dx);
					/*		du = 0.8 * dx;*/
					etoq( du, ru );
					ru[0].sign = 0;
					qfloor( n, m );
					qsub( m, n, m );
					qfloor( ru, ru );
					qadd( ru, m, ru );
					if( n[0].sign == 0 )
					{
						qrecur( n, x, ru, rans );
					}
					else
						{
						qmov( ru, m );
						qmov( n, ru );
						qrecur( m, x, ru, rans );
						qmov( m, ru );
					}
					if( rans[0].exponent == 0 )
					{
						qclear( y );
						goto done;
					}
					qhank( ru, x, y );
					if( n[0].sign == 0 )
						qdiv( rans, y, y );
					else
						qmul( rans, y, y );
					goto done;
				}
			}

	qmul( x, x, hankzz );	/*z = -x * x / 4.0;		*/
	/* x = 0, n = 0 is special case	*/

	if( hankzz[0].exponent < 3 )
	{
		if( n[0].exponent < 3 )
			qmov( qone, y );
		else
			qclear( y );
		goto done;
	}

	hankzz[0].exponent -= 2;
	hankzz[0].sign = -1;

	/*                   inf      2   k
	*                 v  -    (-z /4)
	*  J (z)  =  (z/2)   >  ------------
	*   v                -       -
	*                   k=0  k! | (v+k+1)
	*/
	qmov( n, ans );
	qadd( qone, ans, ans );
	if( (n[0].exponent < 3) || (qcmp( qone, n ) == 0) )
		qmov( qone, u );
	else
	{
		qgamma( ans, u );
		qdiv( u, qone, u );		/*u = 1.0/gamma(n+1);*/
	}
	qmov( u, ans );		/*ans = u;*/
	qmov( qone, k );	/*k = 1.0;*/

	while( u[0].exponent > (QELT) (ans[0].exponent - NBITS)
			|| u[0].exponent > (QELT) (qone[0].exponent - NBITS) )
	{
		qadd( n, k, t );
		qmul( t, k, t );	/*u *= z / (k * (n+k));*/
		qdiv( t, hankzz, t );
		qmul( t, u, u );
		/* remember largest term summed */
#if ERRCK
		if( u[0].exponent > bt )
			bt = u[0].exponent;
		if( ans[0].exponent > bt )
			bt = ans[0].exponent;
#endif
		qadd( u, ans, ans );	/*ans += u;*/
		qadd( qone, k, k );	/*k += 1.0;*/
	}

	/* estimate cancellation error */
#if ERRCK
	i = bt - ans[0].exponent;
	if( i > NBITS/2
			|| DEBUG )
		printf( "qjn pseries: %d bits cancellation\n", i );
#endif

	/* ans *= exp( n * log( x/2.0 ) );*/

	if( n[0].exponent < 3 )
	{
		qmov( ans, y );
	}
	else
	{
		qmov( x, t );
		t[0].exponent -= 1;
		qflog( t, u );
		qmul( u, n, u );
		qfexp( u, t );
		qmul( ans, t, y );
	}

#if DEBUG
	qsub( y, yh, yh );
	du = qtoe( yh, DOROUNDING );
	printf( "qjn - qhank = %.5e\n", du );
#endif
done:
	if( sign < 0 )
		y[0].sign = ~y[0].sign;
}


/* Hankel's asymptotic expansion
* for large x.
* AMS55 #9.2.5.
*/

void qhank(Qfloatp const n,Qfloatp const x,Qfloatp y)
{
	int bt;
	int flag, sign, nsum, i;
	double dconv;

	bt = 0;
	nsum = 0;
	qmul( n, n, m );	/* m = 4.0*n*n;*/
	m[0].exponent += 2;
	qmov( qone, j );	/* j = 1.0;*/
	qmov( x, hankzz );	/* z = 8.0 * x;*/
	hankzz[0].exponent += 3;
	qmov( qone, k );	/* k = 1.0;*/
	qmov( qone, hankc );	/* hankc = 1.0;*/
	qsub( qone, m, u );	/* u = (m - 1.0)/z;*/
	qdiv( hankzz, u, u );
	qmov( u, hanks );	/* hanks = u;*/
	sign = 1;
	qmov( qone, ans );	/* conv = 1.0;*/
	flag = 0;
	qmov( qone, t );	/* t = 1.0;*/

	while( t[0].exponent > (qone[0].exponent - NBITS)
			|| u[0].exponent > (qone[0].exponent - NBITS) )
	{
		qadd( qtwo, k, k );	/* k += 2.0;*/
		qadd( qone, j, j );	/* j += 1.0;*/
		sign = -sign;
		qmul( k, k, t );	/* u *= (m - k * k)/(j * z);*/
		qsub( t, m, t );
		qdiv( j, t, t );
		qdiv( hankzz, t, t );
		qmul( t, u, u );
		if( sign < 0 )		/* hankc += sign * u;*/
			qsub( u, hankc, hankc );
		else
			qadd( u, hankc, hankc );
		/* remember largest term summed */
#if ERRCK
		if( u[0].exponent > bt )
			bt = u[0].exponent;
		if( hankc[0].exponent > bt )
			bt = hankc[0].exponent;
#endif
		/*	printf( "Hank P: %.5E %.5E", u, p ); */
		qadd( qtwo, k, k );	/* k += 2.0;*/
		qadd( qone, j, j );	/* j += 1.0;*/
		qmul( k, k, t );	/* u *= (m - k * k)/(j * z);*/
		qsub( t, m, t );
		qdiv( j, t, t );
		qdiv( hankzz, t, t );
		qmul( t, u, u );
		if( sign < 0 )		/* q += sign * u;*/
			qsub( u, hanks, hanks );
		else
			qadd( u, hanks, hanks );
		/* remember largest term summed */
#if ERRCK
		if( u[0].exponent > bt )
			bt = u[0].exponent;
		if( hanks[0].exponent > bt )
			bt = hanks[0].exponent;
#endif
		/*	printf( " Q: %.5E %.5E\n", u, q ); */
		qdiv( hankc, u, t );	/* t = fabs(u/p);*/
		t[0].sign = 0;
		if( qcmp( t, ans ) < 0 )	/* ( t < conv )*/
		{
			qmov( t, ans );	/* conv = t;*/
			qmov( hanks, hankqq );	/* qq = hanks; */
			qmov( hankc, hankpp );	/* pp = hankc; */
			flag = 1;
			nsum += 1;
		}
		/* stop if the terms start getting larger */
		else
			{
			if( flag != 0 )
			{
				goto hank1;
			}
		}
	}

hank1:

	/* estimate cancellation error */
#if ERRCK
	i = bt - hankpp[0].exponent;
	if(i > NBITS/2
			|| DEBUG )
	{
		dconv = qtoe( n, DOROUNDING );
		printf( "qhank(%.5e,", dconv );
		dconv = qtoe( x, DOROUNDING );
		printf( "%.5e): ", dconv );
		printf( "%d bits cancellation after %d terms\n", i, nsum );
	}
#endif

#if DEBUG
	dconv = qtoe( ans, DOROUNDING );
	printf( "qhank: last term / sum = %.4E\n", dconv );
#endif

	qmov( n, t );	/* u = x - (0.5*n + 0.25) * PI;*/
	qmul( qhalf, t, t );
	qmov( qone, hanks );
	hanks[0].exponent -= 2;
	qadd( hanks, t, t );
	qmul( qpi, t, t );
	qsub( t, x, u );

	/* t = sqrt( 2.0/(PI*x) ) * ( pp * cos(u) - qq * sin(u) ); */
	qmul( qpi, x, t );
	qdiv( t, qtwo, t );
	qfsqrt( t, hankzz );
	qfsin( u, hanks );
	qfcos( u, hankc );
	qmul( hankc, hankpp, k );
	qmul( hanks, hankqq, m );
	qsub( m, k, k );
	qmul( k, hankzz, y );
}



/* Reduce the order by backward recurrence.
* AMS55 #9.1.27 and 9.1.73.
*/
void qrecur(Qfloatp n,Qfloatp x,Qfloatp newn,Qfloatp ans)
{
	int ctr;
#if DEBUG
	double da, db;
#endif
	Qfloat pkm2[1], pkm1[1], pk[1], pkp1[1];
	Qfloat qkm2[1], qkm1[1], qk[1], xk[1];
	Qfloat yk[1], r[1], kf[1], k1[1]; 

	/* continued fraction for Jn(x)/Jn-1(x)  */

	/* fstart: */

#if DEBUG
	da = qtoe( n, DOROUNDING );
	db = qtoe( newn, DOROUNDING );
	printf( "qn = %.6e, qnewn = %.6e, qcfrac = ", da, db );
#endif

	qclear( pkm2 );
	qmov( qone, qkm2 );
	qmov( x, pkm1 );
	qadd( n, n, qkm1 );
	qmul( x, x, xk );   /* xk = -x * x; */
	qneg( xk );
	qmov( qkm1, yk );
	qmov( qone, ans );
	ctr = 0;
    r[0].exponent=0;
	do	{
		qadd( qtwo, yk, yk );
		/*	pk = pkm1 * yk +  pkm2 * xk; */
		qmul( pkm1, yk, pk );
		qmul( pkm2, xk, t );
		qadd( t, pk, pk );
		/*	qk = qkm1 * yk +  qkm2 * xk;*/
		qmul( qkm1, yk, qk );
		qmul( qkm2, xk, t );
		qadd( t, qk, qk );
		qmov( pkm1, pkm2 );
		qmov( pk, pkm1 );
		qmov( qkm1, qkm2 );
		qmov( qk, qkm1 );
		if( qk[0].exponent != 0 ) {
			qdiv( qk, pk, r );
        }
		if( r[0].exponent != 0 ) {
			/*		t = fabs( (ans - r)/r ); */
			qsub( r, ans, t );
			qdiv( r, t, t );
			t[0].sign = 0;
			qmov( r, ans );
		}
		else
			qmov( qone, t );

		if( ++ctr > 1000 ) {
			printf( "qrecur: continued fraction did not converge\n" );
			goto done;
		}
		/*
		if( t < MACHEP )
		goto done;
		*/

		if( pk[0].exponent > (qone[0].exponent+NBITS) )
		{
			pkm2[0].exponent -= NBITS;
			pkm1[0].exponent -= NBITS;
			qkm2[0].exponent -= NBITS;
			qkm1[0].exponent -= NBITS;
		}
	}
	while( t[0].exponent > (qone[0].exponent - NBITS) );

done:

#if DEBUG
	da = qtoe( ans, DOROUNDING );
	printf( "%.6e\n", da );
#endif

	/* Change n to n-1 if n < 0 and the continued fraction is small
	*/
	/*
	if( nflag )
	{
	if( fabs(ans) < 0.125 )
	{
	nflag = -1;
	*n = *n - 1.0;
	goto fstart;
	}
	}
	*/

	qmov( newn, kf );

	/* backward recurrence
	*              2k
	*  J   (x)  =  --- J (x)  -  J   (x)
	*   k-1         x   k         k+1
	*/

	qmov( qone, pk );
	qdiv( ans, qone, pkm1 );
	qsub( qone, n, k1 );
	qadd( k1, k1, r );
	do {
		/*	pkm2 = (pkm1 * r  -  pk * x) / x;*/
		qmul( pkm1, r, pkm2 );
		qmul( pk, x, t );
		qsub( t, pkm2, pkm2 );
		qdiv( x, pkm2, pkm2 );
		qmov( pk, pkp1 );
		qmov( pkm1, pk );
		qmov( pkm2, pkm1 );
		qsub( qtwo, r, r );
		/*
		t = fabs(pkp1) + fabs(pk);
		if( (k > (kf + 2.5)) && (fabs(pkm1) < 0.25*t) ) {
		k1 -= 1.0;
		t = x*x;
		pkm2 = ( (r*(r+2.0)-t)*pk - r*x*pkp1 )/t;
		pkp1 = pk;
		pk = pkm1;
		pkm1 = pkm2;
		r -= 2.0;
		}
		*/
		qsub( qone, k1, k1 );
		qadd( qhalf, kf, t );
	}
	while( qcmp( k1, t ) > 0 );

	/* Take the larger of the last two iterates
	* on the theory that it may have less cancellation error.
	*/
	/*
	if( cancel ) {
	if( (kf >= 0.0) && (fabs(pk) > fabs(pkm1)) ) {
	k += 1.0;
	pkm2 = pk;
	}
	}
	*/
	qmov( k1, newn );
#if DEBUG
	da = qtoe( k1, DOROUNDING );
	db = qtoe( pkm2, DOROUNDING );
	printf( "qnewn %.6e qrans %.6e\n", da, db );
#endif
	qmov( pkm2, ans );
}
