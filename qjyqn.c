/*						qjyqn.c		*/
/* Auxiliary function for Hankel's asymptotic expansion */

/* Jn(x) = sqrt(2/(pi x)) [ P(n,x) cos X  -  Q(n,x) sin X ]
* Yn(x) = sqrt(2/(pi x)) [ P(n,x) sin X  +  Q(n,x) cos X ]
*
 * where arg of sine and cosine = X = x - (0.5n + 0.25)*PI.
* We solve this for Pn(x):
* Jn(x) cos X  +  Yn(x) sin X  =  sqrt( 8/(pi x)) Pn(x)
*
 * Series expansions are set to terminate at less than full
* working precision.
*
 */

#include "qhead.h"

/*extern double MAXNUM, PI;*/
/*#define EUL 0.57721566490153286060*/

int qjyqn(Qfloatp qnn, Qfloatp x, Qfloatp y )
{
	Qfloat nfac[1], nm1fac[1], f[1], a[1], psi1[1];
	Qfloat g[1], h[1], jn[1], yn[1], psin[1], z[1], s[1];
	Qfloat qn[1], p[1], q[1], t[1];
	Qfloat tt[1], u[1], odd[1], tlast[1];
	int temp;
	int i, k, n, sign, kpn;
	double dn;

	dn = qtoe( qnn, DOROUNDING );
	n = dn;

	dn = qtoe( x, DOROUNDING );
	if( dn > 64.0 )		/* use asymptotic expansion if x > 32.5 */
		goto hank;

	if( n < 0 )	{
			n = -n;
			if( (n & 1) == 0 )	/* -1**n */		sign = 1;
			else
				sign = -1;
		}
	else
		sign = 1;


	if( x[0].exponent == 0 )	{
			mtherr("qpn", OVERFLOW );
			qinfin(yn);
			return( 0 );
		}

	/*							qpn.c 2		*/


	/* factorial of n */
	qmov( qone, nfac);	/*nfac = 1.0;*/
	qclear( f );		/*f = 0;*/
	for( i=0; i<n-1; i++ )
	{
		qadd( qone, f, f );	/*f += 1.0;*/
		qmul( f, nfac, nfac );	/*nfac *= f;*/
	}
	qmov( nfac, nm1fac );		/*nm1fac = nfac;*/
	qadd( qone, f, f );		/*f += 1.0;*/
	qmul( nfac, f, nfac );		/*nfac *= f;*/
	/*printf("nfac %.4E\n", nfac );*/
	qdiv( nfac, qone, a );		/*a = 1.0/nfac;*/

	/* psi function */
	qclear( psi1 );			/*psi1 = 0;*/
	qclear( psin );			/*psin = 0;*/
	for( i=1; i<=n; i++ )
	{
		itoq( i, z );
		qdiv( z, qone, z );
		qadd( psin, z, psin );	/*psin += 1.0/i;*/
	}
	/*printf("psin %.4E\n", psin );*/
	qmul( x, x, z );
	z[0].exponent -= 2;
	z[0].sign = -1;			/*z = -x*x/4.0;*/
	qdiv( nfac, psin, s );		/*s = psin/nfac;*/
	kpn = n+1;
	qmov( qone, f );		/* k = 1 */
	itoq( kpn, g );
	do
	    {
		qdiv( f, qone, h );
		qadd( psi1, h, psi1 );	/*psi1 += 1.0/k;*/
		qdiv( g, qone, h );
		qadd( psin, h, psin );	/*psin += 1.0/kpn;*/
		qmul( f, g, h );
		qdiv( h, z, h );
		qmul( h, a, a );	/*a *= z/(k*kpn);*/
		qadd( psi1, psin, h );
		qmul( h, a, h );
		/*
		if( s[0] != h[0] )
		{
		if( s[1] > bt[1] )
		qmov( s, bt );
		if( h[1] > bt[1] )
		qmov( h, bt );
		}
		*/
		qadd( s, h, s );	/*s += a*(psi1+psin);*/
		qadd( qone, f, f );	/*k += 1;*/
		qadd( qone, g, g );	/*kpn += 1;*/
	}
	while( ((int) s[0].exponent - (int) a[0].exponent) < NBITS/2 );
	/* estimate of cancellation error */
	/*
	qdiv( s, bt, bt );
	bt[1] -= 144;
	qtoasc( bt, asy, 4 );
	printf( "yn est error %s\n", asy );
	*/
	/*printf("infinite %.4E\n", s);*/

	/*							qpn.c 3		*/
	/* finite sum */

	qclear( f );			/*f = 0;*/
	if( n > 0 ){
			z[0].sign = 0;			/*z = -z;*/
			qmov( nm1fac, a );		/*a = nm1fac;*/
			kpn = n - 1;
			qmov( a, f );			/*f = a;*/
			itoq( kpn, g );
			qmov( qone, h );		/* k = 1 */
			for( k=1; k<n; k++ )
			{
				qmul( g, h, nm1fac );
				qdiv( nm1fac, z, nm1fac );
				qmul( a, nm1fac, a );	/*a *= z/(k*kpn);*/
				qadd( f, a, f );	/*f += a;*/
				qsub( qone, g, g );	/*kpn -= 1;*/
				qadd( qone, h, h );
			}
		}
	/*printf("finite %.4E\n", f);*/
	/* x/2**n */
	qmov( x, a );
	a[0].exponent -= 1;			/*a = x/2;*/
	qmov( qone, z );		/*z = 1.0;*/
	for( i=1; i<=n; i++ )
		qmul( z, a, z );	/*z *= a;*/
	/*printf("x/2**n %.4E\n", z );*/
	/* combine the terms */
	/*s = 2.0*(log(a)+EUL)*jn(n,x) - z*s - f/z;*/
	qflog( a, g );
	qadd( g, qeul, g );
	i = n;
	itoq( i, h );
	bessel_J( h, x, jn );
	qmul( jn, g, g );
	g[0].exponent += 1;

	qdiv( z, f, f );
	qsub( f, g, g );

	qmul( z, s, f );
	qsub( f, g, g );

	qdiv( qpi, g, yn );		/*s /= PI;*/
	if( sign < 0 )	yn[0].sign = ~yn[0].sign;		/*return(s*sign);*/
	goto findp;

	/*		jvpq.c
	* Hankel's asymptotic expansion for Bessel functions Jv(x)
	* Note: does not converge to 144 bit accuracy
	* for x less than about 51.5.
	*/


hank:

	if( n != 0 )	{
			i = n;
			itoq( i, u );
			qmul( u, u, u );	/* u = 4.0 * n * n */
			u[0].exponent += 2;
		}
	else
		qclear(u);

	qmov( x, z );
	z[0].exponent += 3;		/* z = 8.0 * x */
	qclear(qn);		/* n = 0.0     */
	qmov( qone, t );	/* t = 1.0 */
	qclear(p);		/* p = 0.0 */
	qclear(q);		/* q = 0.0 */
	qmov( qone, odd );	/* odd = -1.0 */
	odd[0].sign = -1;
	sign = 1;
	qmov( qone, tlast );

	do
	    {
		if( sign > 0 )		qadd( t, p, p );	/* p += t */
		else
		    qsub( t, p, p );	/* p -= t */

		qadd( qone, qn, qn );		/* n += 1.0 */
		qadd( qtwo, odd, odd );		/* odd += 2.0 */
		/* t *= (u - odd*odd)/(n * z) */
		qmul( odd, odd, tt );
		qsub( tt, u, tt );
		qmul( t, tt, t );
		qmul( qn, z, tt );
		qdiv( tt, t, t );


		if( sign > 0 )		qadd( t, q, q );	/* q += t */
		else
		    qsub( t, q, q );	/* q -= t */
		qadd( qone, qn, qn );		/* n += 1.0 */
		qadd( qtwo, odd, odd );		/* odd += 2.0 */
		/* t *= (u - odd*odd)/(n * z) */
		qmul( odd, odd, tt );
		qsub( tt, u, tt );
		qmul( tt, t, t );
		qmul( qn, z, tt );	/* /(n*z) */
		qdiv( tt, t, t );
		temp = t[0].sign;
		t[0].sign = 0;
		if( qcmp(t, tlast) > 0 ) /* stop if terms get bigger */		goto done;
		qmov( t, tlast );
		t[0].sign = temp;
		if( qn[0].exponent > (QELT) (EXPONE + 6) )		goto done;
		sign = -sign;
	}
	while( ((int) q[0].exponent - (int) t[0].exponent) < NBITS/2 );

done:

	/*qmov( p, y );*/
	qmov( q, y );
	/*
	qtoasc( tlast, asy, 4 );
	printf( "qjypq tlast %s\n", asy );
	*/
	return(0);
findp:

	/* Jn(x) = sqrt(2/(pi x)) [ P(n,x) cos X  -  Q(n,x) sin X ]
	* Yn(x) = sqrt(2/(pi x)) [ P(n,x) sin X  +  Q(n,x) cos X ]
	*
 * where arg of sine and cosine = X = x - (0.5n + 0.25)*PI.
	* We solve this for Pn(x):
	* Jn(x) cos X  +  Yn(x) sin X  =  sqrt( 8/(pi x)) Pn(x)
	*/

	qmov( qone, tt );
	tt[0].exponent -= 2;		/* 0.25 */
	if( n != 0 )	{
			i = n;
			itoq( i, u );
			u[0].exponent -= 1;		/* 0.5n  */
			qadd( u, tt, tt );
		}
	else
		qclear(u);

	qmul( qpi, tt, tt );	/* times pi */
	qsub( tt, x, tt );	/* subtracted from x */

	/* P(x) = Jn cosX + Yn sinX */
	/*
	qfcos( tt, f );
	qmul( jn, f, f );
	qfsin( tt, g );
	qmul( yn, g, g );
	qadd( f, g, f );
	*/

	/* Q(x) = Yn cosX - Jn sinX */
	qfcos( tt, f );
	qmul( yn, f, f );
	qfsin( tt, g );
	qmul( jn, g, g );
	qsub( g, f, f );

	/* constant factor = sqrt( pi x/2 ) */
	qmul( x, qpi, tt );
	tt[0].exponent -= 1;
	qfsqrt( tt, odd );
	qmul( odd, f, y );
	return( 0 );
}
