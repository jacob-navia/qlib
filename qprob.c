/*  qprob.c  */
/* various probability integrals
* computed via incomplete beta and gamma integrals
*/

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/

#include "qhead.h"

/* binomial distribution */
/*							qbdtr
*
 *	Binomial distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qbdtr( k, n, p, y );
* int k, n;
* QELT *p, *y;
*
 * qbdtr( k, n, p, y );
*
 * DESCRIPTION:
*
 * Returns (in y) the sum of the terms 0 through k of the Binomial
* probability density:
*
 *   k
*   --  ( n )   j      n-j
*   >   (   )  p  (1-p)
*   --  ( j )
*  j=0
*
 * The terms are not summed directly; instead the incomplete
* beta integral is employed, according to the formula
*
 * y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).
*
 * The arguments must be positive, with p ranging from 0 to 1.
*
 */


void qbdtr(Qfloatp const kk,const Qfloatp nn,const Qfloatp p,Qfloatp y)
{
	Qfloat dk[1], dn[1], dp[1];
	int li;
	long long k = qtoll(kk),n=qtoll(nn);

	if( k >= n )
	{
		qmov( qone, y );
		return;
	}

	li = k + 1;
	itoq( li, dk );	/* dk = k */

	li = n - k;
	itoq( li, dn );

	qmov( p, dp );
	qsub( dp, qone, dp );

	qincb( dn, dk, dp, y );
}

/*							qbdtrc
*
 *	Complemented binomial distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qbdtrc( k, n, p, y );
* int k, n;
* QELT *p, *y;
*
 * y = qbdtrc( k, n, p, y );
*
 * DESCRIPTION:
*
 * Returns the sum of the terms k+1 through n of the Binomial
* probability density:
*
 *   n
*   --  ( n )   j      n-j
*   >   (   )  p  (1-p)
*   --  ( j )
*  j=k+1
*
 * The terms are not summed directly; instead the incomplete
* beta integral is employed, according to the formula
*
 * y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).
*
 * The arguments must be positive, with p ranging from 0 to 1.
*
 */

void qbdtrc(const int k,const int n,Qfloatp const p,const Qfloatp y)
{
	Qfloat dk[1], dn[1];
	int li;

	if( k < 0 )
	{
		qmov( qone, y );
		return;
	}
	if( k == n )
	{
		qclear( y );
		return;
	}
	li = k + 1;
	itoq( li, dk );	/* dk = k */

	li = n - k;
	itoq( li, dn );

	qincb( dk, dn, p, y );
}

/*							qbdtri
*
 *	Inverse binomial distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qbdtri( k, n, y, p );
* int k, n;
* QELT *p, *y;
*
 * qbdtri( k, n, y, p );
*
 * DESCRIPTION:
*
 * Finds the event probability p such that the sum of the
* terms 0 through k of the Binomial probability density
* is equal to the given cumulative probability y.
*
 * This is accomplished using the inverse beta integral
* function and the relation
*
 * 1 - p = incbi( n-k, k+1, y ).
*
 */

void qbdtri(int k,int n,Qfloatp y,Qfloatp p)
{
	Qfloat dk[1], dn[1];
	int li;


	if( (n <= k) || (k < 0) )
	{
		qclear( y );
		return;
	}
	li = k + 1;
	itoq( li, dk );	/* dk = k */

	li = n - k;
	itoq( li, dn );

	beta_distribution_invQ( dn, dk, y, p );
	qsub( p, qone, p );
}


/*							qchdtr
*
 *	Chi-square distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qchdtr( df, x, y );
* QELT *df, *x, *y;
*
 * qchdtr( df, x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the area under the left hand tail (from 0 to x)
* of the Chi square probability density function with
* v degrees of freedom.
*
 *
 *                                  inf.
*                                    -
*                        1          | |  v/2-1  -t/2
*  P( x | v )   =   -----------     |   t      e     dt
*                    v/2  -       | |
*                   2    | (v/2)   -
*                                   x
*
 * where x is the Chi-square variable.
*
 * The incomplete gamma integral is used, according to the
* formula
*
 *	y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).
*
 *
 * The arguments must both be positive.
*
 */

void qchdtr( Qfloatp df, Qfloatp x, Qfloatp y )
{
	Qfloat a[1], b[1];

	qmov( df, a );
	qmov( x, b );
	a[0].exponent -= 1;
	b[0].exponent -= 1;
	qigam( a, b, y );
}

/*							qchdtc
*
 *	Complemented Chi-square distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qchdtc( df, x, y );
* QELT df[], x[], y[];
*
 * qchdtc( df, x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the area under the right hand tail (from x to
* infinity) of the Chi square probability density function
* with v degrees of freedom:
*
 *
 *                                  inf.
*                                    -
*                        1          | |  v/2-1  -t/2
*  P( x | v )   =   -----------     |   t      e     dt
*                    v/2  -       | |
*                   2    | (v/2)   -
*                                   x
*
 * where x is the Chi-square variable.
*
 * The incomplete gamma integral is used, according to the
* formula
*
 *	y = chdtr( v, x ) = igamc( v/2.0, x/2.0 ).
*
 *
 * The arguments must both be positive.
*
 */

void qchdtc(Qfloatp df, Qfloatp x, Qfloatp y)
{
	Qfloat a[1], b[1];

	qmov( df, a );
	qmov( x, b );
	a[0].exponent -= 1;
	b[0].exponent -= 1;
	qigamc( a, b, y );
}


/*							qchdti
*
 *	Inverse of complemented Chi-square distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qchdti( df, y, x );
* QELT *df, *x, *y;
*
 * qchdti( df, y, x );
*
 *
 *
 *
 * DESCRIPTION:
*
 * Finds the Chi-square argument x such that the integral
* from x to infinity of the Chi-square density is equal
* to the given cumulative probability y.
*
 * This is accomplished using the inverse gamma integral
* function and the relation
*
 *    x/2 = igami( df/2, y );
*
 *
 *
 *
 * ACCURACY:
*
 * See igami.c.
*
 * ERROR MESSAGES:
*
 *   message         condition      value returned
* chdtri domain   y < 0 or y > 1        0.0
*                     v < 1
*
 */

void qchdti(Qfloatp df,Qfloatp y, Qfloatp x )
{
	Qfloat a[1];

	qmov( df, a );
	a[0].exponent -= 1;
	qigami( a, y, x );
	x[0].exponent += 1;
}



/*							qfdtr
*
 *	F distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qfdtr( ia, ib, x, y );
* int ia, ib;
* QELT *x, *y;
*
 * qfdtr( ia, ib, x, y );
*
 * DESCRIPTION:
*
 * Returns the area from zero to x under the F density
* function (also known as Snedcor's density or the
* variance ratio density).  This is the density
* of x = (u1/df1)/(u2/df2), where u1 and u2 are random
* variables having Chi square distributions with df1
* and df2 degrees of freedom, respectively.
*
 * The incomplete beta integral is used, according to the
* formula
*
 *	P(x) = incbet( df1/2, df2/2, (df1*x/(df2 + df1*x) ).
*
 *
 * The arguments a and b are greater than zero, and x is
* nonnegative.
*
 */

void qfdtr(int ia,int ib,Qfloatp x,Qfloatp y )
{
	Qfloat a[1], b[1], u[1], v[1], w[1];
	int li;

	li = ia;
	itoq( li, a );
	li = ib;
	itoq( li, b );

	qmov( a, w );		/* ax/(b+ax)  */
	qmul( x, w, w );
	qadd( b, w, u );
	qdiv( u, w, w );
	qmov( a, u );
	u[0].exponent -= 1;
	qmov( b, v );
	v[0].exponent -= 1;
	qincb( u, v, w, y ); /* incbet( a/2.0, b/2.0, (ax/(b+ax)) )  */
}


/*							qfdtrc
*
 *	Complemented F distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qfdtrc( ia, ib, x, y );
* int ia, ib;
* QELT x[], y[];
*
 * qfdtrc( ia, ib, x, y );
*
 * DESCRIPTION:
*
 * Returns the area from x to infinity under the F density
* function (also known as Snedcor's density or the
* variance ratio density).
*
 *
 *                      inf.
*                       -
*              1       | |  a-1      b-1
* 1-P(x)  =  ------    |   t    (1-t)    dt
*            B(a,b)  | |
*                     -
*                      x
*
 *
 * The incomplete beta integral is used, according to the
* formula
*
 *	P(x) = incbet( df2/2, df1/2, (df2/(df2 + df1*x) ).
*
 */

void qfdtrc(int ia,int ib,Qfloatp x, Qfloatp y)
{
	Qfloat a[1], b[1], u[1], v[1], w[1];
	int li;

	li = ia;
	itoq( li, a );
	li = ib;
	itoq( li, b );

	/* b/(b + ax)  */
	qmov( a, u );
	qmul( x, u, u );
	qadd( b, u, u );
	qdiv( u, b, w );

	qmov( a, u );
	u[0].exponent -= 1;
	qmov( b, v );
	v[0].exponent -= 1;
	qincb( v, u, w, y ); /* incbet( b/2.0, a/2.0, (b/(b+ax)) )  */
}

/*							qfdtri
*
 *	Inverse of complemented F distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qfdtri( ia, ib, y, x );
* int ia, ib;
* QELT x[], y[];
*
 * qfdtri( ia, ib, y, x );
*
 * DESCRIPTION:
*
 * Finds the F density argument x such that the integral
* from x to infinity of the F density is equal to the
* given probability p.
*
 * This is accomplished using the inverse beta integral
* function and the relations
*
 *      z = incbi( df2/2, df1/2, p )
*      x = df2 (1-z) / (df1 z).
*
 * Note: the following relations hold for the inverse of
* the uncomplemented F distribution:
*
 *      z = incbi( df1/2, df2/2, p )
*      x = df2 z / (df1 (1-z)).
*
 */

void qfdtri(int ia,int ib,Qfloatp y,Qfloatp x )
{
	Qfloat a[1], b[1], u[1], v[1], w[1];
	int li;

	li = ia;
	itoq( li, a );
	li = ib;
	itoq( li, b );

	qmov( a, u );
	u[0].exponent -= 1;
	qmov( b, v );
	v[0].exponent -= 1;
	beta_distribution_invQ( v, u, y, w ); /* incbi( b/2.0, a/2.0, y )  */

	/* x = (b - bw)/aw  */
	qmul( b, w, u );
	qsub( u, b, u );
	qmul( a, w, v );
	qdiv( v, u, x );
}


/*							qgdtr
*
 *	Gamma distribution function
*
 *
 *
 * SYNOPSIS:
*
 * int qgdtr( a, b, x, y );
* QELT *a, *b, *x, *y;
*
 * qgdtr( a, b, x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the integral from zero to x of the gamma probability
* density function:
*
 *
 *                x
*        b       -
*       a       | |   b-1  -at
* y =  -----    |    t    e    dt
*       -     | |
*      | (b)   -
*               0
*
 *  The incomplete gamma integral is used, according to the
* relation
*
 * y = igam( b, ax ).
*
 */

void qgdtr(Qfloatp a,Qfloatp b,Qfloatp x,Qfloatp y )
{
	Qfloat w[1];

	qmul( a, x, w );
	qigam( b, w, y );  /* igam( b, a * x )  */
}

/*							qgdtrc
*
 *	Complemented gamma distribution function
*
 *
 *
 * SYNOPSIS:
*
 * int qgdtrc( a, b, x, y );
* QELT *a, *b, *x, *y;
*
 * qgdtrc( a, b, x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the integral from x to infinity of the gamma
* probability density function:
*
 *
 *               inf.
*        b       -
*       a       | |   b-1  -at
* y =  -----    |    t    e    dt
*       -     | |
*      | (b)   -
*               x
*
 *  The incomplete gamma integral is used, according to the
* relation
*
 * y = igamc( b, ax ).
*
 */

void qgdtrc(Qfloatp a,Qfloatp b,Qfloatp x,Qfloatp y )
{
	Qfloat w[1];

	qmul( a, x, w );
	qigamc( b, w, y );  /* igamc( b, a * x )  */
}

/*							qnbdtr
*
 *	Negative binomial distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qnbdtr( k, n, p, y );
* int k, n;
* QELT *p, *y;
*
 * qnbdtr( k, n, p, y );
*
 * DESCRIPTION:
*
 * Returns the sum of the terms 0 through k of the negative
* binomial distribution:
*
 *   k
*   --  ( n+j-1 )   n      j
*   >   (       )  p  (1-p)
*   --  (   j   )
*  j=0
*
 * In a sequence of Bernoulli trials, this is the probability
* that k or fewer failures precede the nth success.
*
 * The terms are not computed individually; instead the incomplete
* beta integral is employed, according to the formula
*
 * y = nbdtr( k, n, p ) = incbet( n, k+1, p ).
*
 * The arguments must be positive, with p ranging from 0 to 1.
*
 */

void qnbdtr(const int k,const int n,const Qfloatp p,Qfloatp y )
{
	Qfloat dk[1], dn[1];
	int li;

	if( k == 0 )
	{
		qmov( qone, y );
		return;
	}

	li = k + 1;
	itoq( li, dk );

	li = n;
	itoq( li, dn );

	qincb( dn, dk, p, y );
}

/*							qnbdtc
*
 *	Complemented negative binomial distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qnbdtc( k, n, p, y );
* int k, n;
* QELT *p, *y;
*
 * qnbdtc( k, n, p, y );
*
 * DESCRIPTION:
*
 * Returns the sum of the terms k+1 to infinity of the negative
* binomial distribution:
*
 *   inf
*   --  ( n+j-1 )   n      j
*   >   (       )  p  (1-p)
*   --  (   j   )
*  j=k+1
*
 * The terms are not computed individually; instead the incomplete
* beta integral is employed, according to the formula
*
 * y = nbdtrc( k, n, p ) = incbet( k+1, n, 1-p ).
*
 * The arguments must be positive, with p ranging from 0 to 1.
*
 */

void qnbdtc( int k,int n,Qfloatp p,Qfloatp y)
{
	Qfloat dk[1], dn[1], w[1];
	int li;

	if( k == 0 )
	{
		qmov( qone, y );
		return;
	}

	li = k + 1;
	itoq( li, dk );

	li = n;
	itoq( li, dn );

	qsub( p, qone, w );
	qincb( dk, dn, w, y );
}



/*							qpdtr
*
 *	Poisson distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qpdtr( k, m, y );
* int k;
* QELT *m, *y;
*
 * qpdtr( k, m, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the sum of the first k terms of the Poisson
* distribution:
*
 *   k         j
*   --   -m  m
*   >   e    --
*   --       j!
*  j=0
*
 * The terms are not summed directly; instead the incomplete
* gamma integral is employed, according to the relation
*
 * y = pdtr( k, m ) = igamc( k+1, m ).
*
 * The arguments must both be positive.
*
 */

int qPoissonDistribution(int k,Qfloatp m, Qfloatp y)
{
	Qfloat v[1];
	int li;


	li= k+1;
	itoq( li, v );

	qigamc( v, m, y );  /* igamc( k+1, m ) */
	return 0;
}

/*							qpdtrc
*
 *	Complemented poisson distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qpdtrc( k, m, y );
* int k;
* QELT *m, *y;
*
 * qpdtrc( k, m, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the sum of the terms k+1 to infinity of the Poisson
* distribution:
*
 *  inf.       j
*   --   -m  m
*   >   e    --
*   --       j!
*  j=k+1
*
 * The terms are not summed directly; instead the incomplete
* gamma integral is employed, according to the formula
*
 * y = pdtrc( k, m ) = igam( k+1, m ).
*
 * The arguments must both be positive.
*
 */

int qPoissonDistributionComp(const int k,const Qfloatp m,Qfloatp y)
{
	Qfloat v[1];
	int li;


	li= k+1;
	itoq( li, v );

	qigam( v, m, y );  /* igam( k+1, m ) */
	return 0;
}


/*							qpdtri
*
 *	Inverse Poisson distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qpdtri( k, y, m );
* int k;
* QELT *m, *y;
*
 * qpdtri( k, y, m );
*
 *
 *
 *
 * DESCRIPTION:
*
 * Finds the Poisson variable x such that the integral
* from 0 to x of the Poisson density is equal to the
* given probability y.
*
 * This is accomplished using the inverse gamma integral
* function and the relation
*
 *    m = igami( k+1, y ).
*
 */

int qPoissonDistributionInv(int k,Qfloatp y,Qfloatp m)
{
	Qfloat v[1];
	int li;


	li= k+1;
	itoq( li, v );

	qigami( v, y, m );
	return 0;
}

