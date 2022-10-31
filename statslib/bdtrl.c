/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/
#include <errno.h>
extern double incbet ( double, double, double );
extern double incbi ( double, double, double );
extern double pow ( double, double );
extern double expm1 ( double );
extern double log1p ( double );

/*							bdtrcl()
*
*	Complemented binomial distribution
*
*
*
* SYNOPSIS:
*
* int k, n;
* double p, y, bdtrcl();
*
* y = bdtrcl( k, n, p );
*
*
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
*
*
* ACCURACY:
*
* See incbet.c.
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* bdtrcl domain     x<0, x>1, n<k       0.0
*/
double binomial_c(int k,int n,double p )
{
	double dk, dn;

	if( (p < 0.0L) || (p > 1.0L) )
		goto domerr;
	if( k < 0 )
		return( 1.0L );

	if( n < k )
	{
domerr:
		errno = ERANGE;
		return( 0.0L );
	}

	if( k == n )
		return( 0.0L );
	dn = n - k;
	if( k == 0 )
	{
		if( p < .01L )
			dk = -expm1( dn * log1p(-p) );
		else
			dk = 1.0L - pow( 1.0L-p, dn );
	}
	else
	{
		dk = k + 1;
		dk = incbet( dk, dn, p );
	}
	return( dk );
}



/*							bdtrl.c
*
 *	Binomial distribution
*
*
*
* SYNOPSIS:
*
* int k, n;
* double p, y, bdtrl();
*
* y = bdtrl( k, n, p );
*
*
*
* DESCRIPTION:
*
* Returns the sum of the terms 0 through k of the Binomial
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
*
*
* ACCURACY:
*
* Tested at random points (k,n,p) with a and b between 0
* and 10000 and p between 0 and 1.
*    Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE      0,10000      3000       1.6e-14     2.2e-15
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* bdtrl domain        k < 0            0.0
*                     n < k
*                     x < 0, x > 1
*
*/
double binomial(int k,int n,double p )
{
	double dk, dn, q;

	if( (p < 0.0L) || (p > 1.0L) )
		goto domerr;
	if( (k < 0) || (n < k) )
	{
domerr:
		errno = ERANGE;
		return( 0.0L );
	}

	if( k == n )
		return( 1.0L );

	q = 1.0L - p;
	dn = n - k;
	if( k == 0 )
	{
		dk = pow( q, dn );
	}
	else
	{
		dk = k + 1;
		dk = incbet( dn, dk, q );
	}
	return( dk );
}


/*							bdtril()
*
 *	Inverse binomial distribution
*
 *
 *
 * SYNOPSIS:
*
 * int k, n;
* double p, y, bdtril();
*
 * p = bdtril( k, n, y );
*
 *
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
 * ACCURACY:
*
 * See incbi.c.
* Tested at random k, n between 1 and 10000.  The "domain" refers to p:
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE       0,1        3500       2.0e-15     8.2e-17
*
 * ERROR MESSAGES:
*
 *   message         condition      value returned
* bdtril domain     k < 0, n <= k         0.0
*                  x < 0, x > 1
*/

/*								bdtr() */


#ifdef __LCC__
double __declspec(naked) binomial_inv(int k,int n,double y )
{
}
#endif
double bdtri(int k,int n, double y)
{
	double dk, dn, p;

	if( (y < 0.0L) || (y > 1.0L) || (k < 0) || (n <= k) ) {
		errno = ERANGE;
		return(0);
	}

	dn = n - k;
	if( k == 0 )
	{
		if( y > 0.8L )
			p = -expm1( log1p(y-1.0L) / dn );
		else
		    p = 1.0L - pow( y, 1.0L/dn );
	}
	else
	    {
		dk = k + 1;
		p = incbet( dn, dk, y );
		if( p > 0.5 )
			p = incbi( dk, dn, 1.0L-y );
		else
		    p = 1.0 - incbi( dn, dk, y );
	}
	return( p );
}
