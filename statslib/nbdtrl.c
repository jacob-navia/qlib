/*
Cephes Math Library Release 2.3:  January,1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

#include <errno.h>
int mtherr (char *, int);
extern double incbet( double, double, double );
extern double pow( double, double );
extern double incbi( double, double, double );

/*							nbdtrcl.c
*
 *	Complemented negative binomial distribution
*
 *
 *
 * SYNOPSIS:
*
 * int k, n;
* double p, y, nbdtrcl();
*
 * y = nbdtrcl( k, n, p );
*
 *
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
 *
 *
 * ACCURACY:
*
 * See incbetl.c.
*
 */
#ifdef __LCC__
double __declspec(naked) binomial_neg_distribution_c(int k,int n,double p)
{
}
#endif
double nbdtrc(int k,int n,double p )
{
	double dk, dn;

	if( (p < 0.0) || (p > 1.0) )
		goto domerr;
	if( k < 0 )
	{
domerr:
		mtherr( "nbdtrl", EDOM );
		return( 0.0 );
	}
	dn = n;
	if( k == 0 )
		return( 1.0 - pow( p, dn ) );

	dk = k+1;
	return( incbet( dk, dn, 1.0 - p ) );
}



/*							nbdtrl.c
*
 *	Negative binomial distribution
*
 *
 *
 * SYNOPSIS:
*
 * int k, n;
* double p, y, nbdtrl();
*
 * y = nbdtrl( k, n, p );
*
 *
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
 *
 *
 * ACCURACY:
*
 * Tested at random points (k,n,p) with k and n between 1 and 10,000
* and p between 0 and 1.
*
 * arithmetic   domain     # trials      peak         rms
*    Absolute error:
*    IEEE      0,10000     10000       9.8e-15     2.1e-16
*
 */
#ifdef __LCC__
double __declspec(naked) binomial_neg_distribution(int k,int n,double p)
{
}
#endif
double nbdtr(int k,int n, double p )
{
	double dk, dn;

	if( (p < 0.0) || (p > 1.0) )
		goto domerr;
	if( k < 0 )
	{
domerr:
		mtherr( "nbdtrl", EDOM );
		return( 0.0 );
	}
	dn = n;
	if( k == 0 )
		return( pow( p, dn ) );

	dk = k+1;
	return( incbet( dn, dk, p ) );
}


/*							nbdtril
*
 *	Functional inverse of negative binomial distribution
*
 *
 *
 * SYNOPSIS:
*
 * int k, n;
* double p, y, nbdtril();
*
 * p = nbdtril( k, n, y );
*
 *
 *
 * DESCRIPTION:
*
 * Finds the argument p such that nbdtr(k,n,p) is equal to y.
*
 * ACCURACY:
*
 * Tested at random points (a,b,y), with y between 0 and 1.
*
 *               a,b                     Relative error:
* arithmetic  domain     # trials      peak         rms
*    IEEE     0,100
* See also incbil.c.
*/
#ifdef __LCC__
double __declspec(naked) binomial_neg_distribution_inv(int k,int n,double p)
{
}
#endif
double nbdtri(int k,int n,double p )
{
	double dk, dn, w;

	if( (p < 0.0) || (p > 1.0) )
		goto domerr;
	if( k < 0 ) {
domerr:
		mtherr( "nbdtr", EDOM );
		return( 0.0 );
	}
	dk = k+1;
	dn = n;
	w = incbi( dn, dk, p );
	return( w );
}
