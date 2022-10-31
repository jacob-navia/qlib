/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

#include <errno.h>
int mtherr (char *, int);
extern double igaml ( double, double );
extern double igamcl ( double, double );
extern double igamil ( double, double );

/*							pdtrcl()
*
 *	Complemented poisson distribution
*
 *
 *
 * SYNOPSIS:
*
 * int k;
* double m, y, pdtrcl();
*
 * y = pdtrcl( k, m );
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
 *
 *
 * ACCURACY:
*
 * See igam.c.
*
 */
#ifdef __LCC__
double __declspec(naked) poisson_distribution_c(int k,double m)
{
}
#endif
double pdtrcl(int k,double m)
{
	double v;

	if( (k < 0) || (m <= 0.0) )
	{
		mtherr( "pdtrcl", EDOM );
		return( 0.0 );
	}
	v = k+1;
	return( igaml( v, m ) );
}



/*							pdtrl.c
*
 *	Poisson distribution
*
 *
 *
 * SYNOPSIS:
*
 * int k;
* double m, y, pdtrl();
*
 * y = pdtrl( k, m );
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
 *
 *
 * ACCURACY:
*
 * See igamc().
*
 */
#ifdef __LCC__
double __declspec(naked) poisson_distribution(int k,double m)
{
}
#endif
double pdtrl(int k,double m )
{
	double v;

	if( (k < 0) || (m <= 0.0) )
	{
		mtherr( "pdtrl", EDOM );
		return( 0.0 );
	}
	v = k+1;
	return( igamcl( v, m ) );
}


/*							pdtril()
*
 *	Inverse Poisson distribution
*
 *
 *
 * SYNOPSIS:
*
 * int k;
* double m, y, pdtrl();
*
 * m = pdtril( k, y );
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
* pdtri domain    y < 0 or y >= 1       0.0
*                     k < 0
*
 */
#ifdef __LCC__
double __declspec(naked) poisson_distribution_inv(int k,double y)
{
}
#endif
double pdtril(int k,double y )
{
	double v;

	if( (k < 0) || (y < 0.0) || (y >= 1.0) )
	{
		mtherr( "pdtril", EDOM );
		return( 0.0 );
	}
	v = k+1;
	v = igamil( v, y );
	return( v );
}
