
/*								chdtr() */


/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/
#include <errno.h>
int mtherr (char *, int);
extern double igamc( double, double );
extern double igamma( double, double );
extern double igami( double, double );

/*							chdtrc()
*
 *	Complemented Chi-square distribution
*
 *
 *
 * SYNOPSIS:
*
 * double v, x, y, chdtrc();
*
 * y = chdtrc( v, x );
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
 *
 *
 * ACCURACY:
*
 * See igamc().
*
 * ERROR MESSAGES:
*
 *   message         condition      value returned
* chdtrc domain  x < 0 or v < 1        0.0
*/
#ifdef __LCC__
double __declspec(naked) chi_sqr_distribution_c(double df,double x)
{
}
#endif
double chdtrc(double df,double x)
{

	if( (x < 0.0L) || (df < 1.0L) )
	{
		mtherr( "chdtrcl", EDOM );
		return(0.0L);
	}
	return( igamc( 0.5L*df, 0.5L*x ) );
}



/*							chdtrl.c
*
 *	Chi-square distribution
*
 *
 *
 * SYNOPSIS:
*
 * double df, x, y, chdtr();
*
 * y = chdtr( df, x );
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
 *
 *
 * ACCURACY:
*
 * See igam().
*
 * ERROR MESSAGES:
*
 *   message         condition      value returned
* chdtr domain   x < 0 or v < 1        0.0
*/
#ifdef __LCC__
double __declspec(naked) chi_sqr_distribution(double df,double x)
{
}
#endif
double chdtr(double df,double x)
{

	if( (x < 0.0L) || (df < 1.0L) )
	{
		mtherr( "chdtrl", EDOM );
		return(0.0L);
	}
	return( igamma( 0.5L*df, 0.5L*x ) );
}


/*							chdtri()
*
 *	Inverse of complemented Chi-square distribution
*
 *
 *
 * SYNOPSIS:
*
 * double df, x, y, chdtri();
*
 * x = chdtri( df, y );
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
#ifdef __LCC__
double __declspec(naked) chi_sqr_distribution_cinv(double df,double y)
{
}
#endif
double chdtri(double df,double y )
{
	double x;

	if( (y < 0.0L) || (y > 1.0L) || (df < 1.0L) )
	{
		mtherr( "chdtril", EDOM );
		return(0.0L);
	}

	x = igami( 0.5L * df, y );
	return( 2.0L * x );
}
