/*							stdtrl.c
*
 *	Student's t distribution
*
 *
 *
 * SYNOPSIS:
*
 * double p, t, stdtrl();
* int k;
*
 * p = stdtrl( k, t );
*
 *
 * DESCRIPTION:
*
 * Computes the integral from minus infinity to t of the Student
* t distribution with integer k > 0 degrees of freedom:
*
 *                                      t
*                                      -
*                                     | |
*              -                      |         2   -(k+1)/2
*             | ( (k+1)/2 )           |  (     x   )
*       ----------------------        |  ( 1 + --- )        dx
*                     -               |  (      k  )
*       sqrt( k pi ) | ( k/2 )        |
*                                   | |
*                                    -
*                                   -inf.
*
* Relation to incomplete beta integral:
*
 *        1 - stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )
* where
*        z = k/(k + t**2).
*
 * For t < -1.6, this is the method of computation.  For higher t,
* a direct method is derived from integration by parts.
* Since the function is symmetric about t=0, the area under the
* right tail of the density is found by calling the function
* with -t instead of t.
*
* ACCURACY:
*
 * Tested at random 1 <= k <= 100.  The "domain" refers to t.
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE     -100,-1.6    10000       5.7e-18     9.8e-19
*    IEEE     -1.6,100     10000       3.8e-18     1.0e-19
*/

/*							stdtril.c
*
 *	Functional inverse of Student's t distribution
*
 *
 *
 * SYNOPSIS:
*
 * double p, t, stdtril();
* int k;
*
 * t = stdtril( k, p );
*
 *
 * DESCRIPTION:
*
 * Given probability p, finds the argument t such that stdtrl(k,t)
* is equal to p.
*
* ACCURACY:
*
 * Tested at random 1 <= k <= 100.  The "domain" refers to p:
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE       0,1        3500       4.2e-17     4.1e-18
*/


/*
Cephes Math Library Release 2.3:  January, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

#include <errno.h>
#include <math.h>
int mtherr (char *, int);

extern double PI, MACHEPL, MAXNUML;
extern double sqrt( double );
extern double atan( double );
extern double incbet( double, double, double );
extern double incbi( double, double, double );
extern double fabs( double );

#ifdef __LCC__
double __declspec(naked) students_t(int k,double t)
{
}
#endif
double stdtr(int k,double t )
{
	double x, rk, z, f, tz, p, xsqk;
	int j;

	if( k <= 0 )
	{
		mtherr( "stdtrl", EDOM );
		return(0.0);
	}

	if( t == 0.0 )
		return( 0.5 );

	if( t < -1.6 )
	{
		rk = k;
		z = rk / (rk + t * t);
		p = 0.5 * incbet( 0.5*rk, 0.5, z );
		return( p );
	}

	/*	compute integral from -t to + t */

	if( t < 0.0 )
		x = -t;
	else
		x = t;

	rk = k;	/* degrees of freedom */
	z = 1.0 + ( x * x )/rk;

	/* test if k is odd or even */
	if( (k & 1) != 0)
	{

		/*	computation for odd k	*/

		xsqk = x/sqrt(rk);
		p = atan( xsqk );
		if( k > 1 )
		{
			f = 1.0;
			tz = 1.0;
			j = 3;
			while(  (j<=(k-2)) && ( (tz/f) > MACHEPL )  )
			{
				tz *= (j-1)/( z * j );
				f += tz;
				j += 2;
			}
			p += f * xsqk/z;
		}
		p *= 2.0/PI;
	}


	else
	{

		/*	computation for even k	*/

		f = 1.0;
		tz = 1.0;
		j = 2;

		while(  ( j <= (k-2) ) && ( (tz/f) > MACHEPL )  )
		{
			tz *= (j - 1)/( z * j );
			f += tz;
			j += 2;
		}
		p = f * x/sqrt(z*rk);
	}

	/*	common exit	*/


	if( t < 0.0 )
		p = -p;	/* note destruction of relative accuracy */

	p = 0.5 + 0.5 * p;
	return(p);
}

#ifdef __LCC__
double __declspec(naked) students_t_inv(int k,double p)
{
}
#endif
double stdtri(int k,double p )
{
	double t, rk, z;
	int rflg;

	if( k <= 0 || p <= 0.0 || p >= 1.0 )
	{
		mtherr( "stdtril", EDOM );
		return(0.0);
	}

	rk = k;

	if( p > 0.25 && p < 0.75 )
	{
		if( p == 0.5 )
			return( 0.0 );
		z = 1.0 - 2.0 * p;
		z = incbi( 0.5, 0.5*rk, fabs(z) );
		t = sqrt( rk*z/(1.0-z) );
		if( p < 0.5 )
			t = -t;
		return( t );
	}
	rflg = -1;
	if( p >= 0.5)
	{
		p = 1.0 - p;
		rflg = 1;
	}
	z = incbi( 0.5*rk, 0.5, 2.0*p );

	if( MAXNUML * z < rk )
		return(rflg* MAXNUML);
	t = sqrt( rk/z - rk );
	return( rflg * t );
}
