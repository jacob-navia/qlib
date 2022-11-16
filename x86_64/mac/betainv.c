/*							incbi()
*
 *      Inverse of imcomplete beta integral
*
 *
 *
 * SYNOPSIS:
*
 * long double a, b, x, y, beta_distribution_inv();
*
 * x = beta_distribution_inv( a, b, y );
*
 *
 *
 * DESCRIPTION:
*
 * Given y, the function finds x such that
*
 *  beta_distribution( a, b, x ) = y .
*
 * The routine performs interval halving or Newton iterations to find the
* root of incbet(a,b,x) - y = 0.
*
 *
 * ACCURACY:
*
 *                      Relative error:
*                x     a,b
* arithmetic   domain  domain  # trials    peak       rms
*    IEEE      0,1   .25,100     50000    1.4e-13   4.0e-15
*    IEEE      0,1    .5,10000    5000    4.0e-12   2.0e-13
*    IEEE      0,1     0,5       45000    7.0e-13   5.4e-15
* With a and b constrained to half-integer or integer values:
*    IEEE      0,1    .5,100     10000    2.8e-14   1.4e-15
* With a = .5, b constrained to half-integer or integer values:
*    IEEE      0,1    .5,10000    2000    1.1e-10   1.6e-11
*/


/*
Cephes Math Library Release 2.3:  March,1995
Copyright 1984, 1995 by Stephen L. Moshier
Adapted to lcc-win32 by Jacob Navia
*/

#include <math.h>
#include <stats.h>
#include <float.h>
#include "qhead.h"
#define MINLOG  -8.872283911167299960540E1
extern long double normal_distribution_inv(long double y0);
long double beta_distribution_inv(long double aa,long double bb,long double yy0 )
{
	long double a, b, y0, d, y, x, x0, x1, lgm, yp, di, dithresh, yl, yh;
	int i, rflg, dir, nflg;


	if( yy0 <= 0 )
		return(0.0);
	if( yy0 >= 1.0 )
		return(1.0);

	if( aa <= 1.0 || bb <= 1.0 )
	{
		nflg = 1;
		dithresh = 4.0 * LDBL_EPSILON;
		rflg = 0;
		a = aa;
		b = bb;
		y0 = yy0;
		x = a/(a+b);
		y = beta_distribution( a, b, x );
		goto ihalve;
	}
	else
	{
		nflg = 0;
		dithresh = 1.0e-4;
	}
	/* approximation to inverse function */

	yp = -normal_distribution_inv(yy0);

	if( yy0 > 0.5 )
	{
		rflg = 1;
		a = bb;
		b = aa;
		y0 = 1.0 - yy0;
		yp = -yp;
	}
	else
	{
		rflg = 0;
		a = aa;
		b = bb;
		y0 = yy0;
	}

	lgm = (yp * yp - 3.0)/6.0;
	x0 = 2.0/( 1.0/(2.0*a-1.0)  +  1.0/(2.0*b-1.0) );
	y = yp * sqrtl( x0 + lgm ) / x0
		- ( 1.0/(2.0*b-1.0) - 1.0/(2.0*a-1.0) )
		* (lgm + 5.0/6.0 - 2.0/(3.0*x0));
	y = 2.0 * y;
	if( y < MINLOG )
	{
		x0 = 1.0;
		goto under;
	}
	x = a/( a + b * exp(y) );
	y = beta_distribution( a, b, x );
	yp = (y - y0)/y0;
	if( fabsl(yp) < 1.0e-2 )
		goto newt;

ihalve:

	/* Resort to interval halving if not close enough */
	x0 = 0.0;
	yl = 0.0;
	x1 = 1.0;
	yh = 1.0;
	di = 0.5;
	dir = 0;

	for( i=0; i<400; i++ )
	{
		if( i != 0 )
		{
			x = x0  +  di * (x1 - x0);
			if( x == 1.0 )
				x = 1.0 - LDBL_EPSILON;
			y = beta_distribution( a, b, x );
			yp = (x1 - x0)/(x1 + x0);
			if( fabsl(yp) < dithresh )
			{
				x0 = x;
				goto newt;
			}
		}

		if( y < y0 )
		{
			x0 = x;
			yl = y;
			if( dir < 0 )
			{
				dir = 0;
				di = 0.5;
			}
			else if( dir > 1 )
				di = 0.5 * di + 0.5;
			else
				di = (y0 - y)/(yh - yl);
			dir += 1;
			if( x0 > 0.75 )
			{
				if( rflg == 1 )
				{
					rflg = 0;
					a = aa;
					b = bb;
					y0 = yy0;
				}
				else
				{
					rflg = 1;
					a = bb;
					b = aa;
					y0 = 1.0 - yy0;
				}
				x = 1.0 - x;
				y = beta_distribution( a, b, x );
				goto ihalve;
			}
		}
		else
		{
			x1 = x;
			if( rflg == 1 && x1 < LDBL_EPSILON )
			{
				x0 = 0.0;
				goto done;
			}
			yh = y;
			if( dir > 0 )
			{
				dir = 0;
				di = 0.5;
			}
			else if( dir < -1 )
				di = 0.5 * di;
			else
				di = (y - y0)/(yh - yl);
			dir -= 1;
		}
	}
	mtherr( "beta_distribution_inv", PLOSS );
	if( x0 >= 1.0 )
	{
		x0 = 1.0 - LDBL_EPSILON;
		goto done;
	}
	if( x == 0.0 )
	{
under:
		mtherr( "beta_distribution_inv", UNDERFLOW );
		x0 = 0.0;
		goto done;
	}

newt:

	if( nflg )
		goto done;

	x0 = x;
	lgm = lgamma(a+b) - lgamma(a) - lgamma(b);

	for( i=0; i<20; i++ )
	{
		/* Compute the function at this point. */
		if( i != 0 )
			y = beta_distribution(a,b,x0);
		/* Compute the derivative of the function at this point. */
		d = (a - 1.0) * log(x0) + (b - 1.0) * log(1.0-x0) + lgm;
		if( d < MINLOG )
		{
			x0 = 0.0;
			goto under;
		}
		d = expl(d);
		/* compute the step to the next approximation of x */
		d = (y - y0)/d;
		x = x0;
		x0 = x0 - d;
		if( x0 <= 0.0 )
		{
			x0 = 0.0;
			goto under;
		}
		if( x0 >= 1.0 )
		{
			x0 = 1.0 - LDBL_EPSILON;
			goto done;
		}
		if( fabsl(d/x0) < 64.0 * LDBL_EPSILON )
			goto done;
	}

done:
	if( rflg )
	{
		if( x0 <= LDBL_EPSILON )
			x0 = 1.0 - LDBL_EPSILON;
		else
			x0 = 1.0 - x0;
	}
	return( x0 );
}
