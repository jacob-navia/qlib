/*							igamil()
*
 *      Inverse of complemented imcomplete gamma integral
*
 *
 *
 * SYNOPSIS:
*
 * double a, x, y, igamil();
*
 * x = igamil( a, y );
*
 *
 *
 * DESCRIPTION:
*
 * Given y, the function finds x such that
*
 *  igamc( a, x ) = y.
*
 * Starting with the approximate value
*
 *         3
*  x = a t
*
 *  where
*
 *  t = 1 - d - ndtri(y) sqrt(d)
*
* and
*
 *  d = 1/9a,
*
 * the routine performs up to 10 Newton iterations to find the
* root of igamc(a,x) - y = 0.
*
 *
 * ACCURACY:
*
 * Tested for a ranging from 0.5 to 30 and x from 0 to 0.5.
*
 *                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    DEC       0,0.5         3400       8.8e-16     1.3e-16
*    IEEE      0,0.5        10000       1.1e-14     1.0e-15
*
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/
#include <errno.h>
#include <float.h>
int mtherr (char *, int);
extern double normal_distribution_inv( double );
extern double exp( double );
extern double fabs( double );
extern double log( double );
extern double sqrt( double );
extern double lgamma( double );
extern double igamc( double, double );
#define MACHEPL DBL_EPSILON
#define MAXNUML DBL_MAX
#define MAXLOGL 1.1356523406294143949492E4
#define MINLOGL -1.1355137111933024058873E4

#ifdef __LCC__
double __declspec(naked) gamma_incomplete_cinv(double a,double y0)
{
}
#endif
double igami(double a,double y0 )
{
	double x0, x1, x, yl, yh, y, d, lgm, dithresh;
	int i, dir;

	/* bound the solution */
	x0 = DBL_MAX;
	yl = 0.0;
	x1 = 0.0;
	yh = 1.0;
	dithresh = 4.0 * MACHEPL;

	/* approximation to inverse function */
	d = 1.0L/(9.0L*a);
	y = ( 1.0L - d - normal_distribution_inv(y0) * sqrt(d) );
	x = a * y * y * y;

	lgm = lgamma(a);

	for( i=0; i<10; i++ )
	{
		if( x > x0 || x < x1 )
			goto ihalve;
		y = igamc(a,x);
		if( y < yl || y > yh )
			goto ihalve;
		if( y < y0 )
		{
			x0 = x;
			yl = y;
		}
		else
		{
			x1 = x;
			yh = y;
		}
		/* compute the derivative of the function at this point */
		d = (a - 1.0L) * log(x0) - x0 - lgm;
		if( d < -MAXLOGL )
			goto ihalve;
		d = -exp(d);
		/* compute the step to the next approximation of x */
		d = (y - y0)/d;
		x = x - d;
		if( i < 3 )
			continue;
		if( fabs(d/x) < dithresh )
			goto done;
	}

	/* Resort to interval halving if Newton iteration did not converge. */
ihalve:

	d = 0.0625;
	if( x0 == MAXNUML )
	{
		if( x <= 0.0L )
			x = 1.0L;
		while( x0 == MAXNUML )
		{
			x = (1.0 + d) * x;
			y = igamc( a, x );
			if( y < y0 )
			{
				x0 = x;
				yl = y;
				break;
			}
			d = d + d;
		}
	}
	d = 0.5;
	dir = 0;

	for( i=0; i<400; i++ )
	{
		x = x1  +  d * (x0 - x1);
		y = igamc( a, x );
		lgm = (x0 - x1)/(x1 + x0);
		if( fabs(lgm) < dithresh )
			break;
		lgm = (y - y0)/y0;
		if( fabs(lgm) < dithresh )
			break;
		if( x <= 0.0 )
			break;
		if( y > y0 )
		{
			x1 = x;
			yh = y;
			if( dir < 0 )
			{
				dir = 0;
				d = 0.5;
			}
			else if( dir > 1 )
				d = 0.5 * d + 0.5;
			else
				d = (y0 - yl)/(yh - yl);
			dir += 1;
		}
		else
		{
			x0 = x;
			yl = y;
			if( dir > 0 )
			{
				dir = 0;
				d = 0.5;
			}
			else if( dir < -1 )
				d = 0.5L * d;
			else
				d = (y0 - yl)/(yh - yl);
			dir -= 1;
		}
	}
	//if( x == 0.0L )
		//mtherr( "igamil", EUNDFL );

done:
	return( x );
}
