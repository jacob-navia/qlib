
/*							btdtrl.c
 *
 *	Beta distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, btdtrl();
 *
 * y = btdtrl( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area from zero to x under the beta density
 * function:
 *
 *
 *                          x
 *            -             -
 *           | (a+b)       | |  a-1      b-1
 * P(x)  =  ----------     |   t    (1-t)    dt
 *           -     -     | |
 *          | (a) | (b)   -
 *                         0
 *
 *
 * The mean value of this distribution is a/(a+b).  The variance
 * is ab/[(a+b)^2 (a+b+1)].
 *
 * This function is identical to the incomplete beta integral
 * function, incbetl(a, b, x).
 *
 * The complemented function is
 *
 * 1 - P(1-x)  =  incbetl( b, a, x );
 *
 *
 * ACCURACY:
 *
 * See incbetl.c.
 *
 *
								btdtrl() */
/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1995 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/
extern double incbet( double, double, double );
#ifdef __LCC
double __declspec(naked) beta_distribution(double a,double b,double x)
{
}
#endif
double btdtr( double a,double b,double x )
{

	return( incbet( a, b, x ) );
}
