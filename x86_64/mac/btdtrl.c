
/*							btdtrl.c
 *
 *	Beta distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * long double a, b, x, y, btdtrl();
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
extern long double incbetl ( long double, long double, long double );

#ifdef __LCC__
long double __declspec(naked) beta_distribution(long double a,long double b,long double x)
{
}
#endif
long double btdtrl( long double a,long double b,long double x )
{

	return( incbetl( a, b, x ) );
}

#ifndef __LCC__
long double beta_distribution(long double a,long double b,long double x)
{
	return btdtrl(a,b,x);
}
#endif
