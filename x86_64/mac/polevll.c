/*							polevll.c
 *							p1evll.c
 *
 *	Evaluate polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * long double x, y, coef[N+1], polevl[];
 *
 * y = polevll( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 *  The function p1evll() assumes that coef[N] = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as polevll().
 *
 *  This module also contains the following globally declared constants:
 * MAXNUML = 1.189731495357231765021263853E4932L;
 * MACHEPL = 5.42101086242752217003726400434970855712890625E-20L;
 * MAXLOGL =  1.1356523406294143949492E4L;
 * MINLOGL = -1.1355137111933024058873E4L;
 * LOGE2L  = 6.9314718055994530941723E-1L;
 * LOG2EL  = 1.4426950408889634073599E0L;
 * PIL     = 3.1415926535897932384626L;
 * PIO2L   = 1.5707963267948966192313L;
 * PIO4L   = 7.8539816339744830961566E-1L;
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic.  This routine is used by most of
 * the functions in the library.  Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *
 */


/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/
#include <fenv.h>
/* almost 2^16384 */
long double MAXNUML = 1.189731495357231765021263853E4932L;
/* 2^-64 */
long double MACHEPL = 5.42101086242752217003726400434970855712890625E-20L;
/* log( MAXNUML ) */
long double MAXLOGL =  1.1356523406294143949492E4L;
/* log(smallest denormal number = 2^-16446) */
long double MINLOGL = -1.13994985314888605586758E4L;
long double LOGE2L  = 6.9314718055994530941723E-1L;
long double LOG2EL  = 1.4426950408889634073599E0L;
long double PIL     = 3.1415926535897932384626L;
long double PIO2L   = 1.5707963267948966192313L;
long double PIO4L   = 7.8539816339744830961566E-1L;
long double NANL = 0.0L / 0.0L;
long double NEGZEROL = -0.0L;

/* Polynomial evaluator:
 *  P[0] x^n  +  P[1] x^(n-1)  +  ...  +  P[n]
 */
long double polevll(long double x,long double *p,int n )
{
	long double y;

	y = *p++;
	do {
		y = y * x + *p++;
	}
	while( --n );
	return y;
}



/* Polynomial evaluator:
 *  x^n  +  P[0] x^(n-1)  +  P[1] x^(n-2)  +  ...  +  P[n]
 */
long double p1evll(long double x,long double *p,int n )
{
	long double y;

	n -= 1;
	y = x + *p++;
	do	{
		y = y * x + *p++;
	}
	while( --n );
	return y;
}
