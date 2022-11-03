/*							log.c
 *
 *	Natural logarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, log();
 *
 * y = log( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of x.
 *
 * The argument is separated into its exponent and fractional
 * parts.  If the exponent is between -1 and +1, the logarithm
 * of the fraction is approximated by
 *
 *     log(1+x) = x - 0.5 x**2 + x**3 P(x)/Q(x).
 *
 * Otherwise, setting  z = 2(x-1)/x+1),
 * 
 *     log(x) = z + z**3 P(z)/Q(z).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0.5, 2.0    150000      1.44e-16    5.06e-17
 *    IEEE      +-MAXNUM    30000       1.20e-16    4.78e-17
 *    DEC       0.5, 2.0    20000       2.0e-17     7.0e-18
 *    DEC       1, MAXNUM   17700       1.9e-17     6.1e-18
 *
 * In the tests over the interval [+-MAXNUM], the logarithms
 * of the random arguments were uniformly distributed over
 * [0, MAXLOG].
 *
 * ERROR MESSAGES:
 *
 * log singularity:  x = 0; returns MINLOG
 * log domain:       x < 0; returns MINLOG
 */

/*
Cephes Math Library Release 2.2:  December, 1990
Copyright 1984, 1990 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include "mconf.h"
static char fname[] = {"log"};

/* Coefficients for log(1+x) = x - x**2/2 + x**3 P(x)/Q(x)
 * 1/sqrt(2) <= x < sqrt(2)
 */
#ifdef UNK
static double P[] = {
 1.01875663804580931796E-4,
 4.97494994976747001425E-1,
 4.70579119878881725854E0,
 1.44989225341610930846E1,
 1.79368678507819816313E1,
 7.70838733755885391666E0,
};
static double Q[] = {
/* 1.00000000000000000000E0, */
 1.12873587189167450590E1,
 4.52279145837532221105E1,
 8.29875266912776603211E1,
 7.11544750618563894466E1,
 2.31251620126765340583E1,
};
#endif

/* Coefficients for log(x) = z + z**3 P(z)/Q(z),
 * where z = 2(x-1)/(x+1)
 * 1/sqrt(2) <= x < sqrt(2)
 */

static double R[3] = {
-7.89580278884799154124E-1,
 1.63866645699558079767E1,
-6.41409952958715622951E1,
};
static double S[3] = {
/* 1.00000000000000000000E0,*/
-3.56722798256324312549E1,
 3.12093766372244180303E2,
-7.69691943550460008604E2,
};

#define SQRTH 0.70710678118654752440
extern double MINLOG;

double log(x)
double x;
{
int e;
#ifdef DEC
short *q;
#endif
double y, z;
double frexp(), ldexp(), polevl(), p1evl();

/* Test for domain */
if( x <= 0.0 )
	{
	if( x == 0.0 )
		mtherr( fname, SING );
	else
		mtherr( fname, DOMAIN );
	return( MINLOG );
	}

/* separate mantissa from exponent */

#ifdef DEC
q = (short *)&x;
e = *q;			/* short containing exponent */
e = ((e >> 7) & 0377) - 0200;	/* the exponent */
*q &= 0177;	/* strip exponent from x */
*q |= 040000;	/* x now between 0.5 and 1 */
#endif

/* Note, frexp is used so that denormal numbers
 * will be handled properly.
 */
#ifdef IBMPC
x = frexp( x, &e );
/*
q = (short *)&x;
q += 3;
e = *q;
e = ((e >> 4) & 0x0fff) - 0x3fe;
*q &= 0x0f;
*q |= 0x3fe0;
*/
#endif

/* Equivalent C language standard library function: */
#ifdef UNK
x = frexp( x, &e );
#endif

#ifdef MIEEE
x = frexp( x, &e );
#endif



/* logarithm using log(x) = z + z**3 P(z)/Q(z),
 * where z = 2(x-1)/x+1)
 */

if( (e > 2) || (e < -2) )
{
if( x < SQRTH )
	{ /* 2( 2x-1 )/( 2x+1 ) */
	e -= 1;
	z = x - 0.5;
	y = 0.5 * z + 0.5;
	}	
else
	{ /*  2 (x-1)/(x+1)   */
	z = x - 0.5;
	z -= 0.5;
	y = 0.5 * x  + 0.5;
	}

x = z / y;


/* rational form */
z = x*x;
z = x * ( z * polevl( z, R, 2 ) / p1evl( z, S, 3 ) );
y = e;
z = z - y * 2.121944400546905827679e-4;
z = z + x;
z = z + e * 0.693359375;
goto ldone;
}



/* logarithm using log(1+x) = x - .5x**2 + x**3 P(x)/Q(x) */

if( x < SQRTH )
	{
	e -= 1;
	x = ldexp( x, 1 ) - 1.0; /*  2x - 1  */
	}	
else
	{
	x = x - 1.0;
	}


/* rational form */
z = x*x;
#if DEC
y = x * ( z * polevl( x, P, 5 ) / p1evl( x, Q, 6 ) );
#else
y = x * ( z * polevl( x, P, 5 ) / p1evl( x, Q, 5 ) );
#endif
if( e )
	y = y - e * 2.121944400546905827679e-4;
y = y - ldexp( z, -1 );   /*  y - 0.5 * z  */
z = x + y;
if( e )
	z = z + e * 0.693359375;

ldone:

return( z );
}
