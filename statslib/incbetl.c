/*							incbetl.c
*
 *	Incomplete beta integral
*
 *
 * SYNOPSIS:
*
 * double a, b, x, y, incbetl();
*
 * y = incbetl( a, b, x );
*
 *
 * DESCRIPTION:
*
 * Returns incomplete beta integral of the arguments, evaluated
* from zero to x.  The function is defined as
*
 *                  x
*     -            -
*    | (a+b)      | |  a-1     b-1
*  -----------    |   t   (1-t)   dt.
*   -     -     | |
*  | (a) | (b)   -
*                 0
*
 * The domain of definition is 0 <= x <= 1.  In this
* implementation a and b are restricted to positive values.
* The integral from x to 1 may be obtained by the symmetry
* relation
*
 *    1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).
*
 * The integral is evaluated by a continued fraction expansion
* or, when b*x is small, by a power series.
*
 * ACCURACY:
*
 * Tested at random points (a,b,x) with x between 0 and 1.
* arithmetic   domain     # trials      peak         rms
*    IEEE       0,5       20000        4.5e-18     2.4e-19
*    IEEE       0,100    100000        3.9e-17     1.0e-17
* Half-integer a, b:
*    IEEE      .5,10000  100000        3.9e-14     4.4e-15
* Outputs smaller than the IEEE gradual underflow threshold
* were excluded from these statistics.
*
 * ERROR MESSAGES:
*
 *   message         condition      value returned
* incbetl domain     x<0, x>1          0.0
*/
/*
Cephes Math Library, Release 2.3:  January, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

#include <errno.h>
#include <float.h>
int mtherr (char *, int);
#define MAXGAML 1755.455L
static double big = 9.223372036854775808e18L;
static double biginv = 1.084202172485504434007e-19L;
#define MACHEPL DBL_EPSILON
#define MINLOGL -1.1355137111933024058873E4L
#define MAXLOGL 1.1356523406294143949492E4L
extern double tgamma( double );
extern double lgamma( double );
extern double exp( double );
extern double log( double );
extern double fabs( double );
extern double pow( double, double );
static double incbcfl( double, double, double );
static double incbdl( double, double, double );
static double pseriesl( double, double, double );

#ifdef __LCC__
double __declspec(naked) beta_incomplete(double aa, double bb,double x)
{
}
#endif
double incbetl(double aa,double bb,double xx )
{
	double a, b, t, x, w, xc, y;
	int flag;

	if( aa <= 0.0L || bb <= 0.0L )
		goto domerr;

	if( (xx <= 0.0L) || ( xx >= 1.0L) )
	{
		if( xx == 0.0L )
			return( 0.0L );
		if( xx == 1.0L )
			return( 1.0L );
domerr:
		mtherr( "incbetl", EDOM );
		return( 0.0L );
	}

	flag = 0;
	if( (bb * xx) <= 1.0L && xx <= 0.95L)
	{
		t = pseriesl(aa, bb, xx);
		goto done;
	}

	w = 1.0L - xx;

	/* Reverse a and b if x is greater than the mean. */
	if( xx > (aa/(aa+bb)) )
	{
		flag = 1;
		a = bb;
		b = aa;
		xc = xx;
		x = w;
	}
	else
		{
		a = aa;
		b = bb;
		xc = w;
		x = xx;
	}

	if( flag == 1 && (b * x) <= 1.0 && x <= 0.95)
	{
		t = pseriesl(a, b, x);
		goto done;
	}

	/* Choose expansion for optimal convergence */
	y = x * (a+b-2.0) - (a-1.0);
	if( y < 0.0L )
		w = incbcfl( a, b, x );
	else
		w = incbdl( a, b, x ) / xc;

	/* Multiply w by the factor
	a      b   _             _     _
	x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */

	y = a * log(x);
	t = b * log(xc);
	if( (a+b) < MAXGAML && fabs(y) < MAXLOGL && fabs(t) < MAXLOGL )
	{
		t = pow(xc,b);
		t *= pow(x,a);
		t /= a;
		t *= w;
		t *= tgamma(a+b) / (tgamma(a) * tgamma(b));
		goto done;
	}
	else
	{
		/* Resort to logarithms.  */
		y += t + lgamma(a+b) - lgamma(a) - lgamma(b);
		y += log(w/a);
		if( y < MINLOGL )
			t = 0.0;
		else
			t = exp(y);
	}

done:

	if( flag == 1 )
	{
		if( t <= MACHEPL )
			t = 1.0 - MACHEPL;
		else
			t = 1.0 - t;
	}
	return( t );
}

/* Continued fraction expansion #1
* for incomplete beta integral
*/

static double incbcfl(double a,double b,double x )
{
	double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
	double k1, k2, k3, k4, k5, k6, k7, k8;
	double r, t, ans, thresh;
	int n;

	k1 = a;
	k2 = a + b;
	k3 = a;
	k4 = a + 1.0L;
	k5 = 1.0L;
	k6 = b - 1.0L;
	k7 = k4;
	k8 = a + 2.0L;

	pkm2 = 0.0L;
	qkm2 = 1.0L;
	pkm1 = 1.0L;
	qkm1 = 1.0L;
	ans = 1.0L;
	r = 1.0L;
	n = 0;
	thresh = 3.0L * MACHEPL;
	do
		{

		xk = -( x * k1 * k2 )/( k3 * k4 );
		pk = pkm1 +  pkm2 * xk;
		qk = qkm1 +  qkm2 * xk;
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;

		xk = ( x * k5 * k6 )/( k7 * k8 );
		pk = pkm1 +  pkm2 * xk;
		qk = qkm1 +  qkm2 * xk;
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;

		if( qk != 0.0L )
			r = pk/qk;
		if( r != 0.0L )
		{
			t = fabs( (ans - r)/r );
			ans = r;
		}
		else
			t = 1.0L;

		if( t < thresh )
			goto cdone;

		k1 += 1.0L;
		k2 += 1.0L;
		k3 += 2.0L;
		k4 += 2.0L;
		k5 += 1.0L;
		k6 -= 1.0L;
		k7 += 2.0L;
		k8 += 2.0L;

		if( (fabs(qk) + fabs(pk)) > big )
		{
			pkm2 *= biginv;
			pkm1 *= biginv;
			qkm2 *= biginv;
			qkm1 *= biginv;
		}
		if( (fabs(qk) < biginv) || (fabs(pk) < biginv) )
		{
			pkm2 *= big;
			pkm1 *= big;
			qkm2 *= big;
			qkm1 *= big;
		}
	}
	while( ++n < 400 );
	//mtherr( "incbetl", EPLOSS );

cdone:
	return(ans);
}


/* Continued fraction expansion #2
* for incomplete beta integral
*/

static double incbdl(double a,double b,double x )
{
	double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
	double k1, k2, k3, k4, k5, k6, k7, k8;
	double r, t, ans, z, thresh;
	int n;

	k1 = a;
	k2 = b - 1.0L;
	k3 = a;
	k4 = a + 1.0L;
	k5 = 1.0L;
	k6 = a + b;
	k7 = a + 1.0L;
	k8 = a + 2.0L;

	pkm2 = 0.0L;
	qkm2 = 1.0L;
	pkm1 = 1.0L;
	qkm1 = 1.0L;
	z = x / (1.0L-x);
	ans = 1.0L;
	r = 1.0L;
	n = 0;
	thresh = 3.0L * MACHEPL;
	do
		{

		xk = -( z * k1 * k2 )/( k3 * k4 );
		pk = pkm1 +  pkm2 * xk;
		qk = qkm1 +  qkm2 * xk;
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;

		xk = ( z * k5 * k6 )/( k7 * k8 );
		pk = pkm1 +  pkm2 * xk;
		qk = qkm1 +  qkm2 * xk;
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;

		if( qk != 0.0L )
			r = pk/qk;
		if( r != 0.0L )
		{
			t = fabs( (ans - r)/r );
			ans = r;
		}
		else
			t = 1.0L;

		if( t < thresh )
			goto cdone;

		k1 += 1.0L;
		k2 -= 1.0L;
		k3 += 2.0L;
		k4 += 2.0L;
		k5 += 1.0L;
		k6 += 1.0L;
		k7 += 2.0L;
		k8 += 2.0L;

		if( (fabs(qk) + fabs(pk)) > big )
		{
			pkm2 *= biginv;
			pkm1 *= biginv;
			qkm2 *= biginv;
			qkm1 *= biginv;
		}
		if( (fabs(qk) < biginv) || (fabs(pk) < biginv) )
		{
			pkm2 *= big;
			pkm1 *= big;
			qkm2 *= big;
			qkm1 *= big;
		}
	}
	while( ++n < 400 );
	//mtherr( "incbetl", EPLOSS );

cdone:
	return(ans);
}

/* Power series for incomplete gamma integral.
Use when b*x is small.  */

static double pseriesl(double a,double b,double x )
{
	double s, t, u, v, n, t1, z, ai;

	ai = 1.0L / a;
	u = (1.0L - b) * x;
	v = u / (a + 1.0L);
	t1 = v;
	t = u;
	n = 2.0L;
	s = 0.0L;
	z = MACHEPL * ai;
	while( fabs(v) > z )
	{
		u = (n - b) * x / n;
		t *= u;
		v = t / (a + n);
		s += v;
		n += 1.0L;
	}
	s += t1;
	s += ai;

	u = a * log(x);
	if( (a+b) < MAXGAML && fabs(u) < MAXLOGL )
	{
		t = tgamma(a+b)/(tgamma(a)*tgamma(b));
		s = s * t * pow(x,a);
	}
	else
		{
		t = lgamma(a+b) - lgamma(a) - lgamma(b) + u + log(s);
		if( t < MINLOGL )
			s = 0.0L;
		else
			s = exp(t);
	}
	return(s);
}
