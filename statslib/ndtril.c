/*							ndtril.c
*
 *	Inverse of Normal distribution function
*
 *
 *
 * SYNOPSIS:
*
 * long double x, y, ndtril();
*
 * x = ndtril( y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the argument, x, for which the area under the
* Gaussian probability density function (integrated from
* minus infinity to x) is equal to y.
*
 *
 * For small arguments 0 < y < exp(-2), the program computes
* z = sqrt( -2 log(y) );  then the approximation is
* x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z) .
* For larger arguments,  x/sqrt(2 pi) = w + w^3 R(w^2)/S(w^2)) ,
* where w = y - 0.5 .
*
 * ACCURACY:
*
 *                      Relative error:
* arithmetic   domain        # trials      peak         rms
*  Arguments uniformly distributed:
*    IEEE       0, 1           5000       7.8e-19     9.9e-20
*  Arguments exponentially distributed:
*    IEEE     exp(-11355),-1  30000       1.7e-19     4.3e-20
*
 *
 * ERROR MESSAGES:
*
 *   message         condition    value returned
* ndtril domain      x <= 0        -MAXNUML
* ndtril domain      x >= 1         MAXNUML
*
 */


/*
Cephes Math Library Release 2.3:  January, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

#include <errno.h>
int mtherr (char *, int);
#define MAXNUML 1.189731495357231765021263853E4932L

/* ndtri(y+0.5)/sqrt(2 pi) = y + y^3 R(y^2)
0 <= y <= 3/8
Peak relative error 6.8e-21.  */
static long double s2pi = 2.506628274631000502416E0L;
static long double P0[8] = {
 8.779679420055069160496E-3L,
-7.649544967784380691785E-1L,
 2.971493676711545292135E0L,
-4.144980036933753828858E0L,
 2.765359913000830285937E0L,
-9.570456817794268907847E-1L,
 1.659219375097958322098E-1L,
-1.140013969885358273307E-2L,
};
static long double Q0[7] = {
/* 1.000000000000000000000E0L, */
-5.303846964603721860329E0L,
 9.908875375256718220854E0L,
-9.031318655459381388888E0L,
 4.496118508523213950686E0L,
-1.250016921424819972516E0L,
 1.823840725000038842075E-1L,
-1.088633151006419263153E-2L,
};
/* Approximation for interval z = sqrt(-2 log y ) between 2 and 8
*/
/*  ndtri(p) = z - ln(z)/z - 1/z P1(1/z)/Q1(1/z)
z = sqrt(-2 ln(p))
2 <= z <= 8, i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
Peak relative error 5.3e-21  */
static long double P1[10] = {
 4.302849750435552180717E0L,
 4.360209451837096682600E1L,
 9.454613328844768318162E1L,
 9.336735653151873871756E1L,
 5.305046472191852391737E1L,
 1.775851836288460008093E1L,
 3.640308340137013109859E0L,
 3.691354900171224122390E-1L,
 1.403530274998072987187E-2L,
 1.377145111380960566197E-4L,
};
static long double Q1[9] = {
/* 1.000000000000000000000E0L, */
 2.001425109170530136741E1L,
 7.079893963891488254284E1L,
 8.033277265194672063478E1L,
 5.034715121553662712917E1L,
 1.779820137342627204153E1L,
 3.845554944954699547539E0L,
 3.993627390181238962857E-1L,
 1.526870689522191191380E-2L,
 1.498700676286675466900E-4L,
};
/* ndtri(x) = z - ln(z)/z - 1/z P2(1/z)/Q2(1/z)
z = sqrt(-2 ln(y))
8 <= z <= 32
i.e., y between exp(-32) = 1.27e-14 and exp(-512) = 4.38e-223
Peak relative error 1.0e-21  */
static long double P2[8] = {
 3.244525725312906932464E0L,
 6.856256488128415760904E0L,
 3.765479340423144482796E0L,
 1.240893301734538935324E0L,
 1.740282292791367834724E-1L,
 9.082834200993107441750E-3L,
 1.617870121822776093899E-4L,
 7.377405643054504178605E-7L,
};
static long double Q2[7] = {
/* 1.000000000000000000000E0L, */
 6.021509481727510630722E0L,
 3.528463857156936773982E0L,
 1.289185315656302878699E0L,
 1.874290142615703609510E-1L,
 9.867655920899636109122E-3L,
 1.760452434084258930442E-4L,
 8.028288500688538331773E-7L,
};
/*  ndtri(x) = z - ln(z)/z - 1/z P3(1/z)/Q3(1/z)
32 < z < 2048/13
Peak relative error 1.4e-20  */
static long double P3[8] = {
 2.020331091302772535752E0L,
 2.133020661587413053144E0L,
 2.114822217898707063183E-1L,
-6.500909615246067985872E-3L,
-7.279315200737344309241E-4L,
-1.275404675610280787619E-5L,
-6.433966387613344714022E-8L,
-7.772828380948163386917E-11L,
};
static long double Q3[7] = {
/* 1.000000000000000000000E0L, */
 2.278210997153449199574E0L,
 2.345321838870438196534E-1L,
-6.916708899719964982855E-3L,
-7.908542088737858288849E-4L,
-1.387652389480217178984E-5L,
-7.001476867559193780666E-8L,
-8.458494263787680376729E-11L,
};
extern long double polevll ( long double, void *, int );
extern long double p1evll ( long double, void *, int );
extern long double logl ( long double );
extern long double sqrtl ( long double );

long double normal_distribution_inv(long double y0)
{
	long double x, y, z, y2, x0, x1;
	int code;

	if( y0 <= 0.0L )
	{
		mtherr( "ndtril", EDOM );
		return( -MAXNUML );
	}
	if( y0 >= 1.0L )
	{
		mtherr( "ndtri", EDOM );
		return( MAXNUML );
	}
	code = 1;
	y = y0;
	// exp(-2) = 0.135335283236612691893999494972484403407631545909575881468158872654073374101487689937098122490657048755077
	if( y > (1.0L - 0.1353352832366126918939L) ) /* 0.135... = exp(-2) */
	{
		y = 1.0L - y;
		code = 0;
	}

	if( y > 0.1353352832366126918939L )
	{
		y = y - 0.5L;
		y2 = y * y;
		x = y + y * (y2 * polevll( y2, P0, 7 )/p1evll( y2, Q0, 7 ));
		x = x * s2pi;
		return(x);
	}

	x = sqrtl( -2.0L * logl(y) );
	x0 = x - logl(x)/x;
	z = 1.0L/x;
	if( x < 8.0L )
		x1 = z * polevll( z, P1, 9 )/p1evll( z, Q1, 9 );
	else if( x < 32.0L )
		x1 = z * polevll( z, P2, 7 )/p1evll( z, Q2, 7 );
	else
	    x1 = z * polevll( z, P3, 7 )/p1evll( z, Q3, 7 );
	x = x0 - x1;
	if( code != 0 )
		x = -x;
	return( x );
}
