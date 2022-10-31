/*							qfresnl
*
 *	Fresnel integral
*
 *
 *
 * SYNOPSIS:
*
 * int qfresnl( x, s, c );
* QELT *x, *s, *c;
*
 * qfresnl( x, s, c );
*
 *
 * DESCRIPTION:
*
 * Evaluates the Fresnel integrals
*
 *           x
*           -
*          | |
* C(x) =   |   cos(pi/2 t**2) dt,
*        | |
*         -
*          0
*
 *           x
*           -
*          | |
* S(x) =   |   sin(pi/2 t**2) dt.
*        | |
*         -
*          0
*
 *
 * The integrals are evaluated by a power series for x < 1.
* For large x auxiliary functions f(x) and g(x) are employed
* such that
*
 * C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
* S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
*
 * Routine qfresfg computes f and g.
*
 *
 * ACCURACY:
*
 * Series expansions are truncated at less than full working precision.
*/

/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
*/

#include <stdio.h>
#include "qhead.h"

/* asymptotic expansion converges to 144 bits at BRKPT >= 8.2 */
#if NBITS > 144
#define BRKPT 8.2
#else
#define BRKPT 5.7
#endif
static double dtemp = 0.0;
static Qfloat x2[1];
static Qfloat c[1];
static Qfloat q[1];
static Qfloat d[1];
static Qfloat t[1];
static Qfloat x4[1];
static Qfloat x3[1];
static Qfloat t1[1];
static Qfloat cc0[1];
static Qfloat ss0[1];
static Qfloat f[1];
static Qfloat g[1];

int qfresfg( Qfloatp, Qfloatp, Qfloatp );

int qfresnl(Qfloatp const x,Qfloatp const ss,Qfloatp cc)
{
	int max;
	int sign, i, iconv, ic2;
	int il;

	dtemp = qtoe( x, DOROUNDING );
	if( dtemp > BRKPT ) /* 7.0 5.7 */
		goto asymp0;

	max = 0;
	qmul( x, x, x2 );
	qmul( x, x2, x3 );	/* x3 = x**3		*/
	qmul( x2, x2, x4 );	/* x4 = x**4		*/
	sign = -1;
	qmov( qone, c );	/* c=1			*/
	qmov( qpi, q );
	q[0].exponent -= 1;
	qmul( q, x3, x3 );	/* x3 = pi/2 * x**3	*/
	qmul( q, q, q );	/* q = (pi/2)**2	*/

	qmov( qone, t );
	qmov( qone, cc0 );
	qmov( qone, ss0 );
	il = 3;
	itoq( il, d );
	qdiv( d, ss0, ss0 );	/* ss = pi/2 * x**3 / 3	*/
	i = 1;
	do
		{
		qadd( c, qone, c );	/* c += 1		*/
		qdiv( c, t, t );	/* t /= c		*/
		qmul( q, t, t );	/* t *= (pi/2)**2	*/
		qmul( t, x4, t );	/* t *= x**4		*/
		il = 4*i + 1;
		itoq( il, d );
		qmov( t, f );
		qdiv( d, f, f );
		if( f[0].exponent > max )
			max = f[0].exponent;
		if( sign < 0 )
			qsub( f, cc0, cc0 );
		else
			qadd( cc0, f, cc0 );
		iconv = cc0[0].exponent - f[0].exponent;
		qadd( c, qone, c );
		qdiv( c, t, t );
		qmov( t, f );
		il = 4*i + 3;
		itoq( il, d );
		qdiv( d, f, f );
		if( sign < 0 )
			qsub( f, ss0, ss0 );
		else
			qadd( ss0, f, ss0 );
		if( f[0].exponent > max )
			max = f[0].exponent;
		ic2 = (int) ss0[0].exponent - (int) f[0].exponent;
		if( ic2 < iconv )
			iconv = ic2;
		sign = -sign;
		++i;
	}
	while( iconv < 144 );

	max = NBITS - max;
	if( ss0[0].exponent < cc0[0].exponent )
		max += ss0[0].exponent;
	else
		max += cc0[0].exponent;
	if( max < (QELT) (3*NBITS/7) )
		printf( "qfresf: prec %d bits\n", max );

	qmul( x, cc0, cc );		/*	C(x)	*/
	qmul( x3, ss0, ss );		/*	S(x)	*/
	goto done;


asymp0:

	qfresfg( x, f, g );

	qmov( qpi, t );		/*	t = PIO2 * x * x;	*/
	t[0].exponent -= 1;
	qmul( x, t, t );
	qmul( x, t, t );
	qfcos( t, c );		/*	c = cos(t);		*/
	qfsin( t, t );		/*	t = sin(t);		*/
	qmov( qone, cc );	/*	cc = 0.5  +  f * t  -  g * c;	*/
	cc[0].exponent -= 1;
	qmul( f, t, x3 );
	qadd( x3, cc, cc );
	qmul( g, c, x3 );
	qsub( x3, cc, cc );

	qmov( qone, ss );	/*	ss = 0.5  -  f * c  -  g * t;	*/
	ss[0].exponent -= 1;
	qmul( f, c, x3 );
	qsub( x3, ss, ss );
	qmul( g, t, x3 );
	qsub( x3, ss, ss );

done:

	if( x[0].sign != 0 )
	{
		cc[0].sign = -1;
		ss[0].sign = -1;
	}

	return(0);
}

/*		the auxiliary functions		*/
static Qfloat lss[1] = { 0};
static Qfloat lcc[1] = { 0};

int qfresfg(Qfloatp x,Qfloatp f,Qfloatp g )
{
	int sign, max;
#define co c
#define si q

	dtemp = qtoe( x, DOROUNDING );
	if( dtemp <= BRKPT )
		goto nasymp;


	/*		Asymptotic power series auxiliary functions
	*		for large argument
	*/

	qmov( qone, c );	/* 	c = 1.0;		*/
	qmov( qone, t );	/*	t = 1.0;		*/
	qadd( qone, qone, t1 );	/*	t1 = 2.0;		*/
	qmul( qpi, x, q );	/*	q = PIO2 * 2.0 * x;	*/
	qmul( q, x, x3 );	/*	x3 = nn * x;		*/
	qmul( x3, x3, x4 );	/*	x4 = x3 * x3;		*/
	sign = -1;
	qmov( qone, f );	/*	f = 1.0;		*/
	qmov( qone, g );	/*	g = 1.0;		*/

	/*
	if( dtemp > 20000.0 )
	goto done5;
	*/

	while( ((int) qone[0].exponent - (int) t[0].exponent) < NBITS/2 )
	{
		qadd( qtwo, c, c );	/*	c += 2.0;		*/
		qmul( t, c, t );	/*	t *= c / x4;		*/
		qdiv( x4, t, t );
		if( sign < 0 )
			qsub( t, f, f ); /*	f -= t;			*/
		else
			qadd( t, f, f ); /*	f += t;			*/
		qadd( qtwo, c, c );	/*	c += 2.0;		*/
		qmul( t, c, t );	/*	t *= c;			*/
		/*
		qtoe( t, (unsigned short *) &dtemp );
		printf( "%.4E\n", dtemp );
		*/
		if( t[0].exponent > t1[0].exponent )
		{
			max = qone[0].exponent - t1[0].exponent;
			if( max < 3*NBITS/7 )
				printf( "qfres asymp prec %d bits\n", max );
			goto done5;
		}
		qmov( t, t1 );		/*	t1 = t;			*/
		if( sign < 0 )
			qsub( t, g, g ); /*	g -= t;			*/
		else
			qadd( t, g, g ); /*	g += t;			*/
		sign = -sign;
	}


done5:
	qmul( x3, q, d );	/*	g /= x3 * q;		*/
	qdiv( d, g, g );
	qdiv( q, f, f );	/*	f /= q;		*/
	return 0;





nasymp:

	qfresnl( x, lss, lcc );
	qmov( qone, co );
	co[0].exponent -= 1;
	qsub( lcc, co, lcc );
	qsub( lss, co, lss );

	qmul( qpi, x2, x2 );
	x2[0].exponent -= 1;
	qfcos( x2, co );
	qfsin( x2, si );

	/* f = (.5-ss)*cos - (.5-cc)*sin	*/

	qmul( lss, co, f );
	qmul( lcc, si, t );
	qsub( t, f, f );

	/* g = (.5-cc)*cos + (.5-ss)*sin	*/
	qmul( lcc, co, g );
	qmul( lss, si, t );
	qadd( t, g, g );
	return 0;
}

