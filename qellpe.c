/*							qellpe.c
 *
 *	Complete elliptic integral of the second kind
 *
 *
 *
 * SYNOPSIS:
 *
 * int qellpe(x, y);
 * QELT *x, *y;
 *
 * qellpe(x, y);
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *            pi/2
 *             -
 *            | |                 2
 * E(m)  =    |    sqrt( 1 - m sin t ) dt
 *          | |
 *           -
 *            0
 *
 * Where m = 1 - m1, using the arithmetic-geometric mean method.
 *
 *
 * ACCURACY:
 *
 * Original Method terminates at NBITS/2.
 * Now accuracy has been increased to NBITS. jn
 *
 * Cephes Math Library, Release 2.1:  February, 1989
 * Copyright 1984, 1987, 1989 by Stephen L. Moshier
 * Adapted to lcc-win by jacob navia
 * Evaluation points (functions.wolfram.com)
 * 0.5:
 * 1.35064388104767550252017473533872584134952236692435454532\
 *   32537088578778908361273690402360778224915636099470783313\
 *   4660666779461120524
 * 0.25:
 * 1.46746220933942715545979526699091613602536175232723196050\
 *   07906364908242272712906356540385307335046021821754957207\
 *   9406740695677032698
 *
 * Complete elliptic integral of second kind
 *	Argument is m
 */

#include "qhead.h"

void qellpe(Qfloatp x,Qfloatp y)
{
	Qfloat a[1], b[1], c[1], d[1], e[1], temp[1];

	if( qcmpto1( x ) > 0 || x[0].sign != 0 )
	{
		mtherr( "qellpe", DOMAIN );
		return;
	}
	if( x[0].exponent < (QELT) (EXPONE - 2*NBITS) )
	{
		qmov( qone, y );
		return;
	}
	qsub( x, qone, temp );
	if( temp[0].exponent < (QELT) (EXPONE - 129) )
	{
		qmov( qpi, y );
		y->exponent -= 1;
		return;
	}
	// This function was returning the results for 1-x instead 
	// of the correct ones, at least according to Mathematica and
	// MATLAB. I have replaced the code with arrows with exclamation
	// marks
	//                      qmov( temp, e ); <------
	qmov(x,e); // !!!
	qfsqrt( temp, c );		/* c = sqrt(x)		*/
	/*qsub( temp, qone, b );*/	/* b = sqrt( 1 - x )	*/
	//                      qfsqrt( x, b ); <-------
	qfsqrt(temp,b); // !!!
	qmov( qone, a );	/* a = 1	*/
	qmov( qone, d );

	while( ((int) a[0].exponent - (int) c[0].exponent) < (NBITS/2) )
	{
		qsub( b, a, c );	/* c = (a - b)/2.0	*/
		if (c[0].exponent == 0)
			break;
		c[0].exponent -= 1;
		qmul( a, b, temp );	/* temp = ( a * b )	*/
		qadd( a, b, a );	/* a = (a + b)/2.0	*/
		a[0].exponent -= 1;
		qfsqrt( temp, b ); // b = sqrt(a * b)
		//qadd( d, d, d );	/* d += d 		*/
		d->exponent += 1;
		qsquare( c, temp );	/* e += c * c * d	*/
		qmul( temp, d, temp );
		qadd( e, temp, e );
	}

//	qmov( qpi, temp );	/* get pi/2	*/
//	temp[0].exponent -= 1;
	qdiv( a, qPi_Div_2, temp );	/* Integral of first kind */
	qsub( e, qtwo, b );	/* ( 2 - e )/2	*/
	b[0].exponent -= 1;
	qmul( temp, b, y );
}
