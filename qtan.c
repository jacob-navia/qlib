/*							qftan.c
 *      Circular tangent check routine
 *
 *
 *
 * SYNOPSIS:
 *
 * int qftan( x, y );
 * QELT *x, *y;
 *
 * qftan( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Domain of approximation is reduced by the transformation
 *
 * x -> x - pi floor((x + pi/2)/pi)
 *
 *
 * then tan(x) is the continued fraction
 *
 *                  2   2   2
 *             x   x   x   x
 * tan(x)  =  --- --- --- --- ...
 *            1 - 3 - 5 - 7 -
 *
 */
/*							qcot
 *
 *	Circular cotangent check routine
 *
 *
 *
 * SYNOPSIS:
 *
 * int qcot( x, y );
 * QELT *x, *y;
 *
 * qcot( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * cot (x) = 1 / tan (x).
 *
 */

#include "qhead.h"

void qftan(Qfloatp x,Qfloatp y)
{
	Qfloat xxx[1], e[1], r[1], j[1], xx[1], m2[1];
	int i, sign, n,z;
	int li;

	sign = signof(x);
	qmov( x, xxx );
	setpositive(xxx);
	/*	range reduction to +-pi/2		*/
	qadd( xxx, qPi_Div_2, e );	/*  x - pi * int( (x + pi/2)/pi  )	*/
	qmul( qinv_pi, e, e ); // Multiply by inverse instead of dividing by pi
	qfloor( e, e );
	qmul( e, qpi, e );
	qsub( e, xxx, xxx );

	qmov( qtwo, m2 );
	qneg( m2 );

	/* loop count = 17 for convergence to 9*16 bit mantissa if x < 1 */
	/* accuracy better than double precision for x beyond 2 */
	n = NBITS/9; //8;	/* 17 */
	li = 2*n + 1;
	itoq( li, j );
	qmov( j, e );
	qsquare( xxx, xx );
	qneg( xx );

	/* continued fraction expansion */
	for( i=0; i<n; i++) {
		qdiv( e, xx, r );
		qadd( m2, j, j );
		z = qadd( r, j, e );
		if (z == 0) break;
		// printf("[%d],exponent: %d 0x%x\n",
			// i,r[0].exponent-EXPONE,r[0].exponent-EXPONE);
		//printf("[%d],%140.132qf\n",i,r[0]);
	}

	qdiv( e, xxx, y );
	if( sign != 0 )
		setnegative(y);
}



void qcot(Qfloatp x,Qfloatp y)
{
	Qfloat z[1];

	qftan( x, z );
	qinv( z, y );
}
