/*							qerf.c
*
 *	Error function
*
 *
 *
 * SYNOPSIS:
*
 * int qerf( x, y );
* QELT *x, *y;
*
 * qerf( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * The integral is
*
 *                           x
*                            -
*                 2         | |          2
*   erf(x)  =  --------     |    exp( - t  ) dt.
*              sqrt(pi)   | |
*                          -
*                           0
*
 *
 * Cephes Math Library Release 2.2:  June, 1992
 * Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
 * Adapted by jacob navia to lcc-win.
0.25| 0.276326390168236932985068267764815712065353977892311254082471931262685083897302225185184805849741366790049258312413967593202915724630
0.5 | 0.520499877813046537682746653891964528736451575757963700058805725647193521716853570914788218734787757032966124386194391236065414690591
0.75| 0.711155633653515131598937834591410777374205954096537232278133397125036368764049561110793253890098900043530345370418642747913490195346
1.0 | 0.842700792949714869341220635082609259296066997966302908459937897834717254096010841261983325348144888454158261532021694364852339058255

*/

#include "qhead.h"

void EXPORT qerf(Qfloatp x,Qfloatp y )
{
	Qfloat xx[1], z[1], a[1], b[1];
	int i, j;
	int sign;

	qmov( x, xx );
	qmul( xx, xx, z );	/* square of y				*/

	sign = signof(xx);

	if( xx[0].exponent < (QELT) (EXPONE + 2) )
	{
		i = 57; /* 27;*/
		// This was losing precision in a BIG way.
		// Fixed. Jacob
		//if( xx[1] < EXPONE )
			//i = 20; /*10;*/
		j = 4 * i + 1;
		itoq( j, a);		/*  2 * i  +  1			*/
		qmov( a, y );

		/* continued fraction expansion */
		while( j > 1 )
		{
			qsub( qone, a, a );
			qmul( a, z, b );
			qdiv( y, b, y );
			qsub( qone, a, a );
			qadd( a, y, y );

			qsub( qone, a, a );
			qmul( a, z, b );
			qdiv( y, b, y );
			qsub( qone, a, a );
			qsub( y, a, y );
			j -= 4;
		}
		qdiv( y, xx, y );

		qneg(z);		/* exp( -xx**2 )	*/
		qfexp( z, a );
		qmul( a, y, y );
		qmul( oneopi, y, y );
	}

	/*							qerf.c 2	*/

	/* alternate continued fraction expansion for large x */
	else {
		j = 176;
		//if( xx[1] > (QELT) (EXPONE + 1) )
			//j = 120; /*76;*/
		itoq( j, a );
		if( sign != 0 )
			xx[0].sign = 0;
		qmov( xx, y );
		while( j > 0 ) {
			qdiv( y, a, y );
			y[0].exponent -= 1;
			qadd( xx, y, y );
			qsub( qone, a, a );
			j -= 1;
		}
		qdiv( y, qone, y );

		qneg(z);		/* exp( -xx**2 )	*/
		qfexp( z, a );
		qmul( a, y, y );
		qmul( oneopi, y, y );
		decreaseExponent(y);
		qsub( y, qone, y );
		if( sign != 0 )
			qneg(y);
	}
}
