/*							qndtr.c
*
 *	Normal distribution function
*
 *
 *
 * SYNOPSIS:
*
 * int qndtr( x, y );
* QELT *x, *y;
*
 * qndtr( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the area under the Gaussian probability density
* function, integrated from minus infinity to x:
*
 *                            x
*                             -
*                   1        | |          2
*    ndtr(x)  = ---------    |    exp( - t /2 ) dt
*               sqrt(2pi)  | |
*                           -
*                          -inf.
*
 *             =  ( 1 + erf(z) ) / 2
*             =  erfc(z) / 2
*
 * where z = x/sqrt(2).
*
 */

#include "qhead.h"

/* loop counts adjusted for full convergence to 9 word mantissa */


void qndtr(Qfloatp x,Qfloatp y )
{
	int i,j, sign;
	Qfloat xx[1], z[1], a[1], b[1];

	qmul( qinv_sqrt2, x, xx );
	qsquare( xx, z );	/* square of y				*/

	sign = xx[0].sign;
	xx[0].sign = 0;

	if( xx[0].exponent < (QELT) (EXPONE + 2) ) /* xx < 4 */
	{

		i = NBITS * 3 / 8;  /*  168 */

		/*
		* if( xx[0].exponent < EXPONE )
		*	i = NBITS/6;
		*/

		j = 4 * i + 1;
		itoq( j, y);		/*  2 * i  +  1			*/
		//qmov( a, y );

		/* continued fraction expansion
		* Hart et al., p. 137; AMS55 #26.2.15
		*                    1                  x    x**2   2 x**2   3 x**2
		* P(x) - 0.5  =  --------- exp(-x**2)  ---  ------  ------   ------ ...
		*                sqrt(2pi)             1 -    3 +     5 -      7 +
		*/
		while( j > 1) {
			j--;
			itoq(j, a);
			qmuli( a, z, b );
			qdiv( y, b, y );
			j--;
			itoq(j, a);
			qadd( a, y, y );

			j--;
			itoq(j, a);
			qmuli( a, z, b );
			qdiv( y, b, y );
			j--;
			itoq(j, a);
			qsub( y, a, y );
		}
		qdiv( y, xx, y );
		/* *(y+1) += 1; */
		qneg(z);		/* exp( -xx**2 )	*/
		qfexp( z, a );
		qmul( a, y, y );
		qmul( oneopi, y, y );

		if( sign != 0 )
			qneg(y);

		qadd( y, qone, y );
		y[0].exponent -= 1;

	}


	/* alternate continued fraction expansion for large x
	* AMS55 #7.1.14
	*
 *              inf.
	*                -
	*               | |                1   1/2  2/2  3/2  4/2
	* 2 exp(x**2)   |   exp(-t**2) =  ---  ---  ---  ---  ---  ...
	*             | |                 x +  x +  x +  x +  x +
	*              -
	*               x
	*/
	else
	{

		j = NBITS * 6 / 5;  /* 172 */

		/*
		*if( xx[0].exponent > (EXPONE + 1) )
		*	j = NBITS;
		*/

		if( sign != 0 )
			xx[0].sign = 0;
		qmov( xx, y );
		while( j > 0 ) {
			itoq(j, a);
			qdiv( y, a, y );
			y->exponent -= 1;
			qadd( xx, y, y );
			j -= 1;
		}
		qinv( y, y );
		qneg(z);		/* exp( -xx**2 )	*/
		qfexp( z, a );
		qmul( a, y, y );
		qmul( oneopi, y, y );

		y->exponent -= 2;
		if( sign == 0 )
			qsub( y, qone, y );
	}

}
