/*							qerfc.c
 *
 *	Complementary error function
 *
 *
 *
 * SYNOPSIS:
 *
 * int qerfc( x, y );
 * QELT *x, *y;
 *
 * qerfc( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *  1 - erf(x) =
 *
 *                           inf.
 *                             -
 *                  2         | |          2
 *   erfc(x)  =  --------     |    exp( - t  ) dt
 *               sqrt(pi)   | |
 *                           -
 *                            x
 *
 * Cephes Math Library Release 2.2:  June, 1992
 * Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
 * Adapted to lcc-win by jacob navia

							qerfc.c	
0.25| 0.723673609831763067014931732235184287934646022107688745917528068737314916102697774814815194150258633209950741687586032406797084275370
0.5 | 0.479500122186953462317253346108035471263548424242036299941194274352806478283146429085211781265212242967033875613805608763934585309409
0.75| 0.288844366346484868401062165408589222625794045903462767721866602874963631235950438889206746109901099956469654629581357252086509804654
1.0 | 0.157299207050285130658779364917390740703933002033697091540062102165282745903989158738016674651855111545841738467978305635147660941745
 */

#include "qhead.h"

/* exp(x**2) * erfc(x) */

void EXPORT qerfc(Qfloatp x,Qfloatp y )
{
	Qfloat xx[1], z[1], a[1], b[1];
	int i, j;
	int sign;

	qmov( x, xx );
	qmul( xx, xx, z );	/* square of y				*/

	sign = xx[0].sign;

	if( xx[0].exponent < (QELT) (EXPONE + 2) ) /* xx < 4.0 */
	{
		i = 74; /*27;*/
		if( xx[0].exponent < EXPONE ) /* xx < 1.0 */
			i = 40; /*10;*/
		j = 4 * i + 1;
		itoq( j, a);		/*  2 * i  +  1			*/
		qmov( a, y );

		/* continued fraction expansion */
		while( j > 1 )
		{
			qsub( qone, a, a );
			qmuli( a, z, b );
			qdiv( y, b, y );
			qsub( qone, a, a );
			qadd( a, y, y );

			qsub( qone, a, a );
			qmuli( a, z, b );
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
		qsub( y, qone, y );
		/*qdiv( a, y, y );*/
	}


	/* alternate continued fraction expansion for large x */
	else
		{
		j = 192;
		if( xx[0].exponent > (QELT) (EXPONE + 1) )  /* xx > 2.0 */
			j = 220; /*76;*/
		itoq( j, a );
		if( sign != 0 )
			xx[0].sign = 0;
		qmov( xx, y );
		while( j > 0 )
		{
			qdiv( y, a, y );
			decreaseExponent(y); 
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
		if( sign != 0 )
			qsub( y, qtwo, y );
	}
}
