/*							qhyp.c
*
 *	Confluent hypergeometric function
*
 *
 *
 * SYNOPSIS:
*
 * int qhyp( a, b, x, y );
* QELT *a, *b, *x, *y;
*
 * qhyp( a, b, x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Computes the confluent hypergeometric function
*
 *                          1           2
*                       a x    a(a+1) x
*   F ( a,b;x )  =  1 + ---- + --------- + ...
*  1 1                  b 1!   b(b+1) 2!
*
 *
 * ACCURACY:
*
 * Series expansion is truncated at less than full working precision.
*
 */

/*
Cephes Math Library Release 2.1:  November, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
*/

/*	hyperg.c	22 Oct 83  */

/*	confluent hypergeometric function
*
 *                          1           2
*                       a x    a(a+1) x
*   F ( a,b;x )  =  1 + ---- + --------- + ...
*  1 1                  b 1!   b(b+1) 2!
*
 */

#include "qhead.h"

void qhyp(Qfloatp a,Qfloatp b,Qfloatp x,Qfloatp y)
{
	int count=0;
	Qfloat an[1], bn[1], a0[1], sum[1], n[1], t[1];

	qmov( a, an );		/*an = a;*/
	qmov( b, bn );		/*bn = b;*/
	qmov( qone, a0 );	/*a0 = 1.0;*/
	qmov( qone, sum );	/*sum = 1.0;*/
	qmov( qone, n );	/*n = 1.0;*/

	do
		{
		if( an[0].exponent == 0 )
			goto done;
		if( bn[0].exponent == 0 )
		{
			qinfin(sum);
			goto done;
		}
		/*
		if( (a0 > 1.0e34) || (n > 130) )
		goto asymf;
		*/
		qmul( bn, n, t );
		qdiv( t, a0, a0 );
		qmul( an, x, t );
		qmul( t, a0, a0 );	/*a0 *= (an * x) / (bn * n);*/
		qadd( sum, a0, sum );	/*sum += a0;*/
		qincr(an,an,0);         /* an += 1.0 */
		qincr(bn,bn,0);         /* bn += 1.0 */
		qincr(n,n,0);           /* , += 1.0 */
		count++;
		if (count > 60)
			break;
	}
	// Here was -70 originally, but that gives only 20 digits!
	// Added cutoff at 60 iterations and increased the iteration count.
	while( ((int) a0[0].exponent - (int) sum[0].exponent) > -280 ); 

	/*printf("1F1( %.2E %.2E %.5E ) = %.3E %.6E\n", a, b, x,  n, sum);*/
done:
	qmov( sum, y );
}
