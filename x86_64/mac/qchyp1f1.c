/*							qhyp.c		*/

/*							qchyp1f1.c	*/

/*	confluent hypergeometric function
*
 *                          1           2
*                       a x    a(a+1) x
*   F ( a,b;x )  =  1 + ---- + --------- + ...
*  1 1                  b 1!   b(b+1) 2!
*
 *
 * Series summation terminates at 70 bits accuracy.
*
 */

#include "qhead.h"

extern qcmplx qcone;

static qcmplx an, bn, a0, sum, n, t;

int qchyp1f1(qcmplx *a,qcmplx *b,qcmplx *x,qcmplx *y)
{
	int ae, se;

	qmov( qone, &qcone.re[0] );
	qclear( &qcone.im[0] );
	qcmov( a, &an );		/*an = a;*/
	qcmov( b, &bn );		/*bn = b;*/
	qcmov( &qcone, &a0 );	/*a0 = 1.0;*/
	qcmov( &qcone, &sum );	/*sum = 1.0;*/
	qcmov( &qcone, &n );	/*n = 1.0;*/

	do
	    {
		if( an.re[0].exponent == 0 && an.im[0].exponent == 0 )
			goto done;
		if( bn.re[0].exponent == 0 && bn.im[0].exponent == 0 )
		{
			qinfin(sum.re);
			qinfin(sum.im);
			goto done;
		}
		/*
		if( (a0 > 1.0e34) || (n > 130) )
		goto asymf;
		*/
		qcmul( &bn, &n, &t );
		qcdiv( &t, &a0, &a0 );
		qcmul( &an, x, &t );
		qcmul( &t, &a0, &a0 );	/*a0 *= (an * x) / (bn * n);*/
		qcadd( &sum, &a0, &sum );	/*sum += a0;*/
		qcadd( &an, &qcone, &an );	/*an += 1.0;*/
		qcadd( &bn, &qcone, &bn );	/*bn += 1.0;*/
		qcadd( &n, &qcone, &n );	/*n += 1.0;*/
		ae = a0.re[0].exponent;
		if (a0.im[0].exponent > ae)
			ae = a0.im[0].exponent;
		se = sum.re[0].exponent;
		if (sum.im[0].exponent > se)
			se = sum.im[0].exponent;
	}
	while( ((int) ae - (int) se) > -70 );

	/*printf("1F1( %.2E %.2E %.5E ) = %.3E %.6E\n", a, b, x,  n, sum);*/
done:
	qcmov( &sum, y );
	return(0);
}
