/*							qndtri.c
*
 *	Inverse of Normal distribution function
*
 *
 *
 * SYNOPSIS:
*
 * int qndtri(y, x);
* QELT *y, *x;
*
 * qndtri(y, x);
*
 *
 *
 * DESCRIPTION:
*
 * Returns the argument, x, for which the area under the
* Gaussian probability density function (integrated from
* minus infinity to x) is equal to y.
*
 * The routine refines a trial solution computed by the double
* precision function ndtri.
*
 */


/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989, 1998 by Stephen L. Moshier
*/

#include <stats.h>
#include "qhead.h"

extern double PI;

void qndtri(Qfloatp qy1,Qfloatp qx0)
{
	long double y0, x0;
	int i, k, righttail;
	Qfloat temp[1];
	Qfloat qy0[1];
	static Qfloat qcc[1];
	static Qfloat qcl[1];
	Qfloat qd[1];
	Qfloat qy[1];
	Qfloat qx[1];
	static int qcflg = 0;

	if( qcflg == 0 )
	{
		qmov( qpi, temp );
		temp[0].exponent += 1;
		qfsqrt( temp, temp );
		qdiv( temp, qone, qcc );
		qflog( qcc, qcl );
		qcflg = 1;
	}

	qmov( qy1, qy0 );

	if( qcmp(qy0, qzero) <= 0 )
	{
		mtherr( "qndtri", DOMAIN );
		qinfin( qx0 );
		qneg( qx0 );
		return;
	}

	if( qcmp(qy0, qone) >= 0 )
	{
		mtherr( "qndtri", DOMAIN );
		qinfin( qx0 );
		return;
	}

	/* Avoid a convergence problem that happens when y is close to 1.  */
	if( qcmp(qy0, qhalf) >= 0 )
	{
		qsub( qy0, qone, qy0 );
		righttail = 1;
	}
	else
	    {
		righttail = 0;
	}

	if( qy0[0].exponent < (QELT) (EXPONE - 1021) ) /* 4.5e-308 */
	{
		k = 7;
		/* x = sqrt( -2 log y ) */
		qflog( qy0, qd );
		qd[0].exponent += 1;
		qd[0].sign = 0;
		qfsqrt( qd, qx );
		/* fine adjustment:
		* x = x - (log x  +  log sqrt 2pi)/x
		*/
		qflog( qx, temp );
		qsub( qcl, temp, temp );
		qdiv( qx, temp, temp );
		qsub( qx, temp, qx );
		qx[0].sign = -1;
	}
	else
	{
		__qtoe64( qy0, (unsigned short *) &y0 );
		x0 = normal_distribution_inv( y0 );
		__e64toq( (unsigned short *) &x0, qx );
		k = 3;
	}


	for( i=0; i<k; i++ )
	{
		qndtr( qx, qy );
		/* debugging code */
		/*
		qtoasc( qx, s, 5 );
		printf( "%s", s );
		qtoasc( qy, s, 5 );
		printf( " %s", s );
		qsub( qy0, qy, temp );
		qdiv( qy0, temp, temp );
		qtoasc( temp, s, 5 );
		printf( " %s\n", s );
		*/
		/*   */
		qmul( qx, qx, qd );
		qd[0].exponent -= 1;
		qneg(qd);
		qfexp( qd, qd );
		qmul( qcc, qd, qd );
		if( qd[0].exponent > 3 )
		{
			qsub( qy0, qy, temp );
			qdiv( qd, temp, temp );
			qsub( temp, qx, qx );
		}
	}

	qmov( qx, qx0 );
	if( righttail )
		qneg( qx0 );
}
