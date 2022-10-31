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
 * Cephes Math Library Release 2.1:  January, 1989
 * Copyright 1984, 1987, 1989, 1998 by Stephen L. Moshier
 * Reworked 2018 by J.N.
*/

#include "qhead.h"
#include "mconf.h"

// qcl = log(1.0/sqrt(2*pi))
//  -0.918938533204672741780329736405617639861397473637783412817151540482765695927260
static Qfloat qcl[1]= {
{-1,EXPONE-1,{0xeb3f8e4325f5a534ULL,0x94bc900144192023ULL,0xcfb08f8d13458b4dULL,
	0xdec6a3133daa155dULL,0x212f9d7fe00e86bfULL,0x93eabf905c5569bbULL,0xb05cab571b4cda5bULL}}};
// qcc = 1.0/sqrt(2*pi)
// 0.39894228040143267793994605993438186847585863116493465766592582967065792589930183850125
static Qfloat qcc[1] = {
{0,EXPONE-2, {0xcc42299ea1b28468ULL, 0x7e59e2805d5c717fULL, 0xa7053e60a24e0db1ULL,
0x808e369004c73a5dULL, 0x3bf158aa88fa7b6eULL, 0x47eba75766e9d11bULL, 0xd6284ef2e17c2a26ULL}}};

extern double ndtri(double);

void qndtri(Qfloatp qy1,Qfloatp qx0)
{
	double y0, x0;
	int i, k, righttail;
	Qfloat qy0[1], qd[1], qy[1], qx[1], temp[1];
	
	qmov( qy1, qy0 );
	
	if( qy0->sign || qy0->exponent == 0 ) {
		mtherr( "qndtri", DOMAIN );
		qinfin( qx0 );
		qneg( qx0 );
		return;
	}
	
	if( qcmp(qy0, qone) >= 0 ) {
		mtherr( "qndtri", DOMAIN );
		qinfin( qx0 );
		return ;
	}
	
	/* Avoid a convergence problem that happens when y is close to 1.  */
	if( qcmp(qy0, qhalf) >= 0 ) {
		qsub( qy0, qone, qy0 );
		righttail = 1;
	}
	else {
		righttail = 0;
	}
	
	if( qy0[0].exponent < (QELT) (EXPONE - 1021) ) /* 4.5e-308 */ {
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
	else {
		y0 = qtoe( qy0, 1 );
		x0 = ndtri( y0 );
		etoq( x0, qx );
		k = 4;
	}
	
	
	for( i=0; i<k; i++ ) {
		qndtr( qx, qy );
		qsquare( qx, qd );
		qd[0].exponent -= 1;
		qneg(qd);
		qfexp( qd, qd );
		qmul( qcc, qd, qd );
		if( qd[0].exponent > 3 ) {
			qsub( qy0, qy, temp );
			qdiv( qd, temp, temp );
			qsub( temp, qx, qx );
		}
	}
	qmov( qx, qx0 );	
	if( righttail )
		qneg( qx0 );
}
