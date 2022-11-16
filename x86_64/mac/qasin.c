/*							qasin.c
*
 *	Inverse circular sine
*
 *
 *
 * SYNOPSIS:
*
 * int qasin( x, y );
* Qfloatp *x, *y;
*
 * qasin( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns radian angle between -pi/2 and +pi/2 whose sine is x.
*
 *    asin(x) = arctan (x / sqrt(1 - x^2))
*
 * If |x| > 0.5 it is transformed by the identity
*
 *    asin(x) = pi/2 - 2 asin( sqrt( (1-x)/2 ) ).
*
 */
/*							qacos
*
 *	Inverse circular cosine
*
 *
 *
 * SYNOPSIS:
*
 * int qacos( x, y );
* Qfloatp x[], y[];
*
 * qacos( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns radian angle between 0 and pi whose cosine
* is x.
*
 * acos(x) = pi/2 - asin(x)
*
Cephes Math Library Release 2.3:  April, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

/* arc sine */

#include "qhead.h"

void qasin(Qfloatp x,Qfloatp y)
{
	Qfloat a[1], z[1], zz[1], temp[1];
	int sign, flg;

	sign = signof(x);
	qmov( x, a );
	setpositive(a);
	if( qcmp( a, qone ) > 0 ) {
		mtherr( "qasin", DOMAIN );
		qclear(y);
		return;
	}
	if( qcmp( a, qhalf) > 0 ) {
		qsub( a, qhalf, zz );	/*zz = 0.5 -a;*/
		qadd( qhalf, zz, zz );	/*zz = (zz + 0.5)/2.0;*/
		if( exponent(zz) > 0 ) {
			zz->exponent -= 1 ;
		}
		qfsqrt( zz, z );		/* z = sqrt( zz );*/
		flg = 1;
	}
	else
	{
		qmov( a, z );		/*z = a;*/
		qmul( z, z, zz );	/* zz = z * z;*/
		flg = 0;
	}

	qsub( zz, qone, temp );	/* 1 - x**2 */
	qfsqrt( temp, temp );
	qdiv( temp, z, temp );
	qatn( temp, y );
	if( flg != 0 ) {
		y->exponent += 1;
		qsub( y, qPi_Div_2, y );	/* z = PIO2 - z;*/
	}
	y->sign = sign;
}



void qacos(Qfloatp x,Qfloatp y)
{
	Qfloat r[1],temp[1];

	qmov( x, temp );
	setpositive(temp);
	if( qcmp( temp, qone ) > 0 )
	{
		mtherr( "qacos", DOMAIN );
		qclear(y);
		return;
	}
	qasin( x, r );
	qsub( r, qPi_Div_2, y );
}
