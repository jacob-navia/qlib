/*							qatn
 *
 *	Inverse circular tangent
 *      (arctangent)
 *
 *
 *
 * SYNOPSIS:
 *
 * int qatn( x, y );
 * QELT *x, *y;
 *
 * qatn( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns radian angle between -pi/2 and +pi/2 whose tangent
 * is x.
 *
 * Range reduction is from three intervals into the interval
 * from zero to pi/8.
 *
 *                     2     2     2
 *               x    x   4 x   9 x
 * arctan(x) =  ---  ---  ----  ----  ...
 *              1 -  3 -  5 -   7 -
 *
 *
 *	Quadrant correct inverse circular tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * int qatn2( y, x, z );
 * QELT *x, *y, *z;
 *
 * qatn2( y, x, z );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns radian angle -PI < z < PI whose tangent is y/x.
 *
Cephes Math Library Release 2.3:  April, 1995
Copyright 1984, 1995 by Stephen L. Moshier
 */

/* arctangent check routine */

#include "qhead.h"


/* 
  This constants have been verified with  the PARI gp calculator J.N.
  tan(pi/8) = sqrt(2) - 1 =
 * 4.1421356237309504880168872420969807856967187537694807317667973799073248E-1
 */
static Qfloat qtp8[1] = {{
	0,EXPONE-2,{0xd413cccfe7799211ULL,0x65f626cdd52afa7cULL,
	0x75bd82ea24eea133ULL,0xb45eb2160cce6455ULL,0x2bf20c10eae28b0eULL,
	0xa2c7f9bf720f6ce4ULL,0x3dd2a1790e71ed2aULL // rounded
}}};

/* tan(3pi/8) = sqrt(2) + 1 =
 * 2.4142135623730950488016887242096980785696718753769480731766797379907325E0
 */
static Qfloat qt3p8[1] = {{
	0,EXPONE+1,{0x9a827999fcef3242ULL,0x2cbec4d9baa55f4fULL,
	0x8eb7b05d449dd426ULL,0x768bd642c199cc8aULL,0xa57e41821d5c5161ULL,
	0xd458ff37ee41ed9cULL,0x87ba542f21ce3da5ULL
}}};

void qatn(Qfloatp x,Qfloatp y)
{
	int i, j, nsq;
	int sign;
	Qfloat z[1], a[1], b[1], xx[1], qj[1], yy[1];

	qmov( x, xx );
	if( signof(xx) != 0 ) {
		setpositive(xx);
		sign = -1;
	}
	else
		sign = 1;

	/* range reduction */
	if( qcmp(xx, qt3p8) > 0 ) {
		qmov( qpi, yy );
		decreaseExponent(yy);
		//qdiv( xx, qone, xx );
		qinv(xx,xx);
		qneg( xx );
	}
	else if( qcmp(xx, qtp8) > 0 ) {
		qmov( qpi, yy );
		yy[0].exponent -= 2;
		qsub( qone, xx, a );	/* x = (x-1.0)/(x+1.0) */
		qincr( xx, b );
		qdiv( b, a, xx );
	}
	else { qclear( yy ); }

	qsquare( xx, z );	/* square of x				*/
	if( exponent(z) == 0 ) {
		qmov( xx, y );
		goto done;
	}
	/* loop count for full convergence
	* x < sqrt(2)-1: i = 2*NBITS/9
	* x < 1: i = 4*NBITS/5
	*/

	i = 2*NBITS/9;
	j = 2 * i + 1;
	itoq( j, qj );		/*  2 * i  +  1			*/
	qmov( qj, b );

	/* continued fraction expansion */
	while( j > 1 ) {
		nsq = i * i;
		itoq( nsq, a );	/* i**2				*/
		qmuli( a, z, a );	/* i**2 * x**2			*/
		qdiv( b, a, b );	/*  i**2 x**2 / (2*i + 1)	*/
		j -= 2;
		i -= 1;
		qsub( qtwo, qj, qj );	/* 2*i + 1			*/
		qadd( qj, b, b );
	}

	qdiv( b, xx, y );

done:
	qadd( yy, y, y );

	if( sign < 0 )
		qneg(y);
}

/*							qatn2	*/
/* angle whose tangent is y/x */

void qatn2(Qfloatp y,Qfloatp x,Qfloatp z)
{
	Qfloat v[1], w[1];
	int code;


	code = 0;

	if( (signof(x) != 0) && (exponent(x) > 0) )
		code = 2;
	if( (signof(y) != 0) && (exponent(y) > 0 ) )
		code |= 1;

	if( exponent(x) <= 1 ) /* x zero */
	{
		if( code & 1 ) /* y negative */
		{
#if ANSIC
			qmov (qpi, z); /* - pi/2 */
			decreaseExponent(z);
			qneg(z);
#else
			qmov( qpi, z ); /* 3*pi/2 */
			decreaseExponent(z);
			qadd( qpi, z, z );
#endif
			return;
		}
		if( exponent(y) <= 1 ) /* y zero */
		{
			qclear(z);
			return;
		}
		qmov( qpi, z ); /* y positive */
		decreaseExponent(z);		/* PI/2 */
		return;
	}

	if( exponent(y) <= 1 ) /* y zero */
	{
		if( code & 2 ) /* x negative */
		{
			qmov( qpi, z );
			return;
		}
		qclear(z);
		return;
	}


	switch( code )
	{
	default:
	case 0:
	case 1:
		qclear(w);
		break;
	case 2:
		qmov( qpi, w );
		break;
	case 3:
		qmov( qpi, w );
		qneg(w);
		break;
	}

	qdiv( x, y, v );	/* z = w + arctan( y/x ) */
	qatn( v, z );
	qadd( w, z, z );
}


