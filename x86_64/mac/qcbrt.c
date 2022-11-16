/*							qcbrt.c
*
 *	Cube root
*
 *
 *
 * SYNOPSIS:
*
 * int qcbrt( x, y );
* QELT *x, *y;
*
 * qcbrt( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the cube root of the argument, which may be negative.
*
Cephes Math Library Release 2.2:  January, 1991
Copyright 1984, 1991 by Stephen L. Moshier
*/


/* cube root */
/* 21 December 83, SLM */

#include "qhead.h"

#define CBRT2 1.2599210498948731647672106L
#define CBRT4 1.5874010519681994747517056L
#define CBRT2I 0.793700525984099737375852L
#define CBRT4I 0.629960524947436582383605L

void qcbrt(Qfloatp xx,Qfloatp y)
{
	int sign,i;
	long e, rem;
	long double x, z;
	Qfloat *a,tmp[1],b[1],c[1];


	if( exponent(xx) == 0 ) {
		qclear( y );
		return;
	}

	if (xx == y)
		a = tmp;
	else
		a = y;
	sign = signof(xx);
	qmov( xx, c );
	setpositive(c);

	/* extract power of 2, leaving mantissa between 0.5 and 1 */
	qfrexp( c, &e, a );

	__qtoe64( a, (unsigned short *) &z );
	/* Approximate cube root of number between .5 and 1,
	* peak relative error = 9.2e-6
	*/
	x = (((-1.3466110473359520655053e-1  * z
			+ 5.4664601366395524503440e-1) * z
			- 9.5438224771509446525043e-1) * z
			+ 1.1399983354717293273738e0 ) * z
		+ 4.0238979564544752126924e-1;
	/* Newton iterations. */
	x -= ( x - (z/(x*x)) )*0.33333333333333333333333L;
	x -= ( x - (z/(x*x)) )*0.33333333333333333333333L;
	x -= ( x - (z/(x*x)) )*0.33333333333333333333333L;
	x -= ( x - (z/(x*x)) )*0.33333333333333333333333L;

	/* exponent divided by 3 */
	if( e >= 0 ) {
		rem = e;
		e = e / 3;
		rem = rem - 3 * e;
		if( rem == 1 )
			x = x * CBRT2;
		else if( rem == 2 )
			x = x * CBRT4;
	}
	else { /* argument less than 1 */
		e = -e;
		rem = e;
		e /= 3;
		rem -= 3*e;
		if( rem == 1 )
			x *= CBRT2I;
		else if( rem == 2 )
			x *= CBRT4I;
		e = -e;
	}
	__e64toq( (unsigned short *) &x, a );

	/* multiply by power of 2 */
	qldexp( a, e, a );

	for (i=0; i<3;i++) {
		/* More Newton iterations. */
		qmul( a, a, b );	/* 	x * x		*/
		qdiv( b, c, b );	/*	z / x*x	*/
		qsub( b, a, b );	/*	x - z/x*x	*/
		qdiv( qthree, b, b );	/*	.../3		*/
		qsub( b, a, a );	/*	x - ...		*/
	}
	if( sign != 0 )
		setnegative(a);
	if (xx == y)
		qmov( a, y );
}
