/*							qfcos.c
 *
 *	Circular cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * int qfcos( x, y );
 * QELT *x, *y;
 *
 * qfcos( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * cos(x) = sin(pi/2 - x)
 *
 */
/*							qfcos.c		*/

/* cosine check routine */

#include "qhead.h"

void qfcos(Qfloatp x,Qfloatp y)
{
	Qfloat a[1];
	/* cos(x) = sin( pi/2 - x ) */
	qsub( x, qPi_Div_2, a );
	qfsin( a, y );
}


/*							qcosm1
 *
 *	Complemented circular cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * int qcosm1( x, y );
 * QELT *x, *y;
 *
 * qcosm1( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns cos(x) - 1.
 * Computation by Taylor series, useful only for small x.
 *
 */

void qcosm1(Qfloatp x,Qfloatp y )
{
	Qfloat a[1], s[1], x2[1];
	int sign,n=1;

	qsquare( x, x2 );
	qclear(s);
	qmov( qone, a );
	sign = -1;
	do {
		n++;
		qdivi( n, a, a );
		qmul( x2, a, a );
		if( sign > 0 )
			qadd( a, s, s );
		else
			qsub( a, s, s );
		sign = -sign;
		n++;
		qdivi( n, a, a );
	}
	while( ( exponent(s) -  exponent(a) ) < NBITS );

	qmov( s, y );
}
