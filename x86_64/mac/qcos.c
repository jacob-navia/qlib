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
	/* cos(x) = sin( pi/2 - x ) */
	qsub( x, qPi_Div_2, y );
	qfsin( y, y );
}

#if NOTUSED

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

void qcosm1(QELT * x,QELT * y )
{
	QELT a[NQ], s[NQ], x2[NQ], n[NQ];
	int sign;

	qmul( x, x, x2 );
	qclear(s);
	qmov( qone, a );
	qmov( qone, n );
	sign = -1;
	do
	    {
		qadd( qone, n, n );
		qdiv( n, a, a );
		qmul( x2, a, a );
		if( sign > 0 )
			qadd( a, s, s );
		else
			qsub( a, s, s );
		sign = -sign;
		qadd( qone, n, n );
		qdiv( n, a, a );
	}
	while( ( exponent(s) -  exponent(a) ) < NBITS );

	qmov( s, y );
}
#endif
