 /*						qpow
 *
 *	Power function check routine
 *
 *
 *
 * SYNOPSIS:
 *
 * int qpow( x, y, z );
 * QELT *x, *y, *z;
 *
 * qpow( x, y, z );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes x raised to the yth power.
 *
 *       y
 *      x  =  exp( y log(x) ).
 *
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
 */

#include "qhead.h"
#define isinfq(x) (x->exponent >= MAXEXP)

void qfpow(Qfloatp x,Qfloatp y,Qfloatp z)
{
	Qfloat w[1];
	long long li;

	if (x->exponent >= MAXEXP) {
		qmov(x,z);
		return;
	}
	if (y->exponent >= MAXEXP) {
		qmov(y,z);
		return;
	}
	qfloor( y, w );
	if( qequal(y,w) ) {
		qifrac( y, &li, w );
		if( li < 0 )
			li = -li;
		if( li <= MAXEXP ) {
			qpowi( x, y, z );
			return;
		}
	}
	/* z = exp( y * log(x) ) */

	qflog( x, w );
	qmul( y, w, w );
	qfexp( w, z );
}


/*							qpowi
 *
 *	Real raised to integer power
 *
 *
 *
 * SYNOPSIS:
 *
 * int qpowi( x, n, z );
 * QELT *x, *n, *z;
 *
 * qpowi( x, n, z );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns argument x raised to the nth power, where n is integer valued.
 * The routine efficiently decomposes n as a sum of powers of
 * two. The desired power is a product of two-to-the-kth
 * powers of x.
 *
 */

void qpowi(Qfloatp x,Qfloatp y,Qfloatp z)
{
	Qfloat w[1];
	long long li, lx, e;
	int signx, signy;

	qifrac( y, &li, w );
	qmov( x, w );
	signx = signof(w);
	setpositive(w);
	if( li < 0 )
		lx = -li;
	else
		lx = li;

	/* Check for 2^N. */
	if( qequal( qtwo , w ) ) {
		qmov( qtwo, z );
		if( signx && (lx & 1) )
			setnegative(z);
		else
			setpositive(z);
		e = exponent(w) + li - 1;
		if( e > MAXEXP )
			qinfin( z );
		else if( e <= 0 ) {
			qclear(z);
		}
		else {
			z[0].exponent = e;
		}
		return;
	}

	if( lx == 0x7fffffff ) {
		qfpow( x, y, z );
		return ;
	}

	if( exponent(x) == 0 ) {
		if( li == 0 ) 	{ qmov( qone, z ); }
		else if( li < 0 ) { qinfin( z ); }
		else { qclear( z ); }
		return;
	}

	if( li == 0L ) {
		qmov( qone, z );
		return;
	}

	if( li < 0 ) {
		li = -li;
		signy = -1;
	}
	else signy = 0;


	/* First bit of the power */
	if( li & 1 ) {
		qmov( w, z );
	}
	else {
		qmov( qone, z );
		signx = 0;
	}

	li >>= 1;
	while( li != 0L ) {
		qmul( w, w, w );	/* arg to the 2-to-the-kth power */
		if( li & 1L )	/* if that bit is set, then include in product */
			qmul( w, z, z );
		li >>= 1;
	}

	if( signx )
		qneg( z ); /* odd power of negative number */

	if( signy ) {
		if( exponent(z) != 0 ) {
			//qdiv( z, qone, z );
			qinv(z,z);
		}
		else {
			qinfin( z );
			//mtherr( "qpowi", OVERFLOW );
		}
	}
}
