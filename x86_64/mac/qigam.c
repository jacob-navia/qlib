/*							qigam.c
*	Check routine for incomplete gamma integral
*
 *
 *
 * SYNOPSIS:
*
 * For the left tail:
* int qigam( a, x, y );
* QELT *a, *x, *y;
* qigam( a, x, y );
*
 * For the right tail:
* int qigamc( a, x, y );
* QELT *a, *x, *y;
* qigamc( a, x, y );
*
 *
 * DESCRIPTION:
*
 * The function is defined by
*
 *                           x
*                            -
*                   1       | |  -t  a-1
*  igam(a,x)  =   -----     |   e   t   dt.
*                  -      | |
*                 | (a)    -
*                           0
*
 *
 * In this implementation both arguments must be positive.
* The integral is evaluated by either a power series or
* continued fraction expansion, depending on the relative
* values of a and x.
*
 *
 * ACCURACY:
*
 * Expansions terminate at less than full working precision.
*
 */

/*							qigam.c	*/
/*	Check routine for incomplete gamma integral */
/*	SLM, 22 Jan 84	*/



/*
* incomplete gamma integral
*
 *
 *          inf.
*           -
*   -      |   -t  a-1
*  | (a)   |  e   t   dt  =  qigamc(a,x)
*          |
*         -
*          x
*
 *
 */

#include "qhead.h"
#define DELTAEXP 355

void qigamc(Qfloatp a,Qfloatp x,Qfloatp y)
{
	Qfloat ans[1], c[1], yc[1], ax[1], z[1], pk[1], pkm1[1], pkm2[1], qk[1];
	Qfloat qkm1[1], qkm2[1], r[1], t[1];

	if( (x->sign != 0) || ( a->sign != 0) || (x->exponent == 0) || (a->exponent == 0) )
	{
		mtherr( "qigam", DOMAIN );
		return;
	}

	qsub( a, x, z);		/* z = x - a; */

	if( (x->exponent <= (QELT) (EXPONE-1)) || (z[0].sign != 0 ) )
	{
		qigam( a, x, y );
		qsub( y, qone, y );
		return;
	}

	/* 	ax = exp( a * log(x) - x - lgam(a) ); */
	qflog( x, ax );
	qmul( a, ax, ax );
	qsub( x, ax, ax );
	qlgam( a, c );
	qsub( c, ax, c );
	qfexp( c, ax );

	/* continued fraction */
	qsub( a, qone, y); 	/* y = 1.0 - a; */
	qadd( x, y, z );
	qadd( qone, z, z);	/* z = x + y + 1.0; */
	qmov( qzero, c );
	qmov( qone, pkm2 );	/* pkm2 = 1.0; */
	qmov( x, qkm2 );	/* qkm2 = x; */
	qadd( x, qone, pkm1);	/* pkm1 = x + 1.0; */
	qmul( z, x, qkm1);	/* qkm1 = z * x; */
	qdiv( qkm1, pkm1, ans);	/* ans = pkm1/qkm1; */

	do
	    {
		qadd( qone, c, c);		/* c += 1.0; */
		qadd( qone, y, y);		/* y += 1.0; */
		qadd( qtwo, z, z);		/* z += 2.0; */
		qmul( y, c, yc );		/* yc = y * c; */
		qmul( pkm2, yc, r);
		qmul( pkm1, z, pk);
		qsub( r, pk, pk );	/* pk = pkm1 * z  -  pkm2 * yc; */
		qmul( qkm2, yc, r );
		qmul( qkm1, z, qk );
		qsub( r, qk, qk );	/* qk = qkm1 * z  -  qkm2 * yc; */
		if( qk[0].exponent > 0 )
		{
			qdiv( qk, pk, r );	/* r = pk/qk; */
			qsub( r, ans, t );	/* t = ans - r */
			qmov( r, ans );		/* ans = r; */
		}
		else
			qmov( qone, t );		/* t = 1.0; */

		qmov( pkm1, pkm2 );		/* pkm2 = pkm1; */
		qmov( pk, pkm1 );		/* pkm1 = pk; */
		qmov( qkm1, qkm2 );		/* qkm2 = qkm1; */
		qmov( qk, qkm1 );		/* qkm1 = qk; */
	}
	while( ans[0].exponent -  t[0].exponent < DELTAEXP ); /* was 67 10**-20 */

	qmul( ax, ans, y );  /* return ans * ax */
}


void qigam(Qfloatp a,Qfloatp x, Qfloatp y)
{
	Qfloat z[1],ax[1],c[1],r[1],ans[1];
	if( (x->sign != 0) || ( a->sign != 0) || (x->exponent == 0) || (a->exponent == 0) )
	{
		mtherr( "qigam", DOMAIN );
		return;
	}

	qsub( a, x, z);		/* z = x - a; */

	if( (x->exponent > (QELT) (EXPONE-1)) && (z->sign == 0 ) )
	{
		qigamc( a, x, y );
		qsub( y, qone, y );
		return;
	}
	/* 	ax = exp( a * log(x) - x - lgam(a) ); */
	qflog( x, ax );
	qmul( a, ax, ax );
	qsub( x, ax, ax );
	qlgam( a, c );
	qsub( c, ax, c );
	qfexp( c, ax );

	/* power series */
	qmov( a, r );		/* r = a; */
	qmov( qone, c );	/* c = 1.0; */
	qmov( qone, ans );	/* ans = 1.0; */

	do
	    {
		qadd( qone, r, r );		/* r += 1.0; */
		qdiv( r, x, z );
		qmul( z, c, c );		/* c *= x/r; */
		qadd( c, ans, ans );		/* ans += c; */
	}
	while( ans[0].exponent -  c[0].exponent < DELTAEXP ); /* was 67 while( c/ans > stop ); */

	qdiv( a, ax, z );	/* ans * ax / a  */
	qmul( z, ans, y );
}
