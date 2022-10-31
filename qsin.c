/*							qfsin.c
 *      Circular sine check routine

 *
 *
 * SYNOPSIS:
 *
 * int qfsin(Qfloatp x, Qfloatp y );
 *
 * qfsin( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of pi/2.
 * Then
 *
 *               3    5    7
 *              z    z    z
 * sin(z) = z - -- + -- - -- + ...
 *              3!   5!   7!
 *
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1985, 1995, 1996 by Stephen L. Moshier
Adapted to lcc-win by jacob, 2018
 */

#include "qhead.h"
static Qfloat one_sixth[1] = {
{0,0x0007fffe,{0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,
0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,
0xaaaaaaaaaaaaaaabULL}}};
void qfsin(Qfloatp x,Qfloatp y)
{
	int sign,r;
	long long mod;
	unsigned long long A;
	Qfloat a[1], b[1], z[1], xx[1];

	if( signof(x) )
		sign = -1;
	else
		sign = 1;
	qmov( x, xx );
	setpositive(xx);
	/* range reduction to [0, pi/2]	*/
	qdiv( qPi_Div_2, xx, a );
	qfloor( a, xx );
	/* b = xx - 8 * floor(xx/8) */
	if( exponent(xx) >= 3) {
		xx[0].exponent -= 3;
		qfloor( xx, b );
		b[0].exponent += 3;
		xx[0].exponent += 3;
	}
	else
		qclear(b);

	qsub( b, xx, b );
	qifrac( b, &mod, b );

	qsub( xx, a, b );
	qmul( b, qPi_Div_2, xx );

	mod &= 3;
	if( mod > 1 )
		sign = -sign;
	if( mod & 1 )
		qsub( xx, qPi_Div_2, xx );	/* xx = 1 - xx */

	qsquare( xx, z );
	qneg( z );

	//qmov( qthree, a );
	A = 3;
	qmul(z,one_sixth,b);
	qincr(b,y);

	/* power series */
	do {
		//qadd( qone, a, a );	/* a += 1	*/
		A++;
		qdivi( A, b, b );	/* b /= a	*/
		//qadd( qone, a, a );	/* a += 1	*/
		A++;
		qdivi( A, b, b );	/* b /= a	*/
		qmul( z, b, b );	/* b *= z	*/
		r = qadd( b, y, y );	/* y += b	*/
	}
	while( r != 0 );
	/* The loop count abpve never exceeeds 100 
		if (A > 100) printf("%ld\n",A); // Never prints anything
	*/
	qmul( xx, y, y );
	if( sign < 0)
		setnegative(y);
}

/* sin(x) - x */

int qsinmx3(Qfloatp x,Qfloatp y )
{
	Qfloat z[1], b[1], n[1];
	int r;

	qsquare( x, z );
	qneg( z );

	qmov( qone, n );
	qclear( y );

	/* compute the cube term x^3/3! */
	qmov( x, b );
	qincr( n, n );
	qdiv( n, b, b );     /* x / 2 */
	qincr( n, n );
	qdiv( n, b, b );     /* x / 6 */
	qmul( z, b, b );     /* x^3 / 6 */

	/* power series */
	do
		{
		qincr( n, n );	/* n += 1	*/
		qdiv( n, b, b );	/* b /= n	*/
		qincr( n, n );	/* n += 1	*/
		qdiv( n, b, b );	/* b /= n	*/
		qmul( z, b, b );	/* b *= z	*/
		r = qadd( b, y, y );	/* y += b	*/
	}
	while(r != 0 );
	return 0;
}
