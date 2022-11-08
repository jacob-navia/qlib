/*							qflog.c
*
 *	Natural logarithm
*
 *
 *
 * SYNOPSIS:
*
 * int qflog( x, y );
* QELT *x, *y;
*
 * qflog( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the base e (2.718...) logarithm of x.
*
 * After reducing the argument into the interval [1/sqrt(2), sqrt(2)],
* the logarithm is calculated by
*
 *       x-1
* w  =  ---
*       x+1
*                     3     5
*                    w     w
* ln(x) / 2  =  w + --- + --- + ...
*                    3     5
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier


Adapted to lcc-win by jacob navia
*/

#include "qhead.h"

#if 1
// This is a modified sqrt2: sqrt2/2.
static Qfloat sqrt2[1] = {
0,EXPONE-1,0xb504f333f9de6484ULL,0x597d89b3754abe9fULL,
0x1d6f60ba893ba84cULL,0xed17ac8583339915ULL,0x4afc83043ab8a2c3ULL,
0xa8b1fe6fdc83db39ULL,0x0f74a85e439c7b4aULL};
#include <string.h>
void qflog(Qfloatp x,Qfloatp y)
{
	Qfloat xx[1], z[1], a[1], b[1], t[1], qj[1];
	int ex;
	int i;

	if( signof(x) ) {
		qclear(y);
		mtherr( "qflog", DOMAIN );
		return;
	}
	if( exponent(x) == 0 ) {
		qinfin( y );
		y[0].exponent = -1;
		mtherr( "qflog", SING );
		return;
	}
	/* range reduction: log x = log( 2**ex * m ) = ex * log2 + log m */
	qmov(x, xx );
	ex = exponent(xx);
	if( ex == EXPONE )
	{ /* log 1 = 0 */
		if( qequal(x, qone) )
		{
			qclear(y);
			return ;
		}
	}
	ex -= (EXPONE-1);
	xx[0].exponent = (EXPONE-1);
	/* Adjust range to 1/sqrt(2), sqrt(2) */
	if( qcmp( xx, sqrt2 ) < 0 )
	{
		ex -= 1;
		xx[0].exponent += 1;
	}

	qadd( qone, xx, b );
	qsub( qone, xx, a );
	if( exponent(a) == 0 ) {
		qclear(y);
		goto bdone;
	}
	qdiv( b, a, y );	/* store (x-1)/(x+1) in y */

	qmul( y, y, z );

	qmov( z, a );
	qmov( qthree, qj );
	qdiv( qthree, z, t);
	qadd( t, qone,b);
	i=0;
	do {
		qadd_subtract( qtwo, qj, qj, ADDITION );	/* 2 * i + 1		*/
		qmul( z, a, a );
		qdiv( qj, a, t );
		qadd_subtract( t, b, b, ADDITION );
		i++;
	}
	while( ( (int)exponent(b) -  (int)exponent(t) ) <= NBITS );

	qmul( b, y, y );
	y->exponent += 1;

bdone:
	/* now add log of 2**ex */
	if( ex != 0 ) {
		itoq( ex, b );
		qmul(b,qlog2,t);
#if 0
		qmul( C2, b, t );
		qadd( t, y, y );
		qmul( C1, b, t );
#endif
		qadd( t, y, y );
	}
	;
}
#else
extern QELT qpi[];
static QELT C1[NQ] = {
0x00000000,0x00080101,0x00000000,0x80000000,0x00000000};
static QELT C2[NQ] = {
 0x00000000,0x00080008,0x00000000,0xb17217f7,0xd1cf79ab,0xc9e3b398,0x03f2f6af,
0x40f34326,0x7298b62d,0x8a0d175b,0x8baafa2b,0xe7b87620,0x6debac98,0x559552fb,
0x4afa1b10,0xed2eadde
};
static QELT qfour[NQ] = { 0,0x80003,0,0x80000000,0};
void qflog(QELT *a,QELT *y)
{
	int i=0;
	QELT an[NQ],bn[NQ],tmp[NQ],delta[NQ];

qfloat qa,qan,qbn;

	qmov(qone,an);
	qmul(C1,a,bn);
	qdiv(bn,qfour,bn);
	do {
		qmov(an,tmp);
		qadd(an,bn,an);
		an[1] -= 1;
		qmul(tmp,bn,bn);
		qfsqrt(bn,bn);
		qsub(an,bn,delta);
		i++;
		if (i > 30) {
			qmov(a,&qa);
			qmov(an,&qan);
			qmov(bn,&qbn);
printf("not converging!\n%130.123qe\n%130.123qe\n%130.123qe\n%130.123qe\n",qa,qan,qbn,*(qfloat *)delta);
			break;
		}
	} while (i<20);
	qadd(an,an,an);
	qdiv(an,qpi,an);
	qsub(C2,an,y);
}
#endif
