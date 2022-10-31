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
#if 1
// This is a modified sqrt2: sqrt2/2.
static Qfloat sqrt2[1] = {{
0,EXPONE-1,{0xb504f333f9de6484ULL,0x597d89b3754abe9fULL,
0x1d6f60ba893ba84cULL,0xed17ac8583339915ULL,0x4afc83043ab8a2c3ULL,
0xa8b1fe6fdc83db39ULL,0x0f74a85e439c7b4aULL}}};
#if 0
QfloatAccum qlog2Accum[1] = {
{0,EXPONE-1,{0,0xb17217f7d1cf79abULL,0xc9e3b39803f2f6afULL,
0x40f343267298b62dULL,0x8a0d175b8baafa2bULL,0xe7b876206debac98ULL,
0x559552fb4afa1b10ULL,0xed2eae35c1382144ULL,0x27573b291169b825ULL,0x3e96ca16224ae8c5ULL}}};
#endif

#include <string.h>
void qflog(const Qfloatp x,Qfloatp y)
{
	Qfloat xx[1], z[1], a[1], b[1], t[1];
	int ex,r,i=0;
	QELT den;

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
	if( ex == EXPONE ) { /* log 1 = 0 */
		if( qequal(x, qone) ) {
			qclear(y);
			return ;
		}
	}
	ex -= (EXPONE-1);
	xx[0].exponent = (EXPONE-1);
	/* Adjust range to 1/sqrt(2), sqrt(2) */
	if( qcmp( xx, sqrt2 ) < 0 ) {
		ex -= 1;
		xx[0].exponent += 1;
	}

	qincr( xx, b );
	qsub( qone, xx, a );
	if( exponent(a) == 0 ) {
		qclear(y);
		goto bdone;
	}
	qdiv( b, a, y );	/* store (x-1)/(x+1) in y */

	qsquare( y, z );

	qmov( z, a );
	den = 3;
	qmul(oneThird, z, t);
	qincr( t, b);
	do {
		//qadd_subtract( qtwo, qj, qj, ADDITION );	/* 2 * i + 1		*/
		den += 2;
		qmul( z, a, a );
		qdivi(den, a, t);
		i++;
		r = qadd_subtract( t, b, b, ADDITION );
	} while( r != 0 );

	qmul( b, y, y );
	y->exponent += 1;

bdone:
	/* now add log of 2**ex */
	if( ex != 0 ) {
		itoq( ex, b );
		qfma(C2,b,y,y);
		qfma(C1,b,y,y);
	}
}
#else
float128_t logf128(float128_t);
void qflog(const Qfloatp a,Qfloatp y)
{
	float128_t d;
	Qfloat x[1],z[1],w[1];

	if( signof(a) ) {
		qclear(y);
		mtherr( "qflog", DOMAIN );
		return;
	}
	if( exponent(a) == 0 ) {
		qinfin( y );
		y[0].exponent = -1;
		mtherr( "qflog", SING );
		return;
	}
	d = qtoe113(a);
	d = logf128(d);
	e113toq(d,x);

	x->sign = -1;
	qfexp(x,z);
	qmul(a,z,w);
	x->sign = 0;
	qadd(x,w,w);
	qsub(qone,w,x);

	x->sign = -1;
	qfexp(x,z);
	qmul(a,z,w);
	x->sign = 0;
	qadd(x,w,w);
	qsub(qone,w,y);

#if 0
	x->sign = -1;
	qfexp(x,z);
	qmul(a,z,w);
	x->sign = 0;
	qadd(x,w,w);
	qsub(qone,w,x);

	x->sign = -1;
	qfexp(x,z);
	qmul(a,z,w);
	x->sign = 0;
	qadd(x,w,w);
	qsub(qone,w,y);
#endif
}
#endif
#else
// Implementation using the AGM (Arithmetic Geometric Mean)
void qflog(Qfloatp x,Qfloatp y)
{
    Qfloat  z[1], qs[1], t[1], qm[1];
    int exponent,m;

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
    exponent = x->exponent-EXPONE;
    m = (NBITS+3)/2 - exponent;
    itoq(m,qm);

    qmov(x,qs);
    qs[0].exponent += m; // s = x * (2 ^ m)
    qinv(qs,qs);
    qs[0].exponent += 2; // 4 * 1/s --> 4/s
    qagm(qone,qs,z);     // z = AGM(1,s)
    z->exponent++;       // z = 2*AGM(1,4/s)
    qdiv(z,qpi,t);       // t = Pi/(2*AGM(1,4/s))

    qmul(qm,qlog2,qm);   // qm = m * log2
    qsub(qm,t,y);        // result = Pi/(2*AGM(1,4/s)) - m * log2
}

#endif
#if 0
extern QELT qpi[];
static QELT C1[NQ] = {
0x00000000,0x00080101,0x00000000,0x80000000,0x00000000};
// log(2) = 
//0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964
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
