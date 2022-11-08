/*							qfexp.c
 *
 *	Exponential function check routine
 *
 *
 *
 * SYNOPSIS:
 *
 * int qfexp( x, y );
 * QELT *x, *y;
 *
 * qfexp( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns e (2.71828...) raised to the x power.
 *
 * Cephes Math Library Release 2.3:  January, 1995
 * Copyright 1984, 1995 by Stephen L. Moshier
 * Adapted to lcc-win by jacob navia
 */

#include "qhead.h"

/* C1 + C2 = ln 2 */
static Qfloat C1[1] = {
	0,EXPONE-1,0xb172000000000000ULL,0,0,0,0,0,0};
static Qfloat C2[1] = {
0,0x0007ffed,0xbfbe8e7bcd5e4f1dULL,0x9cc01f97b57a079aULL,
0x193394c5b16c5068ULL,0xbadc5d57d15f3dc3ULL,0xb1036f5d64c2acaaULL,
0x97da57d0d8876975ULL,0x71ae09c10a213ab9ULL
};

/* 363408 */
static Qfloat maxExp[1] =
{0,0x00080013,0xb172000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
 0x0000000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,};

void qfexp(Qfloatp x,Qfloatp y )
{
	Qfloat num[1], den[1], x2[1];
	long long i=0;
	int sign,exp;

	/* range reduction theory: x = i + f, 0<=f<1;
	* e^x = e^i * e^f
	* e^i = 2^(i/log 2).
	* Let i/log2 = i1 + f1, 0 <= f1 < 1.
	* Then e^i = 2^i1 * 2^f1, so
	* e^x = 2^i1 * e^(f1 log 2) * e^f.
	*/

	/* Catch overflow that might cause an endless recursion below.  */
	if( exponent(x) >= EXPONE + 12 ) {
		if (qcmp(x,maxExp)> 0) {
			if( signof(x) != 0 )
				goto underf;
			else	 goto overf;
		}
	}
	if( exponent(x) == 0 )
	{
		qmov( qone, y );
		return;
	}
	qmov(x, x2);
	qdiv( qlog2, x2, den );
	qadd( qhalf, den, den );
	qfloor( den, num );
	qifrac( num, &i, den );
	qmul( num, C1, den );
	qsub( den, x2, x2 );
	qmul( num, C2, den );
	qsub( den, x2, x2 );

	exp = x2[0].exponent;
	if (exp != 0)
		exp -= 1;
	x2[0].exponent=exp;		/* x/2				*/
	qtanh( x2, x2 );	/* tanh( x/2 )			*/
	qadd( x2, qone, num );	/* 1 + tanh			*/
	qneg( x2 );
	qadd( x2, qone, den );	/* 1 - tanh			*/
	qdiv( den, num, y );	/* (1 + tanh)/(1 - tanh)	*/

	i += exponent(y);
	sign = signof(y);
	if( i > MAXEXP )
	{
overf:
		//mtherr( "qfexp", OVERFLOW );
		qinfin(y);
		return;
	}
	if( i <= 0 ) {
underf:
		qclear(y);
		return;
	}
	y->sign = sign;
	y->exponent = i;
}
#if 0
QELT q256[NQ] = { 0,EXPONE+8,0,0x80000000,0};
extern QELT qtwo[];
void qfexp(QELT *x,QELT *y )
{
	int i=0;
	QELT r[NQ],ans[NQ],den[NQ],num[NQ],f[NQ],tmp[NQ];

	qremain(x,qlog2,r);
//printf("%130.120qe\n",*(qfloat *)r);
	qdiv(q256,r,r);
	qmov(qone,ans);
	qadd(ans,r,ans);
	qmov(r,num);
	qmov(qtwo,den);
	qmul(r,num,r);
	qmov(qtwo,f);
	for (i=0; i<30; i++) {
		qdiv(den,num,tmp);
		qadd(ans,tmp,ans);
		qmul(num,r,num);
		qadd(f,qone,f);
		qmul(f,den,den);
	}
	qmul(ans,ans,ans);
	qmul(ans,ans,ans);
	qmul(ans,ans,ans);
	qmul(ans,ans,ans);
	qmul(ans,ans,ans);
	qmul(ans,ans,ans);
	qmul(ans,ans,ans);
	qmul(ans,ans,ans);
	qmov(ans,y);
}
#endif
