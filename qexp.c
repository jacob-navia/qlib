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

#include <stdio.h>
#include "qhead.h"
extern Qfloat InverseFact[];
#if MEASURE
unsigned char tabFreq[4096];
unsigned Max;
double tabAvg[100];
static int idx=0,idx1=0;
#endif
/*
	exp(z) = 1 + z/1! + (z^2)/2! + (z^3)/3! + (z^4)/4! + ...
*/
static void Series(Qfloatp d,Qfloatp y)
{
    Qfloat num[1],sum[1],z[1],term[1],*fac;
	int i,r;

    qmov(d,z);
	fac = &InverseFact[0];

	qadd(qone,d,sum);    // 1 + Z/1!
	qmul(d,z,num);
	fac++;

	qmov(num,term);
	term[0].exponent -= 1; // (z^2)/2!
	qadd(sum,term,sum);
	qmul(num,z,num);
	fac++;
	
    for (i=3;i<160;i++) {
		
        qmul(fac,num,term); // 1/i! * z^n --> (z^n)/i!
        r = qadd(sum,term,sum);
		if (r == 0) {
			break;
		}
        qmul(num,z,num); // z^n
		fac++;
	}   
	qmov(sum,y);
#if MEASURE
	tabFreq[idx++]=(unsigned char)i;
	if (i>Max) Max=i;
	if (idx == 4096) {
		double s=0;
		for (i=0; i<4096;i++) {
			s += tabFreq[i];
		}
		idx = 0;
		tabAvg[idx1++]=s;
	}
#endif
}


/* 363408 */
static Qfloat maxExp[1] = {
{0,0x00080013,{0xb172000000000000ULL,0,0,0,0,0,0}}};

void qfexp(Qfloatp x,Qfloatp y )
{
	Qfloat num[1], den[1], x2[1];
	long long i=0;
	int sign;

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
	if( exponent(x) == 0 ) {
		qmov( qone, y );
		return;
	}
	if (exponent(x) <= EXPONE-1) {
		qmov(x,x2);
		Series(x2,y);
		return;
	}
	//qmov(x, x2); // This was unnecessary
	// Instead of dividing by log2, multiply by the precalculated inverse
	qfma(qinv_log2,x,qhalf,den);
	qfloor( den, num );
	qifrac( num, &i, den );
	qmul( num, C1, den );
	qsub( den, x, x2 );
	qmul( num, C2, den );
	qsub( den, x2, x2 );
#if 0

	int exp = x2[0].exponent;
	if (exp != 0)
		exp -= 1;
	x2[0].exponent=exp;		/* x/2				*/
	qtanh( x2, x2 );	/* tanh( x/2 )			*/
	qincr( x2, num );	/* 1 + tanh			*/
	qneg( x2 );
	qincr( x2, den );	/* 1 - tanh			*/
	qdiv( den, num, y );	/* (1 + tanh)/(1 - tanh)	*/
#else
	Series(x2,y);
#endif
	i += exponent(y);
	sign = signof(y);
	if( i > MAXEXP ) {
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
