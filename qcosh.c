/*							qcosh.c
 *
 *	Hyperbolic cosine
 *
 *
 * SYNOPSIS:
 *
 * int qcosh(x, y);
 * QELT *x, *y;
 *
 * qcosh(x, y);
 *
 *
 * DESCRIPTION:
 *
 * cosh(x)  =  ( exp(x) + exp(-x) )/2.
 *
 * Cephes Math Library Release 2.3:  March, 1995
 * Copyright 1985, 1995 by Stephen L. Moshier
*/


#include "qhead.h"
/*cosh(x) = 1 + x2/2! + x4/4! + ... */
extern Qfloat InverseFact[];

void qcosh(Qfloatp d,Qfloatp y)
{
#if 1
	Qfloat h[1], w[1];

	qfexp( d, w );		/* w = exp(x) */
	qinv( w, h );	/* 1/exp(x) */
	qadd( w, h, y );
	decreaseExponent(y);
#else
    Qfloat num[1],sum[1],z[1],term[1],*fac;
    int i,r;

    qmul(d,d,z);
    fac = &InverseFact[0];

    qmov(qone,sum);    // 1 + Z/1!
    fac++;

    qmov(z,term);
    term[0].exponent -= 1; // (z^2)/2!
    qadd(sum,term,sum);
    qmul(z,z,num);         // z^4
    fac += 2;
    
    for (i=3;i<INVERSE_FACTORIALS/2;i++) {
    
        qmul(fac,num,term); // 1/i! * z^n --> (z^n)/i!
        r = qadd(sum,term,sum);
        if (r == 0) {
            break;
        }   
        qmul(num,z,num); // z^n
        fac += 2;
    }   
    qmov(sum,y);
#endif
}
