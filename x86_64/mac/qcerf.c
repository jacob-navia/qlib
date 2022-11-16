/* qcerf.c
Complex error function.

2 z                      2
erf z =  -------- 1F1 ( 1/2, 3/2, -z )
sqrt(pi)
*/

#include "qhead.h"

int qchyp1f1(qcmplx *, qcmplx *, qcmplx *, qcmplx *);

int qcerf(qcmplx *zz,qcmplx *y)
{
	qcmplx a, b, z, z2;
	int neg, cj;

	qcmov(zz, &z);
	neg = 0;
	if (z.re[0].sign)
	{
		qneg(z.re);
		qneg(z.im);
		neg = 1;
	}
	cj = 0;
	if (z.im[0].sign)
	{
		qneg(z.im);
		cj = 1;
	}
	qmov(qhalf,a.re);
	qclear(a.im);

	qmov(qthree,b.re);
	b.re[0].exponent -= 1;
	qclear(b.im);

	qcmul(&z, &z, &z2);
	qneg(z2.re);
	qneg(z2.im);

	qchyp1f1(&a, &b, &z2, &b);

	qmul(oneopi, z.re, a.re);
	qmul(oneopi, z.im, a.im);
	qcmul(&a, &b, y);
	if (neg)
	{
		qneg(y->re);
		qneg(y->im);
	}
	if (cj)
		qneg(y->im);
	return 0;
}

/* qcerfw.c
Faddeeva error function.

2
w(z) =  exp(-z )  erfc(-iz)

*/


int qcerfw(qcmplx *x,qcmplx *y)
{
	qcmplx z, w;

	/* -iz */
	qmov(x->im,z.re);
	qmov(x->re,z.im);
	qneg(z.im);
	/* 1 - erf(-iz) */
	qcerf(&z, &w);
	qsub( w.re, qone, w.re);
	qneg(w.im);
	/* exp(-z^2) */
	qcmul(x, x, &z);
	qneg(z.re);
	qneg(z.im);
	qcexp(&z, y);

	qcmul(y, &w, y);
	return 0;
}
