/* Integral of Planck's radiation formula.
*
 *                                       1
*                                ------------------
*                                 5
*                                t  (exp(1/bw) - 1)
*
 * Set
*   b = T/c2
*   u = exp(1/bw)
*
 *  In terms of polylogarithms Li_n(u)¸ the integral is
*
 *                  (           Li (u)      Li (u)                  )
*     1          4 (              3           2          log(1-u)  )
*    ----  -  6 b  ( Li (u)  -  ------  +  --------  +  ---------- )
*       4          (   4          bw              2             3  )
*    4 w           (                        2 (bw)        6 (bw)   )
*
 *   Since u > 1, the Li_n are complex valued.  This is not
* the best way to calculate the result, which is real, but it
* is adopted as a the priori formula against which other formulas
* can be verified.
*/
#include "qhead.h"
static Qfloat c1[1], c2[1], u[1];

int qplanck (Qfloatp w,Qfloatp T,Qfloatp y)
{
	static Qfloat b[1], d[1], q[1];
	qcmplx cs, cu, cp;
	double dd;

	dd = 3.7417749e-16;
	etoq (dd, c1);

	dd = 0.01438769;
	etoq (dd, c2);

	/* Temperature, degrees K.
	dd = 3000.;
	etoq (&dd, b);  */
	/* b = T / c2 */
	qdiv(c2, T, b);

	/* d = b*w */
	qmul(b, w, d);

	/* p = exp(1/d) */
	qdiv (d, qone, u);
	qfexp (u, cp.re);
	qclear(cp.im);

	/* s = polylog(4,p) */
	qcpolylog (4, &cp, &cs);

	/* s = s - polylog(3,p) / d */
	qcpolylog (3, &cp, &cu);

	qdiv (d, cu.re, cu.re);
	qdiv (d, cu.im, cu.im);
	qcsub (&cu, &cs, &cs);

	/* s = s + polylog(2,p) / (2*d*d) */

	qcpolylog (2, &cp, &cu);
	qmul(d, d, q);
	q[0].exponent += 1;
	qdiv(q, cu.re, cu.re);
	qdiv(q, cu.im, cu.im);
	qcadd(&cu, &cs, &cs);

	/* s = s + log(1-p) / (6*d*d*d) */
	qsub(cp.re, qone, cu.re);
	qclear(cu.im);
	qclog(&cu, &cu);

	qmul(d, q, q);
	qmul(qthree, q, q);
	qdiv(q, cu.re, cu.re);
	qdiv(q, cu.im, cu.im);
	qcadd(&cu, &cs, &cs);

	/* y = .25/(w*w*w*w) - 6 * b*b*b*b * s */
	qmul(b, b, d);
	qmul(d, d, d);
	qmov(d, u);
	u[0].exponent += 1;
	qmul(qthree, u, u);
	qmul(u, cs.re, cs.re);
	qmul(u, cs.im, cs.im);

	qmul(w, w, u);
	qmul(u, u, u);
	qdiv(u, qone, u);
	u[0].exponent -= 2;
	qsub( cs.re, u, cs.re);

	/* r = b*b*b*b*pi*pi*pi*pi/7.5; */
	qmul(qpi, qpi, u);
	qmul(u, u, u);
	qmul(d, u, u);
	dd = 7.5;
	etoq(dd, d);
	qdiv(d, u, u);

	/*  s = s + r */
	qadd(u, cs.re, cs.re);
	/* y=a*s */
	qmul(c1, cs.re, y);
	return 0;
}


/* Right-hand tail of Planck integral.  */
int qplanckc (Qfloatp w,Qfloatp T,Qfloatp y)
{
	Qfloat p[1];
	qplanck(w, T, p);
	qmul(c1, u, u);
	u[0].exponent -= 1;
	qsub (p, u, y);
	return 0;
}

