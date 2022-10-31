#include "mconf.h"
#if 0
 21 // Return arithmetic-geometric mean (positive arguments).
 22 double agm (double a, double g)
 23 {
 24     if (a <= 0 || g <= 0) {    // domain error
 25         return 0;
 26     }
 27     for (;;) {
 28         double a0 = a;
 29         double g0 = g;
 30         a = (a0 + g0) / 2;
 31         g = sqrt(a0 * g0);
 32         if (a == g || fabs(a - g) >= fabs(a0 - g0)) {
 33             return a;
 34         }
 35     }
 36 }
#endif
void qagm(Qfloatp a,Qfloatp g,Qfloat *r)
{
	Qfloat a0[1],g0[1],wa[1],wg[1],delta[1],tmp1[1];
	int d;

	qadd(a,g,wa);
	wa[0].exponent -= 1;
	qmul(a,g,wg);
	qfsqrt(wg,wg);
	qmov(g,wg);
	while (1) {
 // 28         double a0 = a;
		qmov(wa,a0);
 // 29         double g0 = g;
//  30         a = (a0 + g0) / 2;
		qadd(a0,wg,wa);
		wa->exponent -= 1;
//   31         g = sqrt(a0 * g0);
		qmul(a0,wg,tmp1);
		qfsqrt(tmp1,wg);
// Test
//	if (a == g)
		if (wg->exponent != wa->exponent) continue;
		if (qequal(wa,wg))
			break;
		qsub(wg,wa,delta); // delta = a - g
		delta[0].sign = 0;
		qsub(g0,a0,tmp1);  // tmp1 = a0 - g0
		tmp1[0].sign = 0;
		d = qcmp(delta,tmp1); // d --> 1 if delta > tmp1, 0 if equal; -1 if delta < tmp1
		if (d >= 0)
			break;
	}
	qmov(wa,r);
}
