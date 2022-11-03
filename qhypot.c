#include "qhead.h"

void qhypot(Qfloatp xx, Qfloatp yy, Qfloatp z)
{
	Qfloatp big,small;
	Qfloat x[1],y[1];
	int s;

	qmov(xx,x);
	qmov(yy,y);
	x[0].sign=0;
	y[0].sign = 0;
	s = qcmp(x,y);
	if (s > 0) {
		big = x;
		small = y;
	}
	else if (s < 0) {
		big = y;
		small = x;
	}
	else {
		qmul(qsqrt2,x,z);
		return;
	}
	qdiv(big,small,z);
	qmul(z,z,z);
	qadd(qone,z,z);
	qfsqrt(z,z);
	qmul(big,z,z);
}
