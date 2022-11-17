#include <qhead.h>

int qkolmogorov(Qfloatp const x,Qfloatp y)
{
	Qfloat t[1],p[1],r[1],t1[1];
	int res,sign;

	res=qcmp(x,qzero);
	if (res <= 0) 
		return 0;
	// Calculate x=--2.0*x*x
	qmuli(qtwo,x,t); // 2*y
	qmul(x,t,t);
	t[0].sign=-1;

	sign = 1;
	qmov(qzero,p); // p = 0.0
	qmov(qone,r); // r = 1.0

	do { // calculate t = exp(x*r*r);
		
		qmul(r,r,t1);
		qmul(x,t1,t1);
		qfexp(t1,t);
	
		// p += sign * t
		if (sign < 0)
			t[0].sign = ~t[0].sign;
		res=qadd(t,p,p);
		qadd(r,qone,r);
		sign = -sign;
	} while (res );
	p[0].exponent += 1;
	qmov(p,y);
	return 1;
}
