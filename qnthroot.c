#include "qhead.h"

void qnthroot(Qfloatp x,Qfloatp N,Qfloatp r)
{
	double d;
	long long n = qtoll(N);
	int res;
	Qfloat z[1],y[1],nm1[1];

	d = qtoe(x,NOROUNDING);
	d = pow(d, 1.0/n);
	etoq(d,z);
	qsub(qone,N,nm1);
	/*
            1  /       X          \
	delta = -- | --------------   |
            n  \  pow(x[k-1],n-1) /
	x[k+1]= x[k]+delta
	*/

	qfpow(z,nm1,y);
	qdiv(y,x,y);
	qsub(z,y,y);
	qdivi(n,y,y);
	res=qadd(z,y,z);
	if (res == 0) goto doexit;

	qfpow(z,nm1,y);
	qdiv(y,x,y);
	qsub(z,y,y);
	qdivi(n,y,y);
	res=qadd(z,y,z);
	if (res == 0) goto doexit;

	qfpow(z,nm1,y);
	qdiv(y,x,y);
	qsub(z,y,y);
	qdivi(n,y,y);
	res=qadd(z,y,z);
	if (res == 0) goto doexit;

	qfpow(z,nm1,y);
	qdiv(y,x,y);
	qsub(z,y,y);
	qdivi(n,y,y);
	res=qadd(z,y,z);
	if (res == 0) goto doexit;
#if 0
	qfpow(z,nm1,y);
	qdiv(y,x,y);
	qsub(z,y,y);
	qdivi(n,y,y);
	qadd(z,y,r);
#endif
doexit:
	qmov(z,r);
}
