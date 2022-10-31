#include "qhead.h"

void qldexp(Qfloatp x,long n,Qfloatp y)
{
	long k;
	long e = exponent(x);
	long s = signof(x);

	k = e  +  n;
	qmov( x, y );
	y->sign = s;
	y->exponent = k;
	if( (k >= MAXEXP) || (n > (2 * (long)MAXEXP)) )
		qinfin(y);
	if( (k <= 0) || (n < (-2 * (long)MAXEXP)) )
		qclear(y);
}


void qfrexp(Qfloatp x,long *e,Qfloatp y)
{

	if( x[0].exponent == 0 )	{
		*e = 0;
		qclear(y);
	}
	else
	    {
		int expo = exponent(x);
		int sign = signof(x);
		*e = (long) expo - (long) EXPONE + 1;
		qmov( x, y );
		y[0].exponent = (EXPONE - 1);
		y[0].sign = sign;
	}
}
