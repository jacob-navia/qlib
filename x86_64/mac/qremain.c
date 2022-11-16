/* Floating point remainder.
c = remainder after dividing b by a.
If n = integer part of b/a, rounded toward zero,
then qremain(a,b,c) gives c = b - n * a.  */

#include "qhead.h"

void qremain(Qfloatp a,Qfloatp b,Qfloatp c)
{
	QfloatAccum den[1], num[1];
	int j;
	int i;
	QfloatAccum quot[1];

	if( qcmp(a,qzero) == 0 )
	{
		qclear( c );
		return;
	}
	qmovz(a,den);
	qmovz( b, num );
	num[0].sign = 0;

	/* Execute divide steps until num < den.
	Least significant integer quotient bits left in quot[].  */

	quot[0].mantissa[8] = 0;
	memset( quot,0,sizeof(quot) );
	while( cmpm(num,den) >= 0 )
	{
		if( cmpm(den,num) <= 0 )
		{
			__subm(den, num);
			j = 1;
		}
		else
		    {
			j = 0;
		}
		__shup1(quot);
		quot[0].mantissa[8] |= j;
		__shup1(num);
		num[0].exponent -= 1;
	}
	/* Normalize.  */
	while( num[0].mantissa[1] != 0 )
	{
		__shdn1( num );
		num[0].exponent += 1;
	}
	i = 0;
	do
	    {
		if( num[0].mantissa[2] & SIGNBIT )
			goto normok;
		else if( (num[0].mantissa[2] & 0xffff000000000000ULL) == 0 )
		{
			shiftupn(num,16);
			num[0].exponent -= 16;
			i += 16;
		}
		else if( (num[0].mantissa[2] & 0xff00000000000000ULL) == 0 )
		{
			shiftupn(num,8);
			num[0].exponent -= 8;
			i += 8;
		}
		else
		{
			__shup1(num);
			num[0].exponent -= 1;
			i += 1;
		}
	}
	while( i<=NBITS );
	memset(num,0,sizeof(num));

normok:
	/* Sign of remainder = sign of numerator. */
	num[0].sign = b[0].sign;
	__pack( num, c );
}
