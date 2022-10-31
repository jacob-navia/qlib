/* Floating point remainder.
c = remainder after dividing b by a.
If n = integer part of b/a, rounded toward zero,
then qremain(a,b,c) gives c = b - n * a.
Integer return value contains low order bits of the integer quotient n.

According to the C99 standard,
189)   When y != 0, the remainder r = x REM y is defined regardless
of the rounding mode by the mathematical relation r = x - ny,
where n is the integer nearest the exact value of x / y;
whenever | n - x / y | = 1/2, then n is even. Thus, the remainder
is always exact. If r = 0, its sign shall be that of x.
This definition is applicable for all implementations.
*/

#include "qhead.h"


int qremquo(Qfloatp a,Qfloatp b,Qfloatp c)
{
	QELT j, qi;
	QfloatAccum den[1], num[1], quot[1];
	Qfloat t1[1],t2[1],temp[1];

	if( qcmp(a,qzero) == 0 )
	{
		mtherr( "qremquo", SING );
		qclear( c );
		return 0;
	}
	qmovz(a, den);
	setsign(den, 0);
	qmovz(b,num);
	setsign(num, 0);

	/* Execute divide steps until exponent of num < exponent of den.
	Least significant integer quotient bits left in quot[].  */

	quot[0].mantissa[8] = 0;
	memset( quot,0,sizeof(quot));

	while (num[0].exponent >= den[0].exponent)
	{
		shup1(quot);
		if( cmpm(den,num) <= 0 )
		{
			subm(den, num);
			quot[0].mantissa[8] |= 1;
		}
		shup1(num);
		decreaseExponent(num);
	}
	/* Normalize.  */
	while( num[0].mantissa[1] != 0 )
	{
		shdn1( num );
		num[0].exponent += 1;
	}
	j = 0;
	do
		{
		if( num[0].mantissa[1] & SIGNBIT )
			goto normok;
		else if( (num[0].mantissa[1] & 0xffff000000000000ULL) == 0 )
		{
			shiftupn(num,16);
			num[0].exponent -= 16;
			j += 16;
		}
		else if( (num[0].mantissa[1] & 0xff00000000000000ULL) == 0 )
		{
			shiftupn(num,8);
			num[0].exponent -= 8;
			j += 8;
		}
		else
		{
			shup1(num);
			num[0].exponent -= 1;
			j += 1;
		}
	}
	while( j<=NBITS );
	memset(num,0,sizeof(num));

normok:
	qi = quot[0].mantissa[8];
	/* Round integer quotient to nearest or even.  */

	/* Test remainder >= 0.5 * denominator  */
	pack(den,t1);
	qmov (t1, temp);
	temp[0].exponent -= 1;
	pack(num,t1);
	j = qcmp (t1, temp);
	if (j >= 0) {
		if (j == 0) {
			/* Round to even */
			if (qi & 1)
				goto qroundit;
		}
		else
		{
qroundit:
			qi += 1;
			pack(den,t1);
			pack(num,t2);
			qsub (t1, t2, t2);
			qmovz(t2,num);
		}
	}

	/* Sign of remainder = sign of numerator. */
	num[0].sign = b[0].sign;

	/* Sign of quotient.  */
	if (a[0].sign != b[0].sign)
		qi = -qi;

	pack( num, c );
	return qi;
}
