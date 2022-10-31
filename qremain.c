
/* Floating point remainder.
   c = remainder after dividing b by a.
   If n = integer part of b/a, rounded toward zero,
   then qremain(a,b,c) gives c = b - n * a.  */

#include "qhead.h"

int cmpm( QfloatAccump, QfloatAccump );


void qremain(Qfloatp const qa,Qfloatp const qb,Qfloatp qc )
{
	QELT *c = (QELT *)qc;
	QELT *b = (QELT *)qb;
	QELT *a = (QELT *)qa;
	QfloatAccum quot[1];
	QfloatAccum den[1], num[1];
	QELT j;
	int i;

	if( qequal(qa,qzero) ) {
		mtherr( "qremain", SING );
		qclear( (Qfloatp)c );
		return ;
	}
	den[0].mantissa[ACCUM_LENGTH-1] = 0;
	qmov( (Qfloatp)a, (Qfloatp)den );
	den[0].sign = den[0].exponent = 0;
	num[0].mantissa[ACCUM_LENGTH-1] = 0;
	qmov( (Qfloatp)b, (Qfloatp)num );
	num[0].sign = 0;

	/* Execute divide steps until num < den.
	   Least significant integer quotient bits left in quot[].  */

	memset(&quot,0,sizeof(quot));
	while( qcmp((Qfloatp)num,(Qfloatp)den) >= 0 ) {
		if( cmpm(den,num) <= 0 ) {
				subm(den, num);
				j = 1;
		}
		else {
				j = 0;
		}
		shup1(quot);
		quot[0].mantissa[ACCUM_LENGTH-1] |= j;
		shup1(num);
		num[0].exponent -= 1;
	}
	/* Normalize.  */
	while( num[0].mantissa[0] != 0 ) {
		shdn1( num );
		num[0].exponent += 1;
	}
	i = 0;
	do {
		if( num[0].mantissa[1] & SIGNBIT )
				goto normok;
		else if( (num[0].mantissa[0] & 0xffff000000000000ULL) == 0)
		{
				shiftupn(num,16);
				num[0].exponent -= 16;
				i += 16;
		}
		else if( (num[0].mantissa[0] & 0xff00) == 0 ) {
				shiftupn(num,8);
				num[0].exponent -= 8;
				i += 8;
		}
		else {
				shup1(num);
				num[0].exponent -= 1;
				i += 1;
		}
	}
	while( i<=NBITS );
	qclear((Qfloatp)num);

normok:
	/* Sign of remainder = sign of numerator. */
	num[0].sign = qb[0].sign;
	qmov( (Qfloatp)num, qc );
}
