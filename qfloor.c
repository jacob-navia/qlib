
/* qfloor - largest integer not greater than x
* qround - nearest integer to x
*/

/* #include "mconf.h"  */
#include "qhead.h"

static QELT bmask[] = {
	0xffffffffffffffffULL,
	0xfffffffffffffffeULL,
	0xfffffffffffffffcULL,
	0xfffffffffffffff8ULL,
	0xfffffffffffffff0ULL,
	0xffffffffffffffe0ULL,
	0xffffffffffffffc0ULL,
	0xffffffffffffff80ULL,
	0xffffffffffffff00ULL,
	0xfffffffffffffe00ULL,
	0xfffffffffffffc00ULL,
	0xfffffffffffff800ULL,
	0xfffffffffffff000ULL,
	0xffffffffffffe000ULL,
	0xffffffffffffc000ULL,
	0xffffffffffff8000ULL,
	0xffffffffffff0000ULL,
	0xfffffffffffe0000ULL,
	0xfffffffffffc0000ULL,
	0xfffffffffff80000ULL,
	0xfffffffffff00000ULL,
	0xffffffffffe00000ULL,
	0xffffffffffc00000ULL,
	0xffffffffff800000ULL,
	0xffffffffff000000ULL,
	0xfffffffffe000000ULL,
	0xfffffffffc000000ULL,
	0xfffffffff8000000ULL,
	0xfffffffff0000000ULL,
	0xffffffffe0000000ULL,
	0xffffffffc0000000ULL,
	0xffffffff80000000ULL,
	0xffffffff00000000ULL
};

void qfloor(const Qfloatp x,Qfloatp y)
{
	QELT *p,tmp; 
	long e;
	int sign,modified=0;

	if( exponent(x) == 0 ) {
		qclear( y );
		return;
	}
	e = (long)exponent(x) - (EXPONE-1);
	sign = x->sign;
	if( e <= 0 ) {
		if( sign != 0 ) {
			qmov( qminusone, y );
		}
		else {
			qclear( y );
		}
		return;
	}
	if (e >= NBITS) {
		// If the number is too big, the difference between floor(x)
		// and x is almost zero. Return the number unchanged.
		if (x != y) qmov(x,y);
		return;
	}

	/* number of bits to clear out */
	e = NBITS - e;

	if (x != y) qmov(x,y);

	p = (QELT *)(&y->mantissa[MANTISSA_LENGTH-1]);

	while( e >= 64 ) {
		if (*p) {
			*p-- = 0;
			modified += 1;
		}
		else p--;
		e -= 64;
	}
	if (e >= 32) {
		tmp = *p;
		*p &= 0xFFFFFFFF00000000ULL;
		e -= 32;
		/* clear the remaining bits */
		*p &= ((QELT)bmask[e] << 32);
		modified += (*p != tmp);
	}
	else {
		/* clear the remaining bits */
		tmp = *p;
		*p &= (bmask[e&31]);
		modified += (tmp != *p);
	}

	/* truncate negatives toward minus infinity */
	if( sign && modified ) {
		qsub( qone, y, y );
	}
}





void qround(Qfloatp x,Qfloatp y)
{
	Qfloat z[1], f[1];
	int r;

	qfloor( x, z );
	qsub( z, x, f );
	r = qcmp( f, qhalf );
	if( r > 0 )
		goto rndup;
	if( r == 0 )
	{
		setpositive(z);
		/* round to even */
		z[0].exponent -= 1;
		qfloor( z, f );
		z[0].exponent += 1;
		f[0].exponent += 1;
		if( qcmp(z,f) != 0 ) {
rndup:
			qadd( qone, z, z );
		}
	}
	qmov( z, y );
}

