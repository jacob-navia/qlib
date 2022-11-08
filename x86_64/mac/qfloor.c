
/* qfloor - largest integer not greater than x
* qround - nearest integer to x
*/

/* #include "mconf.h"  */
#include "qhead.h"

static int bmask[] = {
	0xffffffff,
	0xfffffffe,
	0xfffffffc,
	0xfffffff8,
	0xfffffff0,
	0xffffffe0,
	0xffffffc0,
	0xffffff80,
	0xffffff00,
	0xfffffe00,
	0xfffffc00,
	0xfffff800,
	0xfffff000,
	0xffffe000,
	0xffffc000,
	0xffff8000,
	0xffff0000,
	0xfffe0000,
	0xfffc0000,
	0xfff80000,
	0xfff00000,
	0xffe00000,
	0xffc00000,
	0xff800000,
	0xff000000,
	0xfe000000,
	0xfc000000,
	0xf8000000,
	0xf0000000,
	0xe0000000,
	0xc0000000,
	0x80000000,
	0x00000000
};

void qfloor(Qfloatp x,Qfloatp y)
{
	Qfloat t[1];
	unsigned long long *p;
	long e;

	if( exponent(x) == 0 ) {
		qclear( y );
		return;
	}
	qmov(x,t);
	e = (long)exponent(t) - (EXPONE-1);

	if( e <= 0 ) {
		if( signof(t) != 0 ) {
			qmov( qminusone, y );
		}
		else {
			qclear( y );
		}
		return;
	}
	if (e > NBITS) {
		qclear(y);
		return;
	}

	/* number of bits to clear out */
	e = NBITS - e;

	qmov( t, y );
	p = (unsigned long long *)(&y[0].mantissa[6]);

	while( e >= 64 ) {
		*p-- = 0;
		e -= 64;
	}
	if (e >= 32) {
		*p &= 0xFFFFFFFF00000000ULL;
		e -= 32;
		/* clear the remaining bits */
		*p &= ((unsigned long long)bmask[e] << 32);
	}
	else
		/* clear the remaining bits */
		*p &= (bmask[e&31]);

	/* truncate negatives toward minus infinity */
	if( signof(t) != 0 ) {
		if( qcmp( t, y ) != 0 )
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
		int sign = signof(z);
		setpositive(z);
		/* round to even */
		z[0].exponent -= 1;
		qfloor( z, f );
		z[0].exponent += 1;
		f[0].exponent += 1;
		if( qcmp(z,f) != 0 )
		{
rndup:
			qadd( qone, z, z );
		}
	}
	qmov( z, y );
}

