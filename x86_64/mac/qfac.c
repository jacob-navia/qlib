/* Factorial function */

#include "qhead.h"

void qfact( Qfloatp x, Qfloatp y)
{
	Qfloat p[1], j[1];
	long i, n;
	long double dn, dp, dj;

	if( x->sign != 0 )	{
		qinfin( y );
		mtherr( "qfac", DOMAIN );
		return;
	}

	__qtoe64( x, (unsigned short *) &dn );
	n = dn;

	if( n == 0 ) {
		qmov( qone, y );
		return;
	}

	if( n == 1 ) {
		qmov( qone, y );
		return;
	}

	if( n > 3210 ) {
		qinfin(y);
		mtherr( "qfac", OVERFLOW );
		return;
	}
	/* Cheat by using normal arithmetic */
	dp = 1.0;
	dj = 1.0;
	i = 1;
	do	{
		if( i > 17 )
			goto fmore;
		i += 1;
		dj += 1.0;
		dp *= dj;
	}
	while( i < n );

	__e64toq( (unsigned short *) &dp, y );
	return;

fmore:

	__e64toq( (unsigned short *) &dj, j );
	__e64toq( (unsigned short *) &dp, p );
	qmov(qone,j);
	qmov(qone,p);
	i=1;
	do
		{
		i += 1;
		qadd( qone, j, j );
		qmuli( j, p, p );
	}
	while( i < n );

	qmov( p, y );
}


