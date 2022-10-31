/* Factorial function */

#include "qhead.h"

void qfact( Qfloatp x, Qfloatp y)
{
	Qfloat p[1], j[1];
	long i, n;
	double dn, dp, dj;

	if( x->sign != 0 )	{
		qinfin( y );
		mtherr( "qfac", DOMAIN );
		return;
	}

	dn = qtoe( x, NOROUNDING );
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

	etoq( dp, y );
	return;

fmore:

	etoq( dj, j );
	etoq( dp, p );
	qmov(qone,j);
	qmov(qone,p);
	i=1;
	do
		{
		i += 1;
		qincr( j, j );
		qmuli( j, p, p );
	}
	while( i < n );

	qmov( p, y );
}


