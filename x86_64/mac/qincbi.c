
/*							qincbi.c */
/* Inverse beta integral */

/*							qincbi()
*
 *      Inverse of imcomplete beta integral
*
 *
 *
 * SYNOPSIS:
*
 * double a, b, x, y, incbi();
*
 * x = incbi( a, b, y );
*
 *
 *
 * DESCRIPTION:
*
 * Given y, the function finds x such that
*
 *  incbet( a, b, x ) = y.
*
 * the routine performs up to 10 Newton iterations to find the
* root of incbet(a,b,x) - y = 0.
*
 *
 * ACCURACY:
*
 *                      Relative error:
* arithmetic   range      # trials      peak         rms
*    DEC       0,0.5         3400       8.8e-16     1.3e-16
*    IEEE
*
 */

#include "qhead.h"
#define BITSPREC 310

#include <stats.h>

void EXPORT beta_distribution_invQ(Qfloatp a,Qfloatp b,Qfloatp yy,Qfloatp ans)
{
	int dir;
	double da, db, dx, dy;
	Qfloat t[1], u[1], d[1], xx[1], y[1];
	Qfloat xl[1], xh[1], yl[1], yh[1];

	/* Estimate starting value by double precision routine.  */
	qtoe( a, (unsigned short *) &da );
	qtoe( b, (unsigned short *) &db );
	qtoe( yy, (unsigned short *) &dy );
	dx = beta_distribution_inv( da, db, dy );
	etoq( (unsigned short *) &dx, xx );

	/* Initialize search delta */
	dx=1e-105;
	etoq( (unsigned short *) &dx, t );
	dir = 0;
	qmov( qhalf, d );

	/* Bracket the solution.  */
	if( qcmp( xx, qone ) >= 0 )
	{
		qmov( qone, xh );
		qmov( qone, yh );
		goto findlow;
	}
	else if( qcmp( xx, qzero ) <= 0 )
	{
		qclear( xl );
		qclear( yl );
		goto findhigh;
	}
	else
	{
		qincb( a, b, xx, y );
		if( qcmp( y, yy ) == 0 )
		{
			qmov( xx, ans );
			return;
		}
		if( qcmp( y, yy ) < 0 )
		{
			qmov( y, yl );
			qmov( xx, xl );
			goto findhigh;
		}
		else
		{
			qmov( y, yh );
			qmov( xx, xh );
			goto findlow;
		}
	}

findhigh:

	while( qcmp( y, yy ) <= 0 )
	{
		t[0].exponent += 1;
		qadd( t, qone, u );
		qmul( xx, u, u );
		if( qcmp (u, qone) >= 0 )
		{
			qmov( qone, xh );
			qmov( qone, yh );
			goto iterate;
		}
		qincb( a, b, u, y );
	}
	qmov( u, xh );
	qmov( y, yh );
	goto iterate;

findlow:

	while( qcmp( y, yy ) >= 0 )
	{
		t[0].exponent += 1;
		qsub( t, qone, u );
		qmul( xx, u, u );
		if( qcmp (u, qzero) <= 0 )
		{
			qmov( qzero, xh );
			qmov( qzero, yh );
			goto iterate;
		}
		qincb( a, b, u, y );
	}
	qmov( u, xl );
	qmov( y, yl );


iterate:

	qsub( xl, xh, t );
	qmul( d, t, t );
	qadd( xl, t, xx );
	if( qcmp( xx, qone ) >= 0 )
	{
		qmov( qone, t );
		t[0].exponent -= NBITS;
		qsub( t, qone, xx );
	}
	qincb( a, b, xx, y );

	if( qcmp( y, yy ) < 0 )
	{
		/* New lower limit.  */
		qmov( xx, xl );
		qmov( y, yl );
		if( dir < 0 )
		{
			qmov( qhalf, d );
			dir = 0;
		}
		else if( dir > 1 )
		{
			/* Accelerate toward high limit.  */
			d[0].exponent -= 1;
			qadd( qhalf, d, d );
		}
		else
		{
			/* New guess at relative location of solution.  */
			qsub( y, yy, t );
			qsub( yl, yh, u );
			qdiv( u, t, d );
			dir += 1;
		}
	}
	else
	{
		/* New upper limit.  */
		qmov( xx, xh );
		qmov( y, yh );
		if( dir > 0 )
		{
			dir = 0;
			qmov( qhalf, d );
		}
		else if( dir < -1 )
		{
			/* Accelerate toward low limit.  */
			d[0].exponent -= 1;
		}
		else
		{
			/* Guess location of solution.  */
			qsub( yy, y, t );
			qsub( yl, yh, u );
			qdiv( u, t, d );
		}
		dir -= 1;
	}

	/* Exit test */
	qsub( xl, xh, t );
	qadd( xl, xh, u );
	if( ((int)u[0].exponent - (int) t[0].exponent) > BITSPREC )
		goto done;
	qsub( y, yy, t  );
	if( ((int)yy[0].exponent - (int) t[0].exponent) > BITSPREC )
		goto done;
	goto iterate;

done:
	qmov( xx, ans );
}
