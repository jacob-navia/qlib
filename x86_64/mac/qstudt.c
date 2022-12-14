/*							stdtr.c
*
 *	Student's t distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qstudt( k, t, y );
* int k;
* QELT *t, *y;
*
 * qstudt( k, t, y );
*
 *
 * DESCRIPTION:
*
 * Computes the integral from minus infinity to t of the Student
* t distribution with integer k > 0 degrees of freedom:
*
 *                                      t
*                                      -
*                                     | |
*              -                      |         2   -(k+1)/2
*             | ( (k+1)/2 )           |  (     x   )
*       ----------------------        |  ( 1 + --- )        dx
*                     -               |  (      k  )
*       sqrt( k pi ) | ( k/2 )        |
*                                   | |
*                                    -
*                                   -inf.
*
* Relation to incomplete beta integral:
*
 *        1 - stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )
* where
*        z = k/(k + t**2).
*
 * For t < -2, this is the method of computation.  For higher t,
* a direct method is derived from integration by parts.
* Since the function is symmetric about t=0, the area under the
* right tail of the density is found by calling the function
* with -t instead of t.
*
* ACCURACY:
*
 */

/*							studnt.c	*/
/*	STUDNT.C	24 NOV 83

C		STUDNT.FOR		LATEST REV: 31 AUG 77
C				SLM, 31 AUG 77
C
C	EVALUTATES INTEGRAL OF STUDENT'S T DISTRIBUTION FROM
C	MINUS INFINITY TO T
C
C		USAGE:
C	CALL STUDNT(K,T,P)
C
C	K = INTEGER NUMBER OF DEGREES OF FREEDOM
C	T = RANDOM VARIABLE ARGUMENT
C	P = OUTPUT AREA
C
C	THE DENSITY FUNCTION IS
C	A*Z**-(K+2)/2,
C	WHERE Z = 1 + (T**2)/K
C	AND A = GAMMA((K+1)/2)/(GAMMA(K/2) * SQRT(K*PI)).
C	THE INTEGRAL IS EVALUATED IN CLOSED FORM BY INTEGRATION BY
C	PARTS.  THE RESULT IS EXACT, TO WITHIN ROUNDOFF ERROR.
C
C	SUBROUTINE LGAM, LOG OF GAMMA FUNCTION, IS NEEDED.
*/

/*							studnt.c 2	*/

#include "qhead.h"

void beta_distribution_invQ(Qfloatp a,Qfloatp b,Qfloatp p,Qfloatp ans);

void qstudt(int k,Qfloatp t,Qfloatp y)
{
	int kl;
	int jj;
	Qfloat x[1], rk[1], z[1], f[1], tz[1], p[1], xsqk[1], j[1];
	if( k <= 0 )
	{
		mtherr( "qstudt", DOMAIN );
		qclear(y);
		return;
	}
	/*	COMPUTE INTEGRAL FROM ZERO TO ABS(T)	*/

	kl = k;
	itoq( kl, rk );	/* degrees of freedom */

	qmov( t, x );
	x[0].sign = 0;
	if( t[0].sign != 0 && qcmp(x, qtwo) > 0 )
	{
		qmul( x, x, z );
		qadd( rk, z, z );
		qdiv( z, rk, z );
		rk[0].exponent -= 1;
		qincb( rk, qhalf, z, p );
		qmul( qhalf, p, y );
		return;
	}

	/*z = 1.0 + ( x * x )/rk;*/
	qmul( x, x, z );
	qdiv( rk, z, z );
	qadd( qone, z, z );
	/* test if k is odd or even */
	if( (k & 1) != 0)
	{

		/*	COMPUTATION FOR ODD K	*/

		qfsqrt( rk, xsqk );
		qdiv( xsqk, x, xsqk );	/*xsqk = x/sqrt(rk);*/

		qatn( xsqk, p );	/*p = arctan( xsqk );*/

		if( k > 1 )
		{
			qmov( qone, f );	/*f = 1.0;*/
			qmov( qone, tz );	/*tz = 1.0;*/
			qmov( qtwo, j );	/*j = 3;*/
			jj = 3;
			while( jj<=(k-2) )
			{
				qdiv( z, tz, tz );
				qmul( j, tz, tz );
				qadd( qone, j, j );
				qdiv( j, tz, tz );	/*tz *= (j-1)/( z * j );*/
				qadd( f, tz, f );	/*f += tz;*/
				jj += 2;
				qadd( qone, j, j );
			}
			qmul( xsqk, f, f);
			qdiv( z, f, f );
			qadd( p, f, p );	/*p += f * xsqk/z;*/
		}
		p[0].exponent += 1;
		qdiv( qpi, p, p );	/*p *= 2.0/PI;*/
	}
    /*							studnt.c 3	*/

	else
	{

		/*	COMPUTATION FOR EVEN K	*/

		qmov( qone, f );	/*f = 1.0;*/
		qmov( qone, tz );	/*tz = 1.0;*/
		qmov( qone, j );
		jj = 2;

		while( jj <= (k-2) )
		{
			qmul( j, tz, tz );
			qdiv( z, tz, tz );
			qadd( qone, j, j );
			qdiv( j, tz, tz );	/*tz *= (j - 1)/( z * j );*/
			qadd( f, tz, f );	/*f += tz;*/
			jj += 2;
			qadd( qone, j, j );
		}
		qmul( f, x, p );
		qmul( z, rk, f );
		qfsqrt( f, f );
		qdiv( f, p, p );	/*p = f * x/sqrt(z*rk);*/
	}

	/*	COMMON EXIT	*/


	if( t[0].sign != 0 )
		p[0].sign = -1;

	p[0].exponent -= 1;
	qadd( qhalf, p, p );
	qmov( p, y );	/*p = 0.5 + 0.5 * p;*/
}

/*							qstdtri
*
 *	Functional inverse of Student's t distribution
*
 *
 *
 * SYNOPSIS:
*
 * int qstdtri( k, p, t );
* int k;
* QELT *p, *t;
*
 * qstdtri( k, p, t );
*
 *
 * DESCRIPTION:
*
 * Given probability p, finds the argument t such that stdtr(k,t)
* is equal to p.
*
*/


void qstdtri(int k,Qfloatp p,Qfloatp t)
{
	int kl;
	int rflg;
	Qfloat x[1], rk[1], z[1], f[1], tz[1], xsqk[1], j[1];
	if( k <= 0 )
	{
		mtherr( "qstdtri", DOMAIN );
		qclear(t);
		return;
	}
	/* z = incbi(k/2, 1/2, 2p)  */
	kl = k;
	itoq( kl, rk );
	qmul( qhalf, rk, j );

	rflg = qcmp(p, qhalf);
	if( rflg == 0 )
	{
		qclear(t);
		return;
	}
	qmov( qhalf, f );
	f[0].exponent -= 1;
	if( qcmp(p, f) > 0 ) /* 0.25 */
	{
		qsub( f, qone, f );  /* 0.75 */
		if( qcmp(p, f) < 0 )
		{
			qmov( p, tz );
			qmul( qtwo, tz, tz );
			qsub( tz, qone, tz );
			tz[0].sign = 0;
			beta_distribution_invQ( qhalf, j, tz, z );
			qsub( z, qone, f );
			qmul(rk, z, z );
			qdiv( f, z, z );
			qfsqrt( z, z );
			if( rflg < 0 )
				qneg(z);
			qmov( z, t );
			return;
		}
	}
	if( rflg > 0 )
		qsub( p, qone, f );
	else
		qmov( p, f );
	qmul( qtwo, f, f );
	beta_distribution_invQ( j, qhalf, f, z );
	/* z = k/(k+t^2)  */
	qdiv( z, rk, tz );
	qsub( rk, tz, tz );
	qfsqrt( tz, tz );
	if( rflg < 0 )
		qneg( tz );
	qmov( tz, t );
}
