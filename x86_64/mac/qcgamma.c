/*							cgamma
*
 *	Complex gamma function
*
 *
 *
 * SYNOPSIS:
*
 * int qcgamma( x, y );
* qcmplx *x, *y;
*
 * qcgamma( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns complex-valued gamma function of the complex argument.
*
 * gamma(x) = exp (log(gamma(x)))
*
 */
/*							qclgam
*
 *	Natural logarithm of complex gamma function
*
 *
 *
 * SYNOPSIS:
*
 * int qclgam( x, y );
* qcmplx *x, *y;
*
 * qclgam( x, y );
*
 *
 *
 * DESCRIPTION:
*
 * Returns the base e (2.718...) logarithm of the complex gamma
* function of the argument.
*
 * The logarithm of the gamma function is approximated by the
* logarithmic version of Stirling's asymptotic formula.
* Arguments of real part less than +32 are increased by recurrence.
* The cosecant reflection formula is employed for arguments
* having real part less than -34.
*
Cephes Math Library Release 2.7:  March, 1998
Copyright 1984, 1998 Stephen L. Moshier
*/

/* Complex variable natural logarithm of gamma function */
#include "qhead.h"
#define qcarg(z,a) qatn2((z)->im, (z)->re, a)

extern qcmplx qcone;
int initqgam(void);
/* See qgamma.c for coefficients. */
extern Qfloat IMPORT qgamcof[NG];
extern Qfloat IMPORT qgam12[];

extern int qgamini;

void qclgam(qcmplx *x, qcmplx *y )
{
	qcmplx v, w, g, xx, t;
	Qfloat a[1], b[1];
	Qfloat *p;
	int i, cj;
	long long il;

	if( qgamini == 0 )
		initqgam();

	qmov( qone, &qcone.re[0] );
	qclear( &qcone.im[0] );
	qcmov( x, &xx );

	cj = 0;
	if (xx.im[0].sign != 0)
	{
		cj = 1;
		xx.im[0].sign = 0;
	}

	if( (xx.re[0].exponent > (QELT) (EXPONE + 8))
	    && (x->re[0].sign != 0) )
		    {
		qmov(&xx.re[0], a);
		qfloor(a,b);
		if (qcmp(a,b) == 0)
			    {
qlgover:
			mtherr("qlgam", SING);
			return;
		}
		qcmov(&qcone, &t);
		qcsub(&xx, &t, &t);
		qclgam( &t, &t );  /* ln gamma(1-z)  */
		qneg(b);
		qifrac(b,&il,a);
		/*
		if (il & 1)
		il += 1;
		*/
		lltoq(&il, b);
		qadd (b, &xx.re[0], &xx.re[0]);
		qmul( &xx.re[0], qpi, &g.re[0] );	/* PI/(sin(PI*z))	*/
		qmul( &xx.im[0], qpi, &g.im[0] );
		qcsin( &g, &g );
		if(g.re[0].exponent == 0 && g.im[0].exponent == 0)
			    goto qlgover;
		qclog(&g, &g);
		qflog( qpi, &v.re[0] );
		qclear( &v.im[0] );
		qcsub( &g, &v, y );
		qcsub( &t, y, y );	/* ... /gamma(x)	*/
		qmul(qpi, b, b);
		qsub(b, &y->im[0], &y->im[0]);
		goto qcldone;
	}

		/* range reduction: transform argument to be greater than 32.
	To satisfy Im {clgam(z)} = arg cgamma(z), accumulate
	arg v during the recurrence.  */
		/*qcmov( x, &xx );*/
		qclear( a );
	qclear( b );
	qcmov( &qcone, &v );
#if 1 == 28
	while( xx.re[1] <= (QELT) (EXPONE + 8) )
	#else
	    while( xx.re[0].exponent <= (QELT) (EXPONE + 5) )
	#endif
	    {
		qcmul( &xx, &v, &v );
		qcarg( &xx, b);
		qadd( b, a, a );
		qcadd( &qcone, &xx, &xx );
	}
	qcabs(&v, b);
	qflog(b,&v.re[0]);
	qmov(a,&v.im[0]);
	qcneg(&v);

		/*  Asymptotic series in 1/x**2 */
		qcmul( &xx, &xx, &w );
	qcdiv( &w, &qcone, &w );

		p = &qgamcof[0];
	qmul( &w.re[0], p, &g.re[0]);
	qmul( &w.im[0], p, &g.im[0] );
	p += 1;
	qadd( &g.re[0], p, &g.re[0] );
	for( i=0; i<NG-2; i++ )
		{
		qcmul( &w, &g, &g );
		p += 1;
		qadd( &g.re[0], p, &g.re[0] );
	}

		qcdiv( &xx, &g, &g );

		/* g + (x - 0.5)*log(x) - x + qgam12	*/
		qsub( qhalf, &xx.re[0], &t.re[0] );
	qmov( &xx.im[0], &t.im[0] );
	qclog( &xx, &w );
	qcmul( &t, &w, &t );
	qcsub( &xx, &t, &t );
	qadd( qgam12, &t.re[0], &t.re[0] );
	qcadd( &g, &t, &g );

		qcadd( &v, &g, y );
qcldone:
	if (cj)
		    {
		y->im[0].sign = ~(y->im[0].sign);
	}
}

/*							qgamma()	*/

/* Complex variable gamma function check routine */


void qcgamma(qcmplx *x,qcmplx *y )
{
	qcmplx xx;

	qcmov( x, &xx );
	qclgam( &xx, y );
	qcexp( y, y );
}

