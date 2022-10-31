/*				qcmplx.c
* Q type complex number arithmetic
*
 * The syntax of arguments in
*
 * cfunc( a, b, c )
*
 * is
* c = b + a
* c = b - a
* c = b * a
* c = b / a.
*/
#include "qhead.h"
#include <math.h>
extern double MAXNUM;

qcmplx qczero = {0};
//qcmplx EXPORT qcone = { {{ 0,EXPONE,{0x80000000}}} ,  0  };
qcmplx qcone[1] = {
	{{ // qcone[0]
		{0,EXPONE, // qcone[0].sign, qcone[0].exponent
			{0x8000000000000000ULL,0,0,0,0,0,0}
		},
		//{0}
	}}
};

//void csqrt(cmplx *, cmplx *);

/*	c = b + a	*/

void EXPORT qcadd(qcmplx *a,qcmplx *b,qcmplx *c)
{

	qadd( &a->re[0], &b->re[0], &c->re[0] );
	qadd( &a->im[0], &b->im[0], &c->im[0] );
}


/*	c = b - a	*/

void EXPORT qcsub(qcmplx *a,qcmplx *b,qcmplx *c )
{

	qsub( &a->re[0], &b->re[0], &c->re[0] );
	qsub( &a->im[0], &b->im[0], &c->im[0] );
}

/*	c = b * a
*/
void EXPORT qcmul(qcmplx *aa, qcmplx *bb, qcmplx *cc )
{
	qcmplx a, b;
	Qfloat y[1], t[1], p[1];
	long e1, e2, e3, ear, eai, ebr, ebi;

	qcmov(aa, &a);
	qcmov(bb, &b);
	if( (eai = a.im[0].exponent) == 0 ) {/* pure real */
		qmul( a.re, b.re, cc->re );
		qmul( a.re, b.im, cc->im );
		return ;
	}
	if( (ebi = b.im[0].exponent) == 0 ) {/* pure real */
		qmul( b.re, a.re, cc->re );
		qmul( b.re, a.im, cc->im );
		return ;
	}
	if( (ear = a.re[0].exponent) == 0 )	{/* pure imaginary */
		qmul( a.im, b.re, y );
		qmul( a.im, b.im, t );
		goto imret;
	}
	if( (ebr = b.re[0].exponent) == 0 )
	{/* pure imaginary */
		qmul( b.im, a.re, y );
		qmul( b.im, a.im, t );
imret:
		if( t[0].exponent != 0 )
			qneg( t );
		qmov( t, cc->re );
		qmov( y, cc->im );
		return;
	}
	/* Overflow proofing: extract all the exponents
	* and operate with values near 1.
	*/
	a.re[0].exponent = EXPONE;
	a.im[0].exponent = EXPONE;
	b.re[0].exponent = EXPONE;
	b.im[0].exponent = EXPONE;
	qmul( b.re, a.re, y );
	qmul( b.im, a.im, t );
	e2 = ebr + ear;
	e3 = ebi + eai;
	e1 = e2 - e3;
	/* Equalize exponents in preparation for subtract. */
	if( e1 >= 0 )
	{
		if( t[0].exponent <= e1 )
			goto mreal;
		t[0].exponent -= e1;
	}
	else {
		e2 = e3;
		e1 = -e1;
		if( y[0].exponent > e1 )
			y[0].exponent -= e1;
		else
			qclear(y);
	}
	qsub( t, y, y );
mreal:
	if( y[0].exponent != 0 )
	{
		e2 = (e2 - (long)EXPONE) + ((long )y[0].exponent - (long)EXPONE);
		if( e2 > (long)MAXEXP )
		{
			mtherr( "qcmul", OVERFLOW );
			qinfin( y );
		}
		else
		    {
			if( e2 <= 0 )
				qclear(y);
			else
			    y[0].exponent = e2;
		}
	}
	qmul( b.re, a.im, t );
	qmul( b.im, a.re, p );
	e2 = ebr + eai;
	e3 = ebi + ear;
	e1 = e2 - e3;
	if( e1 >= 0 )
	{
		if( (long) p[0].exponent <= e1 )
			goto mimag;
		p[0].exponent -= e1;
	}
	else
	{
		e2 = e3;
		e1 = -e1;
		if( (long) t[0].exponent > e1 )
			t[0].exponent -= e1;
		else
			qclear(t);
	}
	qadd( p, t, p );
mimag:
	if( p[0].exponent != 0 )
	{
		e2 = (e2 - (long)EXPONE) + ((long )p[0].exponent - (long)EXPONE);
		if( e2 > (long)MAXEXP )
		{
			mtherr( "qcmul", OVERFLOW );
			qinfin( p );
		}
		else
		    {
			if( e2 <= 0 )
				qclear(p);
			else
			    p[0].exponent = e2;
		}
	}
	a.re[0].exponent = ear;
	a.im[0].exponent = eai;
	b.re[0].exponent = ebr;
	b.im[0].exponent = ebi;
	qmov( y, cc->re );
	qmov( p, cc->im );
}

/*						cmplx.c	*/


/*	c = b / a */

void EXPORT qcdiv( qcmplx *a,qcmplx * b,qcmplx * c )
{
	Qfloat y[1], p[1], t[1], s[1];
	long e1, e2, e3, en, ed, ear, eai, ebr, ebi;

	qmov( &a->re[0], s );
	qmov( &a->im[0], t );
	if( (eai = t[0].exponent) == 0 )
	{ /* pure real */
		if( s[0].exponent == 0 )
			goto overf; /* divide by zero */
		qdiv( s, &b->re[0], &c->re[0] );
		qdiv( s, &b->im[0], &c->im[0] );
		return;
	}
	if( (ear = s[0].exponent) == 0 )
	{ /* pure imaginary */
		qmov( &b->im[0], s );
		qdiv( t, &b->re[0], &c->im[0] );
		qdiv( t, s, &c->re[0] );
		qneg( &c->im[0] );
		return ;
	}
	/* Anti-overflow technique.
	* Extract all the exponents and operate with numbers near 1.
	*/
	a->re[0].exponent = EXPONE;
	a->im[0].exponent = EXPONE;
	if( (ebr = b->re[0].exponent) != 0 )
		b->re[0].exponent = EXPONE;
	if( (ebi = b->im[0].exponent) != 0 )
		b->im[0].exponent = EXPONE;

	/* y = a->re * a->re  +  a->im * a->im, with exponent ed */
	qmul( &a->re[0], &a->re[0], y );
	qmul( &a->im[0], &a->im[0], t );
	e2 = ear + ear;
	e3 = eai + eai;
	e1 = e2 - e3;
	if( e1 >= 0 )
	{
		ed = e2;
		if( (long) t[0].exponent <= e1 )
			goto mdenom;
		t[0].exponent -= e1;
	}
	else
	{
		ed = e3;
		e1 = -e1;
		if( (long) y[0].exponent > e1 )
			y[0].exponent -= e1;
		else
			qclear(y);
	}
	qadd( t, y, y );
mdenom:
	/* p = (b->re * a->re  +  b->im * a->im)/y */
	qmul( &b->re[0], &a->re[0], p );
	qmul( &b->im[0], &a->im[0], t );
	e2 = ebr + ear;
	e3 = ebi + eai;
	e1 = e2 - e3;
	if( e1 >= 0 )
	{
		en = e2;
		if( (long) t[0].exponent <= e1 )
			goto mreal;
		t[0].exponent -= e1;
	}
	else
	{
		en = e3;
		e1 = -e1;
		if( (long) p[0].exponent > e1 )
			p[0].exponent -= e1;
		else
			qclear(p);
	}
	qadd( t, p, p );
mreal:
	if( p[0].exponent != 0 ) {
		qdiv( y, p, p );
		e1 = en - ed + (long )p[0].exponent;
		if( e1 > (long)MAXEXP )
		{
			mtherr( "qcdiv", OVERFLOW );
			qinfin( p );
		}
		else
		    {
			if( e1 <= 0 )
				qclear(p);
			else
			    p[0].exponent = e1;
		}
	}
	/* c->im = (b->im * a->re  -  b->re * a->im)/y */
	qmul( &b->im[0], &a->re[0], s );
	qmul( &b->re[0], &a->im[0], t );
	e2 = ebi + ear;
	e3 = ebr + eai;
	e1 = e2 - e3;
	if( e1 >= 0 )
	{
		en = e2;
		if( (long )t[0].exponent <= e1 )
			goto mimag;
		t[0].exponent -= e1;
	}
	else
	{
		en = e3;
		e1 = -e1;
		if( (long) s[0].exponent > e1 )
			s[0].exponent -= e1;
		else
			qclear(s);
	}
	qsub( t, s, s );
mimag:
	if( s[0].exponent != 0 )
	{
		qdiv( y, s, s );
		e1 = en - ed + s[0].exponent;
		if( e1 > (long)MAXEXP )
		{
			mtherr( "qcdiv", OVERFLOW );
			qinfin( s );
		}
		else
		    {
			if( e1 <= 0 )
				qclear(s);
			else
			    s[0].exponent = e1;
		}
	}
	/* restore input exponents */
	a->re[0].exponent = ear;
	a->im[0].exponent = eai;
	b->re[0].exponent = ebr;
	b->im[0].exponent = ebi;
	qmov( s, &c->im[0] );
	qmov( p, &c->re[0] );
	return ;

overf:
	qinfin( &c->re[0] );
	qclear( &c->im[0] );
	return ;
}

/*	b = a	*/

void EXPORT qcmov( qcmplx *a, qcmplx *b )
{

	qmov( &a->re[0], &b->re[0] );
	qmov( &a->im[0], &b->im[0] );
}

void EXPORT qcneg( qcmplx *a )
{

	if( a->re[0].exponent != 0 ) /* don't produce minus 0 */
		a->re[0].sign = ~a->re[0].sign;
	if( a->im[0].exponent != 0 )
		a->im[0].sign = ~a->im[0].sign;
}


/* Absolute value of complex a, returns real y
*/
void EXPORT qcabs( qcmplx *a, Qfloatp y)
{
	Qfloat b[1], d[1];
	long ea, eb;

	ea = (unsigned int )a->re[0].exponent;
	eb = (unsigned int )a->im[0].exponent;
	if( ((ea - eb) > NBITS) || (eb == 0) )
	{
		qmov( &a->re[0], y );
		y[0].sign = 0;
		return;
	}
	if( ((eb - ea) > NBITS) || (ea == 0) )
	{
		qmov( &a->im[0], y );
		y[0].sign = 0;
		return;
	}
	/* Rescale to make geometric mean of re and im near to 1 */
	ea -= EXPONE;
	eb -= EXPONE;
	ea = (ea + eb)/2;
	qmov( &a->re[0], b );
	b[0].exponent = b[0].exponent - ea;
	qmul( b, b, b );
	qmov( &a->im[0], d );
	d[0].exponent = d[0].exponent - ea;
	qmul( d, d, d );
	/* sqrt( re**2 + im**2 ) */
	qadd( b, d, d );
	qfsqrt( d, y );
	/* restore scale factor */
	y[0].exponent = y[0].exponent + ea;
}

/* complex exponential */

void EXPORT qcexp( qcmplx *a, qcmplx *c )
{
	Qfloat r[1], t[1], u[1];

	if( a->re[0].exponent != 0 )
		qfexp( &a->re[0], r );
	else
	    qmov( qone, r );
	if( a->im[0].exponent != 0 )
	{
		qfcos( &a->im[0], t );
		qfsin( &a->im[0], u );
	}
	else
	    {
		qmov( qone, t );
		qclear( u );
	}
	qmul( r, t, &c->re[0] );
	qmul( r, u, &c->im[0] );
}


/* complex logarithm */

void EXPORT qclog( qcmplx *a, qcmplx *c )
{
	Qfloat x[1], y[1];

	qcabs( a, y );
	qflog( y, x );
#if 1 //ANSIC
	qatn2( &a->im[0], &a->re[0], y );
#else
	qatn2( &a->re[0], &a->im[0], y );
	if( qcmp(y, qpi) > 0 )
	{
		qsub( qpi, y, y );
		qsub( qpi, y, y );
	}
#endif
	qmov( x, &c->re[0] );
	qmov( y, &c->im[0] );
}

void EXPORT qcsin(qcmplx *a,qcmplx *c )
{
	Qfloat e[1], ch[1], sh[1];

	qfexp( &a->im[0], e );
	qdiv( e, qone, ch );
	qsub( ch, e, sh );
	if( sh[0].exponent > 0 )
		sh[0].exponent -= 1;
	qadd( ch, e, ch );
	if( ch[0].exponent > 0 )
		ch[0].exponent -= 1;

	qfsin( &a->re[0], e );
	qmul( e, ch, ch );

	qfcos( &a->re[0], e );
	qmul( e, sh, &c->im[0] );
	qmov( ch, &c->re[0] );
}


void EXPORT qccos(qcmplx *a,qcmplx *c)
{
	Qfloat e[1], ch[1], sh[1];

	qfexp( &a->im[0], e );
	qdiv( e, qone, ch );
	qsub( ch, e, sh );
	if( sh[0].exponent > 0 )
		sh[0].exponent -= 1;
	qadd( ch, e, ch );
	if( ch[0].exponent > 0 )
		ch[0].exponent -= 1;

	qfsin( &a->re[0], e );
	qmul( e, sh, sh );
	qneg( sh );

	qfcos( &a->re[0], e );
	qmul( e, ch, &c->re[0] );
	qmov( sh, &c->im[0] );
}

void EXPORT qcasin(qcmplx *a,qcmplx *w )
{
	qcmplx ca, ct, zz, z2;


	qmov( &a->re[0], &ca.re[0] );
	qmov( &a->im[0], &ca.im[0] );

	qmov( &ca.im[0], &ct.re[0] );	/* ct.re = -ca.im,    iz */
	qneg( &ct.re[0] );
	qmov( &ca.re[0], &ct.im[0] );	/* ct.im = ca.re */
	qcmul( &ca, &ca, &zz );	/* sqrt( 1 - z*z) */
	qsub( &zz.re[0], qone, &zz.re[0] );	/* zz.re = 1 - zz.re */
	qneg( &zz.im[0] );		/* zz.im = -zz.im */
	qcsqrt( &zz, &z2 );

	qcadd( &z2, &ct, &zz );
	qclog( &zz, &zz );
	qmov( &zz.im[0], &w->re[0] );	/* w->re = zz.im	 mult by 1/i = -i */
	qmov( &zz.re[0], &w->im[0] );	/* w->im = -zz.re */
	qneg( &w->im[0] );
}

void EXPORT qcsqrt(qcmplx *z,qcmplx *w)
{
	qcmplx q, s, y;
	double dr, dt, dx, dy;
	/* cmplx dz; */
	long ea, eb, ec;

	qcmov( z, &y );
	ea = (unsigned int )y.re[0].exponent;
	eb = (unsigned int )y.im[0].exponent;
	if( (ea == 0) && (eb == 0) )
	{ /* sqrt(0) */
		qclear( w->re );
		qclear( w->im );
		return;
	}
	/* Rescale to make max(re, im) near to 1 */
	if( eb > ea )
		ec = eb;
	else
		ec = ea;
	ec -= EXPONE;
	ec &= ~1;	/* make scale factor an even power of 2 */
	/*
	ec = ec/2;
	ec = 2*ec;
	*/
	if( ea > ec )
		y.re[0].exponent -= ec;
	else
	    qclear( &y.re[0] );
	if( eb > ec )
		y.im[0].exponent -= ec;
	else
	    qclear( &y.im[0] );
	dx = qtoe( &y.re[0], NOROUNDING );
	dy = qtoe( &y.im[0], NOROUNDING );
	/* csqrt( &dz, &dz ); */
	dr = sqrt(dx*dx + dy*dy);
	if( dx > 0 )
	{
		dt = sqrt( 0.5 * (dr + dx) );
		dr = fabs( 0.5 * dy / dt );
	}
	else
	{
		dr = sqrt( 0.5 * (dr - dx) );
		dt = fabs( 0.5 * dy / dr );
	}
	if (dy < 0)
		dr = -dr;
	etoq( dt, &q.re[0] );
	etoq( dr, &q.im[0] );
	/* Fix signs.  */
	q.re[0].sign = 0;
	q.im[0].sign = z->im[0].sign;
	/* Newton iteration */
	qcdiv( &q, &y, &s );
	qcadd( &q, &s, &q );
	if( q.re[0].exponent > 0 )
		q.re[0].exponent -= 1;
	if( q.im[0].exponent > 0 )
		q.im[0].exponent -= 1;

	qcdiv( &q, &y, &s );
	qcadd( &q, &s, &q );
	if( q.re[0].exponent > 0 )
		q.re[0].exponent -= 1;
	if( q.im[0].exponent > 0 )
		q.im[0].exponent -= 1;

#if NQ > 12
	qcdiv( &q, &y, &s );
	qcadd( &q, &s, &q );
	if( q.re[0].exponent > 0 )
		q.re[0].exponent -= 1;
	if( q.im[0].exponent > 0 )
		q.im[0].exponent -= 1;
#endif
	/* restore half the scale */
	ec >>= 1;
	if( q.re[0].exponent != 0 )
		q.re[0].exponent += ec;
	if( q.im[0].exponent != 0 )
		q.im[0].exponent += ec;
	qcmov( &q, w );
}



void EXPORT qcacos(qcmplx *z,qcmplx *w )
{
	qcmplx t;
	Qfloat p[1];

	qcasin( z, &t );
	qmov( qpi, p );
	p[0].exponent -= 1;
	qsub( &t.re[0], p, &w->re[0] );
	qmov( &t.im[0], &w->im[0] );
	qneg( &w->im[0] );
}


void EXPORT qctan(qcmplx *z,qcmplx *w )
{
	Qfloat d[1], zr[1], zi[1], t[1];

	/* d = cos( 2.0 * z->re ) + cosh( 2.0 * z->im ) */
	qmov( z->re, zr );
	zr[0].exponent += 1;
	qfcos( zr, d );
	qmov( z->im, zi );
	zi[0].exponent += 1;
	qcosh( zi, t );
	qadd( t, d, d );

	/* w->re = sin( 2.0 * z->re ) / d; */
	qfsin( zr, t );
	qdiv( d, t, w->re );
	/* w->im = sinh( 2.0 * z->im ) / d; */
	qsinh( zi, t );
	qdiv( d, t, w->im );
}

void EXPORT qccot(qcmplx *z,qcmplx *w )
{
	Qfloat d[1], zr[1], zi[1], t[1];


	/* d = cosh(2.0 * z->im) - cos(2.0 * z->re) */
	qmov( &z->im[0], zi );
	zi[0].exponent += 1;
	qcosh( zi, d );
	qmov( &z->re[0], zr );
	zr[0].exponent += 1;
	qfcos( zr, t );
	qsub( t, d, d );

	/* w->re = sin( 2.0 * z->re ) / d */
	qfsin( zr, t );
	qdiv( d, t, &w->re[0] );
	/* w->im = -sinh( 2.0 * z->im ) / d */
	qsinh( zi, t );
	qdiv( d, t, &w->im[0] );
	qneg( &w->im[0] );
}


void EXPORT qredpi(Qfloatp x,Qfloatp y)
{
	Qfloat t[1];
	double di;
	int i;

	qdiv( qpi, x, t );	/* t = x/PI */
	t[0].sign = 0;
	qadd( qhalf, t, t );	/* t += 0.5 */

	/*qifrac( t, &i, s );*/	/* i = t */
	di = qtoe( t, DOROUNDING );
	i = di;
	if( i != 0 )
	{
		itoq( i, t );		/* t = i */
		qmul( qpi, t, t );
		t[0] = x[0];
		qsub( t, x, y );
	}
	else
		qmov( x, y );
}

void EXPORT qcatan(qcmplx * z,qcmplx *w )
{
	Qfloat a[1], t[1], x[1], x2[1], y[1], y2[1];

	qmov( z->re, x );	/* x = z->re */
	qmov( z->im, y );	/* y = z->im */


	if( x[0].exponent == 0 ) /* pure imaginary */
	{
		qclear( x2 );
		qclear( w->re );
		if( qcmp(qone, y) == 0 )
		{
			qinfin( w->im );
			qneg( w->im );
			goto overf;
		}
		qneg(y);
		if( qcmp(qone, y) == 0 )
		{
			qinfin( w->im );
overf:
			mtherr( "qcatan", SING );
			return ;
		}
		qneg(y);
		goto imag;
	}
	qmul( x, x, x2 );	/* x2 = x * x */
	/* a = 1.0 - x2 - (y * y) */
	qmul( y, y, y2 );
	qsub( y2, qone, a );
	qsub( x2, a, a );

	/* t = atan2( a, 2.0 * x )/2.0 */
	qmov( x, y2 );
	y2[0].exponent += 1;
#if ANSIC
	qatn2( y2, a, t );
#else
	qatn2( a, y2, t );
#endif
	if( t[0].exponent > 0 )
		t[0].exponent -= 1;
	qredpi( t, w->re );

imag:
	if( y[0].exponent == 0 )
	{
		qclear( w->im );
		return ;
	}
	qsub( qone, y, t );	/* t = y - 1.0 */

	/* a = x2 + (t * t) */
	qmul( t, t, a );
	qadd( x2, a, a );

	qadd( qone, y, t );	/* t = y + 1.0 */
	/* a = (x2 + (t * t))/a */
	qmul( t, t, y2 );
	qadd( x2, y2, y2 );
	qdiv( a, y2, a );

	/* w->im = log(a)/4.0 */
	qflog(a, t);
	if( t[0].exponent > 1 )
		t[0].exponent -= 2;
	qmov( t, w->im );
}


/* Complex hyperbolic sine.  */

void EXPORT qcsinh (qcmplx *z,qcmplx *w)
{
	Qfloat x[1], y[1];

	qmov( z->re, x );
	qmov( z->im, y );
	/* sinh (x) * cos (y); */
	qsinh( x, w->re );
	qfcos( y, w->im );
	qmul( w->im, w->re, w->re );
	/* cosh (x) * sin (y); */
	qcosh( x, w->im );
	qfsin( y, x );
	qmul( x, w->im, w->im );
}



/* Complex inverse hyperbolic sine.  */

void EXPORT qcasinh (qcmplx *z,qcmplx *w)
{
	qcmplx u;

	qclear( u.re );
	qmov( qone, u.im );
	qcmul( z, &u, &u );
	qcasin( &u, w );
	qclear( u.re );
	qmov( qone, u.im );
	qneg( u.im );
	qcmul( &u, w, w );
}


/* Complex hyperbolic cosine.  */

void EXPORT qccosh (qcmplx *z, qcmplx *w)
{
	Qfloat x[1], y[1], u[1];

	qmov( z->re, x );
	qmov( z->im, y );

	/* cosh (x) * cos (y) */
	qcosh( x, w->re );
	qfcos( y, u );
	qmul( u, w->re, w->re);

	/* sinh (x) * sin (y) */
	qsinh( x, w->im );
	qfsin( y, u );
	qmul( u, w->im, w->im);
}


/* Complex inverse hyperbolic cosine.  */

void EXPORT qcacosh (qcmplx *z, qcmplx *w)
{
	qcmplx u;

	qcacos( z, w );
	qclear( u.re );
	qmov( qone, u.im );
	qcmul( &u, w, w );
}


/* Complex hyperbolic tangent.  */

void EXPORT qctanh (qcmplx *z, qcmplx *w)
{
	Qfloat x[1], y[1], d[1];

	qmov (z->re, x);
	x[0].exponent += 1;
	qsinh( x, w->re);
	qmov (z->im, y);
	y[0].exponent += 1;
	qfsin( y, w->im);

	/* cosh 2x + cos 2y  */
	qcosh (x, d);
	qfcos (y, x);
	qadd (x, d, d);

	qdiv (d, w->re, w->re);
	qdiv (d, w->im, w->im);
}


/* Complex inverse hyperbolic tangent.  */

void EXPORT qcatanh (qcmplx *z, qcmplx *w)
{
	qcmplx u;

	qclear( u.re );
	qmov( qone, u.im );
	qcmul (z, &u, &u);  /* i z */
	qcatan (&u, w);
	qclear( u.re );
	qmov( qone, u.im );
	qneg( u.im );
	qcmul (&u, w, w);  /* -i catan iz */
}


/* z = complex x raised to the complex y power */
void EXPORT qcpow( qcmplx *x, qcmplx *y, qcmplx *z )
{
	qcmplx w;

	if( (x->re[0].exponent == 0) && (x->im[0].exponent == 0) )
	{ /* powers of zero */
		if( y->re[0].exponent == 0 )
		{ /* real part of exponent = 0 */
			qmov( qone, &z->re[0] );
			qclear( &z->im[0] );
			if( y->im[0].exponent != 0 ) /* indeterminate angle */
				mtherr( "qcpow", DOMAIN );
			return;
		}
		if( y->re[0].sign != 0 )
		{ /* real part negative -> infinity */
			qinfin( &z->re[0] );
			qclear( &z->im[0] );
			mtherr( "qcpow", DOMAIN );
			return;
		}
		qclear( &z->re[0] ); /* 0**(+x) = 0 */
		qclear( &z->im[0] );
		return;
	}
	qclog( x, &w );
	qcmul( &w, y, &w );
	qcexp( &w, z );
}
