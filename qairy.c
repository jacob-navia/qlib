/*							qairy.c
 *
 *	Airy functions
 *
 *
 *
 * SYNOPSIS:
 *
 * int airyq( x, ai, aip, bi, bip );
 * QELT *x, *ai, *aip, *bi, *bip;
 *
 * airyq( x, ai, aip, bi, bip );
 *
 *
 *
 * DESCRIPTION:
 *
 * Solution of the differential equation
 *
 *	y"(x) = xy.
 *
 * The function returns the two independent solutions Ai, Bi
 * and their first derivatives Ai'(x), Bi'(x).
 *
 * Evaluation is by power series summation for small x,
 * by asymptotic expansion for large x.
 *
 *
 * ACCURACY:
 *
 * The asymptotic expansion is truncated at less than full working precision.
 *
 */

/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
 */

/* Airy functions */

#include "qhead.h"

/* Flag to determine which functions to compute.
 * aiconf = 0, compute all four functions.
 * aiconf = 1, do Ai, Bi only.
 * aiconf = -1, do Ai', Bi' only.
 */
int aiconf = 0;

// There were some errors in the last digits position in the hexadecimal
// constants below. I corrected them but still, the precision goes no further than
// around 105 digits
/*
 * c2 = 1/(Gamma(1/3) * 3**(1/3)) =
 *  2.5881940379280679840518356018920396347909113835493458221000181385610277E-1
  0.258819403792806798405183560189203963479091138354934582210001813856102772676790280654196405827275384313371193211789133381275035952167626014785050989848
 */
static Qfloat qc2[1] = {
	{0,EXPONE-2,{0x8483fa15b87c545dULL,0x74df94a043b3dc3cULL,0xd463e126ba75161fULL,
	0x4d8518753816a8c2ULL,0xc1c78a7d3ce1e0baULL,0xfd77ec4ffe5e9208ULL,0x9f393967f9e21a43LL}}};
/* Gamma(2/3)
 * = 1.3541179394264004169452880281545137855193272660567936983940224679637830E0
 *
 * c1 = 1/(Gamma(2/3) * 3**(2/3)) =
  0.35502805388781723926006318600418317639797917419917724058332651030081004245012671295717424605404027168842044873034949583975829267044616193710504024002
 */
static Qfloat qc1[1] = {{
	0,EXPONE-2,{0xb5c63cb138adc2f5ULL,0x2daf7609cc9edb52ULL,0xf56200d64df3d6c8ULL,
	0x42466728441708aeULL,0x08f972aaf93e7687ULL,0xa49e4b218179dc50ULL,0x67d51ee376c02db1ULL}}};
static Qfloat qthird[1]= {{ /* 1/3 */
0x00000000,0x00007fff,{0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,
0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaaULL,0xaaaaaaaaaaaaaaaa }}};
static Qfloat qsqrt3[1]= {{ // sqrt(1/3)
0x00000000,0x00008000,{0x93cd3a2c8198e269ULL,0x0c7c0f257d92be83ULL,
0x0c9d66eec69e17ddULL,0x97b58cc2cf6c8cf6ULL,0x1859454874fb1f3fULL,0x388658e5 }}};
static Qfloat q216[1]={{
0x00000000,0x00008008,{0xd800000000000000ULL
}}};
static void qainfg(Qfloatp zeta,Qfloatp uf,Qfloatp ug,Qfloatp aisin,Qfloatp aicos);
static void qanfgp(Qfloatp zeta,Qfloatp uf,Qfloatp ug,Qfloatp aisin,Qfloatp aicos);

void EXPORT airyq(Qfloatp x,Qfloatp ai,Qfloatp aip,Qfloatp bi,Qfloatp bip)
{
	double dx;
	int i, k;
	int maxf;
	Qfloat aif[1], aig[1], z[1], t[1], zeta[1], aisin[1], aicos[1];
	Qfloat uf[1], ug[1], temp[1], qden[1], qnum[1], qk[1];

	qmov( x, t );		/* t = abs(x);		*/
	t[0].sign = 0;
	qclear( aisin );
	qclear( aicos );
	qclear( zeta );
	if( t[0].exponent < 3 ) {
		qgamma( qthird, uf ); /* gamma(1/3) */
		qinv( uf, bip );
		qcbrt( qthird, z );   /* cbrt(1/3) */
		qmul( bip, z, aip );
		qneg( aip );
		qmov(qthird,temp);
		temp[0].exponent += 1;
		qgamma( temp, ug ); /* gamma(2/3) */
		qdiv( ug, z, bi );
		qmul( z, bi, ai );
		return;
	}

	qmul( t, t, z );	/* z = t * t * t		*/
	qmul( z, t, z );


	dx = qtoe( x, NOROUNDING );

	if( dx > 20.0 )
		goto asymp;
	if( dx < -20.0 )
		goto aineg;

	/* f, g */
	t[0].sign = x[0].sign;	/* put back the sign */
	z[0].sign = x[0].sign;		/* restore sign of z = x * x * x;*/

	if( aiconf < 0 )
		goto doaip;

	qmov( qone, aif );	/*f = 1.0;*/
	qmov( t, aig );		/* g = x;*/
	qmov( qone, uf );	/* uf = 1.0;*/
	qmov( t, ug );		/* ug = x;*/
	qmov( qone, qk );	/* k = 1.0;*/
	while( ((int) aif[0].exponent - (int) uf[0].exponent) < NBITS )
	{
		qmul( uf, z, uf );	/* uf *= z;*/
		qincr( qk, qk );	/* k += 1.0;*/
		qdiv( qk, uf, uf );	/* uf /=k;*/
		qmul( ug, z, ug );	/* ug *= z;*/
		qincr( qk, qk );	/* k += 1.0;*/
		qdiv( qk, ug, ug );	/* ug /=k;*/
		qdiv( qk, uf, uf );	/* uf /=k;*/
		qadd( aif, uf, aif );	/* f += uf;*/
		qincr( qk, qk );	/* k += 1.0;*/
		qdiv( qk, ug, ug );	/* ug /=k;*/
		qadd( aig, ug, aig );	/* g += ug;*/
	}
	qmul( qc1, aif, uf );		/* uf = c1 * f;*/
	qmul( qc2, aig, ug );		/* ug = c2 * g;*/
	qsub( ug, uf, ai );		/* *ai = uf - ug;*/
	qadd( ug, uf, bi );		/* *bi = sqrt3 * (uf + ug);*/
	qdiv( qsqrt3, bi, bi );
	if( aiconf > 0 )
		goto aidone;

	/* the deriviative of ai */
doaip:

	qmov( qone, qk );		/* k = 4.0;*/
	qk[0].exponent += 2;
	qmul( t, t, uf );		/* uf = x * x/2.0;*/
	uf[0].exponent -= 1;
	qmul( qthird, z, ug );		/* ug = z/3.0;*/
	qmov( uf, aif );		/* f = uf;*/
	qincr(ug, aig );		/* g = 1.0 + ug;*/
	qmul( qthird, uf, uf );		/* uf /= 3.0;*/

	while( ((int) aig[0].exponent - (int) ug[0].exponent) < NBITS )
	{
		qmul( z, uf, uf );	/* uf *= z;*/
		qdiv( qk, ug, ug );	/* ug /=k;*/
		qincr(qk, qk );	/* k += 1.0;*/
		qmul( z, ug, ug );	/* ug *= z;*/
		qdiv( qk, uf, uf );	/* uf /=k;*/
		qadd( uf, aif, aif );	/* f += uf;*/
		qincr( qk, qk );	/* k += 1.0;*/
		qdiv( qk, ug, ug );	/* ug /=k;*/
		qdiv( qk, uf, uf );	/* uf /=k;*/
		qadd( ug, aig, aig );	/* g += ug;*/
		qincr( qk, qk );	/* k += 1.0;*/
	}

	qmul( qc1, aif, uf );		/* uf = c1 * f;*/
	qmul( qc2, aig, ug );		/* ug = c2 * g;*/
	qsub( ug, uf, aip );		/* *aip = uf - ug;*/
	qadd( uf, ug, bip );		/* *bip = sqrt3 * (uf + ug);*/
	qdiv( qsqrt3, bip, bip );
	goto aidone;

	    /* negative x */

aineg:


	t[0].sign = 0;
	z[0].sign = 0;
	qfsqrt( z, zeta );
	zeta[0].exponent += 1;		/* zeta = 2.0 * sqrt(z) / 3.0	*/
	qmul( qthird, zeta, zeta );

	if( aiconf < 0 )
		goto ainegp;

	qainfg(zeta,uf,ug,aisin,aicos);
	qincr( uf, uf );

	qmov( qpi, temp ); /* zeta + pi/4 */
	temp[0].exponent -= 2;
	qadd( zeta, temp, temp );
	qfsin( temp, aisin );
	qfcos( temp, aicos );

	qmul( uf, aisin, temp );
	qmul( ug, aicos, ai );
	qsub( ai, temp, ai );

	qfsqrt( t, temp );
	qmul( qpi, temp, temp );
	qfsqrt( temp, temp );
	qdiv( temp, ai, ai );

	qmul( ug, aisin, qk );
	qmul( uf, aicos, bi );
	qadd( bi, qk, bi );
	qdiv( temp, bi, bi );

	if( aiconf > 0 )
		goto aidone;


ainegp:

	qanfgp(zeta,uf,ug,aisin,aicos);
	qincr( uf, uf );

	qmov( qpi, temp ); /* zeta + pi/4 */
	temp[0].exponent -= 2;
	qadd( zeta, temp, temp );
	qfsin( temp, aisin );
	qfcos( temp, aicos );

	qmul( uf, aicos, temp );
	qmul( ug, aisin, aip );
	qadd( aip, temp, aip );

	qfsqrt( t, temp );
	qdiv( qpi, temp, temp );
	qfsqrt( temp, temp );
	temp[0].sign = -1;
	qmul( temp, aip, aip );

	qmul( ug, aicos, qk );
	qmul( uf, aisin, bip );
	qsub( qk, bip, bip );
	temp[0].sign = 0;
	qmul( temp, bip, bip );
	goto aidone;



	/* large positive x */
asymp:

	qfsqrt( z, zeta );
	zeta[0].exponent += 1;		/* zeta = 2.0 * sqrt(z) / 3.0	*/
	qmul( qthird, zeta, zeta );

	if( aiconf < 0 )
		goto adoaip;

	maxf = MAXEXP;
	k = 1;
	qmov( qone, qk );
	qmov( qone, qden );
	qmov( qone, uf );
	qmov( qone, ug );

	do
	    {
		i = 2*k + 1;
		qmov( qk, temp );
		temp[0].exponent += 1;
		qincr( temp, temp );
		qmov( temp, qnum );

		while( i < (6*k-1) )
		{
			i += 2;
			qadd( qtwo, temp, temp );
			qmuli( temp, qnum, qnum );
		}
		qmuli( q216, qden, qden );
		qmuli( qk, qden, qden );
		qmul( zeta, qden, qden );
		qdiv( qden, qnum, temp );
		if( k & 1 )
			qsub( temp, uf, uf );
		else
			qadd( temp, uf, uf );
		qadd( temp, ug, ug );	/* for Bi */
		k += 1;
		qincr( qk, qk );
		if( temp[0].exponent < maxf )
			maxf = temp[0].exponent;
		if( temp[0].exponent > maxf )
			break;
		/*
		if( k > 40 )
		break;
		*/
	}
	while( ((int) qone[0].exponent - (int) temp[0].exponent) < NBITS/2 );

	qfsqrt( t, temp );
	qmul( qpi, temp, temp );
	qfsqrt( temp, temp );

	qfexp( zeta, qk );
	qmul( qk, temp, ai );
	qdiv( ai, uf, ai );
	ai[0].exponent -= 1;

	qdiv( temp, qk, bi );
	qmul( ug, bi, bi );

	if( aiconf > 0 )
		goto aidone;




adoaip:

	maxf = MAXEXP;
	k = 1;
	qmov( qone, qk );
	qmov( qone, qden );
	qmov( qone, uf );
	qmov( qone, ug );

	do
	    {
		i = 2*k + 1;
		qmov( qk, temp );
		temp[0].exponent += 1;
		qincr( temp, temp );
		qmov( temp, qnum );

		while( i < (6*k-3) )
		{
			i += 2;
			qadd( qtwo, temp, temp );
			qmuli( temp, qnum, qnum );
		}
		/*
		ck * -(6k+1)/(6k-1)
		*/
		qadd( qtwo, temp, temp );	/* *(6k+1) */
		qadd( qtwo, temp, temp );
		qmuli( temp, qnum, qnum );
		if( qnum[0].sign != 0 )
			qnum[0].sign = 0;
		else
			qnum[0].sign = -1;
		qmuli( q216, qden, qden );
		qmuli( qk, qden, qden );
		qmul( zeta, qden, qden );
		qdiv( qden, qnum, temp );
		if( k & 1 )
			qsub( temp, uf, uf );
		else
			qadd( temp, uf, uf );
		qadd( temp, ug, ug );	/* for Bi' */
		if( temp[0].exponent < maxf )
			maxf = temp[0].exponent;
		if( temp[0].exponent > maxf )
			break;
		k += 1;
		qincr( qk, qk );
		/*
		if( k > 40 )
		break;
		*/
	}
	while( ((int) qone[0].exponent - (int) temp[0].exponent) < NBITS/2 );

	qfsqrt( t, temp );
	qdiv( qpi, temp, temp );
	qfsqrt( temp, temp );
	qfexp( zeta, qk );
	qdiv( qk, temp, aip );
	qmul( aip, uf, aip );
	aip[0].exponent -= 1;
	if( aip[0].sign != 0 )
		aip[0].sign = 0;
	else
	    aip[0].sign = -1;
	qmul( qk, temp, bip );
	qmul( bip, ug, bip );

aidone:

	;
}


/* Auxiliary functions */

static Qfloat ai[1];
static Qfloat aip[1];
static Qfloat bi[1];
static Qfloat bip[1];
static Qfloat tt[1];
static void qainfg(Qfloatp zeta,Qfloatp uf,Qfloatp ug,Qfloatp aisin,Qfloatp aicos)
{
	Qfloat temp[1],qden[1],qnum[1],qk[1];
	QELT maxf;
	int k, i;

	maxf = MAXEXP;
	k = 1;
	qmov( qone, qk );
	qmov( qone, qden );
	/*qmov( qone, uf );*/
	qclear( uf );
	qclear( ug );
	do
	    {
		i = 2*k + 1;
		qmov( qk, temp );
		temp[0].exponent += 1;
		qincr( temp, temp );
		qmov( temp, qnum );

		while( i < (6*k-1) )
		{
			i += 2;
			qadd( qtwo, temp, temp );
			qmuli( temp, qnum, qnum );
		}
		qmuli( q216, qden, qden );
		qmuli( qk, qden, qden );
		qmul( zeta, qden, qden );
		qdiv( qden, qnum, temp );
		if( k & 1 )
		{
			if( k & 2 )
				qsub( temp, ug, ug );
			else
				qadd( temp, ug, ug );
		}
		else
		{
			if( k & 2 )
				qsub( temp, uf, uf );
			else
				qadd( temp, uf, uf );
		}
		if( temp[0].exponent < maxf )
			maxf = temp[0].exponent;
		if( temp[0].exponent > maxf )
			break;
		qincr( qk, qk );
		k += 1;
		/*
		if( k > 70 )
		break;
		*/
	}
	while( ((int) qone[0].exponent - (int) temp[0].exponent) < NBITS/2 );
}


static void qanfgp(Qfloatp zeta,Qfloatp uf,Qfloatp ug,Qfloatp aisin,Qfloatp aicos)
{
	Qfloat temp[1],qden[1],qnum[1],qk[1];
	QELT maxf;
	int k, i;

	maxf = MAXEXP;
	k = 1;
	qmov( qone, qk );
	qmov( qone, qden );
	/*qmov( qone, uf );*/
	qclear( uf );
	qclear( ug );

	do
	    {
		i = 2*k + 1;
		qmov( qk, temp );
		temp[0].exponent += 1;
		qincr( temp, temp );
		qmov( temp, qnum );

		while( i < (6*k-3) )
		{
			i += 2;
			qadd( qtwo, temp, temp );
			qmuli( temp, qnum, qnum );
		}
		/*
		ck * -(6k+1)/(6k-1)
		*/
		qadd( qtwo, temp, temp );	/* *(6k+1) */
		qadd( qtwo, temp, temp );
		qmuli( temp, qnum, qnum );
		if( qnum[0].sign != 0 )
			qnum[0].sign = 0;
		else
			qnum[0].sign = -1;
		qmuli( q216, qden, qden );
		qmuli( qk, qden, qden );
		qmul( zeta, qden, qden );
		qdiv( qden, qnum, temp );
		if( k & 1 )
		{
			if( k & 2 )
				qsub( temp, ug, ug );
			else
				qadd( temp, ug, ug );
		}
		else
		{
			if( k & 2 )
				qsub( temp, uf, uf );
			else
				qadd( temp, uf, uf );
		}
		if( temp[0].exponent < maxf )
			maxf = temp[0].exponent;
		if( temp[0].exponent > maxf )
			break;
		k += 1;
		qincr( qk, qk );
		/*
		if( k > 70 )
		break;
		*/
	}
	while( ((int) qone[0].exponent - (int) temp[0].exponent) < NBITS/2 );

}



/* Functions for polynomial approximation */

/* ai * sin(zeta) + bi * cos(zeta) */
int qainf(Qfloatp qx,Qfloatp qy,Qfloatp aisin,Qfloatp aicos)
{
	Qfloat temp[1],uf[1],ug[1],zeta[1];
	double dx;

	aiconf = 1;
	if( qx[0].exponent < 3 )
	{
		qmov( qone, qy );
		goto fdone;
	}

	dx = qtoe( qx, NOROUNDING );
	if( dx < 0.020539595906443729 )
	{
		qinv( qx, zeta );
		qainfg(zeta,uf,ug,aisin,aicos);
		qmov( uf, qy );
		goto fdone;
	}

	/*
	t = 4.0 * x * x / 9.0;
	t = 1.0/cbrt(t);
	airy( t, &ai, &aip, &bi, &bip );
	z = exp( 1.0/x );
	t = sqrt(t);
	t = sqrt(PI*t);
	y = bi * t / z;
	*/

	qmul( qx, qx, tt );
	tt[0].exponent += 2;
	qdivi( 9, tt, tt );
	qcbrt( tt, tt );

	qinv( tt, tt );
	tt[0].sign = -1;
	airyq( tt, ai, aip, bi, bip );
	tt[0].sign = 0;

	qfsqrt( tt, tt );
	qmul( tt, qpi, tt );
	qfsqrt( tt, tt );

	qmul( tt, ai, ai );
	qmul( tt, bi, bi );

	if( (aisin[0].exponent == 0) || (aicos[0].exponent == 0) )
	{
		qdiv( qx, qone, temp );
		qmov( qpi, tt );
		tt[0].exponent -= 2;
		qadd( temp, tt, tt );
		qfsin( tt, aisin );
		qfcos( tt, aicos );
	}

	qmul( aisin, ai, temp );
	qmul( aicos, bi, tt );
	qadd( tt, temp, qy );
	qsub( qone, qy, qy );

fdone:
	;
	return 0;
}

/* -ai * cos(zeta) + bi * sin( zeta) */

int qaing(Qfloatp qx,Qfloatp qy,Qfloatp aisin,Qfloatp aicos)
{
	Qfloat temp[1],zeta[1],uf[1],ug[1];
	double dx;

	aiconf = 1;
	if( qx[0].exponent < 3 ) {
		qclear(qy);
		goto fdone;
	}

	dx = qtoe( qx, DOROUNDING );
	if( dx <= 0.020539595906443729 ) {
		qinv( qx, zeta );
		qainfg(zeta,uf,ug,aisin,aicos);
		qmov( ug, qy );
		goto fdone;
	}

	qmul( qx, qx, tt );
	tt[0].exponent += 2;
	qdivi( 9, tt, tt );
	qcbrt( tt, tt );

	qinv( tt, tt );
	tt[0].sign = -1;
	airyq( tt, ai, aip, bi, bip );
	tt[0].sign = 0;

	qfsqrt( tt, tt );
	qmul( tt, qpi, tt );
	qfsqrt( tt, tt );

	qmul( tt, ai, ai );
	qmul( tt, bi, bi );

	if( (aisin[0].exponent == 0) || (aicos[0].exponent == 0) )
	{
		qinv( qx, temp );
		qmov( qpi, tt );
		tt[0].exponent -= 2;
		qadd( temp, tt, tt );
		qfsin( tt, aisin );
		qfcos( tt, aicos );
	}
	qmul( aisin, bi, temp );
	qmul( aicos, ai, tt );
	qsub( tt, temp, qy );

fdone:
	;
	return 0;
}

#if 0
/* -aip * cos(zeta) + bip * sin(zeta) */
int qaipnf(qx, qy)
QELT qx[], qy[];
{
	QELT temp[NQ];
	double dx;

	aiconf = -1;
	if( qx[1] < 3 )
	{
		qmov( qone, qy );
		goto fdone;
	}

	dx = qtoe( qx, NOROUNDING );
	if( dx <= 0.020539595906443729 )
	{
		qdiv( qx, qone, zeta );
		qanfgp();
		qmov( uf, qy );
		goto fdone;
	}

	/*
	t = 4.0 * x * x / 9.0;
	t = 1.0/cbrt(t);
	airy( t, &ai, &aip, &bi, &bip );
	z = exp( 1.0/x );
	t = sqrt(t);
	t = sqrt(PI*t);
	y = bi * t / z;
	*/

	qmul( qx, qx, tt );
	tt[1] += 2;
	qdivi( 9, tt, tt );
	qcbrt( tt, tt );

	qdiv( tt, qone, tt );
	tt[0] = -1;
	airyq( tt, ai, aip, bi, bip );
	tt[0] = 0;

	qfsqrt( tt, tt );
	qdiv( tt, qpi, tt );
	qfsqrt( tt, tt );

	qmul( tt, aip, aip );
	qmul( tt, bip, bip );

	if( (aisin[1] == 0) || (aicos[1] == 0) )
	{
		qdiv( qx, qone, temp );
		qmov( qpi, tt );
		tt[1] -= 2;
		qadd( temp, tt, tt );
		qfsin( tt, aisin );
		qfcos( tt, aicos );
	}

	qmul( aicos, aip, temp );
	qmul( aisin, bip, tt );
	qsub( temp, tt, qy );
	qsub( qone, qy, qy );

fdone:
	;
	return 0;
}
#endif



#if 0
/* -aip * sin(zeta) - bip * cos(zeta) */
int qaipng(qx, qy)
QELT qx[], qy[];
{
	QELT temp[NQ];
	double dx;

	aiconf = -1;
	if( qx[1] < 3 )
	{
		qmov( qone, qy );
		goto fdone;
	}

	dx = qtoe( qx, NOROUNDING );
	if( dx <= 0.020539595906443729 )
	{
		qdiv( qx, qone, zeta );
		qanfgp();
		qmov( ug, qy );
		goto fdone;
	}

	qmul( qx, qx, tt );
	tt[1] += 2;
	qdivi( 9, tt, tt );
	qcbrt( tt, tt );

	qdiv( tt, qone, tt );
	tt[0] = -1;
	airyq( tt, ai, aip, bi, bip );
	tt[0] = 0;

	qfsqrt( tt, tt );
	qdiv( tt, qpi, tt );
	qfsqrt( tt, tt );

	qmul( tt, aip, aip );
	qmul( tt, bip, bip );

	if( (aisin[1] == 0) || (aicos[1] == 0) )
	{
		qdiv( qx, qone, temp );
		qmov( qpi, tt );
		tt[1] -= 2;
		qadd( temp, tt, tt );
		qfsin( tt, aisin );
		qfcos( tt, aicos );
	}

	qmul( aisin, aip, temp );
	qmul( aicos, bip, tt );
	qadd( temp, tt, qy );
	qneg( qy );

fdone:
	;
	return 0;
}

#endif


/*
 *            - x**(1/4)  R(1/zeta)
 * Ai'(x)  =  ----------------------
 *             2 sqrt(pi) exp(zeta)
 */
void qaipp(Qfloatp qx,Qfloatp qy)
{
	Qfloat temp[1];

	aiconf = -1;
	if( qx[0].exponent < 10 )
	{
		qclear( qy );
		goto fdone;
	}

	/* tt = cbrt(9*zeta/4 */
	qmul( qx, qx, tt );
	tt[0].exponent += 2;
	qdiv( qnine, tt, tt );
	qcbrt( tt, tt );

	qdiv( tt, qone, tt );
	airyq( tt, ai, aip, bi, bip );

	qfsqrt( tt, tt );
	qdiv( tt, qpi, tt );
	qfsqrt( tt, tt );
	qmul( tt, aip, qy );

	qdiv( qx, qone, temp );
	qfexp( temp, temp );
	qmul( temp, qy, qy );
	qy[0].exponent += 1;
	qneg( qy );

fdone:
	;
}



#if 0

/*
 *          x**(1/4) exp(zeta)
 * Bi'(x) = ----------------- [ 1 + 1/zeta R(1/zeta) ]
 *              sqrt(pi)
 */
static int qbipp(qx, qy)
QELT qx[], qy[];
{
	QELT temp[NQ];

	aiconf = -1;
	if( qx[1] < 10 )
	{
		qclear( qy );
		goto fdone;
	}

	qmul( qx, qx, tt );
	tt[1] += 2;
	qdiv( qnine, tt, tt );
	qcbrt( tt, tt );

	qdiv( tt, qone, tt );
	airyq( tt, ai, aip, bi, bip );

	qfsqrt( tt, tt );
	qdiv( tt, qpi, tt );
	qfsqrt( tt, tt );
	qmul( tt, bip, qy );

	qdiv( qx, qone, temp );
	qfexp( temp, temp );
	qdiv( temp, qy, qy );
	qsub( qone, qy, qy );

fdone:
	;
	return 0;
}
#endif



#if 0
/*
 *                   R(1/zeta)
 * Ai(x)  =  ----------------------------
 *           2 sqrt(pi sqrt(x)) exp(zeta)
 */
static int qaiasym(qx, qy)
QELT qx[], qy[];
{
	QELT temp[NQ];

	aiconf = 1;
	if( qx[1] < 10 )
	{
		qclear( qy );
		goto fdone;
	}

	/* tt = cbrt(9*zeta/4 */
	qmul( qx, qx, tt );
	tt[1] += 2;
	qdiv( qnine, tt, tt );
	qcbrt( tt, tt );

	qdiv( tt, qone, tt );
	airyq( tt, ai, aip, bi, bip );

	qfsqrt( tt, tt );
	qmul( tt, qpi, tt );
	qfsqrt( tt, tt );
	qmul( tt, ai, qy );

	qdiv( qx, qone, temp );
	qfexp( temp, temp );
	qmul( temp, qy, qy );
	qy[1] += 1;

fdone:
	;
	return 0;
}
#endif


#if 0
/*
 *             exp(zeta)
 * Bi(x) = ----------------- [ 1 + 1/zeta R(1/zeta) ]
 *          sqrt(pi sqrt(x))
 */
int qbiasym(qx, qy)
QELT qx[], qy[];
{
	QELT temp[NQ];

	aiconf = 1;
	if( qx[1] < 10 )
	{
		qclear( qy );
		goto fdone;
	}

	qmul( qx, qx, tt );
	tt[1] += 2;
	qdiv( qnine, tt, tt );
	qcbrt( tt, tt );

	qdiv( tt, qone, tt );
	airyq( tt, ai, aip, bi, bip );

	qfsqrt( tt, tt );
	qmul( tt, qpi, tt );
	qfsqrt( tt, tt );
	qmul( tt, bi, qy );

	qdiv( qx, qone, temp );
	qfexp( temp, temp );
	qdiv( temp, qy, qy );
	qsub( qone, qy, qy );

fdone:
	;
	return 0;
}
#endif

