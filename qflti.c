/*							qflt.c
 *			QFLOAT
 *
 *	Extended precision floating point routines
 *
 *	asctoq( string, q )	ascii string to q type
 *	eto( &d, q )		IEEE double precision to q type
 *	e24toq( &d, q )		IEEE single precision to q type
 *	qabs(q)			absolute value
 *	qadd( a, b, c )		c = b + a
 *	qclear(q)		q = 0
 *	qcmp( a, b )		compare a to b
 *	qdiv( a, b, c )		c = b / a
 *	qifrac( x, &l, frac )   x to integer part l and q type fraction
 *	qfrexp( x, l, y )	find exponent l and fraction y between .5 and 1
 *	qldexp( x, l, y )	multiply x by 2^l
 *	qinfin( x )		    set x to infinity, leaving its sign alone
 *	qmov( a, b )		b = a
 *	qmul( a, b, c )		c = b * a
 *	qmuli( a, b, c )	c = b * a, a has only 16 significant bits
 *	qisneg(q)		    returns sign of q
 *	qneg(q)			    q = -q
 *	qnrmlz(q)		    adjust exponent and mantissa
 *	qsub( a, b, c )		c = b - a
 *	qtoasc( a, s, n )	q to ASCII string, n digits after decimal
 *	qtod( q, &d )		convert q type to DEC double precision
 *	double qtoe( q, doround )		convert q type to IEEE double precision
 *	qtoe24( q, &d )		convert q type to IEEE single precision
 *
 *   qinv(src,result)    Inverse: result = 1/src. Jacob, 2018
 *
 * Data structure of the number (a "word" is 16 bits)
 *
 *	sign word		(0 for positive, -1 for negative)
 *	exponent		(EXPONE for 1.0)
 *	high guard word		(always zero after normalization)
 *	N-1 mantissa words	(most significant word first,
 *				 most significant bit is set)
 *
 * Numbers are stored in C language as arrays.  All routines
 * use pointers to the arrays as arguments.
 *
 * The result is always normalized after each arithmetic operation.
 * All arithmetic results are chopped. No rounding is performed except
 * on conversion to double precision.
 */

/*
 * Revision history:
 *
 * SLM,  5 Jan 84	PDP-11 assembly language version
 * SLM,  2 Mar 86	fixed bug in asctoq()
 * SLM,  6 Dec 86	C language version
 * JN    1995-2011       Many improvements, modifications, etc
 * JN   2017 Ported to ARM64
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "qhead.h"
void mdnorm(QfloatAccump x);
#define NTEN 17
#define MINNTEN -65536*2
#define MAXNTEN  65536*2
#if 1
#pragma pack(push,1)
static Qfloat qtens[NTEN+1] = {
/* 10**262144 */
//{0, EXPONE+0xD49A8, {0xb75b3786baad6104ULL, 0xe01ab22ec1c080a6ULL, 0xf662839762ea5dc8ULL, 
//0xd391e3b3292701d8ULL, 0x2ec9e02d62d75239ULL, 0xa5943f6f1d015a17ULL, 0x9143c3d1af4579caULL}},
/* 10 ** 131072 */
{0, 0x000ea4d4,{0xd8a7940bfe7769bdULL,0x084da5a21ee8580dULL,0x1cd54215cfa3bc72ULL,
0x741442da159de9aaULL,0xbd165810b2181f41ULL,0x4f411a904733b27cULL,
0x32d368da6306f6a2ULL}}, // rounded since next digit is 9410500921eeea19
/* 65536 */
{0,EXPONE+0x35269,{0xeb81cf19f0160e7fULL,0x71cf55536ae2e412ULL,0x4ba19289c718ebe1ULL,
0x0c03d8b3134c5b9eULL,0x50d4a3d75514e6a8ULL,0x3113f680b5b060fbULL,
0x0b783161dd0b904dULL}},
/* 32768  */
{0,EXPONE+0x1a934,{0xf58a326a5783a749ULL,0xb3c1068fdbaa558bULL,0xa3a3210f13afdf93ULL,
0x00ba376f07afc6c5ULL,0x26ff6a4c606cf3f9ULL,0x962356cdefe0ea0fULL,
0x3dd8ece4f5b52fd5ULL}},
/* 16384 */
{0,EXPONE+0xd49a,{0xb1485471f16603b5ULL,0x6226e1134dc6c888ULL,0xeea4a43b2888144aULL,
0x7c4f561e933c14d4ULL,0xa23cb61d4e761864ULL,0x2b5b0598d235b819ULL,
0x9984733ae85f4836ULL}},
/* 8192 */
{0,EXPONE+0x6a4d,{0x96a3a1d17faf211aULL,0x0c7c2892305f4e12ULL,0x072b211aceb5055eULL,
0x412e0901a7bb6c7dULL,0x98ab274a132284c8ULL,0x0c3f77471b93102aULL,
0xd54759925092b445ULL}},
/* 4096 */
{0,EXPONE+0x3526,{0xc46052028a20979aULL,0xc94c153f804a4a92ULL,0x65761fb2444e2267ULL,
0xdd5cf7c945f22a3fULL,0xfff0d1d2e18469b1ULL,0xea38aac7e917a569ULL,
0x2db0d164ea315c84ULL}},
/* 2048 */
{0,EXPONE+0x1a93,{0x9e8b3b5dc53d5de4ULL,0xa74d28ce329ace52ULL,0x6a3197bbebe3034fULL,
0x77154ce2bcba1964ULL,0x8b21c11eb962b1b6ULL,0x1b93cf2ee5ca6f7eULL,
0x928e61d08e2d6942ULL}},
/* 1024 */
{0,EXPONE+0x0d49,{0xc976758681750c17ULL,0x650d3d28f18b50ceULL,0x526b988275249b0fULL,
0xd6f4b6d27bd1c61cULL,0x253069a5c329f9afULL,0x4a9cf154d24daa37ULL,
0x8187a351927d3cc8ULL}},
/* 512 */
{0,EXPONE+0x06a4,{0xe319a0aea60e91c6ULL,0xcc655c54bc5058f8ULL,0x9c6583981d134cbaULL,
0x422d38ea3584cde4ULL,0x0b9a1d7d634df2d8ULL,0x74a24bbae09b3399ULL,
0x549d5d6f25948477ULL}},
/* 256 */
{0,EXPONE+0x0352,{0xaa7eebfb9df9de8dULL,0xddbb901b98feeab7ULL,0x851e4cbf3de2f98aULL,
0xae780c7fea81c788ULL,0x5a6b43a2a6c495b8ULL,0xccd604d64e2ddab2ULL,
0xbb01f9e94dce0d7bULL}},
/* 128 */
{0,EXPONE+0x01a9,{0x93ba47c980e98cdfULL,0xc66f336c36b10137ULL,0x0234f3fd7b08dd39ULL,
0x0bc3c54e3f40f7e6ULL,0x424ba54f80400000ULL,0x0000000000000000ULL,
0x0000000000000000ULL}},
/* 64 */
{0,EXPONE+0x00d4,{0xc2781f49ffcfa6d5ULL,0x3cbf6b71c76b25fbULL,0x50f8080000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL}},
/* 32 */
{0,EXPONE+0x006a,{0x9dc5ada82b70b59dULL,0xf020000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL}},
/* 16 */
{0,EXPONE+0x0035,{0x8e1bc9bf04000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL}},
/* 8 */
{0,EXPONE+0x001a,{0xbebc200000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL}},
/* 4 */
{0,EXPONE+0x000d,{0x9c40000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL}},
/* 2 */
{0,EXPONE+0x0006,{0xc800000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL}},
/* 1 */
{0,EXPONE+0x0003,{0xa000000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL}},
};

static Qfloat qmtens[NTEN+1] = {
/* -262144 */
//{0,EXPONE-0xd49a8,{0xb2b64ef19d3f7af4ULL, 0xc34ec2f101cf3b1eULL, 0x4b17aa38482448caULL, 
//0xee9191606531f7d5ULL, 0xbdc54d6b86351fd6ULL, 0x3d261d335adec36eULL, 0x6c131d936cdca200}},
/* -131072 */
{0,0x00015b2d,{0x973ecee321687655ULL,0x586003998e195e6dULL,0xfb30d3b4168e9301ULL,
0xfec838f9b5d01eacULL,0x6c03ed24bfa41e54ULL,0x2090f8357dc75bc3ULL,
0x3219388138d6228fULL}},
/* -65536 */
{0,0x0004ad97,{0x8b2358ebdc628061ULL,0x2591506b30320e73ULL,0x99f6b65165d787bfULL,
0x0947fbefe1b7e2e9ULL,0xcaac64f8d86be15aULL,0x1509a10228680295ULL,
0x5d8179ee450f3086ULL}},
/* -32768 */
{0,EXPONE-108853,{0x8573f08e5ca085e5ULL,0x718e8c4ea73abf7cULL,0xd2b17e716f5893a0ULL,
0xd31f7cb6ffc3e1d3ULL,0x74a104cf5e82cd5bULL,0x197a57c855637174ULL,
0xb64206fe2d649b78ULL}},
/* -16384 */
{0,EXPONE-54427,{0xb8d5bbe70e108517ULL,0x456e9e094283bb25ULL,0x538a30f4bfe56662ULL,
0xed3491e867b318b7ULL,0x1a29f86c2138b9fdULL,0x1876005af6efa18bULL,
0xffe9c27e143dc2aaULL}},
/* -8192 */
{0,EXPONE-27214,{0xd986c20b686da869ULL,0x5d1d4fd85b05f4c2ULL,0xeef183849aee88f7ULL,
0x6e7ad89a55416d56ULL,0x354824302cb9572cULL,0x4628bee7c8b0d71aULL,
0xdc43411943b307d6ULL}},
/* -4096 */
{0,EXPONE-13607,{0xa6dd04c8d2ce9fdeULL,0x2de38123a1c3cffcULL,0x20305d0244e091baULL,
0x5e2d7403972f6f2bULL,0x04e663ac16d51215ULL,0x9834265791e1cf28ULL,
0xb1e18621d364461bULL}},
/* -2048 */
{0,EXPONE-6804,{0xceae534f34362de4ULL,0x492512d4f2ead2cbULL,0x8263ca5cbc774bd9ULL,
0x71aad59046c74249ULL,0xddddc0dfbb0d6f2cULL,0x0583c8fe05b57dc3ULL,
0x2cba0316ecf72da6ULL}},
/* -1024 */
{0,EXPONE-3402,{0xa2a682a5da57c0bdULL,0x87a601586bd3f698ULL,0xf53e94d1b2357c32ULL,
0xc0eaff3755a2ddcdULL,0x5191a70f9e339aa9ULL,0x62702807d10d8cfaULL,
0x791d963cf9e43fbcULL}},
/* -512 */
{0,EXPONE-1701,{0x9049ee32db23d21cULL,0x7132d332e3f204d4ULL,0xe7317d62209b6a93ULL,
0xd4c94a9da0693e0cULL,0xfefde89b01e51068ULL,0xe0e47572dfadea99ULL,
0x47d7c1375d35bafcULL}},
/* -256 */
{0,EXPONE-851,{0xc0314325637a1939ULL,0xfa911155fefb5308ULL,0xa23e2ed27766e8ccULL,
0x9b03537708b1648fULL,0x7a730e52e6e74360ULL,0x55f3dff45e1bdf17ULL,
0x996c8440ac86a8a5ULL}},
/* -128 */
{0,EXPONE-426,{0xddd0467c64bce4a0ULL,0xac7cb3f6d05ddbdeULL,0xe26ca6063461fffaULL,
0x4ed775fc49f27952ULL,0xf2f3a45f9009f3c9ULL,0xee49e1c12e611caaULL,
0x30eacfefed54ddebULL}},
/* -64 */
{0,EXPONE-213,{0xa87fea27a539e9a5ULL,0x3f2398d747b36224ULL,0x2a1fee40d90aab31ULL,
0x0e128b5d938cfb3fULL,0x6f203147025d1129ULL,0xe4926e2d6669a035ULL,
0xeeb11beaf4823361ULL}},
/* -32 */
{0,EXPONE-107,{0xcfb11ead453994baULL,0x67de18eda5814af2ULL,0x0b5b1aa028ccd99eULL,
0x59e338e387ad8e28ULL,0x56760892162196afULL,0x3c7ede8038a3b843ULL,
0xd185e627c28efd15ULL}},
/* -16 */
{0,EXPONE-54,{0xe69594bec44de15bULL,0x4c2ebe687989a9b3ULL,0xbf716c1add27f085ULL,
0x23ccd3484db670aaULL,0xd6e012cf7fa7cf42ULL,0x858294ca267ae62dULL,
0xb82b4f8baffdb8f7ULL}},
/* -8 */
{0,EXPONE-27,{0xabcc77118461cefcULL,0xfdc20d2b36ba7c3dULL,0x3d4d3d758161697cULL,
0x7068f3b46d2f8350ULL,0x5705755fd37a3b04ULL,0x8dd947039f9cbde6ULL,
0x55091bdd8125039bULL}},
/* -4 */
{0,EXPONE-14,{0xd1b71758e219652bULL,0xd3c36113404ea4a8ULL,0xc154c985f06f6944ULL,
0x67381d7dbf487fcbULL,0x923a29c779a6b50bULL,0x0f27bb2fec56d5cfULL,
0xaacd9e83e425aee6ULL}},
/* -2 */
{0,EXPONE-7,{0xa3d70a3d70a3d70aULL,0x3d70a3d70a3d70a3ULL,0xd70a3d70a3d70a3dULL,
0x70a3d70a3d70a3d7ULL,0x0a3d70a3d70a3d70ULL,0xa3d70a3d70a3d70aULL,
0x3d70a3d70a3d70a4ULL}},
/* -1 */
{0,EXPONE-4,{0xccccccccccccccccULL,0xccccccccccccccccULL,0xccccccccccccccccULL,
0xccccccccccccccccULL,0xccccccccccccccccULL,0xccccccccccccccccULL,
0xcccccccccccccccdULL}},
};
#pragma pack(pop)
#else
#include "qtens-bigendian.h"
#endif
int normlz(QfloatAccump x,int *SC);
int shift(QfloatAccump x,int SC);
void inverse(QfloatAccump,Qfloatp);
void square(QfloatAccump a,QfloatAccump b,int prec );

extern void pack(QfloatAccump,Qfloatp);
#define qneg(z) ((z)[0].sign ^= -1)

/* Test if negative  */

int qisneg(Qfloatp x)
{
	return (x->sign != 0);
}

void qifrac(Qfloatp x,long long *i,Qfloatp frac)
{
	QfloatAccum ac1[1];
	int SC;

	SC = x->exponent - (EXPONE - 1);
	if( SC <= 0 )	{
		/* if exponent <= 0, integer = 0 and argument is fraction */
		*i = 0L;
		qmov(x,frac);
		return;
	}
	qmovz( x, ac1);
	if( SC > 63 )	{
		/*
		;	long integer overflow: output large integer
		;	and correct fraction
		*/
		*i = 0x7fffffffffffffffULL;
		if (SC < sizeof(QELT) * 8)
			shift( ac1 , SC);
		else {
			// Total failure
			qclear(frac);
			return;
		}
	}
	else {
		shift( ac1 , SC );
		*i = ac1[0].mantissa[0];
	}
	if( signof(x) )
		*i = -(*i);
	ac1[0].sign = 0;
	ac1[0].exponent = EXPONE - 1;
	ac1[0].mantissa[0] = 0;
	if( normlz(ac1,&SC) ) {
		qclear( frac );
		return;
	}
	ac1[0].exponent -= SC;
	pack( ac1, frac );
}

/*
;	subtract
;
;	QELT a[NQ], b[NQ], c[NQ];
;	qsub( a, b, c );	 c = b - a
 */

#undef qadd
#undef qsub

int qsub(const Qfloatp a,const Qfloatp b, Qfloatp c )
{

	return qadd_subtract( a, b, c ,1);
}


/*
;	add
;
;	QELT a[NQ], b[NQ], c[NQ];
;	qadd( a, b, c );	 c = b + a
Basic Floating Point Addition Algorithm
---------------------------------------
Assuming that the operands are already in the right format, performing
floating point addition: Result = X + Y = (Xm x 2Xe) + (Ym x 2Ye)
involves the following steps:
(1) Align binary point:
 Initial result exponent: the larger of Xe, Ye
 Compute exponent difference: Ye - Xe
 If Ye > Xe Right shift Xm that many positions to form Xm 2 Xe-Ye
 If Xe > Ye Right shift Ym that many positions to form Ym 2 Ye-Xe
(2) Compute sum of aligned mantissas:
i.e Xm2 Xe-Ye + Ym or Xm + Xm2 Ye-Xe
(3) If normalization of result is needed, then a normalization step follows:
 Left shift result, decrement result exponent (e.g., if result is 0.001xx) or
 Right shift result, increment result exponent (e.g., if result is 10.1xx)
Continue until MSB of data is 1
(4) Check result exponent:
 If larger than maximum exponent allowed return exponent overflow
 If smaller than minimum exponent allowed return exponent underflow
(5) If result mantissa is 0, may need to set the exponent to zero by a special step
to return a proper zero.

Returns:
0 No changes. One of the numbers is too small
1 Normal return
-1 Underflow
-2 Overflow
 */

int  qadd(const Qfloatp a, const Qfloatp b,Qfloatp c )
{

	return qadd_subtract( a, b, c ,0);
}

#ifndef QADD_SUBTRACT_IN_ASM 
int  qadd_subtract(const Qfloatp a,const Qfloatp b,Qfloatp c ,int subflg)
{
	long lt;
	long long i;
	int lost,SC;
	QfloatAccum ac1[1];
	QfloatAccum ac2[1];
	QfloatAccump plarger,psmaller;

	/* compare exponents */
	lt = (long) a->exponent - (long) b->exponent;
	if( lt > 0 )	{	/* put the larger number in ac2 */
		qmovz( a, ac2);
		qmovz( b, ac1);
		if (subflg)
			ac2[0].sign = ~ac2[0].sign;
		lt = -lt;
	}
	else  {
		qmovz(a,ac1);
		qmovz(b,ac2);
		if (subflg)
			ac1[0].sign = ~ac1[0].sign;
	}
	plarger = ac2;
	psmaller = ac1;
	lost = 0;
	// At this point lt is always zero or less than zero
	if( lt != 0 )	{
		if( lt < -NBITS-1 )	{	/* answer same as larger addend */
			pack(plarger,c);
			return 0;
		}

		lost = shift( psmaller ,lt); /* shift the smaller number down */
	}
	else { /* exponents were the same, so must compare mantissae */
		i = cmpm( psmaller, plarger );
		if( i == 0 )	{	/* the numbers are identical */
			/* if different signs, result is zero */
			if( psmaller->sign != plarger->sign )
				goto underf;
			/* if exponents zero, result is zero */
			if( plarger->exponent == 0 )
				goto underf;
			/* if same sign, result is double */
			if( psmaller->exponent >= MAXEXP )			{
				qclear(c);
				if( plarger->sign != 0 )
					qneg(c);
				goto overf;
			}
			qmov(a,c);
			c->exponent += 1;
			return 1;
		}
		if( i > 0 )		{	/* put the larger number in ac2 */
			QfloatAccump tmpP = psmaller;
			psmaller = plarger;
			plarger = tmpP;
		}
	}

	if( psmaller->sign == plarger->sign )	{
		addm( psmaller, plarger );
		subflg = 0;
	}
	else {
		subm( psmaller, plarger );
		subflg = 1;
	}
	if (normlz(plarger,&SC))
		goto underf;
	lt = (long long)plarger->exponent - SC;
	if (lt > MAXEXP)
		goto overf;
	if (lt < 0)
		goto underf;
	plarger->exponent = lt;
	if (roundAccum(plarger,lost,subflg)) goto overf;
	pack( plarger, c );
	return 1;

underf:
	qclear(c);
	return -1;

overf:
	qinfin(c);
	return -2;
}
#endif
void  qdivi(long long a,const Qfloatp b,Qfloatp c )
{
	long lt;
	int  signOfA,toshift,exponentOfA;
	QfloatAccum ac3[1];
	QELT d;

	if (exponent(b) == 0 ) goto divunderf;
	if (a == 0) goto divoverflow;	
	/* Avoid exponent underflow in mdnorm.  */
	lt = b->exponent;
	if (a > 0)
		signOfA = 0;
	else {
		signOfA = -1;
		a = -a;
	}
	ac3[0].exponent = 4;
	toshift = bsr64(a);
	d = (unsigned long long)a << toshift;
	exponentOfA = EXPONE - 1 + (64-toshift);
	divi(d,b,ac3);
	mdnorm(ac3);
	/* calculate sign */
	if( signOfA == signof(b) )
		ac3[0].sign = 0;
	else
		ac3[0].sign = -1;

	/* calculate exponent */
	lt = lt + (long long )ac3[0].exponent -4L - exponentOfA;
	lt = lt + (long long)EXPONE + 1;
	if( lt >= MAXEXP ) goto divoverflow;
	else if( lt <= 0 ) goto divunderf;
	else ac3[0].exponent = lt;
	pack(ac3,c);
	return;
divunderf:
	qclear(c); 
	return;
divoverflow:
	qinfin(c);

}

#if 0
void qmulAccum(const Qfloatp a,const QfloatAccump ac2,Qfloatp c)
{
	QfloatAccum ac3[1];
	long long lt;

	if (exponent(a) == 0) goto underf;
	qmovz( a,ac3);
	mulm( ac2,a,ac3);
	mdnorm(ac3);
    lt = ((long long)exponent(ac3) - (EXPONE-1)) + ((long long)ac2->exponent - (EXPONE - 1));
	/* calculate sign of product */
	ac3[0].sign = ( signof(a) == signof(ac2) ) ? 0 : -1;
	lt = lt + EXPONE - 1;
	if( lt >= (long long)MAXEXP )	goto overf;
	if( lt <= 0 )	goto underf;
	ac3[0].exponent = lt;
	pack( ac3, c );
	return;
underf:
	qclear(c);
	return ;
overf:
	qinfin(c);
}
#endif

int  qfma(const Qfloatp a,const Qfloatp b,const Qfloatp addend,Qfloatp c )
{
	long long lt,i;
	QfloatAccum ac3[1];
	QfloatAccum ac2[1];
	QfloatAccum ac1[1];
	QfloatAccump plarger,psmaller,tmpP;
	int lost=0,SC,subflg=0;

	// Avoid jumps for the normal case.
	if (exponent(a) == 0 )	goto underf;
	if (exponent(b) == 0) goto underf;
    /* detect multiplication by small integer a */
    if ((a->mantissa[1] == 0) && (a->mantissa[2] == 0) &&
		(a->mantissa[3] == 0) && (a->mantissa[4] == 0) &&
		(a->mantissa[5] == 0) && (a->mantissa[6] == 0)) {
			mulin( a->mantissa[0],b, ac3 );
			lt = ((long long)a->exponent - (EXPONE-1)) + ((long long)ac3[0].exponent - (EXPONE - 1)); 
    }    
    else /* detect multiplication by small integer b */
        if ((b->mantissa[1] == 0) && (b->mantissa[2] == 0) &&
			(b->mantissa[3] == 0) && (b->mantissa[4] == 0) &&
			(b->mantissa[5] == 0) && (b->mantissa[6] == 0)) {
                mulin(  b->mantissa[0], a, ac3 );
                lt = ((long)exponent(b) - (EXPONE-1)) + ((long )ac3[0].exponent - (EXPONE - 1)); 
            }    
    else {
		mulm( a, b, ac3 );
		lt = ((long)exponent(b) - (EXPONE-1)) + ((long )ac3[0].exponent - (EXPONE - 1));
    }

	/* calculate sign of product */
	ac3[0].sign = ( signof(b) == signof(a) ) ? 0 : -1;
	lt = lt + EXPONE - 1;
	if( lt >= (long)MAXEXP )	goto overf;
	if( lt <= 0 )	goto underf;
	ac3[0].exponent = lt;
//void  qadd_subtract(const Qfloatp a,const Qfloatp b,Qfloatp c ,int subflg)
	/* compare exponents */
	lt = (long) lt - (long) addend->exponent;
	if (lt > 0) { /* put the larger number in ac2 */
		memcpy(ac2,ac3,sizeof(QfloatAccum));
		qmovz(addend,ac1);
		lt = -lt;
	}
	else {
		memcpy(ac1,ac3,sizeof(QfloatAccum));
		qmovz(addend,ac2);
	}
	plarger = ac2;
	psmaller = ac1;

	if( lt != 0 )	{
		if( lt < -NBITS-1 )	{	/* answer same as larger addend */
			if (roundAccum(plarger,lost,subflg)) goto overf;
			pack(plarger,c);
			return 0;
		}

		lost = shift( psmaller , lt); /* shift the smaller number down */
	}
	else { /* exponents were the same, so must compare mantissae */
		i = cmpm( psmaller, plarger );
		if( i == 0 )	{	/* the numbers are identical */
			/* if different signs, result is zero */
			if( psmaller->sign != plarger->sign )
				goto underf;
			/* if exponents zero, result is zero */
			if( psmaller->exponent == 0 )
				goto underf;
			/* if same sign, result is double */
			if( psmaller->exponent >= MAXEXP )			{
				qclear(c);
				if( psmaller->sign != 0 )
					qneg(c);
				goto overf;
			}
			pack(ac3,c);
			c->exponent += 1;
			return 1;
		}
		if( i > 0 )		{	/* put the larger number in ac2 */
			tmpP = psmaller;
			psmaller = plarger;
			plarger = tmpP;
		}
		lost = 0;
	}

	if( psmaller->sign == plarger->sign )	{
		addm( psmaller, plarger );
		subflg = 0;
	}
	else {
		subm( psmaller, plarger );
		subflg = 1;
	}
	if (normlz(plarger,&SC))
		goto underf;
	lt = (long long)plarger->exponent - SC;
	if (lt > MAXEXP)
		goto overf;
	if (lt < 0)
		goto underf;
	plarger->exponent = lt;
	if (roundAccum(plarger,lost,subflg)) goto overf;
	pack( plarger, c );
	return 1;
underf:
	qclear(c);
	return -1;

overf:
	qinfin(c);
	return -2;
}
#if 0
void qfm(Qfloatp a,Qfloatp b,Qfloatp c,Qfloatp d)
{
	Qfloat t[1];
//	qmul(a,b,t);
//	qadd_subtract(t,c,d,ADDITION);
	qfma1(a,b,c,d);
//	if (qcmp(t,d) != 0)
//		printf("Wrong\n");
}
#endif
#ifdef x86_64
void  qsquare(const Qfloatp a,Qfloatp c )
{
	long long lt;
	unsigned __int128 u;
	QfloatAccum ac3[1];
	QfloatAccum ac2[1];

	// Avoid jumps for the normal case.
	if (exponent(a) == 0 )	goto underf;
	/* detect multiplication by small integer a */
    if ((a->mantissa[1] == 0) && (a->mantissa[2] == 0) &&
		(a->mantissa[3] == 0) && (a->mantissa[4] == 0) &&
		(a->mantissa[5] == 0) && (a->mantissa[6] == 0)) {
			u = a->mantissa[0];
			u = u*u;
			memset(ac3,0,sizeof(ac3));
			ac3->mantissa[1] = u >> 64;
			ac3->mantissa[2] = (QELT)u;
			ac3->exponent = 1+a->exponent;
    }    
	else {
		qmovz( a, ac2 );
		//memset(ac3,0,sizeof(ac3));
		ac3[0].exponent = a->exponent;
		square( ac2, ac3, 8 );
	}
	mdnorm(ac3);
	lt = ((long)exponent(a) - (EXPONE)) + ((long )ac3[0].exponent - (EXPONE));
	/* product is always positive */
	ac3[0].sign = 0;
	lt = lt + EXPONE;
	if( lt >= (long)MAXEXP )	goto overf;
	if (lt <= 0) goto underf;
	ac3[0].exponent = lt;
	pack( ac3, c );
	return;
underf:
	qclear(c);
	return ;

overf:
	qinfin(c);
}
#endif
/* Multiply, a has at most WORDSIZE significant bits */

void qmuli(Qfloatp a,Qfloatp b,Qfloatp c )
{
	int lt;
	QfloatAccum act[1];

	if (a->exponent == 0 || b->exponent == 0)	{
		goto underf;
	}
	//qmovz(b,act);
	mulin(a->mantissa[0],b, act);
	mdnorm(act);

	/* calculate sign of product */
	act[0].sign = ( signof(b) == signof(a) ) ? 0 : -1;

	/* calculate exponent */
	lt = ((long)act[0].exponent - (EXPONE-1)) + ((long )exponent(a) - (EXPONE - 1));
	if( lt >= MAXEXP )	goto overf;
	lt = lt + EXPONE - 1;
	if( lt >= MAXEXP )	goto overf;
	if( lt < 0 )	goto underf;
	act[0].exponent = lt;
	pack( act, c );
	return ;

underf:
	qclear(c);
	return ;

overf:
	qinfin(c);
	mtherr( "qmuli", OVERFLOW );
}

/*
;	shift mantissa
;
;	Shifts mantissa area up or down by the number of bits
;	given by the argument sc.
*/

#if defined(x86_64) || defined(NOASM)
int shift(QfloatAccump x,int sc)
{
	QELT *p;
	int lost,toshift;

	if( sc == 0 )
		return(0);

	lost = 0;
	p = (QELT *)x;
	if( sc < 0 ) {
	/* This is where 99% of the execution goes.The four operations
	   when calling this routine have always SC < 0 */
		sc = -sc;
		if (sc >= WORDSIZE) {
			if (sc > NBITS_ACC)
				sc = NBITS_ACC;
			toshift = sc/WORDSIZE;
			sc &= 63;
			if (toshift == 1) {
				p[9] = p[8];
				p[8] = p[7];
				p[7] = p[6];
				p[6] = p[5];
				p[5] = p[4];
				p[4] = p[3];
				p[3] = p[2];
				lost |= p[2];
				p[2]=0;
			}
			else if (toshift == 2) {
				p[9] = p[7];
				p[8] = p[6];
				p[7] = p[5];
				p[6] = p[4];
				p[5] = p[3];
				p[4] = p[2];
				lost |= p[2]|p[3];
				p[2] = p[3] = 0;
			}
			else if (toshift == 3) {
				p[9] = p[6];
				p[8] = p[5];
				p[7] = p[4];
				p[6] = p[3];
				p[5] = p[2];
				lost |= p[2]|p[3]|p[4];
				p[2] = p[3] = p[4] = 0;
			}
			else if (toshift == 4) {
				p[9] = p[5];
				p[8] = p[4];
				p[7] = p[3];
				p[6] = p[2];
				lost |= p[2]|p[3]|p[4]|p[5];
				p[2] = p[3] = p[4] = p[5] = 0;
			}
			else if (toshift == 5) {
				p[9] = p[4];
				p[8] = p[3];
				p[7] = p[2];
				lost |= p[2]|p[3]|p[4]|p[5]|p[6];
				p[2] = p[3] = p[4] = p[5] = p[6] = 0;
			}
			else if (toshift == 6) {
				p[9] = p[3];
				p[8] = p[2];
				lost |= p[2]|p[3]|p[4]|p[5]|p[6]|p[7];
				p[2] = p[3] = p[4] = p[5] = p[6] = p[7] = 0;
			}
			else if (toshift == 7) {
				p[9] = p[2];
				lost |= p[2]|p[3]|p[4]|p[5]|p[6]|p[7]|p[8];
				p[2] = p[3] = p[4] = p[5] = p[6] = p[7] = p[8] = 0;
			}
			else {
				memmove(p+2+toshift,p+2,(ACCUM_LENGTH-toshift)*sizeof(QELT));
				sc -= WORDSIZE*toshift;
				while (toshift) {
					toshift--;
					lost |= p[2+toshift];
					p[2+toshift] = 0;
				}
			}
		}
		if (sc > 0) {
			lost |= x->mantissa[ACCUM_LENGTH-1];
			shiftdownn(x,sc);
		}
	}
	else {
	/* This shifts with SC # 0 are called in somehow pathological 
	   cases, like sin(1e80) for instance */
		if (sc >= WORDSIZE) {
			if (sc >= NBITS) {
				memset(x->mantissa,0,sizeof(x->mantissa));
				x->mantissa[1] = 1ULL << (WORDSIZE-1);
				return 1;
			}
			toshift = sc/WORDSIZE;
			memmove(p+2,p+2+toshift,((ACCUM_LENGTH-1)-toshift)*sizeof(QELT));
			sc -= WORDSIZE*toshift;
			while (toshift) {
				toshift--;
				lost |= p[ACCUM_LENGTH-toshift];
				p[ACCUM_LENGTH-toshift] = 0;
			}
		}
		if (sc) {
			shiftupn(x,sc);
		}
	}
	return( lost );
}
#endif

/*
; Fill entire number, including exponent and mantissa, with
; largest possible number.
 */

void qinfin(Qfloatp x)
{
	int i;

	x->exponent = MAXEXP;
	for( i=0; i<MANTISSA_LENGTH; i++ )
		x->mantissa[i] = -1LL;
}

/* normalization program */
void qnrmlz(Qfloatp x)
{
	QfloatAccum ac1[1];
	int SC;

	qmovz( x, ac1 );
	normlz( ac1,&SC );	/* shift normalize the mantissa */
	ac1[0].exponent -= SC;	/* subtract the shift count from the exponent */
	pack( ac1, x );
}

static void TrimTrailingZeroes(char *string)
{
	char *p = string + (strlen(string)-1);

	if (*p == '0') {
		while (*p == '0') {
			p--;
		}
		if (*p != '.')
			p++;
		*p = 0;
	}
}

/*						qtoasc.c	*/
/* Convert q type number to ASCII string */
static Qfloat ten[1] = {{0,EXPONE+3,{0xa000000000000000ULL,0,0,0,0,0,0}}};
/* Get values for powers of ten.  */
int EXPORT qtoasc(Qfloatp q,char *string,int width,int ndigs,int flags )
{
	int i, k, expon;
	char *s, *ss,*pstr;
	long long digit;
	Qfloat x[1], xt[1], ac4[1];
	Qfloatp tenth;
	Qfloat qzero[1];
	int sign;
	Qfloat *p,*r;


	if (exponent(q) == MAXEXP) {
		for (i=0;i<8;i++) {
			if (q->mantissa[i] != -1)
				break;
		}
		if (i == 8) {
			if (signof(q) == 0)
				strcpy(string,"#inf");
			else
				strcpy(string,"-#inf");
			return 0;
		}
	}
	//memset(ten,0,sizeof(ten));
	//memset(qzero,0,sizeof(qzero));
	//ten[0].exponent = EXPONE+3;
	//ten[0].mantissa[0] = 0xa000000000000000ULL;
	qmov( q, x );
	sign = signof(x);
	x[0].sign = 0;
	expon = 0;

	i = qcmp(qone,x);
	if( i == 0 )
		goto isone;
	if( iszero(x) )
	{
		qclear( x );
		goto isone;
	}

	if( i < 0 )
	{
		k = MAXNTEN;
		p = (Qfloat *)&qtens[0];
		qmov( qone, ac4 );
		qmov( x, xt );
		while( qcmp( ten, x ) <= 0 )
		{
			if( qcmp( p, xt ) <= 0 )
			{
				qdiv( p, xt, xt );
				qmul( p, ac4, ac4 );
				expon += k;
			}
			k >>= 1;
			if( k == 0 )
				break;
			p++;
		}
		qdiv( ac4, x, x );
	}
	else // All numbers smaller than 1 go here
	{
		k = MINNTEN;
		p = (Qfloat *)&qmtens[0];
		r = (Qfloat *)&qtens[0];
		tenth = &qmtens[NTEN];
		while( qcmp( tenth, x ) > 0 )
		{
			if( qcmp( p, x ) >= 0 )
			{
				qmul( r, x, x );
				expon += k;
			}
			k /= 2;
			/* Prevent infinite loop due to arithmetic error: */
			if( k == 0 )
				break;
			p++;
			r++;
		}
		qmuli( ten, x, x );
		expon -= 1;
	}

isone:
	qifrac( x, &digit, x );
	/* The following check handles numbers very close to 10**(2**n)
	* when there is a mistake due to arithmetic error.
	*/
	if( digit >= 10 )
	{
		qdivi( 10, x, x );
		expon += 1;
		digit = 1;
	}
	s = string;
	if( sign != 0 )
		*s++ = '-';
	*s++ = (char )digit + 060;
	*s++ = '.';
	if (expon == -1)
		ndigs--;
	if( ndigs < 0 )
		ndigs = 0;
	if( ndigs > NDEC )
		ndigs = NDEC;
	for( k=0; k<ndigs; k++ )
	{
		qmuli( ten, x, x );
		qifrac( x, &digit, x );
		*s++ = (char )digit + 060;
	}

	*s = '\0';	/* mark end of string */

	/* round off the ASCII string */

	qmuli( ten, x, x );
	qifrac( x, &digit, x );
	if( digit > 4 )
	{
		/* Check for critical rounding case */
		if( digit == 5 )
		{
			if( qcmp( x, qzero ) != 0 )
				goto roun;	/* round to nearest */
			if( (*(s-1) & 1) == 0 )
				goto doexp;	/* round to even */
		}
roun:
		--s;
		k = *s & 0x7f;
		/* Carry out to most significant digit? */
		if( k == '.' )
		{
			--s;
			k = *s & 0x7f;
			k += 1;
			*s = k;
			/* Most significant digit rounds to 10? */
			if( k > '9' )
			{
				*s = '1';
				expon += 1;
			}
			goto doexp;
		}
		/* Round up and carry out from less significant digits. */
		k += 1;
		*s = k;
		if( k > '9' )
		{
			*s = '0';
			goto roun;
		}
	}

doexp:
	if (1) {
		if (expon >= 0) {
			pstr = strchr(string,'.');
			int iexp = expon;

			if (expon < strlen(string)-2) {
				for (i=0;i<iexp;i++) {
					*pstr = pstr[1];
					pstr++;
				}
				*pstr = '.';
				if (1)
					TrimTrailingZeroes(string);
				return 0;
			}
			pstr[0] = 'x';
			if (1)
				TrimTrailingZeroes(string);
			pstr[0] = '.';
			ss = string+strlen(string)-1;
			if (*ss == '.') {
				ss++;
				*ss++ = '0';
			}
			else ss++;
		}
		else if (expon < 0) {
			pstr = strchr(string,'.');
			int iexp = -expon;

			if (iexp < ndigs-2) {
				pstr--;
				pstr[1] = pstr[0];
				pstr[0] = '.';
				iexp--;
				if (1)
					TrimTrailingZeroes(string);
				if (iexp > 0) {
					memmove(pstr+iexp+1,pstr+1,1+strlen(string));
					memset(pstr+1,'0',iexp);
					memmove(pstr+1,pstr,strlen(string)+1);
					pstr[0] = '0';
				}
				else if (iexp == 0) {
					memmove(pstr+1,pstr,1+strlen(string));
					pstr[0] = '0';
				}
				return 0;
			}
			else {
				pstr[0] = 'x';
				if (1)
					TrimTrailingZeroes(string);
				pstr[0] = '.';
				ss = string+strlen(string)-1;
				if (*ss == '.') {
					ss++;
					*ss++ = '0';
				}
				else ss++;
			}
		}
	}
	sprintf(string+strlen(string),"E%+05d",expon);
	return 0;
}

int qprint(Qfloatp p)
{
	char buf[500];

	qtoasc(p,buf,200,136,0);
	return printf("%s\n",buf);
}
int qcmpto1(const Qfloatp p)
{
	int i;

	if( signof(p) )	{ /* the signs are different */
		return( -1 );
	}
	/* both are the same sign */
    if (p->exponent != EXPONE) {
	     if (p->exponent > EXPONE)
	             return 1;
	     else
	             return -1;
    }
	if (p->mantissa[0] > 0x8000000000000000ULL)
		return 1;
	else if (p->mantissa[0] < 0x8000000000000000ULL)
		return -1;
	for (i=1; i<MANTISSA_LENGTH;i++) {
		if (p->mantissa[i])
			return 1;
	}
	return(0);	/* equality */
}
/*
;								ASCTOQ
;		ASCTOQ.MAC		LATEST REV: 11 JAN 84
;					SLM, 3 JAN 78
;                                       Adapted to lcc-win by jacob navia 1995-2008
;
;	Convert ASCII string to qfloat precision floating point
;
;       Original comments by SLM below.
;		Numeric input is free field decimal number
;		with max of 15 digits with or without
;		decimal point entered as ASCII from teletype.
;	Entering E after the number followed by a second
;	number causes the second number to be interpreted
;	as a power of 10 to be multiplied by the first number
;	(i.e., "scientific" notation).
;
;	Usage:
;		asctoq( string, q );
*/
int asctoq(char * s, Qfloatp y ,char **pend)
{
	QfloatAccum yy[1];
	Qfloat qt[1];
	int esign, nsign, decflg, nexp, exp, prec, base, k,SC;
	long lexp;
	Qfloatp p;
	QfloatAccum ac2[1];
	char *save = s;

	nsign = 0;
	esign = 1;
	decflg = 0;
	nexp = 0;
	exp = 0;
	prec = 0;
	memset( yy, 0,sizeof(yy) );
	memset( ac2,0,sizeof(ac2));
	base = 10;
	if (s[0] == '-') {
		nsign = -1;
		s++;
	}
	else if (s[0] == '+')
		s++;
	if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X'))	{
		base = 16;
		s += 2;
	}

nxtcom:
	if( (*s >= '0') && (*s <= '9') )
		k = *s - '0';
	else if (*s >= 'a' && *s <= 'f')
		k = 10 + *s - 'a';
	else if (*s >= 'A' && *s <= 'F')
		k = 10 + *s - 'A';
	else
		k = -1;

	if( (k >= 0) && (k < base) ){
		if( (prec == 0) && (decflg == 0) && (*s == '0') )	goto donchr;
		if( prec < NDEC)	{
			if( decflg )		nexp += 1;	/* count digits after decimal point */
			if( base == 16 )		{
				shiftupn(yy,4);
			}
			else {
				shup1( yy );	/* multiply current number by 10 */
				memcpy(ac2, yy, sizeof(ac2)) ;
				shiftupn(ac2,2);
				addm( ac2, yy );
			}
			memset( ac2,0,sizeof(ac2));
			ac2[0].mantissa[MANTISSA_LENGTH] = k;

			addm( ac2, yy );
		}
		prec += 1;
		goto donchr;
	}

	switch( *s ){
	case ' ':
		goto daldone;
	case 'E':
	case 'e':
	case 'P':
	case 'p':
		goto expnt;
	case '.':	/* decimal point */
		if( decflg )		goto error;
		++decflg;
		break;
	case '\0':
#ifndef OSK
	case '\n':
#endif
	case '\r':
		goto daldone;
	default:
		goto daldone;
error:
		printf( "asctoq conversion error: %.10s\n" ,save);
		qclear(y);
		return 0;
	}
donchr:
	++s;
	goto nxtcom;


	/*		EXPONENT INTERPRETATION */
expnt:

	/* 0.0eXXX is zero, regardless of XXX.  Check for the 0.0. */
	for( exp = 0; exp < 8; exp++ )
	{
		if( yy[0].mantissa[exp] != 0 )		goto read_expnt;
	}
	qclear(y);
	goto aexit;

read_expnt:
	exp = 0;
	++s;
	/* check for + or - */
	if( *s == '-' )	{
		esign = -1;
		++s;
	}
	if( *s == '+' )	++s;
	while( (*s >= '0') && (*s <= '9') )	{
		/* Check for oversize decimal exponent.  */
		if( exp >= 32767 || exp < 0 )		{
			if( esign < 0 )			goto zero;
			else
			    goto infinite;
		}
		exp *= 10;
		exp += *s++ - '0';
	}
	if( esign < 0 )	exp = -exp;

daldone:
	if (base == 16)	{
		/* Base 16 hexadecimal floating constant.  */
		normlz (yy,&SC);
		if (SC > NBITS)		{
			memset(yy,0,sizeof(yy));
			goto aexit;
		}
		/* Adjust the exponent.  NEXP is the number of hex digits,
		EXP is a power of 2.  */
		lexp = (EXPONE - 1 + NBITS) - SC + yy[0].exponent + exp - 4 * nexp;
		if (lexp > MAXEXP)		goto infinite;
		if (lexp < 0)		goto zero;
		yy[0].exponent = lexp;
		yy[0].sign = nsign;
		pack( yy, y );
		goto aexit;
	}

	nexp = exp - nexp;

	if( normlz(yy,&SC) )	{
		qclear(y);
		return 0;
	}

	yy[0].exponent = EXPONE - 1 + NBITS - SC;
	yy[0].sign = nsign;
	pack( yy, y ) ;

	/* Escape from excessively large exponent.  */
	if( nexp >= 2 * MAXNTEN )	{
infinite:
		qinfin(y);
		mtherr( "asctoq", OVERFLOW );
		goto aexit;
	}
	if( nexp <= -2 * MAXNTEN )	{
zero:
		qclear(y);
		mtherr( "asctoq", UNDERFLOW );
		printf("exponent=%d, input: %.20s\n",nexp,save);
		goto aexit;
	}

	/* multiply or divide by 10**NEXP */
	if( nexp == 0 )	goto aexit;
	esign = 0;
	if( nexp < 0 )	{
		esign = -1;
		nexp = -nexp;
	}

	p = &qtens[NTEN];
	exp = 1;
	qmov( qone, qt );

	do
	    {
		if( exp & nexp )
			qmul( p, qt, qt );
		exp <<= 1;
		p --;
	}
	while( exp <= MAXNTEN );

	if( esign < 0 )	qdiv( qt, y, y );
	else
	    qmul( qt, y, y );
aexit:
	return 0;
}

/*
; Q type to 80-bit IEEE long double precision
;	long double d;
;	QELT q[N+2];
;	qtoe( q, &d );
 */
void qtoe64(Qfloatp x,unsigned short * e )
{
	int SC;
	int exponent,sign,isExpanded=0;
	unsigned long long *pll = (unsigned long long *)e;
	QELT ac1[ACC_LEN];

	e += 4;

	*e = 0;	/* output high order */
	sign = x->sign;
	if( sign )
		*e = 0x8000;	/* output sign bit */
	if (x->exponent == 0) {
		qmovz( x, (QfloatAccump)&ac1[0] );
		isExpanded = 1;
		if( normlz((QfloatAccump)ac1,&SC) )	goto o64zero;
		ac1[1] -= SC;
	}

	exponent = x->exponent - (EXPONE-0x3fff);	/* adjust exponent for offsets */
	if( exponent >= 0x7fff)	{
		*e-- |= 0x7ffe;
		*e-- = 0xffff;
		*e-- = 0xffff;
		*e-- = 0xffff;
		*e-- = 0xffff;
		return ;
	}
	/* We can't handle denormal numbers if the q-type exponenent is 0. */
	if( exponent <= 0 )	{
		if( exponent > -65 )		{
			SC = exponent;
			qmovz(x,(QfloatAccump)ac1);
			shift( (QfloatAccump)ac1, SC );
			exponent = 0;
		}
		else
		    {
o64zero:
			*(--e) = 0;
			*(--e) = 0;
			*(--e) = 0;
			*(--e) = 0;
			return ;
		}
	}

#if 0
	/* round off to nearest or even */
	k = ac1[2+3];
	if( (k & SIGNBIT) != 0 )	{
		if( (k & ((QELT) SIGNBIT - 1)) == 0 )		{
			/* check all less significant bits */
			/* 0 1 2  3   4   5   6 */
			/* S E M 1,2 3,4 5,6 7,8 */
			for( j=2+4; j<=NQ; j++ )
			{
				if( ac1[j] != 0 )
					goto yesrnd;
			}
			/* round to even */
			if( (ac1[4] & 1) == 0 )
				goto nornd;
		}
yesrnd:
		memset( ac2,0,sizeof(ac2) );
		ac2[2+3] = SIGNBIT;
		addm( ac2, ac1 );
		if( ac1[2] )		{
			shdn1(ac1);
			i += 1;
		}
		if( (i == 0) && (ac1[3] & SIGNBIT) )
			i += 1;
	}
#endif


	*e |= exponent;	/* high order output already has sign bit set */
	if (isExpanded) {
		*pll = *(unsigned long long *)&ac1[4];
	}
	else {
		*pll = x->mantissa[0];
	}
}

double normalize(double ld)
{
	 Qfloatp ten;
	 Qfloat q[1];
	 int i = 0;

	 etoq(ld,q);
	 // Multiply by 10^4096
	 ten = &qtens[1];
	 qmul(q,ten,q);
	 ten = &qtens[10]; // 10**1024
	 do {
	         qmul(q,ten,q);
	                         i++;
	                         if (i > 50)
	                                 break;
	 } while (qcmp(q,qone) < 0);
	 ld = qtoe(q, DOROUNDING);
	 return ld;
}
/*
; Convert 80-bit IEEE long double precision to Q type
;	long double d;
;	QELT q[N+2];
;	e64toq( &d, q );
 */

void e64toq(unsigned short *e,Qfloatp y )
{
	int r,SC;
	unsigned long long *p;
	QfloatAccum yy[1];
	int sign,denorm;


	denorm = 0;	/* flag if denormalized number */

	e += 4;

	r = *e;
	sign = 0;
	if( r & 0x8000 )
		sign = -1;

	r &= 0x7fff;
	/* If zero exponent, then the mantissa is denormalized. */
	if( r == 0 )	{
		denorm = 1;
	}
	r += (EXPONE-0x3fff);
	p = ((unsigned long long *) (e-4)) ;
	if( denorm )	{ /* if zero exponent, then normalize the mantissa */
		memset(yy,0,sizeof(QfloatAccum));
		yy->exponent = r;
		yy->mantissa[1] = *p;
		yy->sign = sign;
		/* For Intel long double, shift denormal significand up 1
		-- but only if the top significand bit is zero.  */
		if( (yy->mantissa[2] & 0x8000000000000000ULL) == 0 )
			shup1( yy );
		if( normlz( yy ,&SC) )
			memset(yy,0,sizeof(yy));
		else
		    yy->exponent -= SC;
		pack(yy,y);
	}
	else {
		y->mantissa[0] = *p;
		y->sign = sign;
		y->exponent = r;
		memset(&y->mantissa[1],0,(MANTISSA_LENGTH-1)*sizeof(QELT));
	}
}

void qldexp(Qfloatp x,long n,Qfloatp y)
{
	 long k;
	 long e = exponent(x);
	 long s = signof(x);

	 k = e  +  n;
	 qmov( x, y );
	 y->sign = s;
	 y->exponent = k;
	 if( (k >= MAXEXP) || (n > (2 * (long)MAXEXP)) )
	         qinfin(y);
	 if( (k <= 0) || (n < (-2 * (long)MAXEXP)) )
	         qclear(y);
}


void qfrexp(Qfloatp x,long *e,Qfloatp y)
{

	 if( x[0].exponent == 0 )        {
	         *e = 0;
	         qclear(y);
	 }
	 else {
	         int expo = exponent(x);
	         int sign = signof(x);
	         *e = (long) expo - (long) EXPONE + 1;
	         qmov( x, y );
	         y[0].exponent = (EXPONE - 1);
	         y[0].sign = sign;
	 }
}

void i128toq(__int128 ll,Qfloatp result)
{
	int SC=0;
	unsigned __int128 longl;
	QfloatAccum ac1[1];

	if (ll == 0) {
		qclear(result);
		return;
	}
	memset(ac1,0,sizeof(ac1));
	if (ll < 0) {
		ll = -ll;
		ac1[0].sign = -1;
	}
	ac1->exponent = EXPONE - 1;
	longl = ll;
	ac1[0].mantissa[1] = longl >> 64;
	ac1[0].mantissa[2] = longl;
	if (!normlz( ac1,&SC)) {
		ac1[0].exponent += (128 - SC);
	}
	pack(ac1,result);
}
/*
; Q type to IEEE double precision
;       double d;
;       QELT q[NQ];
;       double d = qtoe( q, DOROUNDING/NOROUNDING );
 */
double qtoe(Qfloatp x, int doround )
{
	int j;
	long long i;
	QELT *p,k;
	QfloatAccum ac1[1];
	int SC;
	union ud ud;
	QELT p0,p1;
	int denorm = 0;

	/* output sign bit */
	if (x->sign) ud.d.sign =  1;
	else ud.d.sign = 0;
	// suppose numbers are normalized 
	if ((x->mantissa[0]&0x8000000000000000ULL)==0) {
		qmovz( x, ac1);
		if( normlz(ac1,&SC) ) return 0.0;
		ac1->exponent -= SC;
		denorm = 1;
	}
	i = (long) x->exponent - (EXPONE - 1023);      /* adjust exponent for offsets */
	ud.d.exponent = i;
	/* Handle denormalized small numbers.  */
	if( i <= 0 )    {
		if (denorm == 0) qmovz( x, ac1 );
		if( i > -53 ) {
			SC = i - 1;
			shift( ac1,SC );
			i = 0;
			denorm = 1;
		}
		else {
			return 0.0;
		}
	}

	if (doround) {
		 /* round off to nearest or even */
		if (denorm == 0) qmovz(x,ac1);
		denorm = 1;
		 k = ac1[0].mantissa[1];
		 if( (k & 0x8000) != 0 )  {
			QfloatAccum ac2[1];
			if( (k & 0x07fff) == 0x4000 ) {
				if( (k & 0x8000) != 0 ) {
				/* check all less significant bits */
				for( j=2; j<=MANTISSA_LENGTH; j++ ) {
					if( ac1->mantissa[j] != 0 )
						goto yesrnd;
					}
				}
			goto nornd;
			}
yesrnd:
			memset( ac2 ,0, sizeof(ac2));
			ac2->mantissa[1] = 0x8000;
			addm( ac2, ac1 );
			if ( ac1->mantissa[0] ) {
				shdn1(ac1);
				i += 1;
	         }
	         if( (i == 0) && (ac1[0].mantissa[1] & SIGNBIT) )    
				i += 1;
		}
	}
nornd:

	if( i > 2047 )  {       /* Saturate at largest number less than infinity. */
	        mtherr( "qtoe", OVERFLOW );
	        return HUGE_VAL;
	}

	if (denorm) {
    	p = &ac1[0].mantissa[0];
		SC = 5;
		shift( ac1,SC );
		ud.d.mantissa0 = ((p[0]&0xf)<< 16)|((p[1]>> 48)&0xffff);
		ud.d.mantissa1 = (p[1]>> 16) &0xffffffff;
	}
	else {
		// Do a shiftupn of 5 for only two QELTs
		p0 = x->mantissa[0] >> (WORDSIZE-5);
		p1 = (x->mantissa[0] << 5) | (x->mantissa[1] >> (WORDSIZE-5));
		ud.d.mantissa0 = ((p0&0xf)<< 16)|((p1>> 48)&0xffff);
		ud.d.mantissa1 = (p1>> 16) &0xffffffff;
	}
	return ud.dd;
}

double qtoefast(Qfloatp x)
{
	union ud ud;
	QELT p0,p1;

	/* output sign bit */
	ud.d.sign = x->sign ?  1 : 0;
	// suppose numbers are normalized 
	ud.d.exponent = (long) x->exponent - (EXPONE - 1023);      /* adjust exponent for offsets */
	// Do a shiftupn of 5 for only two QELTs
	p0 = x->mantissa[0] >> (WORDSIZE-5);
	p1 = (x->mantissa[0] << 5) | (x->mantissa[1] >> (WORDSIZE-5));
	ud.d.mantissa0 = ((p0&0xf)<< 16)|((p1>> 48)&0xffff);
	ud.d.mantissa1 = (p1>> 16) &0xffffffff;
	return ud.dd;
}

/*
;       Negate
 */

#undef qneg
void qneg(Qfloatp x)
{
	 if (x[0].exponent == 0)
	         return;

	 x[0].sign  = (x[0].sign == 0) ? -1 : 0;
}
#define qneg(z) (z)[0].sign ^= -1;

/*
; Convert IEEE double precision to Q type
;       double d;
;       QELT q[NQ];
;       etoq( &d, q );
 */

int etoq(double d,Qfloatp y )
{
	QfloatAccum yy[1];
	int denorm;
	int SC;
	QELT UpperBits;
	_Double *dd;
	QELT m0,m1;

	if (d == 0) {
		qclear(y);
		return 0;
	}
	denorm = 0;     /* flag if denormalized number */
	dd = (_Double *)&d;
	/*
	r = *e & 0x7fff;
	if( r == 0 )
	return 0;
	*/
	if (dd->sign)
		y->sign = -1;
	else y->sign = 0;

	UpperBits = ((dd->mantissa0 >> 16)&0xf)|0x10;
	y->exponent = dd->exponent +EXPONE - 0x3ff;
	/* If zero exponent, then the mantissa is denormalized.
	* So take back the understood high mantissa bit. */
	if( dd->exponent == 0 )    {
		memset(yy,0,sizeof(yy));
		denorm = 1;
		UpperBits &= ~0x10;
		yy[0].exponent = y->exponent;
		yy[0].sign = y->sign;
	}
	m0 = dd->mantissa0;
	m1 = dd->mantissa1;
	y->mantissa[0] = (m0 << (48-5))|(m1<< (16-5))| (UpperBits << (64-5));
	if( denorm )    { /* if zero exponent, then normalize the mantissa */
		SC = -5;
		shift(yy,SC);
		yy[0].mantissa[1]=y->mantissa[0];
		if( normlz( yy,&SC ) ) {
			qclear(y);
			return 0;
		}
		else
			yy[0].exponent -= SC-1;
	}
	if (denorm) pack( yy, y );
	else memset(y->mantissa+1,0,sizeof(y->mantissa)-sizeof(QELT));
	return 0;
}


int isinfq(Qfloatp a)
{
	 return a->exponent >= MAXEXP;
}

long long qtoll(Qfloatp x)
{
	QfloatAccum ac1[1];
	int SC;

	 qmovz(x,ac1);
	 SC = ac1[0].exponent - (EXPONE - 1);
	 if (SC <= 0)
	         return 0L;
	 if (SC > 63)
	         return 0x8000000000000000LL;
	 shift(ac1, SC);
	 return ac1[0].mantissa[0];
}

__int128 qtoi128(Qfloatp x)
{
	QfloatAccum ac1[1];
	int SC;
	__int128 result;

	qmovz(x,ac1);
	SC = ac1[0].exponent - (EXPONE - 1);
	if (SC <= 0)
		return 0;
	if (SC > 127)
		return ((unsigned __int128)0x8000000000000000ULL << 64);
	shift(ac1,SC);
	result = ac1[0].mantissa[1];
	result <<= 64;
	result |= ac1[0].mantissa[0];
	if (ac1[0].sign)
		result = -result;
	return result;
}

void qabs ( Qfloatp x)
{
	 x->sign = 0;
}

void qinv(const Qfloatp a,Qfloatp r)
{
	long lt;
	int sign = a->sign;
	QfloatAccum ac4[1];
	long exponentOfA;

	if (exponent(a) == 0) goto divoverflow;	
	qmovz( a, ac4 );
	/* We need to save exponent of a since a could be equal to r! */
	exponentOfA = a->exponent;
	lt = EXPONE-1;
	ac4[0].exponent = 4;
	inverse(ac4,r);
	/* calculate sign */
	r->sign = sign;

	/* calculate exponent */
	lt = lt + (long long)r->exponent - 4  - (long long)exponentOfA;
	lt = lt + EXPONE + 1;
	if( lt >= MAXEXP ) goto divoverflow;
	else if( lt <= 0 ) goto divunderf;
	else r->exponent = lt;
	return;
divunderf:
	qclear(r); 
	return;
divoverflow:
	qinfin(r);
}

#if 1
void  qincr(const Qfloatp a,Qfloatp c);
float128_t inversef128(float128_t);
/* Third order newton iteration for inverse. This
   procedure is 40% slower than the "inverse" of qfltbi.c}

   xj+1 = xj + xj*(1 - xj*a) + xj*(1 - xj*a)^2. 

   xj+1 = xj * (1 + (1 - xj*a) + (1 - xj*a)^2). 

   xj+1 = xj * (1 + (1 - xj*a) * (1 + (1 - xj*a))). 
*/
void Inverse(Qfloatp a,Qfloatp y)
{
	float128_t d;
	Qfloat xj[1],z[1],term[1];
	d = qtoe113(a);
	d = inversef128(d);
	e113toq(d,xj);
	
	//qsub(z,qone,term);		//                              1 - xj*a
	if (xj[0].exponent == 0) {
		qmov(xj,y);
		return;
	}
	qfma(xj,a,qminusone,term);			//                                  xj*a
	term->sign = ~term->sign;
	qincr(term,z);		//                         1 + (1 - xj*a)
	qfma(term,z,qone,z);	//           (1 - xj*a) * (1 + (1 - xj*a))
	//qincr(z,z);			//       1 + (1 - xj*a) * (1 + (1 - xj*a))
	qmul(xj,z,xj);			// xj * (1 + (1 - xj*a) * (1 + (1 - xj*a)))

	qfma(xj,a,qminusone,term);			// xj*a
	term->sign = ~term->sign;
	//qsub(z,qone,term);		// 1-xj*a
	qincr(term,z);
	qfma(term,z,qone,z);
	//qincr(z,z);
	qmul(xj,z,y);
}
#endif

void  qincr(const Qfloatp a,Qfloatp c)
{
	long lt;
	long long i;
	int lost,SC;
	QfloatAccum ac1[1];
	QfloatAccum ac2[1];
	QfloatAccump plarger,psmaller,tmpP;

	/* compare exponents */
	lt = (long) a->exponent - (long) EXPONE;
    qmovz(a,ac2);
    memset(ac1,0,sizeof(ac1));
    ac1[0].exponent = EXPONE;
    ac1[0].mantissa[1]=0x8000000000000000ULL;
    if( lt > 0 )    {   /* put the larger number in ac2 */
        lt = -lt;
        plarger = ac2;
        psmaller = ac1;
    }
    else  {
        plarger = ac1;
        psmaller = ac2;
    }

	lost = 0;
	if( lt != 0 )	{
		if( lt < -NBITS-1 )		goto done;	/* answer same as larger addend */

		if (psmaller == ac1 && lt < 0) {
			long n = -lt/64;
			long m = lt&63;
			ac1[0].mantissa[1]=0;
			ac1[0].mantissa[1+n] = (0x8000000000000000ULL) >> (64 - m);
		}

		else lost = shift( psmaller , lt); /* shift the smaller number down */
		
	}
	else { /* exponents were the same, so must compare mantissae */
		i = cmpm( psmaller, plarger );
		if( i == 0 )	{	/* the numbers are identical */
			/* if different signs, result is zero */
			if( a->sign )
				goto underf;
			/* if exponents zero, result is zero */
			if( psmaller->exponent == 0 )
				goto underf;
			/* if same sign, result is double */
			if( psmaller->exponent >= MAXEXP )			{
				qclear(c);
				if( psmaller->sign != 0 )
					qneg(c);
				goto overf;
			}
			qmov(a,c);
			c->exponent += 1;
			return;
		}
		if( i > 0 )		{	/* put the larger number in ac2 */
			tmpP = psmaller;
			psmaller = plarger;
			plarger = tmpP;
		}
	}

	if( a->sign == 0 )	{
		if (lt == 0) {
			__int128 tmp = (__int128)0x8000000000000000ULL +  (__int128)ac2[0].mantissa[1];
			plarger->mantissa[0] = tmp>>64;
			plarger->mantissa[1] = (QELT)tmp;
			if (plarger == ac1)
				memcpy(plarger->mantissa+2,psmaller->mantissa+2,6*sizeof(QELT));
		}
		else
			addm( psmaller, plarger );
	}
	else {
		subm( psmaller, plarger );
	}
	if (normlz(plarger,&SC))
		goto underf;
	lt = (long long)plarger->exponent - SC;
	if (lt > MAXEXP)
		goto overf;
	if (lt < 0)
		goto underf;
	plarger->exponent = lt;
	/* round off */
	if (roundAccum(plarger,lost,0))
		goto overf;
done:
	pack( plarger, c );
	return;

underf:
	qclear(c);
	return;

overf:
	qinfin(c);
}

union ud113 {
	ld113 d;
	float128_t dd;
};
float128_t qtoe113(Qfloatp src)
{
	union ud113 u; 
	unsigned long long m1,m0;

//	u.dd = 3.141592653589793238462643383279502884197169399L;
	u.d.sign = (src->sign >> 31);
	m1 = (src->mantissa[0]) << 1;
	m1 = m1|(src->mantissa[1]>>63);
	m0 = (src->mantissa[1] << 1);
	
	u.d.exponent = src->exponent - (EXPONE-1) + 0x3fff-1;
	u.d.mantissahigh = (m1>>16);
	u.d.mantissalow = (m1 << 48)|(m0 >> 16);
	return u.dd;
}

void e113toq(float128_t ld,Qfloatp r)
{
	union ud113 *u;


	u = (union ud113 *)&ld;
	r->sign = 0 - (u->d.sign);
	r->exponent = u->d.exponent - 0x3fff + 1 + EXPONE-1;
	r->mantissa[0] = ((unsigned long long)u->d.mantissahigh << 15)|0x8000000000000000ULL;
	r->mantissa[0] |= u->d.mantissalow >> 49;
	r->mantissa[1] = (u->d.mantissalow << 15);
	r->mantissa[2] = r->mantissa[3] = r->mantissa[4] = r->mantissa[5] = 0;
	r->mantissa[6] = 0;
}
