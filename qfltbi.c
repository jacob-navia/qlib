/* NOTE, this version of qfltb.c uses the 128-bit `__int128'
type.  It can be used only with WORDSIZE = 64 bits.  */
//#define NOASM
/*
* Utilities for extended precision arithmetic, called by qflt.c.
* These should all be written in machine language for speed.
*
 * addm( x, y )		add significand of x to that of y
* shdn1( x )		shift significand of x down 1 bit
* shdn8( x )		shift significand of x down 8 bits
* shdn16( x )		shift significand of x down 16 bits
* shup1( x )		shift significand of x up 1 bit
* shup8( x )		shift significand of x up 8 bits
* shup16( x )		shift significand of x up 16 bits
* divm( a, b )		divide significand of a into b
* mulm( a, b )		multiply significands, result in b
* mdnorm( x )		normalize and round off
*
 * Copyright (c) 1984 - 1988 by Stephen L. Moshier.  All rights reserved.
   Rewritten by Jacob navia for ARM 64 2018
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "qhead.h"
#if WORDSIZE != 64
Error qfltbi.c works only with WORDSIZE 64.
#endif
#ifdef NOASM
char *__Accelerators = "";
#else
char *__Accelerators = " with asm accelerators";
#endif

void mdnorm(QfloatAccump x);
void subm2(QfloatAccump x,QfloatAccump y);
#if !defined(NOASM) && !defined(x86_64)
void mulv4(QfloatAccump a,QfloatAccump b,QfloatAccump c);
void square4(QfloatAccump a,QfloatAccump b);
void subm4(QfloatAccump x,QfloatAccump y);
#endif
/*
;	Shift mantissa up by 1 bit
*/
#ifdef NOASM
void shup1(QfloatAccump x)
{
	QELT newbits,bits;
	int i;

	bits = x->mantissa[ACCUM_LENGTH] >> (WORDSIZE-1);
	x->mantissa[ACCUM_LENGTH] <<= 1;
	for( i=ACCUM_LENGTH-1; i>0; i-- ) {
		newbits = x->mantissa[i] >> (WORDSIZE - 1);
		x->mantissa[i] <<= 1;
		x->mantissa[i] |= bits;
		bits = newbits;
	}
	x->mantissa[0] <<= 1;
	x->mantissa[0] |= bits;
}
#endif
#ifdef NOTUSED
static void shup32(QfloatAccump qx)
{
	int i;
	QELT newbyt, oldbyt;

	oldbyt = 0;

	for( i=ACCUM_LENGTH; i>=0; i-- ) {
		newbyt = qx->mantissa[i] >> (WORDSIZE-32);
		qx->mantissa[i] <<= 32;
		qx->mantissa[i] |= oldbyt;
		oldbyt = newbyt;
	}
}
#endif

#ifdef NOASM
/*
;	Shift mantissa down by 1 bit
*/

void shdn1(QfloatAccump qx)
{
	QELT newbits;
	QELT bits, u;
	int i;

	bits = 0;
	for( i=0; i< ACCUM_LENGTH; i++ ) {
		u = qx->mantissa[i];
		newbits = u << (WORDSIZE - 1);
		u >>= 1;
		u |= bits;
		bits = newbits;
		qx->mantissa[i] = u;
	}	
}
#endif

#ifdef NOASM 
void shiftdown2(QfloatAccump x)
{
	QELT newbyt, oldbyt;
	int i;

	oldbyt = 0;
	for( i=0; i<=ACCUM_LENGTH; i++ )
	{
		newbyt = x->mantissa[i] << (WORDSIZE-2);
		x->mantissa[i] >>= 2;
		x->mantissa[i] |= oldbyt;
		oldbyt = newbyt;
	}
}
#else
void shiftdown2(QfloatAccump x);
#endif

#ifdef NOASM
void shiftdownn(QfloatAccump x,int n)
{
	QELT newbyt, oldbyt;
	int i;

	if (n >= 64)
		printf("n=%d\n",n);
	oldbyt = 0;
	for( i=0; i<=ACCUM_LENGTH; i++ )
	{
		newbyt = x->mantissa[i] << (WORDSIZE-n);
		x->mantissa[i] >>= n;
		x->mantissa[i] |= oldbyt;
		oldbyt = newbyt;
	}
}

void shiftupn(QfloatAccump x,int n)
{
	int i;
	QELT newbyt, oldbyt;

	oldbyt = 0;
	for( i=ACCUM_LENGTH; i>= 0; i-- )
	{
		newbyt = x->mantissa[i] >> (WORDSIZE - n);
		x->mantissa[i] <<= n;
		x->mantissa[i] |= oldbyt;
		oldbyt = newbyt;
	}
}
#endif
#if defined(NOASM)
void divi(unsigned long long src,Qfloatp qb, QfloatAccump b)
{
	unsigned __int128 u,d,qu;

	if (src == 0) {
		memset(b,0,sizeof(QfloatAccum));
		return;
	}
	d = src;
	qmovz(qb,b);
	/* Do single precision divides if so. */
	shiftdown2(b);
	u = ((unsigned __int128)b->mantissa[1] << 64) | b->mantissa[2];
	// 2
	qu = u/d;
	b->mantissa[1] = (QELT)qu;
	u = ((u - d * qu) << 64) | b->mantissa[3];
	// 3
	qu = u/d;
	b->mantissa[2] = (QELT)qu;
	u = ((u - d * qu) << 64) | b->mantissa[4];
	// 4
	qu = u/d;
	b->mantissa[3] = (QELT)qu;
	u = ((u - d * qu) << 64) | b->mantissa[5];
	// 5
	qu = u/d;
	b->mantissa[4] = (QELT)qu;
	u = ((u - d * qu) << 64) | b->mantissa[6];
	// 6
	qu = u/d;
	b->mantissa[5] = (QELT)qu;
	u = ((u - d * qu) << 64) | b->mantissa[7];
	// 7
	qu = u/d;
	b->mantissa[6] = (QELT)qu;
	u = ((u - d * qu) << 64) | b->mantissa[8];
	// 8
	qu = u/d;
	b->mantissa[7] = (QELT)qu;
	u = ((u - d * qu) << 64);
	// 9
	qu = u/d;
	b->mantissa[8] = (QELT)qu;
}
/*
extern unsigned __int128 udivmodti4(unsigned __int128 a,unsigned __int128 b,unsigned __int128 *r);
void divi(unsigned long long src, QfloatAccump b)
{
	unsigned __int128 u,d,qu,mod;

	d = src;
	shiftdown2(b);
	u = ((unsigned __int128)b->mantissa[1] << 64) | b->mantissa[2];
	// 2
	qu = udivmodti4(u,d,&mod);
	b->mantissa[1] = (QELT)qu;
	u = (mod << 64) | b->mantissa[3];
	// 3
	qu = udivmodti4(u,d,&mod);
	b->mantissa[2] = (QELT)qu;
	u = (mod << 64) | b->mantissa[4];
	// 4
	qu = udivmodti4(u,d,&mod);
	b->mantissa[3] = (QELT)qu;
	u = (mod << 64) | b->mantissa[5];
	// 5
	qu = udivmodti4(u,d,&mod);
	b->mantissa[4] = (QELT)qu;
	u = (mod << 64) | b->mantissa[6];
	// 6
	qu = udivmodti4(u,d,&mod);
	b->mantissa[5] = (QELT)qu;
	u = (mod << 64) | b->mantissa[7];
	// 7
	qu = udivmodti4(u,d,&mod);
	b->mantissa[6] = (QELT)qu;
	u = (mod << 64) | b->mantissa[8];
	// 8
	qu = udivmodti4(u,d,&mod);
	b->mantissa[7] = (QELT)qu;
	u = (mod << 64) | b->mantissa[9];
	// 9
	qu = udivmodti4(u,d,&mod);
	b->mantissa[8] = (QELT)qu;
}
*/
#endif
#ifdef NOASM
/* Variable precision multiply of significands.
* c must not be in the same location as either a or b.
*/
static void mulv(QfloatAccump a,QfloatAccump b,QfloatAccump c,int prec)
{
	QELT *p, *q, *r;
	unsigned __int128 u, lp;
	int k, i;

	k = prec+2;
	memset(c->mantissa,0,k*sizeof(QELT));

	r = &c->mantissa[k];
	for( ; k>=2; k--,r-- ) {
		q = &b->mantissa[1];
		p = &a->mantissa[k-1];
		for( i=k; i>=2; i-- ) {
			if( (*p ) && (*q)) {
				lp = (*p);
			    lp *= (unsigned __int128)(*q);
				u =  (lp&0xffffffffffffffffULL);
				u += *r;
				(*r) = (u&0xffffffffffffffffULL);
				u = u >> 64;
				u += (lp >> 64);
				u += *(r-1);
				*(r-1) = u;
				*(r-2) += (u >> 64);
			}
			p--;
			q++;
		}
	}
}
#else
void mulv(QfloatAccump a,QfloatAccump b,QfloatAccump c,int prec);
#endif

#if defined(NOASM) || defined(x86_64)
static void mulv2(QfloatAccump a,QfloatAccump b,QfloatAccump c)
{
	QELT *p, *q, *r;
	unsigned __int128 u, lp;
	int k, i;

	k = 4;
	memset(c->mantissa,0,k*sizeof(QELT));

	r = &c->mantissa[4];

	// 4
	q = &b->mantissa[1];
	p = &a->mantissa[4-1];
	for( i=4; i>=2; i-- ) {
		if( (*p ) && (*q)) {
			lp = (*p);
		    lp *= (unsigned __int128)(*q);
			u =  (lp&0xffffffffffffffffULL);
			u += *r;
			(*r) = (u&0xffffffffffffffffULL);
			u = u >> 64;
			u += (lp >> 64);
			u += *(r-1);
			*(r-1) = u;
			*(r-2) += (u >> 64);
		}
		p--;
		q++;
	}
	r--;
	// 3
	q = &b->mantissa[1];
	p = &a->mantissa[3-1];
	for( i=3; i>=2; i-- ) {
		if( (*p ) && (*q)) {
			lp = (*p);
		    lp *= (unsigned __int128)(*q);
			u =  (lp&0xffffffffffffffffULL);
			u += *r;
			(*r) = (u&0xffffffffffffffffULL);
			u = u >> 64;
			u += (lp >> 64);
			u += *(r-1);
			*(r-1) = u;
			*(r-2) += (u >> 64);
		}
		p--;
		q++;
	}
	r--;
	// 2
	q = &b->mantissa[1];
	p = &a->mantissa[2-1];
	if( (*p ) && (*q)) {
		lp = (*p);
	    lp *= (unsigned __int128)(*q);
		u =  (lp&0xffffffffffffffffULL);
		u += *r;
		(*r) = (u&0xffffffffffffffffULL);
		u = u >> 64;
		u += (lp >> 64);
		u += *(r-1);
		*(r-1) = u;
		*(r-2) += (u >> 64);
	}
}
#else
void mulv2(QfloatAccump a,QfloatAccump b,QfloatAccump c);
#endif
/* Variable precision square.
* b must be in a different location from a.
*/
//#pragma GCC optimize("-O0")
//#pragma GCC optimize("-fwrapv")
#ifdef NOASM
void square(QfloatAccump a,QfloatAccump b,int prec )
{
	QELT *p, *q, *r;
	unsigned __int128 u, lp;
	int k;

	memset(b->mantissa,0,sizeof(b->mantissa));
	b->sign = 0;
	r = &b->mantissa[prec+2];
	for( k=prec+1; k>=1; k-- ) {
		q = &a->mantissa[1];
		p = &a->mantissa[k];
		while( p >= q ) {
			if( (*p) && (*q) ) {
				/*	printf( "%d %d %d\n", p - &a[3], q - &a[3], r - &b[3] );*/
				lp = ((unsigned __int128)(*p)) * (*q);
				if( p != q ) {
					if( ((lp >> 127) & 1 ) )
						*(r-2) += 1;
					lp <<= 1;
				}
				u =  lp&0xffffffffffffffffULL;;
				u += *r;
				(*r) = (u&0xffffffffffffffffULL);
				u = u >> 64;
				u += (lp >> 64);
				u += *(r-1);
				*(r-1) = u;
				*(r-2) += (u >> 64);
			}
			p--;
			q++;
		}
		--r;
	}
	shup1(b);
}
#else
void square(QfloatAccump a,QfloatAccump b,int prec );
#endif
/*
//	x[n+1] = x[n] * (2 - a * x[n]) -->
	x[n+1] = 2*x[n] - a * x[n] * x[n]
*/
//#pragma GCC optimize("O2")
#if defined(x86_64) || defined(NOASM)
void inverse_internal(QfloatAccump qa,QfloatAccump quot)
{
    QfloatAccum sqr[1], prod[1];
	QELT newbits,bits;
    unsigned __int128 u;
	memset(quot,0,sizeof(QfloatAccum));
	memset(prod,0,sizeof(prod));
	memset(sqr,0,sizeof(sqr));
	u = (unsigned __int128)0x4000000000000000ULL << 64;
	//u = u << 64;
	u = (QELT)(u / (unsigned __int128)qa->mantissa[1]);
	quot[0].mantissa[1]=u;
	quot[0].mantissa[2]=0;
	u = u*u;
	//sqr[0].mantissa[0] = u >> 127;
	u <<= 1;
	sqr[0].mantissa[1]= u >> 64;
	sqr[0].mantissa[2] = (QELT)u;
	sqr[0].mantissa[3] = 0;
	sqr[0].mantissa[4] = 0;
	mulv2(qa,sqr,prod);
	subm2(prod,quot);
	//shup1(quot);
	// inline shup1(quot) This improves division by 0.5%
	bits = quot->mantissa[4] >> (WORDSIZE-1);
	quot->mantissa[4] <<= 1;
	newbits = quot->mantissa[3] >> (WORDSIZE - 1);
	quot->mantissa[3] <<= 1;
	quot->mantissa[3] |= bits;
	bits = newbits;

	newbits = quot->mantissa[2] >> (WORDSIZE - 1);
	quot->mantissa[2] <<= 1;
	quot->mantissa[2] |= bits;
    bits = newbits;

	newbits = quot->mantissa[1] >> (WORDSIZE - 1);
	quot->mantissa[1] <<= 1;
	quot->mantissa[1] |= bits;

	quot->mantissa[0] |= newbits;
	square(quot,sqr,4);
	mulv(qa,sqr,prod,4);
	subm(prod,quot);
	shup1(quot);

	square(quot,sqr,8);
	mulv(qa,sqr,prod,8);
	subm(prod,quot);
	shup1(quot);

	square(quot,sqr,8);
	mulv(qa,sqr,prod,8);
	subm(prod,quot);
	shup1(quot);
}
#endif

#ifdef NOASM 
void divm(Qfloatp a, Qfloatp b,QfloatAccump ac3)
{
	QfloatAccum quot[1],qa[1],prod[1],qb[1];

	qmovz(a,qa);
	qmovz(b,qb);
	inverse_internal(qa,quot);
	mulv( quot, qb, prod, 8 );
	memcpy(qb->mantissa,prod->mantissa,sizeof(qb->mantissa));
}
#endif
#if defined(x86_64) || defined(NOASM)
void inverse(Qfloatp qa,Qfloatp qb)
{
	QfloatAccum quot[1],divisor[1];

	memcpy(divisor,qa,sizeof(QfloatAccum)),
	inverse_internal(divisor,quot);
	quot[0].exponent = 4;
	quot[0].sign = qa->sign;
	mdnorm( quot );
	pack(quot,qb);
}
#endif
#ifdef NOASM
/*void divm(QfloatAccump qa,QfloatAccump qb)
{
	QfloatAccum quot[1];

	inverse_internal(qa,quot);
	mulvCopy(quot,qb,qb);
	StoreProd(qb);
}
#endif
#ifdef x86_64
void inverse(QfloatAccump qa,Qfloatp qb)
{
	QfloatAccum quot[1];

	inverse_internal(qa,quot);
	StoreQuot(quot);
	quot[0].exponent = qa->exponent;
	quot[0].sign = qa->sign;
	mdnorm( quot );
	pack(quot,qb);
}
*/

#endif
#ifdef NOASM
void mulm(Qfloatp qa,Qfloatp qb,QfloatAccump ac3)
{
	QELT *p, *q;
	QfloatAccum act[1],b[1];
	QELT *r;
	unsigned __int128 lp, a;
	int i, k;

	qmovz(qa,ac3);
	qmovz(qb,b);
	memset( act[0].mantissa, 0, sizeof(act[0].mantissa) );
	act[0].sign = ac3->sign;
	act[0].exponent = ac3->exponent;
	r = &act[0].mantissa[10];

	q = &ac3->mantissa[1];
	p = &b->mantissa[8];

	for( i=8; i>=1; i-- ) {
		if( (*p ) && (*q) ) {
			lp = (*p);
		    lp *= (unsigned __int128)(*q);
			a =  (lp&0xffffffffffffffffULL);
			a += *r;
			(*r) = (a&0xffffffffffffffffULL);
			a = a >> 64;
			a += (lp >> 64);
			a += *(r-1);
			*(r-1) = a;
			*(r-2) += (a >> 64);
		}
		p--;
		q++;
	}
	r--;
	for( k=9; k>=2; k-- ) {
		q = &ac3->mantissa[1];
		p = &b->mantissa[k-1];

		for( i=k-1; i>=1; i-- ) {
			if( (*p ) && (*q) ) {
				lp = (*p);
			    lp *= (unsigned __int128)(*q);
				a =  (lp&0xffffffffffffffffULL);
				a += *r;
				(*r) = (a&0xffffffffffffffffULL);
				a = a >> 64;
				a += (lp >> 64);
				a += *(r-1);
				*(r-1) = a;
				*(r-2) += (a >> 64);
			}
			p--;
			q++;
		}
		--r;
	}
	memcpy( ac3, act, sizeof(QfloatAccum) );
}
#endif
#if defined(NOASM)
void mulin(QELT y,Qfloatp b,QfloatAccump ac3)
{
	unsigned long long *p, *r;
	QELT t,tm1,tm2;
	unsigned __int128 lp, a;
	int i;

	qmovz(b,ac3);
	if ( ac3->mantissa[8] ) {
		lp = ac3->mantissa[8];
		lp *= y;
		ac3->mantissa[9] = lp;
		t = ac3->mantissa[8] = (lp >> 64);
	}
	else t = 0;
	tm1 = 0;
	tm2 = 0;
	r = &ac3->mantissa[8];
	p = &ac3->mantissa[7];
	for( i=MANTISSA_LENGTH-1; i>=0; i-- ) {
		if( *p ) {
			lp = *p;
			lp *= y;
			a = ((unsigned __int128)(t)+(unsigned long long)lp);
			t = a;
			a = (unsigned __int128)(tm1) + (unsigned __int128)(lp >> 64) + (a >> 64);
			tm1 = a;
			tm2 += a >> 64;
		}
		*r-- = t;
		t = tm1;
		tm1 = tm2;
		tm2 = 0;
		--p;
	}
	*r = t;
	*(r-1) = tm1;
}
#endif

#if defined(NOASM)
static QfloatAccum rndbit[1]={
{0,0,{0,0,0,0,0,0,0,1,0,0,0}}
};

void addbit(QfloatAccump x)
{
	addm(rndbit,x);
}
#endif

#if defined(NOASM) || defined(x86_64)
void mdnorm(QfloatAccump x)
{
	int i;

	while (x->mantissa[0]) {
	// The "while" will be executed at most once. It is "just in case"
	// shdn1 is more efficient than the more general shiftdownn. That makes
	// a better choice.
		shdn1(x);
		if( x->exponent < MAXEXP )
			x->exponent += 1;
		else {
			x->exponent = MAXEXP;
			break;
			//mtherr("mdnorm", OVERFLOW);
		}
	}
	if (0 == ( x->mantissa[1] & 0x8000000000000000ULL ) && x->mantissa[1]) {
		i = bsr64(x->mantissa[1]);
		if (x->exponent >= i) {
			shiftupn(x,i);
			x->exponent -= i;
		}
	}
	if (x->mantissa[MANTISSA_LENGTH+1] & 0x8000000000000000ULL  ) {
		/*
		if(  ((x[NQ] & SIGNBIT) == SIGNBIT) && ((x[NQ-1] & 1) == 0)  )
		goto nornd;
		*/
		addbit( x );
//		x->mantissa[MANTISSA_LENGTH+1]=0;
	}
	// Same reasoning as above: x->mantissa[0] is never different
	// from one.
	while( x->mantissa[0] ) {
		shdn1(x);
		if( x->exponent < MAXEXP )
			x->exponent += 1;
		else {
			x->exponent = MAXEXP;
			break;
			//mtherr("mdnorm", OVERFLOW);
		}
	}
}
#endif
#ifdef NOASM
void addm(QfloatAccump x,QfloatAccump y )
{
	unsigned __int128 a;
	int i,carry=0;

	for( i=ACCUM_LENGTH; i>= 0; i-- ) {
		a = ((unsigned __int128 )x->mantissa[i]) + ((unsigned __int128 )y->mantissa[i]) + carry;
		if( ((unsigned long long) (a >> 64)) & 0x1 )
			carry = 1;
		else
			carry = 0;
		y->mantissa[i] = a;
	}
}
#else
void addm(QfloatAccump x,QfloatAccump y );
#endif

#ifdef NOASM
void subm(QfloatAccump x,QfloatAccump y )
{
	unsigned __int128 a;
	int i,carry;

	carry = 0;
	for( i=ACCUM_LENGTH; i>= 0; i-- )
	{
		a = ((unsigned __int128 )y->mantissa[i]) - ((unsigned __int128 )x->mantissa[i]) - carry;
		if( ((unsigned long long) (a >> 64)) & 0x1 )
			carry = 1;
		else
			carry = 0;
		y->mantissa[i] = a;
	}
}
#else
void subm(QfloatAccump x,QfloatAccump y );
#endif
#if defined(x86_64) || defined(NOASM)
void subm2(QfloatAccump x,QfloatAccump y)
{
	unsigned __int128 a;
	int carry;

	a = ((unsigned __int128 )y->mantissa[4]) - ((unsigned __int128 )x->mantissa[4]);
	if( ((unsigned long long) (a >> 64)) & 0x1 ) carry = 1; else carry = 0;
	y->mantissa[4] = a;
	a = ((unsigned __int128 )y->mantissa[3]) - ((unsigned __int128 )x->mantissa[3]) - carry;
	if( ((unsigned long long) (a >> 64)) & 0x1 ) carry = 1; else carry = 0;
	y->mantissa[3] = a;
	a = ((unsigned __int128 )y->mantissa[2]) - ((unsigned __int128 )x->mantissa[2]) - carry;
	if( ((unsigned long long) (a >> 64)) & 0x1 ) carry = 1; else carry = 0;
	y->mantissa[2] = a;
	a = ((unsigned __int128 )y->mantissa[1]) - ((unsigned __int128 )x->mantissa[1]) - carry;
	if( ((unsigned long long) (a >> 64)) & 0x1 ) carry = 1; else carry = 0;
	y->mantissa[1] = a;
	a = ((unsigned __int128 )y->mantissa[0]) - ((unsigned __int128 )x->mantissa[0]) - carry;
	y->mantissa[0] = a;
	y->mantissa[5] = y->mantissa[6]=0;
}
#endif
#ifdef NOASM
/*
;	move a to b
*/
void  qmov(const Qfloat *a,Qfloat *b)
{
	memcpy(b,a,sizeof(Qfloat));
}
#endif

#ifdef NOASM
void qmovz(Qfloatp a,QfloatAccump b)
{
	b->mantissa[0]=0;
	memcpy(b->mantissa+1,a->mantissa,sizeof(a->mantissa));
	memset(b->mantissa+1+MANTISSA_LENGTH,0,sizeof(b->mantissa)-sizeof(a->mantissa)-sizeof(QELT));
	b->sign = a->sign;
	b->exponent = a->exponent;
}
#endif
/*
; Clear out entire number, including sign and exponent, pointed
; to by x
*/

#ifdef NOASM
void  qclear(Qfloatp x)
{
	memset(x,0,sizeof(Qfloat));
}
#endif

#if defined(NOASM)
int qequal(Qfloatp x,Qfloatp y)
{
	return 0 == qcmp(x,y);
}
#endif

#ifdef NOASM
void pack(QfloatAccump src,Qfloatp dst)
{
	dst->sign = src->sign;
	dst->exponent = src->exponent;
	memcpy(dst->mantissa,src->mantissa+1, sizeof(dst->mantissa));
}
#endif
#if defined(NOASM)
int cmpm(QfloatAccump a,QfloatAccump b)
{
	for (int i = 0; i<MANTISSA_LENGTH+2;i++) {
		if (a->mantissa[i] != b->mantissa[i]) {
			if (a->mantissa[i] > b->mantissa[i])
				return 1;
			else
				return -1;
		}
	}
	return 0;
}
#endif

void prt(Qfloatp x)
{
	char buf[256];

	qtoasc(x,buf,150,132,0);
	printf("%s\n",buf);
}


#ifdef NOASM
void __swapm(QfloatAccump a,QfloatAccump b)
{
	QELT *pa = (QELT *)a,*pb = (QELT *)b,t1,t2;
	int i;

	for (i=0; i<sizeof(QfloatAccum)/sizeof(QELT);i++) {
		t1 = *pa;
		t2 = *pb;
		*pa++ = t2;
		*pb++ = t1;
	}
}
#endif
#if defined(NOASM)
/*
; convert long integer to q type
;
;       long l;
;       QELT q[NQ];
;       itoq( l, q );
; note &l is the memory address of l
;
; Optimized 2018 by J.N. 
; 1) Avoid any copy of the data. 
; 2) Do not call normalize and just normalize immediately We know we have
;    at most 64 bits.
*/
void itoq(long long ll,Qfloatp y)
{
	int toshift;

	if (ll == 0) {
		qclear(y);
		return;
	}
	memset(&y->mantissa[1],0,sizeof(y->mantissa)-sizeof(QELT));
	if( ll < 0 )    {
		ll = -ll;       /* make it positive */
		y->sign = -1; /* put correct sign in the q type number */
	 }
	else y->sign = 0;

	y->exponent = EXPONE - 1;      /* exponent if normalize shift count were 0 */

	toshift =  bsr64(ll);
	y->mantissa[0] = ll << toshift;
	y->exponent += (64-toshift);
}
#endif

#if defined(NOASM) || defined(x86_64)
/*
;	normalize
;
; shift normalizes the mantissa area pointed to by R1
; shift count (up = positive) returned in SC
*/
int normlz(QfloatAccump x,int *SC)
{
	int toshift,sc;

	if (x->exponent == 0)
		return 0;
	if (x->mantissa[0] == 0) {
		sc = 0;
		while (x->mantissa[1] == 0) {
			x->mantissa[1] = x->mantissa[2];
			x->mantissa[2] = x->mantissa[3];
			x->mantissa[3] = x->mantissa[4];
			x->mantissa[4] = x->mantissa[5];
			x->mantissa[5] = x->mantissa[6];
			x->mantissa[6] = x->mantissa[7];
			x->mantissa[7] = x->mantissa[8];
			x->mantissa[8] = 0;
			sc += 64;
			if (sc > NBITS_ACC) {
				*SC = sc;
				return 1;
			}
		}
		if ((x->mantissa[1]&0x8000000000000000LL)==0) {
			toshift = bsr64(x->mantissa[1]);
			*SC = sc+toshift;
			shiftupn(x,toshift);
		}
		else *SC = sc;
	}
	else {
		/* normalize by shifting down out of the high guard word
		of the mantissa */

		toshift = 64 - bsr64( x->mantissa[0] );
		*SC = -toshift;
		shiftdownn(x,toshift);
	}
	return(0);
}
#endif

#if defined(NOASM) || defined(x86_64)
int roundAccum(QfloatAccump plarger,int lost,int subflg)
{
	int SC;
	long lt;

	/* round off */
	QELT i = plarger->mantissa[MANTISSA_LENGTH+1];
	if( i & SIGNBIT )	{
		if( i == SIGNBIT )		{
			if( lost == 0 )			{
			/* Critical case, round to even */
			if( (plarger->mantissa[MANTISSA_LENGTH+2] & 1) == 0 )
				return 0;
			}
			else {
				if( subflg != 0 )
					return 0;
			}
		}
		addbit(plarger);
		normlz(plarger,&SC);
		if( SC ) {
			lt = (long )plarger->exponent - SC;
			if( lt >= (long)MAXEXP )
				 return QLIB_OVERFLOW;
			plarger->exponent = lt;
		}
	}
	return 0;
}
#endif

#if defined(NOASM)
/* same algorithm than itoq */
void lltoq(long long ll,Qfloatp y)
{
	int toshift;
	long long longl;

	longl = ll;
	if (longl == 0) {
		qclear(y);
		return;
	}
	memset(&y->mantissa[1],0,sizeof(y->mantissa)-sizeof(QELT));
	if( longl < 0 )    {
		longl = -longl;       /* make it positive */
		y->sign = -1; /* put correct sign in the q type number */
	 }
	else y->sign = 0;

	y->exponent = EXPONE - 1;      /* exponent if normalize shift count were 0 */

	toshift =  bsr64( longl );
	longl = longl << toshift;
	y->mantissa[0]=longl;
	y->exponent += (64-toshift);
}
#endif

void etoqfast(double e,Qfloatp y )
{
	_Double *dd = (_Double *)&e;
	QELT m0,m1;

	if (dd->sign)
		y->sign = -1;
	else
		y->sign = 0;

	y->mantissa[0] = ((dd->mantissa0 >> 16)&0xf)|0x10;
	y->exponent = dd->exponent +EXPONE - 0x3ff;
	
	m0 = dd->mantissa0;
	m1 = dd->mantissa1;
	y->mantissa[1] = (m0 << 48)|(m1<< 16);
	y->mantissa[0] = (y->mantissa[0] << (WORDSIZE-5))|(y->mantissa[1] >> 5);
	y->mantissa[1] = y->mantissa[2] = y->mantissa[3] = 0;
	y->mantissa[4] = y->mantissa[5] = y->mantissa[6] = 0;
}
#if defined(x86_64) || defined(NOASM)
/* QELT a[NQ], b[NQ];
 * qcmp( a, b )
 *
 *  returns +1 if a > b
 *           0 if a == b
 *          -1 if a < b
*/

int qcmp(Qfloatp p, Qfloatp q )
{
	int i,msign;
	Qfloat r[1];

	if( ( exponent(p) <= (QELT) NBITS)  && ( exponent(q) <= (QELT) NBITS ) )	{
		qadd_subtract( q, p, r,1 );
		if( exponent(r) == 0 )
			return( 0 );
		if( signof(r) == 0 )
			return( 1 );
		return( -1 );
	}

	if( signof(p) != signof(q) )	{ /* the signs are different */
		if( signof(p) == 0 )
			return( 1 );
		else
		    return( -1 );
	}

	/* both are the same sign */
	if( signof(p) == 0 )
		msign = 1;
	else
		msign = -1;

    if (p->exponent != q->exponent) {
	     if (p->exponent > q->exponent)
	             return msign;
	     else
	             return -msign;
    }

	for (i=0; i<MANTISSA_LENGTH;i++) {
		if (q->mantissa[i] != p->mantissa[i])
			goto diffL;
	}

	return(0);	/* equality */
diffL:

	if( (p->mantissa[i]) > q->mantissa[i] )
		return( msign );		/* p is bigger */
	else
		return( -msign );	/* p is littler */
}
#endif

#if defined( NOASM) || defined(x86_64)
/*
;	multiply
;
;	QELT a[NQ], b[NQ], c[NQ];
;	qmul( a, b, c );	c = b * a
Basic Floating Point Multiplication Algorithm
Assuming that the operands are already in the right format, performing
floating point multiplication:
Result = R = X * Y = (-1)Xs (Xm x 2Xe) * (-1)Ys (Ym x 2Ye)
involves the following steps:
(1) If one or both operands is equal to zero, return the result as zero, otherwise:
(2) Compute the sign of the result Xs XOR Ys
(3) Compute the mantissa of the result:
 Multiply the mantissas: Xm * Ym
 Round the result to the allowed number of mantissa bits
(4) Compute the exponent of the result:
Result exponent = biased exponent (X) + biased exponent (Y) - bias
(5) Normalize if needed, by shifting mantissa right, incrementing result exponent.
(6) Check result exponent for overflow/underflow:
 If larger than maximum exponent allowed return exponent overflow
 If smaller than minimum exponent allowed return exponent underflow
 */
void  qmul(const Qfloatp a,const Qfloatp b,Qfloatp c )
{
	int lt;
	QfloatAccum ac3[1];

	// Avoid jumps for the normal case.
	if (exponent(a) == 0 )	goto underf;
	if (exponent(b) == 0) goto underf;
    /* detect multiplication by small integer a */
    if ((a->mantissa[1] == 0) && (a->mantissa[2] == 0) &&
        (a->mantissa[3] == 0) && (a->mantissa[4] == 0) &&
        (a->mantissa[5] == 0) && (a->mantissa[6] == 0)) {
            //qmovz( b, ac3 );
            mulin( a->mantissa[0],b, ac3 );
			mdnorm(ac3);
            lt = ((long long)a->exponent - (EXPONE-1)) + (ac3[0].exponent - (EXPONE - 1)); 
    }    
    else /* detect multiplication by small integer b */
        if ((b->mantissa[1] == 0) && (b->mantissa[2] == 0) &&
            (b->mantissa[3] == 0) && (b->mantissa[4] == 0) &&
            (b->mantissa[5] == 0) && (b->mantissa[6] == 0)) {
                //qmovz( a, ac3 );
                mulin( b->mantissa[0], a, ac3 );
				mdnorm(ac3);
                lt = ((long long)exponent(b) - (EXPONE-1)) + ((long long)ac3[0].exponent - (EXPONE - 1)); 
            }    
    else {
		mulm(a,b,ac3);
		mdnorm(ac3);
        lt = ((long long)exponent(b) - (EXPONE-1)) + ((long long)ac3[0].exponent - (EXPONE - 1));
    }
	/* calculate sign of product */
	ac3[0].sign = ( signof(b) == signof(a) ) ? 0 : -1;
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

#if defined( NOASM) || defined(x86_64)
/*
;	divide
;
;	QELT a[NQ], b[NQ], c[NQ];
;	qdiv( a, b, c );	c = b / a
Basic Floating Point Division Algorithm
(http://meseec.ce.rit.edu/eecc250-winter99/250-1-27-2000.pdf)
Assuming that the operands are already in the right format, performing
floating point multiplication:
Result = R = X / Y = (-1)Xs (Xm x 2Xe) / (-1)Ys (Ym x 2Ye)
involves the following steps:
(1) If the divisor Y is zero return Infinity, if both are zero return NaN
(2) Compute the sign of the result Xs XOR Ys
(3) Compute the mantissa of the result:
 The dividend mantissa is extended to 48 bits by adding 0's to the right of the least
significant bit.
 When divided by a 24 bit divisor Ym, a 24 bit quotient is produced.
(4) Compute the exponent of the result:
Result exponent = [biased exponent (X) - biased exponent (Y)] + bias
(5) Normalize if needed, by shifting mantissa left, decrementing result exponent.
(6) Check result exponent for overflow/underflow:
 If larger than maximum exponent allowed return exponent overflow
 If smaller than minimum exponent allowed return exponent underflow
 */
/* for Newton iteration version:
 * extern short qtwo[];
 * static short qt[NQ] = {0};
 * static short qu[NQ] = {0};
 */
void  qdiv(const Qfloatp a,const Qfloatp b,Qfloatp c )
{
	long lt;
	QfloatAccum ac3[1];

	if (exponent(b) == 0 ) goto divunderf;
	if (exponent(a) == 0) goto divoverflow;	

	/* Avoid exponent underflow in mdnorm.  */
	lt = b->exponent;
	if (a->mantissa[1] == 0 && a->mantissa[2] == 0 &&
		a->mantissa[3] == 0 && a->mantissa[4] == 0 &&
		a->mantissa[5] == 0 && a->mantissa[6] == 0 ) {
		divi(a->mantissa[0],b,ac3); 
	}
	else {
		divm( a, b, ac3 );
	}
	ac3[0].exponent = 4;
	mdnorm(ac3);

	/* calculate sign */
	if( signof(a) == signof(b) )
		ac3[0].sign = 0;
	else
		ac3[0].sign = -1;

	/* calculate exponent */
	lt = lt + (long )ac3[0].exponent -4L - (long )exponent(a);

	lt = lt + (long)EXPONE + 1;
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
#endif
