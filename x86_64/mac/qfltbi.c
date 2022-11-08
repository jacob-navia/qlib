/* NOTE, this version of qfltb.c uses the 64-bit `long long int'
type available in the GNU C compiler.  It can be used only
with WORDSIZE = 32 bits.  */

#define INT64 long long
#define NQ	16
#define ACC_LEN 20
#define MAXEXP 1048576

/* Define nonzero for processors that can shift by 31 bits quickly. */
#define FASTSHIFT 1

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
*/
#include <stdlib.h>
#include <string.h>
#define WORDSIZE 32
/* #define N (NQ-2) */
#if WORDSIZE != 32
Error qfltbi.c works only with WORDSIZE 32.
#endif
#define QELT int
int __shdn1(QELT *x);
int __shup1(QELT *x);
int __subm(QELT *x,QELT *y);
int __addm(QELT *x,QELT *y);
int __mdnorm(QELT x[]);
/*
;	Shift mantissa down by 8 bits
*/

int __shdn8(QELT *x)
{
	register QELT newbyt, oldbyt;
	int i;

	x += 1;
	oldbyt = 0;
	for( i=0; i<NQ-1; i++ )
	{
		newbyt = *x << (WORDSIZE - 8);
		*x >>= 8;
		*x |= oldbyt;
		oldbyt = newbyt;
		++x;
	}
	return 0;
}

/*
;	Shift mantissa up by 8 bits
*/

int __shup8(QELT *x)
{
	int i;
	register QELT newbyt, oldbyt;

	x += NQ;
	oldbyt = 0;

	for( i=0; i<NQ-1; i++ )
	{
		newbyt = *x >> (WORDSIZE - 8);
		*x <<= 8;
		*x |= oldbyt;
		oldbyt = newbyt;
		--x;
	}
	return 0;
}



/*
;	Shift mantissa up by 16 bits
*/

int __shup16(QELT *x)
{
	int i;
	register QELT newbyt, oldbyt;

	x += NQ;
	oldbyt = 0;

	for( i=0; i<NQ-1; i++ )
	{
		newbyt = *x >> 16;
		*x <<= 16;
		*x |= oldbyt;
		oldbyt = newbyt;
		--x;
	}
	return 0;
}

/*
;	Shift mantissa down by 16 bits
*/

int __shdn16(QELT *x)
{
	register QELT newbyt, oldbyt;
	int i;

	x += 1;
	oldbyt = 0;
	for( i=0; i<NQ-1; i++ )
	{
		newbyt = *x << 16;
		*x >>= 16;
		*x |= oldbyt;
		oldbyt = newbyt;
		++x;
	}
	return 0;
}

#if 1
int __shiftdownn(QELT *x,int n)
{
	int r = n;
	while (n >= 16) {
		__shdn16(x);
		n -= 16;
	}
	while (n >= 8) {
		__shdn8(x);
		n -= 8;
	}
	while (n > 0) {
		__shdn1(x);
		n--;
	}
	return r;
}

int __shiftupn(QELT *x,int n)
{
	int r = n;
	while (n >= 16) {
		__shup16(x);
		n -= 16;
	}
	while (n >= 8) {
		__shup8(x);
		n -= 8;
	}
	while (n > 0) {
		__shup1(x);
		n--;
	}
	return r;
}
#endif
#if 1
int __divi(QELT a[],QELT b[])
{
	unsigned long long sqr[11+1], prod[11+1], quot[11+1],qu;
	unsigned long long *pl=a;
	int i, prec, k;
	unsigned __int128 *p;
	unsigned __int128 u,d;

	{
		memcpy(prod,b,sizeof(QELT)*20);
		/* Do single precision divides if so. */
		prod[10] = 0;
		__shdn1( prod );
		__shdn1( prod );
		d = pl[2];
		//p = (__int128 *)(&prod[1]);
		u = ((unsigned __int128)prod[2] << 64) | prod[3];
		//u = *p;
		for( i=2; i<9; i++ )
		{
			unsigned __int128 t1,t2;
			if (d != 0)
				//qu = u / d;
				qu = u/d;
			else
					qu = 0;
			prod[i] = (unsigned long long)qu;
			//u = ((u - (__int128)d * qu) << 64) | prod[i+2];
			t1 = (u-(unsigned __int128)d*qu);
			t2 = t1 << 64;
			u = t2|prod[i+2];
		}
		if (d != 0)
			prod[9] = u / d;
		else
			prod[9] = 0;
	}
	return 0;
}
#endif
/* Variable precision multiply of significands.
* c must not be in the same location as either a or b.
*/
int __mulv(long long a[],long long b[],long long c[],int prec)
{
	long long *p, *q, *r;
	unsigned __int128 u, lp, u1;
	int k, i;

	k = prec+3;
	p = &c[1];
	do
		*p++ = 0;
	while( --k );

	r = &c[prec+3];
	for( k=prec+2; k>=2; k-- )
	{
		q = &b[2];
		p = &a[k];
		for( i=k; i>=2; i-- )
		{
			if( (*p == 0) || (*q == 0) )
			{
				--p;
				++q;
				continue;
			}
			lp = (unsigned __int128)(*p);
		    lp *= (*q);
			p--;
			q++;
			u = (unsigned __int128)(*r);
		   	u += ((unsigned long long) lp);
			*r = u;
			u1 = (unsigned __int128)(*(r-1));
		    u1 += (lp >> 64);
		    u1 +=  (u >> 64);
			*(r-1) = u1;
			*(r-2) += (unsigned long long)(u1 >> 64);
		}
		--r;
	}
	return 0;
}
/* Variable precision square.
* b must be in a different location from a.
*/
int __squarev(long long a[],long long b[],int prec )
{
	long long *p, *q, *r;
	unsigned __int128 u, lp, u1;
	int k;

	k = prec+2;
	if (k <= 0) {
		printf("error in squarev!\nprecision is 0x%x (%d)\n",prec,prec);
		return 0;
	}
	p = &b[1];
	do
		*p++ = 0;
	while( --k );

	r = &b[prec+3];
	for( k=prec+2; k>=2; k-- )
	{
		q = &a[2];
		p = &a[k];
		while( p >= q )
		{
			if( (*p == 0) || (*q == 0) )
			{
				--p;
				++q;
				continue;
			}
			/*	printf( "%d %d %d\n", p - &a[3], q - &a[3], r - &b[3] );*/
			lp = ((__int128)(*p)) * (*q);
			if( p != q )
			{
				int g = ((lp>>64) &0x8000000000000000ULL);
				if( ((lp >> 64) & 0x8000000000000000ULL ) )
					*(r-2) += 1;
				lp <<= 1;
			}
			--p;
			++q;
			//u = (unsigned __int128)(*r) + ((unsigned long long) lp);
			u = (unsigned __int128)(*r);
			u += ((unsigned long long) lp);
			*r = (unsigned long long)u;

			//u1 = (unsigned __int128)(*(r-1)) + (lp >> 64) + (u >> 64);
			u1 = (unsigned __int128)(*(r-1));
			u1 += (lp >> 64);
			u1 += (u >> 64);
			*(r-1) = u1;
			*(r-2) += (unsigned long long)(u1 >> 64);
		}
		--r;
	}
	__shup1(b);
	return 0;
}

#if 1
int __divm(long long a[],long long b[])
{
    long long sqr[12+2], prod[12+2], quot[12+2],qu;
    unsigned long long *pl=a;
    int i, prec, k;
    unsigned __int128 *p;
    unsigned __int128 u,d;
	/* Slower procedure is required */
	memset(quot,0,sizeof(quot));
	memset(prod,0,sizeof(prod));
	memset(sqr,0,sizeof(sqr));
	pl = (long long *)a;
	pl = (long long *)(&pl[2]);
	u = 0x4000000000000000ULL;
	u = u << 64;
	quot[2] = u / ((unsigned long long)a[2]);
	prec = 2;
	k = 1;
	while( prec < 9 )
	{
		k = 2 * k;
		if( k > 9 )
			prec = 9;
		else
			prec = k;
#if 1
		__squarev( quot, sqr, prec );
		__mulv( a, sqr, prod, prec );
#else
		memcpy(sqr,quot,sizeof(sqr));
		__mulm(quot,sqr);
		__mulm(a,sqr);
		memcpy(prod,sqr,sizeof(sqr));
#endif
		__subm( prod, quot );
		__shup1( quot );
	}
#if 1
	__mulv( quot, b, prod, 9 );
#else
	mulm( b, quot);
	memcpy(prod,quot,sizeof(sqr));
#endif
	prod[0] = b[0];
	__mdnorm( prod );
	memcpy( b, prod, sizeof(QELT)*ACC_LEN );
	return 0;
}
#else
int __divm(QELT *a,QELT *b)
{
	QELT sqr[ACC_LEN+2], prod[ACC_LEN+2], quot[ACC_LEN+2];
	int i, prec, k;
	QELT d, qu, *p;
	unsigned INT64 u;

	memset(quot,0,sizeof(quot));
	memset(prod,0,sizeof(prod));
	memset(sqr,0,sizeof(sqr));
	quot[3] = ((unsigned INT64)0x4000000000000000ULL) / a[6];
	prec = 2;
	k = 1;
	while( prec < ACC_LEN-2 )
	{
		k = 2 * k;
		if( k > ACC_LEN-2 )
			prec = ACC_LEN - 2;
		else
			prec = k;
		__squarev( quot, sqr, prec );
		__mulv( a, sqr, prod, prec );
		__subm( prod, quot );
		__shup1( quot );
	}
	__mulv( quot, b, prod, ACC_LEN-2 );
	prod[0] = b[0];
	prod[1] = b[1];
	__mdnorm( prod );
	memcpy( b, prod, sizeof(QELT)*ACC_LEN );
	return 0;
}
static int __mulv(QELT a[],QELT b[],QELT c[],int prec)
{
	register QELT *p, *q, *r;
	register unsigned INT64 u, lp;
	int k, i;

	k = prec+2;
	p = &c[2];
	do
		*p++ = 0;
	while( --k );

	r = &c[prec+3];
	for( k=prec+2; k>=3; k-- )
	{
		q = &b[3];
		p = &a[k];
		for( i=k; i>=3; i-- )
		{
			if( (*p == 0) || (*q == 0) )
			{
				--p;
				++q;
				continue;
			}
			lp = (unsigned INT64)(*p--) * (*q++);
			u = (unsigned INT64)(*r) + ((unsigned int) lp);
			*r = u;
			u = (unsigned INT64)(*(r-1)) + (lp >> 32) + (u >> 32);
			*(r-1) = u;
			*(r-2) += u >> 32;
		}
		--r;
	}
	return 0;
}
/* Variable precision square.
* b must be in a different location from a.
*/
static int __squarev(QELT a[],QELT b[],int prec )
{
	QELT *p, *q, *r;
	register unsigned INT64 u, lp;
	int k;

	k = prec+2;
	p = &b[2];
	do
		*p++ = 0;
	while( --k );

	r = &b[prec+3];
	for( k=prec+2; k>=3; k-- )
	{
		q = &a[3];
		p = &a[k];
		while( p >= q )
		{
			if( (*p == 0) || (*q == 0) )
			{
				--p;
				++q;
				continue;
			}
			/*	printf( "%d %d %d\n", p - &a[3], q - &a[3], r - &b[3] );*/
			lp = (unsigned INT64)(*p) * (*q);
			if( p != q )
			{
				if( (lp >> 32) & 0x80000000 )
					*(r-2) += 1;
				lp <<= 1;
			}
			--p;
			++q;
			u = (unsigned INT64)(*r) + ((unsigned int) lp);
			*r = u;
			u = (unsigned INT64)(*(r-1)) + (lp >> 32) + (u >> 32);
			*(r-1) = u;
			*(r-2) += u >> 32;
		}
		--r;
	}
	__shup1(b);
	return 0;
}
#endif

int __mulm(long long b[],long long ac3[])
{
	long long *p, *q;
	long long act[13];
	long long *r;
	unsigned __int128 lp, a, a1;
	int i, k, m, o;

	memset( act, 0, sizeof(act) );
	act[0] = ac3[0];
	r = &act[11];
	for( k=10; k>=2; k-- )
	{
		if( k == 10 )
		{
			m = 9;
			o = 3;
		}
		else
			{
			m = k;
			o = 2;
		}
		q = &ac3[o];
		p = &b[m];

		for( i=m; i>=o; i-- )
		{
			if( (*p == 0) || (*q == 0) )
			{
				--p;
				++q;
				continue;
			}
			lp = (unsigned __int128)(*p);
		    lp *= (*q);
			p--;
			q++;
			a = (unsigned __int128)(*r);
			a += ((long long) lp);
			*r = a;
			a1 = (unsigned __int128)(*(r-1));
		    a1 += (lp >> 64);
		    a1 += (a >> 64);
			*(r-1) = a1;
			//*(r-2) += a >> 64;
			a1 = a1 >> 64;
			*(r-2) += (unsigned long long)a1;
		}
		--r;
	}
	__mdnorm( act );
	memcpy( ac3, act, sizeof(QELT)*ACC_LEN );
	return 0;
}
#if 1
int __mulin(unsigned long b[],unsigned long long ac3[])
{
	long long *p, *r;
	long long act[13];
	long long y;
	unsigned __int128 lp, a;
	int i;

	memset( act,0, sizeof(act) );
	act[0] = ac3[0];
	r = &act[11];
	p = (unsigned long long *)(b);
	y = p[1];
	p = &ac3[9];
	for( i=10; i>=2; i-- )
	{
		if( *p == 0 )
		{
			--p;
			--r;
			continue;
		}
		//lp = (__int128)(*p--) * y;
		lp = *p;
		p--;
		lp *= y;
		a = ((unsigned __int128)(*r)+(unsigned long long)lp);
		*r = a;
		a = (unsigned __int128)(*(r-1)) + (lp >> 64) + (a >> 64);
		//*(r-1) += (lp >> 64);
		//*(r-1) += (a >> 64);

		*(r-1) = a;
		*(r-2) += a >> 64;
		--r;
	}
	__mdnorm( act );
	memcpy( ac3, act, 10*sizeof(long long) );
	return 0;
}
#endif
static long long rndbit[10]={
0,0,0,0,
0,0,0,0,
1,0
};

int __mdnorm(QELT x[])
{
	int i;
	long long *p = (long long *)x;

	for( i=0; i<64; i++ )
	{
		if( p[1] == 0 )
			break;
		__shdn1(x);
		if( x[1] < MAXEXP )
			x[1] += 1;
		else
		{
			x[1] = MAXEXP;
			break;
			//mtherr("mdnorm", OVERFLOW);
		}
	}
	for( i=0; i<64; i++ )
	{
		if( p[2] & 0x8000000000000000ULL )
			break;
		/* Prevent exponent underflow.
		Rounding may be incorrect when this happens.  */
		if( x[1] >= 1 )
		{
			__shup1(x);
			x[1] -= 1;
		}
		else break; //mtherr("mdnorm", UNDERFLOW);
	}

	if (p[9] & 0x8000000000000000ULL  )
	{
		/*
		if(  ((x[NQ] & SIGNBIT) == SIGNBIT) && ((x[NQ-1] & 1) == 0)  )
		goto nornd;
		*/
		__addm( rndbit, x );
	}
	while( p[1] )
	{
		__shdn1( x );
		if( x[1] < MAXEXP )
			x[1] += 1;
		else
		{
			x[1] = MAXEXP;
			break;
			//mtherr("mdnorm", OVERFLOW);
		}
	}
	p[9] = 0;
	return 0;
}

/*
;	move a to b
*/
#if 1
int  __qmov(QELT *a,QELT *b)
{
	register int i;

	i = NQ;
	do
		{
		*b++ = *a++;
	}
	while( --i );
	return 0;
}

int  __qmovz(QELT *a,QELT *b)
{
	register int i;

	i = NQ;
	do
		{
		*b++ = *a++;
	}
	while( --i );
	*b++ = 0;
	return 0;
}
#endif
/*
; Clear out entire number, including sign and exponent, pointed
; to by x
*/

#if 1
void  __qclear(QELT *x)
{
	register int i;

	for( i=0; i<NQ; i++ )
		*x++ = 0;
}
int __qequal(QELT *x,QELT *y)
{
	return 0 == __qcmp(x,y);
}

/*
;       Shift mantissa down by 1 bit
*/

int __shdn1(QELT *x)
{
        register QELT newbits;
        register QELT bits, u;
        int i;

        x += 2; /* point to mantissa area */

        bits = 0;
        for( i=0; i<NQ-1; i++ )
        {
                u = *x;
                newbits = u << (WORDSIZE - 1);
                u >>= 1;
                u |= bits;
                bits = newbits;
                *x++ = u;
        }
        return 0;
}

/*
;       Shift mantissa up by 1 bit
*/
int __shup1(QELT *x)
{
        register QELT newbits;
        register QELT bits, u;
        int i;

        x += NQ;
        bits = 0;

        for( i=0; i<NQ-1; i++ )
        {
                u = *x;
                newbits = u >> (WORDSIZE - 1);
                u <<= 1;
                u |= bits;
                bits = newbits;
                *x-- = u;
        }
        return 0;
}

/*
;       Compare mantissas
;       Only accumulator should be passed here,
;       not normal numbers
;
;       QELT a[NQ], b[NQ];
;       cmpm( a, b );
;
;       for the mantissas:
;       returns +1 if a > b
;                0 if a == b
;               -1 if a < b
*/
int cmpm(QELT *a,QELT *b )
{
        int i;
        long long *pa=(long long *)a,*pb=(long long *)b;

        pa++;
        pb++;
        for(i=1 ; i<10; i++ )
        {
                if( *pa++ != *pb++ )            goto difrnt;
        }
        return(0);

difrnt:
        if( (unsigned long long) *(--pa) > (unsigned long long) *(--pb) )       return(1);
        else
            return(-1);
}

/*
;       add mantissas
;       x + y replaces y
*/

int __addm(QELT *x,QELT *y)
{
        register unsigned INT64 a;
        int i;
        QELT carry;

        x += NQ;
        y += NQ;
        carry = 0;
        for( i=0; i<NQ-1; i++ )
        {
                a = (unsigned INT64 )(*x) + (unsigned INT64 )(*y) + carry;
                if( ((unsigned int) (a >> 32)) & 0x1 )
                        carry = 1;
                else
                        carry = 0;
                *y = a;
                --x;
                --y;
        }
        return 0;
}

/*
;       subtract mantissas
;       y - x replaces y
*/

int __subm(QELT *x,QELT *y)
{
        register unsigned INT64 a;
        int i;
        QELT carry;

        x += NQ;
        y += NQ;
        carry = 0;
        for( i=0; i<NQ-1; i++ )
        {
                a = (unsigned INT64 )(*y) - (unsigned INT64 )(*x) - carry;
                if( ((unsigned int) (a >> 32)) & 0x1 )
                        carry = 1;
                else
                        carry = 0;
                *y = a;
                --x;
                --y;
        }
        return 0;
}

#endif

