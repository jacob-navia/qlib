/*							qrand.c
*
 *	Pseudorandom number generator
*
 *
 *
 * SYNOPSIS:
*
 * int qfrand( q );
* QELT q[NQ];
*
 * qfrand( q );
*
 *
 *
 * DESCRIPTION:
*
 * Yields a random number 1.0 <= q < 2.0.
*
 * A three-generator congruential algorithm adapted from Brian
* Wichmann and David Hill (BYTE magazine, March, 1987,
* pp 127-8) is used to generate random 16-bit integers.
* These are copied into the significand area to produce
* a pseudorandom bit pattern.
*/


#include "qhead.h"

/* This function implements the three
congruential generators.  */

static long long sx = 1;
static long long sy = 10000;
static long long sz = 3000;

int ranwh(void)
{
	long long r, s;

	/*  sx = sx * 171 mod 30269 */
	r = sx/177;
	s = sx - 177 * r;
	sx = 171 * s - 2 * r;
	if( sx < 0 )
		sx += 30269;


	/* sy = sy * 172 mod 30307 */
	r = sy/176;
	s = sy - 176 * r;
	sy = 172 * s - 35 * r;
	if( sy < 0 )
		sy += 30307;

	/* sz = 170 * sz mod 30323 */
	r = sz/178;
	s = sz - 178 * r;
	sz = 170 * s - 63 * r;
	if( sz < 0 )
		sz += 30323;
	/* The results are in static sx, sy, sz. */
	return 0;
}


int qfrand(Qfloatp q )
{
	QELT r;
	int i;

	/* Positive sign, exponent of 1.0, clear the high guard word.  */
	q->sign = 0;
	q->exponent = EXPONE;

	/* Fill significand with pseudorandom patterns.  */
	for( i=0; i<MANTISSA_LENGTH; i++ )
	{
		ranwh();
		r = ((sx * sy) + sz) & 0xffff;
		ranwh();
		r = r | (((sx * sy) + sz) << 32L);
		q->mantissa[i] = r;
	}

	/* Ensure the significand is normalized.  */
	q->mantissa[0] |= 0x80000000;
	return 0;
}
