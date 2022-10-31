/*============================================================================

This C source file is part of the SoftFloat IEC/IEEE Floating-point Arithmetic
Package, Release 2b.

Written by John R. Hauser.  This work was made possible in part by the
International Computer Science Institute, located at Suite 600, 1947 Center
Street, Berkeley, California 94704.  Funding was partially provided by the
National Science Foundation under grant MIP-9311980.  The original version
of this code was written as part of a project to build a fixed-point vector
processor in collaboration with the University of California at Berkeley,
overseen by Profs. Nelson Morgan and John Wawrzynek.  More information
is available through the Web page `http://www.cs.berkeley.edu/~jhauser/
arithmetic/SoftFloat.html'.

THIS SOFTWARE IS DISTRIBUTED AS IS, FOR FREE.  Although reasonable effort has
been made to avoid it, THIS SOFTWARE MAY CONTAIN FAULTS THAT WILL AT TIMES
RESULT IN INCORRECT BEHAVIOR.  USE OF THIS SOFTWARE IS RESTRICTED TO PERSONS
AND ORGANIZATIONS WHO CAN AND WILL TAKE FULL RESPONSIBILITY FOR ALL LOSSES,
COSTS, OR OTHER PROBLEMS THEY INCUR DUE TO THE SOFTWARE, AND WHO FURTHERMORE
EFFECTIVELY INDEMNIFY JOHN HAUSER AND THE INTERNATIONAL COMPUTER SCIENCE
INSTITUTE (possibly via similar legal warning) AGAINST ALL LOSSES, COSTS, OR
OTHER PROBLEMS INCURRED BY THEIR CUSTOMERS AND CLIENTS DUE TO THE SOFTWARE.

Derivative works are acceptable, even for commercial purposes, so long as
(1) the source code for the derivative work includes prominent notice that
the work is derivative, and (2) the source code includes prominent notice with
these four paragraphs for those parts of this code that are retained.
---

I have adapted this code to lcc-win, corrected some bugs in the 128B bit part.
Actually only the 128 bit extension is used here.
=============================================================================*/
/*
max normal    1.1897314953572317650857593266280070e+4932
              7ffeffff ffffffff ffffffff ffffffff

min normal    3.3621031431120935062626778173217526e-4932
              00010000 00000000 00000000 00000000

max subnormal 3.3621031431120935062626778173217520e-4932
              0000ffff ffffffff ffffffff ffffffff

min subnormal 6.4751751194380251109244389582276466e-4966
              00000000 00000000 00000000 00000001

Infinity      7fff0000 00000000 00000000 00000000

quiet NaN     7fff8000 00000000 00000000 00000000

signaling NaN 7fff0000 00000000 00000000 00000001
*/
typedef int flag;
typedef char int8;
typedef unsigned short bits16;
typedef short sbits16;
typedef unsigned int bits32;
typedef unsigned long long bits64;
typedef int sbits32;
typedef long long sbits64;
typedef long long int64;
typedef unsigned long long uint64;
typedef int int32;
typedef unsigned int uint32;

#define LIT64( a ) a##LL

typedef unsigned int float32;
typedef unsigned long long float64;

typedef struct {
    unsigned long long low;
    unsigned short high;
} floatx80;

typedef struct _float128_t{
    unsigned long long low, high;
} float128_t;

#define UNLIKELY(x) (__builtin_expect(!!(x), 0))
typedef struct _integer192 {
	unsigned long long a0,a1,a2;
} integer192;
/*----------------------------------------------------------------------------
| Internal canonical NaN format.
*----------------------------------------------------------------------------*/
typedef struct {
    flag sign;
    bits64 high, low;
} commonNaNT;

/*----------------------------------------------------------------------------
| Software IEC/IEEE floating-point underflow tininess-detection mode.
*----------------------------------------------------------------------------*/
enum {
    float_tininess_after_rounding  = 0,
    float_tininess_before_rounding = 1
};

enum {
    float_round_nearest_even = 0,
    float_round_down         = 1,
    float_round_up           = 2,
    float_round_to_zero      = 3
};

enum {
    float_flag_invalid   =  1,
    float_flag_divbyzero =  4,
    float_flag_overflow  =  8,
    float_flag_underflow = 16,
    float_flag_inexact   = 32
};

/*----------------------------------------------------------------------------
| Routine to raise any or all of the software IEC/IEEE floating-point
| exception flags.
*----------------------------------------------------------------------------*/
void float_raise( int8 );

enum {
    FALSE = 0,
    TRUE  = 1
};


/*----------------------------------------------------------------------------
| Floating-point rounding mode, extended double-precision rounding precision,
| and exception flags.
*----------------------------------------------------------------------------*/
int float_rounding_mode = float_round_nearest_even;
int float_exception_flags = 0;

/*----------------------------------------------------------------------------
| Packs the sign `zSign', the exponent `zExp', and the significand formed
| by the concatenation of `zSig0' and `zSig1' into a quadruple-precision
| floating-point value, returning the result.  After being shifted into the
| proper positions, the three fields `zSign', `zExp', and `zSig0' are simply
| added together to form the most significant 32 bits of the result.  This
| means that any integer portion of `zSig0' will be added into the exponent.
| Since a properly normalized significand will have an integer portion equal
| to 1, the `zExp' input should be 1 less than the desired result exponent
| whenever `zSig0' and `zSig1' concatenated form a complete, normalized
| significand.
*----------------------------------------------------------------------------*/

#ifndef ASSEMBLY_MACROS
static float128_t
 packFloat128( flag zSign, int32 zExp, bits64 zSig_high, bits64 zSig_low )
{
    float128_t z;

    z.low = zSig_low;
    z.high = ( ( (bits64) zSign )<<63 ) + ( ( (bits64) zExp )<<48 ) + zSig_high;
    return z;

}
#else
static inline float128_t __declspec(naked) packFloat128( flag zSign, int32 zExp, bits64 zSig0, bits64 zSig1 )
{
	_asm("\tshl\t$63,%rcx");
	_asm("\tshl\t$48,%rdx");
	_asm("\taddq\t%rcx,%r8");
	_asm("\tpinsrq\t$0,%r9,%xmm0");
	_asm("\taddq\t%rdx,%r8");
	_asm("\tpinsrq $1,%r8,%xmm0");
}
#endif

/*----------------------------------------------------------------------------
| Returns 1 if the quadruple-precision floating-point value `a' is a
| signaling NaN; otherwise returns 0.
*----------------------------------------------------------------------------*/
#if 1//ndef ASSEMBLY_MACROS
flag float128_is_signaling_nan( float128_t a )
{

    return
           ( ( ( a.high>>47 ) & 0xFFFF ) == 0xFFFE )
        && ( a.low || ( a.high & LIT64( 0x00007FFFFFFFFFFF ) ) );

}
#else
#define float128_is_signaling_nan(a) \
	( ( ( a.high>>47 ) & 0xFFFF ) == 0xFFFE ) \
        && ( a.low || ( a.high & LIT64( 0x00007FFFFFFFFFFF ) ) )

#endif

/*----------------------------------------------------------------------------
| Shifts `a' right by the number of bits given in `count'.  If any nonzero
| bits are shifted off, they are ``jammed'' into the least significant bit of
| the result by setting the least significant bit to 1.  The value of `count'
| can be arbitrarily large; in particular, if `count' is greater than 32, the
| result will be either 0 or 1, depending on whether `a' is zero or nonzero.
| The result is stored in the location pointed to by `zPtr'.
*----------------------------------------------------------------------------*/

static void shift32RightJamming( bits32 a, int count, bits32 *zPtr )
{
    bits32 z;

    if ( count == 0 ) {
        z = a;
    }
    else if ( count < 32 ) {
        z = ( a>>count ) | ( ( a<<( ( - count ) & 31 ) ) != 0 );
    }
    else {
        z = ( a != 0 );
    }
    *zPtr = z;

}
/*----------------------------------------------------------------------------
| Shifts the 128-bit value formed by concatenating `a0' and `a1' left by the
| number of bits given in `count'.  Any bits shifted off are lost.  The value
| of `count' must be less than 64.  The result is broken into two 64-bit
| pieces which are stored at the locations pointed to by `z0Ptr' and `z1Ptr'.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static float128_t  shortShift128Left(int count, bits64 a0, bits64 a1)
{
	float128_t r;
    r.low = a1<<count;
    r.high =
        ( count == 0 ) ? a0 : ( a0<<count ) | ( a1>>( ( - count ) & 63 ) );
	return r;
}
#else
static inline float128_t  __declspec(naked) shortShift128Left(int count, bits64 a0, bits64 a1)
{
	_asm("\tmovq\t%r8,%rax");
	_asm("\tsal\t%cl,%rax");
	_asm("\tpinsrq\t$0,%rax,%xmm0");
	_asm("\tsal\t%cl,%rdx");
	_asm("\tnegb\t%cl");
	_asm("\tandb\t$63,%cl");
	_asm("\tshrq\t%cl,%r8");
	_asm("\torq\t%rdx,%r8");
	_asm("\tpinsrq	$1,%r8,%xmm0");
}
#endif

#ifndef ASSEMBLY_MACROS
static float128_t  shortShift128Left16(bits64 a0, bits64 a1)
{
	float128_t r;
    r.low = a1<<16;
    r.high = ( a0<<16 ) | ( a1>>( ( -16 ) & 63 ) );
	return r;
}

#else
static inline float128_t  __declspec(naked) shortShift128Left16(bits64 a0, bits64 a1)
{
	_asm("\tmovq\t%rdx,%rax");
	_asm("\tsal\t$16,%rdx");
	_asm("\tpinsrq\t$0,%rdx,%xmm0");
	_asm("\tshr\t$48,%rax");
	_asm("\tsal\t$16,%rcx");
	_asm("\torq\t%rax,%rcx");
	_asm("\tpinsrq\t$1,%rcx,%xmm0");
}

#endif

/*----------------------------------------------------------------------------
| Returns the result of converting the quadruple-precision floating-point NaN
| `a' to the canonical NaN format.  If `a' is a signaling NaN, the invalid
| exception is raised.
*----------------------------------------------------------------------------*/
static commonNaNT float128ToCommonNaN( float128_t a )
{
    commonNaNT z;
	float128_t tmp;

    if ( float128_is_signaling_nan( a ) ) float_raise( float_flag_invalid );
    z.sign = a.high>>63;
    tmp = shortShift128Left16(a.high, a.low);
	z.high = tmp.high, z.low = tmp.low;
    return z;

}

/*----------------------------------------------------------------------------
| Shifts `a' right by the number of bits given in `count'.  If any nonzero
| bits are shifted off, they are ``jammed'' into the least significant bit of
| the result by setting the least significant bit to 1.  The value of `count'
| can be arbitrarily large; in particular, if `count' is greater than 64, the
| result will be either 0 or 1, depending on whether `a' is zero or nonzero.
| The result is stored in the location pointed to by `zPtr'.
*----------------------------------------------------------------------------*/

static void shift64RightJamming( bits64 a, int count, bits64 *zPtr )
{
    bits64 z;

    if ( count == 0 ) {
        z = a;
    }
    else if ( count < 64 ) {
        z = ( a>>count ) | ( ( a<<( ( - count ) & 63 ) ) != 0 );
    }
    else {
        z = ( a != 0 );
    }
    *zPtr = z;

}

/*----------------------------------------------------------------------------
| Shifts the 128-bit value formed by concatenating `a0' and `a1' right by 64
| _plus_ the number of bits given in `count'.  The shifted result is at most
| 64 nonzero bits; this is stored at the location pointed to by `z0Ptr'.  The
| bits shifted off form a second 64-bit result as follows:  The _last_ bit
| shifted off is the most-significant bit of the extra result, and the other
| 63 bits of the extra result are all zero if and only if _all_but_the_last_
| bits shifted off were all zero.  This extra result is stored in the location
| pointed to by `z1Ptr'.  The value of `count' can be arbitrarily large.
|     (This routine makes more sense if `a0' and `a1' are considered to form
| a fixed-point value with binary point between `a0' and `a1'.  This fixed-
| point value is shifted right by the number of bits given in `count', and
| the integer part of the result is returned at the location pointed to by
| `z0Ptr'.  The fractional part of the result may be slightly corrupted as
| described above, and is returned at the location pointed to by `z1Ptr'.)
*----------------------------------------------------------------------------*/

static void shift64ExtraRightJamming(bits64 a0, bits64 a1, int count, bits64 *z0Ptr, bits64 *z1Ptr )
{
    bits64 z0, z1;
    int8 negCount = ( - count ) & 63;

    if ( count == 0 ) {
        z1 = a1;
        z0 = a0;
    }
    else if ( count < 64 ) {
        z1 = ( a0<<negCount ) | ( a1 != 0 );
        z0 = a0>>count;
    }
    else {
        if ( count == 64 ) {
            z1 = a0 | ( a1 != 0 );
        }
        else {
            z1 = ( ( a0 | a1 ) != 0 );
        }
        z0 = 0;
    }
    *z1Ptr = z1;
    *z0Ptr = z0;

}

/*----------------------------------------------------------------------------
| Shifts the 128-bit value formed by concatenating `a0' and `a1' right by the
| number of bits given in `count'.  Any bits shifted off are lost.  The value
| of `count' can be arbitrarily large; in particular, if `count' is greater
| than 128, the result will be 0.  The result is broken into two 64-bit pieces
| which are stored at the locations pointed to by `z0Ptr' and `z1Ptr'.
*----------------------------------------------------------------------------*/
static float128_t shift128Right(bits64 a0, bits64 a1, int count)
{
    float128_t r;
    int8 negCount = ( - count ) & 63;

    if ( count == 0 ) {
        r.low = a1;
        r.high = a0;
    }
    else if ( count < 64 ) {
        r.low = ( a0<<negCount ) | ( a1>>count );
        r.high = a0>>count;
    }
    else {
        ///// BUG!!!!!z1 = ( count < 64 ) ? ( a0>>( count & 63 ) ) : 0;
        r.low = ( count < 128 ) ? ( a0>>( count & 63 ) ) : 0;
        r.high = 0;
    }
    return r;

}
#ifndef ASSEMBLY_MACROS
static void shift128RightSmallCount(bits64 a0, bits64 a1, int count, float128_t *tmp)
{
    int8 negCount = ( - count ) & 63;

    tmp->low = ( a0<<negCount ) | ( a1>>count );
    tmp->high = a0>>count;

}
#else
static inline __declspec(naked) void shift128RightSmallCount(bits64 a0, bits64 a1, int count, float128_t *tmp)
{
    _asm("\tmovq\t%rcx,%r11");
    _asm("\tmovl\t%r8d,%ecx");
    _asm("\tnegl\t%ecx");
    _asm("\tandb\t$63,%cl");
    _asm("\tmovq\t%r11,%rax");
    _asm("\tsal\t%cl,%rax");
    _asm("\tmovb\t%r8b,%cl");
    _asm("\tshrq\t%cl,%rdx");
    _asm("\torq\t%rdx,%rax");
    _asm("\tshrq\t%cl,%r11");
    _asm("\tmovq\t%rax,(%r9)");
    _asm("\tmovq\t%r11,8(%r9)");
}
#endif
/*----------------------------------------------------------------------------
| Shifts the 128-bit value formed by concatenating `a0' and `a1' right by the
| number of bits given in `count'.  If any nonzero bits are shifted off, they
| are ``jammed'' into the least significant bit of the result by setting the
| least significant bit to 1.  The value of `count' can be arbitrarily large;
| in particular, if `count' is greater than 128, the result will be either
| 0 or 1, depending on whether the concatenation of `a0' and `a1' is zero or
| nonzero.  The result is broken into two 64-bit pieces which are stored at
| the locations pointed to by `z0Ptr' and `z1Ptr'.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
void shift128RightJamming(bits64 a0, bits64 a1, int count, bits64 *z0Ptr, bits64 *z1Ptr )
{
    bits64 z0, z1;
    int8 negCount = ( - count ) & 63;

    if ( count == 0 ) {
        z1 = a1;
        z0 = a0;
    }
    else if ( count < 64 ) {
        z1 = ( a0<<negCount ) | ( a1>>count ) | ( ( a1<<negCount ) != 0 );
        z0 = a0>>count;
    }
    else {
        if ( count == 64 ) {
            z1 = a0 | ( a1 != 0 );
        }
        else if ( count < 128 ) {
            z1 = ( a0>>( count & 63 ) ) | ( ( ( a0<<negCount ) | a1 ) != 0 );
        }
        else {
            z1 = ( ( a0 | a1 ) != 0 );
        }
        z0 = 0;
    }
    *z1Ptr = z1;
    *z0Ptr = z0;

}
#else
void shift128RightJamming( bits64 a0, bits64 a1, int count, bits64 *z0Ptr, bits64 *z1Ptr );
#endif
#ifndef ASSEMBLY_MACROS
/*----------------------------------------------------------------------------
| Shifts the 192-bit value formed by concatenating `a0', `a1', and `a2' right
| by 64 _plus_ the number of bits given in `count'.  The shifted result is
| at most 128 nonzero bits; these are broken into two 64-bit pieces which are
| stored at the locations pointed to by `z0Ptr' and `z1Ptr'.  The bits shifted
| off form a third 64-bit result as follows:  The _last_ bit shifted off is
| the most-significant bit of the extra result, and the other 63 bits of the
| extra result are all zero if and only if _all_but_the_last_ bits shifted off
| were all zero.  This extra result is stored in the location pointed to by
| `z2Ptr'.  The value of `count' can be arbitrarily large.
|     (This routine makes more sense if `a0', `a1', and `a2' are considered
| to form a fixed-point value with binary point between `a1' and `a2'.  This
| fixed-point value is shifted right by the number of bits given in `count',
| and the integer part of the result is returned at the locations pointed to
| by `z0Ptr' and `z1Ptr'.  The fractional part of the result may be slightly
| corrupted as described above, and is returned at the location pointed to by
| `z2Ptr'.)
*----------------------------------------------------------------------------*/
void shift128ExtraRightJamming( bits64 a0, bits64 a1, bits64 a2,
     int count, bits64 *z0Ptr, bits64 *z1Ptr, bits64 *z2Ptr )
{
    bits64 z0, z1, z2;
    int8 negCount = ( - count ) & 63;

    if ( count == 0 ) {
        z2 = a2;
        z1 = a1;
        z0 = a0;
    }
    else {
        if ( count < 64 ) {
            z2 = a1<<negCount;
            z1 = ( a0<<negCount ) | ( a1>>count );
            z0 = a0>>count;
        }
        else {
            if ( count == 64 ) {
                z2 = a1;
                z1 = a0;
            }
            else {
                a2 |= a1;
                if ( count < 128 ) {
                    z2 = a0<<negCount;
                    z1 = a0>>( count & 63 );
                }
                else {
                    z2 = ( count == 128 ) ? a0 : ( a0 != 0 );
                    z1 = 0;
                }
            }
            z0 = 0;
        }
        z2 |= ( a2 != 0 );
    }
    *z2Ptr = z2;
    *z1Ptr = z1;
    *z0Ptr = z0;

}
#else

void shift128ExtraRightJamming( bits64 a0, bits64 a1, bits64 a2,
     int count, bits64 *z0Ptr, bits64 *z1Ptr, bits64 *z2Ptr );
#endif
#ifndef ASSEMBLY_MACROS
static float128_t shift128ExtraRightJammingByOne( bits64 a0, bits64 a1, bits64 a2, bits64 *z2Ptr )
{
    float128_t result;

    *z2Ptr = (a1<<63)|(a2 != 0);
    result.low = ( a0<<63 ) | ( a1>>1);
    result.high = a0>>1;
    return result;
}
#else
static inline __declspec(naked) float128_t shift128ExtraRightJammingByOne(
                          bits64 a0, bits64 a1, bits64 a2, bits64 *z2Ptr )
{
	_asm("\torq\t%r8,%r8");
	_asm("\tsetne\t%r8b");
	_asm("\tandl\t$1,%r8d");
	_asm("\tmovq\t%rdx,%r11");
	_asm("\tsal\t$63,%rdx");
	_asm("\torq\t%r8,%rdx");
	_asm("\tmovq\t%rdx,(%r9)");
// result.low = ( a0<<63 ) | ( a1>>1);
	_asm("\tmovq\t%rcx,%r10");
	_asm("\tsal\t$63,%r10");
	_asm("\tshrq\t$1,%r11");
	_asm("\torq\t%r11,%r10");
	_asm("\tpinsrq\t$0,%r10,%xmm0");
// result.high = a0>>1
	_asm("\tshrq\t$1,%rcx");
	_asm("\tpinsrq\t$1,%rcx,%xmm0");
}
#endif
#ifndef ASSEMBLY_MACROS
/*
This is the same procedure as shift128ExtraRightJamming but specialized for a2 zero
and count bigger than zero and smaller than 64. This is the case for two important
calls to the above procedure in the divide and square root operation, where they are
in the critical path.
*/
static float128_t shift128ExtraRightJammingSmallCount(
     int count,
     bits64 a0,
     bits64 a1,
     bits64 *z2Ptr
 )
{
    float128_t z;
    int negCount = ( - count ) & 63;

    *z2Ptr = a1<<negCount;
    z.low = ( a0<<negCount ) | ( a1>>count );
    z.high = a0>>count;
	return z;
}
#else
static inline float128_t __declspec(naked) shift128ExtraRightJammingSmallCount(
     int count, bits64 a0, bits64 a1, bits64 *z2Ptr)
{
	_asm("\tmovq\t%rcx,%r10");
	_asm("\tnegl\t%ecx");
	_asm("\tandl\t$63,%ecx");

	_asm("\tmovq\t%r8,%rax");
	_asm("\tsal\t%cl,%rax");
	_asm("\tmovq\t%rax,(%r9)");

	_asm("\tmovq\t%rdx,%rax");
	_asm("\tsal\t%cl,%rax");
	_asm("\tmovl\t%r10d,%ecx");
	_asm("\tshrq\t%cl,%r8");
	_asm("\torq\t%r8,%rax");
	_asm("\tpinsrq\t$0,%rax,%xmm0");
	_asm("\tshrq\t%cl,%rdx");
	_asm("\tpinsrq\t$1,%rdx,%xmm0");
}
#endif

/*----------------------------------------------------------------------------
| Shifts the 192-bit value formed by concatenating `a0', `a1', and `a2' left
| by the number of bits given in `count'.  Any bits shifted off are lost.
| The value of `count' must be less than 64.  The result is broken into three
| 64-bit pieces which are stored at the locations pointed to by `z0Ptr',
| `z1Ptr', and `z2Ptr'.
*----------------------------------------------------------------------------*/

static void shortShift192Left( bits64 a0, bits64 a1, bits64 a2, int count,
     bits64 *z0Ptr, bits64 *z1Ptr, bits64 *z2Ptr )
{
    bits64 z0, z1, z2;
    int8 negCount;

    z2 = a2<<count;
    z1 = a1<<count;
    z0 = a0<<count;
    if ( 0 < count ) {
        negCount = ( ( - count ) & 63 );
        z1 |= a2>>negCount;
        z0 |= a1>>negCount;
    }
    *z2Ptr = z2;
    *z1Ptr = z1;
    *z0Ptr = z0;

}

/*----------------------------------------------------------------------------
| Adds the 128-bit value formed by concatenating `a0' and `a1' to the 128-bit
| value formed by concatenating `b0' and `b1'.  Addition is modulo 2^128, so
| any carry out is lost.  The result is broken into two 64-bit pieces which
| are stored at the locations pointed to by `z0Ptr' and `z1Ptr'.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static float128_t add128( bits64 a0, bits64 a1, bits64 b0, bits64 b1)
{
    float128_t z;
    z.low = a1 + b1;
    z.high = a0 + b0 + ( z.low < a1 );
    return z;
}
#else
static inline __declspec(naked) float128_t
	add128(bits64 a0, bits64 a1, bits64 b0, bits64 b1)
{
	_asm("\taddq\t%r9,%rdx");
	_asm("\tpinsrq\t$0,%rdx,%xmm0");
	_asm("\tadcq\t%r8,%rcx");
	_asm("\tpinsrq\t$1,%rcx,%xmm0");
}
#endif

/*----------------------------------------------------------------------------
| Adds the 192-bit value formed by concatenating `a0', `a1', and `a2' to the
| 192-bit value formed by concatenating `b0', `b1', and `b2'.  Addition is
| modulo 2^192, so any carry out is lost.  The result is broken into three
| 64-bit pieces which are stored at the locations pointed to by `z0Ptr',
| `z1Ptr', and `z2Ptr'.
| It is used in divide and square root. In all the usages, the argument 'b0'
| is zero. Original code modified to assume this, eliminating one argument
| j.n.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static void add192( bits64 a0, bits64 a1, bits64 a2, bits64 b1, bits64 b2,
     bits64 *z0Ptr,bits64 *z1Ptr,bits64 *z2Ptr)
{
    bits64 z0, z1;
    int8 carry1;

    z0 = a2 + b2;
    carry1 = ( z0 < a2 );
    *z2Ptr = z0;
    z1 = a1 + b1;
    z0 = a0 + (z1 < a1);
    z1 += carry1;
    z0 += ( z1 < carry1 );
    *z1Ptr = z1;
    *z0Ptr = z0;
}
#else
static inline void __declspec(naked) add192( bits64 a0, bits64 a1, bits64 a2, bits64 b1,
bits64 b2, bits64 *z0Ptr,bits64 *z1Ptr,bits64 *z2Ptr)
{
//          z0 = a2 + b2;           z0 --> r11
	_asm("\tmovq	%r8,%r11");
	_asm("\taddq	32(%rsp),%r11");
//          carry1 = ( z0 < a2 );    carry1 --> rax
//	_asm("\tcmpq	%r8,%r11");
	_asm("\tsetc	%al");
//      	*z2Ptr = z0;
	_asm("\tmovq	56(%rsp),%r10");
	_asm("\tmovq	%r11,(%r10)");
//          z1 = a1 + b1;          z1 --> r8
	_asm("\tmovq	%rdx,%r8");
	_asm("\taddq	%r9,%r8");
//          z0 = a0 + (z1 < a1);  z0 --> r10
	_asm("\tcmpq	%rdx,%r8");
	_asm("\tsetb	%r10b");
	_asm("\tmovq	48(%rsp),%r11");
	_asm("\tandl	$1,%r10d");
	_asm("\taddq	%rcx,%r10");
//          z1 += carry1;
	_asm("\tmovsbq	%al,%r9");
	_asm("\taddq	%r9,%r8");
//          z0 += ( z1 < carry1 );
	_asm("\tcmpq	%r9,%r8");
	_asm("\tsetb	%cl");
	_asm("\tmovq	40(%rsp),%rdx");
	_asm("\tandl	$1,%ecx");
	_asm("\taddq	%rcx,%r10");
//          *z1Ptr = z1;
	_asm("\tmovq	%r8,(%r11)");
//          *z0Ptr = z0;
	_asm("\tmovq	%r10,(%rdx)");
}
#endif
/*----------------------------------------------------------------------------
| Subtracts the 128-bit value formed by concatenating `b0' and `b1' from the
| 128-bit value formed by concatenating `a0' and `a1'.  Subtraction is modulo
| 2^128, so any borrow out (carry out) is lost.  The result is broken into two
| 64-bit pieces which are stored at the locations pointed to by `z0Ptr' and
| `z1Ptr'.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static float128_t sub128(bits64 a0, bits64 a1, bits64 b0, bits64 b1)
{
	float128_t result;
    result.low = a1 - b1;
    result.high = a0 - b0 - ( a1 < b1 );
	return result;
}
#else
static inline __declspec(naked) float128_t sub128(bits64 a0, bits64 a1, bits64 b0, bits64 b1)
{
	_asm("\tsubq\t%r9,%rdx");
	_asm("\tpinsrq\t$0,%rdx,%xmm0");
	_asm("\tsbb\t%r8,%rcx");
	_asm("\tpinsrq\t$1,%rcx,%xmm0");
}
#endif

/*----------------------------------------------------------------------------
| Subtracts the 192-bit value formed by concatenating `b0', `b1', and `b2'
| from the 192-bit value formed by concatenating `a0', `a1', and `a2'.
| Subtraction is modulo 2^192, so any borrow out (carry out) is lost.  The
| result is broken into three 64-bit pieces which are stored at the locations
| pointed to by `z0Ptr', `z1Ptr', and `z2Ptr'.
| In all calls the third argument (a2) is always zero. Eliminated
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static void sub192( bits64 a0, bits64 a1, bits64 b0, bits64 b1,
     bits64 b2, bits64 *z0Ptr, bits64 *z1Ptr, bits64 *z2Ptr )
{
    bits64 z0, z1, z2;
    int8 borrow0, borrow1;

    z2 = - (long long)b2;
    borrow1 = ( 0 < b2 );
    z1 = a1 - b1;
    borrow0 = ( a1 < b1 );
    z0 = a0 - b0;
    z0 -= ( z1 < borrow1 );
    z1 -= borrow1;
    z0 -= borrow0;
    *z2Ptr = z2;
    *z1Ptr = z1;
    *z0Ptr = z0;
}
#else
static void inline __declspec(naked) sub192( bits64 a0, bits64 a1, bits64 b0, bits64 b1,
     bits64 b2, bits64 *z0Ptr, bits64 *z1Ptr, bits64 *z2Ptr )
{
	_asm("\tmovq\t32(%rsp),%r10");
	_asm("\tmovq\t%r10,%r11");
	_asm("\tmovq\t56(%rsp),%rax");
	_asm("\tnegq\t%r10");

	_asm("\tmovq\t%r10,(%rax)");

	_asm("\torq\t%r11,%r11");
	_asm("\tseta\t%al");

	_asm("\tmovq\t%r8,%r10");
	_asm("\tmovq\t%rdx,%r8");
	_asm("\tsubq\t%r9,%r8");

	_asm("\tcmpq\t%r9,%rdx");
	_asm("\tsetb\t%r11b");

	_asm("\tsubq\t%r10,%rcx");

	_asm("\tmovsbq\t%al,%r10");
	_asm("\tcmpq\t%r10,%r8");
	_asm("\tsetb\t%dl");
	_asm("\tandl\t$1,%edx");
	_asm("\tsubq\t%rdx,%rcx");

	_asm("\tmovsbq\t%al,%r10");
	_asm("\tsubq\t%r10,%r8");

	_asm("\tmovsbq\t%r11b,%r10");
	_asm("\tsubq\t%r10,%rcx");

	_asm("\tmovq\t48(%rsp),%r10");
	_asm("\tmovq\t%r8,(%r10)");
	_asm("\tmovq\t40(%rsp),%r11");
	_asm("\tmovq\t%rcx,(%r11)");
}
#endif
/*----------------------------------------------------------------------------
| Multiplies `a' by `b' to obtain a 128-bit product.  The product is broken
| into two 64-bit pieces which are stored at the locations pointed to by
| `z0Ptr' and `z1Ptr'.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static void mul64To128( bits64 a, bits64 b, bits64 *z0Ptr, bits64 *z1Ptr )
{
    bits32 aHigh, aLow, bHigh, bLow;
    bits64 z0, zMiddleA, zMiddleB, z1;

    aLow = (bits32)a;
    aHigh = a>>32;
    bLow = (bits32)b;
    bHigh = b>>32;
    z1 = ( (bits64) aLow ) * bLow;
    zMiddleA = ( (bits64) aLow ) * bHigh;
    zMiddleB = ( (bits64) aHigh ) * bLow;
    z0 = ( (bits64) aHigh ) * bHigh;
    zMiddleA += zMiddleB;
    z0 += ( ( (bits64) ( zMiddleA < zMiddleB ) )<<32 ) + ( zMiddleA>>32 );
    zMiddleA <<= 32;
    z1 += zMiddleA;
    z0 += ( z1 < zMiddleA );
    *z1Ptr = z1;
    *z0Ptr = z0;

}
#else
static void inline __declspec(naked) mul64To128( bits64 a, bits64 b, bits64 *z0Ptr, bits64 *z1Ptr )
{
	_asm("\tmovq\t%rcx,%rax");
	_asm("\tmulq\t%rdx");
	_asm("\tmovq\t%rdx,(%r8)");
	_asm("\tmovq\t%rax,(%r9)");
}
#endif

/*----------------------------------------------------------------------------
| Multiplies the 128-bit value formed by concatenating `a0' and `a1' by
| `b' to obtain a 192-bit product.  The product is broken into three 64-bit
| pieces which are stored at the locations pointed to by `z0Ptr', `z1Ptr', and
| `z2Ptr'.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static float128_t mul128By64To192( bits64 a0, bits64 a1, bits64 b, bits64 *z2Ptr )
{
    bits64 z0, z1, more1;

    mul64To128( a1, b, &z1, z2Ptr );
    mul64To128( a0, b, &z0, &more1 );
    return add128( z0, more1, 0, z1);
}
#else

extern inline __declspec(naked) float128_t mul128By64To192(bits64 a0,bits64 a1,bits64 b,bits64 *z2Ptr)
{
	// mul64To128( a1, b, &z1, &z2 );
	_asm("\tmovq\t%rdx,%rax");
	_asm("\tmulq\t%r8");
	_asm("\tmovq\t%rdx,%r11");
	_asm("\tmovq\t%rax,(%r9)");
	// mul64To128( a0, b, &z0, &more1 );
	_asm("\tmovq\t%rcx,%rax");
	_asm("\tmulq\t%r8");
	// add128( z0, more1, 0, z1);
	_asm("\taddq\t%r11,%rax");
	_asm("\tadcq\t$0,%rdx");
	_asm("\tpinsrq\t$0,%rax,%xmm0");
	_asm("\tpinsrq\t$1,%rdx,%xmm0");
}
#endif
#ifndef ASSEMBLY_MACROS
/*----------------------------------------------------------------------------
| Multiplies the 128-bit value formed by concatenating `a0' and `a1' to the
| 128-bit value formed by concatenating `b0' and `b1' to obtain a 256-bit
| product.  The product is broken into four 64-bit pieces which are stored at
| the locations pointed to by `z0Ptr', `z1Ptr', `z2Ptr', and `z3Ptr'.
*----------------------------------------------------------------------------*/

static void mul128To256( bits64 a0, bits64 a1, bits64 b0, bits64 b1, bits64 *z0Ptr,
     bits64 *z1Ptr, bits64 *z2Ptr, bits64 *z3Ptr )
{
    bits64 z0, z1, z2, z3;
    bits64 more1, more2;
    float128_t z, y, x, m;

    mul64To128( a1, b1, &z2, &z3 );
    mul64To128( a1, b0, &z1, &more2 );

    z = add128( z1, more2, 0, z2);

    mul64To128( a0, b0, &z0, &more1 );
    y = add128( z0, more1, 0, z.high );

    mul64To128( a0, b1, &more1, &more2 );

    x = add128( more1, more2, 0, z.low );

    m = add128( y.high, y.low, 0, x.high);

    *z3Ptr = z3;
    *z2Ptr = x.low;
    *z1Ptr = m.low;
    *z0Ptr = m.high;

}
#else
extern void __mul128To256(bits64 a0, bits64 a1, bits64 b0, bits64 b1, bits64 *z0Ptr, bits64 *z1Ptr, bits64 *z2Ptr, bits64 *z3Ptr );
#define mul128To256 __mul128To256
#endif
/*----------------------------------------------------------------------------
| Returns an approximation to the 64-bit integer quotient obtained by dividing
| `b' into the 128-bit value formed by concatenating `a0' and `a1'.  The
| divisor `b' must be at least 2^63.  If q is the exact quotient truncated
| toward zero, the approximation returned lies between q and q + 2 inclusive.
| If the exact quotient q is larger than 64 bits, the maximum positive 64-bit
| unsigned integer is returned.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static bits64 estimateDiv128To64( bits64 a0, bits64 a1, bits64 b )
{
    bits64 b0, b1;
    bits64 rem0, rem1, term0, term1;
    bits64 z;
    float128_t y;

    if ( b <= a0 ) return LIT64( 0xFFFFFFFFFFFFFFFF );
    b0 = b>>32;
    z = ( b0<<32 <= a0 ) ? LIT64( 0xFFFFFFFF00000000 ) : ( a0 / b0 )<<32;
    mul64To128( b, z, &term0, &term1 );
    y = sub128( a0, a1, term0, term1);
	rem0 = y.high, rem1 = y.low;
    while ( ( (sbits64) rem0 ) < 0 ) {
        z -= LIT64( 0x100000000 );
        b1 = b<<32;
        y = add128( rem0, rem1, b0, b1);
        rem0 = y.high, rem1 = y.low;
    }
    rem0 = ( rem0<<32 ) | ( rem1>>32 );
    z |= ( b0<<32 <= rem0 ) ? 0xFFFFFFFF : rem0 / b0;
    return z;

}
#else
static inline bits64 __declspec(naked) estimateDiv128To64( bits64 a0, bits64 a1, bits64 b )
{
	_asm("\tmovq\t%rdx,%rax");
	_asm("\tmovq\t%rcx,%rdx");
	_asm("\tdivq\t%r8");
}
#endif

/*----------------------------------------------------------------------------
| Returns an approximation to the square root of the 32-bit significand given
| by `a'.  Considered as an integer, `a' must be at least 2^31.  If bit 0 of
| `aExp' (the least significant bit) is 1, the integer returned approximates
| 2^31*sqrt(`a'/2^31), where `a' is considered an integer.  If bit 0 of `aExp'
| is 0, the integer returned approximates 2^31*sqrt(`a'/2^30).  In either
| case, the approximation returned lies strictly within +/-2 of the exact
| value.
*----------------------------------------------------------------------------*/

static bits32 estimateSqrt32( int aExp, bits32 a )
{
    static const bits16 sqrtOddAdjustments[] = {
        0x0004, 0x0022, 0x005D, 0x00B1, 0x011D, 0x019F, 0x0236, 0x02E0,
        0x039C, 0x0468, 0x0545, 0x0631, 0x072B, 0x0832, 0x0946, 0x0A67
    };
    static const bits16 sqrtEvenAdjustments[] = {
        0x0A2D, 0x08AF, 0x075A, 0x0629, 0x051A, 0x0429, 0x0356, 0x029E,
        0x0200, 0x0179, 0x0109, 0x00AF, 0x0068, 0x0034, 0x0012, 0x0002
    };
    unsigned char index;
    bits32 z;

    index = ( a>>27 ) & 15;
    if ( aExp & 1 ) {
        z = 0x4000 + ( a>>17 ) - sqrtOddAdjustments[ index ];
        z = ( ( a / z )<<14 ) + ( z<<15 );
        a >>= 1;
    }
    else {
        z = 0x8000 + ( a>>17 ) - sqrtEvenAdjustments[ index ];
        z = a / z + z;
        z = ( 0x20000 <= z ) ? 0xFFFF8000 : ( z<<15 );
        if ( z <= a ) return (bits32) ( ( (sbits32) a )>>1 );
    }
    return ( (bits32) ( ( ( (bits64) a )<<31 ) / z ) ) + ( z>>1 );

}

/*----------------------------------------------------------------------------
| Returns the number of leading 0 bits before the most-significant 1 bit of
| `a'.  If `a' is zero, 32 is returned.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static int8 countLeadingZeros32( bits32 a )
{
    static const int8 countLeadingZerosHigh[] = {
        8, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    int8 shiftCount;

    shiftCount = 0;
    if ( a < 0x10000 ) {
        shiftCount += 16;
        a <<= 16;
    }
    if ( a < 0x1000000 ) {
        shiftCount += 8;
        a <<= 8;
    }
    shiftCount += countLeadingZerosHigh[ a>>24 ];
    return shiftCount;
}
#else
static inline int __declspec(naked) countLeadingZeros32( bits32 a )
{
	_asm("\tmovl\t$32,%eax");
	_asm("\txorl\t%edx,%edx");
	_asm("\tbsrl\t%ecx,%ecx");
	_asm("\tcmovne\t%ecx,%edx");
	_asm("\tsetne\t%cl");
	_asm("\tandl\t$1,%ecx");
	_asm("\tsubl\t%edx,%eax");
	_asm("\tsubl\t%ecx,%eax");
}
#endif
#ifndef ASSEMBLY_MACROS
static int countLeadingZeros32NoChecks(bits32 a)
{
	return countLeadingZeros32(a);
}
#else
static inline int __declspec(naked) countLeadingZeros32NoChecks(bits32 a)
{
	_asm("\tmovl\t$31,%eax");
	_asm("\tbsrl\t%ecx,%edx");
	_asm("\tsubl\t%edx,%eax");
}
#endif
/*----------------------------------------------------------------------------
| Returns the number of leading 0 bits before the most-significant 1 bit of
| `a'.  If `a' is zero, 64 is returned.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static int8 countLeadingZeros64( bits64 a )
{
    int8 shiftCount;

    shiftCount = 0;
    if ( a < ( (bits64) 1 )<<32 ) {
        shiftCount += 32;
    }
    else {
        a >>= 32;
    }
    shiftCount += countLeadingZeros32( a );
    return shiftCount;

}
#else
static int inline __declspec(naked) countLeadingZeros64( bits64 a )
{
	_asm("\tbsrq\t%rcx,%rdx");
	_asm("\tmovl\t$63,%eax");
	_asm("\tsubl\t%edx,%eax");
}
#endif

/*----------------------------------------------------------------------------
| Returns 1 if the 128-bit value formed by concatenating `a0' and `a1'
| is equal to the 128-bit value formed by concatenating `b0' and `b1'.
| Otherwise, returns 0.
*----------------------------------------------------------------------------*/

#ifndef ASSEMBLY_MACROS
static flag eq128( bits64 a0, bits64 a1, bits64 b0, bits64 b1 )
{

    return ( a0 == b0 ) && ( a1 == b1 );

}
#else
static int inline __declspec(naked) eq128( bits64 a0, bits64 a1, bits64 b0, bits64 b1 )
{
	_asm("\txorl\t%eax,%eax");
	_asm("\tcmpq\t%r8,%rcx");
	_asm("\tsete\t%cl");
	_asm("\tcmpq	%r9,%rdx");
	_asm("\tsete\t%al");
	_asm("\tandb\t%cl,%al");
}
#endif

#ifndef ASSEMBLY_MACROS
static float128_t packFloat128Zero(void)
{
	return packFloat128(0,0,0,0);
}
#else
static inline float128_t __declspec(naked) packFloat128Zero(void)
{
	_asm("\tpxor\t%xmm0,%xmm0");
}
#endif

/*----------------------------------------------------------------------------
| Returns 1 if the 128-bit value formed by concatenating `a0' and `a1' is less
| than or equal to the 128-bit value formed by concatenating `b0' and `b1'.
| Otherwise, returns 0.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static int le128( bits64 a0, bits64 a1, bits64 b0, bits64 b1 )
{

    return ( a0 < b0 ) || ( ( a0 == b0 ) && ( a1 <= b1 ) );

}
#else
static int __declspec(naked) inline le128( bits64 a0, bits64 a1, bits64 b0, bits64 b1 )
{
	// a0=rcx,a1=rdx,b0=r8,b1=r9
	_asm("\txorl\t%eax,%eax");
	_asm("\tcmpq\t%r8,%rcx");
	_asm("\tsetb\t%al");
	_asm("\tsete\t%cl");
	_asm("\tcmpq\t%r9,%rdx");
	_asm("\tsetbe\t%dl");
	_asm("\tandb\t%dl,%cl");
	_asm("\torb\t%cl,%al");
}
#endif

/*----------------------------------------------------------------------------
| Returns 1 if the 128-bit value formed by concatenating `a0' and `a1' is less
| than the 128-bit value formed by concatenating `b0' and `b1'.  Otherwise,
| returns 0.
*----------------------------------------------------------------------------*/

flag lt128( bits64 a0, bits64 a1, bits64 b0, bits64 b1 )
{

    return ( a0 < b0 ) || ( ( a0 == b0 ) && ( a1 < b1 ) );

}
#ifndef ASSEMBLY_MACROS
flag ne128( bits64 a0, bits64 a1, bits64 b0, bits64 b1 )
{

    return ( a0 != b0 ) || ( a1 != b1 );

}
#else
int __declspec(naked) ne128( bits64 a0, bits64 a1, bits64 b0, bits64 b1 )
{
	_asm("\txorl\t%eax,%eax");
	_asm("\tcmpq	%r8,%rcx");
	_asm("\tsetne	%al");
	_asm("\tcmpq	%r9,%rcx");
	_asm("\tsetne	%cl");
	_asm("\torb	%cl,%al");
	_asm("\tret");
}
#endif



/*----------------------------------------------------------------------------
| Functions and definitions to determine:  (1) whether tininess for underflow
| is detected before or after rounding by default, (2) what (if anything)
| happens when exceptions are raised, (3) how signaling NaNs are distinguished
| from quiet NaNs, (4) the default generated quiet NaNs, and (5) how NaNs
| are propagated from function inputs to output.  These details are target-
| specific.
*----------------------------------------------------------------------------*/

/*============================================================================

This C source fragment is part of the SoftFloat IEC/IEEE Floating-point
Arithmetic Package, Release 2b.

Written by John R. Hauser.  This work was made possible in part by the
International Computer Science Institute, located at Suite 600, 1947 Center
Street, Berkeley, California 94704.  Funding was partially provided by the
National Science Foundation under grant MIP-9311980.  The original version
of this code was written as part of a project to build a fixed-point vector
processor in collaboration with the University of California at Berkeley,
overseen by Profs. Nelson Morgan and John Wawrzynek.  More information
is available through the Web page `http://www.cs.berkeley.edu/~jhauser/
arithmetic/SoftFloat.html'.

THIS SOFTWARE IS DISTRIBUTED AS IS, FOR FREE.  Although reasonable effort has
been made to avoid it, THIS SOFTWARE MAY CONTAIN FAULTS THAT WILL AT TIMES
RESULT IN INCORRECT BEHAVIOR.  USE OF THIS SOFTWARE IS RESTRICTED TO PERSONS
AND ORGANIZATIONS WHO CAN AND WILL TAKE FULL RESPONSIBILITY FOR ALL LOSSES,
COSTS, OR OTHER PROBLEMS THEY INCUR DUE TO THE SOFTWARE, AND WHO FURTHERMORE
EFFECTIVELY INDEMNIFY JOHN HAUSER AND THE INTERNATIONAL COMPUTER SCIENCE
INSTITUTE (possibly via similar legal warning) AGAINST ALL LOSSES, COSTS, OR
OTHER PROBLEMS INCURRED BY THEIR CUSTOMERS AND CLIENTS DUE TO THE SOFTWARE.

Derivative works are acceptable, even for commercial purposes, so long as
(1) the source code for the derivative work includes prominent notice that
the work is derivative, and (2) the source code includes prominent notice with
these four paragraphs for those parts of this code that are retained.

=============================================================================*/
/*----------------------------------------------------------------------------
| Returns the fraction bits of the extended double-precision floating-point
| value `a'.
*----------------------------------------------------------------------------*/

static bits64 extractFloatx80Frac( floatx80 a )
{

    return a.low;

}

/*----------------------------------------------------------------------------
| Returns the exponent bits of the extended double-precision floating-point
| value `a'.
*----------------------------------------------------------------------------*/

static int32 extractFloatx80Exp( floatx80 a )
{

    return a.high & 0x7FFF;

}

/*----------------------------------------------------------------------------
| Returns the sign bit of the extended double-precision floating-point value
| `a'.
*----------------------------------------------------------------------------*/

static flag extractFloatx80Sign( floatx80 a )
{

    return a.high>>15;

}

/*----------------------------------------------------------------------------
| Returns the fraction bits of the single-precision floating-point value `a'.
*----------------------------------------------------------------------------*/

static bits32 extractFloat32Frac( float32 a )
{

    return a & 0x007FFFFF;

}
/*----------------------------------------------------------------------------
| Returns the exponent bits of the single-precision floating-point value `a'.
*----------------------------------------------------------------------------*/

static int extractFloat32Exp( float32 a )
{

    return ( a>>23 ) & 0xFF;

}

/*----------------------------------------------------------------------------
| Returns the sign bit of the single-precision floating-point value `a'.
*----------------------------------------------------------------------------*/

static flag extractFloat32Sign( float32 a )
{

    return a>>31;

}
/*----------------------------------------------------------------------------
| Returns 1 if the single-precision floating-point value `a' is a signaling
| NaN; otherwise returns 0.
*----------------------------------------------------------------------------*/

flag float32_is_signaling_nan( float32 a )
{

    return ( ( ( a>>22 ) & 0x1FF ) == 0x1FE ) && ( a & 0x003FFFFF );

}

/*----------------------------------------------------------------------------
| Returns the result of converting the single-precision floating-point NaN
| `a' to the canonical NaN format.  If `a' is a signaling NaN, the invalid
| exception is raised.
*----------------------------------------------------------------------------*/

static commonNaNT float32ToCommonNaN( float32 a )
{
    commonNaNT z;

    if ( float32_is_signaling_nan( a ) ) float_raise( float_flag_invalid );
    z.sign = a>>31;
    z.low = 0;
    z.high = ( (bits64) a )<<41;
    return z;

}
/*----------------------------------------------------------------------------
| Returns the fraction bits of the double-precision floating-point value `a'.
*----------------------------------------------------------------------------*/

static bits64 extractFloat64Frac( float64 a )
{

    return a & LIT64( 0x000FFFFFFFFFFFFF );

}

/*----------------------------------------------------------------------------
| Returns the exponent bits of the double-precision floating-point value `a'.
*----------------------------------------------------------------------------*/

static int extractFloat64Exp( float64 a )
{

    return ( a>>52 ) & 0x7FF;

}

/*----------------------------------------------------------------------------
| Returns the sign bit of the double-precision floating-point value `a'.
*----------------------------------------------------------------------------*/

static flag extractFloat64Sign( float64 a )
{

    return a>>63;

}

/*----------------------------------------------------------------------------
| Normalizes the subnormal single-precision floating-point value represented
| by the denormalized significand `aSig'.  The normalized exponent and
| significand are stored at the locations pointed to by `zExpPtr' and
| `zSigPtr', respectively.
*----------------------------------------------------------------------------*/

static void
 normalizeFloat32Subnormal( bits32 aSig, int *zExpPtr, bits32 *zSigPtr )
{
    int8 shiftCount;

    shiftCount = countLeadingZeros32( aSig ) - 8;
    *zSigPtr = aSig<<shiftCount;
    *zExpPtr = 1 - shiftCount;

}
/*----------------------------------------------------------------------------
| Normalizes the subnormal double-precision floating-point value represented
| by the denormalized significand `aSig'.  The normalized exponent and
| significand are stored at the locations pointed to by `zExpPtr' and
| `zSigPtr', respectively.
*----------------------------------------------------------------------------*/

static void
 normalizeFloat64Subnormal( bits64 aSig, int *zExpPtr, bits64 *zSigPtr )
{
    int8 shiftCount;

    shiftCount = countLeadingZeros64( aSig ) - 11;
    *zSigPtr = aSig<<shiftCount;
    *zExpPtr = 1 - shiftCount;

}

/*----------------------------------------------------------------------------
| Underflow tininess-detection mode, statically initialized to default value.
| (The declaration in `softfloat.h' must match the `int8' type here.)
*----------------------------------------------------------------------------*/
int8 float_detect_tininess = float_tininess_after_rounding;

/*----------------------------------------------------------------------------
| Raises the exceptions specified by `flags'.  Floating-point traps can be
| defined here if desired.  It is currently not possible for such a trap
| to substitute a result value.  If traps are not implemented, this routine
| should be simply `float_exception_flags |= flags;'.
*----------------------------------------------------------------------------*/

void float_raise( int8 flags )
{

    float_exception_flags |= flags;

}
/*----------------------------------------------------------------------------
| Returns 1 if the double-precision floating-point value `a' is a signaling
| NaN; otherwise returns 0.
*----------------------------------------------------------------------------*/

flag float64_is_signaling_nan( float64 a )
{

    return
           ( ( ( a>>51 ) & 0xFFF ) == 0xFFE )
        && ( a & LIT64( 0x0007FFFFFFFFFFFF ) );

}
/*----------------------------------------------------------------------------
| Returns the result of converting the double-precision floating-point NaN
| `a' to the canonical NaN format.  If `a' is a signaling NaN, the invalid
| exception is raised.
*----------------------------------------------------------------------------*/

static commonNaNT float64ToCommonNaN( float64 a )
{
    commonNaNT z;

    if ( float64_is_signaling_nan( a ) ) float_raise( float_flag_invalid );
    z.sign = a>>63;
    z.low = 0;
    z.high = a<<12;
    return z;

}

/*----------------------------------------------------------------------------
| The pattern for a default generated single-precision NaN.
*----------------------------------------------------------------------------*/
#define float32_default_nan 0xFFC00000


/*----------------------------------------------------------------------------
| Returns 1 if the double-precision floating-point value `a' is a NaN;
| otherwise returns 0.
*----------------------------------------------------------------------------*/

flag float64_is_nan( float64 a )
{

    return ( LIT64( 0xFFE0000000000000 ) < (bits64) ( a<<1 ) );

}


/*----------------------------------------------------------------------------
| The pattern for a default generated quadruple-precision NaN.  The `high' and
| `low' values hold the most- and least-significant bits, respectively.
*----------------------------------------------------------------------------*/
#define float128_default_nan_high LIT64( 0xFFFF800000000000 )
#define float128_default_nan_low  LIT64( 0x0000000000000000 )

/*----------------------------------------------------------------------------
| Returns 1 if the quadruple-precision floating-point value `a' is a NaN;
| otherwise returns 0.
*----------------------------------------------------------------------------*/

flag isnanf128( float128_t a )
{

    return
           ( LIT64( 0xFFFE000000000000 ) <= (bits64) ( a.high<<1 ) )
        && ( a.low || ( a.high & LIT64( 0x0000FFFFFFFFFFFF ) ) );

}


/*----------------------------------------------------------------------------
| Returns the result of converting the canonical NaN `a' to the quadruple-
| precision floating-point format.
*----------------------------------------------------------------------------*/

static float128_t commonNaNToFloat128( commonNaNT a )
{
    float128_t z;

    shift128RightSmallCount( a.high, a.low, 16,&z);
    z.high |= ( ( (bits64) a.sign )<<63 ) | LIT64( 0x7FFF800000000000 );
    return z;

}

/*----------------------------------------------------------------------------
| Takes two quadruple-precision floating-point values `a' and `b', one of
| which is a NaN, and returns the appropriate NaN result.  If either `a' or
| `b' is a signaling NaN, the invalid exception is raised.
*----------------------------------------------------------------------------*/

float128_t propagateFloat128NaN( float128_t a, float128_t b )
{
    flag aIsNaN, aIsSignalingNaN, bIsNaN, bIsSignalingNaN;

    aIsNaN = isnanf128( a );
    aIsSignalingNaN = float128_is_signaling_nan( a );
    bIsNaN = isnanf128( b );
    bIsSignalingNaN = float128_is_signaling_nan( b );
    a.high |= LIT64( 0x0000800000000000 );
    b.high |= LIT64( 0x0000800000000000 );
    if ( aIsSignalingNaN | bIsSignalingNaN ) float_raise( float_flag_invalid );
    if ( aIsSignalingNaN ) {
        if ( bIsSignalingNaN ) goto returnLargerSignificand;
        return bIsNaN ? b : a;
    }
    else if ( aIsNaN ) {
        if ( bIsSignalingNaN | ! bIsNaN ) return a;
 returnLargerSignificand:
        if ( lt128( a.high<<1, a.low, b.high<<1, b.low ) ) return b;
        if ( lt128( b.high<<1, b.low, a.high<<1, a.low ) ) return a;
        return ( a.high < b.high ) ? a : b;
    }
    else {
        return b;
    }

}

/*----------------------------------------------------------------------------
| Takes a 64-bit fixed-point value `absZ' with binary point between bits 6
| and 7, and returns the properly rounded 32-bit integer corresponding to the
| input.  If `zSign' is 1, the input is negated before being converted to an
| integer.  Bit 63 of `absZ' must be zero.  Ordinarily, the fixed-point input
| is simply rounded to an integer, with the inexact exception raised if the
| input cannot be represented exactly as an integer.  However, if the fixed-
| point input is too large, the invalid exception is raised and the largest
| positive or negative integer is returned.
*----------------------------------------------------------------------------*/
static int32 roundAndPackInt32( flag zSign, bits64 absZ )
{
    int8 roundingMode;
    flag roundNearestEven;
    int8 roundIncrement, roundBits;
    int32 z;

    roundingMode = float_rounding_mode;
    roundNearestEven = ( roundingMode == float_round_nearest_even );
    roundIncrement = 0x40;
    if ( ! roundNearestEven ) {
        if ( roundingMode == float_round_to_zero ) {
            roundIncrement = 0;
        }
        else {
            roundIncrement = 0x7F;
            if ( zSign ) {
                if ( roundingMode == float_round_up ) roundIncrement = 0;
            }
            else {
                if ( roundingMode == float_round_down ) roundIncrement = 0;
            }
        }
    }
    roundBits = absZ & 0x7F;
    absZ = ( absZ + roundIncrement )>>7;
    absZ &= ~ ( ( ( roundBits ^ 0x40 ) == 0 ) & roundNearestEven );
    z = absZ;
    if ( zSign ) z = - z;
    if ( ( absZ>>32 ) || ( z && ( ( z < 0 ) ^ zSign ) ) ) {
        float_raise( float_flag_invalid );
        return zSign ? (sbits32) 0x80000000 : 0x7FFFFFFF;
    }
    if ( roundBits ) float_exception_flags |= float_flag_inexact;
    return z;

}

/*----------------------------------------------------------------------------
| Takes the 128-bit fixed-point value formed by concatenating `absZ0' and
| `absZ1', with binary point between bits 63 and 64 (between the input words),
| and returns the properly rounded 64-bit integer corresponding to the input.
| If `zSign' is 1, the input is negated before being converted to an integer.
| Ordinarily, the fixed-point input is simply rounded to an integer, with
| the inexact exception raised if the input cannot be represented exactly as
| an integer.  However, if the fixed-point input is too large, the invalid
| exception is raised and the largest positive or negative integer is
| returned.
*----------------------------------------------------------------------------*/

static int64 roundAndPackInt64( flag zSign, bits64 absZ0, bits64 absZ1 )
{
    int8 roundingMode;
    flag roundNearestEven, increment;
    int64 z;

    roundingMode = float_rounding_mode;
    roundNearestEven = ( roundingMode == float_round_nearest_even );
    increment = ( (sbits64) absZ1 < 0 );
    if ( ! roundNearestEven ) {
        if ( roundingMode == float_round_to_zero ) {
            increment = 0;
        }
        else {
            if ( zSign ) {
                increment = ( roundingMode == float_round_down ) && absZ1;
            }
            else {
                increment = ( roundingMode == float_round_up ) && absZ1;
            }
        }
    }
    if ( increment ) {
        ++absZ0;
        if ( absZ0 == 0 ) goto overflow;
        absZ0 &= ~ ( ( (bits64) ( absZ1<<1 ) == 0 ) & roundNearestEven );
    }
    z = absZ0;
    if ( zSign ) z = - z;
    if ( z && ( ( z < 0 ) ^ zSign ) ) {
 overflow:
        float_raise( float_flag_invalid );
        return
              zSign ? (sbits64) LIT64( 0x8000000000000000 )
            : LIT64( 0x7FFFFFFFFFFFFFFF );
    }
    if ( absZ1 ) float_exception_flags |= float_flag_inexact;
    return z;

}



/*----------------------------------------------------------------------------
| Returns the least-significant 64 fraction bits of the quadruple-precision
| floating-point value `a'.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static bits64 extractFloat128Frac1( float128_t a )
{

    return a.low;

}
#else
#define extractFloat128Frac1(a) (a.low)
/*
static inline __declspec(naked) bits64 extractFloat128Frac1(float128_t a)
{
	_asm("\tpextrq\t$0,%xmm0,%rax");
}
*/
#endif

/*----------------------------------------------------------------------------
| Returns the most-significant 48 fraction bits of the quadruple-precision
| floating-point value `a'.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static bits64 extractFloat128Frac0( float128_t a )
{

    return a.high & LIT64( 0x0000FFFFFFFFFFFF );

}
#else
#define extractFloat128Frac0(a) (a.high & LIT64( 0x0000FFFFFFFFFFFF ))
/*
static inline __declspec(naked) bits64 extractFloat128Frac0(float128_t a)
{
	_asm("\tpextrq\t$1,%xmm0,%rax");
	_asm("\tshlq\t$16,%rax");
	_asm("\tshrq\t$16,%rax");
}
*/
#endif
/*----------------------------------------------------------------------------
| Returns the exponent bits of the quadruple-precision floating-point value
| `a'.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static int32 extractFloat128Exp( float128_t a )
{

    return ( a.high>>48 ) & 0x7FFF;

}
#else
#define extractFloat128Exp(a) (( a.high>>48 ) & 0x7FFF)
/*
static int inline __declspec(naked) extractFloat128Exp(float128_t a)
{
	_asm("\tpextrw\t$7,%xmm0,%eax");
	_asm("\tandl\t$0x7fff,%eax");
}
*/
#endif

/*----------------------------------------------------------------------------
| Returns the sign bit of the quadruple-precision floating-point value `a'.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
static flag extractFloat128Sign( float128_t a )
{

    return a.high>>63;

}
#else
#define extractFloat128Sign(a) (a.high >> 63)
/*
static int inline __declspec(naked) extractFloat128Sign(float128_t a)
{
	_asm("\tpextrw\t$7,%xmm0,%rax");
	_asm("\tshrl\t$15,%eax");
}
*/
#endif

/*----------------------------------------------------------------------------
| Normalizes the subnormal quadruple-precision floating-point value
| represented by the denormalized significand formed by the concatenation of
| `aSig0' and `aSig1'.  The normalized exponent is stored at the location
| pointed to by `zExpPtr'.  The most significant 49 bits of the normalized
| significand are stored at the location pointed to by `zSig0Ptr', and the
| least significant 64 bits of the normalized significand are stored at the
| location pointed to by `zSig1Ptr'.
*----------------------------------------------------------------------------*/

void normalizeFloat128Subnormal( bits64 aSig0, bits64 aSig1, int32 *zExpPtr,
     bits64 *zSig0Ptr, bits64 *zSig1Ptr )
{
    int8 shiftCount;
	float128_t tmp;

    if ( aSig0 == 0 ) {
        shiftCount = countLeadingZeros64( aSig1 ) - 15;
        if ( shiftCount < 0 ) {
            *zSig0Ptr = aSig1>>( - shiftCount );
            *zSig1Ptr = aSig1<<( shiftCount & 63 );
        }
        else {
            *zSig0Ptr = aSig1<<shiftCount;
            *zSig1Ptr = 0;
        }
        *zExpPtr = - shiftCount - 63;
    }
    else {
        shiftCount = countLeadingZeros64( aSig0 ) - 15;
        if (shiftCount) {
        	tmp = shortShift128Left(shiftCount, aSig0, aSig1);
		*zSig0Ptr = tmp.high; *zSig1Ptr = tmp.low;
	}
	else *zSig0Ptr = aSig0,*zSig1Ptr = aSig1;
        *zExpPtr = 1 - shiftCount;
    }

}


/*----------------------------------------------------------------------------
| Takes an abstract floating-point value having sign `zSign', exponent `zExp',
| and extended significand formed by the concatenation of `zSig0', `zSig1',
| and `zSig2', and returns the proper quadruple-precision floating-point value
| corresponding to the abstract input.  Ordinarily, the abstract value is
| simply rounded and packed into the quadruple-precision format, with the
| inexact exception raised if the abstract input cannot be represented
| exactly.  However, if the abstract value is too large, the overflow and
| inexact exceptions are raised and an infinity or maximal finite value is
| returned.  If the abstract value is too small, the input value is rounded to
| a subnormal number, and the underflow and inexact exceptions are raised if
| the abstract input cannot be represented exactly as a subnormal quadruple-
| precision floating-point number.
|     The input significand must be normalized or smaller.  If the input
| significand is not normalized, `zExp' must be 0; in that case, the result
| returned is a subnormal number, and it must not require rounding.  In the
| usual case that the input significand is normalized, `zExp' must be 1 less
| than the ``true'' floating-point exponent.  The handling of underflow and
| overflow follows the IEC/IEEE Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
#ifndef ASSEMBLY_MACROS
float128_t roundAndPackFloat128( int zSign, int32 zExp, bits64 zSig0, bits64 zSig1, bits64 zSig2 )
{
    int roundingMode;
    int roundNearestEven, increment, isTiny;

    roundingMode = float_rounding_mode;
    roundNearestEven = ( roundingMode == float_round_nearest_even );
    increment = ( (sbits64) zSig2 < 0 );
    if ( ! roundNearestEven ) goto setRounding;
    if ( 0x7FFD <= (bits32) zExp ) goto specialConditions;
continuation:
    if ( zSig2 ) float_exception_flags |= float_flag_inexact;
    if ( increment ) {
        float128_t tmp = add128( zSig0, zSig1,0,1);
        tmp.low &= ~ ( ( zSig2 + zSig2 == 0 ) & roundNearestEven );
        return packFloat128( zSign, zExp, tmp.high, tmp.low );
    }
    if ( ( zSig0 | zSig1 ) == 0 ) zExp = 0;
    return packFloat128( zSign, zExp, zSig0, zSig1 );
// -------------------------------------------------rounding settings
setRounding:
    if ( roundingMode == float_round_to_zero ) {
        increment = 0;
    }
    else {
        if ( zSign ) {
            increment = ( roundingMode == float_round_down ) && zSig2;
        }
        else {
            increment = ( roundingMode == float_round_up ) && zSig2;
        }
    }
    if ( 0x7FFD > (bits32) zExp ) goto continuation;
// --------------------------------------------Special conditions section
specialConditions:
    if (    ( 0x7FFD < zExp )
         || (    ( zExp == 0x7FFD )
              && eq128(
                     LIT64( 0x0001FFFFFFFFFFFF ),
                     LIT64( 0xFFFFFFFFFFFFFFFF ),
                     zSig0,
                     zSig1
                 )
              && increment
            )
       ) {
        float_raise( float_flag_overflow | float_flag_inexact );
        if (    ( roundingMode == float_round_to_zero )
             || ( zSign && ( roundingMode == float_round_up ) )
             || ( ! zSign && ( roundingMode == float_round_down ) )
           ) {
            return
                packFloat128(
                    zSign,
                    0x7FFE,
                    LIT64( 0x0000FFFFFFFFFFFF ),
                    LIT64( 0xFFFFFFFFFFFFFFFF )
                );
        }
        return packFloat128( zSign, 0x7FFF, 0, 0 );
    }
    if ( zExp < 0 ) {
        isTiny =
               ( float_detect_tininess == float_tininess_before_rounding )
            || ( zExp < -1 )
            || ! increment
            || lt128(
                   zSig0,
                   zSig1,
                   LIT64( 0x0001FFFFFFFFFFFF ),
                   LIT64( 0xFFFFFFFFFFFFFFFF )
               );
        shift128ExtraRightJamming(
            zSig0, zSig1, zSig2, - zExp, &zSig0, &zSig1, &zSig2 );
        zExp = 0;
        if ( isTiny && zSig2 ) float_raise( float_flag_underflow );
        if ( roundNearestEven ) {
            increment = ( (sbits64) zSig2 < 0 );
        }
        else {
            if ( zSign ) {
                increment = ( roundingMode == float_round_down ) && zSig2;
            }
            else {
                increment = ( roundingMode == float_round_up ) && zSig2;
            }
        }
    }
    goto continuation;
}
#else
float128_t roundAndPackFloat128( int zSign, int32 zExp, bits64 zSig0, bits64 zSig1, bits64 zSig2 );
#endif

/*----------------------------------------------------------------------------
| Takes an abstract floating-point value having sign `zSign', exponent `zExp',
| and significand formed by the concatenation of `zSig0' and `zSig1', and
| returns the proper quadruple-precision floating-point value corresponding
| to the abstract input.  This routine is just like `roundAndPackFloat128'
| except that the input significand has fewer bits and does not have to be
| normalized.  In all cases, `zExp' must be 1 less than the ``true'' floating-
| point exponent.
*----------------------------------------------------------------------------*/
float128_t normalizeRoundAndPackFloat128(flag zSign,int32 zExp,bits64 zSig0,bits64 zSig1)
{
    int8 shiftCount;
	float128_t tmp;
    bits64 zSig2;

    if ( zSig0 == 0 ) {
        zSig0 = zSig1;
        zSig1 = 0;
        zExp -= 64;
    }
    shiftCount = countLeadingZeros64( zSig0 ) - 15;
    if ( 0 <= shiftCount ) {
        zSig2 = 0;
	if (shiftCount) {
        	tmp = shortShift128Left(shiftCount, zSig0, zSig1);
		zSig0 = tmp.high, zSig1 = tmp.low;
	}
    }
    else {
        shift128ExtraRightJamming(
            zSig0, zSig1, 0, - shiftCount, &zSig0, &zSig1, &zSig2 );
    }
    zExp -= shiftCount;
    return roundAndPackFloat128( zSign, zExp, zSig0, zSig1, zSig2 );

}

#ifndef ASSEMBLY_MACROS
/*----------------------------------------------------------------------------
| Returns the result of converting the 32-bit two's complement integer `a' to
| the quadruple-precision floating-point format.  The conversion is performed
| according to the IEC/IEEE Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
float128_t _op_cast_int_float128_t( int a )
{
    flag zSign;
    uint32 absA;
    int8 shiftCount;
    bits64 zSig0;

    if ( a == 0 ) return packFloat128Zero();
    zSign = ( a < 0 );
    absA = zSign ? - a : a;
    shiftCount = countLeadingZeros32NoChecks( absA ) + 17;
    zSig0 = absA;
    return packFloat128( zSign, 0x402E - shiftCount, zSig0<<shiftCount, 0 );
}
#else
float128_t __declspec(naked) _op_cast_int_float128_t(int a)
{
//     if ( a == 0 ) return packFloat128Zero();
//     zSign = ( a < 0 );");
	_asm("\torl\t%ecx,%ecx");
	_asm("\tje\t_$CastOfZero");
	_asm("\tsetl	%r10b");
	_asm("\tshlq	$63,%r10");
//     absA = zSign ? - a : a;");
	_asm("\tmovl	%ecx,%eax");
	_asm("\tcdq");
	_asm("\txorl	%edx,%eax");
	_asm("\tsubl	%edx,%eax");
	_asm("\tmovl	%eax,%r8d");
//     shiftCount = countLeadingZeros32NoChecks( absA ) + 17;
	_asm("\tmovl	%eax,%ecx");
	_asm("\tmovl	$31,%eax");
	_asm("\tbsrl	%ecx,%edx");
	_asm("\tsubl	%edx,%eax");
	_asm("\taddl	$17,%eax");
//     zSig0 = absA;
//     return packFloat128( zSign, 0x402E - shiftCount, zSig0<<shiftCount, 0 );
	_asm("\txorl	%r9d,%r9d");
	_asm("\tpinsrq	$0,%r9,%xmm0");
	_asm("\tmovb	%al,%cl");
	_asm("\tsalq	%cl,%r8");
	_asm("\tmovl	$16430,%edx");
	_asm("\tsubl	%eax,%edx");
	_asm("\tshl	$48,%rdx");
	_asm("\taddq	%rdx,%r8");
	_asm("\taddq	%r10,%r8");
	_asm("\tpinsrq	$1,%r8,%xmm0");
	_asm("\tret");
	_asm("_$CastOfZero:");
	_asm("\tpxor\t%xmm0,%xmm0");
	_asm("\tret");
}
#endif

float128_t _op_cast_unsigned_int_float128_t( unsigned a )
{
    int shiftCount;
    bits64 zSig0;

    if ( a == 0 ) return packFloat128Zero();
    shiftCount = countLeadingZeros32NoChecks( a ) + 17;
    zSig0 = a;
    return packFloat128( 0, 0x402E - shiftCount, zSig0<<shiftCount, 0 );
}

/*----------------------------------------------------------------------------
|  int64_to_float128
| Returns the result of converting the 64-bit two's complement integer `a' to
| the quadruple-precision floating-point format.  The conversion is performed
| according to the IEC/IEEE Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
float128_t _op_cast_long_long_float128_t(long long a)
{
    flag zSign;
    uint64 absA;
    int8 shiftCount;
    int32 zExp;
    bits64 zSig0, zSig1;
	float128_t tmp;

    if ( a == 0 ) return packFloat128Zero();
    zSign = ( a < 0 );
    absA = zSign ? - a : a;
    shiftCount = countLeadingZeros64( absA ) + 49;
    zExp = 0x406E - shiftCount;
    if ( 64 <= shiftCount ) {
        zSig1 = 0;
        zSig0 = absA;
        shiftCount -= 64;
    }
    else {
        zSig1 = absA;
        zSig0 = 0;
    }
    if (shiftCount) {
        tmp = shortShift128Left(shiftCount, zSig0, zSig1);
        return packFloat128( zSign, zExp, tmp.high, tmp.low );
    }
    return packFloat128(zSign, zExp, zSig0, zSig1);

}
float128_t _op_cast_unsigned_long_long_float128_t(unsigned long long a)
{
    int8 shiftCount;
    int32 zExp;
    bits64 zSig0, zSig1;
    float128_t tmp;

    if ( a == 0 ) return packFloat128Zero();
    shiftCount = countLeadingZeros64( a ) + 49;
    zExp = 0x406E - shiftCount;
    if ( 64 <= shiftCount ) {
        zSig1 = 0;
        zSig0 = a;
        shiftCount -= 64;
    }
    else {
        zSig1 = a;
        zSig0 = 0;
    }
    if (shiftCount) {
        tmp = shortShift128Left(shiftCount, zSig0, zSig1);
        return packFloat128( 0, zExp, tmp.high, tmp.low );
    }
    return packFloat128(0, zExp, zSig0, zSig1);
}

/*----------------------------------------------------------------------------
| Returns the result of converting the canonical NaN `a' to the extended
| double-precision floating-point format.
*----------------------------------------------------------------------------*/

static long double commonNaNToFloatx80( commonNaNT a )
{
    floatx80 z;

    z.low = LIT64( 0xC000000000000000 ) | ( a.high>>1 );
    z.high = ( ( (bits16) a.sign )<<15 ) | 0x7FFF;
    return *(long double *)&z;

}


/*----------------------------------------------------------------------------
| Returns the result of converting the single-precision floating-point value
| `a' to the double-precision floating-point format.  The conversion is
| performed according to the IEC/IEEE Standard for Binary Floating-Point
| Arithmetic.
*----------------------------------------------------------------------------*/
float128_t _op_cast_float_float128_t(float float_a)
{
    flag aSign;
    int aExp;
    bits32 aSig;
	float32 a = *(float32 *)&float_a;

    aSig = extractFloat32Frac( a );
    aExp = extractFloat32Exp( a );
    aSign = extractFloat32Sign( a );
    if ( aExp == 0xFF ) {
        if ( aSig ) return commonNaNToFloat128( float32ToCommonNaN( a ) );
        return packFloat128( aSign, 0x7FFF, 0, 0 );
    }
    if ( aExp == 0 ) {
        if ( aSig == 0 ) return packFloat128( aSign, 0, 0, 0 );
        normalizeFloat32Subnormal( aSig, &aExp, &aSig );
        --aExp;
    }
    return packFloat128( aSign, aExp + 0x3F80, ( (bits64) aSig )<<25, 0 );

}


/*----------------------------------------------------------------------------
| Returns the result of converting the double-precision floating-point value
| `a' to the quadruple-precision floating-point format.  The conversion is
| performed according to the IEC/IEEE Standard for Binary Floating-Point
| Arithmetic.
*----------------------------------------------------------------------------*/

float128_t _op_cast_double_float128_t(double da)
{
    flag aSign;
    int aExp;
	float64 a = *(float64 *)&da;
    bits64 aSig;
    float128_t tmp;

    aSig = extractFloat64Frac( a );
    aExp = extractFloat64Exp( a );
    aSign = extractFloat64Sign( a );
    if ( aExp == 0x7FF ) {
        if ( aSig ) return commonNaNToFloat128( float64ToCommonNaN( a ) );
        return packFloat128( aSign, 0x7FFF, 0, 0 );
    }
    if ( aExp == 0 ) {
        if ( aSig == 0 ) return packFloat128( aSign, 0, 0, 0 );
        normalizeFloat64Subnormal( aSig, &aExp, &aSig );
        --aExp;
    }
    shift128RightSmallCount( aSig, 0, 4,&tmp);
    return packFloat128( aSign, aExp + 0x3C00, tmp.high, tmp.low );

}


/*----------------------------------------------------------------------------
| Returns 1 if the extended double-precision floating-point value `a' is a
| signaling NaN; otherwise returns 0.
*----------------------------------------------------------------------------*/

flag floatx80_is_signaling_nan( floatx80 a )
{
    bits64 aLow;

    aLow = a.low & ~ LIT64( 0x4000000000000000 );
    return
           ( ( a.high & 0x7FFF ) == 0x7FFF )
        && (bits64) ( aLow<<1 )
        && ( a.low == aLow );

}

/*----------------------------------------------------------------------------
| Returns the result of converting the extended double-precision floating-
| point NaN `a' to the canonical NaN format.  If `a' is a signaling NaN, the
| invalid exception is raised.
*----------------------------------------------------------------------------*/
static commonNaNT floatx80ToCommonNaN( floatx80 a )
{
    commonNaNT z;

    if ( floatx80_is_signaling_nan( a ) ) float_raise( float_flag_invalid );
    z.sign = a.high>>15;
    z.low = 0;
    z.high = a.low<<1;
    return z;

}


/*----------------------------------------------------------------------------
| Returns the result of converting the extended double-precision floating-
| point value `a' to the quadruple-precision floating-point format.  The
| conversion is performed according to the IEC/IEEE Standard for Binary
| Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
float128_t  _op_cast_long_double_float128_t(long double A)
{
    flag aSign;
    int aExp;
    bits64 aSig;
    floatx80 a = *(floatx80 *)&A;
    float128_t tmp;

    aSig = extractFloatx80Frac( a );
    aExp = extractFloatx80Exp( a );
    aSign = extractFloatx80Sign( a );
    if ( ( aExp == 0x7FFF ) && (bits64) ( aSig<<1 ) ) {
        return commonNaNToFloat128( floatx80ToCommonNaN( a ) );
    }
    shift128RightSmallCount( aSig<<1, 0, 16 ,&tmp);
    return packFloat128( aSign, aExp, tmp.high, tmp.low );

}


/*----------------------------------------------------------------------------
| Returns the result of converting the quadruple-precision floating-point
| value `a' to the 32-bit two's complement integer format.  The conversion
| is performed according to the IEC/IEEE Standard for Binary Floating-Point
| Arithmetic---which means in particular that the conversion is rounded
| according to the current rounding mode.  If `a' is a NaN, the largest
| positive integer is returned.  Otherwise, if the conversion overflows, the
| largest integer with the same sign as `a' is returned.
*----------------------------------------------------------------------------*/
int roundToInt( float128_t a )
{
    flag aSign;
    int32 aExp, shiftCount;
    bits64 aSig0, aSig1;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    if ( ( aExp == 0x7FFF ) && ( aSig0 | aSig1 ) ) aSign = 0;
    if ( aExp ) aSig0 |= LIT64( 0x0001000000000000 );
    aSig0 |= ( aSig1 != 0 );
    shiftCount = 0x4028 - aExp;
    if ( 0 < shiftCount ) shift64RightJamming( aSig0, shiftCount, &aSig0 );
    return roundAndPackInt32( aSign, aSig0 );

}

/*----------------------------------------------------------------------------
| Returns the result of converting the quadruple-precision floating-point
| value `a' to the 32-bit two's complement integer format.  The conversion
| is performed according to the IEC/IEEE Standard for Binary Floating-Point
| Arithmetic, except that the conversion is always rounded toward zero.  If
| `a' is a NaN, the largest positive integer is returned.  Otherwise, if the
| conversion overflows, the largest integer with the same sign as `a' is
| returned.
*----------------------------------------------------------------------------*/
int32 _op_cast_float128_t_int( float128_t a )
{
    flag aSign;
    int32 aExp, shiftCount;
    bits64 aSig0, aSig1, savedASig;
    int32 z;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    aSig0 |= ( aSig1 != 0 );
    if ( 0x401E < aExp ) {
        if ( ( aExp == 0x7FFF ) && aSig0 ) aSign = 0;
        goto invalid;
    }
    else if ( aExp < 0x3FFF ) {
        if ( aExp || aSig0 ) float_exception_flags |= float_flag_inexact;
        return 0;
    }
    aSig0 |= LIT64( 0x0001000000000000 );
    shiftCount = 0x402F - aExp;
    savedASig = aSig0;
    aSig0 >>= shiftCount;
    z = aSig0;
    if ( aSign ) z = - z;
    if ( ( z < 0 ) ^ aSign ) {
 invalid:
        float_raise( float_flag_invalid );
        return aSign ? (sbits32) 0x80000000 : 0x7FFFFFFF;
    }
    if ( ( aSig0<<shiftCount ) != savedASig ) {
        float_exception_flags |= float_flag_inexact;
    }
    return z;

}

/*----------------------------------------------------------------------------
| Returns the result of converting the quadruple-precision floating-point
| value `a' to the 64-bit two's complement integer format.  The conversion
| is performed according to the IEC/IEEE Standard for Binary Floating-Point
| Arithmetic---which means in particular that the conversion is rounded
| according to the current rounding mode.  If `a' is a NaN, the largest
| positive integer is returned.  Otherwise, if the conversion overflows, the
| largest integer with the same sign as `a' is returned.
*----------------------------------------------------------------------------*/
int64 float128_to_int64( float128_t a )
{
    flag aSign;
    int32 aExp, shiftCount;
    bits64 aSig0, aSig1;
	float128_t tmp;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    if ( aExp ) aSig0 |= LIT64( 0x0001000000000000 );
    shiftCount = 0x402F - aExp;
    if ( shiftCount <= 0 ) {
        if ( 0x403E < aExp ) {
            float_raise( float_flag_invalid );
            if (    ! aSign
                 || (    ( aExp == 0x7FFF )
                      && ( aSig1 || ( aSig0 != LIT64( 0x0001000000000000 ) ) )
                    )
               ) {
                return LIT64( 0x7FFFFFFFFFFFFFFF );
            }
            return (sbits64) LIT64( 0x8000000000000000 );
        }
        if (shiftCount) {
            tmp = shortShift128Left(- shiftCount, aSig0, aSig1);
            aSig0 = tmp.high, aSig1 = tmp.low;
        }
    }
    else {
        shift64ExtraRightJamming( aSig0, aSig1, shiftCount, &aSig0, &aSig1 );
    }
    return roundAndPackInt64( aSign, aSig0, aSig1 );

}


/*----------------------------------------------------------------------------
*float128_to_int64_round_to_zero
| Returns the result of converting the quadruple-precision floating-point
| value `a' to the 64-bit two's complement integer format.  The conversion
| is performed according to the IEC/IEEE Standard for Binary Floating-Point
| Arithmetic, except that the conversion is always rounded toward zero.
| If `a' is a NaN, the largest positive integer is returned.  Otherwise, if
| the conversion overflows, the largest integer with the same sign as `a' is
| returned.
*----------------------------------------------------------------------------*/
long long  _op_cast_float128_t_long_long(float128_t a)
{
    flag aSign;
    int32 aExp, shiftCount;
    bits64 aSig0, aSig1;
    int64 z;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    if ( aExp ) aSig0 |= LIT64( 0x0001000000000000 );
    shiftCount = aExp - 0x402F;
    if ( 0 < shiftCount ) {
        if ( 0x403E <= aExp ) {
            aSig0 &= LIT64( 0x0000FFFFFFFFFFFF );
            if (    ( a.high == LIT64( 0xC03E000000000000 ) )
                 && ( aSig1 < LIT64( 0x0002000000000000 ) ) ) {
                if ( aSig1 ) float_exception_flags |= float_flag_inexact;
            }
            else {
                float_raise( float_flag_invalid );
                if ( ! aSign || ( ( aExp == 0x7FFF ) && ( aSig0 | aSig1 ) ) ) {
                    return LIT64( 0x7FFFFFFFFFFFFFFF );
                }
            }
            return (sbits64) LIT64( 0x8000000000000000 );
        }
        z = ( aSig0<<shiftCount ) | ( aSig1>>( ( - shiftCount ) & 63 ) );
        if ( (bits64) ( aSig1<<shiftCount ) ) {
            float_exception_flags |= float_flag_inexact;
        }
    }
    else {
        if ( aExp < 0x3FFF ) {
            if ( aExp | aSig0 | aSig1 ) {
                float_exception_flags |= float_flag_inexact;
            }
            return 0;
        }
        z = aSig0>>( - shiftCount );
        if (    aSig1
             || ( shiftCount && (bits64) ( aSig0<<( shiftCount & 63 ) ) ) ) {
            float_exception_flags |= float_flag_inexact;
        }
    }
    if ( aSign ) z = - z;
    return z;

}

unsigned _op_cast_float128_t_unsigned_int(float128_t a)
{
	return _op_cast_float128_t_long_long(a);
}

/*----------------------------------------------------------------------------
| Returns the result of converting the canonical NaN `a' to the single-
| precision floating-point format.
*----------------------------------------------------------------------------*/

static float commonNaNToFloat32( commonNaNT a )
{

    float32 s = ( ( (bits32) a.sign )<<31 ) | 0x7FC00000 | ( a.high>>41 );
	return *(float *)&s;

}

/*----------------------------------------------------------------------------
| Packs the sign `zSign', exponent `zExp', and significand `zSig' into a
| single-precision floating-point value, returning the result.  After being
| shifted into the proper positions, the three fields are simply added
| together to form the result.  This means that any integer portion of `zSig'
| will be added into the exponent.  Since a properly normalized significand
| will have an integer portion equal to 1, the `zExp' input should be 1 less
| than the desired result exponent whenever `zSig' is a complete, normalized
| significand.
*----------------------------------------------------------------------------*/

static float32 packFloat32( flag zSign, int zExp, bits32 zSig )
{

    return ( ( (bits32) zSign )<<31 ) + ( ( (bits32) zExp )<<23 ) + zSig;

}
/*----------------------------------------------------------------------------
| Takes an abstract floating-point value having sign `zSign', exponent `zExp',
| and significand `zSig', and returns the proper single-precision floating-
| point value corresponding to the abstract input.  Ordinarily, the abstract
| value is simply rounded and packed into the single-precision format, with
| the inexact exception raised if the abstract input cannot be represented
| exactly.  However, if the abstract value is too large, the overflow and
| inexact exceptions are raised and an infinity or maximal finite value is
| returned.  If the abstract value is too small, the input value is rounded to
| a subnormal number, and the underflow and inexact exceptions are raised if
| the abstract input cannot be represented exactly as a subnormal single-
| precision floating-point number.
|     The input significand `zSig' has its binary point between bits 30
| and 29, which is 7 bits to the left of the usual location.  This shifted
| significand must be normalized or smaller.  If `zSig' is not normalized,
| `zExp' must be 0; in that case, the result returned is a subnormal number,
| and it must not require rounding.  In the usual case that `zSig' is
| normalized, `zExp' must be 1 less than the ``true'' floating-point exponent.
| The handling of underflow and overflow follows the IEC/IEEE Standard for
| Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/

static float32 roundAndPackFloat32( flag zSign, int zExp, bits32 zSig )
{
    int8 roundingMode;
    flag roundNearestEven;
    int8 roundIncrement, roundBits;
    flag isTiny;

    roundingMode = float_rounding_mode;
    roundNearestEven = ( roundingMode == float_round_nearest_even );
    roundIncrement = 0x40;
    if ( ! roundNearestEven ) {
        if ( roundingMode == float_round_to_zero ) {
            roundIncrement = 0;
        }
        else {
            roundIncrement = 0x7F;
            if ( zSign ) {
                if ( roundingMode == float_round_up ) roundIncrement = 0;
            }
            else {
                if ( roundingMode == float_round_down ) roundIncrement = 0;
            }
        }
    }
    roundBits = zSig & 0x7F;
    if ( 0xFD <= (bits16) zExp ) {
        if (    ( 0xFD < zExp )
             || (    ( zExp == 0xFD )
                  && ( (sbits32) ( zSig + roundIncrement ) < 0 ) )
           ) {
            float_raise( float_flag_overflow | float_flag_inexact );
            return packFloat32( zSign, 0xFF, 0 ) - ( roundIncrement == 0 );
        }
        if ( zExp < 0 ) {
            isTiny =
                   ( float_detect_tininess == float_tininess_before_rounding )
                || ( zExp < -1 )
                || ( zSig + roundIncrement < 0x80000000 );
            shift32RightJamming( zSig, - zExp, &zSig );
            zExp = 0;
            roundBits = zSig & 0x7F;
            if ( isTiny && roundBits ) float_raise( float_flag_underflow );
        }
    }
    if ( roundBits ) float_exception_flags |= float_flag_inexact;
    zSig = ( zSig + roundIncrement )>>7;
    zSig &= ~ ( ( ( roundBits ^ 0x40 ) == 0 ) & roundNearestEven );
    if ( zSig == 0 ) zExp = 0;
    return packFloat32( zSign, zExp, zSig );

}


/*----------------------------------------------------------------------------
| float128_to_float32
| Returns the result of converting the quadruple-precision floating-point
| value `a' to the single-precision floating-point format.  The conversion
| is performed according to the IEC/IEEE Standard for Binary Floating-Point
| Arithmetic.
*----------------------------------------------------------------------------*/
float _op_cast_float128_t_float(float128_t a)
{
    flag aSign;
    int32 aExp;
    bits64 aSig0, aSig1;
    bits32 zSig;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    if ( aExp == 0x7FFF ) {
        if ( aSig0 | aSig1 ) {
            return commonNaNToFloat32( float128ToCommonNaN( a ) );
        }
        return packFloat32( aSign, 0xFF, 0 );
    }
    aSig0 |= ( aSig1 != 0 );
    shift64RightJamming( aSig0, 18, &aSig0 );
    zSig = aSig0;
    if ( aExp || zSig ) {
        zSig |= 0x40000000;
        aExp -= 0x3F81;
    }
    return roundAndPackFloat32( aSign, aExp, zSig );

}
/*----------------------------------------------------------------------------
| Returns the result of converting the canonical NaN `a' to the double-
| precision floating-point format.
*----------------------------------------------------------------------------*/

static float64 commonNaNToFloat64( commonNaNT a )
{

    return
          ( ( (bits64) a.sign )<<63 )
        | LIT64( 0x7FF8000000000000 )
        | ( a.high>>12 );

}

/*----------------------------------------------------------------------------
| Packs the sign `zSign', exponent `zExp', and significand `zSig' into a
| double-precision floating-point value, returning the result.  After being
| shifted into the proper positions, the three fields are simply added
| together to form the result.  This means that any integer portion of `zSig'
| will be added into the exponent.  Since a properly normalized significand
| will have an integer portion equal to 1, the `zExp' input should be 1 less
| than the desired result exponent whenever `zSig' is a complete, normalized
| significand.
*----------------------------------------------------------------------------*/

static float64 packFloat64( flag zSign, int zExp, bits64 zSig )
{

    return ( ( (bits64) zSign )<<63 ) + ( ( (bits64) zExp )<<52 ) + zSig;

}

/*----------------------------------------------------------------------------
| Takes an abstract floating-point value having sign `zSign', exponent `zExp',
| and significand `zSig', and returns the proper double-precision floating-
| point value corresponding to the abstract input.  Ordinarily, the abstract
| value is simply rounded and packed into the double-precision format, with
| the inexact exception raised if the abstract input cannot be represented
| exactly.  However, if the abstract value is too large, the overflow and
| inexact exceptions are raised and an infinity or maximal finite value is
| returned.  If the abstract value is too small, the input value is rounded
| to a subnormal number, and the underflow and inexact exceptions are raised
| if the abstract input cannot be represented exactly as a subnormal double-
| precision floating-point number.
|     The input significand `zSig' has its binary point between bits 62
| and 61, which is 10 bits to the left of the usual location.  This shifted
| significand must be normalized or smaller.  If `zSig' is not normalized,
| `zExp' must be 0; in that case, the result returned is a subnormal number,
| and it must not require rounding.  In the usual case that `zSig' is
| normalized, `zExp' must be 1 less than the ``true'' floating-point exponent.
| The handling of underflow and overflow follows the IEC/IEEE Standard for
| Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/

static float64 roundAndPackFloat64( flag zSign, int zExp, bits64 zSig )
{
    int8 roundingMode;
    flag roundNearestEven;
    int roundIncrement, roundBits;
    flag isTiny;

    roundingMode = float_rounding_mode;
    roundNearestEven = ( roundingMode == float_round_nearest_even );
    roundIncrement = 0x200;
    if ( ! roundNearestEven ) {
        if ( roundingMode == float_round_to_zero ) {
            roundIncrement = 0;
        }
        else {
            roundIncrement = 0x3FF;
            if ( zSign ) {
                if ( roundingMode == float_round_up ) roundIncrement = 0;
            }
            else {
                if ( roundingMode == float_round_down ) roundIncrement = 0;
            }
        }
    }
    roundBits = zSig & 0x3FF;
    if ( 0x7FD <= (bits16) zExp ) {
        if (    ( 0x7FD < zExp )
             || (    ( zExp == 0x7FD )
                  && ( (sbits64) ( zSig + roundIncrement ) < 0 ) )
           ) {
            float_raise( float_flag_overflow | float_flag_inexact );
            return packFloat64( zSign, 0x7FF, 0 ) - ( roundIncrement == 0 );
        }
        if ( zExp < 0 ) {
            isTiny =
                   ( float_detect_tininess == float_tininess_before_rounding )
                || ( zExp < -1 )
                || ( zSig + roundIncrement < LIT64( 0x8000000000000000 ) );
            shift64RightJamming( zSig, - zExp, &zSig );
            zExp = 0;
            roundBits = zSig & 0x3FF;
            if ( isTiny && roundBits ) float_raise( float_flag_underflow );
        }
    }
    if ( roundBits ) float_exception_flags |= float_flag_inexact;
    zSig = ( zSig + roundIncrement )>>10;
    zSig &= ~ ( ( ( roundBits ^ 0x200 ) == 0 ) & roundNearestEven );
    if ( zSig == 0 ) zExp = 0;
    return packFloat64( zSign, zExp, zSig );

}

/*----------------------------------------------------------------------------
| Returns the result of converting the quadruple-precision floating-point
| value `a' to the double-precision floating-point format.  The conversion
| is performed according to the IEC/IEEE Standard for Binary Floating-Point
| Arithmetic.
*----------------------------------------------------------------------------*/
double _op_cast_float128_t_double( float128_t a )
{
    flag aSign;
    int32 aExp;
    bits64 aSig0, aSig1;
	float128_t tmp;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    if ( aExp == 0x7FFF ) {
        if ( aSig0 | aSig1 ) {
            return commonNaNToFloat64( float128ToCommonNaN( a ) );
        }
        return packFloat64( aSign, 0x7FF, 0 );
    }
    tmp = shortShift128Left(14, aSig0, aSig1);
    aSig0 = tmp.high, aSig1 = tmp.low;
    aSig0 |= ( aSig1 != 0 );
    if ( aExp || aSig0 ) {
        aSig0 |= LIT64( 0x4000000000000000 );
        aExp -= 0x3C01;
    }
    aSig1 = roundAndPackFloat64( aSign, aExp, aSig0 );
	return *(double *)&aSig1;

}

/*----------------------------------------------------------------------------
| Packs the sign `zSign', exponent `zExp', and significand `zSig' into an
| extended double-precision floating-point value, returning the result.
*----------------------------------------------------------------------------*/

static long double packFloatx80( flag zSign, int32 zExp, bits64 zSig )
{
    floatx80 z;

    z.low = zSig;
    z.high = ( ( (bits16) zSign )<<15 ) + zExp;
    return *(long double *)&z;

}

/*----------------------------------------------------------------------------
| Takes an abstract floating-point value having sign `zSign', exponent `zExp',
| and extended significand formed by the concatenation of `zSig0' and `zSig1',
| and returns the proper extended double-precision floating-point value
| corresponding to the abstract input.  Ordinarily, the abstract value is
| rounded and packed into the extended double-precision format, with the
| inexact exception raised if the abstract input cannot be represented
| exactly.  However, if the abstract value is too large, the overflow and
| inexact exceptions are raised and an infinity or maximal finite value is
| returned.  If the abstract value is too small, the input value is rounded to
| a subnormal number, and the underflow and inexact exceptions are raised if
| the abstract input cannot be represented exactly as a subnormal extended
| double-precision floating-point number.
|     If `roundingPrecision' is 32 or 64, the result is rounded to the same
| number of bits as single or double precision, respectively.  Otherwise, the
| result is rounded to the full precision of the extended double-precision
| format.
|     The input significand must be normalized or smaller.  If the input
| significand is not normalized, `zExp' must be 0; in that case, the result
| returned is a subnormal number, and it must not require rounding.  The
| handling of underflow and overflow follows the IEC/IEEE Standard for Binary
| Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
static long double  roundAndPackFloatx80(
     int8 roundingPrecision, flag zSign, int32 zExp, bits64 zSig0, bits64 zSig1
 )
{
    int8 roundingMode;
    flag roundNearestEven, increment, isTiny;
    int64 roundIncrement, roundMask, roundBits;

    roundingMode = float_rounding_mode;
    roundNearestEven = ( roundingMode == float_round_nearest_even );
    if ( roundingPrecision == 80 ) goto precision80;
    if ( roundingPrecision == 64 ) {
        roundIncrement = LIT64( 0x0000000000000400 );
        roundMask = LIT64( 0x00000000000007FF );
    }
    else if ( roundingPrecision == 32 ) {
        roundIncrement = LIT64( 0x0000008000000000 );
        roundMask = LIT64( 0x000000FFFFFFFFFF );
    }
    else {
        goto precision80;
    }
    zSig0 |= ( zSig1 != 0 );
    if ( ! roundNearestEven ) {
        if ( roundingMode == float_round_to_zero ) {
            roundIncrement = 0;
        }
        else {
            roundIncrement = roundMask;
            if ( zSign ) {
                if ( roundingMode == float_round_up ) roundIncrement = 0;
            }
            else {
                if ( roundingMode == float_round_down ) roundIncrement = 0;
            }
        }
    }
    roundBits = zSig0 & roundMask;
    if ( 0x7FFD <= (bits32) ( zExp - 1 ) ) {
        if (    ( 0x7FFE < zExp )
             || ( ( zExp == 0x7FFE ) && ( zSig0 + roundIncrement < zSig0 ) )
           ) {
            goto overflow;
        }
        if ( zExp <= 0 ) {
            isTiny =
                   ( float_detect_tininess == float_tininess_before_rounding )
                || ( zExp < 0 )
                || ( zSig0 <= zSig0 + roundIncrement );
            shift64RightJamming( zSig0, 1 - zExp, &zSig0 );
            zExp = 0;
            roundBits = zSig0 & roundMask;
            if ( isTiny && roundBits ) float_raise( float_flag_underflow );
            if ( roundBits ) float_exception_flags |= float_flag_inexact;
            zSig0 += roundIncrement;
            if ( (sbits64) zSig0 < 0 ) zExp = 1;
            roundIncrement = roundMask + 1;
            if ( roundNearestEven && ( roundBits<<1 == roundIncrement ) ) {
                roundMask |= roundIncrement;
            }
            zSig0 &= ~ roundMask;
            return packFloatx80( zSign, zExp, zSig0 );
        }
    }
    if ( roundBits ) float_exception_flags |= float_flag_inexact;
    zSig0 += roundIncrement;
    if ( zSig0 < roundIncrement ) {
        ++zExp;
        zSig0 = LIT64( 0x8000000000000000 );
    }
    roundIncrement = roundMask + 1;
    if ( roundNearestEven && ( roundBits<<1 == roundIncrement ) ) {
        roundMask |= roundIncrement;
    }
    zSig0 &= ~ roundMask;
    if ( zSig0 == 0 ) zExp = 0;
    return packFloatx80( zSign, zExp, zSig0 );
 precision80:
    increment = ( (sbits64) zSig1 < 0 );
    if ( ! roundNearestEven ) {
        if ( roundingMode == float_round_to_zero ) {
            increment = 0;
        }
        else {
            if ( zSign ) {
                increment = ( roundingMode == float_round_down ) && zSig1;
            }
            else {
                increment = ( roundingMode == float_round_up ) && zSig1;
            }
        }
    }
    if ( 0x7FFD <= (bits32) ( zExp - 1 ) ) {
        if (    ( 0x7FFE < zExp )
             || (    ( zExp == 0x7FFE )
                  && ( zSig0 == LIT64( 0xFFFFFFFFFFFFFFFF ) )
                  && increment
                )
           ) {
            roundMask = 0;
 overflow:
            float_raise( float_flag_overflow | float_flag_inexact );
            if (    ( roundingMode == float_round_to_zero )
                 || ( zSign && ( roundingMode == float_round_up ) )
                 || ( ! zSign && ( roundingMode == float_round_down ) )
               ) {
                return packFloatx80( zSign, 0x7FFE, ~ roundMask );
            }
            return packFloatx80( zSign, 0x7FFF, LIT64( 0x8000000000000000 ) );
        }
        if ( zExp <= 0 ) {
            isTiny =
                   ( float_detect_tininess == float_tininess_before_rounding )
                || ( zExp < 0 )
                || ! increment
                || ( zSig0 < LIT64( 0xFFFFFFFFFFFFFFFF ) );
            shift64ExtraRightJamming( zSig0, zSig1, 1 - zExp, &zSig0, &zSig1 );
            zExp = 0;
            if ( isTiny && zSig1 ) float_raise( float_flag_underflow );
            if ( zSig1 ) float_exception_flags |= float_flag_inexact;
            if ( roundNearestEven ) {
                increment = ( (sbits64) zSig1 < 0 );
            }
            else {
                if ( zSign ) {
                    increment = ( roundingMode == float_round_down ) && zSig1;
                }
                else {
                    increment = ( roundingMode == float_round_up ) && zSig1;
                }
            }
            if ( increment ) {
                ++zSig0;
                zSig0 &=
                    ~ ( ( (bits64) ( zSig1<<1 ) == 0 ) & roundNearestEven );
                if ( (sbits64) zSig0 < 0 ) zExp = 1;
            }
            return packFloatx80( zSign, zExp, zSig0 );
        }
    }
    if ( zSig1 ) float_exception_flags |= float_flag_inexact;
    if ( increment ) {
        ++zSig0;
        if ( zSig0 == 0 ) {
            ++zExp;
            zSig0 = LIT64( 0x8000000000000000 );
        }
        else {
            zSig0 &= ~ ( ( (bits64) ( zSig1<<1 ) == 0 ) & roundNearestEven );
        }
    }
    else {
        if ( zSig0 == 0 ) zExp = 0;
    }
    return packFloatx80( zSign, zExp, zSig0 );

}


/*----------------------------------------------------------------------------
| float128_to_floatx80
| Returns the result of converting the quadruple-precision floating-point
| value `a' to the extended double-precision floating-point format.  The
| conversion is performed according to the IEC/IEEE Standard for Binary
| Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
long double _op_cast_float128_t_long_double(float128_t a)
{
    flag aSign;
    int32 aExp;
    bits64 aSig0, aSig1;
	float128_t tmp;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    if ( aExp == 0x7FFF ) {
        if ( aSig0 | aSig1 ) {
            return commonNaNToFloatx80( float128ToCommonNaN( a ) );
        }
        return packFloatx80( aSign, 0x7FFF, LIT64( 0x8000000000000000 ) );
    }
    if ( aExp == 0 ) {
        if ( ( aSig0 | aSig1 ) == 0 ) return packFloatx80( aSign, 0, 0 );
        normalizeFloat128Subnormal( aSig0, aSig1, &aExp, &aSig0, &aSig1 );
    }
    else {
        aSig0 |= LIT64( 0x0001000000000000 );
    }
    tmp = shortShift128Left(15, aSig0, aSig1);
    return roundAndPackFloatx80( 80, aSign, aExp, tmp.high, tmp.low );

}


/*----------------------------------------------------------------------------
| Rounds the quadruple-precision floating-point value `a' to an integer, and
| returns the result as a quadruple-precision floating-point value.  The
| operation is performed according to the IEC/IEEE Standard for Binary
| Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
float128_t float128_round_to_int( float128_t a )
{
    flag aSign;
    int32 aExp;
    bits64 lastBitMask, roundBitsMask;
    int8 roundingMode;
    float128_t z;

    aExp = extractFloat128Exp( a );
    if ( 0x402F <= aExp ) {
        if ( 0x406F <= aExp ) {
            if (    ( aExp == 0x7FFF )
                 && ( extractFloat128Frac0( a ) | extractFloat128Frac1( a ) )
               ) {
                return propagateFloat128NaN( a, a );
            }
            return a;
        }
        lastBitMask = 1;
        lastBitMask = ( lastBitMask<<( 0x406E - aExp ) )<<1;
        roundBitsMask = lastBitMask - 1;
        z = a;
        roundingMode = float_rounding_mode;
        if ( roundingMode == float_round_nearest_even ) {
            if ( lastBitMask ) {
                z = add128( z.high, z.low, 0, lastBitMask>>1);
                if ( ( z.low & roundBitsMask ) == 0 ) z.low &= ~ lastBitMask;
            }
            else {
                if ( (sbits64) z.low < 0 ) {
                    ++z.high;
                    if ( (bits64) ( z.low<<1 ) == 0 ) z.high &= ~1;
                }
            }
        }
        else if ( roundingMode != float_round_to_zero ) {
            if (   extractFloat128Sign( z )
                 ^ ( roundingMode == float_round_up ) ) {
                z = add128( z.high, z.low, 0, roundBitsMask);
            }
        }
        z.low &= ~ roundBitsMask;
    }
    else {
        if ( aExp < 0x3FFF ) {
            if ( ( ( (bits64) ( a.high<<1 ) ) | a.low ) == 0 ) return a;
            float_exception_flags |= float_flag_inexact;
            aSign = extractFloat128Sign( a );
            switch ( float_rounding_mode ) {
             case float_round_nearest_even:
                if (    ( aExp == 0x3FFE )
                     && (   extractFloat128Frac0( a )
                          | extractFloat128Frac1( a ) )
                   ) {
                    return packFloat128( aSign, 0x3FFF, 0, 0 );
                }
                break;
             case float_round_down:
                return
                      aSign ? packFloat128( 1, 0x3FFF, 0, 0 )
                    : packFloat128Zero();
             case float_round_up:
                return
                      aSign ? packFloat128( 1, 0, 0, 0 )
                    : packFloat128( 0, 0x3FFF, 0, 0 );
            }
            return packFloat128( aSign, 0, 0, 0 );
        }
        lastBitMask = 1;
        lastBitMask <<= 0x402F - aExp;
        roundBitsMask = lastBitMask - 1;
        z.low = 0;
        z.high = a.high;
        roundingMode = float_rounding_mode;
        if ( roundingMode == float_round_nearest_even ) {
            z.high += lastBitMask>>1;
            if ( ( ( z.high & roundBitsMask ) | a.low ) == 0 ) {
                z.high &= ~ lastBitMask;
            }
        }
        else if ( roundingMode != float_round_to_zero ) {
            if (   extractFloat128Sign( z )
                 ^ ( roundingMode == float_round_up ) ) {
                z.high |= ( a.low != 0 );
                z.high += roundBitsMask;
            }
        }
        z.high &= ~ roundBitsMask;
    }
    if ( ( z.low != a.low ) || ( z.high != a.high ) ) {
        float_exception_flags |= float_flag_inexact;
    }
    return z;

}
#ifndef ASSEMBLY_MACROS
/*----------------------------------------------------------------------------
| Returns the result of adding the absolute values of the quadruple-precision
| floating-point values `a' and `b'.  If `zSign' is 1, the sum is negated
| before being returned.  `zSign' is ignored if the result is a NaN.
| The addition is performed according to the IEC/IEEE Standard for Binary
| Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
static float128_t addFloat128Sigs( float128_t a, float128_t b, flag zSign )
{
    int32 aExp, bExp, zExp;
    bits64 aSig0, aSig1, bSig0, bSig1, zSig0, zSig1, zSig2;
    int32 expDiff;
    float128_t tmp;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    bSig1 = extractFloat128Frac1( b );
    bSig0 = extractFloat128Frac0( b );
    bExp = extractFloat128Exp( b );
    expDiff = aExp - bExp;
    if ( 0 < expDiff ) {
        if ( aExp == 0x7FFF ) {
            if ( aSig0 | aSig1 ) return propagateFloat128NaN( a, b );
            return a;
        }
        if ( bExp == 0 ) {
            --expDiff;
        }
        else {
            bSig0 |= LIT64( 0x0001000000000000 );
        }
        shift128ExtraRightJamming(
            bSig0, bSig1, 0, expDiff, &bSig0, &bSig1, &zSig2 );
        zExp = aExp;
    }
    else if ( expDiff < 0 ) {
        if ( bExp == 0x7FFF ) {
            if ( bSig0 | bSig1 ) return propagateFloat128NaN( a, b );
            return packFloat128( zSign, 0x7FFF, 0, 0 );
        }
        if ( aExp == 0 ) {
            ++expDiff;
        }
        else {
            aSig0 |= LIT64( 0x0001000000000000 );
        }
        shift128ExtraRightJamming(
            aSig0, aSig1, 0, - expDiff, &aSig0, &aSig1, &zSig2 );
        zExp = bExp;
    }
    else {
        if ( aExp == 0x7FFF ) {
            if ( aSig0 | aSig1 | bSig0 | bSig1 ) {
                return propagateFloat128NaN( a, b );
            }
            return a;
        }
        tmp = add128( aSig0, aSig1, bSig0, bSig1);
        zSig0 = tmp.high, zSig1 = tmp.low;
        if ( aExp == 0 ) return packFloat128( zSign, 0, zSig0, zSig1 );
        zSig2 = 0;
        zSig0 |= LIT64( 0x0002000000000000 );
        zExp = aExp;
        goto shiftRight1;
    }
    aSig0 |= LIT64( 0x0001000000000000 );
    tmp = add128( aSig0, aSig1, bSig0, bSig1);
    zSig0 = tmp.high, zSig1 = tmp.low;
    --zExp;
    if ( tmp.high < LIT64( 0x0002000000000000 ) ) goto roundAndPack;
    ++zExp;
 shiftRight1:
    tmp = shift128ExtraRightJammingByOne( zSig0, zSig1, zSig2, &zSig2 );
    zSig0 = tmp.high, zSig1 = tmp.low;
 roundAndPack:
    return roundAndPackFloat128( zSign, zExp, zSig0, zSig1, zSig2 );

}

/*----------------------------------------------------------------------------
| Returns the result of subtracting the absolute values of the quadruple-
| precision floating-point values `a' and `b'.  If `zSign' is 1, the
| difference is negated before being returned.  `zSign' is ignored if the
| result is a NaN.  The subtraction is performed according to the IEC/IEEE
| Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
static float128_t subFloat128Sigs( float128_t a, float128_t b, flag zSign )
{
    int32 aExp, bExp, zExp;
    bits64 aSig0, aSig1, bSig0, bSig1, zSig0, zSig1;
    int32 expDiff;
    float128_t z,tmp;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    bSig1 = extractFloat128Frac1( b );
    bSig0 = extractFloat128Frac0( b );
    bExp = extractFloat128Exp( b );
    expDiff = aExp - bExp;
    tmp = shortShift128Left(14, aSig0, aSig1);
	aSig0 = tmp.high,aSig1 = tmp.low;
    tmp = shortShift128Left(14, bSig0, bSig1);
	bSig0 = tmp.high, bSig1 = tmp.low;
    if ( 0 < expDiff ) goto aExpBigger;
    if ( expDiff < 0 ) goto bExpBigger;
    if ( aExp == 0x7FFF ) {
        if ( aSig0 | aSig1 | bSig0 | bSig1 ) {
            return propagateFloat128NaN( a, b );
        }
        float_raise( float_flag_invalid );
        z.low = float128_default_nan_low;
        z.high = float128_default_nan_high;
        return z;
    }
    if ( aExp == 0 ) {
        aExp = 1;
        bExp = 1;
    }
    if ( bSig0 < aSig0 ) goto aBigger;
    if ( aSig0 < bSig0 ) goto bBigger;
    if ( bSig1 < aSig1 ) goto aBigger;
    if ( aSig1 < bSig1 ) goto bBigger;
    return packFloat128( float_rounding_mode == float_round_down, 0, 0, 0 );
 bExpBigger:
    if ( bExp == 0x7FFF ) {
        if ( bSig0 | bSig1 ) return propagateFloat128NaN( a, b );
        return packFloat128( zSign ^ 1, 0x7FFF, 0, 0 );
    }
    if ( aExp == 0 ) {
        ++expDiff;
    }
    else {
        aSig0 |= LIT64( 0x4000000000000000 );
    }
    shift128RightJamming( aSig0, aSig1, - expDiff, &aSig0, &aSig1);
    bSig0 |= LIT64( 0x4000000000000000 );
 bBigger:
    tmp = sub128( bSig0, bSig1, aSig0, aSig1);
	zSig0 = tmp.high, zSig1 = tmp.low;
    zExp = bExp;
    zSign ^= 1;
    goto normalizeRoundAndPack;
 aExpBigger:
    if ( aExp == 0x7FFF ) {
        if ( aSig0 | aSig1 ) return propagateFloat128NaN( a, b );
        return a;
    }
    if ( bExp == 0 ) {
        --expDiff;
    }
    else {
        bSig0 |= LIT64( 0x4000000000000000 );
    }
    shift128RightJamming( bSig0, bSig1, expDiff, &bSig0, &bSig1 );
    aSig0 |= LIT64( 0x4000000000000000 );
 aBigger:
    tmp = sub128( aSig0, aSig1, bSig0, bSig1);
    zSig0 = tmp.high, zSig1 = tmp.low;
    zExp = aExp;
 normalizeRoundAndPack:
    --zExp;
    return normalizeRoundAndPackFloat128( zSign, zExp - 14, zSig0, zSig1 );

}
/*----------------------------------------------------------------------------
| Returns the result of adding the quadruple-precision floating-point values
| `a' and `b'.  The operation is performed according to the IEC/IEEE Standard
| for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
float128_t _op_plus_float128_t_float128_t(float128_t a,float128_t b)
{
    int aSign, bSign;

    aSign = extractFloat128Sign( a );
    bSign = extractFloat128Sign( b );
    if ( aSign == bSign ) {
        return addFloat128Sigs( a, b, aSign );
    }
    else {
        return subFloat128Sigs( a, b, aSign );
    }
}
float128_t  _op_minus_float128_t_float128_t(float128_t a,float128_t b)
{
    int aSign, bSign;

    aSign = extractFloat128Sign( a );
    bSign = extractFloat128Sign( b );
    if ( aSign == bSign ) {
        return subFloat128Sigs( a, b, aSign );
    }
    else {
        return addFloat128Sigs( a, b, aSign );
    }

}
#else
float128_t _op_plus_float128_t_float128_t(float128_t a,float128_t b);
float128_t subFloat128Sigs( float128_t a, float128_t b, flag zSign );
float128_t addFloat128Sigs( float128_t a, float128_t b, flag zSign );
float128_t  _op_minus_float128_t_float128_t(float128_t a,float128_t b);
#endif

float128_t _op_plus_double_float128_t(double a,float128_t b)
{
	float128_t tmp = _op_cast_double_float128_t(a);
	return _op_plus_float128_t_float128_t(tmp,b);
}

float128_t _op_plus_long_double_float128_t(long double a,float128_t b)
{
	float128_t tmp = _op_cast_long_double_float128_t(a);
	return _op_plus_float128_t_float128_t(tmp,b);
}


float128_t _op_plus_int_float128_t(int a,float128_t b)
{
	float128_t tmp = _op_cast_int_float128_t(a);
	return _op_plus_float128_t_float128_t(tmp,b);
}
float128_t _op_plus_long_long_float128_t(long long a,float128_t b)
{
	float128_t tmp = _op_cast_long_long_float128_t(a);
	return _op_plus_float128_t_float128_t(tmp,b);
}

float128_t _op_plus_float128_t_double(float128_t a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	return _op_plus_float128_t_float128_t(a,tmp);
}

float128_t _op_plus_float128_t_int(float128_t a,int b)
{
	float128_t tmp = _op_cast_int_float128_t(b);
	return _op_plus_float128_t_float128_t(a,tmp);
}

float128_t _op_plus_float128_t_long_long(float128_t a,long long b)
{
	float128_t tmp = _op_cast_long_long_float128_t(b);
	return _op_plus_float128_t_float128_t(a,tmp);
}
float128_t _op_plus_float128_t_long_double(float128_t a,long double b)
{
	float128_t tmp = _op_cast_long_double_float128_t(b);
	return _op_plus_float128_t_float128_t(a,tmp);
}

float128_t  _op_minus_double_float128_t(double a,float128_t b)
{
	float128_t tmp = _op_cast_double_float128_t(a);
	return _op_minus_float128_t_float128_t(tmp,b);
}

float128_t  _op_minus_int_float128_t(int a,float128_t b)
{
	float128_t tmp = _op_cast_int_float128_t(a);
	return _op_minus_float128_t_float128_t(tmp,b);
}

float128_t  _op_minus_float128_t_int(float128_t a,int b)
{
	float128_t tmp = _op_cast_int_float128_t(b);
	return _op_minus_float128_t_float128_t(a,tmp);
}

float128_t  _op_minus_float128_t_double(float128_t a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	return _op_minus_float128_t_float128_t(a,tmp);
}

float128_t  _op_minus_float128_t_long_double(float128_t a,long double b)
{
	float128_t tmp = _op_cast_long_double_float128_t(b);
	return _op_minus_float128_t_float128_t(a,tmp);
}

float128_t *_op_minusasgn_pfloat128_t_int(float128_t *pf,int b)
{
	float128_t a = _op_cast_int_float128_t(b);
	*pf = _op_minus_float128_t_float128_t(*pf,a);
	return pf;
}

float128_t *_op_minusasgn_pfloat128_t_long_double(float128_t *pf,long double b)
{
	float128_t tmp = _op_cast_long_double_float128_t(b);
	*pf = _op_minus_float128_t_float128_t(*pf,tmp);
	return pf;
}
float128_t *_op_minusasgn_pfloat128_t_double(float128_t *pf, double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	*pf = _op_minus_float128_t_float128_t(*pf,tmp);
	return pf;
}


float128_t *_op_minusasgn_pfloat128_t_float128_t(float128_t *pf,float128_t b)
{
	*pf = _op_minus_float128_t_float128_t(*pf,b);
	return pf;
}
#ifndef ASSEMBLY_MACROS
/*----------------------------------------------------------------------------
| Returns the result of multiplying the quadruple-precision floating-point
| values `a' and `b'.  The operation is performed according to the IEC/IEEE
| Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
float128_t _op_multiply_float128_t_float128_t(float128_t a,float128_t b)
{
    int aSign, bSign, zSign;
    int32 aExp, bExp, zExp;
    bits64 aSig0, aSig1, bSig0, bSig1, zSig0, zSig1, zSig2, zSig3;
    float128_t z,tmp;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    bSig1 = extractFloat128Frac1( b );
    bSig0 = extractFloat128Frac0( b );
    bExp = extractFloat128Exp( b );
    bSign = extractFloat128Sign( b );
    zSign = aSign ^ bSign;
    if ( aExp == 0x7FFF ) {
        if (    ( aSig0 | aSig1 )
             || ( ( bExp == 0x7FFF ) && ( bSig0 | bSig1 ) ) ) {
            return propagateFloat128NaN( a, b );
        }
        if ( ( bExp | bSig0 | bSig1 ) == 0 ) goto invalid;
        return packFloat128( zSign, 0x7FFF, 0, 0 );
    }
    if ( bExp == 0x7FFF ) {
        if ( bSig0 | bSig1 ) return propagateFloat128NaN( a, b );
        if ( ( aExp | aSig0 | aSig1 ) == 0 ) {
 invalid:
            float_raise( float_flag_invalid );
            z.low = float128_default_nan_low;
            z.high = float128_default_nan_high;
            return z;
        }
        return packFloat128( zSign, 0x7FFF, 0, 0 );
    }
    if ( aExp == 0 ) {
        if ( ( aSig0 | aSig1 ) == 0 ) return packFloat128( zSign, 0, 0, 0 );
        normalizeFloat128Subnormal( aSig0, aSig1, &aExp, &aSig0, &aSig1 );
    }
    if ( bExp == 0 ) {
        if ( ( bSig0 | bSig1 ) == 0 ) return packFloat128( zSign, 0, 0, 0 );
        normalizeFloat128Subnormal( bSig0, bSig1, &bExp, &bSig0, &bSig1 );
    }
    zExp = aExp + bExp - 0x4000;
    aSig0 |= LIT64( 0x0001000000000000 );
    tmp = shortShift128Left16(bSig0, bSig1);
    mul128To256( aSig0, aSig1, tmp.high, tmp.low, &zSig0, &zSig1, &zSig2, &zSig3 );

    tmp = add128( zSig0, zSig1, aSig0, aSig1);

    zSig2 |= ( zSig3 != 0 );
    if ( LIT64( 0x0002000000000000 ) <= tmp.high ) {
        tmp = shift128ExtraRightJammingByOne(tmp.high, tmp.low, zSig2, &zSig2);
        ++zExp;
    }
    return roundAndPackFloat128( zSign, zExp, tmp.high, tmp.low, zSig2 );

}
#else
float128_t _op_multiply_float128_t_float128_t(float128_t a,float128_t b);
#endif
/* This had a parameter to specify a "payload" within the nan.  Dropped */
float128_t nanf128(void)
{
	float128_t z;
	float_raise( float_flag_invalid );
	z.low = float128_default_nan_low;
	z.high = float128_default_nan_high;
	return z;
}
#ifndef ASSEMBLY_MACROS
#define NOT_DIV_INLINING
/*----------------------------------------------------------------------------
| Returns the result of dividing the quadruple-precision floating-point value
| `a' by the corresponding value `b'.  The operation is performed according to
| the IEC/IEEE Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
float128_t _op_divide_float128_t_float128_t(float128_t a,float128_t b)
{
    flag aSign, bSign, zSign;
    int32 aExp, bExp, zExp,carry1;
    bits64 aSig0, aSig1, bSig0, bSig1, zSig0, zSig1, zSig2, z0, z1;
    bits64 rem0, rem1, rem2, rem3,  term2, term3;
    float128_t tmp;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    bSig1 = extractFloat128Frac1( b );
    bSig0 = extractFloat128Frac0( b );
    bExp = extractFloat128Exp( b );
    bSign = extractFloat128Sign( b );
    zSign = aSign ^ bSign;
    if (aExp == 0x7FFF ) goto specialConditionsA;
    if ( bExp == 0x7FFF ) goto specialConditionsB;
    if ( bExp == 0 ) goto specialConditionsC;
    if ( aExp == 0 ) goto specialConditionsD;
continuation:
    zExp = aExp - bExp + 0x3FFD;
    tmp = shortShift128Left(15,
        aSig0 | LIT64( 0x0001000000000000 ), aSig1);
    aSig0 = tmp.high, aSig1 = tmp.low;
    tmp = shortShift128Left(15,
        bSig0 | LIT64( 0x0001000000000000 ), bSig1);
    bSig0 = tmp.high, bSig1 = tmp.low;
    if ( le128( bSig0, bSig1, aSig0, aSig1 ) ) {
        shift128RightSmallCount( aSig0, aSig1, 1, &tmp);
        aSig0 = tmp.high, aSig1 = tmp.low;
        ++zExp;
    }
    if (bSig0 > aSig0)
        zSig0 = estimateDiv128To64( aSig0, aSig1, bSig0 );
    else
        zSig0 = (bits64)-1LL;
    tmp = mul128By64To192( bSig0, bSig1, zSig0, &term2 );
    /* inlining of this call */
    sub192( aSig0, aSig1, tmp.high, tmp.low, term2, &rem0, &rem1, &rem2 );
    while ( (sbits64) rem0 < 0 ) {
        --zSig0;
        /* This call is inlined below */
        z0 = rem2 + bSig1;
        carry1 = ( z0 < rem2 );
        rem2 = z0;
        z1 = rem1 + bSig0;
        z0 = rem0 + (z1 < rem1);
        rem1 = z1 + carry1;
        rem0 = z0 + ( rem1 < carry1 );

    }
    if (bSig0 > rem1)
        zSig1 = estimateDiv128To64( rem1, rem2, bSig0 );
    else
        zSig1 = (bits64)-1LL;
    if ( ( zSig1 & 0x3FFF ) <= 4 ) {
        tmp = mul128By64To192( bSig0, bSig1, zSig1, &term3 );
        sub192( rem1, rem2, tmp.high, tmp.low, term3, &rem1, &rem2, &rem3 );
        while ( (sbits64) rem1 < 0 ) {
            --zSig1;
            /* Inlined call to add192 */
            z0 = rem3 + bSig1;
            carry1 = ( z0 < rem3 );
            rem3 = z0;
            z1 = rem2 + bSig0;
            z0 = rem1 + (z1 < rem2);
            z1 += carry1;
            z0 += ( z1 < carry1 );
            rem2 = z1;
            rem1 = z0;
        }
        zSig1 |= ( ( rem1 | rem2 | rem3 ) != 0 );
    }
    tmp = shift128ExtraRightJammingSmallCount(15, zSig0, zSig1,  &zSig2 );
    return roundAndPackFloat128( zSign, zExp, tmp.high, tmp.low, zSig2 );
/*----------------------------------------------------------------
                Special conditions
---------------------------------------------------------------
*/
specialConditionsA: // We arrive here if aExp is 0x7fff
    if ( aSig0 | aSig1 ) return propagateFloat128NaN( a, b );
    if ( bExp == 0x7FFF ) {
        if ( bSig0 | bSig1 ) return propagateFloat128NaN( a, b );
        goto invalid;
    }
    return packFloat128( zSign, 0x7FFF, 0, 0 );
specialConditionsB: // We arrive here if bExp is 0x7fff
    if ( bSig0 | bSig1 ) return propagateFloat128NaN( a, b );
    return packFloat128( zSign, 0, 0, 0 );

specialConditionsC: // We arrive here if bExp is zero
    if ( ( bSig0 | bSig1 ) == 0 ) {
        if ( ( aExp | aSig0 | aSig1 ) == 0 ) {
     invalid:
             return nanf128();
        }
        float_raise( float_flag_divbyzero );
        return packFloat128( zSign, 0x7FFF, 0, 0 );
    }
    normalizeFloat128Subnormal( bSig0, bSig1, &bExp, &bSig0, &bSig1 );
    goto continuation;
specialConditionsD: // We arrive here if aExp is zero
    if ( ( aSig0 | aSig1 ) == 0 ) return packFloat128( zSign, 0, 0, 0 );
    normalizeFloat128Subnormal( aSig0, aSig1, &aExp, &aSig0, &aSig1 );
    goto continuation;
}
#else
float128_t _op_divide_float128_t_float128_t(float128_t a,float128_t b);
#endif

/*----------------------------------------------------------------------------
| Returns the remainder of the quadruple-precision floating-point value `a'
| with respect to the corresponding value `b'.  The operation is performed
| according to the IEC/IEEE Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/

float128_t float128_rem( float128_t a, float128_t b )
{
    flag aSign,  zSign;
    int32 aExp, bExp, expDiff;
    bits64 aSig0, aSig1, bSig0, bSig1, q, term0, term1, term2;
    bits64 allZero, alternateASig0, alternateASig1, sigMean1;
    sbits64 sigMean0;
    float128_t z,tmp;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    bSig1 = extractFloat128Frac1( b );
    bSig0 = extractFloat128Frac0( b );
    bExp = extractFloat128Exp( b );
    if ( aExp == 0x7FFF ) {
        if (    ( aSig0 | aSig1 )
             || ( ( bExp == 0x7FFF ) && ( bSig0 | bSig1 ) ) ) {
            return propagateFloat128NaN( a, b );
        }
        goto invalid;
    }
    if ( bExp == 0x7FFF ) {
        if ( bSig0 | bSig1 ) return propagateFloat128NaN( a, b );
        return a;
    }
    if ( bExp == 0 ) {
        if ( ( bSig0 | bSig1 ) == 0 ) {
 invalid:
            float_raise( float_flag_invalid );
            z.low = float128_default_nan_low;
            z.high = float128_default_nan_high;
            return z;
        }
        normalizeFloat128Subnormal( bSig0, bSig1, &bExp, &bSig0, &bSig1 );
    }
    if ( aExp == 0 ) {
        if ( ( aSig0 | aSig1 ) == 0 ) return a;
        normalizeFloat128Subnormal( aSig0, aSig1, &aExp, &aSig0, &aSig1 );
    }
    expDiff = aExp - bExp;
    if ( expDiff < -1 ) return a;
    tmp = shortShift128Left(15 - ( expDiff < 0 ),
        aSig0 | LIT64( 0x0001000000000000 ),
        aSig1);
	aSig0 = tmp.high, aSig1 = tmp.low;
    tmp = shortShift128Left(15,
        bSig0 | LIT64( 0x0001000000000000 ), bSig1);
	bSig0 = tmp.high, bSig1 = tmp.low;
    q = le128( bSig0, bSig1, aSig0, aSig1 );
    if ( q ) {
		tmp = sub128( aSig0, aSig1, bSig0, bSig1);
		aSig0 = tmp.high, aSig1 = tmp.low;
	}
    expDiff -= 64;
    while ( 0 < expDiff ) {
        if (bSig0 > aSig0)
            q = estimateDiv128To64( aSig0, aSig1, bSig0 );
        else
            q = (bits64)-1LL;
        q = ( 4 < q ) ? q - 4 : 0;
        tmp = mul128By64To192( bSig0, bSig1, q, &term2 );
		term0 = tmp.high, term1 = tmp.low;
        shortShift192Left( term0, term1, term2, 61, &term1, &term2, &allZero );
        tmp = shortShift128Left(61, aSig0, aSig1);
        aSig0 = tmp.high, aSig1 = tmp.low;
        tmp = sub128( aSig0, 0, term1, term2);
		aSig0 = tmp.high, aSig1 = tmp.low;
        expDiff -= 61;
    }
    if ( -64 < expDiff ) {
        if (bSig0 > aSig0)
            q = estimateDiv128To64( aSig0, aSig1, bSig0 );
        else
            q = (bits64)-1LL;
        q = ( 4 < q ) ? q - 4 : 0;
        q >>= - expDiff;
        shift128RightSmallCount( bSig0, bSig1, 12, &tmp);
        bSig0 = tmp.high, bSig1 = tmp.low;
        expDiff += 52;
        if ( expDiff < 0 ) {
            tmp = shift128Right( aSig0, aSig1, - expDiff);
            aSig0 = tmp.high, aSig1 = tmp.low;
        }
        else if (expDiff > 0) {
            tmp = shortShift128Left(expDiff, aSig0, aSig1);
            aSig0 = tmp.high, aSig1 = tmp.low;
        }
        tmp = mul128By64To192( bSig0, bSig1, q, &term2 );
		term0 = tmp.high, term1 = tmp.low;
        tmp = sub128( aSig0, aSig1, term1, term2 );
        aSig0 = tmp.high, aSig1 = tmp.low;
    }
    else {
        shift128RightSmallCount( aSig0, aSig1, 12, &tmp);
        aSig0 = tmp.high, aSig1 = tmp.low;
        shift128RightSmallCount( bSig0, bSig1, 12, &tmp);
        bSig0 = tmp.high, bSig1 = tmp.low;
    }
    do {
        alternateASig0 = aSig0;
        alternateASig1 = aSig1;
        ++q;
        tmp = sub128( aSig0, aSig1, bSig0, bSig1);
		aSig0 = tmp.high, aSig1 = tmp.low;
    } while ( 0 <= (sbits64) aSig0 );
    tmp = add128( aSig0, aSig1, alternateASig0, alternateASig1);
    sigMean0 = tmp.high, sigMean1 = tmp.low;

    if (    ( sigMean0 < 0 )
         || ( ( ( sigMean0 | sigMean1 ) == 0 ) && ( q & 1 ) ) ) {
        aSig0 = alternateASig0;
        aSig1 = alternateASig1;
    }
    zSign = ( (sbits64) aSig0 < 0 );
    if ( zSign ) {
		tmp = sub128( 0, 0, aSig0, aSig1);
		aSig0 = tmp.high, aSig1 = tmp.low;
	}
    return
        normalizeRoundAndPackFloat128( aSign ^ zSign, bExp - 4, aSig0, aSig1 );

}

/*----------------------------------------------------------------------------
| Returns the square root of the quadruple-precision floating-point value `a'.
| The operation is performed according to the IEC/IEEE Standard for Binary
| Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
float128_t sqrtf128(float128_t a)
{
    flag aSign;
    int32 aExp, zExp;
    bits64 aSig0, aSig1, zSig0, zSig1, zSig2, doubleZSig0;
    bits64 rem0, rem1, rem2, rem3, term0, term1, term2, term3;
    float128_t z, tmp;

    aSig1 = extractFloat128Frac1( a );
    aSig0 = extractFloat128Frac0( a );
    aExp = extractFloat128Exp( a );
    aSign = extractFloat128Sign( a );
    if ( aExp == 0x7FFF ) {
        if ( aSig0 | aSig1 ) return propagateFloat128NaN( a, a );
        if ( ! aSign ) return a;
        goto invalid;
    }
    if ( aSign ) {
        if ( ( aExp | aSig0 | aSig1 ) == 0 ) return a;
 invalid:
        float_raise( float_flag_invalid );
        z.low = float128_default_nan_low;
        z.high = float128_default_nan_high;
        return z;
    }
    if ( aExp == 0 ) {
        if ( ( aSig0 | aSig1 ) == 0 ) return packFloat128Zero();
        normalizeFloat128Subnormal( aSig0, aSig1, &aExp, &aSig0, &aSig1 );
    }
    zExp = ( ( aExp - 0x3FFF )>>1 ) + 0x3FFE;
    aSig0 |= LIT64( 0x0001000000000000 );
    zSig0 = estimateSqrt32( aExp, aSig0>>17 );
    tmp = shortShift128Left(13 - ( aExp & 1 ), aSig0, aSig1);
    aSig0 = tmp.high, aSig1 = tmp.low;
    if ((zSig0 << 32) > aSig0)
        zSig0 = estimateDiv128To64( aSig0, aSig1, zSig0<<32 ) + ( zSig0<<30 );
    else
        zSig0 = (bits64)-1LL +(zSig0<<30);
    doubleZSig0 = zSig0<<1;
    mul64To128( zSig0, zSig0, &term0, &term1 );
    tmp = sub128( aSig0, aSig1, term0, term1);
	rem0 = tmp.high, rem1 = tmp.low;
    while ( (sbits64) rem0 < 0 ) {
        --zSig0;
        doubleZSig0 -= 2;
        tmp = add128( rem0, rem1, zSig0>>63, doubleZSig0 | 1);
        rem0 = tmp.high, rem1 = tmp.low;
    }
    if (doubleZSig0 > rem1)
        zSig1 = estimateDiv128To64( rem1, 0, doubleZSig0 );
    else
        zSig1 = (bits64)-1LL;
    if ( ( zSig1 & 0x1FFF ) <= 5 ) {
        if ( zSig1 == 0 ) zSig1 = 1;
        mul64To128( doubleZSig0, zSig1, &term1, &term2 );
        tmp = sub128( rem1, 0, term1, term2);
		rem1 = tmp.high, rem2 = tmp.low;
        mul64To128( zSig1, zSig1, &term2, &term3 );
        sub192( rem1, rem2, 0, term2, term3, &rem1, &rem2, &rem3 );
        while ( (sbits64) rem1 < 0 ) {
            --zSig1;
            tmp = shortShift128Left(1, 0, zSig1);
			term2 = tmp.high, term3 = tmp.low;
            term3 |= 1;
            term2 |= doubleZSig0;
            add192( rem1, rem2, rem3, term2, term3, &rem1, &rem2, &rem3 );
        }
        zSig1 |= ( ( rem1 | rem2 | rem3 ) != 0 );
    }
    tmp = shift128ExtraRightJammingSmallCount(14, zSig0, zSig1,  &zSig2 );
    return roundAndPackFloat128( 0, zExp, tmp.high, tmp.low, zSig2 );

}

/*----------------------------------------------------------------------------
| Returns 1 if the quadruple-precision floating-point value `a' is equal to
| the corresponding value `b', and 0 otherwise.  The comparison is performed
| according to the IEC/IEEE Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
int _op_equal_float128_t_float128_t(float128_t a,float128_t b)
{

    if (    (    ( extractFloat128Exp( a ) == 0x7FFF )
              && ( extractFloat128Frac0( a ) | extractFloat128Frac1( a ) ) )
         || (    ( extractFloat128Exp( b ) == 0x7FFF )
              && ( extractFloat128Frac0( b ) | extractFloat128Frac1( b ) ) )
       ) {
        if (    float128_is_signaling_nan( a )
             || float128_is_signaling_nan( b ) ) {
            float_raise( float_flag_invalid );
        }
        return 0;
    }
    return
           ( a.low == b.low )
        && (    ( a.high == b.high )
             || (    ( a.low == 0 )
                  && ( (bits64) ( ( a.high | b.high )<<1 ) == 0 ) )
           );

}

int _op_equal_double_float128_t(double a,float128_t b)
{
	float128_t tmp = _op_cast_double_float128_t(a);
	return _op_equal_float128_t_float128_t(tmp,b);
}

int _op_equal_float128_t_int(int a,float128_t b)
{
	float128_t tmp = _op_cast_int_float128_t(a);
	return _op_equal_float128_t_float128_t(tmp,b);
}

int _op_equal_long_long_float128_t(long long a,float128_t b)
{
	float128_t tmp = _op_cast_long_long_float128_t(a);
	return _op_equal_float128_t_float128_t(tmp,b);
}


int _op_equal_float128_t_double(float128_t a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	return _op_equal_float128_t_float128_t(a,tmp);
}

int _op_notEqual_float128_t_float128_t(float128_t a,float128_t b)
{
	return !_op_equal_float128_t_float128_t(a,b);
}

int _op_notEqual_float128_t_double(float128_t a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	return !_op_equal_float128_t_float128_t(a,tmp);
}

int _op_notEqual_double_float128_t(double a,float128_t b)
{
	float128_t tmp = _op_cast_double_float128_t(a);
	return !_op_equal_float128_t_float128_t(tmp,b);
}

int _op_notEqual_int_float128_t(int a,float128_t b)
{
	float128_t tmp = _op_cast_double_float128_t(a);
	return !_op_equal_float128_t_float128_t(tmp,b);
}

int _op_notEqual_float128_t_int(float128_t a,int b)
{
	float128_t tmp = _op_cast_int_float128_t(b);
	return !_op_equal_float128_t_float128_t(a,tmp);
}

/*----------------------------------------------------------------------------
| Returns 1 if the quadruple-precision floating-point value `a' is less than
| or equal to the corresponding value `b', and 0 otherwise.  The comparison
| is performed according to the IEC/IEEE Standard for Binary Floating-Point
| Arithmetic.
*----------------------------------------------------------------------------*/
int  _op_lessequal_float128_t_float128_t(float128_t a,float128_t b)
{
    int aSign, bSign;

    if (    (    ( extractFloat128Exp( a ) == 0x7FFF )
              && ( extractFloat128Frac0( a ) | extractFloat128Frac1( a ) ) )
         || (    ( extractFloat128Exp( b ) == 0x7FFF )
              && ( extractFloat128Frac0( b ) | extractFloat128Frac1( b ) ) )
       ) {
        float_raise( float_flag_invalid );
        return 0;
    }
    aSign = extractFloat128Sign( a );
    bSign = extractFloat128Sign( b );
    if ( aSign != bSign ) {
        return
               aSign
            || (    ( ( (bits64) ( ( a.high | b.high )<<1 ) ) | a.low | b.low )
                 == 0 );
    }
    return
          aSign ? le128( b.high, b.low, a.high, a.low )
        : le128( a.high, a.low, b.high, b.low );

}
int  _op_greater_float128_t_float128_t(float128_t a,float128_t b)
{
	return !_op_lessequal_float128_t_float128_t(a,b);
}
/*----------------------------------------------------------------------------
| Returns 1 if the quadruple-precision floating-point value `a' is less than
| the corresponding value `b', and 0 otherwise.  The comparison is performed
| according to the IEC/IEEE Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/
int _op_less_float128_t_float128_t(float128_t a,float128_t b)
{
    int aSign, bSign;

    if (    (    ( extractFloat128Exp( a ) == 0x7FFF )
              && ( extractFloat128Frac0( a ) | extractFloat128Frac1( a ) ) )
         || (    ( extractFloat128Exp( b ) == 0x7FFF )
              && ( extractFloat128Frac0( b ) | extractFloat128Frac1( b ) ) )
       ) {
        float_raise( float_flag_invalid );
        return 0;
    }
    aSign = extractFloat128Sign( a );
    bSign = extractFloat128Sign( b );
    if ( aSign != bSign ) {
        return
               aSign
            && (    ( ( (bits64) ( ( a.high | b.high )<<1 ) ) | a.low | b.low )
                 != 0 );
    }
    return
          aSign ? lt128( b.high, b.low, a.high, a.low )
        : lt128( a.high, a.low, b.high, b.low );

}

int _op_less_float128_t_int(float128_t a,int b)
{
	float128_t tmp = _op_cast_int_float128_t(b);
	return _op_less_float128_t_float128_t(a,tmp);
}

int _op_less_float128_t_long_double(float128_t a,long double b)
{
	float128_t tmp = _op_cast_long_double_float128_t(b);
	return _op_less_float128_t_float128_t(a,tmp);
}

int _op_less_float128_t_double(float128_t a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	return _op_less_float128_t_float128_t(a,tmp);
}


int _op_greaterequal_float128_t_float128_t(float128_t a,float128_t b)
{
	return !_op_less_float128_t_float128_t(a,b);
}

int _op_greaterequal_float128_t_double(float128_t a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	return !_op_less_float128_t_float128_t(a,tmp);
}

int _op_greaterequal_float128_t_int(float128_t a,int b)
{
	float128_t tmp = _op_cast_int_float128_t(b);
	return !_op_less_float128_t_float128_t(a,tmp);
}

int _op_greaterequal_float128_t_long_long(float128_t a,long long b)
{
	float128_t tmp = _op_cast_long_long_float128_t(b);
	return !_op_less_float128_t_float128_t(a,tmp);
}

int _op_greaterequal_double_float128_t(double a,float128_t b)
{
	float128_t tmp = _op_cast_double_float128_t(a);
	return !_op_less_float128_t_float128_t(tmp,b);
}

int _op_greaterequal_int_float128_t(int a,float128_t b)
{
	float128_t tmp = _op_cast_int_float128_t(a);
	return !_op_less_float128_t_float128_t(tmp,b);
}

int _op_greaterequal_long_double_float128_t(long double a,float128_t b)
{
	float128_t tmp = _op_cast_long_double_float128_t(a);
	return !_op_less_float128_t_float128_t(tmp,b);
}

int _op_greaterequal_long_long_float128_t(long long a,float128_t b)
{
	float128_t tmp = _op_cast_long_long_float128_t(a);
	return !_op_less_float128_t_float128_t(tmp,b);
}


/*----------------------------------------------------------------------------
| Returns 1 if the quadruple-precision floating-point value `a' is equal to
| the corresponding value `b', and 0 otherwise.  The invalid exception is
| raised if either operand is a NaN.  Otherwise, the comparison is performed
| according to the IEC/IEEE Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/

flag float128_eq_signaling( float128_t a, float128_t b )
{

    if (    (    ( extractFloat128Exp( a ) == 0x7FFF )
              && ( extractFloat128Frac0( a ) | extractFloat128Frac1( a ) ) )
         || (    ( extractFloat128Exp( b ) == 0x7FFF )
              && ( extractFloat128Frac0( b ) | extractFloat128Frac1( b ) ) )
       ) {
        float_raise( float_flag_invalid );
        return 0;
    }
    return
           ( a.low == b.low )
        && (    ( a.high == b.high )
             || (    ( a.low == 0 )
                  && ( (bits64) ( ( a.high | b.high )<<1 ) == 0 ) )
           );

}

/*----------------------------------------------------------------------------
| Returns 1 if the quadruple-precision floating-point value `a' is less than
| or equal to the corresponding value `b', and 0 otherwise.  Quiet NaNs do not
| cause an exception.  Otherwise, the comparison is performed according to the
| IEC/IEEE Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/

flag float128_le_quiet( float128_t a, float128_t b )
{
    flag aSign, bSign;

    if (    (    ( extractFloat128Exp( a ) == 0x7FFF )
              && ( extractFloat128Frac0( a ) | extractFloat128Frac1( a ) ) )
         || (    ( extractFloat128Exp( b ) == 0x7FFF )
              && ( extractFloat128Frac0( b ) | extractFloat128Frac1( b ) ) )
       ) {
        if (    float128_is_signaling_nan( a )
             || float128_is_signaling_nan( b ) ) {
            float_raise( float_flag_invalid );
        }
        return 0;
    }
    aSign = extractFloat128Sign( a );
    bSign = extractFloat128Sign( b );
    if ( aSign != bSign ) {
        return
               aSign
            || (    ( ( (bits64) ( ( a.high | b.high )<<1 ) ) | a.low | b.low )
                 == 0 );
    }
    return
          aSign ? le128( b.high, b.low, a.high, a.low )
        : le128( a.high, a.low, b.high, b.low );

}

/*----------------------------------------------------------------------------
| Returns 1 if the quadruple-precision floating-point value `a' is less than
| the corresponding value `b', and 0 otherwise.  Quiet NaNs do not cause an
| exception.  Otherwise, the comparison is performed according to the IEC/IEEE
| Standard for Binary Floating-Point Arithmetic.
*----------------------------------------------------------------------------*/

flag float128_lt_quiet( float128_t a, float128_t b )
{
    flag aSign, bSign;

    if (    (    ( extractFloat128Exp( a ) == 0x7FFF )
              && ( extractFloat128Frac0( a ) | extractFloat128Frac1( a ) ) )
         || (    ( extractFloat128Exp( b ) == 0x7FFF )
              && ( extractFloat128Frac0( b ) | extractFloat128Frac1( b ) ) )
       ) {
        if (    float128_is_signaling_nan( a )
             || float128_is_signaling_nan( b ) ) {
            float_raise( float_flag_invalid );
        }
        return 0;
    }
    aSign = extractFloat128Sign( a );
    bSign = extractFloat128Sign( b );
    if ( aSign != bSign ) {
        return
               aSign
            && (    ( ( (bits64) ( ( a.high | b.high )<<1 ) ) | a.low | b.low )
                 != 0 );
    }
    return
          aSign ? lt128( b.high, b.low, a.high, a.low )
        : lt128( a.high, a.low, b.high, b.low );

}
float128_t _op_multiply_float128_t_int(float128_t a,int b)
{
	float128_t tmp = _op_cast_int_float128_t(b);
	return _op_multiply_float128_t_float128_t(tmp,a);
}

float128_t _op_multiply_float128_t_double(float128_t a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	return _op_multiply_float128_t_float128_t(tmp,a);
}

float128_t _op_multiply_float128_t_long_double(float128_t a,long double b)
{
	float128_t tmp = _op_cast_long_double_float128_t(b);
	return _op_multiply_float128_t_float128_t(tmp,a);
}

float128_t _op_multiply_long_double_float128_t(long double a,float128_t b)
{
	float128_t tmp = _op_cast_long_double_float128_t(a);
	return _op_multiply_float128_t_float128_t(tmp,b);
}

float128_t _op_multiply_int_float128_t(int a,float128_t b)
{
	float128_t tmp = _op_cast_int_float128_t(a);
	return _op_multiply_float128_t_float128_t(tmp,b);
}


float128_t _op_multiply_double_float128_t(double a,float128_t b)
{
	float128_t tmp = _op_cast_double_float128_t(a);
	return _op_multiply_float128_t_float128_t(tmp,b);
}

float128_t _op_divide_float128_t_int(float128_t a,int b)
{
	float128_t tmp = _op_cast_int_float128_t(b);
	return _op_divide_float128_t_float128_t(a,tmp);
}
float128_t _op_divide_int_float128_t(int a,float128_t b)
{
	float128_t tmp = _op_cast_int_float128_t(a);
	return _op_divide_float128_t_float128_t(tmp,b);
}
float128_t _op_divide_float128_t_double(float128_t a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	return _op_divide_float128_t_float128_t(a,tmp);
}

float128_t _op_divide_float128_t_long_double(float128_t a,long double b)
{
	float128_t tmp = _op_cast_long_double_float128_t(b);
	return _op_divide_float128_t_float128_t(a,tmp);
}


float128_t _op_divide_double_float128_t(double a,float128_t b)
{
	float128_t tmp = _op_cast_double_float128_t(a);
	return _op_divide_float128_t_float128_t(tmp,b);
}

float128_t _op_divide_long_double_float128_t(long double a,float128_t b)
{
	float128_t tmp = _op_cast_long_double_float128_t(a);
	return _op_divide_float128_t_float128_t(tmp,b);
}


float128_t *_op_asgn_pfloat128_t_double(float128_t *pfloat128_t,double x)
{
	float128_t r = _op_cast_long_double_float128_t((long double)x);
	*pfloat128_t = r;
	return pfloat128_t;
}

float128_t *_op_asgn_pfloat128_t_float(float128_t *pfloat128_t,float x)
{
	float128_t r = _op_cast_long_double_float128_t((long double)x);
	*pfloat128_t = r;
	return pfloat128_t;
}


long double *_op_asgn_plong_double_float128_t(long double *plong_double,float128_t x)
{
	long double r = _op_cast_float128_t_long_double(x);
	*plong_double = r;
	return plong_double;
}

int *_op_asgn_pint_float128_t(int *pint,float128_t x)
{
	int r = _op_cast_float128_t_int(x);
	*pint = r;
	return pint;
}

double *_op_asgn_pdouble_float128_t(double *pdouble,float128_t x)
{
	double r = _op_cast_float128_t_double(x);
	*pdouble = r;
	return pdouble;
}


float128_t _op_minus_long_double_float128_t(long double a,float128_t b)
{
	float128_t tmp = _op_cast_long_double_float128_t(a);
	return  _op_minus_float128_t_float128_t(tmp,b);
}

int _op_greater_float128_t_long_double(float128_t a,long double b)
{
	float128_t tmp = _op_cast_long_double_float128_t(b);
	return  _op_greater_float128_t_float128_t(a,tmp);
}

int _op_greater_float128_t_int(float128_t a,int b)
{
	float128_t tmp = _op_cast_int_float128_t(b);
	return  _op_greater_float128_t_float128_t(a,tmp);
}


int _op_greater_float128_t_double(float128_t a,double b)
{
	return  _op_greater_float128_t_float128_t(a,
	                _op_cast_double_float128_t(b));
}

float128_t _op_minus_float128_t(float128_t a)
{
	float128_t tmp = a;
#if 0
	return _op_minus_float128_t_float128_t(tmp,a);
#else
	if (tmp.high&0x8000000000000000ULL) {
		tmp.high &= ~0x8000000000000000ULL;
	}
	else tmp.high |= 0x8000000000000000ULL;
	return tmp;
#endif
}

float128_t *_op_asgn_pfloat128_t_long_double(float128_t *a,long double b)
{
	float128_t tmp = _op_cast_long_double_float128_t(b);
	*a = tmp;
	return a;
}
float128_t *_op_asgn_pfloat128_t_int(float128_t *a,int b)
{
	*a = _op_cast_int_float128_t(b);
	return a;
}

float128_t *_op_plusasgn_pfloat128_t_float128_t(float128_t *a,float128_t b)
{
	*a = _op_plus_float128_t_float128_t(*a,b);
	return a;
}

float128_t *_op_plusasgn_pfloat128_t_int(float128_t *a,int b)
{
	float128_t tmp = _op_cast_int_float128_t(b);
	tmp = _op_plus_float128_t_float128_t(tmp,*a);
	*a = tmp;
	return a;
}

int *_op_plusasgn_pint_float128_t(int *a,float128_t b)
{
	int tmp = _op_cast_float128_t_int(b);
	*a += tmp;
	return a;
}

float128_t *_op_plusasgn_pfloat128_t_double(float128_t *a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	tmp = _op_plus_float128_t_float128_t(tmp,*a);
	*a = tmp;
	return a;
}


float128_t *_op_plusasgn_pfloat128_t_long_double(float128_t *a,long double b)
{
	float128_t tmp = _op_cast_long_double_float128_t(b);
	tmp = _op_plus_float128_t_float128_t(*a,tmp);
	*a = tmp;
	return a;
}


float128_t *_op_multasgn_pfloat128_t_long_double(float128_t *a,long double b)
{
	float128_t tmp = _op_cast_long_double_float128_t(b);
	tmp = _op_multiply_float128_t_float128_t(tmp,*a);
	*a = tmp;
	return a;
}

float128_t *_op_multasgn_pfloat128_t_float128_t(float128_t *a,float128_t b)
{
	*a = _op_multiply_float128_t_float128_t(b,*a);
	return a;
}

float128_t *_op_multasgn_pfloat128_t_double(float128_t *a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	tmp = _op_multiply_float128_t_float128_t(tmp,*a);
	*a = tmp;
	return a;
}

float128_t *_op_divasgn_pfloat128_t_double(float128_t *a,double b)
{
	float128_t tmp = _op_cast_double_float128_t(b);
	*a = _op_divide_float128_t_float128_t(*a,tmp);
	return a;
}

float128_t *_op_divasgn_pfloat128_t_float128_t(float128_t *a,float128_t b)
{
	*a = _op_divide_float128_t_float128_t(*a,b);
	return a;
}

float128_t *_op_multasgn_pfloat128_t_int(float128_t *a,int b)
{
	float128_t tmp = _op_cast_int_float128_t(b);
	tmp = _op_multiply_float128_t_float128_t(tmp,*a);
	*a = tmp;
	return a;
}

float128_t fabsf128(float128_t x)
{
	float128_t r = x;
	r.high &= 0x7fffffffffffffffULL;
	return r;
}

#ifdef __LCC__
qfloat _op_cast_float128_t_qfloat(float128_t a)
{
	qfloat r;
	int exp = extractFloat128Exp(a);
	int sign = extractFloat128Sign(a);
	long long low = a.low<<15;
	long long high;

	if (a.high|a.low) {
		high  = ((a.high << 15)|0x8000000000000000ULL);
		high |= (((a.high&1)<<15)|((a.low >> 49)&0xffff));
	}
	else {
		r.sign = r.exponent = 0;
		r.mantissa[0] = r.mantissa[1] = r.mantissa[2] = 0;
		r.mantissa[3] = r.mantissa[4] = r.mantissa[5] = 0;
		r.mantissa[6] = 0;
		return r;
	}
	exp -= 0x3FFF;
	exp += 0x80001;
	r.sign = sign?-1:0;
	r.exponent = exp;
	r.mantissa[0] = high;
	r.mantissa[1] = low;
	return r;
}
#endif

float128_t frexpf128(float128_t x,int *m)
{
	float128_t result;
	int sign = extractFloat128Sign(x);
	int exp = extractFloat128Exp(x);
	long long f0 = extractFloat128Frac0(x);
	long long f1 = extractFloat128Frac1(x);
	int new_exp;
	if ((x.high|x.low)==0) { *m = 0; result.high = result.low = 0; return result; }
	new_exp = 0x3fff - 1;
	result = packFloat128(sign, new_exp,f0,f1);
	*m = exp-new_exp ;
	return result;
}
/* ______________________________________________________________________________
//
                                     C O N S T A N T S
*/
// 0.6931471805599453094172321214581765680755001343602552541
float128_t __xLn2 = {0xf35793c7673007e5ULL, 0x3ffe62e42fefa39eULL};
//1.41421356237309504880168872420969807856967187537694807317667973
float128_t __xSqrt2 = {0xc908b2fb1366ea95ULL, 0x3fff6a09e667f3bcULL};
float128_t __float128Exp2Min = {0,0xc00cfff600000000ULL}; // -16382.75
float128_t __float128Exp2Max = {0,0x400cfff600000000ULL}; //  16382.75
// 1 /log(2) --> 1.442695040888963407359924681001892137426645954152985934135449
float128_t __xLog2e = {0xe1777d0ffda0d23a,0x3fff71547652b82f};
// 2.30258509299404568401799145468436420760110148862877297
float128_t __xLog2_10 = {0x582dd4adac5705a6,0x400026bb1bbb5551};
float128_t __xOne = {0,0x3fff000000000000ULL};
// 1.1897314953572317650857593266280070e+4932
float128_t __MAXNUMf128 = {0xffffffffffffffffULL, 0x7ffeffffffffffffULL };
// 3.3621031431120935062626778173217526e-4932
float128_t __MINNUMf128 = {0,0x1000000000000ULL };
// __MAXLOGf128 = 1.135583025911358400418251384584930671458833e4L;
float128_t __MAXLOGf128 = {0xF35793C7673007E6ULL, 0x400C62E42FEFA39EULL };
// MINLOG = -1.143276959615573793352782661133116431383730e4L;
float128_t __MINLOGf128 = {0x2C89D24D65E96273ULL, 0xC00C654628220780ULL};
// 0.5
float128_t __xHalf = {0, 0x3ffe000000000000};
// PI
float128_t __M_PIF128 = {0x8469898CC51701B8, 0x4000921FB54442D1} ;
float128_t __xZero = {0};
// PI/4 0.785398163397448309615660845819875721049292349843776
float128_t __xPiO4 = {0x8469898CC51701B8, 0x3FFE921FB54442D1};
// 1 / sqrt(pi)
float128_t __xOneOSqrtPi = {0xD11AE3A914FED7FEULL,0x3FFE20DD750429B6ULL};
// sqrt(2)/2 0.70710678118654752440084436210484903928483593
float128_t __xSqrtH = { 0xC908B2FB1366EA95ULL,0x3FFE6A09E667F3BCULL};
// 1/sqrt(pi)
//float128_t __xOneOSqrtPi = {0xD11AE3A914FED7FEULL, 0x3FFE20DD750429B6ULL};
/* __________________________________________________________________________________ */
float128_t inversef128(float128_t in)
{
	return _op_divide_float128_t_float128_t(__xOne,in);
}
#ifndef __LCC__
float128_t logf128(float128_t z)
{
	float128_t f, h, a, b, xOne;
	int k, m, exp;

	if ((z.low|z.high)==0) {
		float_exception_flags |= float_flag_divbyzero;
		// return minus infinity
		h.low = 0;
		h.high = 0xffff000000000000ULL;
		return h;
	}
	if (z.high&0x8000000000000000ULL )  {
		float_exception_flags |=float_flag_invalid;
		h.low = float128_default_nan_low;
		h.high = float128_default_nan_high;
		return h;
	}
	xOne.high =0x3fff000000000000ULL;
	xOne.low = 0;
	if (z.high == xOne.high&& z.low == 0) {
		h.high = h.low = 0;
		return h;
	}
	z = frexpf128(z, &m);

	z = _op_multiply_float128_t_float128_t(z, __xSqrt2);
	a = _op_minus_float128_t_float128_t(z, xOne);
	b = _op_plus_float128_t_float128_t(z, xOne);
	z = _op_divide_float128_t_float128_t(a, b);
        //h = xpr2 (z, 1);
	h = _op_multiply_float128_t_int(z,2);
	z = _op_multiply_float128_t_float128_t(z, z);
	f = h;
	k = 1;
	exp = extractFloat128Exp(h);
	for (; (exp-0x3fff) > -112;) {
		h = _op_multiply_float128_t_float128_t(h, z);
		exp = extractFloat128Exp(h);
		a = _op_divide_float128_t_int(h, k += 2);
		f = _op_plus_float128_t_float128_t(f,a);
	}
	a = _op_cast_double_float128_t(m-0.5);
	b = _op_multiply_float128_t_float128_t(a,__xLn2);
	return _op_plus_float128_t_float128_t(f,b);
}
#endif

#if 0
// Inverse of log 10 1/log(10= 0.434294481903251827651128918916605082294397005803666566114453783165
static float128_t invLog10={0xE32A6AB7555F5A68ULL, 0x3FFDBCB7B1526E50ULL};
extern float128_t logf128(float128_t);
float128_t log10f128(float128_t x)
{
	float128_t result = logf128(x);
	return _op_multiply_float128_t_float128_t(result,invLog10);
}
#endif
static void ClearBlock(void* destination, int bitCount)
{
	unsigned char* dst = destination;

	if (bitCount >= 64) {
		long long *pll = (long long *)dst;
		*pll = 0;
		dst += 8;
		bitCount -= 64;
	}
	if (bitCount >= 32) {
		int *pi = (int *)dst;
		*pi = 0;
		dst += 4;
		bitCount -= 32;
	}
	if (bitCount >= 16) {
		short *ps = (short *)dst;
		*ps = 0;
		dst += 2;
		bitCount -= 16;
	}
	if (bitCount >= 8) {
		*dst++ = 0;
		bitCount -= 8;
	}
	//do trailing bits
	if (bitCount > 0) {
		*dst &= ~((1 << bitCount) - 1);
	}

}


float128_t floorf128(float128_t x)
{
	float128_t result;
	int unbiasedExponent = extractFloat128Exp(x) - 0x3fff;
	if (unbiasedExponent == 0x7fff)
		result= x; // NAN Inf
	else if (unbiasedExponent < 0) {
		if (extractFloat128Sign(x)) {
			result = _op_minus_float128_t(__xOne);
		}
		else result.high = result.low = 0;

	}
	else {
		result = x;
		if (unbiasedExponent < 112)
			ClearBlock(&result,112-unbiasedExponent);
		if (extractFloat128Sign(x))
			result = _op_minus_float128_t_float128_t(result,__xOne);
	}
	return result;
}

float128_t modff128(float128_t x,float128_t *p)
{
	float128_t t1,intPart = floorf128(x);
	if (extractFloat128Sign(x)) {
		intPart = _op_plus_float128_t_int(intPart,1);
	}
	*p = intPart;
	t1= _op_minus_float128_t_float128_t(x,intPart);
	return t1;
}

float128_t ldexpf128(float128_t x,int p)
{
	float128_t result;
//	int sign = extractFloat128Sign(x);
	int exp = extractFloat128Exp(x);
	// ACtually, w0 is not used, it is better not to extract it!
	//bits64 w0 = extractFloat128Frac0(x);
	result.low =  extractFloat128Frac1(x);
	if (exp == 0) return x;
	result.high = (x.high&0x8000FFFFFFFFFFFFULL)|((long long)exp+p)<<48;
	return result;
	//return normalizeRoundAndPackFloat128(sign,exp+p+1,w0,w1);
}
#if 0
static float128_t exp2f128(float128_t x)
{
	float128_t result,d,t1,t2,s,f;
	int k;
	int sign,m;
	int exp;
	long long w0;
	long long w1;

	if (_op_less_float128_t_float128_t(x,__float128Exp2Min)) {
		result.high = result.low = 0;
		return result;
	}
	if (_op_greater_float128_t_float128_t(x,__float128Exp2Max)) {
		return packFloat128(0,0x7FFE,
                        LIT64( 0x0000FFFFFFFFFFFF ),
                        LIT64( 0xFFFFFFFFFFFFFFFF )
                    );
	}
	sign = extractFloat128Sign(x);
	x = modff128(x,&t1);
	k = _op_cast_float128_t_int(t1);
//	if (sign)
//		k = -k;
	x = _op_multiply_float128_t_float128_t(x,__xLn2);
	x = _op_divide_float128_t_int(x,2);
	s = _op_multiply_float128_t_float128_t(x,x);
	f.high = f.low = 0;
	m = 21;
	d = _op_cast_int_float128_t(m);
	for ( ; m > 1; m -= 2) {
		t1 = _op_plus_float128_t_float128_t(d,f);
		f = _op_divide_float128_t_float128_t(s,t1);
		d = _op_cast_int_float128_t(m-2);
	}
	t1 = _op_plus_float128_t_float128_t(d,f);
	f = _op_divide_float128_t_float128_t(x,t1);
	t1 = _op_plus_float128_t_float128_t(d,f);
	t2 = _op_minus_float128_t_float128_t(d,f);
	f = _op_divide_float128_t_float128_t(t1,t2);
	exp = extractFloat128Exp(f);
	exp += k;
	w0 = extractFloat128Frac0(f);
	w1 = extractFloat128Frac1(f);
	return packFloat128(sign,exp,w0,w1);
}

float128_t expf128(float128_t x)
{
	float128_t tmp = _op_multiply_float128_t_float128_t(x,__xLog2e);
	tmp = exp2f128(tmp);
	return tmp;
}

float128_t exp10f128(float128_t x)
{
	float128_t tmp = _op_multiply_float128_t_float128_t(x,__xLog2_10);
	return exp2f128(tmp);
}
#endif

float128_t copysignf128(float128_t x,float128_t y)
{
	float128_t result;

	result.low = x.low;
	result.high = x.high;
	if (y.high>>63) {
		result.high |= 0x8000000000000000ULL;
	}
	else result.high &= ~0x8000000000000000ULL;
	return result;
}

int isinff128(float128_t x)
{
	int exp = extractFloat128Exp(x);
	bits64 sig0,sig1;

	sig0 = extractFloat128Frac0(x);
	sig1 = extractFloat128Frac1(x);
	if ((exp == 0x7fff) && (sig0|sig1)== 0) return 1;
	return 0;
}

//int __declspec(naked) isfinitef128(float128_t x) {}
int finitef128(float128_t x)
{
	int exp = extractFloat128Exp(x);
	/* Bug in the original code fixed. */
	/*if (exp <= 0x7ffe && exp >= 0x10000)  
	This is wrong! exp can't be <= 32766 and bigger than 65536! 
	The comparison is done as UNSIGNED
	*/
	if (exp <= 32766 && exp >= -32768) // This is wrong! exp can't be <= 32766 and bigger than 65536!
		return 1;
	return 0;
}
#if 0
// This needs reworking
int fpclassifyf128(float128_t x)
{
	int exp = extractFloat128Exp(x);
	if (exp <= 0x7ffe && exp >= 0x10000) {
		if (x.high|x.low)
			return 0; // FP_NORMAL
		return 64; // FP_ZERO
	}
	if (isnanf128(x)) return 1; // FP_NAN
	if (isinff128(x)) return 5; // FP_INFINITE
	return 69; // FP_SUBNORMAL
}
#endif

static const float128_t
// 2.0769187434139310514121985316880384E+34F128
two114 = { 0, 0x4071000000000000ULL },
// 4.8148248609680896326399448564623183E-35F128 
twom114 = {0, 0x3F8D000000000000ULL },
// huge   = 1.0E+4900F128
huge = {0x498151922505CB5FULL, 0x7F945D24084EB26FULL},
// tiny   = 1.0E-4900F128;
tiny = {0x52E4D25544B1042EULL, 0x00697769BEAD75ECULL}; 

float128_t scalbnf128 (float128_t x, int n)
{
	bits64 k,hx,lx;

	hx = x.high,lx = x.low;
        k = extractFloat128Exp(x);
        if (k==0) {				/* 0 or subnormal x */
            if ((lx|(hx&0x7fffffffffffffffULL))==0) return x; /* +-0 */
	    x = _op_multiply_float128_t_float128_t(x,two114);
	    hx = x.high;
	    k = ((hx>>48)&0x7fff) - 114;
	}
        if (k==0x7fff) return _op_plus_float128_t_float128_t(x,x);		/* NaN or Inf */
        k = k+n;
        if (n> 50000 || k > 0x7ffe)
	  return _op_multiply_float128_t_float128_t(huge,copysignf128(huge,x)); /* overflow  */
	if (n< -50000) return  /* underflow! */
                   _op_multiply_float128_t_float128_t(tiny,copysignf128(tiny,x));
        if (k > 0) 				/* normal result */
	    {x.high = (hx&0x8000ffffffffffffULL)|(k<<48); return x;}
        if (k <= -114) /* underflow */
	  return _op_multiply_float128_t_float128_t(tiny,copysignf128(tiny,x)); 
        k += 114;				/* subnormal result */
	x.high = (hx&0x8000ffffffffffffULL)|(k<<48);
        return _op_multiply_float128_t_float128_t(x,twom114);
}
#if 0
/*
https://gist.github.com/publik-void/067f7f2fef32dbe5c27d6e215f824c91
*/
float128_t sinf128(float128_t x)
{
	float128_t x2 = _op_multiply_float128_t_float128_t(x,x);
	return x*(1.0L + 
	x2*(-0.166666666666666666666666666666666667L + 
	x2*(0.00833333333333333333333333333333333069L + 
	x2*(-0.000198412698412698412698412698412671319L +
	x2*(2.75573192239858906525573192223995808e-6L + 
	x2*(-2.50521083854417187750521077962123682e-8L + 
	x2*(1.60590438368216145993922289621550506e-10L + 
	x2*(-7.64716373181981647587481187300831335e-13L + 
	x2*(2.81145725434552075980975905006999319e-15L + 
	x2*(-8.22063524662432650297086257962703293e-18L + 
	x2*(1.95729410633890026175367390152305383e-20L + 
	x2*(-3.86817017051340241224838720319634797e-23L + 
	x2*(6.44695023999222092772271073593727141e-26L + 
	x2*(-9.1836779606017064087088551595474321e-29L + 
	x2*(1.13078207057779775850779192271873238e-31L - 
	1.19290046424220296937971101373203567e-34L * x2)))))))))))))));
}
#endif
#ifdef __LCC__
typedef __int128 int128;
float128_t _op_cast___int128_float128_t(int128 a)
{
	float128_t r,b;
	r = _op_cast_long_long_float128_t(a>>64);
	b = _op_cast_long_long_float128_t((long long)a);
	r = ldexpf128(r,64);
	r = _op_plus_float128_t_float128_t(r,b);
	return r;
}
float128_t *_op_asgn_pfloat128_t___int128(float128_t *pfloat128_t,int128 x)
{
	float128_t r,b;
	b = _op_cast___int128_float128_t(x);
	*pfloat128_t = b;
	return pfloat128_t;
}
#endif
