	.text
# --------------------------------------------------------------shiftupn-------
	.globl	shiftupn
shiftupn:
# This function should never use x14,x11 since it is called from mdnorm
#            n --> x1
#            x --> x0
# (WORDSIZE-n) --> x5
#
#void shiftupn(QfloatAccump x,int n)
#{
#	int i;
#	QELT newbyt, oldbyt;
#
#	oldbyt = 0;
#	for( i=ACCUM_LENGTH; i>= 0; i-- ) {
#		newbyt = x->mantissa[i] >> (WORDSIZE - n);
#		x->mantissa[i] <<= n;
#		x->mantissa[i] |= oldbyt;
#		oldbyt = newbyt;
#	}
#}
	ldp x9,x8,[x0,64]
	mov x6,64
	sub x5,x6,x1
	lsr x4,x8,x5

	lsl x8,x8,x1
	lsr x7,x9,x5
	lsl x9,x9,x1
	orr x9,x9,x4
	stp x9,x8,[x0,64]
# Alternate between the two registers loaded, making the same operations
# in both to maximize parallelism
	ldp x3,x2,[x0,48]
	lsr x4,x2,x5
	lsr x6,x3,x5
	lsl x2,x2,x1
	lsl x3,x3,x1
	orr x2,x2,x7
	orr x3,x3,x4
	stp x3,x2,[x0,48]

	ldp x9,x8,[x0,32]
	lsr x4,x8,x5
	lsr x7,x9,x5
	lsl x8,x8,x1
	lsl x9,x9,x1
	orr x8,x8,x6
	orr x9,x9,x4
	stp x9,x8,[x0,32]

	ldp x3,x2,[x0,16]
	lsr x4,x2,x5
	lsr x6,x3,x5
	lsl x2,x2,x1
	lsl x3,x3,x1
	orr x2,x2,x7
	orr x3,x3,x4
	stp x3,x2,[x0,16]

	ldr x8,[x0,8]
	lsl x8,x8,x1
	orr x8,x8,x6
	str x8,[x0,8]
	ret
	.balign 8
# -----------------------------------------------------------------------------
#  void shiftdownn(QfloatAccump x,int n)
#  {
#  	QELT newbyt, oldbyt;
#  	int i;
#  
#  	oldbyt = 0;
#  	for( i=0; i<=ACCUM_LENGTH; i++ ) {
#  		newbyt = x->mantissa[i] << (WORDSIZE-n);
#  		x->mantissa[i] >>= n;
#  		x->mantissa[i] |= oldbyt;
#  		oldbyt = newbyt;
#  	}
#  }
shiftdownn:
	mov x4,64
	sub x5,x4,x1

	ldp x2,x3,[x0,8]
	lsl x4,x2,x5
	lsr x2,x2,x1
	lsl x7,x3,x5
	lsr x3,x3,x1
	orr x3,x3,x4
	stp x2,x3,[x0,8]

	ldp x2,x3,[x0,24]
	lsl x4,x2,x5
	lsr x2,x2,x1
	lsl x7,x3,x5
	orr x2,x2,x7

	lsr x3,x3,x1
	orr x3,x3,x4
	stp x2,x3,[x0,24]

	ldp x2,x3,[x0,40]
	lsl x4,x2,x5
	lsr x2,x2,x1
	orr x2,x2,x7
	lsl x7,x3,x5
	lsr x3,x3,x1
	orr x3,x3,x4
	stp x2,x3,[x0,40]

	ldp x2,x3,[x0,56]
	lsl x4,x2,x5
	lsr x2,x2,x1
	orr x2,x2,x7
	lsl x7,x3,x5
	lsr x3,x3,x1
	orr x3,x3,x4
	stp x2,x3,[x0,56]

	ldr x2,[x0,72]
	lsr x2,x2,x1
	orr x2,x2,x7
	str x2,[x0,72]
	ret
	.balign 8
	.globl	shiftdownn
# -----------------------------------------------------------------------------
	.globl	shdn1
# void shdn1(QfloatAccump x)
# This function should never use x14,x11, since it is called from mdnorm
shdn1:
    ldp x2,x3,[x0,8]
    ldp x4,x5,[x0,24]
    extr x6,xzr,x2,1
    extr x7,x2,x3,1
    stp x6,x7,[x0,8]

    ldp x8,x1,[x0,40]
    extr x6,x3,x4,1
    ldp x2,x3,[x0,56]
    extr x7,x4,x5,1
    stp x6,x7,[x0,24]

    extr x6,x5,x8,1
    extr x7,x8,x1,1
    stp x6,x7,[x0,40]

    ldr x4,[x0,72]
    extr x6,x1,x2,1
    extr x7,x2,x3,1
    stp  x6,x7,[x0,56]

    extr x2,x3,x4,1
    str x2,[x0,72]
	ret
	.balign 8
# -----------------------------------------------------------------------------
	.globl shup1
shup1:
#   9,8
	ldr x3,[x0,72]
	ldp x4,x1,[x0,56]
	adds x12,x3,x3
	str x12,[x0,72]
#   7,6
	adcs x6,x1,x1
	ldp x3,x10,[x0,40]
	adcs x7,x4,x4
	stp x7,x6,[x0,56]
#   5,4
	adcs x13,x10,x10
	ldp x4,x1,[x0,24]
	adcs x12,x3,x3
	stp x12,x13,[x0,40]
#   3,2
	adcs x9,x1,x1
	ldp x3,x10,[x0,8]
	adcs x12,x4,x4
	adcs x13,x10,x10
	stp x12,x9,[x0,24]
#   1,0
	adc x12,x3,x3
#   0
	stp x12,x13,[x0,8]
#   }
	ret
# ----------------------------------------------------------------------mulv---
	.balign 8
# Variable precision multiply of significands.
# static void mulv(QfloatAccump a,QfloatAccump b,QfloatAccump c,int prec)
# x0 first input
# x1 second input
# x2 output
# x3 --> x1,16
# x4 --> x1,24
# x9 *r
# x10 --> x0,16
# x11 --> x0,24
# x12 *(r-1)
# x13 *(r-2)
# x14 --> x1,32
# x0 --> x0,32
# x1 --> x0,40
# x2  --> x1,40
# x16 --> x1,48
# Data: Arg1 in q16,q20
#       Arg2 in v26,v30
# Result in v3,v7
# Result is also present in x0 --> v3.d[0], x1 --> v3.d[1] and x2 --> v4.d[0]
mulv:
# 						8
	mov x6,v19.d[1]
	mov x3,v26.d[0]
	mul x9,x6,x3
	umulh x12,x6,x3
mulvEntry1:
#	mov x8,v20.d[0]
	
# 						7	x0[72] * x1[16] ( x0->mantissa[8] * x1->mantissa[1] )
	mov x10,v16.d[0]
	mov x11,v16.d[1]
	mov x0,v17.d[0]
	mov x1,v17.d[1]
	mov x6,v19.d[0]
	mov x8,v18.d[1]
#						6	x0[64] * x1[24] ( x0->mantissa[7] * x1->mantissa[2] )
mulvEntry2:
	mov x4,v26.d[1]
	mov x14,v27.d[0]
	mul x5,x6,x4
	umulh x7,x6,x4
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,xzr,xzr
#						5	x0[56] * x1[32] ( x0->mantissa[6] * x1->mantissa[3] )
	mul x5,x8,x14
	umulh x7,x8,x14
	adds x9,x9,x5
	mov x2,v27.d[1]
	mov x16,v28.d[0]
	adcs x12,x12,x7
	mov x8,v18.d[0]
	adc x13,x13,xzr
#						4	x0[48] * x1[40] ( x0->mantissa[5] * x1->mantissa[4] )
	mul x5,x8,x2
	umulh x7,x8,x2
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						3	x0[40] * x1[48] ( x0->mantissa[4] * x1->mantissa[5] )
	mul x5,x1,x16
	umulh x7,x1,x16
	adds x9,x9,x5
	mov x15,v28.d[1]
	adcs x12,x12,x7
	adc x13,x13,xzr
#						2	x0[32] * x1[56] ( x0->mantissa[3] * x1->mantissa[6] )
	mul x5,x0,x15
	umulh x7,x0,x15
	adds x9,x9,x5
	mov x15,v29.d[0]
	adcs x12,x12,x7
	adc x13,x13,xzr
#						1	x0[24] * x1[64] ( x0->mantissa[2] * x1->mantissa[7] )
	mul x5,x11,x15
	umulh x7,x11,x15
	adds x9,x9,x5
	mov x15,v29.d[1]
	adcs x12,x12,x7
	adc x13,x13,xzr
#						0	x0[16] * x1[72] ( x0->mantissa[1] * x1->mantissa[8] )
#	mul x5,x10,x15
	umulh x7,x10,x15
#	adds x9,x9,x5
	adds x12,x12,x7
	adc x13,x13,xzr

#  Write out *r, we are done with that memory location
#	mov v7.d[0],x9
#						6	x0[64] * x1[16] ( x0->mantissa[7] * x1->mantissa[1] )
#	                         x6    *   x3
	mul x5,x6,x3
	umulh x7,x6,x3
	adds x9,x12,x5
	adcs x12,x13,x7
	mov x6,v18.d[1]
	adc x13,xzr,xzr
#						5	x0[56] * x1[24] ( x0->mantissa[6] * x1->mantissa[2] )
	mul x5,x6,x4
	umulh x7,x6,x4
	adds x9,x9,x5
	adcs x12,x12,x7
#	x8 assigned above
	adc x13,x13,xzr
#						4	x0[48] * x1[32] ( x0->mantissa[5] * x1->mantissa[3] )
	mul x5,x8,x14
	umulh x7,x8,x14
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						3	x0[40] * x1[40] ( x0->mantissa[4] * x1->mantissa[4] )
	mul x5,x1,x2
	umulh x7,x1,x2
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						2	x0[32] * x1[48] ( x0->mantissa[3] * x1->mantissa[5] )
	mul x5,x0,x16
	umulh x7,x0,x16
	adds x9,x9,x5
	mov x15,v28.d[1]
	adcs x12,x12,x7
	adc x13,x13,xzr
#						1	x0[24] * x1[56] ( x0->mantissa[2] * x1->mantissa[6] )
	mul x5,x11,x15
	umulh x7,x11,x15
	adds x9,x9,x5
	mov x15,v29.d[0]
	adcs x12,x12,x7
	adc x13,x13,xzr
#						0	x0[16] * x1[64] ( x0->mantissa[1] * x1->mantissa[7] )
	mul x5,x10,x15
	umulh x7,x10,x15
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr

#  Write out *r, we are done with that memory location
	mov v6.d[1],x9
#						5	x0[56] * x1[16] ( x0->mantissa[6] * x1->mantissa[1] )
	mul x5,x6,x3
	umulh x7,x6,x3
	adds x9,x12,x5
	adcs x12,x13,x7
	mov x6,v18.d[0]
	adc x13,xzr,xzr
#						4	x0[48] * x1[24] ( x0->mantissa[5] * x1->mantissa[2] )
	mul x5,x6,x4
	umulh x7,x6,x4
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						3	x0[40] * x1[32] ( x0->mantissa[4] * x1->mantissa[3] )
#							 x1   *   x14
	mul x5,x1,x14
	umulh x7,x1,x14
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						2	x0[32] * x1[40] ( x0->mantissa[3] * x1->mantissa[4] )
	mul x5,x0,x2
	umulh x7,x0,x2
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						1	x0[24] * x1[48] ( x0->mantissa[2] * x1->mantissa[5] )
	mul x5,x11,x16
	umulh x7,x11,x16
	adds x9,x9,x5
	mov x15,v28.d[1]
	adcs x12,x12,x7
	adc x13,x13,xzr
#						0	x0[16] * x1[56] ( x0->mantissa[1] * x1->mantissa[6] )
	mul x5,x10,x15
	umulh x7,x10,x15
	adds x17,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr

#  Write out *r, we are done with that memory location
	mov v6.d[0],x17
#						4	x0[48] * x1[16] ( x0->mantissa[5] * x1->mantissa[1] )
#							  x6   *   x3
	mul x5,x6,x3
	umulh x7,x6,x3
	adds x9,x12,x5
	adcs x12,x13,x7
	adc x13,xzr,xzr
#						3	x0[40] * x1[24] ( x0->mantissa[4] * x1->mantissa[2] )
	mul x5,x1,x4
	umulh x7,x1,x4
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						2	x0[32] * x1[32] ( x0->mantissa[3] * x1->mantissa[3] )
	mul x5,x0,x14
	umulh x7,x0,x14
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						1	x0[24] * x1[40] ( x0->mantissa[2] * x1->mantissa[4] )
#							  x11  *   x6
	mul x5,x11,x2
	umulh x7,x11,x2
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						0	x0[16] * x1[48] ( x0->mantissa[1] * x1->mantissa[5] )
	mul x5,x10,x16
	umulh x7,x10,x16
	adds x16,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr

#  Write out *r, we are done with that memory location
	mov v5.d[1],x16
#						3	x0[40] * x1[16] ( x0->mantissa[4] * x1->mantissa[1] )
	mul x5,x1,x3
	umulh x7,x1,x3
	adds x9,x12,x5
	adcs x12,x13,x7
	adc x13,xzr,xzr
#						2	x0[32] * x1[24] ( x0->mantissa[3] * x1->mantissa[2] )
	mul x5,x0,x4
	umulh x7,x0,x4
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						1	x0[24] * x1[32] ( x0->mantissa[2] * x1->mantissa[3] )
	mul x5,x11,x14
	umulh x7,x11,x14
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						0	x0[16] * x1[40] ( x0->mantissa[1] * x1->mantissa[4] )
#							  x10  *  x6
	mul x5,x10,x2
	umulh x7,x10,x2
	adds x15,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr

#  Write out *r, we are done with that memory location
	mov v5.d[0],x15
#						2	x0[32] * x1[16] ( x0->mantissa[3] * x1->mantissa[1] )
	mul x5,x0,x3
	umulh x7,x0,x3
	adds x9,x12,x5
	adcs x12,x13,x7
	adc x13,xzr,xzr
#						1	x0[24] * x1[24] ( x0->mantissa[2] * x1->mantissa[2] )
	mul x5,x11,x4
	umulh x7,x11,x4
	adds x9,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr
#						0	x0[16] * x1[32] ( x0->mantissa[1] * x1->mantissa[3] )
	mul x5,x10,x14
	umulh x7,x10,x14
	adds x14,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr

#  Write out *r, we are done with that memory location
	mov v4.d[1],x14
#						1	x0[24] * x1[16] ( x0->mantissa[2] * x1->mantissa[1] )
	mul x5,x11,x3
	umulh x7,x11,x3
	adds x9,x12,x5
	adcs x12,x13,x7
	adc x13,xzr,xzr
#						0	x0[16] * x1[24] ( x0->mantissa[1] * x1->mantissa[2] )
	mul x5,x10,x4
	umulh x7,x10,x4
	adds x2,x9,x5
	adcs x12,x12,x7
	adc x13,x13,xzr

#  Write out *r, we are done with that memory location
	mov v4.d[0],x2
#						0	x0[16] * x1[16] ( x0->mantissa[1] * x1->mantissa[1] )
	mul x5,x10,x3
	umulh x7,x10,x3
	adds x1,x12,x5
	adcs x0,x13,x7

#  Write out *r, we are done with that memory location
	mov v3.d[1],x1
# Now write out the last two registers
	mov v3.d[0],x0
	ret
	.balign 8
# ----------------------------------------------------------------------Step2--
# At this point sqr(x1) has only two positions filled: mantissa[1] and mantissa[2].
Step2:
# Input data: qa in v16-v20 
#             sqr in registers x14 and x11
# Output:     prod, stored in v3,v7
#  160 	r = &c->mantissa[4];
	mov x3,v16.d[0]
#                                       4
#  This multiplies only the high part, that generates carry and initial values
#  for the operations on lower values 
	mov x9,v17.d[0]
	umulh x12,x9,x11
	mov x8,v16.d[1]
#  165 			if( (*p) && (*q)) {
#	mul x5,x8,x14
	umulh x7,x8,x14
	adds x12,x12,x7
# this is loading zero
#	ldr x15,[x1,32]
	adc x13,xzr,xzr
#  The multiplication of x1, 32 gives always zero since x1,32 is zero
#	mul x5,x3,x15
#	umulh x7,x3,x15
#	adds x4,x4,x5
#	adcs x12,x12,x7
#	str x4,[x2,40]
#	mov x13,xzr
#  Shift temps. First write out *r, we are done with that memory location
#                                       3
#  165 			if( (*p) && (*q)) {
	mul x5,x8,x11
	umulh x7,x8,x11
	adds x4,x12,x5
	adcs x15,x13,x7
	mul x5,x3,x14
	adc x13,xzr,xzr
#  165 			if( (*p) && (*q)) 
	umulh x7,x3,x14
	adds x10,x4,x5
	mul x5,x3,x11
	adcs x12,x15,x7
	umulh x7,x3,x11
	adc x13,x13,xzr
#  Shift temps. First write out *r, we are done with that memory location
#                                       2
#  165 			if( (*p) && (*q)) {
	adds x4,x12,x5
	adcs x15,x13,x7

	subs x23,xzr,x10

	sbcs x5, xzr,x4

	sbc x6,x21,x15
# shift up 1 bit.
	adds x4,x4,x4
	adcs x22,x5,x5
	adc  x21,x6,x6
# store result
#	mov x24,xzr
# -----------------------------------------------------------------------------
#   Step 4: mulv4(qa,sqr,prod)
# inputs: qa in v16-v20
#         sqr in v26-v30
# output: prod in v3,v7
#Step4:
#First: square, i.e. square(quot,sqr,4), sqr = quot*quot
#      quot is in v21

#	mov x6,x21
#	mov x5,x22
	mul x14,x22,x22
	umulh x7,x22,x22
	mul x8,x21,x21
	umulh x9,x21,x21

	subs x3,x22,x21
	sbcs x4,x21,x21

	eor x3,x3,x4
	and x4,x4,1
	add x3,x3,x4

	mul x10,x3,x3
	umulh x11,x3,x3

	adds x12,x14,x8
	adcs x13,x7,x9
	adc x9,x9,xzr

	adds x7,x7,x12
	adcs x8,x8,x13
	adc x9,x9,xzr

	subs x7,x7,x10
	sbcs x8,x8,x11
	sbc x9,x9,xzr

#inline shup1
	adds x14,x14,x14
	adcs x7,x7,x7
	adcs x13,x8,x8
	adc x3,x9,x9

#Second: multiply mulv4(qa,sqr,prod); prod = qa * sqr
# Input sqr in x3 and x8
# input qa in v16

	mov x5,v16.d[0]
	mov x4,v16.d[1]

	MUL x6, x13, x4
	UMULH x7, x13, x4 
	MUL x8, x3, x5
	UMULH x9, x3, x5
	ADDS x10, x6, x8
	ADCS x11, x7, x9 
	ADC x12, xzr, xzr 
	ADDS x7, x7, x10 
	ADCS x8, x8, x11 
	ADC x9, x9, x12
	SUBS x13, x13, x3 
	SBC x3, x3, x3 
	EOR x13, x13, x3 
	AND x3, x3, 1 
	ADD x13, x13, x3 
	SUBS x4, x4, x5 
	SBC x5, x5, x5 
	EOR x4, x4, x5 
	AND x5, x5, 1 
	ADD x4, x4, x5 
	EOR x3, x3, x5
	SUB x3, x3, 1
	MUL x10, x13, x4 
	UMULH x11, x13, x4
	EOR x10, x10, x3 
	EOR x11, x11, x3 
	AND x4, x3, #1
	ADDS x10, x10, x4 
	ADCS x11, x11, xzr
	ADC x3, x3, xzr
	ADDS x7, x7, x10 
	ADCS x8, x8, x11 
	ADC x9, x9, x3
# Third: subtract quot = quot -prod. quot is in v21-v22
	subs x2,xzr,x6

	sbcs x0,x23,x7
	sbcs x4,x22,x8

	sbc x14,x21,x9
#   5,4
	adds x24,x2,x2
#   3,2
	adcs x23,x0,x0
	adcs x22,x4,x4
	adcs x21,x14,x14
#   1,0
#   0

#	ret
#	mov x2,xzr
#	mov x5,xzr
#	mov x12,xzr
	mul x14,x24,x24
	umulh x2,x24,x24

	mul x6,x24,x23
	umulh x7,x24,x23
	add x5,xzr,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,xzr
	adc x5,x5,xzr
	mov v29.d[0],x2
	mov x2,x12
	mov x12,x5
	mov x5,xzr
	b	squareAsmEntry1
	.balign 8
# -------------------------------------------------------------------square----
#  193 static void square(QfloatAccump a,QfloatAccump b,int prec )
#	.globl	squareAsm
	.globl square
square:
	mov v21.d[0],x21
	mov v21.d[1],x22
	mov v22.d[0],x23
	mov v22.d[1],x24
	mov v23.d[0],x25
	mov v23.d[1],x26
	ldp x21,x22,[x0,16]
	mov v24.d[0],x27
	mov v24.d[1],x28
	mov v0.d[0],x30
	mov v0.d[1],x1
	ldp x23,x24,[x0,32]
	ldp x25,x26,[x0,48]
	ldp x27,x28,[x0,64]
	bl squareAsm
	mov x0,v0.d[1]
	mov x21,v21.d[0]
	mov x22,v21.d[1]
	mov x23,v22.d[0]
	mov x24,v22.d[1]
	mov x25,v23.d[0]
	mov x26,v23.d[1]
	mov x27,v24.d[0]
	mov x28,v24.d[1]
	str x6,[x0,8]
	stp q26,q27,[x0,16]
	stp q28,q29,[x0,48]
	mov x30,v0.d[0]
	ret
# ---------------------------------------------------------------squareAsm-----
#  188 /* Variable precision square. */
#  193 static void square(QfloatAccump a,QfloatAccump b,int prec )
	.balign 8
#  188 /* Variable precision square. square(quot,sqr,8) */
squareAsm:
#  Registers at entry:
# input: quot (x21-x28)
# output: sqr (v26-v30)
# carry returned in x6
#								9
#								16,80
#	ldr x15,[x0,72]
#	ldp x10,x11,[x0,48]
#	mov x10,x25
#	mov x11,x26
#								8
#								16,72
	mul x2,x28,x21
	umulh x12,x28,x21
	add x5,xzr,x12,lsr 63
	adds x2,x2,x2
	adc x12,x12,x12
#								24,64
	mul x6,x27,x22
	umulh x7,x27,x22
	add x5,x5,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	adc x5,x5,xzr
#								32,56
	mul x6,x26,x23
	umulh x7,x26,x23
	add x5,x5,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	adc x5,x5,xzr
#								40,48
	mul x6,x25,x24
	umulh x7,x25,x24
	add x5,x5,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x4,x6,x2
	adcs x2,x7,x12
	adc x12,x5,xzr
#	mov x2,x12
#	mov x12,x5
#								7
#								16,64
	mul x6,x27,x21
	umulh x7,x27,x21
	add x5,xzr,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	adc x5,x5,xzr
#								24,56
	mul x6,x26,x22
	umulh x7,x26,x22
	add x5,x5,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	adc x5,x5,xzr
#								32,48
	mul x6,x25,x23
	umulh x7,x25,x23
	add x5,x5,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	adc x5,x5,xzr
#								40,40
	mul x6,x24,x24
	umulh x7,x24,x24
	adds x14,x6,x2
	adcs x2,x7,x12
	adc x12,x5,xzr
#								6
#								16,56
	mul x6,x26,x21
	umulh x7,x26,x21
	add x5,xzr,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	adc x5,x5,xzr
#								24,48
	mul x6,x25,x22
	umulh x7,x25,x22
	add x5,x5,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	adc x5,x5,xzr
#								32,40
	mul x6,x24,x23
	umulh x7,x24,x23
	add x5,x5,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	mov v29.d[0],x2
	adc x5,x5,xzr
	mov x2,x12
	mov x12,x5
#								5
#								16,48
	mul x6,x25,x21
	umulh x7,x25,x21
	add x5,xzr,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	adc x5,x5,xzr
#								24,40
squareAsmEntry1:
	mul x6,x24,x22
	umulh x7,x24,x22
	add x5,x5,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	adc x5,x5,xzr
#								32,32
	mul x6,x23,x23
	umulh x7,x23,x23
# Save in x9 v28.d[1]
	adds x9,x6,x2
#	adcs x12,x7,x12
	adcs x2,x7,x12
#	adc x5,x5,xzr
	adc x12,x5,xzr
#								4
#								16,40
	mul x6,x24,x21
	umulh x7,x24,x21
	add x5,xzr,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
	adc x5,x5,xzr
#								24,32
	mul x6,x23,x22
	umulh x7,x23,x22
	add x5,x5,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
# Save in x10 v28.d[0]
	adds x10,x6,x2
	adcs x2,x7,x12
	adc x12,x5,xzr
#								3
#								16,32
	mul x6,x23,x21
	umulh x7,x23,x21
	add x5,xzr,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
	adds x2,x6,x2
	adcs x12,x7,x12
#								24,24
	mul x6,x22,x22
	umulh x7,x22,x22
# Save in x8 v27.d[1]
	adds x8,x6,x2
	adcs x2,x7,x12
	adc x12,x5,xzr
#								2
#								16,24
	mul x6,x22,x21
	umulh x7,x22,x21
	add x5,xzr,x7,lsr 63
	adds x6,x6,x6
	adc x7,x7,x7
#save in x1 v27.d[0]
	adds x1,x6,x2
	adcs x2,x7,x12
	adc x12,x5,xzr
#								1
#								16,16
	mul x6,x21,x21
	umulh x7,x21,x21
# Save in x0 the value of v26.d[1] instead of writing it down
	adds x0,x6,x2
# Shift the sliding window by storing in the shifted registers
# instead of mov ing them later
	adcs x2,x7,x12
	adc x12,xzr,xzr

#  The last iteration of the loop does NOT write the contents of x12 and x2 
#  ( *(r-1) and *(r-2) ) to memory. Do it now since the loop is finished
#  we leave the data in x2 and x12 instead of writing to memory. Later below
#  that will be used and THEN written
#inline of shup1
#	mov x6,x4
	adds x7,x4,x4
# This value is never used
#	mov v30.d[0],x7
	mov x5,v29.d[0]
# Use the value in x14 saved above
	adcs x7,x14,x14
	mov v29.d[1],x7
#   7,6
	adcs x13,x5,x5
	mov x14,v27.d[1]
#	mov x10,v28.d[0]
# Use the value in x10 saved in 24,32 above, and in x9 in 32,32
	adcs x7,x9,x9
	mov v28.d[1],x7
	mov v29.d[0],x13
#   5,4
	adcs x13,x10,x10
#	mov x6,v26.d[1]
#	mov x6,x0
#	mov x5,v27.d[0]
# Use the value saved in x8 at 24,24 above
	adcs x7,x8,x8
	mov v27.d[1],x7
	mov v28.d[0],x13
#   3,2
# Use the value saved in x1 at 16,24 above
	adcs x5,x1,x1
	adcs x7,x0,x0
	mov v26.d[1],x7
	mov v27.d[0],x5
#   1,0
#	Add now x2 and x12. This spares 16 bytes (write)+16 bytes (read) above
	adcs x5,x2,x2
	adc x6,x12,x12
#   0
	mov v26.d[0],x5
#  231 }
	ret
# -----------------------------------------------------------end-of-squareAsm--
# --------------------------------------------------------------------addm-----
	.balign	8
	.globl	addm
#   
#  void addm(QfloatAccump x,QfloatAccump y )
addm:
#
	ldr x15,[x0,72]
	ldr x12,[x1,72]
	adds x13,x15,x12
	str x13,[x1,72]

	ldp x15,x10,[x0,56]
	ldp x12,x11,[x1,56]
	adcs x14,x10,x11
	ldp x5,x3,[x0,40]
	adcs x13,x15,x12
	ldp x2,x4,[x1,40]
	stp x13,x14,[x1,56]

	adcs x14,x3,x4
	ldp x7,x9,[x0,24]
	adcs x13,x5,x2
	ldp x8,x6,[x1,24]
	stp x13,x14,[x1,40]

	adcs x11,x9,x6
	ldp x15,x10,[x0,8]
	adcs x13,x7,x8
	ldp x12,x6,[x1,8]
	stp x13,x11,[x1,24]

	adcs x14,x10,x6
	adc x13,x12,x15
	stp x13,x14,[x1,8]
	ret
	.balign	8
# -----------------------------------------------------------------------subm--
#  void subm(QfloatAccump x,QfloatAccump y )
subm:
#	ldr x15,[x0,80]
#	ldr x12,[x1,80]
#	subs x13,x12,x15
#	str x13,[x1,80]
	ldr x15,[x0,72]
	ldr x12,[x1,72]
	subs x13,x12,x15
	str x13,[x1,72]

	ldp x15,x10,[x0,56]
	ldp x12,x11,[x1,56]
	sbcs x14,x11,x10
	ldp x5,x8,[x0,40]
	sbcs x13,x12,x15
	ldp x2,x9,[x1,40]
	stp x13,x14,[x1,56]

	sbcs x7,x9,x8
	ldp x15,x10,[x0,24]
	sbcs x13,x2,x5
	ldp x12,x11,[x1,24]
	stp x13,x7,[x1,40]

	sbcs x14,x11,x10
	sbcs x13,x12,x15
	stp x13,x14,[x1,24]

	ldp x15,x10,[x0,8]
	ldp x12,x11,[x1,8]
	sbcs x14,x11,x10
	sbcs x13,x12,x15
	stp x13,x14,[x1,8]
	ret
	.globl	subm
# -------------------------------------------------------------submProdQuot----
	.balign 8
# Input: x0  --> prod->mantissa[1] 
#        x1 -->  prod->mantissa[2] 
#        x2 -->  prod->mantissa[3]
#        x14 --> prod->mantissa[4]
#        x15 --> prod->mantissa[5]
#        x16 --> prod->mantissa[6]
#		 x17 --> prod->mantissa[7]
#    v6.d[1] --> prod->mantissa[8]
#       x21  --> quot->mantissa[1]
#       x22 -->  quot->mantissa[2]
#             ...
#       x28 -->  quot->mantissa[8]
# Output: Modified registers x21-x28
#         output = 2 * (quot - prod)
submProdQuot:
# 1: subtract
	mov x11,v6.d[1]
	subs x28,x28,x11
	sbcs x27,x27,x17
	sbcs x26,x26,x16
	sbcs x25,x25,x15
	sbcs x24,x24,x14
	sbcs x23,x23,x2
	sbcs x22,x22,x1
	sbc x21,x21,x0
# 2: multiply by 2.
	adds x28,x28,x28
	adcs x27,x27,x27
	adcs x26,x26,x26
	adcs x25,x25,x25
	adcs x24,x24,x24
	adcs x23,x23,x23
	adcs x22,x22,x22
	adc x21,x21,x21
# Done
	ret
	.balign 8
# -----------------------------------------------------------------------mulm--
#  291 void mulm(Qfloatp b,Qfloatp b, QfloatAccump ac3)
mulm:
# map
# Input: x0,x1 point to inputs (Qfloat)
# Output x2 points to QfloatAccum
# x1
# x1  x1->mantissa[1]  8
# x5  x1->mantissa[2] 16
# x9  x1->mantissa[3] 24
# x11 x1->mantissa[4] 32
# x20 x1->mantissa[5] 40
# x22 x1->mantissa[6] 48
# x15 x1->mantissa[7] 56
# x0
# x4  x0->mantissa[1]  8
# x8  x0->mantissa[2] 16
# x10 x0->mantissa[3] 24
# x17 x0->mantissa[4] 32
# x21 x0->mantissa[5] 40
# x16 x0->mantissa[6] 48
# x0  x0->mantissa[7] 56
# Sliding window
# x12 sliding window  0
# x13 sliding window -1
# x14 sliding window -2
	mov v0.d[0],x20
	mov v0.d[1],x21
	mov v1.d[0],x22
	ldr x6,[x0]
	str x6,[x2]
	ldp x9,x11,[x1,24]
	ldp x20,x22,[x1,40]
	ldr x15,[x1,56]
	ldp x1,x5,[x1,8]
	ldp x4,x8,[x0,8]
	ldp x10,x17,[x0,24]
	ldp x21,x16,[x0,40]
	ldr x0,[x0,56]

#	7
	mul x12,x0,x5
	umulh x13,x0,x5
#	6
	mul x3,x16,x9
	umulh x7,x16,x9
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,xzr,xzr
#	5
	mul x3,x21,x11
	umulh x7,x21,x11
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr
#	4
	mul x3,x17,x20
	umulh x7,x17,x20
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr
#	3
	mul x3,x10,x22
	umulh x7,x10,x22
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr
#	2
	mul x3,x8,x15
	umulh x7,x8,x15
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr
#	1
# Store result
#	str x12,[x2,88]
# ------------------------------------------------------9
#  ------------------------------------8
#  ------------------------------------7
	mul x3,x0,x5
	umulh x7,x0,x5
	adds x12,x3,x13
	adcs x13,x7,x14
	adc x14,xzr,xzr

#  ------------------------------------6
	mul x3,x16,x9
	umulh x7,x16,x9
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------5
	mul x3,x21,x11
	umulh x7,x21,x11
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------4
	mul x3,x17,x20
	umulh x7,x17,x20
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------3
	mul x3,x10,x22
	umulh x7,x10,x22
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------2
	mul x3,x8,x15
	umulh x7,x8,x15
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------1
# Store result
	str x12,[x2,80]
# ------------------------------------------------------8
#  ------------------------------------7
	mul x3,x0,x1
	umulh x7,x0,x1
	adds x12,x3,x13
	adcs x13,x7,x14
	adc x14,xzr,xzr

#  ------------------------------------6
	mul x3,x16,x5
	umulh x7,x16,x5
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------5
	mul x3,x21,x9
	umulh x7,x21,x9
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------4
	mul x3,x17,x11
	umulh x7,x17,x11
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------3
	mul x3,x10,x20
	umulh x7,x10,x20
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------2
	mul x3,x8,x22
	umulh x7,x8,x22
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------1
	mul x3,x4,x15
	umulh x7,x4,x15
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

# Store result
	str x12,[x2,72]
# ------------------------------------------------------7
#  ------------------------------------6
	mul x3,x16,x1
	umulh x7,x16,x1
	adds x12,x3,x13
	adcs x13,x7,x14
	adc x14,xzr,xzr

#  ------------------------------------5
	mul x3,x21,x5
	umulh x7,x21,x5
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------4
	mul x3,x17,x9
	umulh x7,x17,x9
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------3
	mul x3,x10,x11
	umulh x7,x10,x11
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------2
	mul x3,x8,x20
	umulh x7,x8,x20
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------1
	mul x3,x4,x22
	umulh x7,x4,x22
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

# Store result
	str x12,[x2,64]
# ------------------------------------------------------6
#  ------------------------------------5
	mul x3,x21,x1
	umulh x7,x21,x1
	adds x12,x3,x13
	adcs x13,x7,x14
	adc x14,xzr,xzr

#  ------------------------------------4
	mul x3,x17,x5
	umulh x7,x17,x5
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------3
	mul x3,x10,x9
	umulh x7,x10,x9
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------2
	mul x3,x8,x11
	umulh x7,x8,x11
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------1
	mul x3,x4,x20
	umulh x7,x4,x20
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

# Store result
	str x12,[x2,56]
# ------------------------------------------------------5
#  ------------------------------------4
	mul x3,x17,x1
	umulh x7,x17,x1
	adds x12,x3,x13
	adcs x13,x7,x14
	adc x14,xzr,xzr

#  ------------------------------------3
	mul x3,x10,x5
	umulh x7,x10,x5
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------2
	mul x3,x8,x9
	umulh x7,x8,x9
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------1
	mul x3,x4,x11
	umulh x7,x4,x11
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

# Store result
	str x12,[x2,48]
# ------------------------------------------------------4
#  ------------------------------------3
	mul x3,x10,x1
	umulh x7,x10,x1
	adds x12,x3,x13
	adcs x13,x7,x14
	adc x14,xzr,xzr

#  ------------------------------------2
	mul x3,x8,x5
	umulh x7,x8,x5
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

#  ------------------------------------1
	mul x3,x4,x9
	umulh x7,x4,x9
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

# Store result
	str x12,[x2,40]
# ------------------------------------------------------3
#  ------------------------------------2
	mul x3,x8,x1
	umulh x7,x8,x1
	adds x12,x3,x13
	adcs x13,x7,x14
	adc x14,xzr,xzr

#  ------------------------------------1
	mul x3,x4,x5
	umulh x7,x4,x5
	adds x12,x3,x12
	adcs x13,x7,x13
	adc x14,x14,xzr

# Store result
#	str x12,[x2,32]
# ------------------------------------------------------2
#  ------------------------------------1
	mul x3,x4,x1
	umulh x7,x4,x1
	mov x20,v0.d[0]
	adds x15,x3,x13
	stp x15,x12,[x2,24]
	mov x21,v0.d[1]
	adcs x13,x7,x14
	adc x14,xzr,xzr

# Store result
#	str x12,[x2,24]
	stp x14,x13,[x2,8]
	mov x22,v1.d[0]
	ret
# --------------------------------------------------------------end of mulm----
	.balign	8
	.globl	mulm
# ------------------------------------------------------------------qclear-----
qclear:
	dup v0.16b,wzr
#	stp xzr,xzr,[x0]
#	stp xzr,xzr,[x0,16]
#	stp xzr,xzr,[x0,32]
#	stp xzr,xzr,[x0,48]
	stp q0,q0,[x0]
	stp q0,q0,[x0,32]
	ret
	.balign	8
	.globl qclear
# -----------------------------------------------------------------------------
qmov:
	ldp q0,q1,[x0]
	ldp q2,q3,[x0,32]
	stp q0,q1,[x1]
	stp q2,q3,[x1,32]
	ret
	.balign	8
	.globl qmov
# -----------------------------------------------------------------------------
qmovz:
# copy sign and exponent
	ldr x5,[x0]
# first position is zero
#	str xzr,[x1,8]
#copy the rest
	ldp x7,x8,[x0,8]
	stp x5,xzr,[x1]
	ldp x2,x3,[x0,24]
	stp x7,x8,[x1,16]
	stp x2,x3,[x1,32]
	ldp x5,x6,[x0,40]
	ldr x2,[x0,56]
	stp x5,x6,[x1,48]
	str x2,[x1,64]
	stp xzr,xzr,[x1,72]
#	str xzr,[x1,88]
	ret
	.balign	8
	.globl qmovz
# -----------------------------------------------------------------------------
	.globl pack
pack:
	ldr x2,[x0]
	ldp x3,x4,[x0,16]
	str x2,[x1]
	ldp x5,x6,[x0,32]
	stp x3,x4,[x1,8]
	stp x5,x6,[x1,24]
	ldp x3,x4,[x0,48]
	ldr x5,[x0,64]
	stp x3,x4,[x1,40]
	str x5,[x1,56]
	ret
	.balign	8
	.extern mdnorm
# -----------------------------------------------------------------------------
	.globl	cmpm
cmpm:
	ldp x4,x5,[x0,8]
	ldp x3,x6,[x1,8]
	cmp x4,x3
	bne	.LcmpNotEqual
	cmp x5,x6
	bne	.LcmpNotEqual
	ldp x4,x5,[x0,24]
	ldp x3,x6,[x1,24]
	cmp x4,x3
	bne	.LcmpNotEqual
	cmp x5,x6
	bne	.LcmpNotEqual
	ldp x4,x5,[x0,40]
	ldp x3,x6,[x1,40]
	cmp x4,x3
	bne	.LcmpNotEqual
	cmp x5,x6
	bne	.LcmpNotEqual
	ldp x4,x5,[x0,56]
	ldp x3,x6,[x1,56]
	cmp x4,x3
	bne	.LcmpNotEqual
	cmp	x5,x6
	bne	.LcmpNotEqual
	ldr x4,[x0,72]
	ldr x3,[x1,72]
	cmp x4,x3
	bne	.LcmpNotEqual
	mov	w0, 0
	ret
	.p2align 3
.LcmpNotEqual:
	mov	w0, -1
	csneg	w0, w0, w0, ls
	ret
	.balign	8
# -----------------------------------------------------------------------------
	.globl qequal
qequal:
	ldp x4,x5,[x0]
	ldp x3,x6,[x1]
	cmp x4,x3
	bne	.LqequalNotEqual
	cmp x5,x6
	bne	.LqequalNotEqual
	ldp x4,x5,[x0,16]
	ldp x3,x6,[x1,16]
	cmp x4,x3
	bne	.LqequalNotEqual
	cmp x5,x6
	bne	.LqequalNotEqual
	ldp x4,x5,[x0,32]
	ldp x3,x6,[x1,32]
	cmp x4,x3
	bne	.LqequalNotEqual
	cmp x5,x6
	bne	.LqequalNotEqual
	ldp x4,x5,[x0,48]
	ldp x3,x6,[x1,48]
	cmp x4,x3
	bne	.LqequalNotEqual
	cmp	x5,x6
.LqequalNotEqual:
	cset x0,eq
	ret
	.balign	8
# ---------------------------------------------------------------------mulin---
# void mulin(QELT y,Qfloatp b,QfloatAccump ac3);
# Multiply qfloat by a single number
# x2 = x1 * x0, ac3 = y*b.
# x0 number
# x1 Qfloatp b
# x2 QfloatAccump ac3
	.global	mulin
mulin:
	ldp	x9, x3,[x1, 48]
	mul	x5, x3, x0
	umulh	x7, x3, x0
	adds	x8, x5, xzr
	str	x8, [x2, 72]
	adc x8,x7,xzr

	mul	x5, x9, x0
	umulh	x7, x9, x0
	adds	x8, x5, x8
	str	x8, [x2, 64]
	ldp	x9,x3, [x1, 32]
	adc x8,x7,xzr

	mul	x5, x3, x0
	umulh	x7, x3, x0
	adds	x8, x5, x8
	str	x8, [x2, 56]
	adc x8,x7,xzr

	mul	x5, x9, x0
	umulh	x7, x9, x0
	adds	x8, x5, x8
	str	x8, [x2, 48]
	ldp	x9,x3, [x1, 16]
	adc x8,x7,xzr

	mul	x5, x3, x0
	umulh	x7, x3, x0
	adds	x8, x5, x8
	str	x8, [x2, 40]
	adc x8,x7,xzr

	mul	x5, x9, x0
	umulh	x7, x9, x0
	adds	x8, x5, x8
	str	x8, [x2, 32]
	ldp	x9,x3, [x1]
	adc x8,x7,xzr

	mul	x5, x3, x0
	umulh	x7, x3, x0
	adds	x8, x5, x8
	str	x8, [x2, 24]
	adc x8,x7,xzr
	stp	xzr, x8, [x2, 8]
	str x9,[x2]
	ret
	.balign	8
# -----------------------------------------------------------------------------
# This function should never use x14,x11 since it is called by mdnorm
	.globl	addbit
addbit:
	ldr x10,[x0,64]
	adds x13,x10,1
	str x13,[x0,64]
	b.cc .LaddbitExit

	ldp x15,x10,[x0,48]
	adcs x4,x10,xzr
	adcs x13,x15,xzr
	stp x13,x4,[x0,48]
	b.cc .LaddbitExit

	ldp x15,x10,[x0,32]
	adcs x4,x10,xzr
	adcs x13,x15,xzr
	stp x13,x4,[x0,32]
	b.cc .LaddbitExit

	ldp x15,x10,[x0,16]
	adcs x4,x10,xzr
	adcs x13,x15,xzr
	stp x13,x4,[x0,16]

	ldr x10,[x0,8]
	adcs x13,x10,xzr
	str x13,[x0,8]
.LaddbitExit:
	ret
	.balign	8
# -----------------------------------------------------------------------------
	.global	itoq
	.global lltoq
itoq:
lltoq:
	cmp	x0, 0
	beq	.LitoqIsZero
	stp	xzr, xzr, [x1, 16]
	stp	xzr, xzr, [x1, 32]
	stp	xzr, xzr, [x1, 48]
	csetm w2,lt
	cneg x0,x0,lt
	str w2,[x1]
	mov	w3, 0x80000
	clz x2,x0
	lsl	x0, x0, x2
	add w3,w3,64
	str	x0, [x1,8]
	sub	x0, x3, x2
	str	w0, [x1, 4]
	ret
	.p2align 3
.LitoqIsZero:
	mov	x0, x1
	b	qclear
# -----------------------------------------------------------------------------
# Uses v0 as storage for x30.
	.balign 8
	.global	normlz
normlz:
	mov	x9, x0
	ldr	x0, [x0, 8]
	mov	x15, x1
	cbnz	x0, .LnormlzAdjustDown
	ldr	x0, [x9, 16]
	mov	w14, 0
	.p2align 2
.LnormlzLoop:
	cmp	x0, 0
	bne	.LnormlzAdjustMSW
	ldp	x0, x6, [x9, 24]
	add	w14, w14, 64
	ldp	x5, x4, [x9, 40]
	cmp	w14, 640
	ldp	x3, x2, [x9, 56]
	stp	x0, x6, [x9, 16]
	ldr	x1, [x9, 72]
	stp	x5, x4, [x9, 32]
	stp	x3, x2, [x9, 48]
	stp	x1, xzr, [x9, 64]
	bne	.LnormlzLoop
	str	w14, [x15]
	mov	w0, 1
	ret
	.p2align 3
.LnormlzAdjustDown:
	clz x0,x0
	mov	w13, 64
	sub	w1, w13, w0
	mov	x0, x9
	neg	w2, w1
	str	w2, [x15]
#	mov v0.d[0],x30
#	bl	shiftdownn
	sub x5,x13,x1

	ldp x2,x3,[x0,8]
	lsl x4,x2,x5
	lsr x2,x2,x1
	lsl x7,x3,x5
	lsr x3,x3,x1
	orr x3,x3,x4
	stp x2,x3,[x0,8]

	ldp x2,x3,[x0,24]
	lsl x4,x2,x5
	lsr x2,x2,x1
	orr x2,x2,x7

	lsl x7,x3,x5
	lsr x3,x3,x1
	orr x3,x3,x4
	stp x2,x3,[x0,24]

	ldp x2,x3,[x0,40]
	lsl x4,x2,x5
	lsr x2,x2,x1
	orr x2,x2,x7
	lsl x7,x3,x5
	lsr x3,x3,x1
	orr x3,x3,x4
	stp x2,x3,[x0,40]

	ldp x2,x3,[x0,56]
	lsl x4,x2,x5
	lsr x2,x2,x1
	orr x2,x2,x7
	lsl x7,x3,x5
	lsr x3,x3,x1
	orr x3,x3,x4
	stp x2,x3,[x0,56]

	ldr x2,[x0,72]
	lsr x2,x2,x1
	orr x2,x2,x7
	str x2,[x0,72]
	mov	w0,wzr 
#	mov	x30,v0.d[0] 
	ret
	.p2align 3
.LnormlzAdjustMSW:
	blt	.LnormlzExit0
	clz x1,x0
	add	w14, w14, w1
	str	w14, [x15]

#	bl	shiftupn
	ldp x15,x8,[x9,64]
	mov x6,64
	sub x5,x6,x1
	lsr x4,x8,x5

	lsl x8,x8,x1
	lsr x7,x15,x5
	lsl x15,x15,x1
	orr x15,x15,x4
	stp x15,x8,[x9,64]
# Alternate between the two registers loaded, making the same operations
# in both to maximize parallelism
	ldp x3,x2,[x9,48]
	lsr x4,x2,x5
	lsr x6,x3,x5
	lsl x2,x2,x1
	lsl x3,x3,x1
	orr x2,x2,x7
	orr x3,x3,x4
	stp x3,x2,[x9,48]

	ldp x15,x8,[x9,32]
	lsr x4,x8,x5
	lsr x7,x15,x5
	lsl x8,x8,x1
	lsl x15,x15,x1
	orr x8,x8,x6
	orr x15,x15,x4
	stp x15,x8,[x9,32]

	ldp x3,x2,[x9,16]
	lsr x4,x2,x5
	lsr x6,x3,x5
	lsl x2,x2,x1
	lsl x3,x3,x1
	orr x2,x2,x7
	orr x3,x3,x4
	stp x3,x2,[x9,16]

	ldr x8,[x9,8]
	lsl x8,x8,x1
	orr x8,x8,x6
	str x8,[x9,8]
	mov	w0, 0
#	mov	x30,v0.d[0] 
	ret
	.p2align 3
.LnormlzExit0:
	str	w14, [x15]
	mov	w0, 0
	ret
	.balign 8
# -----------------------------------------------------------------------------
# Uses v1 as store for x30. Calls normlz that uses v0.
	.global	roundAccum
roundAccum:
	ldr	x3, [x0, 72]
	tbnz	x3, 63, .L430
	mov	w0, 0
	ret
	.p2align 3
.L430:
# SIGNBIT is 0x8000000000000000
#	mov	x4, -9223372036854775808
#	cmp	x3, x4
#	beq	.L431
	lsl x3,x3,1
	cbz x3,.L431
.L416:
	mov	v1.d[0], x0
	mov v1.d[1], x30
	bl	addbit
	sub sp,sp,16
	add	x1, sp,xzr 
	mov	x0, v1.d[0]
	bl	normlz
	ldr	x1, [sp],16
	cbz	w1, .L414
	mov	x0, v1.d[0]
	ldr	w3, [x0, 4]
	mov	x2, 1048575
	sub	x3, x3, x1, sxtw
	mov	w1, 1
	cmp	x3, x2
	bgt	.L414
	mov	w1, 0
	str	w3, [x0, 4]
.L414:
	mov x30,v1.d[1]
	mov	w0, w1
	ret
	.p2align 3
.L431:
	cbnz	w1, .L417
	ldr	x2, [x0, 80]
	tbnz	x2, 0, .L416
	mov w0,0
	ret
	.p2align 3
.L417:
	cbz	w2, .L416
	mov w0,0
	ret
# ----------------------------------------------------inverse_internal---------
	.balign 8
	.global inverse_internal
# Register map
#qa    --> v16, 16
#          v17, 32
#          v18, 48
#          v19  64
#          v20  80
#
#sqr   --> v26, 16
#          v27, 32
#          v28, 48
#          v29, 64
#          v30, 80
#
#prod  --> v3,  16
#          v4,  32
#          v5,  48
#          v6,  64
#          v7,  80
#
#quot:
#	v21.d[0] x21
#	v21.d[1] x22
#	v22.d[0] x23
#	v22.d[1] x24
#	v23.d[0] x25
#	v23.d[1] x26
#	v24.d[0] x27
#	v24.d[1] x28
inverse_internal:
	mov x1,1
	mov x22,xzr
	mov v31.d[0],x30
	lsl x1,x1,62
	mov	x2, v16.d[0]
	bl	udiv128by64
#	mov x21,x0
	mul	x4, x21, x21
	mov x23,xzr
	mov x24,xzr
	umulh	x0, x21, x21
	mov x25,xzr
	mov x26,xzr
	extr x11, x0, x4, 63
	mov x27,xzr
	mov x28,xzr
	lsl	x14, x4, 1
# input argument "sqr" is in x14 (lsb) and x11(msb)
	bl	Step2
	mov x6,xzr
	mov x9,xzr
	mov x3,v26.d[0]
	bl	mulvEntry1
	bl	submProdQuot
	bl	squareAsm
	bl	mulv
	bl	submProdQuot
	mov x30,v31.d[0]
	ret
# -----------------------------------------------------------------------------
	.balign 8
	.global	inverse
inverse:
# 1: Save registers in q21-q24, x1,x30 in q0
	sub sp,sp,96
	mov v21.d[0],x21
	mov v21.d[1],x22
	mov v22.d[0],x23
	mov v22.d[1],x24
	mov v23.d[0],x25
	mov v23.d[1],x26
	mov v24.d[0],x27
	mov v24.d[1],x28
	mov v0.d[0],x30
	mov v0.d[1],x1

# 2: Read data
	ldp q16,q17,[x0,16]
	ldp q18,q19,[x0,48]
# 3: Copy sign and exponent
	ldr x5,[x0]
	stp x5,xzr, [sp]
# 4: Calculate the inverse
	bl	inverse_internal
# 5: Save the result
	stp x21,x22,[sp,16]
	stp x23,x24,[sp,32]
	stp x25,x26,[sp,48]
	stp x27,x28,[sp,64]
# 6: Normalize and pack the result
	mov	x0, sp
	bl	mdnorm
	mov	x0, sp
	mov x1,v0.d[1]
	bl	pack
# 7: Restore the saved registers
	mov x21,v21.d[0]
	mov x22,v21.d[1]
	mov x23,v22.d[0]
	mov x24,v22.d[1]
	mov x25,v23.d[0]
	mov x26,v23.d[1]
	mov x27,v24.d[0]
	mov x28,v24.d[1]
	add sp,sp,96
	mov x30,v0.d[0]
	ret
# --------------------------------------------------------------udiv128by64----
	.balign 8
udiv128by64:
# Specialized routine for dividing 0x4000000000000000/a->mantissa[1]
# x0 is zero
# x1 contains 0x4000000000000000
# x2 a->mantissa[1]
# x3 is zero
# Result is returned in x21!
# x4 contains the high 32 bits of a->mantissa[1],x8 the lower 32 bits
	lsr x4,x2,32
	and x8,x2,0xffffffff
# Divide 0x4000000000000000 by high 32 bits of input a->mantissa[1]
	udiv x9,x1,x4
# x9 contains already the high 32 bits of result
	msub x11,x9,x4,x1
	mul x3,x8,x9
	lsl x10,x11,32
	cmp x3,x10
	b.ls .Ludiv128by64Suite
	add x10,x10,x2
	cmp x2,x10
	ccmp x3,x10,0,ls
	b.ls .L5
	sub x9,x9,2
	add x10,x10,x2
.Ludiv128by64Suite:
	sub x10,x10,x3
	udiv x1,x10,x4
	msub x10,x1,x4,x10
	mul x8,x8,x1
	lsl x10,x10,32
	cmp x8,x10
	b.ls .Ludiv128by64Exit
	add x10,x2,x10
	cmp x2,x10
	ccmp x8,x10,0,ls
	cinc x11,x1,ls
	sub x11,x11,2
.Ludiv128by64Exit:
	orr x21,x11,x9,lsl 32
	ret
	.p2align 2
.L5:
	sub x9,x9,1
	b .Ludiv128by64Suite
# ----------------------------------------------------------------mdnorm--------
	.balign 8
	.global	mdnorm
mdnorm:
# None of the functions called here uses x14, so we can use it instead of using
# expensive registers that need to be saved and restored
# x0 is used but never modified by all functions called, so we do not need to
# preserve it in calls
# All called function use x0 as an implicit argument
	mov v31.d[0],x30
	ldr	x5, [x0, 8]
	cbz	x5, .LmdnormNormalize
	mov	w14, 1048575
	b	.LfirstLoopEntryPoint
	.p2align 3
.LmdnormFirstLoop:
	ldr	x1, [x0, 8]
	str	w5, [x0, 4]
	cbz	x1, .LmdnormNormalize
.LfirstLoopEntryPoint:
	bl	shdn1
	ldr	w5, [x0, 4]
	cmp	w5, w14
	add	w5, w5, 1
	bls	.LmdnormFirstLoop
	mov	w5, 1048576
	str	w5, [x0, 4]
.LmdnormNormalize:
	ldr	x5, [x0, 16]
	cmp	x5, 0
	ble	.LtestIfRoundingNeeded
	clz x14,x5
	ldr	w5, [x0, 4]
	cmp	w5, w14
	bcs	.L23
.LtestIfRoundingNeeded:
	ldr	x5, [x0, 72]
	tbnz	x5, 63, .LdoAddBit
.L11:
	ldr	x5, [x0, 8]
	mov	w14, 1048575
	cbnz	x5, .L16
	mov x30,v31.d[0]
	ret
	.p2align 3
.L25:
	ldr	x1, [x0, 8]
	str	w5, [x0, 4]
	cbz	x1, .LmdnormExit
.L16:
	bl	shdn1
	ldr	w5, [x0, 4]
	cmp	w5, w14
	add	w5, w5, 1
	bls	.L25
	mov	w5, 1048576
	str	w5, [x0, 4]
.LmdnormExit:
	mov x30,v31.d[0]
	ret
	.p2align 3
.L23:
	mov	w1, w14
	bl	shiftupn
	ldr	w5, [x0, 4]
	sub	w14, w5, w14
	ldr	x5, [x0, 72]
	str	w14, [x0, 4]
	tbz	x5, 63, .L11
	.p2align 2
.LdoAddBit:
	bl	addbit
	b	.L11
# ----------------------------------------------------------------udiv128x64General----
# divides a 128 bit number by a 64 bit one
# x0,x1 contain the 64 bit number (high,low)
# x1 is reused in the next call
# x2 contains the 64 bit number
# Result (64 bit) of division is in x6
# Modulo in x1

	.balign 8
	.globl udiv128x64General
udiv128x64General:
	mov	x7, x2
	mov	x8, x0
	mov	x10, x1
# x3 is always zero
#	cbnz	x3, .LqmulOverflow 
# Since x2 is normalized, it is only possible to take this branch if the low 
# 64 bits of the number start also with  0x8...
	cmp	x2, x1
	b.ls .L2 
	clz	x3, x2
	sxtw	x9, w3
	cbz	x3, .L64 
	mov	w10, 64
	lsl	x1, x1, x9
	sub	w10, w10, w9
	lsl	x8, x0, x9
	lsl	x7, x2, x9
	lsr	x0, x0, x10
	orr	x10, x0, x1
.L64:
	lsr	x2, x7, 32
	and	x5, x7, 0xffffffff
	udiv	x1, x10, x2
	mov	x3, x1
	msub	x1, x1, x2, x10
	mul	x0, x5, x3
	extr	x1, x1, x8, 32
	cmp	x0, x1
	b.ls .L124 
	add	x1, x1, x7
	cmp	x0, x1
	ccmp	x7, x1, 2, hi
	b.hi .L984 
	sub	x3, x3, 2
	add	x1, x1, x7
.L124:
	sub	x1, x1, x0
	and	x0, x8, 0xffffffff
	udiv	x6, x1, x2
	msub	x2, x6, x2, x1
	mul	x5, x5, x6
	orr	x0, x0, x2, lsl 32
	cmp	x5, x0
	b.ls .L180 
	add	x0, x0, x7
	cmp	x7, x0
	ccmp	x5, x0, 0, ls
	b.ls .L960 
	sub	x6, x6, 2
	add	x0, x0, x7
.L180:
	sub	x5, x0, x5
# Set the result
	orr	x6, x6, x3, lsl 32
	lsr	x1, x5, x9
	ret
	.balign 8
.L2:
	cbnz	x2, .L324 
	mov	x7, 1                   	
	udiv	x7, x7, x2
.L324:
	clz x9,x7
	cbnz x9,.L472
	and	x10, x7, 0xffffffff
	sub	x2, x1, x7
	lsr	x5, x7, 32
	mov	x1, 1                   	
.L352:
	udiv	x6, x2, x5
	msub	x2, x6, x5, x2
	mov	x3, x6
	extr	x0, x2, x8, 32
	mul	x2, x6, x10
	cmp	x2, x0
	b.ls .L404 
	add	x0, x0, x7
	cmp	x2, x0
	ccmp	x7, x0, 2, hi
	b.hi .L976 
	sub	x3, x6, 2
	add	x0, x0, x7
.L404:
	sub	x0, x0, x2
	and	x2, x8, 0xffffffff
	udiv	x6, x0, x5
	msub	x0, x6, x5, x0
	mul	x5, x6, x10
	orr	x0, x2, x0, lsl 32
	cmp	x5, x0
	b.ls .L460 
	add	x0, x0, x7
	cmp	x7, x0
	ccmp	x5, x0, 0, ls
	b.ls .L968 
	sub	x6, x6, 2
	add	x0, x0, x7
.L460:
	sub	x5, x0, x5
	orr	x6, x6, x3, lsl 32
	lsr	x1, x5, x9
	ret
.L472:
	mov	x2, 64
	lsl	x7, x7, x9
	sub	x2, x2, x9
	lsl	x3, x1, x9
	lsr	x5, x7, 32
	and	x10, x7, 0xffffffff
	lsl	x8, x0, x9
	lsr	x1, x1, x2
	lsr	x2, x0, x2
	orr	x2, x2, x3
	udiv	x6, x1, x5
	msub	x1, x6, x5, x1
	mul	x0, x10, x6
	extr	x3, x1, x2, 32
	cmp	x0, x3
	b.ls .L560 
	add	x3, x3, x7
	cmp	x0, x3
	ccmp	x7, x3, 2, hi
	b.hi .L1016 
	sub	x6, x6, 2
	add	x3, x3, x7
.L560:
	sub	x3, x3, x0
	and	x2, x2, 0xffffffff
	udiv	x1, x3, x5
	msub	x3, x1, x5, x3
	mov	x0, x1
	mul	x1, x10, x1
	orr	x2, x2, x3, lsl 32
	cmp	x1, x2
	b.ls .L620 
	add	x2, x2, x7
	cmp	x1, x2
	ccmp	x7, x2, 2, hi
	b.hi .L992 
	sub	x0, x0, 2
	add	x2, x2, x7
.L620:
	sub	x2, x2, x1
	orr	x1, x0, x6, lsl 32
	b .L352 
.L632:
	mov	x8, 64
	lsl	x3, x3, x9
	sub	x8, x8, x9
	lsl	x11, x2, x9
	lsl	x5, x1, x9
	lsl	x10, x0, x9
	lsr	x2, x2, x8
	lsr	x1, x1, x8
	orr	x3, x2, x3
	lsr	x0, x0, x8
	orr	x5, x0, x5
	and	x12, x3, 0xffffffff
	lsr	x7, x3, 32
	udiv	x6, x1, x7
	msub	x1, x6, x7, x1
	mov	x2, x6
	extr	x0, x1, x5, 32
	mul	x1, x12, x6
	cmp	x1, x0
	b.ls .L736 
	add	x0, x0, x3
	cmp	x1, x0
	ccmp	x3, x0, 2, hi
	b.hi .L1008 
	sub	x2, x6, 2
	add	x0, x0, x3
.L736:
	sub	x0, x0, x1
	and	x5, x5, 0xffffffff
	udiv	x6, x0, x7
	msub	x0, x6, x7, x0
	orr	x5, x5, x0, lsl 32
	mul	x0, x12, x6
	cmp	x0, x5
	b.ls .L792 
	add	x5, x5, x3
	cmp	x3, x5
	ccmp	x0, x5, 0, ls
	b.ls .L1000 
	sub	x6, x6, 2
	add	x5, x5, x3
.L792:
	orr	x6, x6, x2, lsl 32
	and	x12, x11, 0xffffffff
	and	x1, x6, 0xffffffff
	lsr	x2, x11, 32
	sub	x5, x5, x0
	mov x11, 0x100000000
	lsr	x7, x6, 32
	mul	x2, x7, x2
	mul	x7, x7, x12
	mul	x12, x1, x12
	madd	x1, x1, x2, x7
	add	x0, x1, x12, lsr 32
	and	x12, x12, 0xffffffff
	cmp	x7, x0
	add	x1, x2, x11
	csel	x2, x1, x2, hi
	add	x2, x2, x0, lsr 32
	add	x0, x12, x0, lsl 32
	cmp	x5, x2
	b.cc .L936 
	ccmp	x10, x0, 2, eq
	mov	x1, x0
	b.cc .L936 
.L884:
# Prepare x1 for the next call to this function
	sub	x0, x10, x1
	cmp	x10, x0
	sbc	x5, x5, x2
	lsr	x0, x0, x9
	lsl	x8, x5, x8
	orr	x1, x8, x0
	ret
	.balign 8
.L936:
	sub	x1, x0, x11
	sub	x6, x6, 1
	cmp	x0, x1
	sbc	x2, x2, x3
	b .L884 
	.balign 8
.L960:
	sub	x6, x6, 1
	b .L180 
.L968:
	sub	x6, x6, 1
	b .L460 
.L976:
	sub	x3, x6, 1
	b .L404 
.L984:
	sub	x3, x3, 1
	b .L124 
.L992:
	sub	x0, x0, 1
	b .L620 
.L1000:
	sub	x6, x6, 1
	b .L792 
.L1008:
	sub	x2, x6, 1
	b .L736 
.L1016:
	sub	x6, x6, 1
	b .L560 
# ----------------------------------------------------------------------divi----
# Divide a Qfloat by a single number of 64 bits
# void divi(unsigned long long src, QfloatAccump b)
# b = b/src
# x0 contains the number
# x1 points to the qfloat
	.balign 8
	.global	divi
divi:
	mov v31.d[0],x30
	mov	x4, x0
	mov	x15, x2

# Start of inlined version of shiftdown2
	ldp x3,x5,[x1,8]
	ldp x6,x7,[x1,24]
	ldp x8,x9,[x1,40]
	ldp x10,x11,[x1,56]

	extr x1,xzr,x3,2
	extr x0,x3,x5,2

	extr x14,x5,x6,2
	extr x16,x6,x7,2
	mov v22.d[0],x14
	mov v22.d[1],x16
# Registers x4, x13, x14, x15, x16, and x17  are NOT used by udiv128x64General
# x4 is used to hold the 64 bit number, x1 for the modulo, x15 for the
# output pointer
# We can use x14, x17, x16 and x13  to store the results of the shifted accumulator
# accumulator instead of saving them somewhere

	extr x14,x7,x8,2
	extr x17,x8,x9,2

	extr x16,x9,x10,2
	extr x13,x10,x11,2
# end of shiftdown2
	mov	x2, x4
	bl	udiv128x64General
	stp	xzr,x6, [x15, 8]

	mov x0,v22.d[0]
	mov	x2, x4
	bl	udiv128x64General
	str	x6, [x15, 24]

	mov x0,v22.d[1]
	mov	x2, x4
	bl	udiv128x64General
	str	x6, [x15, 32]

	mov x0,x14
	mov	x2, x4
	bl	udiv128x64General
	str	x6, [x15, 40]

	mov x0,x17
	mov	x2, x4
	bl	udiv128x64General
	str	x6, [x15, 48]

	mov x0,x16
	mov	x2, x4
	bl	udiv128x64General
	str	x6, [x15, 56]

	mov x0,x13
	mov	x2, x4
	bl	udiv128x64General
	str	x6, [x15, 64]

	mov x0,xzr
	mov	x2, x4
	bl	udiv128x64General
	str	x6, [x15, 72]

	mov x30,v31.d[0]
	ret
# ----------------------------------------------------------------------qcmp----
	.balign 8
	.global	qcmp
qcmp:
	ldr	w4, [x0, 4]
	cmp	w4, 448
	bhi	.L232
	ldr	w2, [x1, 4]
	cmp	w2, 448
	bls	.L250
.L232:
	ldr	w3, [x0]
	ldr	w2, [x1]
	cmp	w3, w2
	beq	.L235
	mov	w2, -1
	cmp	w3, 0
	csneg	w0, w2, w2, ne
	ret
	.p2align 3
.L235:
	ldr	w5, [x1, 4]
	cmp	w3, 0
	mov	w2, -1
	csneg	w2, w2, w2, ne
	cmp	w4, w5
	beq	.L237
	csneg	w0, w2, w2, hi
	ret
	.p2align 3
.L250:
	str	x30, [sp, -80]!
	mov	x3, x1
	mov	x1, x0
	mov	x0, x3
	mov	w3, 1
	add	x2, sp, 16
	bl	qadd_subtract
	ldr	w0, [sp, 20]
	mov	w2, 0
	cbz	w0, .L231
	ldr	w0, [sp, 16]
	mov	w2, -1
	cmp	w0, 0
	csneg	w2, w2, w2, ne
.L231:
	mov	w0, w2
	ldr	x30, [sp], 80
	ret
	.p2align 3
.L237:
	ldp	x5, x6, [x1, 8]
	ldp	x4, x7, [x0, 8]
	cmp	x5, x4
	bne	.L238
	cmp x6,x7
	bne .L238

	ldp	x5, x6, [x1, 24]
	ldp	x4, x7, [x0, 24]
	cmp	x5, x4
	bne	.L238
	cmp x6,x7
	bne .L238

	ldp	x5, x6, [x1, 40]
	ldp	x4, x7, [x0, 40]
	cmp	x5, x4
	bne	.L238
	cmp x6,x7
	bne .L238

	ldr	x5, [x1, 56]
	ldr	x4, [x0, 56]
	cmp	x5, x4
	bne	.L238

	mov	w0, 0
	ret
	.p2align 3
.L238:
	csneg	w0, w2, w2, cc
	ret
# ----------------------------------------------------------------------qmul----
	.balign 8
	.global	qmul
qmul:
	str	x30, [sp, -144]!
	mov v20.d[1],x20
	mov v20.d[0],x19
	mov v25.d[0],x2
	mov	x20, x0
	ldr	w3, [x0, 4]
	cbz	w3, .LqmulReturnZero
	ldr	w3, [x1, 4]
	mov	x19, x1
	cbnz	w3, .LqmulIsNotZero
.LqmulReturnZero:
	mov x0,v25.d[0]
	bl	qclear
	b .LqmulExit
	.p2align 3
.LqmulIsNotZero:
	ldp	x3,x4, [x0, 16]
	add x2,sp,48
	orr x3,x3,x4
	cbnz	x3, .LqmulTestSecondOperandSmall
	ldp	x3,x4, [x0, 32]
	orr x3,x3,x4
	cbnz	x3, .LqmulTestSecondOperandSmall
	ldp	x3, x4, [x0, 48]
	orr x3,x3,x4
	cbnz	x3, .LqmulTestSecondOperandSmall
#	ldr	x3, [x0, 40]
#	cbnz	x3, .LqmulTestSecondOperandSmall
#	ldr	x3, [x0, 48]
#	cbnz	x3, .LqmulTestSecondOperandSmall
#	ldr	x3, [x0, 56]
#	cbnz	x3, .LqmulTestSecondOperandSmall
	ldr	x0, [x0, 8]
	bl	mulin
	add	x0, sp,48
	bl	mdnorm
	ldr	w2, [sp, 52]
	ldr	w0, [x20, 4]
	sub	w2, w2, #1048576
	add	w2, w2, w0
	b	.LqmulSetSignAndExponent
	.p2align 3
.LqmulTestSecondOperandSmall:
	ldr	x3, [x19, 16]
	cbnz	x3, .LqmulDoLongMultiplication
	ldr	x3, [x19, 24]
	cbz	x3, .L227
.LqmulDoLongMultiplication:
	bl	mulm
.LqmulNormalize:
	add	x0, sp,48
	bl	mdnorm
	ldr	w2, [sp, 52]
	ldr	w0, [x19, 4]
	sub	w2, w2, #1048576
	add	w2, w2, w0
.LqmulSetSignAndExponent:
	ldr	w0, [x20]
	add	w2, w2, 524288
	ldr	w3, [x19]
	mov	w1, 1048575
	cmp	w3, w0
	csetm	w0, ne
	cmp	w2, w1
	str	w0, [sp, 48]
	bgt	.LqmulOverflow
	cmp	w2, 0
	ble	.LqmulReturnZero
	mov	x1, v25.d[0]
	add	x0, sp,48
	str	w2, [sp, 52]
	bl	pack
.LqmulExit:
	mov x19,v20.d[0]
	mov x20,v20.d[1]
#	mov x22,v25.d[0]
	ldr	x30, [sp], 144
	ret
	.p2align 3
.LqmulOverflow:
	mov	x0, v25.d[0]
	bl qinfin
	b .LqmulExit
	.p2align 3
.L227:
	ldr	x3, [x19, 32]
	cbnz	x3, .LqmulDoLongMultiplication
	ldr	x3, [x19, 40]
	cbnz	x3, .LqmulDoLongMultiplication
	ldr	x3, [x19, 48]
	cbnz	x3, .LqmulDoLongMultiplication
	ldr	x3, [x19, 56]
	cbnz	x3, .LqmulDoLongMultiplication
	ldr	x0, [x19, 8]
	mov	x1, x20
	bl	mulin
	b	.LqmulNormalize
# -----------------------------------------------------------------------------
# void divm(Qfloatp a, Qfloatp b,QfloatAccump ac3)
# Input: x0, x1
# Output: sp+64 ... sp+128
divm:
# 1: Save registers in q21-q24
	ldr q16,[x0,8]
	mov v21.d[0],x21
	mov v21.d[1],x22
	mov v22.d[0],x23
	mov v22.d[1],x24
	ldr q17,[x0,24]
	mov v23.d[0],x25
	mov v23.d[1],x26
	mov v24.d[0],x27
	mov v24.d[1],x28
	ldr q18,[x0,40]
	str x1,[sp,16]

# 2: Read the data
	ldr x3,[x0,56]
	mov v19.d[0],x3
	mov v19.d[1],xzr
# 3: Calculate the inverse
	bl	inverse_internal
# 4: Move the inverse from x21...x28 to q16-q19, and read the second arg into q26-q29
	ldr x1,[sp,16]
	mov v16.d[0],x21
	mov v16.d[1],x22
	mov v17.d[0],x23
	ldr q27,[x1,24]
	mov v17.d[1],x24
	ldr q28,[x1,40]
	mov v29.d[1],xzr
	mov v18.d[0],x25
	ldr q29,[x1,56]
	mov v18.d[1],x26
	mov v19.d[0],x27
	ldr q26,[x1,8]
#	mov v19.d[1],x28
# 5 multiply by the inverse
#	mov x6,v19.d[1]
	mov x3,v26.d[0]
	mul x9,x28,x3
	umulh x12,x28,x3

	mov x10,x21
	mov x11,x22
	
	mov x0,x23
	mov x1,x24
	mov x6,x27

	bl	mulvEntry1
# 6 Store the result and restore the saved registers
#	add x2,sp,64
	mov x21,v21.d[0]
	mov x22,v21.d[1]
	mov x23,v22.d[0]
	stp q3,q4,[sp,64+16]
	mov x24,v22.d[1]
	mov x25,v23.d[0]
	mov x26,v23.d[1]
	mov x27,v24.d[0]
	stp q5,q6,[sp,64+48]
	mov x28,v24.d[1]
	str xzr,[sp,64+8]
# Done
#	bl	divm
.qdivSetupSignAndExponent:
	mov	w1, 4
	add	x0, sp,64
	str	w1, [sp, 68]
	bl	mdnorm
	mov x4,v20.d[0]
	ldr	w2, [sp, 68]
	ldr	w0, [x21]
	ldp	w3, w1, [x4]
	ldr w9,[x21,4]
	add	x2, x2, x9
	cmp	w3, w0
	sub	x2, x2, x1
	mov	x3, 524286
	mov	x1, 1048575
	add	x2, x2, x3
	csinv	w0, wzr, wzr, eq
	cmp	x2, x1
	str	w0, [sp, 64]
	bgt	.LqdivOverflow
	cmp	x2, 0
	ble	.LqdivUnderflow
	str	w2, [sp, 68]
	mov	x1, v20.d[1]
	add	x0, sp,64
	bl	pack
.LqdivExit:
	ldp	x30, x21, [sp], 160
	ret
.LqdivUnderflow:
	mov	x0, v20.d[1]
	bl	qclear
	b .LqdivExit
# ----------------------------------------------------------------------qdiv----
	.balign 8
	.global qdiv
qdiv:
	stp	x30,x21, [sp, -160]!
	mov v20.d[0],x0
	mov v20.d[1],x2
#	mov	x20, x2
	ldr	w9, [x1, 4]
	cbz	w9, .LqdivUnderflow
	ldr	w6, [x0, 4]
	cbz	w6, .LqdivOverflow
.LqdivTestForSmallDivisor:
	add	x2, sp,64
	ldp	x6,x5, [x0, 16]
	mov	x21, x1
	orr x6,x6,x5
	cbnz	x6, divm
	ldp	x6,x5, [x0, 32]
	orr x6,x6,x5
	cbnz	x6, divm
	ldp	x5,x6, [x0, 48]
	orr x5,x6,x5
	cbnz	x5, divm
	ldr	x0, [x0, 8]
	bl	divi
	b	.qdivSetupSignAndExponent
.LqdivOverflow:
	mov x0,v20.d[1]
	bl qinfin
	b .LqdivExit
