	.section    __TEXT,__text,regular,pure_instructions
	.p2align 3
	.globl	_shift
_shift:
#	dup v0.2d,xzr
# w20 toshift
	tbz w1, #31, L_SC_MORE_THAN_ZERO
# -------------------------------------------------sc < 0
	cmn w1, #63
	neg w8, w1
	bge L27
	cmn w1, #576
	blt L_ERROR
	ldp q0,q1,[x0,16]
	asr w10, w8, 6
	and w8,w8,63
	ldp q2,q3,[x0,48]
	cmp w10,1
	ldp q2,q3,[x0,48]
	beq L_DO_1
	cmp w10,2
	beq L_DO_2
	cmp w10,3
	beq L_DO_3
	cmp w10,4
	beq L_DO_4
	cmp w10,5
	beq L_DO_5
	cmp w10,6
	beq L_DO_6
	cmp w10,7
	beq L_DO_7
	b L_ERROR
	.p2align 3
L_DO_1:
	mov x10,v3.d[1]
	ext v3.16b,v2.16b,v3.16b,8
	ext v2.16b,v1.16b,v2.16b,8
	ext v1.16b,v0.16b,v1.16b,8
	ins v0.d[1],v0.d[0]
	mov v0.d[0],xzr
	b L16
	.p2align 3
L_DO_7:
	orr v0.16b,v0.16b,v1.16b
	orr v0.16b,v0.16b,v2.16b
	orr v0.16b,v0.16b,v3.16b
# Finished storing the 64 bits shifts. Look if we have to
# make a shift by less than 64 bits
	mov x2,v0.d[0]

	dup v0.2d,xzr
	dup v1.2d,xzr
	dup v2.2d,xzr
	dup v3.2d,xzr

	stp q0,q1,[x0,16]
	stp q2,q3,[x0,48]
	cbz w8,L100
# Only the last position needs to be shifted since all others are
# zero: we have just written those zeroes above
#	ldr x2, [x0, 72]
	lsr x9,x2,x8
	stp xzr, x9, [x0, 64]
	orr x0,x10,x2
	ret
	.p2align 3
L100:
	stp xzr, x2, [x0, 64]
	cmp x10,xzr
	cset w0, ne
	ret
	.p2align 3
L_DO_2:
	mov v3.16b,v2.16b
	mov v2.16b,v1.16b
	mov v1.16b,v0.16b
	addp d0,v0.2d
	mov x10,v0.d[0]
	dup v0.2d,xzr
	b L16
	.p2align 3
L_DO_3:
	addp d3,v3.2d
	mov x10,v2.d[1]
	mov x11,v3.d[0]
	orr x10,x10,x11

	ext v3.16b,v1.16b,v2.16b,8

	ext v2.16b,v0.16b,v1.16b,8

	ins v1.d[1],v0.d[0]

	mov v1.d[0],xzr
	dup v0.2d,xzr

	b L16
	.p2align 3
L_DO_4:
	mov v3.16b,v1.16b
	mov v2.16b,v0.16b

	orr  v0.16b,v0.16b,v1.16b
	cmeq d0,d0,0
	not v0.16b,v0.16b
	mov x10,v0.d[0]

	dup v0.2d,xzr
	dup v1.2d,xzr

	b L16
	.p2align 3
L_DO_5:
	addp d3,v3.2d
	addp d2,v2.2d
	mov x10,v3.d[0]
	mov x11,v2.d[0]
	orr x10,x10,x11

	ext v3.16b,v0.16b,v1.16b,8

	ins v2.d[1],v0.d[0]
	mov v2.d[0],xzr

	dup v1.2d,xzr
	dup v0.2d,xzr

	b L16
	.p2align 3
L_DO_6:
	mov v3.16b,v0.16b

	orr v0.16b,v0.16b,v1.16b
	orr v0.16b,v0.16b,v2.16b
	orr v0.16b,v0.16b,v3.16b

	cmeq d0,d0,0
	not v0.16b,v0.16b
	mov x10,v0.d[0]

	dup v0.2d,xzr
	dup v1.2d,xzr
	dup v2.2d,xzr

L16:
	stp q0,q1,[x0,16]
	stp q2,q3,[x0,48]
	cbnz w8,L11
	cmp x10,xzr
	cset w0, ne
	ret
	.p2align 3
L27:
	ldp q0,q1,[x0,16]
	ldp q2,q3,[x0,48]
	mov w10, 0
L11:
	eor x5,x8,63
	ldr x2, [x0, 72]
# shiftdownn

	add x5,x5,1
	mov x3,v0.d[0]
	lsl x7,x3,x5
	lsr x3,x3,x8
	stp xzr,x3,[x0,8]

#	ldp x12,x14,[x0,24]
	mov x12,v0.d[1]
	lsl x4,x12,x5
	lsr x12,x12,x8
	mov x14,v1.d[0]
	orr x12,x12,x7

	lsl x7,x14,x5
	lsr x14,x14,x8
	orr x14,x14,x4
	stp x12,x14,[x0,24]

#	ldp x9,x3,[x0,40]
	mov x9,v1.d[1]
	mov x3,v2.d[0]
	lsl x4,x9,x5
	lsr x9,x9,x8
	orr x9,x9,x7
	lsl x7,x3,x5
	lsr x3,x3,x8
	orr x3,x3,x4
	stp x9,x3,[x0,40]

#	ldp x12,x14,[x0,56]
	mov x12,v2.d[1]
	mov x14,v3.d[0]
	lsl x4,x12,x5
	lsr x12,x12,x8
	orr x12,x12,x7
	lsl x7,x14,x5
	lsr x14,x14,x8
	orr x14,x14,x4
	stp x12,x14,[x0,56]

#	ldr x9,[x0,72]
	lsr x9,x2,x8
	orr x9,x9,x7
	str x9,[x0,72]
	orr x0,x10,x2
	ret
	.p2align 3
L_ERROR:
# sc >= NBITS. This is an error condition
#	add x0, x0, 8
	mov x1, -9223372036854775808
#	mov w10, 1
	dup v0.2d,xzr
	dup v1.2d,xzr
	dup v2.2d,xzr
	dup v3.2d,xzr
	stp q0,q1,[x0]
	stp q0,q1,[x0,32]
	stp q0,q1,[x0,64]
	str x1, [x0, 16]
	mov w0,1
	ret
	.p2align 3
Lret0:
	mov w0, 0
	ret
	.p2align 3
#--------------------------------------------------------- sc > 0
L_SC_MORE_THAN_ZERO:
	cbz w1, Lret0
	stp x29, x30, [sp, -64]!
	add x29, sp, 0
	stp x21, x22, [sp, 32]
	mov w22, w1
	str x23, [sp, 48]
	mov x23, x0
	stp x19, x20, [sp, 16]
	cmp w1, 63
	ble L29
	cmp w1, 447
	bgt L_ERROR
	asr w21, w1, 6
	mov w2, 8
	sub w2, w2, w21
	add x0, x0, 16
	sxtw    x19, w21
	sbfiz   x2, x2, 3, 32
	add x1, x19, 2
	and w22, w22, 63
	mov w20, 0
	add x1, x23, x1, lsl 3
	bl  _memmove
	neg x0, x19, lsl 3
	add x0, x0, 80
	add x0, x23, x0
L26:
	subs    w21, w21, #1
	ldr x1, [x0]
	str xzr, [x0], 8
	orr w20, w20, w1
	bne L26
L24:
	cbz w22, L9
	mov w1, w22
	mov x0, x23
	bl  _shiftupn
L9:
	mov w0, w20
	ldr x23, [sp, 48]
	ldp x19, x20, [sp, 16]
	ldp x21, x22, [sp, 32]
	ldp x29, x30, [sp], 64
	ret
	.p2align 3
L29:
	mov w20, 0
	b   L24
	.pool
#	.section __TEXT,__cstring,cstring_literals
#	.p2align 3
#Lswitch:
#	.quad L_ERROR
#	.quad L_DO_1
#	.quad L_DO_2
#	.quad L_DO_3
#	.quad L_DO_4
#	.quad L_DO_5
#	.quad L_DO_6
#	.quad L_DO_7
#	.quad L_ERROR
#	.quad L_ERROR
#	.quad L_ERROR
#	.quad L_ERROR
#	.quad L_ERROR
#	
