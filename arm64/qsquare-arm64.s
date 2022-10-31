	.section    __TEXT,__text,regular,pure_instructions
	.p2align 3
	.globl _qsquare
_qsquare:                               
	sub	sp, sp, #240            
	stp	x20, x19, [sp, #208]    
	stp	x29, x30, [sp, #224]    
	add	x29, sp, #224           
	mov	x19, x1
	ldr	w7, [x0, #4]
	cbz	w7, LBB_ReturnZero
	mov	x20, x0
	ldp	x8,x9, [x0, #16]
	orr x9,x8,x9
	cbnz	x9, LBBNotASmallNumber
	ldp	x8,x9, [x20, #32]
	orr x9,x8,x9
	cbnz	x9, LBBNotASmallNumber
	ldp	x8,x9, [x20, #48]
	orr x9,x8,x9
	cbz	x9, LBB_IsSmallNumber
LBBNotASmallNumber:
	mov	x1, sp                 
	bl	_qmovz                 
	ldr	w8, [x20, #4]
	str	w8, [sp, #100]
	mov	x0, sp
	add	x1, sp, #96             
	mov	w2, #8
	bl	_square
LBB_Normalize:
	add	x0, sp, #96             
	bl	_mdnorm
	ldr	w8, [x20, #4]
	ldr	w9, [sp, #100]
	add	x8, x9, x8              
	str	wzr, [sp, #96]          
	cmp	x8, #384, lsl #12       
	b.ge	LBB_Overflow
LBB_NotOverflow:
	mov	w9, #2
	movk	w9, #8, lsl #16
	cmp	x8, x9
	b.hs	LBB_NotUnderflow
LBB_ReturnZero:
	mov	x0, x19
	bl	_qclear
	b	LBBSquareExitProc
LBB_NotUnderflow:
	mov	w9, #-524289
	add	w8, w8, w9
	str	w8, [sp, #100]
	add	x0, sp, #96             
	mov	x1, x19
	bl	_pack
LBBSquareExitProc:
	ldp	x29, x30, [sp, #224]    
	ldp	x20, x19, [sp, #208]    
	add	sp, sp, #240
	ret
	.p2align 3
LBB_IsSmallNumber:
	ldr	x9, [x20, #8]
	mul	x10, x9, x9
	mov v0.d[0],xzr
	mov v0.d[1],xzr
#	movi.2d	v0, #0000000000000000
	stp	q0, q0, [sp, #96]
	umulh	x9, x9, x9
	stp	q0, q0, [sp, #160]
	stp	q0, q0, [sp, #128]
	stp	x9, x10, [sp, #112]     
	add	w8, w7, #1              
	str	w8, [sp, #100]          
	b	LBB_Normalize
LBB_Overflow:
	bl _qinfin
	b	LBBSquareExitProc
