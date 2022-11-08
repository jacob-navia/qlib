



	.file	"qasm.asm"
	.text




	.globl ___shdn1
___shdn1:
	clc
	rcrq $1,8(%rdi)
	rcrq $1,16(%rdi)
	rcrq $1,24(%rdi)
	rcrq $1,32(%rdi)
	rcrq $1,40(%rdi)
	rcrq $1,48(%rdi)
	rcrq $1,56(%rdi)
	rcrq $1,64(%rdi)
	rcrq $1,72(%rdi)
	ret
	.align	2



	.globl ___shup1
___shup1:
	clc
	rclq $1,72(%rdi)
	rclq $1,64(%rdi)
	rclq $1,56(%rdi)
	rclq $1,48(%rdi)
	rclq $1,40(%rdi)
	rclq $1,32(%rdi)
	rclq $1,24(%rdi)
	rclq $1,16(%rdi)
	rclq $1,8(%rdi)
	ret
	.align	2





	.align	2
	.globl	___shiftupn
___shiftupn:
	movq	%rsi,%rcx

	movq 8(%rdi),%rsi
	movq 16(%rdi),%rax
	shldq %cl,%rax,%rsi
	movq %rsi,8(%rdi)

	movq 24(%rdi),%rsi
	shldq %cl,%rsi,%rax
	movq %rax,16(%rdi)

	movq 32(%rdi),%rax
	shldq %cl,%rax,%rsi
	movq %rsi,24(%rdi)

	movq 40(%rdi),%rsi
	shldq %cl,%rsi,%rax
	movq %rax,32(%rdi)

	movq 48(%rdi),%rax
	shldq %cl,%rax,%rsi
	movq %rsi,40(%rdi)

	movq 56(%rdi),%rsi
	shldq %cl,%rsi,%rax
	movq %rax,48(%rdi)

	movq 64(%rdi),%rax
	shldq %cl,%rax,%rsi
	movq %rsi,56(%rdi)

	movq 72(%rdi),%rsi
	shldq %cl,%rsi,%rax
	movq %rax,64(%rdi)

	shlq %cl,%rsi
	movq %rsi,72(%rdi)

	movq	%rcx,%rax
	ret
	.align	2



	.align	2
	.globl	___shiftdownn
___shiftdownn:
	movq	%rsi,%rcx
	movq	72(%rdi),%rax

	movq 64(%rdi),%rsi
	shrd %cl,%rsi,%rax
	movq %rax,72(%rdi)

	movq 56(%rdi),%rax
	shrd %cl,%rax,%rsi
	movq %rsi,64(%rdi)

	movq 48(%rdi),%rsi
	shrd %cl,%rsi,%rax
	movq %rax,56(%rdi)

	movq 40(%rdi),%rax
	shrd %cl,%rax,%rsi
	movq %rsi,48(%rdi)

	movq 32(%rdi),%rsi
	shrd %cl,%rsi,%rax
	movq %rax,40(%rdi)

	movq 24(%rdi),%rax
	shrd %cl,%rax,%rsi
	movq %rsi,32(%rdi)

	movq 16(%rdi),%rsi
	shrd %cl,%rsi,%rax
	movq %rax,24(%rdi)

	movq 8(%rdi),%rax
	shrd %cl,%rax,%rsi
	movq %rsi,16(%rdi)

	shrq %cl,%rax
	movq %rax,8(%rdi)
	movl %ecx,%eax
	ret
	.align	2



	.globl	_bsr64
_bsr64:
	bsrq	%rdi,%rax
	ret

	.align	2




	.globl ___addm
___addm:
	movq 72(%rsi),%rax
	addq 72(%rdi),%rax
	movq %rax,72(%rsi)

	movq 64(%rsi),%rax
	adcq 64(%rdi),%rax
	movq %rax,64(%rsi)

	movq 56(%rsi),%rax
	adcq 56(%rdi),%rax
	movq %rax,56(%rsi)

	movq 48(%rsi),%rax
	adcq 48(%rdi),%rax
	movq %rax,48(%rsi)

	movq 40(%rsi),%rax
	adcq 40(%rdi),%rax
	movq %rax,40(%rsi)

	movq 32(%rsi),%rax
	adcq 32(%rdi),%rax
	movq %rax,32(%rsi)

	movq 24(%rsi),%rax
	adcq 24(%rdi),%rax
	movq %rax,24(%rsi)

	movq 16(%rsi),%rax
	adcq 16(%rdi),%rax
	movq %rax,16(%rsi)

	movq 8(%rsi),%rax
	adcq 8(%rdi),%rax
	movq %rax,8(%rsi)

	ret
	.align	2




	.globl ___subm
___subm:
	movq 72(%rsi),%rax
	subq 72(%rdi),%rax
	movq %rax,72(%rsi)

	movq 64(%rsi),%rax
	sbbq 64(%rdi),%rax
	movq %rax,64(%rsi)

	movq 56(%rsi),%rax
	sbbq 56(%rdi),%rax
	movq %rax,56(%rsi)

	movq 48(%rsi),%rax
	sbbq 48(%rdi),%rax
	movq %rax,48(%rsi)

	movq 40(%rsi),%rax
	sbbq 40(%rdi),%rax
	movq %rax,40(%rsi)

	movq 32(%rsi),%rax
	sbbq 32(%rdi),%rax
	movq %rax,32(%rsi)

	movq 24(%rsi),%rax
	sbbq 24(%rdi),%rax
	movq %rax,24(%rsi)

	movq 16(%rsi),%rax
	sbbq 16(%rdi),%rax
	movq %rax,16(%rsi)

	movq 8(%rsi),%rax
	sbbq 8(%rdi),%rax
	movq %rax,8(%rsi)

	ret
	.align	2
	.globl	___pack




___pack:
	movq	(%rdi),%rax
	movq	%rax,(%rsi)

	movupd	16(%rdi),%xmm0
	movupd	%xmm0,8(%rsi)


	movupd	32(%rdi),%xmm0
	movupd	%xmm0,24(%rsi)

	movupd	48(%rdi),%xmm0
	movupd	%xmm0,40(%rsi)

	movq	64(%rdi),%rax
	movq	%rax,56(%rsi)

	ret
	.align	2
	.globl	___qmovz




___qmovz:
	movq	(%rdi),%rax
	movq	%rax,(%rsi)
	movq	$0,8(%rsi)
	movupd	8(%rdi),%xmm0
	movupd	%xmm0,16(%rsi)
	movupd	24(%rdi),%xmm0
	movupd	%xmm0,32(%rsi)
	movupd	40(%rdi),%xmm0
	movupd	%xmm0,48(%rsi)

	movq	56(%rdi),%rax
	movq	%rax,64(%rsi)
	movq	$0,72(%rsi)
	ret
	.align	2




	.globl	___qmov
___qmov:
	movupd	(%rdi),%xmm0
	movupd	16(%rdi),%xmm1
	movupd	%xmm0,(%rsi)
	movupd	%xmm1,16(%rsi)
	movupd	32(%rdi),%xmm0
	movupd	48(%rdi),%xmm1
	movupd	%xmm0,32(%rsi)
	movupd	%xmm1,48(%rsi)
	ret
	.align	2
	.globl ___qclear




___qclear:
	xorpd	%xmm0,%xmm0
	movupd	%xmm0,(%rdi)
	movupd	%xmm0,16(%rdi)
	movupd	%xmm0,32(%rdi)
	movupd	%xmm0,48(%rdi)
	ret
	.globl ___mulin
	.align	2




___mulin:

	xorpd	%xmm1,%xmm1
	subq	$80,%rsp
	movupd	%xmm1,(%rsp)
	movupd	%xmm1,16(%rsp)
	movupd	%xmm1,32(%rsp)
	movupd	%xmm1,48(%rsp)
	movupd	%xmm1,64(%rsp)


	movq	(%rsi),%rax
	movq	%rax,(%rsp)
	movq	8(%rdi),%rdi
	xorl	%ecx,%ecx

	movq	64(%rsi),%rax
	orq	%rax,%rax
	je	L_skip1
	mulq	%rdi
	movq	%rax,72(%rsp)
	movq	%rdx,64(%rsp)

L_skip1:
	movq	56(%rsi),%rax
	orq	%rax,%rax
	je	L_skip2
	mulq	%rdi
	addq	%rax,64(%rsp)
	adcq	%rdx,56(%rsp)
	adcq	%rcx,48(%rsp)

L_skip2:
	movq	48(%rsi),%rax
	orq	%rax,%rax
	je	L_skip3
	mulq	%rdi
	addq	%rax,56(%rsp)
	adcq	%rdx,48(%rsp)
	adcq	%rcx,40(%rsp)

L_skip3:
	movq	40(%rsi),%rax
	orq	%rax,%rax
	je	L_skip4
	mulq	%rdi
	addq	%rax,48(%rsp)
	adcq	%rdx,40(%rsp)
	adcq	%rcx,32(%rsp)

L_skip4:
	movq	32(%rsi),%rax
	orq	%rax,%rax
	je	L_skip5
	mulq	%rdi
	addq	%rax,40(%rsp)
	adcq	%rdx,32(%rsp)
	adcq	%rcx,24(%rsp)

L_skip5:
	movq	24(%rsi),%rax
	orq	%rax,%rax
	je	L_skip6
	mulq	%rdi
	addq	%rax,32(%rsp)
	adcq	%rdx,24(%rsp)
	adcq	%rcx,16(%rsp)

L_skip6:
	movq	16(%rsi),%rax
	orq	%rax,%rax
	je	L_skip7
	mulq	%rdi
	addq	%rax,24(%rsp)
	adcq	%rdx,16(%rsp)
	adcq	%rcx,8(%rsp)

L_skip7:
	movq	8(%rsi),%rax
	orq	%rax,%rax
	je	L_skip8
	mulq	%rdi
	addq	%rax,16(%rsp)
	adcq	%rdx,8(%rsp)
L_skip8:
	movq	%rsi,%r10
	call	mdnorm
	addq	$80,%rsp
	ret
	.align	2












































	.globl	___mulm
___mulm:
	xorpd	%xmm0,%xmm0
	subq	$88,%rsp
	movupd	%xmm0,(%rsp)

	movupd	%xmm0,16(%rsp)
	movupd	%xmm0,32(%rsp)
	movupd	%xmm0,48(%rsp)
	movupd	%xmm0,64(%rsp)
	movq	$0,80(%rsp)

	movq	(%rsi),%rax
	movq	%rax,(%rsp)

	leaq	80(%rsp),%r11

	leaq	64(%rdi),%r9

	movd	%rbx,%xmm2
	leaq	32(%rsi),%r10

	movl	$7,%ebx
L_mulmFirstLoop:
	movq	(%r9),%rax
	movq	(%r10),%rdx
	mulq	%rdx
	xorq	%rcx,%rcx
	addq	%rax,(%r11)
	adcq	%rdx,-8(%r11)
	adcq	%rcx,-16(%r11)

	subq	$8,%r9
	addq	$8,%r10
	decl	%ebx
	cmpl	$3,%ebx
	jge	L_mulmFirstLoop

	subq	$8,%r11

	movl	$8,%r8d
L_mulmSecondMainLoop:
	movl	%r8d,%ebx
	leaq	(%rdi,%rbx,8),%r9
	leaq	16(%rsi),%r10
	jmp	L_mulmTestInnerLoop
	.align	2
L_mulmInnerLoop:
	movq	(%r9),%rax
	movq	(%r10),%rdx
	mulq	%rdx
	subq	$8,%r9
	xorq	%rcx,%rcx
	addq	%rax,(%r11)
	adcq	%rdx,-8(%r11)
	adcq	%rcx,-16(%r11)


	addq	$8,%r10
	subl	$1,%ebx
L_mulmTestInnerLoop:
	cmpl	$2,%ebx
	jge	L_mulmInnerLoop
	addq	$-8,%r11
	subl	$1,%r8d
	cmpl	$2,%r8d
	jge	L_mulmSecondMainLoop

	movq	%rsi,%r10
	call	mdnorm
	movd	%xmm2,%rbx
	addq	$88,%rsp
	ret
	.align	2







squarev:
	xorpd	%xmm0,%xmm0
	movq	%rsi,%r10
	movq	$0,8(%rsi)
	movupd	%xmm0,16(%rsi)
	movupd	%xmm0,32(%rsi)
	movupd	%xmm0,48(%rsi)
	movupd	%xmm0,64(%rsi)

	leaq	24(%rsi,%rdx,8),%r9

	movl	%edx,%esi
	addl	$2,%esi
	.align	2,0x90
L_$squarevStart:

	leaq	16(%rdi),%r11

	leaq	(%rdi,%rsi,8),%r8
	jmp	L_squarevIncrementOuterLoop
	.align	2
L_$squarevOuterLoop:


	movq	(%r8),%rdx
	orq	%rdx,%rdx
	je	L_squarevSkipMul
	movq	(%r11),%rax
	orq	%rax,%rax
	je	L_squarevSkipMul

	mulq	%rdx
	xorl	%ecx,%ecx

	cmpq	%r11,%r8
	je	L_squarevSkipShift



	clc
	rclq	%rax
	rclq	%rdx
	rclq	%rcx
L_squarevSkipShift:
	subq	$8,%r8
	addq	%rax,(%r9)
	adcq	$0,%rdx
	adcq	$0,%rcx
	addq	%rdx,-8(%r9)
	adcq	%rcx,-16(%r9)
	addq	$8,%r11


L_squarevIncrementOuterLoop:
	cmpq	%r11,%r8
	jae	L_$squarevOuterLoop

	subq	$8,%r9
	subl	$1,%esi
L_$squarev64:
	cmpl	$1,%esi
	jge	L_$squarevStart

	clc
	rclq $1,72(%r10)
	rclq $1,64(%r10)
	rclq $1,56(%r10)
	rclq $1,48(%r10)
	rclq $1,40(%r10)
	rclq $1,32(%r10)
	rclq $1,24(%r10)
	rclq $1,16(%r10)
	rclq $1,8(%r10)
	ret
	.align	2
L_squarevSkipMul:
	subq	$8,%r8
	addq	$8,%r11
	jmp	L_squarevIncrementOuterLoop
	.align	2













































mulv:
	movupd	%xmm0,16(%rdx)
	movupd	%xmm0,32(%rdx)
	movupd	%xmm0,48(%rdx)
	movupd	%xmm0,64(%rdx)
	movq	$0,8(%rdx)

	leaq	24(%rdx,%rbx,8),%r11

	addl	$2,%ebx
	.align	2,0x90
L_mulvOuterLoop:

	leaq	16(%rsi),%r8

	leaq	(%rdi,%rbx,8),%r9

	movl	%ebx,%r10d
	subl	$2,%r10d
	.align	2




L_mulvInnerLoop:

	movq	(%r9),%rax
	movq	(%r8),%rdx

	mulq	%rdx

	xorl	%ecx,%ecx
	subq	$8,%r9
	addq	%rax,(%r11)
	adcq	$0,%rdx
	adcl	$0,%ecx
	addq	%rdx,-8(%r11)
	adcq	%rcx,-16(%r11)



	addq	$8,%r8




	subl	$1,%r10d
	jge	L_mulvInnerLoop

	subq	$8,%r11
	subl	$1,%ebx
	cmpl	$2,%ebx
	jge	L_mulvOuterLoop
	ret
	.align	2
	.globl	___divi





___divi:
	subq	$96,%rsp



	movupd	(%rsi),%xmm0
	movupd	%xmm0,(%rsp)


	movq	%rdi,%r9

	movq	%rsi,%r10




	movq	%r10,%rdi
	movl	$2,%esi
	call	___shiftdownn



	movq	8(%r9),%rcx



	movq	16(%r10),%rdx
	movq	24(%r10),%rax



	divq	%rcx

	movq	%rax,16(%rsp)



	movq	32(%r10),%rax

	divq	%rcx
	movq	%rax,24(%rsp)
	movq	40(%r10),%rax

	divq	%rcx
	movq	%rax,32(%rsp)
	movq	48(%r10),%rax

	divq	%rcx
	movq	%rax,40(%rsp)
	movq	56(%r10),%rax

	divq	%rcx
	movq	%rax,48(%rsp)
	movq	64(%r10),%rax

	divq	%rcx
	movq	%rax,56(%rsp)
	movq	72(%r10),%rax

	divq	%rcx
	movq	%rax,64(%rsp)
	xorq	%rax,%rax



	divq	%rcx
	movq	%rax,72(%rsp)




	call	mdnorm
	addq	$96,%rsp
	ret
	.align	2



































	.globl	___divm
___divm:




	subq	$336,%rsp

	movq	%r12,312(%rsp)
	movq	%r13,320(%rsp)
	movq	%rbx,328(%rsp)
	xorpd	%xmm0,%xmm0
	movq	%rdi,%r12
	movq	%rsi,%r13

	movq	%rsp,%rdi
	xorl	%eax,%eax
	movl	$39,%ecx
	rep
	stosq

	movq	16(%r12),%rcx
	xorl	%eax,%eax
	movabsq $0x4000000000000000,%rdx
	divq	%rcx
	movq	%rax,120(%rsp)



	movq	%rax,%rcx
	mulq	%rcx
	xorq	%rcx,%rcx
	clc
	rclq	%rax
	rclq	%rdx
	rclq	%rcx
	movq	%rax,24+208(%rsp)
	movq	%rdx,16+208(%rsp)
	movq	%rcx,8+208(%rsp)


	movl	$2,%ebx
	movq	%rsp,%rdx
	lea    208(%rsp),%rsi
	movq	%r12,%rdi
	call	mulv

	lea	104(%rsp),%rdi

	movq 72(%rdi),%rax
	subq 72(%rsp),%rax
	movq %rax,72(%rdi)

	movq 64(%rdi),%rax
	sbbq 64(%rsp),%rax
	movq %rax,64(%rdi)

	movq 56(%rdi),%rax
	sbbq 56(%rsp),%rax
	movq %rax,56(%rdi)

	movq 48(%rdi),%rax
	sbbq 48(%rsp),%rax
	movq %rax,48(%rdi)

	movq 40(%rdi),%rax
	sbbq 40(%rsp),%rax
	movq %rax,40(%rdi)

	movq 32(%rdi),%rax
	sbbq 32(%rsp),%rax
	movq %rax,32(%rdi)

	movq 24(%rdi),%rax
	sbbq 24(%rsp),%rax
	movq %rax,24(%rdi)

	movq 16(%rdi),%rax
	sbbq 16(%rsp),%rax
	movq %rax,16(%rdi)

	movq 8(%rdi),%rax
	sbbq 8(%rsp),%rax
	movq %rax,8(%rdi)



	clc
	rclq $1,72(%rdi)
	rclq $1,64(%rdi)
	rclq $1,56(%rdi)
	rclq $1,48(%rdi)
	rclq $1,40(%rdi)
	rclq $1,32(%rdi)
	rclq $1,24(%rdi)
	rclq $1,16(%rdi)
	rclq $1,8(%rdi)


	movl	$4,%edx
	lea     208(%rsp),%rsi

	call	squarev

	movl	$4,%ebx
	movq	%rsp,%rdx
	lea    208(%rsp),%rsi
	movq	%r12,%rdi
	call	mulv


	lea	104(%rsp),%rdi

	movq 72(%rdi),%rax
	subq 72(%rsp),%rax
	movq %rax,72(%rdi)

	movq 64(%rdi),%rax
	sbbq 64(%rsp),%rax
	movq %rax,64(%rdi)

	movq 56(%rdi),%rax
	sbbq 56(%rsp),%rax
	movq %rax,56(%rdi)

	movq 48(%rdi),%rax
	sbbq 48(%rsp),%rax
	movq %rax,48(%rdi)

	movq 40(%rdi),%rax
	sbbq 40(%rsp),%rax
	movq %rax,40(%rdi)

	movq 32(%rdi),%rax
	sbbq 32(%rsp),%rax
	movq %rax,32(%rdi)

	movq 24(%rdi),%rax
	sbbq 24(%rsp),%rax
	movq %rax,24(%rdi)

	movq 16(%rdi),%rax
	sbbq 16(%rsp),%rax
	movq %rax,16(%rdi)

	movq 8(%rdi),%rax
	sbbq 8(%rsp),%rax
	movq %rax,8(%rdi)



	clc
	rclq $1,72(%rdi)
	rclq $1,64(%rdi)
	rclq $1,56(%rdi)
	rclq $1,48(%rdi)
	rclq $1,40(%rdi)
	rclq $1,32(%rdi)
	rclq $1,24(%rdi)
	rclq $1,16(%rdi)
	rclq $1,8(%rdi)

	movl	$7,%edx
	lea     208(%rsp),%rsi
	lea     104(%rsp),%rdi
	call	squarev

	movl	$7,%ebx
	movq	%rsp,%rdx
	lea    208(%rsp),%rsi
	movq	%r12,%rdi
	call	mulv

	lea	104(%rsp),%rsi
	movq	%rsp,%rdi
	call	___subm

	leaq	104(%rsp),%rdi
	clc
	rclq $1,72(%rdi)
	rclq $1,64(%rdi)
	rclq $1,56(%rdi)
	rclq $1,48(%rdi)
	rclq $1,40(%rdi)
	rclq $1,32(%rdi)
	rclq $1,24(%rdi)
	rclq $1,16(%rdi)
	rclq $1,8(%rdi)

	mov     $7,%rbx
	movq	%rsp,%rdx
	movq	%r13,%rsi

	call	mulv

	movq	(%r13),%rax
	movq	%rax,(%rsp)


	movq	%r13,%r10
	call	mdnorm


	movq	312(%rsp),%r12
	movq	320(%rsp),%r13
	movq	328(%rsp),%rbx
	addq	$336,%rsp

	ret
	.align	2



	.globl	_qequal
_qequal:
	movq	(%rdi),%rax
	cmpq	(%rsi),%rax
	jne	L_retzero
	movq	8(%rdi),%rax
	cmpq	8(%rsi),%rax
	jne	L_retzero
	movq	16(%rdi),%rax
	cmpq	16(%rsi),%rax
	jne	L_retzero
	movq	24(%rdi),%rax
	cmpq	24(%rsi),%rax
	jne	L_retzero
	movq	32(%rdi),%rax
	cmpq	32(%rsi),%rax
	jne	L_retzero
	movq	40(%rdi),%rax
	cmpq	40(%rsi),%rax
	jne	L_retzero
	movq	48(%rdi),%rax
	cmpq	48(%rsi),%rax
	jne	L_retzero
	movq	56(%rdi),%rax
	cmpq	56(%rsi),%rax
	jne	L_retzero
	movq	$1,%rax
	ret
	.align	2
L_retzero:
	xorq	%rax,%rax
	ret
	.align	2






mdnorm:
	leaq	8(%rsp),%r11
	movabs	$0x8000000000000000,%r8

	orq	$0,8(%r11)
	je	L_doRounding
	bsrq	8(%r11),%rsi
	addq	$1,%rsi
	movq	%r11,%rdi
	call	___shiftdownn


	addl	%ecx,4(%r11)
	cmpl	$0x100000,4(%r11)
	ja	L_mdnormOverflow


	testq	%r8,16(%r11)
	jne	L_doRounding

	movq	16(%r11),%rdx
	bsrq	%rdx,%rcx
	movl	4(%r11),%eax
	movl	$64,%esi
	subq	%rcx,%rsi
	subq	%rsi,%rax
	jl	L_doRounding
	movl	%eax,4(%r11)
	movq	%r11,%rdi
	subl	%esi,4(%r11)
	call	___shiftupn

L_doRounding:

	testq	%r8,72(%r11)
	je	L_shiftdown



	xorl	%eax,%eax
	addq	$1,64(%r11)
	adcq	%rax,56(%r11)
	adcq	%rax,48(%r11)
	adcq	%rax,40(%r11)
	adcq	%rax,32(%r11)
	adcq	%rax,24(%r11)
	adcq	%rax,16(%r11)
	adcq	%rax,8(%r11)

L_shiftdown:
	orq	$0,8(%r11)
	je	L_mdnormExit
	bsrq	8(%r11),%rcx
	addq	$1,%rcx
	addl	%ecx,4(%r11)
	cmpl	$0x100000,4(%r11)
	jae	L_mdnormOverflow
	movq	%r11,%rdi
	movq	%rcx,%rsi
	call	___shiftdownn



L_mdnormExit:
	movupd	(%r11),%xmm0
	movupd	%xmm0,(%r10)
	movupd	16(%r11),%xmm0
	movupd	%xmm0,16(%r10)
	movupd	32(%r11),%xmm0
	movupd	%xmm0,32(%r10)
	movupd	48(%r11),%xmm0
	movupd	%xmm0,48(%r10)
	movq	64(%r11),%rax
	movq	%rax,64(%r10)

	movq	$0,72(%r10)
	ret

L_mdnormOverflow:
	movl	$0x100000,4(%r11)
	jmp	L_mdnormExit
	.align	2














	.globl	_cmpm
_cmpm:


	movq	8(%rsi),%rax
	cmpq	%rax,8(%rdi)
	jne	L_Different

	movq	16(%rsi),%rax
	cmpq	%rax,16(%rdi)
	jne	L_Different

	movq	24(%rsi),%rax
	cmpq	%rax,24(%rdi)
	jne	L_Different

	movq	32(%rsi),%rax
	cmpq	%rax,32(%rdi)
	jne	L_Different

	movq	40(%rsi),%rax
	cmpq	%rax,40(%rdi)
	jne	L_Different

	movq	48(%rsi),%rax
	cmpq	%rax,48(%rdi)
	jne	L_Different

	movq	56(%rsi),%rax
	cmpq	%rax,56(%rdi)
	jne	L_Different

	movq	64(%rsi),%rax
	cmpq	%rax,64(%rdi)
	jne	L_Different

	movq	72(%rsi),%rax
	cmpq	%rax,72(%rdi)
	jne	L_Different


	xorl	%eax,%eax
	ret
	.align	2
L_Different:

	jb	L_cmpmRetMinus1
	movq	$1,%rax
	ret
L_cmpmRetMinus1:

	movq	$-1,%rax
	ret
	.align	2
	.globl	___addbit
___addbit:
	xorl	%ecx,%ecx
	addq	$1,64(%rdi)
	adcq	%rcx,56(%rdi)
	adcq	%rcx,48(%rdi)
	adcq	%rcx,40(%rdi)
	adcq	%rcx,32(%rdi)
	adcq	%rcx,24(%rdi)
	adcq	%rcx,16(%rdi)
	adcq	%rcx,8(%rdi)
	ret
