.globl nearly_equal
.globl compute_angle
.data
epsilonm7:
	.long	-444973056
	.long	-692087595
	.long	16359
ld40:
	.long	0
	.long	-2147483648
	.long	16385
.text
nearly_equal:
	fldt 4(%esp)
	fldt 16(%esp)
	fsubp %st, %st(1)
	fabs
	fldt 28(%esp)
	fucompp
	fnstsw	%ax
	testb	$69, %ah
	sete	%al
	movzbl	%al, %eax
	ret
compute_angle:
	movl 4(%esp), %edx
	movl 8(%esp), %ebx
	movl 12(%esp), %ecx
	fldt (%edx)
	fldt (%ebx)
	fldt (%ecx)
	fsub %st(2)
	fxch %st(1)
	fsub %st(2)
	fxch %st(2)
	fstp %st
	fldt 12(%edx)
	fldt 12(%ebx)
	fldt 12(%ecx)
	fsub %st(2)
	fxch %st(1)
	fsub %st(2)
	fxch %st(2)
	fstp %st
	fldt 24(%edx)
	fldt 24(%ebx)
	fldt 24(%ecx)
	fsub %st(2)
	fxch %st(1)
	fsub %st(2)
	fxch %st(2)
	fstp %st
	fld %st
	fmul %st
	fld %st(3)
	fmul %st
	faddp %st(1)
	fld %st(5)
	fmul %st
	faddp %st(1)
	fsqrt
	fldt	epsilonm7
	fucomp
	fnstsw	%ax
	testb	$69, %ah
	je	compute_angle_sqrt1_fail
	jmp compute_angle_sqrt1
compute_angle_sqrt1_fail:
	fstp %st
	fstp %st
	fstp %st
	fstp %st
	fstp %st
	fstp %st
	fstp %st
	fldt	ld40
	ret
compute_angle_sqrt1:
	fxch %st(1)
	fdiv %st(1)
	fxch %st(3)
	fdiv %st(1)
	fxch %st(5)
	fdiv %st(1)
	fxch %st(1)
	fstp %st
	fld %st(1)
	fmul %st
	fld %st(4)
	fmul %st
	faddp %st(1)
	fld %st(6)
	fmul %st
	faddp %st(1)
	fsqrt
	fldt	epsilonm7
	fucomp
	fnstsw	%ax
	testb	$69, %ah
	je	compute_angle_sqrt1_fail
	jmp compute_angle_sqrt2
compute_angle_sqrt2:
	fxch %st(2)
	fdiv %st(2)
	fxch %st(4)
	fdiv %st(2)
	fxch %st(6)
	fdiv %st(2)
	fxch %st(2)
	fstp %st
	fmulp %st(1)
	fxch %st(1)
	fmulp %st(2)
	fxch %st(2)
	fmulp %st(3)
	faddp %st(1)
	faddp %st(1)
	ret
