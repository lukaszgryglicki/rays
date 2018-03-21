.globl RGB24
.globl RGB16
.globl RR16
.globl RG16
.globl RB16
.globl RR24
.globl RG24
.globl RB24
.globl asm_flush_zbuff
.text
RGB24:
	pushl %edx	
	pushl %ecx
	movl  0x0C(%esp), %eax
	movl  0x10(%esp), %edx
	movl  0x14(%esp), %ecx
	#andl $0xFF, %eax
	#andl $0xFF, %edx
	#andl $0xFF, %ecx
	sall $0x10, %eax
	sall $0x08, %edx
	addl %edx, %eax	
	addl %ecx, %eax	
	popl %ecx
	popl %edx
	ret	
RGB16:
	pushl %edx	
	pushl %ecx	
	movl  0x0C(%esp), %eax
	movl  0x10(%esp), %edx
	movl  0x14(%esp), %ecx
	sarl $0x03, %eax
	sarl $0x02, %edx
	sarl $0x03, %ecx
	sall $0x0b, %eax
	sall $0x05, %edx
	addl %edx, %eax
	addl %ecx, %eax
	popl %ecx
	popl %edx
	ret
RR16:
	movl 0x04(%esp), %eax
	andl $0xF800, %eax
	sarl $0x0B, %eax
	sall $0x03, %eax
	ret
RG16:
	movl 0x04(%esp), %eax
	andl $0x07E0, %eax
	sarl $0x05, %eax
	sall $0x02, %eax
	ret
RB16:
	movl 0x04(%esp), %eax
	andl $0x1F, %eax
	sall $0x03, %eax
	ret
RR24:
	movl 0x04(%esp), %eax
	andl $0xFF0000, %eax
	sarl $0x10, %eax
	ret
RG24:
	movl 0x04(%esp), %eax
	andl $0x00FF00, %eax
	sarl $0x08, %eax
	ret
RB24:
	movl 0x04(%esp), %eax
	andl $0x0000FF, %eax
	ret
asm_flush_zbuff:
	pushl %ebx
	pushl %edx
	pushl %ecx
	
	xorl %edx, %edx
_label_i:
	xorl %ecx, %ecx
_label_j:
	movl %edx, %ebx
	sall $0x2, %ebx
	addl zbuff, %ebx
	movl (%ebx), %eax
	movl %ecx, %ebx		
	sall $0x2, %ebx
	addl %eax, %ebx
	cmpl $0x6FFFFFFF, (%ebx)
	jz _skip_point
	
	pushl %edx
	pushl %ecx
	pushl %ebx
	pushl (%ebx)
	pushl gc
	pushl dsp
	call XSetForeground
	addl $0x0C, %esp
	popl %ebx
	popl %ecx
	popl %edx

	pushl %edx
	pushl %ecx
	pushl %ebx
	pushl %ecx
	pushl %edx
	pushl gc
	pushl win
	pushl dsp
	call XDrawPoint
	addl $0x14, %esp
	popl %ebx
	popl %ecx
	popl %edx

_skip_point:
	movl $0x6FFFFFFF, (%ebx)
	incl %ecx
	cmpl $0x200, %ecx
	jl _label_j
	incl %edx
	cmpl $0x200, %edx
	jl _label_i
	
	popl %ecx
	popl %edx
	popl %ebx
	ret

