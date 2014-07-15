global start
extern _printf
extern _exit

SECTION .data
s0:     db "CPU INFO:",10
l0:     equ $-s0
s1:     db "Lagest number: %i",10,0
s2:     db "Vendor: %.12s",10,0
test:   db "Dies ist ein test",10,0
argc:   db "Number of argc: %i",10,0

SECTION .text

start:
	mov rax,0x2000004 		; write
	mov rdi,1			; stdout
	mov rsi,s0			; string
	mov rdx,l0			; length
	syscall

	mov rdi,argc
	mov rsi,rax
	call _printf


	push rbp			; save base stack pointer
	push rbx			; send to stack
	mov rbp,rsp			; save stack pointer
	sub rsp,16			; calc position on stack

	mov eax,0x0			; vendor 
	cpuid				; request

	mov [rsp],ebx
	mov [rsp+4],edx
	mov [rsp+8],ecx

	mov rdi,s1
	mov esi,eax
	call _printf

	mov rdi,s2
	mov rsi,rsp
	call _printf

	mov rdi,test
	call _printf

	mov rsp,rbp			; get correct stack pointer
	pop rbx				; reload variables
	pop rbp

	mov rsi,0x0			; call exit to flush io
	call _exit

	mov rax,0x2000001		; exit syscall (x86_64)
	mov rdi,0			; status = 0 (exit normally)
	syscall
