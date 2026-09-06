/* Linux x86-64 syscall-only helper entry. Bindings are appended by the owner.
 * No libc, interpreter, dynamic tags, environment additions or shell parsing.
 * Exit 125: invalid argc/role. Exit 126: execve failed. */
.text
.global _start
.type _start,@function
_start:
    mov %rsp, %r12
    mov (%r12), %r13
    cmp $1, %r13
    jb .Linvalid
    cmp $512, %r13
    ja .Linvalid
    mov 8(%r12), %rbx
    test %rbx, %rbx
    jz .Linvalid
    mov %rbx, %r8
    mov %rbx, %rsi
    mov $4096, %ecx
.Lbasename:
    movzbl (%rsi), %eax
    test %al, %al
    jz .Lselect
    cmp $47, %al
    jne .Ladvance
    lea 1(%rsi), %r8
.Ladvance:
    inc %rsi
    dec %ecx
    jnz .Lbasename
    jmp .Linvalid
.Lselect:
    lea helper_table(%rip), %r9
    mov $4, %ecx
.Lrole:
    mov %r8, %rsi
    mov (%r9), %rdi
.Lcompare_role:
    movzbl (%rsi), %eax
    cmp (%rdi), %al
    jne .Lnext_role
    test %al, %al
    jz .Lmatched
    inc %rsi
    inc %rdi
    jmp .Lcompare_role
.Lnext_role:
    add $16, %r9
    dec %ecx
    jnz .Lrole
    jmp .Linvalid
.Lmatched:
    mov 8(%r9), %r14
    lea 16(%r12,%r13,8), %rdx /* original envp */
    lea 72(,%r13,8), %rax   /* argc + 9 pointer slots, including NULL */
    sub %rax, %rsp
    and $-16, %rsp
    mov %rsp, %rsi
    lea loader_label(%rip), %rax
    mov %rax, 0(%rsi)
    lea no_cache(%rip), %rax
    mov %rax, 8(%rsi)
    lea no_hwcaps(%rip), %rax
    mov %rax, 16(%rsi)
    lea empty_string(%rip), %rax
    mov %rax, 24(%rsi)
    lea library_option(%rip), %rax
    mov %rax, 32(%rsi)
    lea library_prefix(%rip), %rax
    mov %rax, 40(%rsi)
    lea argv0_option(%rip), %rax
    mov %rax, 48(%rsi)
    mov %rbx, 56(%rsi)
    mov %r14, 64(%rsi)
    mov $1, %ecx
.Lcopy:
    cmp %r13, %rcx
    jae .Lexec
    mov 8(%r12,%rcx,8), %r15
    mov %r15, %r10
    lea plugin_original(%rip), %rdi
.Lcompare_plugin:
    movzbl (%rdi), %eax
    cmp (%r10), %al
    jne .Lstore
    test %al, %al
    jz .Lplugin
    inc %rdi
    inc %r10
    jmp .Lcompare_plugin
.Lplugin:
    lea plugin_descriptor(%rip), %r15
.Lstore:
    mov %r15, 64(%rsi,%rcx,8)
    inc %rcx
    jmp .Lcopy
.Lexec:
    movq $0, 64(%rsi,%r13,8)
    lea loader_descriptor(%rip), %rdi
    mov $59, %eax /* __NR_execve */
    syscall
    mov $126, %edi
    jmp .Lexit
.Linvalid:
    mov $125, %edi
.Lexit:
    mov $60, %eax /* __NR_exit */
    syscall
    ud2
.size _start, .-_start
.section .rodata
no_cache: .asciz "--inhibit-cache"
no_hwcaps: .asciz "--glibc-hwcaps-mask"
empty_string: .asciz ""
library_option: .asciz "--library-path"
argv0_option: .asciz "--argv0"
.section .note.GNU-stack,"",@progbits
