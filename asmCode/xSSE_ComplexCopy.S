.file	"xSSE_ComplexCopy.S"
.data
ausgabe: .ascii "The routine xSSE_ComplexCopy has not been implemented for 32-bit yet!\n"
ausgabeEnde: .ascii "EXIT\n"
.text
.align 2
.globl xSSE_ComplexCopy
.type	xSSE_ComplexCopy, @function
xSSE_ComplexCopy:

  #Setting up Stack-Frame
    push %ebp
    mov %esp, %ebp
  
  #Save all registers being used in this routine
    push %eax
    push %ebx
    push %ecx
    push %edx
    push %esi
    push %edi


  #Actual Assember-Code
    mov $4, %eax
    mov $1, %ebx
    mov $ausgabe, %ecx
    mov $ausgabeEnde, %edx
    sub %ecx, %edx
    add $5, %edx          //Laenge des Ausgabe-Strings
    int $0x80
    
    mov $1, %eax
    mov $0, %ebx
    int $0x80
    
  #Restoring saved registers
    pop %edi
    pop %esi
    pop %edx
    pop %ecx
    pop %ebx
    pop %eax
    
  #Desinitializing Stack-Frame
    pop %ebp
  
  #Leaving the routine
    ret
