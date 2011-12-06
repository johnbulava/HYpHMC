/**  void xSSE_ReadLongIntFromMemAddr(char* parameters);
*    -->Parameters:                           Offset
*          long int* memAddr                  0
*          long int value                     8
*          long int error                     16
*
*    -->Target: 64-Bit AMD64 architecture with GNU g++
*
*    -->Workspace Demand: 0
*
*    -->Description:
*          This routine reads the long integer value from the 
*          memory address memAddr. Error is always set to 0.
*
**/
.file	"xSSE_ReadLongIntFromMemAddr.S"
.text
.align 2
.globl xSSE_ReadLongIntFromMemAddr
.type	xSSE_ReadLongIntFromMemAddr, @function
xSSE_ReadLongIntFromMemAddr:
  #Guarantee all memory accesses be complete
    mfence
    
  #Save all registers being used in this routine
    push %rbp
    push %rax
    push %rbx

  #Actual Code of Routine    
    mov %rdi, %rbp
    mov 0(%rbp), %rax
    mov (%rax), %rbx    
    mov %rbx, 8(%rbp)
    	    

  #Set Error to Zero 
    mov $0, %rax
    mov %rax, 16(%rbp)
    
  #Restoring saved registers
  ClearStack:
    pop %rbx
    pop %rax
    pop %rbp
    
  #Leaving the routine and guarantee all memory accesses be complete
    mfence
    ret
