/**  void xSSE_WriteLongIntToMemAddr(char* parameters);
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
*          This routine writes the long integer value to the 
*          memory address memAddr. Error is always set to 0.
*
**/
.file	"xSSE_WriteLongIntToMemAddr.S"
.text
.align 2
.globl xSSE_WriteLongIntToMemAddr
.type	xSSE_WriteLongIntToMemAddr, @function
xSSE_WriteLongIntToMemAddr:
  #Guarantee all memory accesses be complete
    mfence
    
  #Save all registers being used in this routine
    push %rbp
    push %rax
    push %rbx

  #Actual Code of Routine    
    mov %rdi, %rbp
    mov 0(%rbp), %rax
    mov 8(%rbp), %rbx
    movnti %rbx, (%rax)
    	    

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
