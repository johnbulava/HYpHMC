/**  void xSSE_PerformNopLoop(char* parameters);
*    -->Parameters:                           Offset
*          long int count                     0
*          long int error                     8
*
*    -->Target: 64-Bit AMD64 architecture with GNU g++
*
*    -->Workspace Demand: 0
*
*    -->Description:
*          This routine performs count iterations of nop.
*          Error is always set to 0.
*
**/
.file	"xSSE_PerformNopLoop.S"
.text
.align 2
.globl xSSE_PerformNopLoop
.type	xSSE_PerformNopLoop, @function
xSSE_PerformNopLoop:
  #Guarantee all memory accesses be complete
    mfence
    
  #Save all registers being used in this routine
    push %rbp
    push %rax

  #Actual Code of Routine    
    mov %rdi, %rbp
    mov 0(%rbp), %rax
    MainLoop:
      nop
      dec %rax
    jg MainLoop
    	    

  #Set Error to Zero 
    mov $0, %rax
    mov %rax, 8(%rbp)
    
  #Restoring saved registers
  ClearStack:
    pop %rax
    pop %rbp
    
  #Leaving the routine and guarantee all memory accesses be complete
    mfence
    ret
