/**  void xSSE_ExtremeFFT4D_Write2DSlice(char* parameters);
*    -->Parameters:                                         Offset       Remark
*          char* workSpace                                  0
*          Complex* slice                                   8
*          Complex* output                                  16
*          long int Size1                                   24
*          long int Size2                                   32
*          long int localIndexCount                         40
*          long int add1                                    48
*          long int add2                                    56
*          long int embeddingOne                            64           Every line xxx cache lines
*          long int embeddingTwo                            72           Every 4 lines xxx cache lines
*          long int error                                   80
*
*
*    -->Target: 64-Bit AMD64 architecture with GNU g++
*
*    -->Workspace Demand: 0
*
*    -->Description:
*          This routine writes the data in the 2D-slice-array
*          back to main memory (output).
*          Size1, size2 are the sizes of the 2D-slice.
*          LocalIndexCount is the number of nested lattices. EmbeddingOne
*          denotes number of inserted cache lines at each
*          end of line. EmbeddingTwo 
*          gives number of inserted dummy cache-lines in output-slide every 4 lines.
*          Error is set to 0 on success or to an error code if the 
*          routine failed.
**/
#include "xSSE_Global.h"
.file	"xSSE_ExtremeFFT4D_Write2DSlice.S"
.text
.align 2
.globl xSSE_ExtremeFFT4D_Write2DSlice
.type  xSSE_ExtremeFFT4D_Write2DSlice, @function
xSSE_ExtremeFFT4D_Write2DSlice:
  #Guarantee all memory accesses be complete
    mfence
    
  #Save all registers being used in this routine
    push %rbp
    push %rax
    push %rbx
    push %rcx
    push %rdx
    push %rsi
    push %rdi
    push %r8
    push %r9
    push %r10
    push %r11
    push %r12
    push %r13
    push %r14
    push %r15

  #Getting Pointer to Parameters
    mov %rdi, %rbp

  #Actual Code of Routine    
    //Load Slice and Output - Addresses
    mov 8(%rbp), %rsi
    mov 16(%rbp), %rdi
  
    //Preparing Constants for Main loop And Error-Check
    mov 24(%rbp), %r8     # Check Size1 > 0
    cmp $0, %r8
    jle Error1
    
    mov 32(%rbp), %r9     # Check Size2 > 0
    cmp $0, %r9
    jle Error2

    mov 40(%rbp), %r10    # Check localIndexCount > 0
    cmp $0, %r10
    jle Error3

    mov 48(%rbp), %r11     # add1   
    mov 56(%rbp), %r12     # add2   
    mov 64(%rbp), %rax
    shl $6, %rax
    mov %rax, 64(%rbp)     
    mov 72(%rbp), %r13
    shl $6, %r13           # 64*embeddingTwo
    
    
    //Main-Loop    
    mov $0, %rax
    MainLoop1:
//      mov %r9, %rbx
      mov 32(%rbp), %rbx

      MainLoop2:
//        mov %r10, %rcx
        mov 40(%rbp), %rcx

	InnerLoop1:
	  mov (%rsi), %rdx
	  mov 8(%rsi), %r15
	  mov 16(%rsi), %r9
	  mov 24(%rsi), %r10
	  xSSE_ExtremeFFT4D_Write2DSlice_SlideWriteCMD %rdx, (%rdi)
	  xSSE_ExtremeFFT4D_Write2DSlice_SlideWriteCMD %r15, 8(%rdi)


	  xSSE_ExtremeFFT4D_Write2DSlice_SlideWriteCMD %r9, 16(%rdi)
	  xSSE_ExtremeFFT4D_Write2DSlice_SlideWriteCMD %r10, 24(%rdi)


          add $32, %rsi
          add $32, %rdi
	
	  sub $2, %rcx
	jg InnerLoop1
	
	add %r12, %rdi
      
        dec %rbx
      jg MainLoop2
    
      add %r11, %rdi
      mov $0, %rdx
      inc %rax
      test $3, %rax
      cmove %r13, %rdx
      add %rdx, %rsi
      add 64(%rbp), %rsi            // 64*embeddingOne

      cmp %rax, %r8
    jg MainLoop1


  #Set Error to Zero 
    mov $0, %rax
    mov %rax, 80(%rbp)
    
  #Restoring saved registers
  ClearStack:
    pop %r15
    pop %r14
    pop %r13
    pop %r12
    pop %r11
    pop %r10
    pop %r9
    pop %r8
    pop %rdi
    pop %rsi
    pop %rdx
    pop %rcx
    pop %rbx
    pop %rax
    pop %rbp
    
  #Leaving the routine and guarantee all memory accesses be complete
    mfence
    ret
    
  #Error - Handling
  Error1:
    //Set Error to 1         # Size1 <= 0
    mov $1, %rax
    mov %rax, 80(%rbp)
    jmp ClearStack
  Error2:
    //Set Error to 2         # Size2 <= 0
    mov $2, %rax
    mov %rax, 80(%rbp)
    jmp ClearStack
  Error3:                    # localIndexCount <= 0
    //Set Error to 3 
    mov $3, %rax
    mov %rax, 80(%rbp)
    jmp ClearStack
