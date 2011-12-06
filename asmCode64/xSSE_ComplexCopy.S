/**  void xSSE_ComplexCopy(char* parameters);
*    -->Parameters:                           Offset
*          char* workSpace                    0
*          Complex* input                     8
*          Complex* output                    16
*          long int blockSize                 24
*          long int L0                        32
*          long int L1                        40
*          long int L2                        48
*          long int L3                        56
*          long int xtrAdd1                   64
*          long int xtrAdd2                   72
*          long int xtrAdd3                   80
*          long int L3LineCombineCount        88
*          long int PrefetchDelayCycles       96
*          long int error                     104
*          long int TEMPORAL                  112
*          long int TEMPORAL                  120
*          long int TEMPORAL                  128
*          long int TEMPORAL                  136
*
*
*    -->Target: 64-Bit AMD64 architecture with GNU g++
*
*    -->Workspace Demand: 
*
*    -->Description:
*          This routine copies complex numbers (double, double)
*          given in a lattice-structure (L0,L1,L2,L3 ) from input 
*          to output. There are blockSize complex numbers per
*          lattice site. XtrAdd? gives additional increments
*          according to a possible embedded lattice-structure.
*          L3LineCombineCount is the number of L3-Lattice lines
*          that are read in one stride. PrefetchDelayCycles is
*          the number of cycles waiting after one prefetch.
*          Error is set to 0 on success or to an error code if
*          the routine failed.
**/

.file	"xSSE_ComplexCopy.S"
.text
.align 2
.globl xSSE_ComplexCopy
.type	xSSE_ComplexCopy, @function
xSSE_ComplexCopy:
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
    mov 8(%rbp), %rsi      # Source
    mov 16(%rbp), %rdi     # Destination        
    mov 0(%rbp), %r15      # Work-Space

    //Preparing Constants for Main loop 
    mov 88(%rbp), %rbx     # L3LineCombineCount
    cmp $0, %rbx
    jle Error2
    mov 48(%rbp), %rax 
    xor %rdx, %rdx
    div %rbx
    mov %rax, 112(%rbp)    # Number of L3-Line-Blocks
    mov %rdx, 120(%rbp)    # Number of Extra L3-Lines    
    mov 56(%rbp), %r8      # L3
    mov 24(%rbp), %r9      # BlockSize
    imul %r9, %r8
    mov %r8, %r9
    shr $2, %r8            # Inner-Block-Loop Count
    jle Error1
    and $3, %r9           
    mov %r9, 128(%rbp)     # Inner-Rest-Loop Count (0 if no rest loop)
    shl $4, %r9
    mov %r9, 136(%rbp)     # Inner-Rest-Loop Increment for pointer (0 if no rest loop)
    mov 96(%rbp), %r10      # PrefetchDelayCycles
    

    //Main-Loop    
    mov 32(%rbp), %rax      # L0
    loopL0:
      mov 40(%rbp), %rbx      # L1
      loopL1:
        mov 112(%rbp), %rcx     # Number of L3-Line-Blocks
        loopL2CombineBlocks:

          //Reading Data
          mov %rsi, %r14          # Save Source - Pointer
          mov 88(%rbp), %r9       # L3LineCombineCount
          loopLineCombineReadA:
            mov %r8, %rdx          # Inner-Block-Loop Count
            InnerBlockLoopReadA:
  	      prefetchnta 0(%rsi)
	    
 	      mov %r10, %r11
	      DelayLoopA:
	        dec %r11
	      jg DelayLoopA
	    
              add $64, %rsi
  	      dec %rdx
            jg InnerBlockLoopReadA

            cmp $0, 136(%rbp)      # Inner-Rest-Loop 
	    jle NoRestLoopReadA
	      prefetchnta 0(%rsi)
              add 136(%rbp), %rsi
            NoRestLoopReadA:
	    
  	    add 80(%rbp), %rsi     # Extra-Add for Source Pointer
	    dec %r9
	  jg loopLineCombineReadA


          //Writing Data
	  mov %r14, %rsi          # Restoring Source-Pointer
          mov 88(%rbp), %r9       # L3LineCombineCount
          loopLineCombineWriteA:
            mov %r8, %rdx          # Inner-Block-Loop Count	  
            InnerBlockLoopWriteA:
	      mov 0(%rsi), %r11
  	      mov 8(%rsi), %r12
	      mov 16(%rsi), %r13
	      mov 24(%rsi), %r14

	      movnti %r11, 0(%rdi)
	      movnti %r12, 8(%rdi)
	      movnti %r13, 16(%rdi)
	      movnti %r14, 24(%rdi)
	  
	      mov 32(%rsi), %r11
	      mov 40(%rsi), %r12
	      mov 48(%rsi), %r13
	      mov 56(%rsi), %r14

	      movnti %r11, 32(%rdi)
	      movnti %r12, 40(%rdi)
	      movnti %r13, 48(%rdi)
	      movnti %r14, 56(%rdi)
	    
              add $64, %rsi
              add $64, %rdi
  	      dec %rdx
            jg InnerBlockLoopWriteA
	    
            mov 128(%rbp), %rdx
	    cmp $0, %rdx
	    jle NoRestLoopWriteA	    
            RestLoopWriteA:	   
	      mov 0(%rsi), %r11
  	      mov 8(%rsi), %r12

	      movnti %r11, 0(%rdi)
	      movnti %r12, 8(%rdi)
	      
              add $16, %rsi
              add $16, %rdi
	      dec %rdx
            jg RestLoopWriteA
	    NoRestLoopWriteA: 

  	    add 80(%rbp), %rsi      # Extra-Add for Source Pointer
	    add 80(%rbp), %rdi      # Extra-Add for Destination Pointer	  
            dec %r9
          jg loopLineCombineWriteA
	  
	  dec %rcx
	jg loopL2CombineBlocks	  

	
        //Extra L3 - Lines
	cmp $0, 120(%rbp)
	jle NoExtraL3Lines	
          //Reading Data
          mov %rsi, %r14           # Save Source - Pointer
   	  mov 120(%rbp), %r9       # Number of Extra L3-Lines
          loopLineCombineReadB:
            mov %r8, %rdx          # Inner-Block-Loop Count
            InnerBlockLoopReadB:
  	      prefetchnta 0(%rsi)
	    
 	      mov %r10, %r11
	      DelayLoopB:
	        dec %r11
	      jg DelayLoopB
	    
              add $64, %rsi
  	      dec %rdx
            jg InnerBlockLoopReadB

            cmp $0, 136(%rbp)      # Inner-Rest-Loop 
	    jle NoRestLoopReadB
	      prefetchnta 0(%rsi)
              add 136(%rbp), %rsi
            NoRestLoopReadB:
	    
  	    add 80(%rbp), %rsi     # Extra-Add for Source Pointer
	    dec %r9
	  jg loopLineCombineReadB

          //Writing Data
	  mov %r14, %rsi           # Restoring Source-Pointer
          mov 120(%rbp), %r9       # Number of Extra L3-Lines
          loopLineCombineWriteB:
            mov %r8, %rdx          # Inner-Block-Loop Count	  
            InnerBlockLoopWriteB:
	      mov 0(%rsi), %r11
  	      mov 8(%rsi), %r12
	      mov 16(%rsi), %r13
	      mov 24(%rsi), %r14

	      movnti %r11, 0(%rdi)
	      movnti %r12, 8(%rdi)
	      movnti %r13, 16(%rdi)
	      movnti %r14, 24(%rdi)
	  
	      mov 32(%rsi), %r11
	      mov 40(%rsi), %r12
	      mov 48(%rsi), %r13
	      mov 56(%rsi), %r14

	      movnti %r11, 32(%rdi)
	      movnti %r12, 40(%rdi)
	      movnti %r13, 48(%rdi)
	      movnti %r14, 56(%rdi)
	    
              add $64, %rsi
              add $64, %rdi
  	      dec %rdx
            jg InnerBlockLoopWriteB
	    
            mov 128(%rbp), %rdx
	    cmp $0, %rdx
	    jle NoRestLoopWriteB    
            RestLoopWriteB:	   
	      mov 0(%rsi), %r11
  	      mov 8(%rsi), %r12

	      movnti %r11, 0(%rdi)
	      movnti %r12, 8(%rdi)
	      
              add $16, %rsi
              add $16, %rdi
	      dec %rdx
            jg RestLoopWriteB
	    NoRestLoopWriteB: 

  	    add 80(%rbp), %rsi      # Extra-Add for Source Pointer
	    add 80(%rbp), %rdi      # Extra-Add for Destination Pointer	  
            dec %r9
          jg loopLineCombineWriteB
	NoExtraL3Lines:
	  
	add 72(%rbp), %rsi
	add 72(%rbp), %rdi	  	
        dec %rbx
      jg loopL1
      add 64(%rbp), %rsi
      add 64(%rbp), %rdi	  
      dec %rax
    jg loopL0
	    

  #Set Error to Zero 
    mov $0, %rax
    mov %rax, 104(%rbp)
    
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
    //Set Error to 1         # L3*blockSize < 4
    mov $1, %rax
    mov %rax, 104(%rbp)
    jmp ClearStack
  Error2:                       # L3LineCombineCount <= 0
    //Set Error to 2 
    mov $2, %rax
    mov %rax, 104(%rbp)
    jmp ClearStack
    
