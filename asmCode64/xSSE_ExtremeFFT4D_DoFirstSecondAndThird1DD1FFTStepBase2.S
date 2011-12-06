/**  void xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2(char* parameters);
*    -->Parameters:                                         Offset       Remark
*          char* workSpace                                  0
*          Complex* slice                                   8
*          long int Size1                                   16
*          long int Size2                                   24
*          long int localIndexCount                         32
*          long int embeddingOne                            40           Every line xxx cache lines
*          long int embeddingTwo                            48           Every 4 lines xxx cache lines
*          Complex* FFTcomplexFacs1                         56           Pointer to complex numbers for thirf FFT in 1-direction
*          Complex facpmI (bitmask)                         64           Bitmask for multiplication with plus/minus I
*          int** prefetchAddressControlList                 80           Pointer to list of prefetchAddressControls
*          long int error                                   120
*
*
*          long int TEMPORAL                                128          Counter: Outer-Loop 
*          long int TEMPORAL                                136          Counter: Inner-Loop
*          long int TEMPORAL                                144          Size1 / 8
*          long int TEMPORAL                                152          Size2 * localIndexCount
*          long int TEMPORAL                                160          1*(Size2*localIndexCount*16 + 64*embeddingOne)
*          long int TEMPORAL                                168          2*(Size2*localIndexCount*16 + 64*embeddingOne)
*          long int TEMPORAL                                176          3*(Size2*localIndexCount*16 + 64*embeddingOne)
*          long int TEMPORAL                                184          4*(Size2*localIndexCount*16 + 64*embeddingOne) + 1*64*embeddingTwo
*          long int TEMPORAL                                192          8*(Size2*localIndexCount*16 + 64*embeddingOne) + 2*64*embeddingTwo - Size2*localIndexCount*16
*          
*
*
*    -->Target: 64-Bit AMD64 architecture with GNU g++
*
*    -->Workspace Demand: 0
*
*    -->Description:
*          This routine performs the first, second, and third FFT-step
*          in 1-direction with base power 2.
*          Input and Output is the 2D-slice-array.
*          Size1, size2 are the sizes of the 2D-slice and must be
*          multiples of 8. PrefetchAddressControl contains 
*          prefetch-address information for efficient prefetching
*          and must be provided by the calling routine. FFTcomplexFacs1
*          contains the complex-numbers for the third FFT-step in form (Re, Re), (-Im, Im).
*          LocalIndexCount is the number of nested lattices. EmbeddingOne
*          gives number of inserted cach line at end of each line. EmbeddingTwo 
*          gives number of inserted dummy cache-lines in output-slide every 4 lines.
*          Error is set to 0 on success or to an error code if the 
*          routine failed.
**/
#include "xSSE_Global.h"
.file	"xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2.S"
.text
.align 2
.globl xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2
.type  xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2, @function
xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2:
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
    //Preparing Constants for Main loop And Error-Check
    mov 16(%rbp), %rax     # Check Size1 > 0
    cmp $0, %rax
    jle Error1
    test $7, %rax          # Check Size1 multiple of 8
    jg Error2
    shr $3, %rax
    mov %rax, 144(%rbp)    # Size1 / 8
    
    mov 24(%rbp), %rax     # Check Size2 > 0
    cmp $0, %rax
    jle Error3

    mov 32(%rbp), %rcx     # Check localIndexCount > 0
    cmp $0, %rcx
    jle Error4
    imul %rcx, %rax   
    mov %rax, 152(%rbp)    # Size2 * localIndexCount
    
    
    mov 40(%rbp), %rbx     # embeddingOne
    shl $6, %rbx
    shl $4, %rax
    add %rbx, %rax
    mov %rax, 160(%rbp)     # 1*(Size2*localIndexCount*16 + 64*embeddingOne)   
    mov %rax, %rdx
    add %rdx, %rax
    mov %rax, 168(%rbp)     # 2*(Size2*localIndexCount*16 + 64*embeddingOne)   
    add %rdx, %rax
    mov %rax, 176(%rbp)     # 3*(Size2*localIndexCount*16 + 64*embeddingOne)   
    add %rdx, %rax
    mov 48(%rbp), %rcx     # embeddingTwo
    shl $6, %rcx
    add %rcx, %rax
    mov %rax, 184(%rbp)     # 4*(Size2*localIndexCount*16 + 64*embeddingOne) + 1*64*embeddingTwo 
    add %rax, %rax
    mov 152(%rbp), %rdx
    shl $4, %rdx
    sub %rdx, %rax
    mov %rax, 192(%rbp)     # 8*(Size2*localIndexCount*16 + 64*embeddingOne) + 2*64*embeddingTwo - Size2*localIndexCount*16


    //Load Slice-Data-Pointers
    mov 8(%rbp), %rsi              # Slice-Data
    mov 160(%rbp), %r8
    mov 168(%rbp), %r9
    mov 176(%rbp), %r10
    mov 184(%rbp), %r11
    add %rsi, %r11
    
    //Load constants
    movapd 64(%rbp), %xmm15	
    mov 56(%rbp), %rax
    movapd 0(%rax), %xmm13
    movapd 16(%rax), %xmm14    
    
    
    //Main-Loop    
    mov 144(%rbp), %rax
    mov %rax, 128(%rbp)            # Counter: Outer-Loop     
    OuterLoop:
      //Load Prefetch-Data
      #ifdef xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2_SlidePrefetchUSE
        mov 80(%rbp), %rax
        mov (%rax), %rdx
        add $8, %rax
        mov %rax, 80(%rbp)
      #endif
    
      mov 152(%rbp), %rax
      mov %rax, 136(%rbp)          # Counter: Inner-Loop 
      InnerLoop:
        movapd (%rsi), %xmm0
	movapd (%rsi, %r8), %xmm1
	movapd (%rsi, %r9), %xmm2
	movapd (%rsi, %r10), %xmm3	
	movapd (%r11), %xmm4
	movapd (%r11, %r8), %xmm5
	movapd (%r11, %r9), %xmm6
	movapd (%r11, %r10), %xmm7

        #ifdef xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2_SlidePrefetchUSE
          mov (%rdx), %rax
          movsx %eax, %rbx
          xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2_SlidePrefetchCMD (%rsi, %rbx)
          shr $32, %rax
        #endif

        //First FFT-step
        subpd %xmm1, %xmm0
        subpd %xmm3, %xmm2
        subpd %xmm5, %xmm4
        subpd %xmm7, %xmm6
	  
	addpd (%rsi), %xmm1
	addpd (%rsi, %r9), %xmm3
	addpd (%r11), %xmm5
	addpd (%r11, %r9), %xmm7
	  
        #ifdef xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2_SlidePrefetchUSE
          movsx %eax, %rbx
          xSSE_ExtremeFFT4D_DoFirstSecondAndThird1DD1FFTStepBase2_SlidePrefetchCMD (%rsi, %rbx)
        #endif
	  
	//Second FFT-step
	movhlps %xmm2, %xmm8
	movlhps %xmm2, %xmm8
	movhlps %xmm6, %xmm9
	movlhps %xmm6, %xmm9
	  
	xorpd %xmm15, %xmm8
	xorpd %xmm15, %xmm9
	  
	movapd %xmm0, %xmm2
	movapd %xmm1, %xmm6
	movapd %xmm4, %xmm10
	movapd %xmm5, %xmm11

	  
	subpd %xmm8, %xmm0
	subpd %xmm3, %xmm1
	subpd %xmm9, %xmm4
	subpd %xmm7, %xmm5
	  
	addpd %xmm2, %xmm8
	addpd %xmm6, %xmm3
	addpd %xmm10, %xmm9
	addpd %xmm11, %xmm7
	  
	//Third FFT-step
	movhlps %xmm4, %xmm6
	movlhps %xmm4, %xmm6
	movhlps %xmm5, %xmm10
	movlhps %xmm5, %xmm10
	movhlps %xmm9, %xmm11
	movlhps %xmm9, %xmm11
	  
        mulpd %xmm13, %xmm9
        mulpd %xmm14, %xmm11
	xorpd %xmm15, %xmm10
        mulpd %xmm13, %xmm4
        mulpd %xmm14, %xmm6

	addpd %xmm9, %xmm11
	subpd %xmm4, %xmm6
	  
	movapd %xmm0, %xmm2
	movapd %xmm1, %xmm4
	movapd %xmm8, %xmm5
	movapd %xmm3, %xmm9
	  
	subpd %xmm6, %xmm0
	subpd %xmm10, %xmm1
	subpd %xmm11, %xmm8
	subpd %xmm7, %xmm3
	  
	addpd %xmm2, %xmm6
	addpd %xmm4, %xmm10
	addpd %xmm5, %xmm11
	addpd %xmm9, %xmm7
	  
	movapd %xmm7, (%rsi)
	movapd %xmm11, (%rsi, %r8)
	movapd %xmm10, (%rsi, %r9)
	movapd %xmm6, (%rsi, %r10)
	movapd %xmm3, (%r11)
	movapd %xmm8, (%r11, %r8)
	movapd %xmm1, (%r11, %r9)
	movapd %xmm0, (%r11, %r10)
	      
        add $16, %rsi
	add $16, %r11
	add $8, %rdx
        subq $1, 136(%rbp)
      jg InnerLoop
      
      //Slice-Pointer-Increment
      add 192(%rbp), %rsi
      add 192(%rbp), %r11
	
      subq $1, 128(%rbp)
    jg OuterLoop


  #Set Error to Zero 
    mov $0, %rax
    mov %rax, 120(%rbp)
    
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
    mov %rax, 120(%rbp)
    jmp ClearStack
  Error2:
    //Set Error to 2         # Size1 not multiple of 8
    mov $2, %rax
    mov %rax, 120(%rbp)
    jmp ClearStack
  Error3:
    //Set Error to 3         # Size2 <= 0
    mov $3, %rax
    mov %rax, 120(%rbp)
    jmp ClearStack
  Error4:                    # localIndexCount <= 0
    //Set Error to 4 
    mov $4, %rax
    mov %rax, 120(%rbp)
    jmp ClearStack
