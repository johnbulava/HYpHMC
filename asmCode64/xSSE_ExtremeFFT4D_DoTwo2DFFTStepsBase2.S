/**  void xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2(char* parameters);
*    -->Parameters:                                         Offset       Remark
*          char* workSpace                                  0
*          Complex* slice                                   8
*          long int Size1                                   16
*          long int Size2                                   24
*          long int localIndexCount                         32
*          long int smallerHalfBlockSize1                   40           Half of Edge-length 1 of smaller FFT
*          long int smallerHalfBlockSize2                   48           Half of Edge-length 2 of smaller FFT
*          long int embeddingOne                            56           Every 4 lines xxx cache lines
*          long int embeddingTwo                            64           Every 4 lines xxx cache lines
*          Complex* FFTcomplexFacs1                         72           Pointer to complex numbers for FFT 1 (embedded: Small, Large, Large)
*          Complex* FFTcomplexFacs2                         80           Pointer to complex numbers for FFT 2 (embedded: Small, Large, Large)
*          int** prefetchAddressControlList                 88           Pointer to list of prefetchAddressControls
*          Complex fac2                                     96           = (2,2)
*          long int error                                   112


*          long int TEMPORAL                                128          Counter: Outer-Loop 1
*          long int TEMPORAL                                136          Counter: Outer-Loop 2
*          long int TEMPORAL                                144          Counter: Inner-Loop 1
*          long int TEMPORAL                                152          Counter: Inner-Loop 2
*          long int TEMPORAL                                160          Counter: Local-Index-Loop


*          long int TEMPORAL                                320          Size1 / (4*smallerHalfBlockSize1)
*          long int TEMPORAL                                328          Size2 / (4*smallerHalfBlockSize2)

*          long int TEMPORAL                                336          1*localIndexCount*smallerHalfBlockSize2*16
*          long int TEMPORAL                                344          2*localIndexCount*smallerHalfBlockSize2*16
*          long int TEMPORAL                                352          3*localIndexCount*smallerHalfBlockSize2*16

*          long int TEMPORAL                                360          embeddingTwo*64
*          long int TEMPORAL                                368          1*(Size2*localIndexCount*16 + embeddingOne*64)
*          long int TEMPORAL                                376          1*smallerHalfBlockSize1*(Size2*localIndexCount*16 + embeddingOne*64)
*          long int TEMPORAL                                384          4*smallerHalfBlockSize1*(Size2*localIndexCount*16 + embeddingOne*64) + smallerHalfBlockSize1*embeddingTwo*64
*          long int TEMPORAL                                392          2*smallerHalfBlockSize1*(Size2*localIndexCount*16 + embeddingOne*64)
*          long int TEMPORAL                                400          3*smallerHalfBlockSize1*(Size2*localIndexCount*16 + embeddingOne*64)
*          
*
*
*    -->Target: 64-Bit AMD64 architecture with GNU g++
*
*    -->Workspace Demand: 0
*
*    -->Description:
*          This routine performs two FFT-steps with base power 2.
*          Input and Output is the 2D-slice-array.
*          Size1, size2 are the sizes of the 2D-slice and must be
*          multiples of 4*smallerHalfBlockSize1 and 4*smallerHalfBlockSize2,
*          respectively. PrefetchAddressControl contains 
*          prefetch-address information for efficient prefetching
*          and must be provided by the calling routine. FFTcomplexFacs
*          contains the complex-numbers for the FFT-step in form (Re, Re), (-Im, Im).
*          LocalIndexCount is the number of nested lattices. EmbeddingOne
*          gives number of inserted cach line at end of each line. EmbeddingTwo 
*          gives number of inserted dummy cache-lines in output-slide every 4 lines.
*          Error is set to 0 on success or to an error code if the 
*          routine failed.
**/
#include "xSSE_Global.h"
.file	"xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2.S"
.text
.align 2
.globl xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2
.type  xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2, @function
xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2:
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
    mov 40(%rbp), %rbx      # Check smallerHalfBlockSize1 > 0
    cmp $0, %rbx
    jle Error1
    
    mov 48(%rbp), %rcx      # Check smallerHalfBlockSize2 > 0
    cmp $0, %rcx
    jle Error2
    
    mov 16(%rbp), %rax      # Check Size1 > 0
    cmp $0, %rax
    jle Error3
    test $3, %rax           # Check Size1 multiple of 4
    jg Error4
    shr $2, %rax
    mov $0, %rdx
    div %rbx
    cmp $0, %rdx
    jne Error5              # Check Size1 multiple of 4*smallerHalfBlockSize1
    mov %rax, 320(%rbp)     # Size1 / (4*smallerHalfBlockSize1)
    
    mov 24(%rbp), %rax      # Check Size2 > 0
    cmp $0, %rax
    jle Error6
    test $3, %rax           # Check Size2 multiple of 4
    jg Error7
    shr $2, %rax
    mov $0, %rdx
    div %rcx
    cmp $0, %rdx
    jne Error8             # Check Size2 multiple of 4*smallerHalfBlockSize2
    mov %rax, 328(%rbp)    # Size2 / (4*smallerHalfBlockSize2) 

    mov 32(%rbp), %rcx     # Check localIndexCount > 0
    cmp $0, %rcx
    jle Error9

    mov 48(%rbp), %rax     # smallerHalfBlockSize2     
    mov %rax, %rdx 
    imul %rcx, %rdx
    shl $4, %rdx
    mov %rdx, 336(%rbp)    # 1*localIndexCount*smallerHalfBlockSize2*16
    mov %rdx, %rbx
    add %rdx, %rbx
    mov %rbx, 344(%rbp)    # 2*localIndexCount*smallerHalfBlockSize2*16
    add %rdx, %rbx
    mov %rbx, 352(%rbp)    # 3*localIndexCount*smallerHalfBlockSize2*16

    mov 40(%rbp), %rax     # smallerHalfBlockSize1     
    mov 64(%rbp), %rbx     # embeddingTwo
    shl $6, %rbx
    mov %rbx, 360(%rbp)    # embeddingTwo*64
    shl $4, %rcx
    mov 24(%rbp), %rdx
    imul %rcx, %rdx
    mov 56(%rbp), %rdi     # embeddingOne
    shl $6, %rdi           # embeddingOne*64
    add %rdi, %rdx
    mov %rdx, 368(%rbp)    # 1*(Size2*localIndexCount*16 + embeddingOne*64)
    imul %rax, %rdx
    mov %rdx, 376(%rbp)    # 1*smallerHalfBlockSize1*(Size2*localIndexCount*16 + embeddingOne*64)
    mov %rdx, %rcx
    add %rdx, %rcx
    mov %rcx, 392(%rbp)    # 2*smallerHalfBlockSize1*(Size2*localIndexCount*16 + embeddingOne*64)
    add %rdx, %rcx
    mov %rcx, 400(%rbp)    # 3*smallerHalfBlockSize1*(Size2*localIndexCount*16 + embeddingOne*64)
    
    shl $2, %rdx
    imul %rax, %rbx
    add %rbx, %rdx
    mov %rdx, 384(%rbp)    # 4*smallerHalfBlockSize1*(Size2*localIndexCount*16 + embeddingOne*64) + smallerHalfBlockSize1*embeddingTwo*64


    
    //Main-Loop    
    mov 320(%rbp), %rdi
    mov %rdi, 128(%rbp)            # Counter: Outer-Loop 1    
    OuterLoop1:
      //Load 1-Direction Pointers
      mov 320(%rbp), %rdi
      sub 128(%rbp), %rdi
      mov 384(%rbp), %r8           # 4*smallerHalfBlockSize1*(Size2*localIndexCount*16 + embeddingOne*64) + smallerHalfBlockSize1*embeddingTwo*64
      imul %rdi, %r8
      add 8(%rbp), %r8             # Slice-Start-Address
      mov %r8, %r9
      mov %r8, %r10
      mov %r8, %r11
      add 376(%rbp), %r9           # 1*smallerHalfBlockSize1*Size2*localIndexCount*16
      add 392(%rbp), %r10          # 2*smallerHalfBlockSize1*Size2*localIndexCount*16
      add 400(%rbp), %r11          # 3*smallerHalfBlockSize1*Size2*localIndexCount*16
      

      //Treatment of embeddingTwo
      mov 40(%rbp), %rax          # smallerHalfBlockSize1
      mov %rax, %rbx
      mov %rbx, %rcx
      shr $2, %rcx
      mov 360(%rbp), %rdx         # embeddingTwo*64
      imul %rdx, %rcx
      add %rcx, %r9
      add %rax, %rbx
      mov %rbx, %rcx
      shr $2, %rcx
      imul %rdx, %rcx
      add %rcx, %r10
      add %rax, %rbx
      mov %rbx, %rcx
      shr $2, %rcx
      imul %rdx, %rcx
      add %rcx, %r11
      

      //Load Pointer to Complex FFT-Facs Direction 1
      mov 72(%rbp), %r12       # Pointer to complex numbers for FFT (1-direction)
    
      mov 40(%rbp), %rdi
      mov %rdi, 144(%rbp)          # Counter: Inner-Loop 1
      InnerLoop1:
        //Load Prefetch-Data
	mov 88(%rbp), %rdx
	mov (%rdx), %r14
	add $8, %rdx
	mov %rdx, 88(%rbp)
      
        //Reset 2-Direction Pointers
        mov $0, %rax
	mov 336(%rbp), %rbx        #  1*localIndexCount*smallerHalfBlockSize2*16
	mov 344(%rbp), %rcx        #  2*localIndexCount*smallerHalfBlockSize2*16
	mov 352(%rbp), %rdx        #  3*localIndexCount*smallerHalfBlockSize2*16

        mov 328(%rbp), %rdi
        mov %rdi, 136(%rbp)        # Counter: Outer-Loop 2
        OuterLoop2:
	  //Load Pointer to Complex FFT-Facs Direction 2
	  mov 80(%rbp), %r13       # Pointer to complex numbers for FFT (2-direction)
	
	  mov 48(%rbp), %rdi
	  mov %rdi, 152(%rbp)      # Counter: Inner-Loop 2
          InnerLoop2:
	    mov 32(%rbp), %rdi
	    mov %rdi, 160(%rbp)    # Counter: Local-Index-Loop
            LocalIndexLoop:
              //Load and S2-Trafo
              movapd (%r8, %rdx), %xmm3
              movapd (%r9, %rdx), %xmm7
              movhlps %xmm3, %xmm0
              movlhps %xmm3, %xmm0
              movapd (%r10, %rdx), %xmm11
              movhlps %xmm7, %xmm1
              movlhps %xmm7, %xmm1
              movapd (%r11, %rdx), %xmm15
              movhlps %xmm11, %xmm4
              movlhps %xmm11, %xmm4
              movapd 0(%r13), %xmm8         //FFT-Factor
              movhlps %xmm15, %xmm5
              movlhps %xmm15, %xmm5
	      
	      
              movapd 16(%r13), %xmm9        //FFT-Factor
	      mulpd %xmm8, %xmm3
              movapd (%r8, %rcx), %xmm2
	      mulpd %xmm8, %xmm7
              movapd (%r9, %rcx), %xmm6
	      mulpd %xmm8, %xmm11
              movapd (%r10, %rcx), %xmm10
	      mulpd %xmm8, %xmm15
              movapd (%r11, %rcx), %xmm14
	      
              #ifdef xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_SlidePrefetchUSE
                mov 0(%r14), %rsi
                movsx %esi, %r15
                xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_SlidePrefetchCMD (%r8, %r15)
                shr $32, %rsi
              #endif

	      mulpd %xmm9, %xmm0
              movapd (%r10, %rax), %xmm8
	      mulpd %xmm9, %xmm1
              movapd (%r11, %rax), %xmm12
	      mulpd %xmm9, %xmm4
              movapd (%r11, %rbx), %xmm13
	      mulpd %xmm9, %xmm5
              movapd (%r10, %rbx), %xmm9	      
	      
	      addpd %xmm0, %xmm3
	      addpd %xmm1, %xmm7
	      addpd %xmm4, %xmm11
	      addpd %xmm5, %xmm15
	      
              movhlps %xmm9, %xmm0
              movhlps %xmm13, %xmm1
              movlhps %xmm9, %xmm0
              movlhps %xmm13, %xmm1
	      
	      mulpd 0(%r13), %xmm9
	      mulpd 0(%r13), %xmm13
	      mulpd 16(%r13), %xmm0
	      mulpd 16(%r13), %xmm1
	      
	      addpd %xmm0, %xmm9
	      addpd %xmm1, %xmm13
	      
	      subpd %xmm3, %xmm2
	      subpd %xmm7, %xmm6
	      subpd %xmm11, %xmm10
	      subpd %xmm15, %xmm14
	      subpd %xmm9, %xmm8
	      subpd %xmm13, %xmm12
	      
	      addpd (%r8, %rcx), %xmm3
	      addpd (%r9, %rcx), %xmm7
	      addpd (%r10, %rcx), %xmm11
	      addpd (%r11, %rcx), %xmm15
	      addpd (%r10, %rax), %xmm9
	      addpd (%r11, %rax), %xmm13
	      
	      
              #ifdef xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_SlidePrefetchUSE
                movsx %esi, %r15
                xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_SlidePrefetchCMD (%r8, %r15)
              #endif


	      //S1-Trafo
	      movhlps %xmm6, %xmm0
	      movhlps %xmm7, %xmm1
	      movhlps %xmm14, %xmm4
	      movhlps %xmm15, %xmm5
	      movlhps %xmm6, %xmm0
	      movlhps %xmm7, %xmm1
	      movlhps %xmm14, %xmm4
	      movlhps %xmm15, %xmm5
	      
	      mulpd 0(%r12), %xmm6
	      mulpd 0(%r12), %xmm7
	      mulpd 0(%r12), %xmm14
	      mulpd 0(%r12), %xmm15
	      mulpd 16(%r12), %xmm0
	      mulpd 16(%r12), %xmm1
	      mulpd 16(%r12), %xmm4
	      mulpd 16(%r12), %xmm5
	      
	      addpd %xmm0, %xmm6
	      addpd %xmm1, %xmm7
	      addpd %xmm4, %xmm14
	      addpd %xmm5, %xmm15
	      
	      movhlps %xmm12, %xmm0
	      movhlps %xmm13, %xmm1
	      movlhps %xmm12, %xmm0
	      movlhps %xmm13, %xmm1
	      
	      mulpd 0(%r12), %xmm12
	      mulpd 0(%r12), %xmm13
	      mulpd 16(%r12), %xmm0
	      mulpd 16(%r12), %xmm1
	      
	      addpd %xmm0, %xmm12
	      addpd %xmm1, %xmm13
	      
	      
	      addpd %xmm2, %xmm6
	      addpd %xmm3, %xmm7
	      addpd %xmm10, %xmm14
	      addpd %xmm11, %xmm15
	      addpd %xmm8, %xmm12
	      addpd %xmm9, %xmm13
	      
	      mulpd 96(%rbp), %xmm2
	      mulpd 96(%rbp), %xmm3
	      mulpd 96(%rbp), %xmm10
	      mulpd 96(%rbp), %xmm11
	      mulpd 96(%rbp), %xmm8
	      mulpd 96(%rbp), %xmm9
	      
	      subpd %xmm6, %xmm2
	      subpd %xmm7, %xmm3
	      subpd %xmm14, %xmm10
	      subpd %xmm15, %xmm11
	      subpd %xmm12, %xmm8
	      subpd %xmm13, %xmm9
	      
              #ifdef xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_SlidePrefetchUSE
                mov 8(%r14), %rsi
                movsx %esi, %r15
                xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_SlidePrefetchCMD (%r8, %r15)
                shr $32, %rsi
              #endif

	      
	      //L2-Trafo
	      movhlps %xmm2, %xmm0
	      movhlps %xmm6, %xmm1
	      movhlps %xmm10, %xmm4
	      movhlps %xmm14, %xmm5
	      movlhps %xmm2, %xmm0
	      movlhps %xmm6, %xmm1
	      movlhps %xmm10, %xmm4
	      movlhps %xmm14, %xmm5
	      
	      mulpd 64(%r13), %xmm2
	      mulpd 64(%r13), %xmm6
	      mulpd 64(%r13), %xmm10
	      mulpd 64(%r13), %xmm14
	      mulpd 80(%r13), %xmm0
	      mulpd 80(%r13), %xmm1
	      mulpd 80(%r13), %xmm4
	      mulpd 80(%r13), %xmm5
	      
	      addpd %xmm0, %xmm2
	      addpd %xmm1, %xmm6
	      addpd %xmm4, %xmm10
	      addpd %xmm5, %xmm14
	      
	      movhlps %xmm3, %xmm0
	      movhlps %xmm7, %xmm1
	      movhlps %xmm11, %xmm4
	      movhlps %xmm15, %xmm5
	      movlhps %xmm3, %xmm0
	      movlhps %xmm7, %xmm1
	      movlhps %xmm11, %xmm4
	      movlhps %xmm15, %xmm5
	      
	      mulpd 32(%r13), %xmm3
	      mulpd 32(%r13), %xmm7
	      mulpd 32(%r13), %xmm11
	      mulpd 32(%r13), %xmm15
	      mulpd 48(%r13), %xmm0
	      mulpd 48(%r13), %xmm1
	      mulpd 48(%r13), %xmm4
	      mulpd 48(%r13), %xmm5
	      
	      addpd %xmm0, %xmm3
	      addpd %xmm1, %xmm7
	      addpd %xmm4, %xmm11
	      addpd %xmm5, %xmm15
	      
	      
	      addpd %xmm8, %xmm10
	      addpd %xmm9, %xmm11
	      addpd %xmm12, %xmm14
	      addpd %xmm13, %xmm15
	      
	      mulpd 96(%rbp), %xmm8
	      mulpd 96(%rbp), %xmm9
	      mulpd 96(%rbp), %xmm12
	      mulpd 96(%rbp), %xmm13
	      
	      subpd %xmm10, %xmm8
	      subpd %xmm11, %xmm9
	      subpd %xmm14, %xmm12
	      subpd %xmm15, %xmm13


              #ifdef xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_SlidePrefetchUSE
                movsx %esi, %r15
                xSSE_ExtremeFFT4D_DoTwo2DFFTStepsBase2_SlidePrefetchCMD (%r8, %r15)
              #endif

	      
	      //L1-Trafo
	      movhlps %xmm8, %xmm0
	      movhlps %xmm9, %xmm1
	      movhlps %xmm10, %xmm4
	      movhlps %xmm11, %xmm5
	      movlhps %xmm8, %xmm0
	      movlhps %xmm9, %xmm1
	      movlhps %xmm10, %xmm4
	      movlhps %xmm11, %xmm5
	      
	      mulpd 64(%r12), %xmm8
	      mulpd 64(%r12), %xmm9
	      mulpd 64(%r12), %xmm10
	      mulpd 64(%r12), %xmm11
	      mulpd 80(%r12), %xmm0
	      mulpd 80(%r12), %xmm1
	      mulpd 80(%r12), %xmm4
	      mulpd 80(%r12), %xmm5
	      
	      addpd %xmm0, %xmm8
	      addpd %xmm1, %xmm9
	      addpd %xmm4, %xmm10
	      addpd %xmm5, %xmm11
	      
	      movhlps %xmm12, %xmm0
	      movhlps %xmm13, %xmm1
	      movhlps %xmm14, %xmm4
	      movhlps %xmm15, %xmm5
	      movlhps %xmm12, %xmm0
	      movlhps %xmm13, %xmm1
	      movlhps %xmm14, %xmm4
	      movlhps %xmm15, %xmm5
	      
	      mulpd 32(%r12), %xmm12
	      mulpd 32(%r12), %xmm13
	      mulpd 32(%r12), %xmm14
	      mulpd 32(%r12), %xmm15
	      mulpd 48(%r12), %xmm0
	      mulpd 48(%r12), %xmm1
	      mulpd 48(%r12), %xmm4
	      mulpd 48(%r12), %xmm5
	      
	      addpd %xmm0, %xmm12
	      addpd %xmm1, %xmm13
	      addpd %xmm4, %xmm14
	      addpd %xmm5, %xmm15
	      
	      //Rest-Trafo of 0,1,4,5 - Block
              movapd (%r8, %rbx), %xmm1
              movapd (%r9, %rbx), %xmm5
	      
	      movhlps %xmm1, %xmm0
	      movhlps %xmm5, %xmm4
	      movlhps %xmm1, %xmm0
	      movlhps %xmm5, %xmm4
	      
	      mulpd 0(%r13), %xmm1
	      mulpd 0(%r13), %xmm5
	      mulpd 16(%r13), %xmm0
	      mulpd 16(%r13), %xmm4
	      
	      addpd %xmm0, %xmm1
	      addpd %xmm4, %xmm5
	      
              movapd (%r9, %rax), %xmm4
	      
	      subpd %xmm5, %xmm4
	      addpd (%r9, %rax), %xmm5
	      
	      movhlps %xmm4, %xmm0
	      movlhps %xmm4, %xmm0
	      
	      mulpd 0(%r12), %xmm4
	      mulpd 16(%r12), %xmm0
	      
	      addpd %xmm0, %xmm4
	      
	      movhlps %xmm5, %xmm0
	      movlhps %xmm5, %xmm0
	      
	      mulpd 0(%r12), %xmm5
	      mulpd 16(%r12), %xmm0
	      
	      addpd %xmm0, %xmm5
	      
              movapd (%r8, %rax), %xmm0
	      
	      subpd %xmm1, %xmm0
	      addpd (%r8, %rax), %xmm1
	      
	      addpd %xmm0, %xmm4
	      addpd %xmm1, %xmm5
	      mulpd 96(%rbp), %xmm0
	      mulpd 96(%rbp), %xmm1
	      subpd %xmm4, %xmm0
	      subpd %xmm5, %xmm1
	      
	      addpd %xmm0, %xmm2
	      addpd %xmm1, %xmm3
	      addpd %xmm4, %xmm6
	      addpd %xmm5, %xmm7
	      
	      mulpd 96(%rbp), %xmm0
	      mulpd 96(%rbp), %xmm1
	      mulpd 96(%rbp), %xmm4
	      mulpd 96(%rbp), %xmm5
	      
	      subpd %xmm2, %xmm0
	      subpd %xmm3, %xmm1
	      subpd %xmm6, %xmm4
	      subpd %xmm7, %xmm5

	      
	      addpd %xmm0, %xmm8
	      mulpd 96(%rbp), %xmm0
	      addpd %xmm1, %xmm9
	      mulpd 96(%rbp), %xmm1
	      addpd %xmm2, %xmm10
	      mulpd 96(%rbp), %xmm2
	      addpd %xmm3, %xmm11
	      mulpd 96(%rbp), %xmm3

	      movapd %xmm8, (%r9, %rdx)
	      addpd %xmm4, %xmm12
	      movapd %xmm9, (%r9, %rcx)
	      addpd %xmm5, %xmm13
	      movapd %xmm10, (%r9, %rbx)
	      addpd %xmm6, %xmm14
	      movapd %xmm11, (%r9, %rax)	      
	      addpd %xmm7, %xmm15

	      mulpd 96(%rbp), %xmm4
              subpd %xmm8, %xmm0	      
	      mulpd 96(%rbp), %xmm5
              subpd %xmm9, %xmm1     
	      mulpd 96(%rbp), %xmm6
              subpd %xmm10, %xmm2	      
	      mulpd 96(%rbp), %xmm7
              subpd %xmm11, %xmm3

	      movapd %xmm12, (%r8, %rdx)
	      movapd %xmm13, (%r8, %rcx)
	      movapd %xmm14, (%r8, %rbx)
	      movapd %xmm15, (%r8, %rax)

	      movapd %xmm0, (%r11, %rdx)
              subpd %xmm12, %xmm4	      
	      movapd %xmm1, (%r11, %rcx)
              subpd %xmm13, %xmm5	      
	      movapd %xmm2, (%r11, %rbx)
              subpd %xmm14, %xmm6	      
	      movapd %xmm3, (%r11, %rax)
              subpd %xmm15, %xmm7	
	      
	      movapd %xmm4, (%r10, %rdx)
	      movapd %xmm5, (%r10, %rcx)
	      movapd %xmm6, (%r10, %rbx)
	      movapd %xmm7, (%r10, %rax)
	      
	      add $16, %r14
              add $16, %rax
              add $16, %rbx
              add $16, %rcx
              add $16, %rdx
              sub $1, 160(%rbp)
            jg LocalIndexLoop
            //Update Fac-Pointers Direction 2
	    add $96, %r13
	    
            sub $1, 152(%rbp)
          jg InnerLoop2
	  //Update 2-Direction Pointers
	  mov 352(%rbp), %rdi        # 3*localIndexCount*smallerHalfBlockSize2*16
	  add %rdi, %rax
	  add %rdi, %rbx
	  add %rdi, %rcx
	  add %rdi, %rdx
	  
          sub $1, 136(%rbp)
        jg OuterLoop2
	//Update Fac-Pointers Direction 1
        add $96, %r12
	
	//Update 1-Direction Pointers
	mov 368(%rbp), %rdi          # 1*Size2*localIndexCount*16
	add %rdi, %r8
	add %rdi, %r9
	add %rdi, %r10
	add %rdi, %r11

	
        //Treatment of embeddingTwo*64
        mov 40(%rbp), %rdi           # smallerHalfBlockSize1
        sub 144(%rbp), %rdi          # Counter: Inner-Loop 1
        inc %rdi
	mov $0, %rdx
	test $3, %rdi
	cmove 360(%rbp), %rdx        # embeddingTwo*64
	add %rdx, %r8
	add 40(%rbp), %rdi           # smallerHalfBlockSize1
	mov $0, %rdx
	test $3, %rdi
	cmove 360(%rbp), %rdx
	add %rdx, %r9
	add 40(%rbp), %rdi           # smallerHalfBlockSize1
	mov $0, %rdx
	test $3, %rdi
	cmove 360(%rbp), %rdx
	add %rdx, %r10
	add 40(%rbp), %rdi           # smallerHalfBlockSize1
	mov $0, %rdx
	test $3, %rdi
	cmove 360(%rbp), %rdx
	add %rdx, %r11
      
        sub $1, 144(%rbp)
      jg InnerLoop1
	
      sub $1, 128(%rbp)
    jg OuterLoop1



  #Set Error to Zero 
    mov $0, %rax
    mov %rax, 112(%rbp)
    
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
    //Set Error to 1         # smallerHalfBlockSize1 <= 0
    mov $1, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error2:
    //Set Error to 2         # smallerHalfBlockSize2 <= 0
    mov $2, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error3:
    //Set Error to 3         # Size1 <= 0
    mov $3, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error4:                    # Size1 not multiple of 4
    //Set Error to 4 
    mov $4, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error5:                    # Size1 not multiple of 4*smallerHalfBlockSize1
    //Set Error to 5 
    mov $5, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error6:                    # Size2 <= 0
    //Set Error to 6 
    mov $6, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error7:                    # Size2 not multiple of 4
    //Set Error to 7 
    mov $7, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error8:                    # Size2 not multiple of 4*smallerHalfBlockSize2
    //Set Error to 8 
    mov $8, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error9:                    # localIndexCount <= 0
    //Set Error to 9 
    mov $9, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
