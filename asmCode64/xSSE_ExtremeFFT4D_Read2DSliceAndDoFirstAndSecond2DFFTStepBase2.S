/**  void xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2(char* parameters);
*    -->Parameters:                                         Offset       Remark
*          char* workSpace                                  0
*          Complex* input                                   8
*          Complex* slice                                   16
*          long int Size1                                   24
*          long int Size2                                   32
*          long int localIndexCount                         40
*          long int* bitInvertedInputAddressAuxData1        48           For input slices
*          long int* bitInvertedOutputAddressAuxData2       56           For output slice
*          int** prefetchAddressControlList                 64           Pointer to list of prefetchAddressControls
*          long int add2                                    72           Increment in 2-direction
*          long int embeddingOne                            80           Every line xxx cache lines
*          long int embeddingTwo                            88           Every 4 lines xxx cache lines
*          Complex facpmI (bitmask)                         96           Bitmask for multiplication with plus/minus I
*          Complex fac2                                     112          = (2,2)
*          long int error                                   128
*
*          long int TEMPORAL                                136          1*Size/4 * add2
*          long int TEMPORAL                                144          2*Size/4 * add2
*          long int TEMPORAL                                152          3*Size/4 * add2
*          long int TEMPORAL                                160          loop counter: inner-loop
*          long int TEMPORAL                                168          loop counter: Size1-loop
*          long int TEMPORAL                                176          loop counter: Size2-loop
*          long int TEMPORAL                                184          output-increment = 1*localIndex*16
*          long int TEMPORAL                                192          output-increment = 2*localIndex*16
*          long int TEMPORAL                                200          output-increment = 3*localIndex*16 
*          long int TEMPORAL                                208          output-increment = 1*(Size2*localIndex*16 + embeddingOne*64)
*          long int TEMPORAL                                216          output-increment = 2*(Size2*localIndex*16 + embeddingOne*64)
*          long int TEMPORAL                                224          output-increment = 3*(Size2*localIndex*16 + embeddingOne*64)
*          long int TEMPORAL                                232          output-increment = 4*(Size2*localIndex*16 + embeddingOne*64) + embeddingTwo*64
*          long int TEMPORAL                                240          add2 - 16*localIndexCount
*          
*
*
*    -->Target: 64-Bit AMD64 architecture with GNU g++
*
*    -->Workspace Demand: 0
*
*    -->Description:
*          This routine copies a 2D-sub-lattice structure of 
*          complex numbers (double, double) from a higher-
*          dimensional lattice structure into the slice-array.
*          Simulataneously, it performs the first two FFT-steps.
*          FacpmI is a bitmask for +-i, depending on the FFT-mode.
*          Size1, size2 are the sizes of the 2D-slice and must be
*          multiples of 4. PrefetchAddressControl contains 
*          prefetch-address information for efficient prefetching
*          and must be provided by the calling routine.
*          LocalIndexCount is the number of nested lattices.
*          Furthermore, the
*          data addresses are bit-inverted by using the lookup
*          tables bitInvertedInputAddressAuxData1, bitInvertedOutputAddressAuxData2.
*          First one refers to first index of input, second to second
*          index of output data. EmbeddingOne gives number of inserted
*          cache lines per line. EmbeddingTwo gives number of inserted
*          dummy cache-lines in output-slide every 4 lines.
*          Error is set to 0 on success or to an error code if the 
*          routine failed.
**/
#include "xSSE_Global.h"
.file	"xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2.S"
.text
.align 2
.globl xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2
.type  xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2, @function
xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2:
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
    //Load Output Pointer
    mov 16(%rbp), %rdi      # Slices        

    //Preparing Constants for Main loop And Error-Check
    mov 24(%rbp), %rax      # Check Size1 > 0
    cmp $0, %rax
    jle Error1
    test $3, %rax           # Check Size1 multiple of 4
    jg Error2
    dec %rax
    shl $3, %rax
    mov %rax, 24(%rbp)
    
    mov 32(%rbp), %rbx      # Check Size2 > 0
    cmp $0, %rbx
    jle Error3
    test $3, %rbx           # Check Size2 multiple of 4
    jg Error4
    dec %rbx
    shl $3, %rbx
    mov %rbx, 32(%rbp)
    
    mov 40(%rbp), %rcx      # Check localIndexCount > 0
    cmp $0, %rcx
    jle Error5
    dec %rcx
    mov %rcx, 40(%rbp)
    
    mov 40(%rbp), %rcx
    inc %rcx
    shl $4, %rcx
    mov %rcx, 184(%rbp)     # = 1*localIndex*16
    mov %rcx, %rdx
    shl $1, %rdx
    mov %rdx, 192(%rbp)     # = 2*localIndex*16    
    add %rcx, %rdx
    mov %rdx, 200(%rbp)     # = 3*localIndex*16    
    
    mov 32(%rbp), %rbx
    shr $3, %rbx
    inc %rbx
    imul %rbx, %rcx
    mov 80(%rbp), %rax
    shl $6, %rax
    add %rax, %rcx
    mov %rcx, 208(%rbp)      # = 1*(Size2*localIndex*16 + embeddingOne*64)
    mov %rcx, %rdx
    shl $1, %rcx
    mov %rcx, 216(%rbp)      # = 2*(Size2*localIndex*16 + embeddingOne*64)
    add %rcx, %rdx
    mov %rdx, 224(%rbp)      # = 3*(Size2*localIndex*16 + embeddingOne*64)   
    shl $1, %rcx    
    mov 88(%rbp), %rbx
    shl $6, %rbx
    add %rbx, %rcx
    mov %rcx, 232(%rbp)      # = 4*(Size2*localIndex*16 + embeddingOne*64) + embeddingTwo*64

    mov 32(%rbp), %rcx
    add $8, %rcx
    imul 72(%rbp), %rcx
    shr $5, %rcx
    mov %rcx, 136(%rbp)     # 1*Size/4 * add2
    mov %rcx, %rdx
    add %rdx, %rcx
    mov %rcx, 144(%rbp)     # 2*Size/4 * add2
    add %rdx, %rcx    
    mov %rcx, 152(%rbp)     # 3*Size/4 * add2
    
    mov 72(%rbp), %rbx
    mov 40(%rbp), %rcx
    inc %rcx
    shl $4, %rcx
    sub %rcx, %rbx
    mov %rbx, 240(%rbp)     # add2 - 16*localIndexCount
    
    
    
    //Main-Loop    
    mov 24(%rbp), %rax
    mov %rax, 168(%rbp)     # Loop-Counter
    MainLoop1:
      //Load Bit-Inverted Read-Addresses 1 and Calculate Prefetch-Startposition
      mov 48(%rbp), %r13 
      add 168(%rbp), %r13
      mov 0(%r13),   %r8
      mov -8(%r13),  %r9
      mov -16(%r13), %r10
      mov -24(%r13), %r11
      mov 8(%rbp), %rsi      # Source
      add %rsi, %r8
      add %rsi, %r9
      add %rsi, %r10
      add %rsi, %r11
      
      
      //Load Prefetch-Control-Data and set loop-counter            
      mov 64(%rbp), %rax
      mov 168(%rbp), %rbx      
      #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_InputSlidePrefetchUSE
        mov (%rax, %rbx), %r15    # PrefetchControlOffsets 
      #endif
      #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_OutputSlidePrefetchUSE
        mov -8(%rax, %rbx), %r14  # PrefetchControlOffsets      
      #endif      
      
      mov 32(%rbp), %rax
      mov %rax, 176(%rbp)       # Loop-Counter
      //Initialize Read-Addresses 2       
      mov $0, %rax
      mov 144(%rbp), %rbx
      mov 136(%rbp), %rcx
      mov 152(%rbp), %rdx  
      mov $0, %rsi    
      MainLoop2:
        mov 40(%rbp), %r13 
        mov %r13, 160(%rbp)    # Inner Loop Counter
        mov 184(%rbp), %r13    # output-increment = 1*localIndex*16
        MainLoop3:
	  movapd 0(%r8, %rax), %xmm0
	  movapd 0(%r8, %rbx), %xmm1
          addpd %xmm0, %xmm1	  	  
	  movapd 0(%r8, %rcx), %xmm2
	  movapd 0(%r8, %rdx), %xmm3
          addpd %xmm2, %xmm3	  


	  movapd 0(%r9, %rax), %xmm4
	  mulpd 112(%rbp), %xmm0
	  movapd 0(%r9, %rbx), %xmm5	  
	  mulpd 112(%rbp), %xmm2
	  movapd 0(%r9, %rcx), %xmm6
	  movapd 0(%r9, %rdx), %xmm7
          addpd %xmm4, %xmm5


          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_OutputSlidePrefetchUSE
            mov 0(%r14), %rsi
            movsx %esi, %r12
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_OutputSlidePrefetchCMD (%rdi, %r12)
            shr $32, %rsi
          #endif

	  
	  movapd 0(%r10, %rax), %xmm8
          addpd %xmm6, %xmm7
	  movapd 0(%r10, %rbx), %xmm9
	  movapd 0(%r10, %rcx), %xmm10
	  movapd 0(%r10, %rdx), %xmm11
          subpd %xmm1, %xmm0
          subpd %xmm3, %xmm2
	  mulpd 112(%rbp), %xmm4
	  mulpd 112(%rbp), %xmm6
          addpd %xmm8, %xmm9
          addpd %xmm10, %xmm11


          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_OutputSlidePrefetchUSE
            movsx %esi, %r12
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_OutputSlidePrefetchCMD (%rdi, %r12)
          #endif

	  
	  movapd 0(%r11, %rax), %xmm12
	  movapd 0(%r11, %rbx), %xmm13
	  movapd 0(%r11, %rcx), %xmm14
	  movapd 0(%r11, %rdx), %xmm15
          subpd %xmm5, %xmm4
          subpd %xmm7, %xmm6
          addpd %xmm12, %xmm13
          addpd %xmm14, %xmm15
	  mulpd 112(%rbp), %xmm8
	  mulpd 112(%rbp), %xmm10
	  mulpd 112(%rbp), %xmm12
	  mulpd 112(%rbp), %xmm14
          subpd %xmm9, %xmm8
          subpd %xmm11, %xmm10
          subpd %xmm13, %xmm12
          subpd %xmm15, %xmm14
	  

          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_InputSlidePrefetchUSE
            mov 0(%r15), %rsi
            movsx %esi, %r12
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_InputSlidePrefetchCMD (%r8, %r12)
            shr $32, %rsi
          #endif


          addpd %xmm0, %xmm4
          addpd %xmm1, %xmm5
	  mulpd 112(%rbp), %xmm0
	  mulpd 112(%rbp), %xmm1
          addpd %xmm2, %xmm6
          addpd %xmm3, %xmm7
	  mulpd 112(%rbp), %xmm2
	  mulpd 112(%rbp), %xmm3
          addpd %xmm8, %xmm12
          addpd %xmm9, %xmm13
	  mulpd 112(%rbp), %xmm8
	  mulpd 112(%rbp), %xmm9
          addpd %xmm10, %xmm14
          addpd %xmm11, %xmm15
	  mulpd 112(%rbp), %xmm10
	  mulpd 112(%rbp), %xmm11
          subpd %xmm4, %xmm0
          subpd %xmm5, %xmm1
          subpd %xmm6, %xmm2
          subpd %xmm7, %xmm3
          subpd %xmm12, %xmm8
          subpd %xmm13, %xmm9
          subpd %xmm14, %xmm10
          subpd %xmm15, %xmm11

          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_InputSlidePrefetchUSE
            movsx %esi, %r12
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_InputSlidePrefetchCMD (%r8, %r12)
          #endif

          addpd %xmm1, %xmm3
          addpd %xmm5, %xmm7
	  mulpd 112(%rbp), %xmm1
	  mulpd 112(%rbp), %xmm5
          addpd %xmm9, %xmm11
          addpd %xmm13, %xmm15
	  mulpd 112(%rbp), %xmm9
	  mulpd 112(%rbp), %xmm13
          subpd %xmm3, %xmm1
          subpd %xmm7, %xmm5
          subpd %xmm11, %xmm9
          subpd %xmm15, %xmm13

          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_OutputSlidePrefetchUSE
            mov 8(%r14), %rsi
            movsx %esi, %r12
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_OutputSlidePrefetchCMD (%rdi, %r12)
            shr $32, %rsi
          #endif


          addpd %xmm4, %xmm12
          addpd %xmm5, %xmm13
	  mulpd 112(%rbp), %xmm4
	  mulpd 112(%rbp), %xmm5
          addpd %xmm6, %xmm14
          addpd %xmm7, %xmm15
	  mulpd 112(%rbp), %xmm6
	  mulpd 112(%rbp), %xmm7
          subpd %xmm12, %xmm4
          subpd %xmm13, %xmm5
          subpd %xmm14, %xmm6
          subpd %xmm15, %xmm7

          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_OutputSlidePrefetchUSE
            movsx %esi, %r12
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_OutputSlidePrefetchCMD (%rdi, %r12)
          #endif
	  

          //Writeback 15, 13, 7, 5
	  mov 216(%rbp),%rsi
	  mov 192(%rbp),%r12     # output-increment = 2*localIndex*16
	  movapd %xmm15, (%rdi)
	  movapd %xmm13, (%rdi,%r12)
	  add %rdi, %r12
	  movapd %xmm7, (%rdi,%rsi)
	  movapd %xmm5, (%r12,%rsi)

	  
          movlhps %xmm2, %xmm5
          movhlps %xmm2, %xmm5
          movlhps %xmm6, %xmm7
          movhlps %xmm6, %xmm7
          movlhps %xmm10, %xmm13
          movhlps %xmm10, %xmm13
          movlhps %xmm14, %xmm15
          movhlps %xmm14, %xmm15

          xorpd 96(%rbp), %xmm5
          xorpd 96(%rbp), %xmm7
          xorpd 96(%rbp), %xmm13
	  xorpd 96(%rbp), %xmm15

          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_InputSlidePrefetchUSE
            mov 8(%r15), %rsi
            movsx %esi, %r12
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_InputSlidePrefetchCMD (%r8, %r12)
            shr $32, %rsi
          #endif
	  
          addpd %xmm0, %xmm5
          addpd %xmm4, %xmm7
	  mulpd 112(%rbp), %xmm0
	  mulpd 112(%rbp), %xmm4
          addpd %xmm8, %xmm13
          addpd %xmm12, %xmm15
	  mulpd 112(%rbp), %xmm8
	  mulpd 112(%rbp), %xmm12
          subpd %xmm5, %xmm0
          subpd %xmm7, %xmm4
          subpd %xmm13, %xmm8
          subpd %xmm15, %xmm12

          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_InputSlidePrefetchUSE
            movsx %esi, %r12
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstAndSecond2DFFTStepBase2_InputSlidePrefetchCMD (%r8, %r12)
          #endif


          //Writeback 15, 12, 7, 4
	  mov 216(%rbp),%rsi
	  mov 200(%rbp),%r12
          movapd %xmm15, (%rdi, %r13)	  
	  movapd %xmm12, (%rdi, %r12)
	  add %rdi, %r12
	  movapd %xmm4, (%r12, %rsi)
	  add %rdi, %rsi
	  movapd %xmm7, (%rsi, %r13)
	  

          movlhps %xmm8, %xmm4
          movhlps %xmm8, %xmm4
          movlhps %xmm9, %xmm7
          movhlps %xmm9, %xmm7
          movlhps %xmm13, %xmm12
          movhlps %xmm13, %xmm12
          movlhps %xmm11, %xmm15
          movhlps %xmm11, %xmm15

          xorpd 96(%rbp), %xmm4
          xorpd 96(%rbp), %xmm7
          xorpd 96(%rbp), %xmm12
	  xorpd 96(%rbp), %xmm15
	  
          addpd %xmm0, %xmm4
          addpd %xmm1, %xmm7
          addpd %xmm5, %xmm12
          addpd %xmm3, %xmm15
	  
	  mulpd 112(%rbp), %xmm0
	  mulpd 112(%rbp), %xmm1
	  mulpd 112(%rbp), %xmm5
	  mulpd 112(%rbp), %xmm3
	  
          subpd %xmm4, %xmm0
          subpd %xmm7, %xmm1
          subpd %xmm12, %xmm5
          subpd %xmm15, %xmm3


          //Writeback rest
	  mov 208(%rbp),%r12
	  mov 216(%rbp),%rsi
	  add %rdi, %r12	  
          movapd %xmm15, (%r12)
	  movapd %xmm3, (%r12, %rsi)
	  add %r13, %r12	  
	  movapd %xmm12, (%r12)
	  movapd %xmm5, (%r12, %rsi)
	  add %r13, %r12	  
	  movapd %xmm7, (%r12)
	  movapd %xmm1, (%r12, %rsi)
	  add %r13, %r12	  
	  movapd %xmm4, (%r12)
	  movapd %xmm0, (%r12, %rsi)
	  
	  
	  add $16, %r14
	  add $16, %r15
	  add $16, %rdi
          add $16, %rax
          add $16, %rbx
          add $16, %rcx
          add $16, %rdx	  
          sub $1, 160(%rbp)
        jge MainLoop3
	
	mov 176(%rbp), %r12
	mov 56(%rbp), %r13
        add (%r13, %r12), %rdi

        mov 240(%rbp), %r12
        add %r12, %rax
        add %r12, %rbx
        add %r12, %rcx
        add %r12, %rdx

        sub $32, 176(%rbp)
      jge MainLoop2
      add 232(%rbp), %rdi
      
      sub $32, 168(%rbp)
    jge MainLoop1


  #Set Error to Zero 
    mov $0, %rax
    mov %rax, 128(%rbp)
    
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
    mov %rax, 128(%rbp)
    jmp ClearStack
  Error2:
    //Set Error to 2         # Size1 not multiple of 4
    mov $2, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack
  Error3:                    # Size2 <= 0
    //Set Error to 3 
    mov $3, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack
  Error4:                    # Size2 not multiple of 4
    //Set Error to 4 
    mov $4, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack
  Error5:                    # localIndexCount <= 0
    //Set Error to 5 
    mov $5, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack
