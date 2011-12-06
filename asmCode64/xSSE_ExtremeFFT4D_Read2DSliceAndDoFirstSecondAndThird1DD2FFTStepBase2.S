/**  void xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2(char* parameters);
*    -->Parameters:                                         Offset       Remark
*          char* workSpace                                  0
*          Complex* input                                   8
*          Complex* slice                                   16
*          long int Size1                                   24
*          long int Size2                                   32
*          long int localIndexCount                         40
*          long int* bitInvertedOutputAddressAuxData1       48           For output slices
*          long int* bitInvertedOutputAddressAuxData2       56           For output slice
*          int** prefetchAddressControlList                 64           Pointer to list of prefetchAddressControls
*          long int add1                                    72           Increment in 1-direction
*          long int add2                                    80           Increment in 2-direction
*          Complex* FFTcomplexFacs                          88           Pointer to complex numbers for third FFT 
*          Complex facpmI (bitmask)                         96           Bitmask for multiplication with plus/minus I
*          long int error                                   112
*
*          long int TEMPORAL                                136          1*Size2/8 * add2
*          long int TEMPORAL                                144          2*Size2/8 * add2
*          long int TEMPORAL                                152          3*Size2/8 * add2
*          long int TEMPORAL                                160          4*Size2/8 * add2
*          long int TEMPORAL                                168          loop counter: inner-loop
*          long int TEMPORAL                                176          loop counter: Size1-loop
*          long int TEMPORAL                                184          loop counter: Size2-loop
*          long int TEMPORAL                                192          1*localIndex*16
*          long int TEMPORAL                                200          2*localIndex*16
*          long int TEMPORAL                                208          3*localIndex*16 
*          long int TEMPORAL                                216          4*localIndex*16 
*          long int TEMPORAL                                224          input-increment =  add2 - 16*localIndexCount
*          long int TEMPORAL                                232          input-increment =  add1 - (Size2/8)*add2
*
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
*          Simulataneously, it performs the first, second, and third
*          1D-FFT-steps in direction 2.
*          FacpmI is a bitmask for +-i, depending on the FFT-mode.
*          Size1, size2 are the sizes of the 2D-slice and must be
*          multiples of 8. PrefetchAddressControl contains 
*          prefetch-address information for efficient prefetching
*          and must be provided by the calling routine.
*          LocalIndexCount is the number of nested lattices.
*          Furthermore, the
*          data addresses are bit-inverted by using the lookup
*          tables bitInvertedOutputAddressAuxData1, bitInvertedOutputAddressAuxData2.
*          First one refers to first index of output, second to second
*          index of output data. FFTcomplexFacs
*          contains complex numbers for third FFT-step.
*          Error is set to 0 on success or to an error code if the 
*          routine failed.
**/
#include "xSSE_Global.h"
.file	"xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2.S"
.text
.align 2
.globl xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2
.type  xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2, @function
xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2:
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
    mov 8(%rbp), %rsi       # Input
    mov 16(%rbp), %rdi      # Slices        

    //Preparing Constants for Main loop And Error-Check
    mov 24(%rbp), %rax      # Check Size1 > 0
    cmp $0, %rax
    jle Error1
    test $0, %rax           # Check Size1 multiple of 1
    jg Error2
    dec %rax
    shl $3, %rax
    mov %rax, 24(%rbp)
    
    mov 32(%rbp), %rbx      # Check Size2 > 0
    cmp $0, %rbx
    jle Error3
    test $7, %rbx           # Check Size2 multiple of 8
    jg Error4
    dec %rbx
    shl $3, %rbx
    mov %rbx, 32(%rbp)
    
    mov 40(%rbp), %rcx      # Check localIndexCount > 0
    cmp $0, %rcx
    jle Error5
    dec %rcx
    mov %rcx, 40(%rbp)


    mov 32(%rbp), %rbx
    shr $3, %rbx
    inc %rbx
    shr $3, %rbx
    mov 80(%rbp), %rcx    
    imul %rbx, %rcx
    mov %rcx, 136(%rbp)      # = 1*Size2/8 * add2
    mov %rcx, %rdx
    add %rdx, %rcx
    mov %rcx, 144(%rbp)      # = 2*Size2/8 * add2
    add %rdx, %rcx
    mov %rcx, 152(%rbp)      # = 3*Size2/8 * add2
    add %rdx, %rcx
    mov %rcx, 160(%rbp)      # = 4*Size2/8 * add2

    mov 40(%rbp), %rcx
    inc %rcx
    shl $4, %rcx
    mov %rcx, 192(%rbp)     # = 1*localIndex*16    
    mov %rcx, %rdx
    add %rdx, %rcx
    mov %rcx, 200(%rbp)     # = 2*localIndex*16    
    add %rdx, %rcx
    mov %rcx, 208(%rbp)     # = 3*localIndex*16    
    add %rdx, %rcx
    mov %rcx, 216(%rbp)     # = 4*localIndex*16    
    
    mov 80(%rbp), %rcx
    sub %rdx, %rcx
    mov %rcx, 224(%rbp)     # = add2 - 16*localIndexCount 
    mov 72(%rbp), %rcx    
    sub 136(%rbp), %rcx
    mov %rcx, 232(%rbp)     # = add1 - (Size2/8)*add2
    

    //Load constant address-offsets and constants
    mov 136(%rbp), %r8
    mov 144(%rbp), %r9
    mov 152(%rbp), %r10
    mov 192(%rbp), %r12
    mov 200(%rbp), %r13
    mov 208(%rbp), %r14
    movapd 96(%rbp), %xmm15	
    mov 88(%rbp), %rax
    movapd 0(%rax), %xmm13
    movapd 16(%rax), %xmm14
    
    
    
    //Main-Loop    
    mov 24(%rbp), %rax
    mov %rax, 176(%rbp)     # Loop-Counter
    MainLoop1:
      //Load input-addresses
      mov 160(%rbp), %r11
      add %rsi, %r11
    
      //Load Bit-Inverted Write-Addresses 1
      mov 48(%rbp),  %rax 
      add 176(%rbp), %rax
      mov 16(%rbp), %rdi      # Slices              
      add 0(%rax),   %rdi
      mov 216(%rbp), %r15
      add %rdi, %r15
      
      //Load Prefetch-Control-Data and set loop-counter            
      mov 64(%rbp), %rax
      mov 176(%rbp), %rbx
      #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_InputSlidePrefetchUSE
        mov 0(%rax, %rbx, 2), %rcx    # Read-PrefetchControlOffsets 
      #endif
      #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_OutputSlidePrefetchUSE
        mov 8(%rax, %rbx, 2), %rdx     # Write-PrefetchControlOffsets      
      #endif      
      
      mov 32(%rbp), %rax
      mov %rax, 184(%rbp)       # Loop-Counter
      MainLoop2:
        mov 40(%rbp), %rax 
        mov %rax, 168(%rbp)    # Inner Loop Counter
        MainLoop3:
	  movapd (%rsi), %xmm0
	  movapd (%r11), %xmm1
	  movapd (%rsi, %r9), %xmm2
	  movapd (%r11, %r9), %xmm3
	  movapd (%rsi, %r8), %xmm4
	  movapd (%r11, %r8), %xmm5
	  movapd (%rsi, %r10), %xmm6
	  movapd (%r11, %r10), %xmm7

          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_OutputSlidePrefetchUSE
            mov (%rdx), %rax
            movsx %eax, %rbx
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_OutputSlidePrefetchCMD (%rdi, %rbx)
            shr $32, %rax
          #endif

          //First FFT-step
          subpd %xmm1, %xmm0
          subpd %xmm3, %xmm2
          subpd %xmm5, %xmm4
          subpd %xmm7, %xmm6
	  
	  addpd (%rsi), %xmm1
	  addpd (%rsi, %r9), %xmm3
	  addpd (%rsi, %r8), %xmm5
	  addpd (%rsi, %r10), %xmm7
	  
          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_OutputSlidePrefetchUSE
            movsx %eax, %rbx
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_OutputSlidePrefetchCMD (%rdi, %rbx)
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

          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_InputSlidePrefetchUSE
            mov (%rcx), %rax
            movsx %eax, %rbx
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_InputSlidePrefetchCMD (%rsi, %rbx)
            shr $32, %rax
          #endif
	  
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

          #ifdef xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_InputSlidePrefetchUSE
            movsx %eax, %rbx
            xSSE_ExtremeFFT4D_Read2DSliceAndDoFirstSecondAndThird1DD2FFTStepBase2_InputSlidePrefetchCMD (%rsi, %rbx)
          #endif
	  
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
	  
	  movapd %xmm7, (%rdi)
	  movapd %xmm11, (%rdi, %r12)
	  movapd %xmm10, (%rdi, %r13)
	  movapd %xmm6, (%rdi, %r14)
	  movapd %xmm3, (%r15)
	  movapd %xmm8, (%r15, %r12)
	  movapd %xmm1, (%r15, %r13)
	  movapd %xmm0, (%r15, %r14)
	  
	  add $16, %rsi
	  add $16, %r11
	  add $16, %rdi
          add $16, %r15
	  add $8, %rcx
	  add $8, %rdx
          subq $1, 168(%rbp)
        jge MainLoop3
	
	add 224(%rbp), %rsi
	add 224(%rbp), %r11
	
	
        //Increment Write-Addresses
	mov 184(%rbp), %rax
	mov 56(%rbp), %rbx
	add (%rax, %rbx), %rdi
	add (%rax, %rbx), %r15

        subq $64, 184(%rbp)
      jge MainLoop2
      
      add 232(%rbp), %rsi
      add 232(%rbp), %r11
      
      subq $8, 176(%rbp)
    jge MainLoop1


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
    //Set Error to 1         # Size1 <= 0
    mov $1, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error2:
    //Set Error to 2         # Size1 not multiple of 1
    mov $2, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error3:                    # Size2 <= 0
    //Set Error to 3 
    mov $3, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error4:                    # Size2 not multiple of 8
    //Set Error to 4 
    mov $4, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
  Error5:                    # localIndexCount <= 0
    //Set Error to 5 
    mov $5, %rax
    mov %rax, 112(%rbp)
    jmp ClearStack
