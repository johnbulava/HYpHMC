/**  void xSSE_ComplexVectorAddition(char* parameters);
*    -->Parameters:                                  Offset       Comment
*          char* workSpace                           0
*          Complex* input                            8
*          Complex* output                           16
*          long int blockSize                        24
*          Complex alpha                             32
*          long int L0                               48
*          long int L1                               56
*          long int L2                               64
*          long int L3                               72
*          long int xtrAdd1                          80
*          long int xtrAdd2                          88
*          long int xtrAdd3                          96
*          long int L3LineCombineCount               104
*          long int PrefetchReadAheadCacheLines      112
*          long int L1CacheSizeInBytesPerWay         120
*          long int error                            128
*
*          long int TEMPORAL                         136           Number of 8-blocks = L3*blockSize/8
*          long int TEMPORAL                         144           Number of rest Complex numbers = L3*blockSize%8
*          long int TEMPORAL                         176           L0-Loop Counter
*          long int TEMPORAL                         184           L1-Loop Counter
*          long int TEMPORAL                         200           L3-Combined block Loop Counter
*          long int TEMPORAL                         208           L3-Combined block count = L2 / L3LineCombineCount
*          long int TEMPORAL                         216           Prefetch-Ahead in Bytes = 64*PrefetchReadAheadCacheLines
*          long int TEMPORAL                         224           Last-Vector-Index
*          long int TEMPORAL                         232           Prefetch-Extra-Add1 = xtrAdd1 + Prefetch Prefetch-Extra-Add2
*          long int TEMPORAL                         240           Prefetch-Extra-Add2 = xtrAdd2 + Prefetch Prefetch-Extra-Add3
*          long int TEMPORAL                         248           Prefetch-Extra-Add3 = xtrAdd3 + if(L3*blockSize%4==0){64}else{16*L3*blockSize%4}
*          long int TEMPORAL                         256           Prefetch-Extra-Add-Counter = L3*blockSize/4 - PrefetchReadAheadCacheLines + if(L3*blockSize%4==0){0}else{1}
*          Complex                                   272           Bit-Mask for (-1, 1)
*
*
*
*    -->Target: 64-Bit AMD64 architecture with GNU g++
*
*    -->Workspace Demand: L3*L3LineCombineCount*blockSize
*
*    -->Description:
*          This routine performs the addition of two vectors of complex 
*          numbers (double, double) given in a lattice-structure (L0,L1,L2,L3) 
*          in the input-Vector and writes the result
*          to output. There are blockSize complex numbers per
*          lattice site. XtrAdd? gives additional increments
*          according to a possible embedded lattice-structure.
*          L3LineCombineCount is the number of L3-Lattice lines
*          that are read in one stride. PrefetchReadAheadCacheLines is
*          the number of cache-lines to prefetch data ahead.
*          L1CacheSizeInBytesPerWay is the number of bytes in the L1-Cache
*          per way. Error is set to 0 on success or to an error code if
*          the routine failed.
**/
.file	"xSSE_ComplexVectorAddition.S"
.text
.align 2
.globl xSSE_ComplexVectorAddition
.type	xSSE_ComplexVectorAddition, @function
xSSE_ComplexVectorAddition:
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
    //Preparing Constants for Main loop 
    mov 24(%rbp), %rax     # blockSize
    cmp $0, %rax
    jle Error1
    mov 48(%rbp), %rax     # L0
    cmp $0, %rax
    jle Error2
    mov 56(%rbp), %rax     # L1
    cmp $0, %rax
    jle Error2
    mov 64(%rbp), %rax     # L2
    cmp $0, %rax
    jle Error2
    mov 72(%rbp), %rax     # L3
    cmp $0, %rax
    jle Error2
    mov 24(%rbp), %rbx     # blockSize
    imul %rbx, %rax
    cmp $8, %rax
    jl Error3
    mov %rax, %rcx
    shr $2, %rax
    mov %rax, %rdx
    sub 112(%rbp), %rdx    
    mov %rdx, 256(%rbp)    
    shr $1, %rax
    mov %rax, 136(%rbp)    # Number of 8-blocks = L3*blockSize/8
    and $7, %rcx
    mov %rcx, %r8
    mov %rcx, 144(%rbp)    # Number of rest Complex numbers = L3*blockSize%8
    shl $4, %r8    
    mov $1, %rax
    mov $0, %r15
    mov $64, %r14
    test $3, %rcx
    cmove %r14, %r8         # if(L3*blockSize%4==0){64}else{16*L3*blockSize%4}
    test $3, %rcx    
    cmove %r15, %rax         # if(L3*blockSize%4==0){0}else{1}
    add %rax, 256(%rbp)    # Prefetch-Extra-Add-Counter = L3*blockSize/4 - PrefetchReadAheadCacheLines + if(L3*blockSize%4==0){0}else{1}   
    
    
    mov 96(%rbp), %rax     # xtrAdd3
    cmp $0, %rax
    jl Error4
    add %rax, %r8
    mov %r8,248(%rbp)      # Prefetch-Extra-Add3 = xtrAdd3 + if(L3*blockSize%4==0){64}else{16*L3*blockSize%4}
    mov 88(%rbp), %rax     # xtrAdd2
    cmp $0, %rax
    jl Error4
    add %rax, %r8
    mov %r8, 240(%rbp)     # Prefetch-Extra-Add2 = xtrAdd2 + Prefetch Prefetch-Extra-Add3
    mov 80(%rbp), %rax     # xtrAdd1
    cmp $0, %rax
    jl Error4
    add %rax, %r8
    mov %r8, 232(%rbp)     # Prefetch-Extra-Add1 = xtrAdd1 + Prefetch Prefetch-Extra-Add2
    mov 104(%rbp), %rbx    # L3LineCombineCount
    cmp $0, %rbx
    jle Error5    
    mov 64(%rbp), %rax     # L2
    mov $0, %rdx
    div %rbx
    cmp $0, %rdx
    jne Error6             # Check L2 multiple of L3LineCombineCount
    mov %rax, 208(%rbp)    # L3-Combined block count
    mov 112(%rbp), %rax    # PrefetchReadAheadCacheLines
    cmp $0, %rax
    jle Error7
    shl $6, %rax
    mov %rax, 216(%rbp)
    mov 120(%rbp), %rax    # L1CacheSizeInBytesPerWay
    cmp $0, %rax
    jle Error8
    mov $0, %rax
    mov %rax, 224(%rbp)
    mov $0x8000000000000000, %rax
    mov %rax, 272(%rbp)
    mov $0, %rax
    mov %rax, 280(%rbp)
    
    movapd 32(%rbp), %xmm14
    movhlps %xmm14, %xmm15
    movlhps %xmm14, %xmm14
    movlhps %xmm15, %xmm15
    xorpd 272(%rbp), %xmm15

    mov 248(%rbp), %r12     # Prefetch-Extra-Add3 = xtrAdd3 + if(L3*blockSize%4==0){64}else{16*L3*blockSize%4}
    mov 256(%rbp),%r14      # Prefetch-Extra-Add-Counter = L3*blockSize/4 - PrefetchReadAheadCacheLines + if(L3*blockSize%4==0){0}else{1}
    mov 136(%rbp),%r10      # Number of 8-blocks = L3*blockSize/8
    mov 144(%rbp),%rdx      # Number of rest Complex numbers = L3*blockSize%8	    
    mov 16(%rbp), %rdi      # Output
    mov 8(%rbp), %rsi       # Input


    //Main-Loop    
    mov 48(%rbp), %rax      # L0
    mov %rax, 176(%rbp)    
    loopL0:
      mov 56(%rbp), %rax      # L1
      mov %rax, 184(%rbp)    
      loopL1:
        mov 208(%rbp), %rax     # L3-Combined block count
        mov %rax, 200(%rbp)  
        loopL2:
          mov $0, %rax
	  sub 216(%rbp), %rax   # Prefetch-Ahead in Bytes = 64*PrefetchReadAheadCacheLines
          mov 232(%rbp), %rbx   # Prefetch-Extra-Add1 = xtrAdd1 + Prefetch Prefetch-Extra-Add2
          mov 240(%rbp), %rcx   # Prefetch-Extra-Add2 = xtrAdd2 + Prefetch Prefetch-Extra-Add3
          mov 248(%rbp), %r11   # Prefetch-Extra-Add3 = xtrAdd3 + if(L3*blockSize%4==0){64}else{16*L3*blockSize%4}
          cmp $1, 176(%rbp)
	  cmove %rax, %rbx
          cmp $1, 184(%rbp)
	  cmove %rbx, %rcx
          cmp $1, 200(%rbp)
	  cmove %rcx, %r11
	  sub 8(%rbp), %rsi
	  mov %rsi, 224(%rbp)    # Last Vector Index


          //Read input and multiply with alpha
          mov 0(%rbp), %rdi      # Work-Space
          mov 8(%rbp), %rsi      # Input
          mov 224(%rbp), %rax    # Last Vector Index
          add %rax, %rsi
          mov 216(%rbp),%r15      # Prefetch-Ahead in Bytes = 64*PrefetchReadAheadCacheLines
          mov 104(%rbp),%r8       # L3LineCombineCount
          loopL3CombiningA:
            mov %r14, %r13
            mov %r10, %r9
		
	    loopL3BlocksA:
	      movapd 0(%rsi), %xmm0
	      
	      prefetchnta 256(%rsi)
	      prefetchnta 320(%rsi)
	      
	      
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      nop
	      
	      
	      
//	      movapd 16(%rsi), %xmm1
//	      movapd 32(%rsi), %xmm2
//	      movapd 48(%rsi), %xmm3
	      movapd 64(%rsi), %xmm4
//	      movapd 80(%rsi), %xmm5
//	      movapd 96(%rsi), %xmm6
//	      movapd 112(%rsi), %xmm7	      
		  
/*	      movhlps %xmm0, %xmm8  
	      movlhps %xmm0, %xmm8  
	      movhlps %xmm1, %xmm9  
	      movlhps %xmm1, %xmm9  
	      movhlps %xmm2, %xmm10  
	      movlhps %xmm2, %xmm10  
	      movhlps %xmm3, %xmm11  
	      movlhps %xmm3, %xmm11  	      
	      
	      mulpd %xmm14, %xmm0
	      mulpd %xmm15, %xmm8
	      mulpd %xmm14, %xmm1
	      mulpd %xmm15, %xmm9
	      mulpd %xmm14, %xmm2
	      mulpd %xmm15, %xmm10
	      mulpd %xmm14, %xmm3
	      mulpd %xmm15, %xmm11
	      
	      addpd %xmm8, %xmm0
	      addpd %xmm9, %xmm1
	      addpd %xmm10, %xmm2
	      addpd %xmm11, %xmm3
	      	      
	      movapd %xmm0, 0(%rdi)
	      movapd %xmm1, 16(%rdi)
	      movapd %xmm2, 32(%rdi)
	      movapd %xmm3, 48(%rdi)	      
	      
	      movhlps %xmm4, %xmm12  
	      movlhps %xmm4, %xmm12  
	      movhlps %xmm5, %xmm13  
	      movlhps %xmm5, %xmm13  
	      movhlps %xmm6, %xmm8  
	      movlhps %xmm6, %xmm8  
	      movhlps %xmm7, %xmm9  
	      movlhps %xmm7, %xmm9  		  

	      mulpd %xmm14, %xmm4
	      mulpd %xmm15, %xmm12
	      mulpd %xmm14, %xmm5
	      mulpd %xmm15, %xmm13
	      mulpd %xmm14, %xmm6
	      mulpd %xmm15, %xmm8
	      mulpd %xmm14, %xmm7
	      mulpd %xmm15, %xmm9
	      
	      addpd %xmm12, %xmm4
	      addpd %xmm13, %xmm5
	      addpd %xmm8, %xmm6
	      addpd %xmm9, %xmm7
	      	      
	      movapd %xmm4, 64(%rdi)
	      movapd %xmm5, 80(%rdi)
	      movapd %xmm6, 96(%rdi)
	      movapd %xmm7, 112(%rdi) */

/*              prefetchnta (%rsi, %r15)
              xor %rax, %rax
	      cmp $1, %r8
	      cmove %r11, %rax
	      dec %r13
	      cmove %r12, %rax
	      add %rax, %r15

              prefetchnta 64(%rsi, %r15)
              xor %rax, %rax
	      cmp $1, %r8
	      cmove %r11, %rax
	      dec %r13
	      cmove %r12, %rax
	      add %rax, %r15*/

	      add $128, %rsi
    	      add $128, %rdi
	      dec %r9
	    jg loopL3BlocksA
		
	    cmp $0, %rdx
	    jle NoExtraL3A
	      mov %rdx, %rcx
	      cmp $4, %rdx
	      jle LessOrEqualFourExtraL3A
	        movapd 1(%rsi), %xmm0
	        movapd 16(%rsi), %xmm1
	        movapd 32(%rsi), %xmm2
	        movapd 48(%rsi), %xmm3
		  
	        movhlps %xmm0, %xmm8  
	        movlhps %xmm0, %xmm8  
	        movhlps %xmm1, %xmm9  
	        movlhps %xmm1, %xmm9  
	        movhlps %xmm2, %xmm10  
	        movlhps %xmm2, %xmm10  
	        movhlps %xmm3, %xmm11  
	        movlhps %xmm3, %xmm11  
	      
	        mulpd %xmm14, %xmm0
	        mulpd %xmm15, %xmm8
	        mulpd %xmm14, %xmm1
	        mulpd %xmm15, %xmm9
	        mulpd %xmm14, %xmm2
	        mulpd %xmm15, %xmm10
	        mulpd %xmm14, %xmm3
	        mulpd %xmm15, %xmm11
	      
	        addpd %xmm8, %xmm0
	        addpd %xmm9, %xmm1
	        addpd %xmm10, %xmm2
	        addpd %xmm11, %xmm3
	      	      
	        movapd %xmm0, 0(%rdi)
	        movapd %xmm1, 16(%rdi)
	        movapd %xmm2, 32(%rdi)
	        movapd %xmm3, 48(%rdi)	      

                prefetchnta (%rsi, %r15)
                xor %rax, %rax
	        cmp $1, %r8
	        cmove %r11, %rax
	        dec %r13
	        cmove %r12, %rax
	        add %rax, %r15

	        add $64, %rsi
    	        add $64, %rdi
		sub $4, %rcx
	      LessOrEqualFourExtraL3A:
           
              prefetchnta (%rsi, %r15)
              xor %rax, %rax
	      cmp $1, %r8
	      cmove %r11, %rax
	      dec %r13
	      cmove %r12, %rax
	      add %rax, %r15
	      
	      ExtraL3InnerLoopA:
	        movapd 0(%rsi), %xmm0
	        movhlps %xmm0, %xmm8  
	        movlhps %xmm0, %xmm8  
	        mulpd %xmm14, %xmm0
	        mulpd %xmm15, %xmm8
	        addpd %xmm8, %xmm0
	        movapd %xmm0, 0(%rdi)
	      
	        add $16, %rsi
    	        add $16, %rdi              
                dec %rcx
	      jg ExtraL3InnerLoopA
	    NoExtraL3A:
	      
	    dec %r8
          jg loopL3CombiningA
	  
	  
	  
	  
	  
	  


          add 96(%rbp), %rsi
          sub $1, 200(%rbp)
	jg loopL2
        add 88(%rbp), %rsi
        sub $1, 184(%rbp)
      jg loopL1
      add 80(%rbp), %rsi
      sub $1, 176(%rbp)
    jg loopL0
    
	 

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
  Error1:                       # blockSize <= 0
    //Set Error to 1
    mov $1, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack    
  Error2:                       # L0, L1, L2, or L3 <= 0
    //Set Error to 2 
    mov $2, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack    
  Error3:                       # L3*blockSize <= 0
    //Set Error to 3
    mov $3, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack    
  Error4:                       # xtrAdd1, xtrAdd2, xtrAdd3 < 0
    //Set Error to 4
    mov $4, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack    
  Error5:                       # L3LineCombineCount <= 0
    //Set Error to 5
    mov $5, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack    
  Error6:                       # L2 not multiple of L3LineCombineCount
    //Set Error to 6
    mov $6, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack
  Error7:                       # PrefetchReadAheadCacheLines <= 0
    //Set Error to 7
    mov $7, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack
  Error8:                       # PrefetchReadAheadCacheLines <= 0
    //Set Error to 8 
    mov $8, %rax
    mov %rax, 128(%rbp)
    jmp ClearStack
