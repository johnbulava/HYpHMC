                  "mov 104(%6), %7\n\t"
                  "movapd 208(%6), %%xmm15\n\t"
  
                  "movapd 512(%0), %%xmm1\n\t"
                  "movhlps %%xmm1, %%xmm0\n\t"
                  "movlhps %%xmm1, %%xmm0\n\t"
                  "mulpd (%7,%5), %%xmm1\n\t"
                  "mulpd 16(%7,%5), %%xmm0\n\t"
                  "addpd %%xmm0, %%xmm1\n\t"
		  
                  "movapd (%0,%11), %%xmm2\n\t"
                  "movhlps %%xmm2, %%xmm8\n\t"
                  "movlhps %%xmm2, %%xmm8\n\t"
                  "mulpd (%7,%4), %%xmm2\n\t"
                  "mulpd 16(%7,%4), %%xmm8\n\t"
                  "addpd %%xmm8, %%xmm2\n\t"
		  
                  "movapd 512(%0,%11), %%xmm3\n\t"
                  "add %4, %7\n\t"
                  "movhlps %%xmm3, %%xmm9\n\t"
                  "movlhps %%xmm3, %%xmm9\n\t"
                  "mulpd (%7,%5), %%xmm3\n\t"
                  "mulpd 16(%7,%5), %%xmm9\n\t"
                  "addpd %%xmm9, %%xmm3\n\t"

                  "movapd (%0,%12), %%xmm4\n\t"
                  "sub %4, %7\n\t"
                  "movhlps %%xmm4, %%xmm10\n\t"
                  "movlhps %%xmm4, %%xmm10\n\t"
                  "mulpd (%7,%3), %%xmm4\n\t"
                  "mulpd 16(%7,%3), %%xmm10\n\t"
                  "addpd %%xmm10, %%xmm4\n\t"

                  "addpd %%xmm2, %%xmm3\n\t"
                  "mulpd %%xmm15, %%xmm2\n\t"
                  "subpd %%xmm3, %%xmm2\n\t"
		  
                  "movapd 512(%0,%12), %%xmm5\n\t"
                  "add %3, %7\n\t"
                  "movhlps %%xmm5, %%xmm11\n\t"
                  "movlhps %%xmm5, %%xmm11\n\t"
                  "mulpd (%7,%5), %%xmm5\n\t"
                  "mulpd 16(%7,%5), %%xmm11\n\t"
                  "addpd %%xmm11, %%xmm5\n\t"
		  
                  "movapd (%0,%13), %%xmm6\n\t"
                  "movhlps %%xmm6, %%xmm12\n\t"
                  "movlhps %%xmm6, %%xmm12\n\t"
                  "mulpd (%7,%4), %%xmm6\n\t"
                  "mulpd 16(%7,%4), %%xmm12\n\t"
                  "addpd %%xmm12, %%xmm6\n\t"

                  "addpd %%xmm4, %%xmm5\n\t"
                  "mulpd %%xmm15, %%xmm4\n\t"
                  "subpd %%xmm5, %%xmm4\n\t"
		  
                  "movapd 512(%0,%13), %%xmm7\n\t"
                  "add %4, %7\n\t"
                  "movhlps %%xmm7, %%xmm13\n\t"
                  "movlhps %%xmm7, %%xmm13\n\t"
                  "mulpd (%7,%5), %%xmm7\n\t"
                  "mulpd 16(%7,%5), %%xmm13\n\t"
                  "addpd %%xmm13, %%xmm7\n\t"
		  
                  "movapd (%14), %%xmm8\n\t"
                  "sub %3, %7\n\t"
                  "sub %4, %7\n\t"
                  "movhlps %%xmm8, %%xmm14\n\t"
                  "movlhps %%xmm8, %%xmm14\n\t"
                  "mulpd (%7,%2), %%xmm8\n\t"
                  "mulpd 16(%7,%2), %%xmm14\n\t"
                  "addpd %%xmm14, %%xmm8\n\t"

                  "addpd %%xmm6, %%xmm7\n\t"
                  "mulpd %%xmm15, %%xmm6\n\t"
                  "subpd %%xmm7, %%xmm6\n\t"
		  
                  "movapd 512(%14), %%xmm9\n\t"
                  "add %2, %7\n\t"
                  "movhlps %%xmm9, %%xmm0\n\t"
                  "movlhps %%xmm9, %%xmm0\n\t"
                  "mulpd (%7,%5), %%xmm9\n\t"
                  "mulpd 16(%7,%5), %%xmm0\n\t"
                  "addpd %%xmm0, %%xmm9\n\t"

                  "movapd (%14,%11), %%xmm10\n\t"
                  "movhlps %%xmm10, %%xmm13\n\t"
                  "movlhps %%xmm10, %%xmm13\n\t"
                  "mulpd (%7,%4), %%xmm10\n\t"
                  "mulpd 16(%7,%4), %%xmm13\n\t"
                  "addpd %%xmm13, %%xmm10\n\t"

                  "addpd %%xmm8, %%xmm9\n\t"
                  "mulpd %%xmm15, %%xmm8\n\t"
                  "subpd %%xmm9, %%xmm8\n\t"
		  
                  "movapd 512(%14,%11), %%xmm11\n\t"
                  "add %4, %7\n\t"
                  "movhlps %%xmm11, %%xmm14\n\t"
                  "movlhps %%xmm11, %%xmm14\n\t"
                  "mulpd (%7,%5), %%xmm11\n\t"
                  "mulpd 16(%7,%5), %%xmm14\n\t"
                  "addpd %%xmm14, %%xmm11\n\t"
		  
                  "movapd (%14,%12), %%xmm12\n\t"
                  "sub %4, %7\n\t"
                  "movhlps %%xmm12, %%xmm0\n\t"
                  "movlhps %%xmm12, %%xmm0\n\t"
                  "mulpd (%7,%3), %%xmm12\n\t"
                  "mulpd 16(%7,%3), %%xmm0\n\t"
                  "addpd %%xmm0, %%xmm12\n\t"

                  "addpd %%xmm10, %%xmm11\n\t"
                  "mulpd %%xmm15, %%xmm10\n\t"
                  "subpd %%xmm11, %%xmm10\n\t"
		  
                  "movapd 512(%14,%12), %%xmm13\n\t"
                  "add %3, %7\n\t"
                  "movhlps %%xmm13, %%xmm0\n\t"
                  "movlhps %%xmm13, %%xmm0\n\t"
                  "mulpd (%7,%5), %%xmm13\n\t"
                  "mulpd 16(%7,%5), %%xmm0\n\t"
                  "addpd %%xmm0, %%xmm13\n\t"
		  
                  "movapd (%14,%13), %%xmm14\n\t"
                  "movhlps %%xmm14, %%xmm0\n\t"
                  "movlhps %%xmm14, %%xmm0\n\t"
                  "mulpd (%7,%4), %%xmm14\n\t"
                  "mulpd 16(%7,%4), %%xmm0\n\t"
                  "addpd %%xmm0, %%xmm14\n\t"

                  "addpd %%xmm12, %%xmm13\n\t"
                  "mulpd %%xmm15, %%xmm12\n\t"
                  "subpd %%xmm13, %%xmm12\n\t"
		  
                  "movapd 512(%14,%13), %%xmm15\n\t"
                  "add %4, %7\n\t"
                  "movhlps %%xmm15, %%xmm0\n\t"
                  "movlhps %%xmm15, %%xmm0\n\t"
                  "mulpd (%7,%5), %%xmm15\n\t"
                  "mulpd 16(%7,%5), %%xmm0\n\t"
                  "addpd %%xmm0, %%xmm15\n\t"
	    
                  "movapd (%0), %%xmm0\n\t"


  	          //Rest of First STEP
                  "addpd %%xmm0, %%xmm1\n\t"
                  "mulpd 208(%6), %%xmm0\n\t"
                  "subpd %%xmm1, %%xmm0\n\t"
                  "addpd %%xmm14, %%xmm15\n\t"
                  "mulpd 208(%6), %%xmm14\n\t"
                  "subpd %%xmm15, %%xmm14\n\t"

