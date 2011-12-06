#include "SSEroutines.h"





/**** SSE_FourierTrafoForthFFTstep ***
* Performs the third FFT step. 
* ATTENTION: Restriction to oneDim = 128. Input and Output X3-Increments
* are fixed. Input and Output are exprected to have same size.
* Input and Output must not be the same!!!
*
* Parameters: The following data are expected in ASMParameterData
* OneDim                : Lattice size in one dimension
* data                  : Input-/Output-Vector
* Increments           : Increments for output base and output pi-shifted
* ASMParameterData[2576...] : Relative Lese-, Schreibe-Adressen
*
**/
void SSE_FourierTrafoForthFFTstep() {
//inline void SSE_FourierTrafoForthFFTstep() {
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;
  int* Dummy11;
  int* Dummy12;
  int* Dummy13;
  int* Dummy14;
  int* Dummy15;  
  long int* para_int = (long int*) ASMParameterData;
  Complex* para_Comp = ASMParameterData;
  
  
  
__asm__ volatile(
        "mov %15, %6\n\t"
        "mov 16(%6), %0\n\t"
        "mov 24(%6), %1\n\t"
        "mov 16(%6), %8\n\t"
	
        "mov $0, %2\n\t"
        "mov $0, %3\n\t"
        "mov $0, %4\n\t"
        "mov $0, %5\n\t"
        "mov %2, 40(%6)\n\t"
        "mov %3, 48(%6)\n\t"
        "mov %4, 56(%6)\n\t"
	
        "mov %5, 64(%6)\n\t"
        "mov %5, 72(%6)\n\t"
        "mov %5, 80(%6)\n\t"
	


        "mov 288(%6), %11\n\t"
        "mov 296(%6), %12\n\t"
        "mov 304(%6), %13\n\t"
        "mov 312(%6), %14\n\t"

        "mov $16384, %10\n\t"
        "mov %10, 144(%6)\n\t"
        "mov 136(%6), %10\n\t"
        "add %10, 144(%6)\n\t"
        "add %10, 144(%6)\n\t"


        "SSE_FourierTrafoForthFFTstep_MainLoop1: \n\t"
          "SSE_FourierTrafoForthFFTstep_MainLoop2: \n\t"
            "SSE_FourierTrafoForthFFTstep_MainLoop3: \n\t"
              "mov $0, %9\n\t"
	      
	      //Prefetch
              "SSE_FourierTrafoForthFFTstep_schleife4READ: \n\t"
                "SSE_FourierTrafoForthFFTstep_InnerLoopREAD: \n\t"

                  "add %0, %14\n\t"
	
                  // LOAD + Multiply with Weights
                  #include "ExtremeFFTasmBlock3.h"
	    
                  "sub %0, %14\n\t"
	    
	    
  	          //Prefetching
                  "mov 2560(%6,%9,8), %7\n\t"		  
                  "prefetch (%8,%7)\n\t"	
                  "prefetch 128(%8,%7)\n\t"	
                  "prefetch 256(%8,%7)\n\t"	
                  "prefetch 384(%8,%7)\n\t"	
                  "prefetch 512(%8,%7)\n\t"	
                  "prefetch 640(%8,%7)\n\t"	
                  "prefetch 768(%8,%7)\n\t"	
                  "prefetch 896(%8,%7)\n\t"	
  	          //Prefetching END



	          // Second + third Step are external
                  #include "ExtremeFFTasmBlock.h"

                  "mov %9, %7\n\t"
                  "and $1984, %7\n\t"
                  "add %7, %9\n\t"
		  
                  "add 144(%6), %9\n\t"
                  "sub %10, %9\n\t"		  

                  // Forth STEP + SAVE
                  #include "ExtremeFFTasmBlock4.h"


                  "sub 144(%6), %9\n\t"
                  "add %10, %9\n\t"		  
                  "sub %7, %9\n\t"


                  "add $16, %0\n\t"
                  "add $16, %9\n\t"
                  "mov %9, %7\n\t"
                  "and $63, %7\n\t"
                "jne SSE_FourierTrafoForthFFTstep_InnerLoopREAD\n\t"
                "sub $64, %0\n\t"


                "add $128, %0\n\t"
                "add $32, %5\n\t"
                "cmp $128,%5\n\t"
              "jl SSE_FourierTrafoForthFFTstep_schleife4READ\n\t"
	      
	      
	      
	      //Write Buffer
              "SSE_FourierTrafoForthFFTstep_schleife4WRITE: \n\t"
                "SSE_FourierTrafoForthFFTstep_InnerLoopWRITE: \n\t"
	
                  // LOAD + Multiply with Weights
                  "add %0, %14\n\t"
                  #include "ExtremeFFTasmBlock3.h"
                  "sub %0, %14\n\t"


	          // Second + third Step are external
                  #include "ExtremeFFTasmBlock.h"


                  // Forth STEP + SAVE
                  "mov %9, %7\n\t"
                  "and $1984, %7\n\t"
                  "add %7, %9\n\t"
                  "add 144(%6), %9\n\t"
                  "sub %10, %9\n\t"		  
                  #include "ExtremeFFTasmBlock4.h"
                  "sub 144(%6), %9\n\t"
                  "add %10, %9\n\t"		  
                  "sub %7, %9\n\t"


                  //Write Buffer - Part 1
                  "sub $256, %9\n\t"		  
                  "shl $6, %9\n\t"		  
                  "add %9, %10\n\t"		  
                  "shr $3, %9\n\t"	
		  
                  "mov 2560(%6,%9), %7\n\t"		  
                  "movapd 18944(%6,%10), %%xmm0\n\t"		  
                  "movapd 18960(%6,%10), %%xmm1\n\t"		  
                  "movapd 18976(%6,%10), %%xmm2\n\t"		  
                  "movapd 18992(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"		  
		  
                  "mov 2576(%6,%9), %7\n\t"		  
                  "movapd 19072(%6,%10), %%xmm0\n\t"		  
                  "movapd 19088(%6,%10), %%xmm1\n\t"		  
                  "movapd 19104(%6,%10), %%xmm2\n\t"		  
                  "movapd 19120(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"		  

		  
                  "mov 2592(%6,%9), %7\n\t"		  
                  "movapd 19200(%6,%10), %%xmm0\n\t"		  
                  "movapd 19216(%6,%10), %%xmm1\n\t"		  
                  "movapd 19232(%6,%10), %%xmm2\n\t"		  
                  "movapd 19248(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"	

		  
                  "mov 2608(%6,%9), %7\n\t"		  
                  "movapd 19328(%6,%10), %%xmm0\n\t"		  
                  "movapd 19344(%6,%10), %%xmm1\n\t"		  
                  "movapd 19360(%6,%10), %%xmm2\n\t"		  
                  "movapd 19376(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"		  


                  "mov 2624(%6,%9), %7\n\t"		  
                  "movapd 19456(%6,%10), %%xmm0\n\t"		  
                  "movapd 19472(%6,%10), %%xmm1\n\t"		  
                  "movapd 19488(%6,%10), %%xmm2\n\t"		  
                  "movapd 19504(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"		  

                  "mov 2640(%6,%9), %7\n\t"		  
                  "movapd 19584(%6,%10), %%xmm0\n\t"		  
                  "movapd 19600(%6,%10), %%xmm1\n\t"		  
                  "movapd 19616(%6,%10), %%xmm2\n\t"		  
                  "movapd 19632(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"		  

                  "mov 2656(%6,%9), %7\n\t"		  
                  "movapd 19712(%6,%10), %%xmm0\n\t"		  
                  "movapd 19728(%6,%10), %%xmm1\n\t"		  
                  "movapd 19744(%6,%10), %%xmm2\n\t"		  
                  "movapd 19760(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"		  

                  "mov 2672(%6,%9), %7\n\t"		  
                  "movapd 19840(%6,%10), %%xmm0\n\t"		  
                  "movapd 19856(%6,%10), %%xmm1\n\t"		  
                  "movapd 19872(%6,%10), %%xmm2\n\t"		  
                  "movapd 19888(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"		  
		  
		  
                  "shl $3, %9\n\t"		  
                  "sub %9, %10\n\t"		  
                  "shr $6, %9\n\t"		  
                  "add $256, %9\n\t"		  
                  //Write Buffer END - Part 1




                  "add $16, %0\n\t"
                  "add $16, %9\n\t"
                  "mov %9, %7\n\t"
                  "and $63, %7\n\t"
                "jne SSE_FourierTrafoForthFFTstep_InnerLoopWRITE\n\t"
                "sub $64, %0\n\t"


                "add $128, %0\n\t"
                "add $32, %5\n\t"
                "cmp 0(%6),%5\n\t"
              "jl SSE_FourierTrafoForthFFTstep_schleife4WRITE\n\t"	     
              "mov $0, %5\n\t"


              //New Write-Buffer-Index
              "mov 144(%6), %9\n\t"
              "sub %10, %9\n\t"
              "mov %9, %10\n\t"

              //New Output-Base
              "mov 72(%6), %9\n\t"
              "mov %9, 80(%6)\n\t"
              "mov %9, %1\n\t"
              "add 24(%6), %1\n\t"
	      
              //New Input-Base
              "mov 64(%6), %9\n\t"
              "mov %9, 72(%6)\n\t"
              "mov %9, %0\n\t"
              "add 16(%6), %0\n\t"
	      
              //New Prefetch-Base
              "add 128(%6), %9\n\t"
              "mov %9, 64(%6)\n\t"
              "mov %9, %8\n\t"
              "add 16(%6), %8\n\t"

              //New Counter - Values
              "mov 56(%6), %4\n\t"
              "mov 48(%6), %3\n\t"
              "mov 40(%6), %2\n\t"
	      
              "add $32, 56(%6)\n\t"
              "mov 0(%6), %9\n\t"	      
              "cmp %9,56(%6)\n\t"
            "jl SSE_FourierTrafoForthFFTstep_MainLoop3\n\t"
            "mov $0, %9\n\t"	      
            "mov %9, 56(%6)\n\t"
	    
            //New Prefetch-Base
            "mov 64(%6), %9\n\t"
            "add 120(%6), %9\n\t"
            "mov %9, 64(%6)\n\t"
            "mov %9, %8\n\t"
            "add 16(%6), %8\n\t"

            "add $32, 48(%6)\n\t"
            "mov 0(%6), %9\n\t"	      
            "cmp %9,48(%6)\n\t"
          "jl SSE_FourierTrafoForthFFTstep_MainLoop2\n\t"
          "mov $0, %9\n\t"	      
          "mov %9, 48(%6)\n\t"

          //New Prefetch-Base
          "mov 64(%6), %9\n\t"
          "add 112(%6), %9\n\t"
          "mov %9, 64(%6)\n\t"
          "mov %9, %8\n\t"
          "add 16(%6), %8\n\t"

          "add $32, 40(%6)\n\t"
          "mov 0(%6), %9\n\t"	      
          "cmp %9,40(%6)\n\t"
        "jl SSE_FourierTrafoForthFFTstep_MainLoop1\n\t"  
	


        //Lead Out 1
        "mov $0, %9\n\t"
        "SSE_FourierTrafoForthFFTstep_schleife4LEADOUT1: \n\t"
          "SSE_FourierTrafoForthFFTstep_InnerLoopLEADOUT1: \n\t"
            "add %0, %14\n\t"
	
            // LOAD + Multiply with Weights
            #include "ExtremeFFTasmBlock3.h"
	    
            "sub %0, %14\n\t"
	    
            // Second + third Step are external
            #include "ExtremeFFTasmBlock.h"

            "mov %9, %7\n\t"
            "and $1984, %7\n\t"
            "add %7, %9\n\t"

		  
            "add 144(%6), %9\n\t"
            "sub %10, %9\n\t"		  

            // Forth STEP + SAVE
            #include "ExtremeFFTasmBlock4.h"

            "sub 144(%6), %9\n\t"
            "add %10, %9\n\t"		  
            "sub %7, %9\n\t"


            //Write Buffer
            "shl $5, %9\n\t"		  
            "add %9, %10\n\t"		  
            "shr $3, %9\n\t"		  
		  
            "mov 2560(%6,%9), %7\n\t"		  
            "movapd 18944(%6,%10), %%xmm0\n\t"		  
            "movapd 18960(%6,%10), %%xmm1\n\t"		  
            "movapd 18976(%6,%10), %%xmm2\n\t"		  
            "movapd 18992(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  
		  
            "mov 2576(%6,%9), %7\n\t"		  
            "movapd 19072(%6,%10), %%xmm0\n\t"		  
            "movapd 19088(%6,%10), %%xmm1\n\t"		  
            "movapd 19104(%6,%10), %%xmm2\n\t"		  
            "movapd 19120(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  

            "mov 2592(%6,%9), %7\n\t"		  
            "movapd 19200(%6,%10), %%xmm0\n\t"		  
            "movapd 19216(%6,%10), %%xmm1\n\t"		  
            "movapd 19232(%6,%10), %%xmm2\n\t"		  
            "movapd 19248(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  

            "mov 2608(%6,%9), %7\n\t"		  
            "movapd 19328(%6,%10), %%xmm0\n\t"		  
            "movapd 19344(%6,%10), %%xmm1\n\t"		  
            "movapd 19360(%6,%10), %%xmm2\n\t"		  
            "movapd 19376(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  

            "shl $3, %9\n\t"		  
            "sub %9, %10\n\t"		  
            "shr $5, %9\n\t"		  
            //Write Buffer END

            "add $16, %0\n\t"
            "add $16, %9\n\t"
            "mov %9, %7\n\t"
            "and $63, %7\n\t"
          "jne SSE_FourierTrafoForthFFTstep_InnerLoopLEADOUT1\n\t"
          "sub $64, %0\n\t"


          "add $128, %0\n\t"
          "add $32, %5\n\t"
          "cmp 0(%6),%5\n\t"
        "jl SSE_FourierTrafoForthFFTstep_schleife4LEADOUT1\n\t"	     
        "mov $0, %5\n\t"	
	
	
        //New Write-Buffer-Index
        "mov 144(%6), %9\n\t"
        "sub %10, %9\n\t"
        "mov %9, %10\n\t"

        //New Output-Base
        "mov 72(%6), %9\n\t"
        "mov %9, 80(%6)\n\t"
        "mov %9, %1\n\t"
        "add 24(%6), %1\n\t"
	      


        //Lead Out 2
        "mov $0, %9\n\t"
        "SSE_FourierTrafoForthFFTstep_schleife4LEADOUT2: \n\t"
          "SSE_FourierTrafoForthFFTstep_InnerLoopLEADOUT2: \n\t"
            //Write Buffer
            "shl $5, %9\n\t"		  
            "add %9, %10\n\t"		  
            "shr $3, %9\n\t"		  
		  
            "mov 2560(%6,%9), %7\n\t"		  
            "movapd 18944(%6,%10), %%xmm0\n\t"		  
            "movapd 18960(%6,%10), %%xmm1\n\t"		  
            "movapd 18976(%6,%10), %%xmm2\n\t"		  
            "movapd 18992(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  
		  
            "mov 2576(%6,%9), %7\n\t"		  
            "movapd 19072(%6,%10), %%xmm0\n\t"		  
            "movapd 19088(%6,%10), %%xmm1\n\t"		  
            "movapd 19104(%6,%10), %%xmm2\n\t"		  
            "movapd 19120(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  

            "mov 2592(%6,%9), %7\n\t"		  
            "movapd 19200(%6,%10), %%xmm0\n\t"		  
            "movapd 19216(%6,%10), %%xmm1\n\t"		  
            "movapd 19232(%6,%10), %%xmm2\n\t"		  
            "movapd 19248(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  

            "mov 2608(%6,%9), %7\n\t"		  
            "movapd 19328(%6,%10), %%xmm0\n\t"		  
            "movapd 19344(%6,%10), %%xmm1\n\t"		  
            "movapd 19360(%6,%10), %%xmm2\n\t"		  
            "movapd 19376(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  

            "shl $3, %9\n\t"		  
            "sub %9, %10\n\t"		  
            "shr $5, %9\n\t"		  
            //Write Buffer END


            "add $16, %0\n\t"
            "add $16, %9\n\t"
            "mov %9, %7\n\t"
            "and $63, %7\n\t"
          "jne SSE_FourierTrafoForthFFTstep_InnerLoopLEADOUT2\n\t"
          "sub $64, %0\n\t"


          "add $128, %0\n\t"
          "add $32, %5\n\t"
          "cmp 0(%6),%5\n\t"
        "jl SSE_FourierTrafoForthFFTstep_schleife4LEADOUT2\n\t"	     
        "mov $0, %5\n\t"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10), "=r" (Dummy11),  "=r" (Dummy12), "=r" (Dummy13), "=r" (Dummy14), "=r" (Dummy15)
          : "r" (para_Comp));

/*printf("reg0 %d\n",(long int) Dummy1);
printf("reg1 %d\n",(long int) Dummy2);
printf("reg2 %d\n",(long int) Dummy3);
printf("reg3 %d\n",(long int) Dummy4);
printf("reg4 %d\n",(long int) Dummy5);
printf("reg5 %d\n",(long int) Dummy6);
printf("reg6 %d\n",(long int) Dummy7);
printf("reg7 %d\n",(long int) Dummy8);
printf("reg8 %d\n",(long int) Dummy9);
printf("reg9 %d\n",(long int) Dummy10);
printf("reg10 %d\n",(long int) Dummy11);
printf("reg11 %d\n",(long int) Dummy12);
printf("reg12 %d\n",(long int) Dummy13);
printf("reg13 %d\n",(long int) Dummy14);
printf("reg14 %d\n",(long int) Dummy15);
printf("reg15 %d\n",(long int) para_Comp);

printBits((long int) Dummy12);
printBits((long int) Dummy13);
printBits((long int) Dummy14);
printBits((long int) Dummy15);
exit(0);*/

  //Always false
  if (para_int[2]!=para_int[2]) { 
    printf("%d\n",*Dummy1);
  }
}


/**** SSE_FourierTrafoThirdFFTstep ***
* Performs the third FFT step. 
* ATTENTION: Restriction to oneDim = 128. Input and Output X3-Increments
* are fixed.
*
* Parameters: The following data are expected in ASMParameterData
* OneDim                : Lattice size in one dimension
* data                  : Input-/Output-Vector
* Incerements           : Increments for output base and output pi-shifted
*
**/
void SSE_FourierTrafoThirdFFTstep() {
//inline void SSE_FourierTrafoThirdFFTstep() {
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;
  int* Dummy11;
  int* Dummy12;
  int* Dummy13;
  int* Dummy14;
  int* Dummy15;  
  long int* para_int = (long int*) ASMParameterData;
  Complex* para_Comp = ASMParameterData;
  
  
  
__asm__ volatile(
        "mov %15, %6\n\t"
        "mov 16(%6), %0\n\t"
        "mov 24(%6), %1\n\t"
        "mov 16(%6), %8\n\t"
	
        "mov $0, %2\n\t"
        "mov $0, %3\n\t"
        "mov $0, %4\n\t"
        "mov $0, %5\n\t"
        "mov %2, 40(%6)\n\t"
        "mov %3, 48(%6)\n\t"
        "mov %4, 56(%6)\n\t"
	
        "mov %5, 64(%6)\n\t"
        "mov %5, 72(%6)\n\t"
        "mov %5, 80(%6)\n\t"
	
        "mov %5, 416(%6)\n\t"
        "mov %5, 424(%6)\n\t"
        "mov %5, 432(%6)\n\t"


        "mov $16384, %10\n\t"
        "mov %10, 408(%6)\n\t"
        "mov 400(%6), %10\n\t"
        "add %10, 408(%6)\n\t"
        "add %10, 408(%6)\n\t"

        "mov 288(%6), %11\n\t"
        "mov 296(%6), %12\n\t"
        "mov 304(%6), %13\n\t"
        "mov 312(%6), %14\n\t"

        "SSE_FourierTrafoThirdFFTstep_MainLoop1: \n\t"
          "SSE_FourierTrafoThirdFFTstep_MainLoop2: \n\t"
            "SSE_FourierTrafoThirdFFTstep_MainLoop3: \n\t"
              "mov $0, %9\n\t"
	      
	      //Prefetch
              "SSE_FourierTrafoThirdFFTstep_schleife4: \n\t"
                "SSE_FourierTrafoThirdFFTstep_InnerLoop: \n\t"
	
                  // LOAD + Multiply with Weights
                  "add %0, %14\n\t"
                  #include "ExtremeFFTasmBlock5.h"
                  "sub %0, %14\n\t"
	    
	    
  	          //Prefetching
                  "mov 2560(%6,%9,4), %7\n\t"		  
                  "prefetch (%8,%7)\n\t"	
                  "prefetch 128(%8,%7)\n\t"	
                  "prefetch 256(%8,%7)\n\t"	
                  "prefetch 384(%8,%7)\n\t"	
  	          //Prefetching END


	          // Second + third Step are external
                  #include "ExtremeFFTasmBlock.h"


                  // Forth STEP + SAVE
                  "mov %9, %7\n\t"
                  "and $1984, %7\n\t"
                  "add %7, %9\n\t"
                  "add 408(%6), %9\n\t"
                  "sub %10, %9\n\t"		  
                  #include "ExtremeFFTasmBlock6.h"
                  "sub 408(%6), %9\n\t"
                  "add %10, %9\n\t"		  
                  "sub %7, %9\n\t"


                  //Write Buffer - Part 1
                  "shl $5, %9\n\t"		  
                  "add %9, %10\n\t"		  
                  "shr $3, %9\n\t"	
		  
                  "mov 3584(%6,%9), %7\n\t"		  
                  "movapd 18944(%6,%10), %%xmm0\n\t"		  
                  "movapd 18960(%6,%10), %%xmm1\n\t"		  
                  "movapd 18976(%6,%10), %%xmm2\n\t"		  
                  "movapd 18992(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"		  
		  
                  "mov 3600(%6,%9), %7\n\t"		  
                  "movapd 19072(%6,%10), %%xmm0\n\t"		  
                  "movapd 19088(%6,%10), %%xmm1\n\t"		  
                  "movapd 19104(%6,%10), %%xmm2\n\t"		  
                  "movapd 19120(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"		  
		  
                  "mov 3616(%6,%9), %7\n\t"		  
                  "movapd 19200(%6,%10), %%xmm0\n\t"		  
                  "movapd 19216(%6,%10), %%xmm1\n\t"		  
                  "movapd 19232(%6,%10), %%xmm2\n\t"		  
                  "movapd 19248(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"	
		  
                  "mov 3632(%6,%9), %7\n\t"		  
                  "movapd 19328(%6,%10), %%xmm0\n\t"		  
                  "movapd 19344(%6,%10), %%xmm1\n\t"		  
                  "movapd 19360(%6,%10), %%xmm2\n\t"		  
                  "movapd 19376(%6,%10), %%xmm3\n\t"		  
                  "movntpd %%xmm0, (%1,%7)\n\t"		  
                  "movntpd %%xmm1, 16(%1,%7)\n\t"		  
                  "movntpd %%xmm2, 32(%1,%7)\n\t"		  
                  "movntpd %%xmm3, 48(%1,%7)\n\t"		  

                  "shl $3, %9\n\t"		  
                  "sub %9, %10\n\t"		  
                  "shr $5, %9\n\t"		  
                  //Write Buffer END - Part 1


                  "add $16, %0\n\t"
                  "add $16, %9\n\t"
                  "mov %9, %7\n\t"
                  "and $63, %7\n\t"
                "jne SSE_FourierTrafoThirdFFTstep_InnerLoop\n\t"
                "sub $64, %0\n\t"


                "add $128, %0\n\t"
                "add $32, %5\n\t"
                "cmp $128,%5\n\t"
              "jl SSE_FourierTrafoThirdFFTstep_schleife4\n\t"
              "mov $0, %5\n\t"


              //New Write-Buffer-Index
              "mov 408(%6), %9\n\t"
              "sub %10, %9\n\t"
              "mov %9, %10\n\t"

              //New Output-Base
              "mov 72(%6), %9\n\t"
              "mov %9, 80(%6)\n\t"
              "mov 424(%6), %9\n\t"
              "mov %9, 432(%6)\n\t"
              "mov %9, %1\n\t"
              "add 24(%6), %1\n\t"
	      
              //New Input-Base
              "mov 416(%6), %9\n\t"
              "mov %9, 424(%6)\n\t"
              "mov 64(%6), %9\n\t"
              "mov %9, 72(%6)\n\t"
              "mov %9, %0\n\t"
              "add 16(%6), %0\n\t"
	      
              //New Prefetch-Base
              "mov 416(%6), %9\n\t"
              "add 152(%6), %9\n\t"
              "mov %9, 416(%6)\n\t"
              "mov 64(%6), %9\n\t"
              "add 128(%6), %9\n\t"
              "mov %9, 64(%6)\n\t"
              "mov %9, %8\n\t"
              "add 16(%6), %8\n\t"

              //New Counter - Values
              "mov 56(%6), %4\n\t"
              "mov 48(%6), %3\n\t"
              "mov 40(%6), %2\n\t"
	      
              "add $32, 56(%6)\n\t"
              "mov 0(%6), %9\n\t"	      
              "cmp %9,56(%6)\n\t"
            "jl SSE_FourierTrafoThirdFFTstep_MainLoop3\n\t"
            "mov $0, %9\n\t"	      
            "mov %9, 56(%6)\n\t"
	    
            //New Prefetch-Base
            "mov 416(%6), %9\n\t"
            "add 144(%6), %9\n\t"
            "mov %9, 416(%6)\n\t"
            "mov 64(%6), %9\n\t"
            "add 120(%6), %9\n\t"
            "mov %9, 64(%6)\n\t"
            "mov %9, %8\n\t"
            "add 16(%6), %8\n\t"

            "add $32, 48(%6)\n\t"
            "mov 0(%6), %9\n\t"	      
            "cmp %9,48(%6)\n\t"
          "jl SSE_FourierTrafoThirdFFTstep_MainLoop2\n\t"
          "mov $0, %9\n\t"	      
          "mov %9, 48(%6)\n\t"

          //New Prefetch-Base
          "mov 416(%6), %9\n\t"
          "add 136(%6), %9\n\t"
          "mov %9, 416(%6)\n\t"
          "mov 64(%6), %9\n\t"
          "add 112(%6), %9\n\t"
          "mov %9, 64(%6)\n\t"
          "mov %9, %8\n\t"
          "add 16(%6), %8\n\t"

          "add $32, 40(%6)\n\t"
          "mov 0(%6), %9\n\t"	      
          "cmp %9,40(%6)\n\t"
        "jl SSE_FourierTrafoThirdFFTstep_MainLoop1\n\t"  
	


        //Lead Out 1
        "mov $0, %9\n\t"
        "SSE_FourierTrafoThirdFFTstep_schleife4LEADOUT1: \n\t"
          "SSE_FourierTrafoThirdFFTstep_InnerLoopLEADOUT1: \n\t"
	
            // LOAD + Multiply with Weights
            "add %0, %14\n\t"
            #include "ExtremeFFTasmBlock5.h"
            "sub %0, %14\n\t"
	    

            // Second + third Step are external
            #include "ExtremeFFTasmBlock.h"


            // Forth STEP + SAVE
            "mov %9, %7\n\t"
            "and $1984, %7\n\t"
            "add %7, %9\n\t"
            "add 408(%6), %9\n\t"
            "sub %10, %9\n\t"		  
            #include "ExtremeFFTasmBlock6.h"
            "sub 408(%6), %9\n\t"
            "add %10, %9\n\t"		  
            "sub %7, %9\n\t"


            //Write Buffer - Part 1
            "shl $5, %9\n\t"		  
            "add %9, %10\n\t"		  
            "shr $3, %9\n\t"	
		  
            "mov 3584(%6,%9), %7\n\t"		  
            "movapd 18944(%6,%10), %%xmm0\n\t"		  
            "movapd 18960(%6,%10), %%xmm1\n\t"		  
            "movapd 18976(%6,%10), %%xmm2\n\t"		  
            "movapd 18992(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  
		  
            "mov 3600(%6,%9), %7\n\t"		  
            "movapd 19072(%6,%10), %%xmm0\n\t"		  
            "movapd 19088(%6,%10), %%xmm1\n\t"		  
            "movapd 19104(%6,%10), %%xmm2\n\t"		  
            "movapd 19120(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  
		  
            "mov 3616(%6,%9), %7\n\t"		  
            "movapd 19200(%6,%10), %%xmm0\n\t"		  
            "movapd 19216(%6,%10), %%xmm1\n\t"		  
            "movapd 19232(%6,%10), %%xmm2\n\t"		  
            "movapd 19248(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"	
		  
            "mov 3632(%6,%9), %7\n\t"		  
            "movapd 19328(%6,%10), %%xmm0\n\t"		  
            "movapd 19344(%6,%10), %%xmm1\n\t"		  
            "movapd 19360(%6,%10), %%xmm2\n\t"		  
            "movapd 19376(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  

            "shl $3, %9\n\t"		  
            "sub %9, %10\n\t"		  
            "shr $5, %9\n\t"		  
            //Write Buffer END - Part 1


            "add $16, %0\n\t"
            "add $16, %9\n\t"
            "mov %9, %7\n\t"
            "and $63, %7\n\t"
          "jne SSE_FourierTrafoThirdFFTstep_InnerLoopLEADOUT1\n\t"
          "sub $64, %0\n\t"


          "add $128, %0\n\t"
          "add $32, %5\n\t"
          "cmp 0(%6),%5\n\t"
        "jl SSE_FourierTrafoThirdFFTstep_schleife4LEADOUT1\n\t"	     
        "mov $0, %5\n\t"	
	

        //New Write-Buffer-Index
        "mov 408(%6), %9\n\t"
        "sub %10, %9\n\t"
        "mov %9, %10\n\t"

        //New Output-Base
        "mov 72(%6), %9\n\t"
        "mov %9, 80(%6)\n\t"
        "mov 424(%6), %9\n\t"
        "mov %9, 432(%6)\n\t"
        "mov %9, %1\n\t"
        "add 24(%6), %1\n\t"
	      
        //Lead Out 2
        "mov $0, %9\n\t"
        "SSE_FourierTrafoThirdFFTstep_schleife4LEADOUT2: \n\t"
          "SSE_FourierTrafoThirdFFTstep_InnerLoopLEADOUT2: \n\t"
 
            //Write Buffer - Part 1
            "shl $5, %9\n\t"		  
            "add %9, %10\n\t"		  
            "shr $3, %9\n\t"	
		  
            "mov 3584(%6,%9), %7\n\t"		  
            "movapd 18944(%6,%10), %%xmm0\n\t"		  
            "movapd 18960(%6,%10), %%xmm1\n\t"		  
            "movapd 18976(%6,%10), %%xmm2\n\t"		  
            "movapd 18992(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  
		  
            "mov 3600(%6,%9), %7\n\t"		  
            "movapd 19072(%6,%10), %%xmm0\n\t"		  
            "movapd 19088(%6,%10), %%xmm1\n\t"		  
            "movapd 19104(%6,%10), %%xmm2\n\t"		  
            "movapd 19120(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  
		  
            "mov 3616(%6,%9), %7\n\t"		  
            "movapd 19200(%6,%10), %%xmm0\n\t"		  
            "movapd 19216(%6,%10), %%xmm1\n\t"		  
            "movapd 19232(%6,%10), %%xmm2\n\t"		  
            "movapd 19248(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"	
		  
            "mov 3632(%6,%9), %7\n\t"		  
            "movapd 19328(%6,%10), %%xmm0\n\t"		  
            "movapd 19344(%6,%10), %%xmm1\n\t"		  
            "movapd 19360(%6,%10), %%xmm2\n\t"		  
            "movapd 19376(%6,%10), %%xmm3\n\t"		  
            "movntpd %%xmm0, (%1,%7)\n\t"		  
            "movntpd %%xmm1, 16(%1,%7)\n\t"		  
            "movntpd %%xmm2, 32(%1,%7)\n\t"		  
            "movntpd %%xmm3, 48(%1,%7)\n\t"		  

            "shl $3, %9\n\t"		  
            "sub %9, %10\n\t"		  
            "shr $5, %9\n\t"		  
            //Write Buffer END - Part 1 
 

            "add $16, %0\n\t"
            "add $16, %9\n\t"
            "mov %9, %7\n\t"
            "and $63, %7\n\t"
          "jne SSE_FourierTrafoThirdFFTstep_InnerLoopLEADOUT2\n\t"
          "sub $64, %0\n\t"


          "add $128, %0\n\t"
          "add $32, %5\n\t"
          "cmp 0(%6),%5\n\t"
        "jl SSE_FourierTrafoThirdFFTstep_schleife4LEADOUT2\n\t"	     
        "mov $0, %5\n\t"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10), "=r" (Dummy11),  "=r" (Dummy12), "=r" (Dummy13), "=r" (Dummy14), "=r" (Dummy15)
          : "r" (para_Comp));  
  
  

/*printf("reg0 %ld\n",(long int) Dummy1);
printf("reg1 %ld\n",(long int) Dummy2);
printf("reg2 %ld\n",(long int) Dummy3);
printf("reg3 %ld\n",(long int) Dummy4);
printf("reg4 %ld\n",(long int) Dummy5);
printf("reg5 %ld\n",(long int) Dummy6);
printf("reg6 %ld\n",(long int) Dummy7);
printf("reg7 %ld\n",(long int) Dummy8);
printf("reg8 %ld\n",(long int) Dummy9);
printf("reg9 %ld\n",(long int) Dummy10);
printf("reg10 %ld\n",(long int) Dummy11);
printf("reg11 %ld\n",(long int) Dummy12);
printf("reg12 %ld\n",(long int) Dummy13);
printf("reg13 %ld\n",(long int) Dummy14);
printf("reg14 %ld\n",(long int) Dummy15);
printf("reg15 %ld\n",(long int) para_Comp);

printBits((long int) Dummy12);
printBits((long int) Dummy13);
printBits((long int) Dummy14);
printBits((long int) Dummy15);
exit(0);*/

  //Always false
  if (para_int[2]!=para_int[2]) { 
    printf("%d\n",*Dummy1);
  }
}


/**** SSE_FourierTrafoSecondFFTstep ***
* Performs the second FFT step. Input data are expected to be in fixed
* size array.
* ATTENTION: Restriction to oneDim = 128
*
* Parameters: The following data are expected in ASMParameterData
* OneDim                : Lattice size in one dimension
* data                  : Input-/Output-Vector
* Incerements           : Increments for output base and output pi-shifted
*
**/
void SSE_FourierTrafoSecondFFTstep() {
//inline void SSE_FourierTrafoSecondFFTstep() {
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;
  int* Dummy11;
  int* Dummy12;
  int* Dummy13;
  int* Dummy14;
  int* Dummy15;  
  long int* para_int = (long int*) ASMParameterData;
  Complex* para_Comp = ASMParameterData;
  
  
/*printf("input %d\n",(long int) input);
printf("output %d\n",(long int) output);
printf("bit %d\n",(long int) bitInverse);
printf("para %d\n",(long int) para_Comp);*/
 
  
__asm__ volatile(
        "mov %15, %6\n\t"
        "mov 16(%6), %0\n\t"
        "mov 24(%6), %1\n\t"
        "mov 40(%6), %2\n\t"
        "mov 48(%6), %3\n\t"
        "mov 56(%6), %4\n\t"
        "mov 64(%6), %5\n\t"

        "mov 72(%6), %8\n\t"
        "mov 80(%6), %9\n\t"
        "mov 88(%6), %10\n\t"
        "mov 96(%6), %11\n\t"
	
        "mov 8(%6), %13\n\t"    //    <=== Hiermit koennen auch weniger Fourier-Trafos gemacht werden!!!
	
	
	//Addresses for prefetching
        "add $128,%1\n\t"
        "mov %1,%14\n\t"
        "mov %14, 2832(%6)\n\t"
	
        "add $256,%14\n\t"
        "mov %14, 2840(%6)\n\t"
	
        "sub $256,%14\n\t"
        "add %8,%14\n\t"
        "mov %14, 2848(%6)\n\t"
	
        "add $256,%14\n\t"
        "mov %14, 2856(%6)\n\t"
	
        "sub $256,%14\n\t"
        "sub %8,%14\n\t"
        "add %9,%14\n\t"
        "mov %14, 2864(%6)\n\t"
	
        "add $256,%14\n\t"
        "mov %14, 2872(%6)\n\t"
	
        "sub $256,%14\n\t"
        "add %8,%14\n\t"
        "mov %14, 2880(%6)\n\t"

        "add $256,%14\n\t"
        "mov %14, 2888(%6)\n\t"

	
        "sub $256,%14\n\t"
        "sub %8,%14\n\t"
        "sub %9,%14\n\t"
        "add %11,%14\n\t"
        "mov %14, 2896(%6)\n\t"
	
        "add $256,%14\n\t"
        "mov %14, 2904(%6)\n\t"
	
        "sub $256,%14\n\t"
        "add %8,%14\n\t"
        "mov %14, 2912(%6)\n\t"
	
        "add $256,%14\n\t"
        "mov %14, 2920(%6)\n\t"
	
        "sub $256,%14\n\t"
        "sub %8,%14\n\t"
        "add %9,%14\n\t"
        "mov %14, 2928(%6)\n\t"
	
        "add $256,%14\n\t"
        "mov %14, 2936(%6)\n\t"
	
        "sub $256,%14\n\t"
        "add %8,%14\n\t"
        "mov %14, 2944(%6)\n\t"

        "add $256,%14\n\t"
        "mov %14, 2952(%6)\n\t"
	
        "sub $128,%1\n\t"
	
	
	
    // *** MODE 0,0,0,0 ***	
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop0000: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
		
          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movapd 256(%0,%7), %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movapd 1024(%0,%7), %%xmm2\n\t"
          "movapd 1280(%0,%7), %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movapd 4096(%0,%7), %%xmm4\n\t"
          "movapd 4352(%0,%7), %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movapd 5120(%0,%7), %%xmm6\n\t"
          "movapd 5376(%0,%7), %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movapd 16384(%0,%7), %%xmm8\n\t"
          "movapd 16640(%0,%7), %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movapd 17408(%0,%7), %%xmm10\n\t"
          "movapd 17664(%0,%7), %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movapd 20480(%0,%7), %%xmm12\n\t"
          "movapd 20736(%0,%7), %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movapd 21504(%0,%7), %%xmm14\n\t"
          "movapd 21760(%0,%7), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %4,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %4,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop0000\n\t"
        "sub %13, %1\n\t"
	
	
    // *** MODE 0,0,0,1 ***	
        "add $128, %0\n\t"
        "add %5, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop0001: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm14\n\t"
	  
		
          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movhpd 256(%0,%7), %%xmm1\n\t"
          "movlpd 264(%0,%7), %%xmm1\n\t"
          "xorpd %%xmm14, %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movapd 1024(%0,%7), %%xmm2\n\t"
          "movhpd 1280(%0,%7), %%xmm3\n\t"
          "movlpd 1288(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm14, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movapd 4096(%0,%7), %%xmm4\n\t"
          "movhpd 4352(%0,%7), %%xmm5\n\t"
          "movlpd 4360(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm14, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movapd 5120(%0,%7), %%xmm6\n\t"
          "movhpd 5376(%0,%7), %%xmm7\n\t"
          "movlpd 5384(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm14, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movapd 16384(%0,%7), %%xmm8\n\t"
          "movhpd 16640(%0,%7), %%xmm9\n\t"
          "movlpd 16648(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm14, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movapd 17408(%0,%7), %%xmm10\n\t"
          "movhpd 17664(%0,%7), %%xmm11\n\t"
          "movlpd 17672(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm14, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movapd 20480(%0,%7), %%xmm12\n\t"
          "movhpd 20736(%0,%7), %%xmm13\n\t"
          "movlpd 20744(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movapd 21504(%0,%7), %%xmm14\n\t"
          "movhpd 21760(%0,%7), %%xmm15\n\t"
          "movlpd 21768(%0,%7), %%xmm15\n\t"
          "xorpd 192(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop0001\n\t"
        "sub %13, %1\n\t"
	
	
    // *** MODE 0,0,1,0 ***	
        "add $384, %0\n\t"
        "add %4, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop0010: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm14\n\t"
	  
          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movapd 256(%0,%7), %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movhpd 1024(%0,%7), %%xmm2\n\t"
          "movlpd 1032(%0,%7), %%xmm2\n\t"
          "xorpd %%xmm14, %%xmm2\n\t"
          "movhpd 1280(%0,%7), %%xmm3\n\t"
          "movlpd 1288(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm14, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movapd 4096(%0,%7), %%xmm4\n\t"
          "movapd 4352(%0,%7), %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movhpd 5120(%0,%7), %%xmm6\n\t"
          "movlpd 5128(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm14, %%xmm6\n\t"
          "movhpd 5376(%0,%7), %%xmm7\n\t"
          "movlpd 5384(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm14, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movapd 16384(%0,%7), %%xmm8\n\t"
          "movapd 16640(%0,%7), %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movhpd 17408(%0,%7), %%xmm10\n\t"
          "movlpd 17416(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm14, %%xmm10\n\t"
          "movhpd 17664(%0,%7), %%xmm11\n\t"
          "movlpd 17672(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm14, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movapd 20480(%0,%7), %%xmm12\n\t"
          "movapd 20736(%0,%7), %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movhpd 21504(%0,%7), %%xmm14\n\t"
          "movlpd 21512(%0,%7), %%xmm14\n\t"
          "xorpd 192(%6), %%xmm14\n\t"
          "movhpd 21760(%0,%7), %%xmm15\n\t"
          "movlpd 21768(%0,%7), %%xmm15\n\t"
          "xorpd 192(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %3,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %3,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop0010\n\t"
        "sub %13, %1\n\t"	
	
	
    // *** MODE 0,0,1,1 ***	
        "add $128, %0\n\t"
        "add %5, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop0011: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm14\n\t"
          "movapd 160(%6), %%xmm13\n\t"

	  
          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movhpd 256(%0,%7), %%xmm1\n\t"
          "movlpd 264(%0,%7), %%xmm1\n\t"
          "xorpd %%xmm14, %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movhpd 1024(%0,%7), %%xmm2\n\t"
          "movlpd 1032(%0,%7), %%xmm2\n\t"
          "xorpd %%xmm14, %%xmm2\n\t"
          "movapd 1280(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm13, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movapd 4096(%0,%7), %%xmm4\n\t"
          "movhpd 4352(%0,%7), %%xmm5\n\t"
          "movlpd 4360(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm14, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movhpd 5120(%0,%7), %%xmm6\n\t"
          "movlpd 5128(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm14, %%xmm6\n\t"
          "movapd 5376(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm13, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movapd 16384(%0,%7), %%xmm8\n\t"
          "movhpd 16640(%0,%7), %%xmm9\n\t"
          "movlpd 16648(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm14, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movhpd 17408(%0,%7), %%xmm10\n\t"
          "movlpd 17416(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm14, %%xmm10\n\t"
          "movapd 17664(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm13, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movapd 20480(%0,%7), %%xmm12\n\t"
          "movhpd 20736(%0,%7), %%xmm13\n\t"
          "movlpd 20744(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movhpd 21504(%0,%7), %%xmm14\n\t"
          "movlpd 21512(%0,%7), %%xmm14\n\t"
          "xorpd 192(%6), %%xmm14\n\t"
          "movapd 21760(%0,%7), %%xmm15\n\t"
          "xorpd 160(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop0011\n\t"
        "sub %13, %1\n\t"	


    // *** MODE 0,1,0,0 ***	
        "add $1408, %0\n\t"
        "add %3, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop0100: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm14\n\t"

	  
          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movapd 256(%0,%7), %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movapd 1024(%0,%7), %%xmm2\n\t"
          "movapd 1280(%0,%7), %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movhpd 4096(%0,%7), %%xmm4\n\t"
          "movlpd 4104(%0,%7), %%xmm4\n\t"
          "xorpd %%xmm14, %%xmm4\n\t"
          "movhpd 4352(%0,%7), %%xmm5\n\t"
          "movlpd 4360(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm14, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movhpd 5120(%0,%7), %%xmm6\n\t"
          "movlpd 5128(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm14, %%xmm6\n\t"
          "movhpd 5376(%0,%7), %%xmm7\n\t"
          "movlpd 5384(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm14, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movapd 16384(%0,%7), %%xmm8\n\t"
          "movapd 16640(%0,%7), %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movapd 17408(%0,%7), %%xmm10\n\t"
          "movapd 17664(%0,%7), %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movhpd 20480(%0,%7), %%xmm12\n\t"
          "movlpd 20488(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm14, %%xmm12\n\t"
          "movhpd 20736(%0,%7), %%xmm13\n\t"
          "movlpd 20744(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movhpd 21504(%0,%7), %%xmm14\n\t"
          "movlpd 21512(%0,%7), %%xmm14\n\t"
          "xorpd 192(%6), %%xmm14\n\t"
          "movhpd 21760(%0,%7), %%xmm15\n\t"
          "movlpd 21768(%0,%7), %%xmm15\n\t"
          "xorpd 192(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %4,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %4,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop0100\n\t"
        "sub %13, %1\n\t"	


    // *** MODE 0,1,0,1 ***	
        "add $128, %0\n\t"
        "add %5, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop0101: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 160(%6), %%xmm14\n\t"
          "movapd 192(%6), %%xmm13\n\t"

          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movhpd 256(%0,%7), %%xmm1\n\t"
          "movlpd 264(%0,%7), %%xmm1\n\t"
          "xorpd %%xmm13, %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movapd 1024(%0,%7), %%xmm2\n\t"
          "movhpd 1280(%0,%7), %%xmm3\n\t"
          "movlpd 1288(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm13, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movhpd 4096(%0,%7), %%xmm4\n\t"
          "movlpd 4104(%0,%7), %%xmm4\n\t"
          "xorpd %%xmm13, %%xmm4\n\t"
          "movapd 4352(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm14, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movhpd 5120(%0,%7), %%xmm6\n\t"
          "movlpd 5128(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm13, %%xmm6\n\t"
          "movapd 5376(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm14, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movapd 16384(%0,%7), %%xmm8\n\t"
          "movhpd 16640(%0,%7), %%xmm9\n\t"
          "movlpd 16648(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm13, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movapd 17408(%0,%7), %%xmm10\n\t"
          "movhpd 17664(%0,%7), %%xmm11\n\t"
          "movlpd 17672(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm13, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movhpd 20480(%0,%7), %%xmm12\n\t"
          "movlpd 20488(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm13, %%xmm12\n\t"
          "movapd 20736(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movhpd 21504(%0,%7), %%xmm14\n\t"
          "movlpd 21512(%0,%7), %%xmm14\n\t"
          "xorpd 192(%6), %%xmm14\n\t"
          "movapd 21760(%0,%7), %%xmm15\n\t"
          "xorpd 160(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"	  

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop0101\n\t"
        "sub %13, %1\n\t"


    // *** MODE 0,1,1,0 ***	
        "add $384, %0\n\t"
        "add %4, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop0110: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm14\n\t"

          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movapd 256(%0,%7), %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movhpd 1024(%0,%7), %%xmm2\n\t"
          "movlpd 1032(%0,%7), %%xmm2\n\t"
          "xorpd %%xmm14, %%xmm2\n\t"
          "movhpd 1280(%0,%7), %%xmm3\n\t"
          "movlpd 1288(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm14, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movhpd 4096(%0,%7), %%xmm4\n\t"
          "movlpd 4104(%0,%7), %%xmm4\n\t"
          "xorpd %%xmm14, %%xmm4\n\t"
          "movhpd 4352(%0,%7), %%xmm5\n\t"
          "movlpd 4360(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm14, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movapd 5120(%0,%7), %%xmm6\n\t"
          "xorpd 160(%6), %%xmm6\n\t"
          "movapd 5376(%0,%7), %%xmm7\n\t"
          "xorpd 160(%6), %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movapd 16384(%0,%7), %%xmm8\n\t"
          "movapd 16640(%0,%7), %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movhpd 17408(%0,%7), %%xmm10\n\t"
          "movlpd 17416(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm14, %%xmm10\n\t"
          "movhpd 17664(%0,%7), %%xmm11\n\t"
          "movlpd 17672(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm14, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movhpd 20480(%0,%7), %%xmm12\n\t"
          "movlpd 20488(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm14, %%xmm12\n\t"
          "movhpd 20736(%0,%7), %%xmm13\n\t"
          "movlpd 20744(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movapd 21504(%0,%7), %%xmm14\n\t"
          "xorpd 160(%6), %%xmm14\n\t"
          "movapd 21760(%0,%7), %%xmm15\n\t"
          "xorpd 160(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %2,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %2,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop0110\n\t"
        "sub %13, %1\n\t"
	

    // *** MODE 0,1,1,1 ***	
        "add $128, %0\n\t"
        "add %5, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop0111: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 160(%6), %%xmm14\n\t"	  
          "movapd 192(%6), %%xmm13\n\t"

          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movhpd 256(%0,%7), %%xmm1\n\t"
          "movlpd 264(%0,%7), %%xmm1\n\t"
          "xorpd %%xmm13, %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movhpd 1024(%0,%7), %%xmm2\n\t"
          "movlpd 1032(%0,%7), %%xmm2\n\t"
          "xorpd %%xmm13, %%xmm2\n\t"
          "movapd 1280(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm14, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movhpd 4096(%0,%7), %%xmm4\n\t"
          "movlpd 4104(%0,%7), %%xmm4\n\t"
          "xorpd %%xmm13, %%xmm4\n\t"
          "movapd 4352(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm14, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movapd 5120(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm14, %%xmm6\n\t"
          "movhpd 5376(%0,%7), %%xmm7\n\t"
          "movlpd 5384(%0,%7), %%xmm7\n\t"
          "xorpd 176(%6), %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movapd 16384(%0,%7), %%xmm8\n\t"
          "movhpd 16640(%0,%7), %%xmm9\n\t"
          "movlpd 16648(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm13, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movhpd 17408(%0,%7), %%xmm10\n\t"
          "movlpd 17416(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm13, %%xmm10\n\t"
          "movapd 17664(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm14, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movhpd 20480(%0,%7), %%xmm12\n\t"
          "movlpd 20488(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm13, %%xmm12\n\t"
          "movapd 20736(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movapd 21504(%0,%7), %%xmm14\n\t"
          "xorpd 160(%6), %%xmm14\n\t"
          "movhpd 21760(%0,%7), %%xmm15\n\t"
          "movlpd 21768(%0,%7), %%xmm15\n\t"
          "xorpd 176(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop0111\n\t"
        "sub %13, %1\n\t"


    // *** MODE 1,0,0,0 ***	
        "add $5504, %0\n\t"
        "add %2, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop1000: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm14\n\t"

          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movapd 256(%0,%7), %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movapd 1024(%0,%7), %%xmm2\n\t"
          "movapd 1280(%0,%7), %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movapd 4096(%0,%7), %%xmm4\n\t"
          "movapd 4352(%0,%7), %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movapd 5120(%0,%7), %%xmm6\n\t"
          "movapd 5376(%0,%7), %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movhpd 16384(%0,%7), %%xmm8\n\t"
          "movlpd 16392(%0,%7), %%xmm8\n\t"
          "xorpd %%xmm14, %%xmm8\n\t"
          "movhpd 16640(%0,%7), %%xmm9\n\t"
          "movlpd 16648(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm14, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movhpd 17408(%0,%7), %%xmm10\n\t"
          "movlpd 17416(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm14, %%xmm10\n\t"
          "movhpd 17664(%0,%7), %%xmm11\n\t"
          "movlpd 17672(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm14, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movhpd 20480(%0,%7), %%xmm12\n\t"
          "movlpd 20488(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm14, %%xmm12\n\t"
          "movhpd 20736(%0,%7), %%xmm13\n\t"
          "movlpd 20744(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movhpd 21504(%0,%7), %%xmm14\n\t"
          "movlpd 21512(%0,%7), %%xmm14\n\t"
          "xorpd 192(%6), %%xmm14\n\t"
          "movhpd 21760(%0,%7), %%xmm15\n\t"
          "movlpd 21768(%0,%7), %%xmm15\n\t"
          "xorpd 192(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %4,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %4,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop1000\n\t"
        "sub %13, %1\n\t"

	
    // *** MODE 1,0,0,1 ***	
        "add $128, %0\n\t"
        "add %5, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop1001: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm13\n\t"
          "movapd 160(%6), %%xmm14\n\t"

          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movhpd 256(%0,%7), %%xmm1\n\t"
          "movlpd 264(%0,%7), %%xmm1\n\t"
          "xorpd %%xmm13, %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movapd 1024(%0,%7), %%xmm2\n\t"
          "movhpd 1280(%0,%7), %%xmm3\n\t"
          "movlpd 1288(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm13, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movapd 4096(%0,%7), %%xmm4\n\t"
          "movhpd 4352(%0,%7), %%xmm5\n\t"
          "movlpd 4360(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm13, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movapd 5120(%0,%7), %%xmm6\n\t"
          "movhpd 5376(%0,%7), %%xmm7\n\t"
          "movlpd 5384(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm13, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movhpd 16384(%0,%7), %%xmm8\n\t"
          "movlpd 16392(%0,%7), %%xmm8\n\t"
          "xorpd %%xmm13, %%xmm8\n\t"
          "movapd 16640(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm14, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movhpd 17408(%0,%7), %%xmm10\n\t"
          "movlpd 17416(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm13, %%xmm10\n\t"
          "movapd 17664(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm14, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movhpd 20480(%0,%7), %%xmm12\n\t"
          "movlpd 20488(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm13, %%xmm12\n\t"
          "movapd 20736(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movhpd 21504(%0,%7), %%xmm14\n\t"
          "movlpd 21512(%0,%7), %%xmm14\n\t"
          "xorpd 192(%6), %%xmm14\n\t"
          "movapd 21760(%0,%7), %%xmm15\n\t"
          "xorpd 160(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop1001\n\t"
        "sub %13, %1\n\t"

    // *** MODE 1,0,1,0 ***	
        "add $384, %0\n\t"
        "add %4, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop1010: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm14\n\t"

          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movapd 256(%0,%7), %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movhpd 1024(%0,%7), %%xmm2\n\t"
          "movlpd 1032(%0,%7), %%xmm2\n\t"
          "xorpd %%xmm14, %%xmm2\n\t"
          "movhpd 1280(%0,%7), %%xmm3\n\t"
          "movlpd 1288(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm14, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movapd 4096(%0,%7), %%xmm4\n\t"
          "movapd 4352(%0,%7), %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movhpd 5120(%0,%7), %%xmm6\n\t"
          "movlpd 5128(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm14, %%xmm6\n\t"
          "movhpd 5376(%0,%7), %%xmm7\n\t"
          "movlpd 5384(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm14, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movhpd 16384(%0,%7), %%xmm8\n\t"
          "movlpd 16392(%0,%7), %%xmm8\n\t"
          "xorpd %%xmm14, %%xmm8\n\t"
          "movhpd 16640(%0,%7), %%xmm9\n\t"
          "movlpd 16648(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm14, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movapd 17408(%0,%7), %%xmm10\n\t"
          "xorpd 160(%6), %%xmm10\n\t"
          "movapd 17664(%0,%7), %%xmm11\n\t"
          "xorpd 160(%6), %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movhpd 20480(%0,%7), %%xmm12\n\t"
          "movlpd 20488(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm14, %%xmm12\n\t"
          "movhpd 20736(%0,%7), %%xmm13\n\t"
          "movlpd 20744(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movapd 21504(%0,%7), %%xmm14\n\t"
          "xorpd 160(%6), %%xmm14\n\t"
          "movapd 21760(%0,%7), %%xmm15\n\t"
          "xorpd 160(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %3,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %3,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop1010\n\t"
        "sub %13, %1\n\t"

	
    // *** MODE 1,0,1,1 ***	
        "add $128, %0\n\t"
        "add %5, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop1011: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm13\n\t"
          "movapd 160(%6), %%xmm14\n\t"

          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movhpd 256(%0,%7), %%xmm1\n\t"
          "movlpd 264(%0,%7), %%xmm1\n\t"
          "xorpd %%xmm13, %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movhpd 1024(%0,%7), %%xmm2\n\t"
          "movlpd 1032(%0,%7), %%xmm2\n\t"
          "xorpd %%xmm13, %%xmm2\n\t"
          "movapd 1280(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm14, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movapd 4096(%0,%7), %%xmm4\n\t"
          "movhpd 4352(%0,%7), %%xmm5\n\t"
          "movlpd 4360(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm13, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movhpd 5120(%0,%7), %%xmm6\n\t"
          "movlpd 5128(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm13, %%xmm6\n\t"
          "movapd 5376(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm14, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movhpd 16384(%0,%7), %%xmm8\n\t"
          "movlpd 16392(%0,%7), %%xmm8\n\t"
          "xorpd %%xmm13, %%xmm8\n\t"
          "movapd 16640(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm14, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movapd 17408(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm14, %%xmm10\n\t"
          "movhpd 17664(%0,%7), %%xmm11\n\t"
          "movlpd 17672(%0,%7), %%xmm11\n\t"
          "xorpd 176(%6), %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movhpd 20480(%0,%7), %%xmm12\n\t"
          "movlpd 20488(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm13, %%xmm12\n\t"
          "movapd 20736(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movapd 21504(%0,%7), %%xmm14\n\t"
          "xorpd 160(%6), %%xmm14\n\t"
          "movhpd 21760(%0,%7), %%xmm15\n\t"
          "movlpd 21768(%0,%7), %%xmm15\n\t"
          "xorpd 176(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop1011\n\t"
        "sub %13, %1\n\t"	
	
	
    // *** MODE 1,1,0,0 ***	
        "add $1408, %0\n\t"
        "add %3, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop1100: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm13\n\t"

          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movapd 256(%0,%7), %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movapd 1024(%0,%7), %%xmm2\n\t"
          "movapd 1280(%0,%7), %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movhpd 4096(%0,%7), %%xmm4\n\t"
          "movlpd 4104(%0,%7), %%xmm4\n\t"
          "xorpd %%xmm13, %%xmm4\n\t"
          "movhpd 4352(%0,%7), %%xmm5\n\t"
          "movlpd 4360(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm13, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movhpd 5120(%0,%7), %%xmm6\n\t"
          "movlpd 5128(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm13, %%xmm6\n\t"
          "movhpd 5376(%0,%7), %%xmm7\n\t"
          "movlpd 5384(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm13, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movhpd 16384(%0,%7), %%xmm8\n\t"
          "movlpd 16392(%0,%7), %%xmm8\n\t"
          "xorpd %%xmm13, %%xmm8\n\t"
          "movhpd 16640(%0,%7), %%xmm9\n\t"
          "movlpd 16648(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm13, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movhpd 17408(%0,%7), %%xmm10\n\t"
          "movlpd 17416(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm13, %%xmm10\n\t"
          "movhpd 17664(%0,%7), %%xmm11\n\t"
          "movlpd 17672(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm13, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movapd 20480(%0,%7), %%xmm12\n\t"
          "xorpd 160(%6), %%xmm12\n\t"
          "movapd 20736(%0,%7), %%xmm13\n\t"
          "xorpd 160(%6), %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movapd 21504(%0,%7), %%xmm14\n\t"
          "xorpd 160(%6), %%xmm14\n\t"
          "movapd 21760(%0,%7), %%xmm15\n\t"
          "xorpd 160(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %4,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %4,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop1100\n\t"
        "sub %13, %1\n\t"	
	
	
    // *** MODE 1,1,0,1 ***	
        "add $128, %0\n\t"
        "add %5, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop1101: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm13\n\t"
          "movapd 160(%6), %%xmm14\n\t"

          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movhpd 256(%0,%7), %%xmm1\n\t"
          "movlpd 264(%0,%7), %%xmm1\n\t"
          "xorpd %%xmm13, %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movapd 1024(%0,%7), %%xmm2\n\t"
          "movhpd 1280(%0,%7), %%xmm3\n\t"
          "movlpd 1288(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm13, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movhpd 4096(%0,%7), %%xmm4\n\t"
          "movlpd 4104(%0,%7), %%xmm4\n\t"
          "xorpd %%xmm13, %%xmm4\n\t"
          "movapd 4352(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm14, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movhpd 5120(%0,%7), %%xmm6\n\t"
          "movlpd 5128(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm13, %%xmm6\n\t"
          "movapd 5376(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm14, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movhpd 16384(%0,%7), %%xmm8\n\t"
          "movlpd 16392(%0,%7), %%xmm8\n\t"
          "xorpd %%xmm13, %%xmm8\n\t"
          "movapd 16640(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm14, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movhpd 17408(%0,%7), %%xmm10\n\t"
          "movlpd 17416(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm13, %%xmm10\n\t"
          "movapd 17664(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm14, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movapd 20480(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm14, %%xmm12\n\t"
          "movhpd 20736(%0,%7), %%xmm13\n\t"
          "movlpd 20744(%0,%7), %%xmm13\n\t"
          "xorpd 176(%6), %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movapd 21504(%0,%7), %%xmm14\n\t"
          "xorpd 160(%6), %%xmm14\n\t"
          "movhpd 21760(%0,%7), %%xmm15\n\t"
          "movlpd 21768(%0,%7), %%xmm15\n\t"
          "xorpd 176(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
          "add %5,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop1101\n\t"
        "sub %13, %1\n\t"		
	
	
    // *** MODE 1,1,1,0 ***	
        "add $384, %0\n\t"
        "add %4, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop1110: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm13\n\t"
          "movapd 160(%6), %%xmm14\n\t"

          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movapd 256(%0,%7), %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movhpd 1024(%0,%7), %%xmm2\n\t"
          "movlpd 1032(%0,%7), %%xmm2\n\t"
          "xorpd %%xmm13, %%xmm2\n\t"
          "movhpd 1280(%0,%7), %%xmm3\n\t"
          "movlpd 1288(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm13, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movhpd 4096(%0,%7), %%xmm4\n\t"
          "movlpd 4104(%0,%7), %%xmm4\n\t"
          "xorpd %%xmm13, %%xmm4\n\t"
          "movhpd 4352(%0,%7), %%xmm5\n\t"
          "movlpd 4360(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm13, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movapd 5120(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm14, %%xmm6\n\t"
          "movapd 5376(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm14, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movhpd 16384(%0,%7), %%xmm8\n\t"
          "movlpd 16392(%0,%7), %%xmm8\n\t"
          "xorpd %%xmm13, %%xmm8\n\t"
          "movhpd 16640(%0,%7), %%xmm9\n\t"
          "movlpd 16648(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm13, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movapd 17408(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm14, %%xmm10\n\t"
          "movapd 17664(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm14, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movapd 20480(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm14, %%xmm12\n\t"
          "movapd 20736(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movhpd 21504(%0,%7), %%xmm14\n\t"
          "movlpd 21512(%0,%7), %%xmm14\n\t"
          "xorpd 176(%6), %%xmm14\n\t"
          "movhpd 21760(%0,%7), %%xmm15\n\t"
          "movlpd 21768(%0,%7), %%xmm15\n\t"
          "xorpd 176(%6), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

          //Prefetching
          "mov 2832(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
//          "add %4,%12\n\t"
          "mov %12,2832(%6,%14)\n\t"
          "mov 2840(%6,%14),%12\n\t"
          "prefetchw (%12)\n\t"
          "prefetchw 64(%12)\n\t"
//          "add %4,%12\n\t"
          "mov %12,2840(%6,%14)\n\t"
          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop1110\n\t"
        "sub %13, %1\n\t"	
	
    // *** MODE 1,1,1,1 ***	
        "add $128, %0\n\t"
        "add %5, %1\n\t"
        "mov  $0, %7\n\t"
        "mov  $0, %14\n\t"
        "SSE_FourierTrafoSecondFFTstep_InnerLoop1111: \n\t"
          "movapd 208(%6), %%xmm15\n\t"
          "movapd 192(%6), %%xmm12\n\t"
          "movapd 176(%6), %%xmm14\n\t"

          "movapd 160(%6), %%xmm13\n\t"


          // LOAD + First STEP
          "movapd (%0,%7), %%xmm0\n\t"
          "movhpd 256(%0,%7), %%xmm1\n\t"
          "movlpd 264(%0,%7), %%xmm1\n\t"
          "xorpd %%xmm12, %%xmm1\n\t"
          "addpd %%xmm0, %%xmm1\n\t"
          "mulpd %%xmm15, %%xmm0\n\t"
          "subpd %%xmm1, %%xmm0\n\t"

          "movhpd 1024(%0,%7), %%xmm2\n\t"
          "movlpd 1032(%0,%7), %%xmm2\n\t"
          "xorpd %%xmm12, %%xmm2\n\t"
          "movapd 1280(%0,%7), %%xmm3\n\t"
          "xorpd %%xmm13, %%xmm3\n\t"
          "addpd %%xmm2, %%xmm3\n\t"
          "mulpd %%xmm15, %%xmm2\n\t"
          "subpd %%xmm3, %%xmm2\n\t"

          "movhpd 4096(%0,%7), %%xmm4\n\t"
          "movlpd 4104(%0,%7), %%xmm4\n\t"
          "xorpd %%xmm12, %%xmm4\n\t"
          "movapd 4352(%0,%7), %%xmm5\n\t"
          "xorpd %%xmm13, %%xmm5\n\t"
          "addpd %%xmm4, %%xmm5\n\t"
          "mulpd %%xmm15, %%xmm4\n\t"
          "subpd %%xmm5, %%xmm4\n\t"

          "movapd 5120(%0,%7), %%xmm6\n\t"
          "xorpd %%xmm13, %%xmm6\n\t"
          "movhpd 5376(%0,%7), %%xmm7\n\t"
          "movlpd 5384(%0,%7), %%xmm7\n\t"
          "xorpd %%xmm14, %%xmm7\n\t"
          "addpd %%xmm6, %%xmm7\n\t"
          "mulpd %%xmm15, %%xmm6\n\t"
          "subpd %%xmm7, %%xmm6\n\t"

          "movhpd 16384(%0,%7), %%xmm8\n\t"
          "movlpd 16392(%0,%7), %%xmm8\n\t"
          "xorpd %%xmm12, %%xmm8\n\t"
          "movapd 16640(%0,%7), %%xmm9\n\t"
          "xorpd %%xmm13, %%xmm9\n\t"
          "addpd %%xmm8, %%xmm9\n\t"
          "mulpd %%xmm15, %%xmm8\n\t"
          "subpd %%xmm9, %%xmm8\n\t"

          "movapd 17408(%0,%7), %%xmm10\n\t"
          "xorpd %%xmm13, %%xmm10\n\t"
          "movhpd 17664(%0,%7), %%xmm11\n\t"
          "movlpd 17672(%0,%7), %%xmm11\n\t"
          "xorpd %%xmm14, %%xmm11\n\t"
          "addpd %%xmm10, %%xmm11\n\t"
          "mulpd %%xmm15, %%xmm10\n\t"
          "subpd %%xmm11, %%xmm10\n\t"

          "movapd 20480(%0,%7), %%xmm12\n\t"
          "xorpd %%xmm13, %%xmm12\n\t"
          "movhpd 20736(%0,%7), %%xmm13\n\t"
          "movlpd 20744(%0,%7), %%xmm13\n\t"
          "xorpd %%xmm14, %%xmm13\n\t"
          "addpd %%xmm12, %%xmm13\n\t"
          "mulpd %%xmm15, %%xmm12\n\t"
          "subpd %%xmm13, %%xmm12\n\t"

          "movhpd 21504(%0,%7), %%xmm14\n\t"
          "movlpd 21512(%0,%7), %%xmm14\n\t"
          "xorpd 176(%6), %%xmm14\n\t"
          "movapd 21760(%0,%7), %%xmm15\n\t"
          "addpd %%xmm14, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm14\n\t"
          "subpd %%xmm15, %%xmm14\n\t"

//          //Prefetching
//          "mov 2832(%6,%14),%12\n\t"
//          "prefetchw (%12)\n\t"
//          "prefetchw 64(%12)\n\t"
//          "add %4,%12\n\t"
//          "mov %12,2832(%6,%14)\n\t"
//          "mov 2840(%6,%14),%12\n\t"
//          "prefetchw (%12)\n\t"
//          "prefetchw 64(%12)\n\t"
//          "add %4,%12\n\t"
//          "mov %12,2840(%6,%14)\n\t"
//          "add $16,%14\n\t"

          // Second + third Step are external
          #include "ExtremeFFTasmBlock.h"
          "mov %1, %12\n\t"
          "add %11, %12\n\t"

          // Forth STEP + SAVE is external
          #include "ExtremeFFTasmBlock2.h"
	
          "add $16, %1\n\t"
          "add $16, %7\n\t"
          "cmp %13, %7\n\t"
        "jl SSE_FourierTrafoSecondFFTstep_InnerLoop1111\n\t"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10), "=r" (Dummy11),  "=r" (Dummy12), "=r" (Dummy13), "=r" (Dummy14), "=r" (Dummy15)
          : "r" (para_Comp));

/*printf("reg0 %d\n",(long int) Dummy1);
printf("reg1 %d\n",(long int) Dummy2);
printf("reg2 %d\n",(long int) Dummy3);
printf("reg3 %d\n",(long int) Dummy4);
printf("reg4 %d\n",(long int) Dummy5);
printf("reg5 %d\n",(long int) Dummy6);
printf("reg6 %d\n",(long int) Dummy7);
printf("reg7 %d\n",(long int) Dummy8);
printf("reg8 %d\n",(long int) Dummy9);
printf("reg9 %d\n",(long int) Dummy10);
printf("reg10 %d\n",(long int) Dummy11);
printf("reg11 %d\n",(long int) Dummy12);
printf("reg12 %d\n",(long int) Dummy13);
printf("reg13 %d\n",(long int) Dummy14);
printf("reg14 %d\n",(long int) Dummy15);
printf("reg15 %d\n",(long int) para_Comp);

printBits((long int) Dummy12);
printBits((long int) Dummy13);
printBits((long int) Dummy14);
printBits((long int) Dummy15);
exit(0);*/

  //Always false
  if (para_int[2]!=para_int[2]) { 
    printf("%d\n",*Dummy1);
  }
}


/**** SSE_FourierTrafoRearrangerWithFirstFFTstep ***
* Rearranges the input-vector as prescribed by the FFT-algorithm.
* ATTENTION: Restriction to oneDim = 128
*
* Parameters: The following data are expected in ASMParameterData
* OneDim                : Lattice size in one dimension
* data                  : Input-/Output-Vector
* BitInverse            : Bit-Inversion data
* b0,b1,b2,b3           : Integer position of this block (beginning)
*
**/
void SSE_FourierTrafoRearrangerWithFirstFFTstep() {
//inline void SSE_FourierTrafoRearrangerWithFirstFFTstep() {
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;
  int* Dummy11;
  int* Dummy12;
  int* Dummy13;
  int* Dummy14;
  int* Dummy15;  
  long int* para_int = (long int*) ASMParameterData;
  Complex* para_Comp = ASMParameterData;
  
  
/*printf("input %d\n",(long int) input);
printf("output %d\n",(long int) output);
printf("bit %d\n",(long int) bitInverse);
printf("para %d\n",(long int) para_Comp);*/
 
  
__asm__ volatile(
        "mov %15, %6\n\t"
        "mov 32(%6), %7\n\t"
        "mov  280(%6), %2\n\t"
        "mov  272(%6), %3\n\t"
        "mov  264(%6), %4\n\t"
        "mov  256(%6), %5\n\t"
        "mov  $0, %8\n\t"
        "mov  $0, %9\n\t"
	
        "SSE_FourierTrafoRearrangerWithFirstFFTstep_schleife1: \n\t"
          "SSE_FourierTrafoRearrangerWithFirstFFTstep_schleife2: \n\t"
            "SSE_FourierTrafoRearrangerWithFirstFFTstep_schleife3: \n\t"
              "SSE_FourierTrafoRearrangerWithFirstFFTstep_schleife4: \n\t"
                //Berechne 16 Input - Base Adressen -- Erste 4 Stück
                "mov  3072(%7,%5), %10\n\t"
                "add   2048(%7,%4), %10\n\t"
                "mov  1024(%7,%3), %11\n\t"
                "mov  1032(%7,%3), %12\n\t"
                "mov      (%7,%2), %13\n\t"
                "mov     8(%7,%2), %14\n\t"
		
                "add %10, %13\n\t"
                "add %10, %14\n\t"
                "add %11, %13\n\t"
                "add %12, %14\n\t"

                "mov %13, 528(%6, %8)\n\t"
                "mov %14, 552(%6, %8)\n\t"
		
                "sub %11, %13\n\t"
                "sub %12, %14\n\t"
                "add %12, %13\n\t"
                "add %11, %14\n\t"
		
                "mov %13, 544(%6, %8)\n\t"
                "mov %14, 536(%6, %8)\n\t"
		
                "sub %10, %13\n\t"
                "sub %10, %14\n\t"
		
		
                //Berechne 16 Input - Base Adressen -- Zweite 4 Stück
                "mov  3072(%7,%5), %10\n\t"
                "add  2056(%7,%4), %10\n\t"
		
                "add %10, %13\n\t"
                "add %10, %14\n\t"
		
                "mov %13, 576(%6, %8)\n\t"
                "mov %14, 568(%6, %8)\n\t"
		
                "sub %12, %13\n\t"
                "sub %11, %14\n\t"
                "add %11, %13\n\t"
                "add %12, %14\n\t"
		
                "mov %13, 560(%6, %8)\n\t"
                "mov %14, 584(%6, %8)\n\t"
		
                "sub %10, %13\n\t"
                "sub %10, %14\n\t"
		
                //Berechne 16 Input - Base Adressen -- Dritte 4 Stück
                "mov  3080(%7,%5), %10\n\t"
                "add  2048(%7,%4), %10\n\t"
		
                "add %10, %13\n\t"
                "add %10, %14\n\t"
		
                "mov %13, 592(%6, %8)\n\t"
                "mov %14, 616(%6, %8)\n\t"
		
                "sub %11, %13\n\t"
                "sub %12, %14\n\t"
                "add %12, %13\n\t"
                "add %11, %14\n\t"
		
                "mov %13, 608(%6, %8)\n\t"
                "mov %14, 600(%6, %8)\n\t"
		
                "sub %10, %13\n\t"
                "sub %10, %14\n\t"
		
		
                //Berechne 16 Input - Base Adressen -- Vierte 4 Stück
                "mov  3080(%7,%5), %10\n\t"
                "add  2056(%7,%4), %10\n\t"
		
                "add %10, %13\n\t"
                "add %10, %14\n\t"
		
                "mov %13, 640(%6, %8)\n\t"
                "mov %14, 632(%6, %8)\n\t"
		
                "sub %12, %13\n\t"
                "sub %11, %14\n\t"
                "add %11, %13\n\t"
                "add %12, %14\n\t"
		
                "mov %13, 624(%6, %8)\n\t"
                "mov %14, 648(%6, %8)\n\t"
		
                //Berechne Target - Base - Adresse
/*                "mov %5, %1\n\t"
                "shl $2, %1\n\t"    //SHIFT hängt von Groesse oneDim ab!!! (N=2->1, N=4->2, N=8->3, usw...)
                "or  %4, %1\n\t"
                "shl $2, %1\n\t"    //SHIFT hängt von Groesse oneDim ab!!! (N=2->1, N=4->2, N=8->3, usw...)
                "or  %3, %1\n\t"
                "shl $2, %1\n\t"    //SHIFT hängt von Groesse oneDim ab!!! (N=2->1, N=4->2, N=8->3, usw...)
                "or  %2, %1\n\t"
                "shl $4, %1\n\t"
                "mov %1, 2576(%6, %9)\n\t"*/
		
		
                "add $128, %8\n\t"
                "add $8, %9\n\t"
                "sub $16, %2\n\t"
                "cmp 248(%6), %2\n\t"
              "jge SSE_FourierTrafoRearrangerWithFirstFFTstep_schleife4\n\t"
              "mov  280(%6), %2\n\t"

              "sub $16, %3\n\t"
              "cmp 240(%6), %3\n\t"
            "jge SSE_FourierTrafoRearrangerWithFirstFFTstep_schleife3\n\t"
            "mov  272(%6), %3\n\t"

            "sub $16, %4\n\t"
            "cmp 232(%6), %4\n\t"
          "jge SSE_FourierTrafoRearrangerWithFirstFFTstep_schleife2\n\t"
          "mov  264(%6), %4\n\t"

          "sub $16, %5\n\t"
          "cmp 224(%6), %5\n\t"
        "jge SSE_FourierTrafoRearrangerWithFirstFFTstep_schleife1\n\t"
		
		
	//Lade Daten für Main-Compute-Loop
        "mov 16(%6), %0\n\t"
        "mov 24(%6), %1\n\t"
        "mov  $1920, %2\n\t"
        "mov  $120, %3\n\t"
        "mov  %0, %8\n\t"
	
        "SSE_FourierTrafoRearrangerWithFirstFFTstep_MainComputeLoop: \n\t"
          "mov  2576(%6,%3), %5\n\t"
          "mov  8(%6), %4\n\t"
          "mov  %2, %7\n\t"
	  
	  
          "SSE_FourierTrafoRearrangerWithFirstFFTstep_InnerLoop: \n\t"
            "movapd 208(%6), %%xmm15\n\t"
		
            // LOAD + First STEP
            "mov  528(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm0\n\t"
            "mov  536(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm1\n\t"
            "addpd %%xmm0, %%xmm1\n\t"
            "mulpd %%xmm15, %%xmm0\n\t"
            "subpd %%xmm1, %%xmm0\n\t"

            "mov  544(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm2\n\t"
            "mov 552(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm3\n\t"
            "addpd %%xmm2, %%xmm3\n\t"
            "mulpd %%xmm15, %%xmm2\n\t"
            "subpd %%xmm3, %%xmm2\n\t"

            "mov 560(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm4\n\t"
            "mov 568(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm5\n\t"
            "addpd %%xmm4, %%xmm5\n\t"
            "mulpd %%xmm15, %%xmm4\n\t"
            "subpd %%xmm5, %%xmm4\n\t"

            "mov 576(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm6\n\t"
            "mov 584(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm7\n\t"
            "addpd %%xmm6, %%xmm7\n\t"
            "mulpd %%xmm15, %%xmm6\n\t"
            "subpd %%xmm7, %%xmm6\n\t"

            "mov 592(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm8\n\t"
            "mov 600(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm9\n\t"
            "addpd %%xmm8, %%xmm9\n\t"
            "mulpd %%xmm15, %%xmm8\n\t"
            "subpd %%xmm9, %%xmm8\n\t"

            "mov 608(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm10\n\t"
            "mov 616(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm11\n\t"
            "addpd %%xmm10, %%xmm11\n\t"
            "mulpd %%xmm15, %%xmm10\n\t"
            "subpd %%xmm11, %%xmm10\n\t"

            "mov 624(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm12\n\t"
            "mov 632(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm13\n\t"
            "addpd %%xmm12, %%xmm13\n\t"
            "mulpd %%xmm15, %%xmm12\n\t"
            "subpd %%xmm13, %%xmm12\n\t"

            "mov 640(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm14\n\t"
            "mov 648(%6, %2), %10\n\t"
            "movapd (%0, %10), %%xmm15\n\t"
            "addpd %%xmm14, %%xmm15\n\t"
            "mulpd 208(%6), %%xmm14\n\t"
            "subpd %%xmm15, %%xmm14\n\t"
	    
	    //Prefetching
            "mov  400(%6, %7), %11\n\t"
            "mov  408(%6, %7), %12\n\t"
            "prefetchnta (%8, %11)\n\t"
            "prefetchnta 64(%8, %11)\n\t"
            "prefetchnta (%8, %12)\n\t"
            "prefetchnta 64(%8, %12)\n\t"
            "add $16, %7\n\t"
	    
	    

	    // Second + third Step are external
           #include "ExtremeFFTasmBlock.h"

            // Forth STEP + SAVE
            "addpd %%xmm0, %%xmm8\n\t"
            "mulpd 208(%6), %%xmm0\n\t"
            "movapd %%xmm8,  2688(%1, %5)\n\t"  // <- Write
            "subpd %%xmm8, %%xmm0\n\t"
            "movapd %%xmm0,  10880(%1, %5)\n\t"  // <- Write
		  
            "addpd %%xmm1, %%xmm9\n\t"
            "mulpd 208(%6), %%xmm1\n\t"
            "movapd %%xmm9,  2560(%1, %5)\n\t"  // <- Write
            "subpd %%xmm9, %%xmm1\n\t"
            "movapd %%xmm1,  10752(%1, %5)\n\t"  // <- Write
		  
            "addpd %%xmm2, %%xmm10\n\t"
            "mulpd 208(%6), %%xmm2\n\t"
            "movapd %%xmm10,  2176(%1, %5)\n\t"  // <- Write
            "subpd %%xmm10, %%xmm2\n\t"
            "movapd %%xmm2,  10368(%1, %5)\n\t"  // <- Write

            "addpd %%xmm3, %%xmm11\n\t"
            "mulpd 208(%6), %%xmm3\n\t"
            "movapd %%xmm11,  2048(%1, %5)\n\t"  // <- Write
            "subpd %%xmm11, %%xmm3\n\t"
            "movapd %%xmm3,  10240(%1, %5)\n\t" // <- Writ

            "addpd %%xmm4, %%xmm12\n\t"
            "mulpd 208(%6), %%xmm4\n\t"
            "movapd %%xmm12, 640(%1, %5)\n\t"  // <- Write
            "subpd %%xmm12, %%xmm4\n\t"
            "movapd %%xmm4,  8832(%1, %5)\n\t"  // <- Write
		  
            "addpd %%xmm5, %%xmm13\n\t"
            "mulpd 208(%6), %%xmm5\n\t"
            "movapd %%xmm13, 512(%1, %5)\n\t"  // <- Write
            "subpd %%xmm13, %%xmm5\n\t"
            "movapd %%xmm5,  8704(%1, %5)\n\t"  // <- Write

            "addpd %%xmm6, %%xmm14\n\t"
            "mulpd 208(%6), %%xmm6\n\t"
            "movapd %%xmm14, 128(%1, %5)\n\t"  // <- Write
            "subpd %%xmm14, %%xmm6\n\t"
            "movapd %%xmm6, 8320(%1, %5)\n\t"  // <- Write

            "addpd %%xmm7, %%xmm15\n\t"
            "mulpd 208(%6), %%xmm7\n\t"
            "movapd %%xmm15, (%1, %5)\n\t"  // <- Write
            "subpd %%xmm15, %%xmm7\n\t"
            "movapd %%xmm7, 8192(%1, %5)\n\t"  // <- Write
	
		
            "add $16, %0\n\t"
            "add $16, %1\n\t"
            "sub $16, %4\n\t"
          "jg SSE_FourierTrafoRearrangerWithFirstFFTstep_InnerLoop\n\t"

          "sub 8(%6), %0\n\t"
          "sub 8(%6), %1\n\t"
          "sub $128, %2\n\t"
          "sub $8, %3\n\t"
        "jge SSE_FourierTrafoRearrangerWithFirstFFTstep_MainComputeLoop \n\t"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10), "=r" (Dummy11),  "=r" (Dummy12), "=r" (Dummy13), "=r" (Dummy14), "=r" (Dummy15)
          : "r" (para_Comp));

/*printf("reg0 %d\n",(long int) Dummy1);
printf("reg1 %d\n",(long int) Dummy2);
printf("reg2 %d\n",(long int) Dummy3);
printf("reg3 %d\n",(long int) Dummy4);
printf("reg4 %d\n",(long int) Dummy5);
printf("reg5 %d\n",(long int) Dummy6);
printf("reg6 %d\n",(long int) Dummy7);
printf("reg7 %d\n",(long int) Dummy8);
printf("reg8 %d\n",(long int) Dummy9);
printf("reg9 %d\n",(long int) Dummy10);
printf("reg10 %d\n",(long int) Dummy11);
printf("reg11 %d\n",(long int) Dummy12);
printf("reg12 %d\n",(long int) Dummy13);
printf("reg13 %d\n",(long int) Dummy14);
printf("reg14 %d\n",(long int) Dummy15);
printf("reg15 %d\n",(long int) para_Comp);

printBits((long int) Dummy12);
printBits((long int) Dummy13);
printBits((long int) Dummy14);
printBits((long int) Dummy15);
exit(0);*/

  //Always false
  if (para_int[2]!=para_int[2]) { 
    printf("%d\n",*Dummy1);
  }
}


/**** SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace ***
* Performs the application of the inverse of MdagM for static Phi in the Fourier-Space.
*
* Parameters:
* oneDim                   : Lattice size in one dimension
* input                    : Input Vector
* output                   : Output Vector
* sinP                     : Table, containing tilde p
* auxData                  : Table, containing eigenvalues and eigenvalues divided by sqrt(p^2)...
*
**/
void SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace(int oneDimL0, int oneDimL1, int oneDimL2, int oneDimL3, Complex* input, Complex* output, Complex* sinP, Complex* auxData, long int* PiPiPiPiIndices) {
  if (oneDimL0<=0) return;
  if (oneDimL1<=0) return;
  if (oneDimL2<=0) return;
  if (oneDimL3<=0) return;
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;
  int* Dummy11;
  int* Dummy12;
  int* Dummy13;
  Complex* para_Comp = ASMParameterData;
  long int* para_int = (long int*) para_Comp;
  para_int[0] = (long int) sinP;
  para_Comp[1] = Complex(0, 0);
  para_int[4] = (long int) input;
  para_int[5] = (long int) output;
  para_int[6] = (long int) PiPiPiPiIndices;
  para_int[7] = (long int) auxData;
  para_int[8] = (long int) xtraSize1*128*(oneDimL3+xtraSize3)*(oneDimL2+xtraSize2);
  para_int[9] = (long int) xtraSize2*128*(oneDimL3+xtraSize3);
  para_int[10] = (long int) xtraSize3*128;
  para_int[11] = 32*(oneDimL0-1);
  para_int[12] = 32*(oneDimL0/2);
  para_int[13] = 32*(oneDimL1-1);
  para_int[14] = 32*(oneDimL1/2);
  para_int[15] = 32*(oneDimL2-1);
  para_int[16] = 32*(oneDimL2/2);
  para_int[17] = 32*(oneDimL3-1);
  para_int[18] = 32*(oneDimL3/2);
  
  
  //Always false
  if (LogLevel>10000) { 
    printf("Compiler-Optimization-Prevention: Dummy = %d\n",(int) para_int[1]);  
  }

/*  int I;
  for (I=0; I<18; I++) {
    printf("para %d: %ld\n",I,para_int[I]);
  }
  printf("input: %ld\n",(long int)input);
  printf("output: %ld\n",(long int)output);
  printf("pipipi %ld\n",(long int)PiPiPiPiIndices);
  printf("para %ld\n",(long int)para_int);
 printf("sinP: %ld\n",(long int)sinP);
  printf("auxData: %ld\n",(long int)auxData);
  printf("------------------\n");
*/


__asm__ volatile(
        "mov %13, %6\n\t"
        "mov 0(%6), %12\n\t"
        "mov 88(%6), %2\n\t"
        "mov 104(%6), %3\n\t"
        "mov 120(%6), %4\n\t"
        "mov 136(%6), %5\n\t"
	
        "mov 32(%6), %0\n\t"
        "mov 40(%6), %1\n\t"
        "mov 48(%6), %9\n\t"
        "mov 56(%6), %7\n\t"
        "mov $0, %10\n\t"
        "mov $0, %11\n\t"
	
        "SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace_schleife1: \n\t"
          "movapd 16(%12,%2), %%xmm8\n\t"
          "add (%9,%2), %11\n\t"

          "SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace_schleife2: \n\t"
            "movapd 4096(%12,%3), %%xmm9\n\t"
            "add 4096(%9,%3), %11\n\t"
	  
            "SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace_schleife3: \n\t"
              "movapd 8208(%12,%4), %%xmm10\n\t"
              "add 8192(%9,%4), %11\n\t"
	    
              "SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace_schleife4: \n\t"
                "movapd 12288(%12,%5), %%xmm11\n\t"
                "add 12288(%9,%5), %11\n\t"

                "mov $1, %8\n\t"
                "SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace_InnerLoop:\n\t"

                  //Lade AuxData
                  "movapd 16(%7), %%xmm2\n\t"
                  "movhlps %%xmm2, %%xmm3\n\t"
                  "movlhps %%xmm2, %%xmm2\n\t"
                  "movlhps %%xmm3, %%xmm3\n\t"

		  //Collective Read
                  "movapd  0(%0,%10), %%xmm4\n\t"
                  "movapd 16(%0,%10), %%xmm5\n\t"
                  "movapd 32(%0,%10), %%xmm6\n\t"
                  "movapd 48(%0,%10), %%xmm7\n\t"
                  "prefetchnta 1024(%0)\n\t"
                  "prefetchnta 512(%7)\n\t"

                  "movapd %%xmm4,%%xmm12\n\t"
                  "movapd %%xmm5,%%xmm13\n\t"
                  "movapd %%xmm6,%%xmm14\n\t"
                  "movapd %%xmm7,%%xmm15\n\t"
	     

                  "mulpd %%xmm3,%%xmm4\n\t"
                  "mulpd %%xmm3,%%xmm5\n\t"
                  "mulpd %%xmm3,%%xmm6\n\t"
                  "mulpd %%xmm3,%%xmm7\n\t"
                  "mulpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm2,%%xmm13\n\t"
                  "mulpd %%xmm2,%%xmm14\n\t"
                  "mulpd %%xmm2,%%xmm15\n\t"

                  "movapd %%xmm4,%%xmm0\n\t"
                  "movapd %%xmm5,%%xmm1\n\t"
                  "movapd %%xmm6,%%xmm2\n\t"
                  "movapd %%xmm7,%%xmm3\n\t"

                  //Multiplikationen mit Realteil 1
                  "mulpd %%xmm11,%%xmm0\n\t"
                  "subpd %%xmm0,%%xmm14\n\t"
                  "mulpd %%xmm11,%%xmm1\n\t"
                  "addpd %%xmm1,%%xmm15\n\t"
                  "mulpd %%xmm11,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm11,%%xmm3\n\t"
                  "subpd %%xmm3,%%xmm13\n\t"


                  "movapd %%xmm4,%%xmm0\n\t"
                  "movapd %%xmm5,%%xmm1\n\t"
                  "movapd %%xmm6,%%xmm2\n\t"
                  "movapd %%xmm7,%%xmm3\n\t"
	      
                  //Multiplikationen mit Realteil 2
                  "mulpd %%xmm9,%%xmm0\n\t"
                  "subpd %%xmm0,%%xmm15\n\t"
                  "mulpd %%xmm9,%%xmm1\n\t"
                  "subpd %%xmm1,%%xmm14\n\t"
                  "mulpd %%xmm9,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm13\n\t"
                  "mulpd %%xmm9,%%xmm3\n\t"
                  "addpd %%xmm3,%%xmm12\n\t"

                  //Vertauschung Real- und Imaginärteil
                  "movhlps %%xmm4,%%xmm0\n\t"
                  "movlhps %%xmm4,%%xmm0\n\t"
                  "movhlps %%xmm5,%%xmm1\n\t"
                  "movlhps %%xmm5,%%xmm1\n\t"
                  "movhlps %%xmm6,%%xmm2\n\t"
                  "movlhps %%xmm6,%%xmm2\n\t"
                  "movhlps %%xmm7,%%xmm3\n\t"
                  "movlhps %%xmm7,%%xmm3\n\t"

                  "movapd %%xmm0,%%xmm4\n\t"
                  "movapd %%xmm1,%%xmm5\n\t"
                  "movapd %%xmm2,%%xmm6\n\t"
                  "movapd %%xmm3,%%xmm7\n\t"


                  //Multiplikationen mit Imaginärteil 1
                  "mulpd %%xmm8,%%xmm0\n\t"
                  "addpd %%xmm0,%%xmm14\n\t"
                  "mulpd %%xmm8,%%xmm1\n\t"
                  "addpd %%xmm1,%%xmm15\n\t"
                  "mulpd %%xmm8,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm8,%%xmm3\n\t"
                  "addpd %%xmm3,%%xmm13\n\t"


                  //Multiplikationen mit Imaginärteil 2
                  "mulpd %%xmm10,%%xmm4\n\t"
                  "subpd %%xmm4,%%xmm15\n\t"
                  "mulpd %%xmm10,%%xmm5\n\t"
                  "addpd %%xmm5,%%xmm14\n\t"
                  "mulpd %%xmm10,%%xmm6\n\t"
                  "addpd %%xmm6,%%xmm13\n\t"
                  "mulpd %%xmm10,%%xmm7\n\t"
                  "subpd %%xmm7,%%xmm12\n\t"


                  //Addiere Hauptdiagnonale
		  //Collective Read
                  "movapd  0(%0,%11), %%xmm4\n\t"
                  "movapd 16(%0,%11), %%xmm5\n\t"
                  "movapd 32(%0,%11), %%xmm6\n\t"
                  "movapd 48(%0,%11), %%xmm7\n\t"
                  "movapd 32(%7), %%xmm3\n\t"
                  "movapd   %%xmm4, %%xmm0\n\t"
                  "movapd   %%xmm5, %%xmm1\n\t"
                  "movapd   %%xmm6, %%xmm2\n\t"
                  "mulpd %%xmm3,%%xmm0\n\t"
                  "addpd %%xmm0,%%xmm12\n\t"
                  "movapd   %%xmm7, %%xmm0\n\t"
                  "mulpd %%xmm3,%%xmm1\n\t"
                  "addpd %%xmm1,%%xmm13\n\t"
                  "mulpd %%xmm3,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm14\n\t"
                  "mulpd %%xmm3,%%xmm0\n\t"
                  "addpd %%xmm0,%%xmm15\n\t"


                  //Collective Write
   	          #ifdef SmallLattice
                    "movapd %%xmm12, 0(%1,%11) \n\t"
                    "movapd %%xmm13,16(%1,%11) \n\t"
                    "movapd %%xmm14,32(%1,%11) \n\t"
                    "movapd %%xmm15,48(%1,%11) \n\t"
		  #else
                    "movntpd %%xmm12, 0(%1,%11) \n\t"
                    "movntpd %%xmm13,16(%1,%11) \n\t"
                    "movntpd %%xmm14,32(%1,%11) \n\t"
                    "movntpd %%xmm15,48(%1,%11) \n\t"
		  #endif
		

                  //Gleiches nochmal für die um pipipipi verschobenen moden
                  //Lade AuxData
                  "movapd 48(%7), %%xmm2\n\t"
                  "movhlps %%xmm2, %%xmm3\n\t"
                  "movlhps %%xmm2, %%xmm2\n\t"
                  "movlhps %%xmm3, %%xmm3\n\t"

		  //Daten sind geladen
                  "movapd %%xmm4,%%xmm12\n\t"
                  "movapd %%xmm5,%%xmm13\n\t"
                  "movapd %%xmm6,%%xmm14\n\t"
                  "movapd %%xmm7,%%xmm15\n\t"
	     

                  "mulpd %%xmm3,%%xmm4\n\t"
                  "mulpd %%xmm3,%%xmm5\n\t"
                  "mulpd %%xmm3,%%xmm6\n\t"
                  "mulpd %%xmm3,%%xmm7\n\t"
                  "mulpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm2,%%xmm13\n\t"
                  "mulpd %%xmm2,%%xmm14\n\t"
                  "mulpd %%xmm2,%%xmm15\n\t"

                  "movapd %%xmm4,%%xmm0\n\t"
                  "movapd %%xmm5,%%xmm1\n\t"
                  "movapd %%xmm6,%%xmm2\n\t"
                  "movapd %%xmm7,%%xmm3\n\t"

                  //Multiplikationen mit Realteil 1
                  "mulpd %%xmm11,%%xmm0\n\t"
                  "subpd %%xmm0,%%xmm14\n\t"
                  "mulpd %%xmm11,%%xmm1\n\t"
                  "addpd %%xmm1,%%xmm15\n\t"
                  "mulpd %%xmm11,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm11,%%xmm3\n\t"
                  "subpd %%xmm3,%%xmm13\n\t"


                  "movapd %%xmm4,%%xmm0\n\t"
                  "movapd %%xmm5,%%xmm1\n\t"
                  "movapd %%xmm6,%%xmm2\n\t"
                  "movapd %%xmm7,%%xmm3\n\t"
	      
                  //Multiplikationen mit Realteil 2
                  "mulpd %%xmm9,%%xmm0\n\t"
                  "subpd %%xmm0,%%xmm15\n\t"
                  "mulpd %%xmm9,%%xmm1\n\t"
                  "subpd %%xmm1,%%xmm14\n\t"
                  "mulpd %%xmm9,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm13\n\t"
                  "mulpd %%xmm9,%%xmm3\n\t"
                  "addpd %%xmm3,%%xmm12\n\t"

                  //Vertauschung Real- und Imaginärteil
                  "movhlps %%xmm4,%%xmm0\n\t"
                  "movlhps %%xmm4,%%xmm0\n\t"
                  "movhlps %%xmm5,%%xmm1\n\t"
                  "movlhps %%xmm5,%%xmm1\n\t"
                  "movhlps %%xmm6,%%xmm2\n\t"
                  "movlhps %%xmm6,%%xmm2\n\t"
                  "movhlps %%xmm7,%%xmm3\n\t"
                  "movlhps %%xmm7,%%xmm3\n\t"

                  "movapd %%xmm0,%%xmm4\n\t"
                  "movapd %%xmm1,%%xmm5\n\t"
                  "movapd %%xmm2,%%xmm6\n\t"
                  "movapd %%xmm3,%%xmm7\n\t"


                  //Multiplikationen mit Imaginärteil 1
                  "mulpd %%xmm8,%%xmm0\n\t"
                  "addpd %%xmm0,%%xmm14\n\t"
                  "mulpd %%xmm8,%%xmm1\n\t"
                  "addpd %%xmm1,%%xmm15\n\t"
                  "mulpd %%xmm8,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm8,%%xmm3\n\t"
                  "addpd %%xmm3,%%xmm13\n\t"


                  //Multiplikationen mit Imaginärteil 2
                  "mulpd %%xmm10,%%xmm4\n\t"
                  "subpd %%xmm4,%%xmm15\n\t"
                  "mulpd %%xmm10,%%xmm5\n\t"
                  "addpd %%xmm5,%%xmm14\n\t"
                  "mulpd %%xmm10,%%xmm6\n\t"
                  "addpd %%xmm6,%%xmm13\n\t"
                  "mulpd %%xmm10,%%xmm7\n\t"
                  "subpd %%xmm7,%%xmm12\n\t"


                  //Addiere Hauptdiagnonale
		  //Collective Read
                  "movapd  0(%0,%10), %%xmm4\n\t"
                  "movapd 16(%0,%10), %%xmm5\n\t"
                  "movapd 32(%0,%10), %%xmm6\n\t"
                  "movapd 48(%0,%10), %%xmm7\n\t"
                  "movapd (%7), %%xmm3\n\t"
                  "mulpd %%xmm3,%%xmm4\n\t"
                  "addpd %%xmm4,%%xmm12\n\t"
                  "mulpd %%xmm3,%%xmm5\n\t"
                  "addpd %%xmm5,%%xmm13\n\t"
                  "mulpd %%xmm3,%%xmm6\n\t"
                  "addpd %%xmm6,%%xmm14\n\t"
                  "mulpd %%xmm3,%%xmm7\n\t"
                  "addpd %%xmm7,%%xmm15\n\t"


                  //Collective Write
   	          #ifdef SmallLattice
                    "movapd %%xmm12, 0(%1,%10) \n\t"
                    "movapd %%xmm13,16(%1,%10) \n\t"
                    "movapd %%xmm14,32(%1,%10) \n\t"
                    "movapd %%xmm15,48(%1,%10) \n\t"
		  #else
                    "movntpd %%xmm12, 0(%1,%10) \n\t"
                    "movntpd %%xmm13,16(%1,%10) \n\t"
                    "movntpd %%xmm14,32(%1,%10) \n\t"
                    "movntpd %%xmm15,48(%1,%10) \n\t"
		  #endif


                  "add $64, %10\n\t"
                  "add $64, %11\n\t"
                  "dec %8\n\t"
                "jge SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace_InnerLoop\n\t"
		  
                "sub $128, %11\n\t"
                "sub 12288(%9,%5), %11\n\t"
                "add $64, %7\n\t"
                "sub $32,%5\n\t"
              "jge SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace_schleife4\n\t"
              "sub 8192(%9,%4), %11\n\t"
       "add 80(%6), %10\n\t"
	      
              "mov 136(%6), %5\n\t"
              "sub $32,%4\n\t"
            "jge SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace_schleife3\n\t"
            "sub 4096(%9,%3), %11\n\t"
       "add 72(%6), %10\n\t"
    
            "mov 120(%6), %4\n\t"
            "sub $32,%3\n\t"
          "jge SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace_schleife2\n\t"
          "sub (%9,%2), %11\n\t"
      "add 64(%6), %10\n\t"
	      
          "mov 104(%6), %3\n\t"
          "sub $32,%2\n\t"
          "cmp 96(%6),%2\n\t"
        "jge SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace_schleife1"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10), "=r" (Dummy11), "=r" (Dummy12), "=r" (Dummy13)
          : "r" (para_int) );

  //Always false
  if (para_int[11]!=32*(oneDimL0-1)) { 
    int* xxx1 = (int*)sinP;
    int* xxx2 = (int*)auxData;
    printf("%d\n", *xxx1);
    printf("%d\n", *xxx2);
    printf("%d\n",*Dummy1);
    printf("%d\n",*Dummy2);
    printf("%d\n",*Dummy3);
    printf("%d\n",*Dummy4);
    printf("%d\n",*Dummy5);
    printf("%d\n",*Dummy6);
    printf("%d\n",*Dummy7);
    printf("%d\n",*Dummy8);
    printf("%d\n",*Dummy9);
    printf("%d\n",*Dummy10);
    printf("%d\n",*Dummy11);
    printf("%d\n",*Dummy12);
    printf("%d\n",*Dummy13);
  }
  

/*    printf("Reg0: %ld\n", (long int) Dummy1);
    printf("Reg1: %ld\n", (long int) Dummy2);
    printf("Reg2: %ld\n", (long int) Dummy3);
    printf("Reg3: %ld\n", (long int) Dummy4);
    printf("Reg4: %ld\n", (long int) Dummy5);
    printf("Reg5: %ld\n", (long int) Dummy6);
    printf("Reg6: %ld\n", (long int) Dummy7);
    printf("Reg7: %ld\n", (long int) Dummy8);
    printf("Reg8: %ld\n", (long int) Dummy9);
    printf("Reg9: %ld\n", (long int) Dummy10);
    printf("Reg10: %ld\n", (long int) Dummy11);
    printf("Reg11: %ld\n", (long int) Dummy12);
    printf("Reg12: %ld\n", (long int) Dummy13);
//    exit(0);*/
}


/**** SSE_Performf_YB_2rho_AndCopyToOutput ***
* Applies yNB^+/2rho -1 to x and stores result in interim and stores yN B^+ applied on x in output.
* ATTENTION: Watch out for factors 2rho!!!
*
* Parameters:
* Lto4                  : Lattice size in one dimension
* yN                    : Parameter yN
* twoRho                : 2 * Parameter Rho
* x                     : Input Vector
* phi                   : Phi - Feld
* interim               : Zwischenspeichervektor
* output                : Output Vector
*
**/
void SSE_Performf_YB_2rho_AndCopyToOutput(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double twoRho, Complex* x, double* phi, Complex* interim, Complex* output) {
  if ((L0<=0) ||(L1<=0) ||(L2<=0) ||(L3<=0)) return;
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;
  
//  phi = (double*) createSuperAlignedComplex(6);
  
  Complex* para_Comp = createSuperAlignedComplex(10);
  long int* para_int = (long int*) para_Comp;
  para_int[0] = 0;
  para_int[1] = 0;
  para_Comp[1] = Complex(1.0 / twoRho, 1.0 / twoRho);
  para_int[4] = (long int) x;
  para_int[5] = (long int) output;
  para_int[6] = (long int) interim;
  para_int[7] = (long int) &(phi[0]);  
  para_Comp[4] = Complex(yN, yN);
  para_int[10] = (long int) 0x8000000000000000;
  para_int[11] = (long int) 0;
  para_int[12] = L0;
  para_int[13] = L1;
  para_int[14] = L2;
  para_int[15] = L3;
  para_int[16] = (long int) xtrS1*128*(L3+xtrS3)*(L2+xtrS2);
  para_int[17] = (long int) xtrS2*128*(L3+xtrS3);
  para_int[18] = (long int) xtrS3*128;
  
  

  //Always false
  if (LogLevel>10000) { 
    printf("Compiler-Optimization-Prevention: Phi = %f\n",phi[0]);  
  }

__asm__ volatile("mov 96(%10), %4\n\t"
        "mov 104(%10), %5\n\t"
        "mov 112(%10), %6\n\t"
        "mov 120(%10), %7\n\t"
        "mov 32(%10), %0\n\t"
        "mov 40(%10), %1\n\t"
        "mov 48(%10), %2\n\t"
        "mov 56(%10), %3\n\t"
	
        "SSE_Performf_YB_2rho_AndCopyToOutput_schleife1: \n\t"
          "SSE_Performf_YB_2rho_AndCopyToOutput_schleife2: \n\t"
            "SSE_Performf_YB_2rho_AndCopyToOutput_schleife3: \n\t"
              "SSE_Performf_YB_2rho_AndCopyToOutput_schleife4: \n\t"
    	        //Laden der Phi - Werte
                "movapd 64(%10), %%xmm0\n\t"
                "movapd 80(%10), %%xmm1\n\t"    //Konstante zum Negieren der ersten Double-Zahl
                "movapd 0(%3), %%xmm4\n\t"
                "movapd 16(%3), %%xmm6\n\t"
                "prefetchnta 1024(%3)\n\t"
                "mulpd %%xmm0, %%xmm4\n\t"
                "mulpd %%xmm0, %%xmm6\n\t"
                "movhlps %%xmm4, %%xmm5\n\t"
                "movhlps %%xmm6, %%xmm7\n\t"
                "movlhps %%xmm4, %%xmm4\n\t"
                "movlhps %%xmm5, %%xmm5\n\t"
                "movlhps %%xmm6, %%xmm6\n\t"
                "movlhps %%xmm7, %%xmm7\n\t"
                "xorpd %%xmm1, %%xmm5\n\t"
                "xorpd %%xmm1, %%xmm7\n\t"
	  
                //Lade ungerade Input-Indices
                "movapd  0(%0), %%xmm8\n\t"
                "movapd 32(%0), %%xmm9\n\t"
                "movapd 64(%0), %%xmm10\n\t"
                "movapd 96(%0), %%xmm11\n\t"
                "prefetchnta 1024(%0)\n\t"

                "movhlps %%xmm8, %%xmm0\n\t"
                "movlhps %%xmm8, %%xmm0\n\t"
                "movhlps %%xmm9, %%xmm1\n\t"
                "movlhps %%xmm9, %%xmm1\n\t"
                "movhlps %%xmm10, %%xmm2\n\t"
                "movlhps %%xmm10, %%xmm2\n\t"
                "movhlps %%xmm11, %%xmm3\n\t"
                "movlhps %%xmm11, %%xmm3\n\t"
	  
                "movapd  %%xmm8, %%xmm12\n\t"
                "movapd  %%xmm9, %%xmm13\n\t"
                "movapd  %%xmm10, %%xmm14\n\t"
                "movapd  %%xmm11, %%xmm15\n\t"
	  
 	        //Phi_0
                "mulpd %%xmm4, %%xmm12\n\t"
                "mulpd %%xmm4, %%xmm13\n\t"
                "mulpd %%xmm4, %%xmm14\n\t"
                "mulpd %%xmm4, %%xmm15\n\t"

	        //Phi_2
                "mulpd %%xmm6, %%xmm10\n\t"
                "subpd %%xmm10, %%xmm12\n\t"
                "mulpd %%xmm6, %%xmm11\n\t"
                "addpd %%xmm11, %%xmm13\n\t"
                "mulpd %%xmm6, %%xmm8\n\t"
                "addpd %%xmm8, %%xmm14\n\t"
                "mulpd %%xmm6, %%xmm9\n\t"
                "subpd %%xmm9, %%xmm15\n\t"


                "movapd  %%xmm0, %%xmm8\n\t"
                "movapd  %%xmm1, %%xmm9\n\t"
                "movapd  %%xmm2, %%xmm10\n\t"
                "movapd  %%xmm3, %%xmm11\n\t"

                //Phi_1
                "mulpd %%xmm5, %%xmm2\n\t"
                "subpd %%xmm2, %%xmm12\n\t"
                "mulpd %%xmm5, %%xmm3\n\t"
                "addpd %%xmm3, %%xmm13\n\t"
                "mulpd %%xmm5, %%xmm0\n\t"
                "subpd %%xmm0, %%xmm14\n\t"
                "mulpd %%xmm5, %%xmm1\n\t"
                "addpd %%xmm1, %%xmm15\n\t"
	  
	        //Phi_3
                "mulpd %%xmm7, %%xmm8\n\t"
                "subpd %%xmm8, %%xmm12\n\t"
                "mulpd %%xmm7, %%xmm9\n\t"
                "addpd %%xmm9, %%xmm13\n\t"
                "mulpd %%xmm7, %%xmm10\n\t"
                "addpd %%xmm10, %%xmm14\n\t"
                "mulpd %%xmm7, %%xmm11\n\t"
                "subpd %%xmm11, %%xmm15\n\t"
	  
	        //Collective Write nach Output
                #ifdef SmallLattice
                  "movapd  %%xmm12,   0(%1)\n\t"
                  "movapd  %%xmm13,  32(%1)\n\t"
                  "movapd  %%xmm14, 64(%1)\n\t"
                  "movapd  %%xmm15, 96(%1)\n\t"
  	        #else
                  "movntpd  %%xmm12,   0(%1)\n\t"
                  "movntpd  %%xmm13,  32(%1)\n\t"
                  "movntpd  %%xmm14, 64(%1)\n\t"
                  "movntpd  %%xmm15, 96(%1)\n\t"
  	        #endif

   	        //Division durch 2rho und sub x
                "movapd 16(%10), %%xmm0\n\t"
                "mulpd %%xmm0, %%xmm12\n\t"	  
                "subpd  0(%0), %%xmm12\n\t"
                "mulpd %%xmm0, %%xmm13\n\t"	  
                "subpd 32(%0), %%xmm13\n\t"
                "mulpd %%xmm0, %%xmm14\n\t"	  
                "subpd 64(%0), %%xmm14\n\t"
                "mulpd %%xmm0, %%xmm15\n\t"	  
                "subpd 96(%0), %%xmm15\n\t"
	  
	        //Collective Write nach interim
                #ifdef SmallLattice
                  "movapd  %%xmm12,   0(%2)\n\t"
                  "movapd  %%xmm13,  32(%2)\n\t"
                  "movapd  %%xmm14, 64(%2)\n\t"
                  "movapd  %%xmm15, 96(%2)\n\t"
	        #else
                  "movntpd  %%xmm12,   0(%2)\n\t"
                  "movntpd  %%xmm13,  32(%2)\n\t"
                  "movntpd  %%xmm14, 64(%2)\n\t"
                  "movntpd  %%xmm15, 96(%2)\n\t"
	        #endif
	  
	  
  	        //Lade gerade Input-Indices
                "movapd  16(%0), %%xmm8\n\t"
                "movapd  48(%0), %%xmm9\n\t"
                "movapd  80(%0), %%xmm10\n\t"
                "movapd 112(%0), %%xmm11\n\t"
                "prefetchnta 1088(%0)\n\t"

                "movhlps %%xmm8, %%xmm0\n\t"
                "movlhps %%xmm8, %%xmm0\n\t"
                "movhlps %%xmm9, %%xmm1\n\t"
                "movlhps %%xmm9, %%xmm1\n\t"
                "movhlps %%xmm10, %%xmm2\n\t"
                "movlhps %%xmm10, %%xmm2\n\t"
                "movhlps %%xmm11, %%xmm3\n\t"
                "movlhps %%xmm11, %%xmm3\n\t"
	  
                "movapd  %%xmm8, %%xmm12\n\t"
                "movapd  %%xmm9, %%xmm13\n\t"
                "movapd  %%xmm10, %%xmm14\n\t"
                "movapd  %%xmm11, %%xmm15\n\t"
	  
	        //Phi_0
                "mulpd %%xmm4, %%xmm12\n\t"
                "mulpd %%xmm4, %%xmm13\n\t"
                "mulpd %%xmm4, %%xmm14\n\t"
                "mulpd %%xmm4, %%xmm15\n\t"

      	        //Phi_2
                "mulpd %%xmm6, %%xmm10\n\t"
                "subpd %%xmm10, %%xmm12\n\t"
                "mulpd %%xmm6, %%xmm11\n\t"
                "addpd %%xmm11, %%xmm13\n\t"
                "mulpd %%xmm6, %%xmm8\n\t"
                "addpd %%xmm8, %%xmm14\n\t"
                "mulpd %%xmm6, %%xmm9\n\t"
                "subpd %%xmm9, %%xmm15\n\t"


                "movapd  %%xmm0, %%xmm8\n\t"
                "movapd  %%xmm1, %%xmm9\n\t"
                "movapd  %%xmm2, %%xmm10\n\t"
                "movapd  %%xmm3, %%xmm11\n\t"

	        //Phi_1
                "mulpd %%xmm5, %%xmm2\n\t"
                "subpd %%xmm2, %%xmm12\n\t"
                "mulpd %%xmm5, %%xmm3\n\t"
                "addpd %%xmm3, %%xmm13\n\t"
                "mulpd %%xmm5, %%xmm0\n\t"
                "subpd %%xmm0, %%xmm14\n\t"
                "mulpd %%xmm5, %%xmm1\n\t"
                "addpd %%xmm1, %%xmm15\n\t"
	  
	        //Phi_3
                "mulpd %%xmm7, %%xmm8\n\t"
                "subpd %%xmm8, %%xmm12\n\t"
                "mulpd %%xmm7, %%xmm9\n\t"
                "addpd %%xmm9, %%xmm13\n\t"
                "mulpd %%xmm7, %%xmm10\n\t"
                "addpd %%xmm10, %%xmm14\n\t"
                "mulpd %%xmm7, %%xmm11\n\t"
                "subpd %%xmm11, %%xmm15\n\t"
	  
 	        //Collective Write nach Output
                #ifdef SmallLattice
                  "movapd  %%xmm12,  16(%1)\n\t"
                  "movapd  %%xmm13,  48(%1)\n\t"
                  "movapd  %%xmm14,  80(%1)\n\t"
                  "movapd  %%xmm15, 112(%1)\n\t"
  	        #else 
                  "movntpd  %%xmm12,  16(%1)\n\t"
                  "movntpd  %%xmm13,  48(%1)\n\t"
                  "movntpd  %%xmm14,  80(%1)\n\t"
                  "movntpd  %%xmm15, 112(%1)\n\t"
 	        #endif


	        //Division durch 2rho und sub x
                "movapd 16(%10), %%xmm0\n\t"
                "mulpd %%xmm0, %%xmm12\n\t"	  
                "subpd  16(%0), %%xmm12\n\t"
                "mulpd %%xmm0, %%xmm13\n\t"	  
                "subpd  48(%0), %%xmm13\n\t"
                "mulpd %%xmm0, %%xmm14\n\t"	  
                "subpd  80(%0), %%xmm14\n\t"
                "mulpd %%xmm0, %%xmm15\n\t"	  
                "subpd 112(%0), %%xmm15\n\t"
	  
 	        //Collective Write nach interim
     	        #ifdef SmallLattice
                  "movapd  %%xmm12,  16(%2)\n\t"
                  "movapd  %%xmm13,  48(%2)\n\t"
                  "movapd  %%xmm14,  80(%2)\n\t"
                  "movapd  %%xmm15, 112(%2)\n\t"	  
	        #else
                  "movntpd  %%xmm12,  16(%2)\n\t"
                  "movntpd  %%xmm13,  48(%2)\n\t"
                  "movntpd  %%xmm14,  80(%2)\n\t"
                 "movntpd  %%xmm15, 112(%2)\n\t"	  
	        #endif
	  
		
                "add $128, %0\n\t"
                "add $128, %1\n\t"
                "add $128, %2\n\t"
                "add $32, %3\n\t"
                "dec %7\n\t"
              "jg SSE_Performf_YB_2rho_AndCopyToOutput_schleife4\n\t"
              "mov 120(%10), %7\n\t"
              "add 144(%10), %0\n\t"
              "add 144(%10), %1\n\t"
              "add 144(%10), %2\n\t"
	      
              "dec %6\n\t"
            "jg SSE_Performf_YB_2rho_AndCopyToOutput_schleife3\n\t"
            "mov 112(%10), %6\n\t"
            "add 136(%10), %0\n\t"
            "add 136(%10), %1\n\t"
            "add 136(%10), %2\n\t"

            "dec %5\n\t"
          "jg SSE_Performf_YB_2rho_AndCopyToOutput_schleife2\n\t"
          "mov 104(%10), %5\n\t"
          "add 128(%10), %0\n\t"
          "add 128(%10), %1\n\t"
          "add 128(%10), %2\n\t"

          "dec %4\n\t"
        "jg SSE_Performf_YB_2rho_AndCopyToOutput_schleife1\n\t"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
          : "r" (para_Comp) );


  //Always false
  if (para_int[12]!=(L0)) { 
    printf("%d\n",*Dummy1);
  }
  destroySuperAlignedComplex(para_Comp);
}


/**** SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace ***
* Performs the application of the Dirac-Operator in the Fourier-Space.
*
*            work[0] = work[0] + auxData[countAux].y*(Complex(p3,p0)*input[count+2] + Complex(p1,-p2)*input[count+3]);
*            work[1] = work[1] + auxData[countAux].y*(Complex(p1,p2)*input[count+2] + Complex(-p3,p0)*input[count+3]);
*            work[2] = work[2] + auxData[countAux].y*(Complex(-p3,p0)*input[count+0] + Complex(-p1,p2)*input[count+1]);
*            work[3] = work[3] + auxData[countAux].y*(Complex(-p1,-p2)*input[count+0] + Complex(p3,p0)*input[count+1]);
*
* Parameters:
* oneDim                   : Lattice size in one dimension
* input                    : Input Vector
* output                   : Output Vector
* sinP                     : Table, containing tilde p
* auxData                  : Table, containing eigenvalues and eigenvalues divided by sqrt(p^2)...
*
**/
void SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(int oneDimL0, int oneDimL1, int oneDimL2, int oneDimL3, Complex* input, Complex* output, Complex* sinP, Complex* auxData) {
  if (oneDimL0<=0) return;
  if (oneDimL1<=0) return;
  if (oneDimL2<=0) return;
  if (oneDimL3<=0) return;
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;
  Complex* para_Comp = new Complex[7];
  long int* para_int = (long int*) para_Comp;
  para_Comp[1] = Complex(0, 0);
  para_int[4] = (long int) input;
  para_int[5] = (long int) output;
  para_int[6] = (long int) xtraSize1*128*(oneDimL3+xtraSize3)*(oneDimL2+xtraSize2);
  para_int[7] = (long int) xtraSize2*128*(oneDimL3+xtraSize3);
  para_int[8] = (long int) xtraSize3*128;
  para_int[9] = 32*(oneDimL0-1);
  para_int[10] = 32*(oneDimL1-1);
  para_int[11] = 32*(oneDimL2-1);
  para_int[12] = 32*(oneDimL3-1);


__asm__ volatile(
        "mov 72(%10), %2\n\t"
        "mov 80(%10), %3\n\t"
        "mov 88(%10), %4\n\t"
        "mov 96(%10), %5\n\t"
        "mov 96(%10), %6\n\t"
        "mov %11, %7\n\t"
	
        "mov 32(%10), %0\n\t"
        "mov 40(%10), %1\n\t"
	
        "SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace_schleife1: \n\t"
          "mov %2, %8\n\t"
          "add %12, %8\n\t"
          "movapd 16(%8), %%xmm8\n\t"

          "SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace_schleife2: \n\t"
            "mov %3, %8\n\t"
            "add %12, %8\n\t"
            "movapd 4096(%8), %%xmm9\n\t"
	  
            "SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace_schleife3: \n\t"
              "mov %4, %8\n\t"
              "mov %6, %9\n\t"
              "add %12, %8\n\t"
              "add %12, %9\n\t"
              "movapd 8208(%8), %%xmm10\n\t"
	    
              "SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace_schleife4: \n\t"
                "movapd 12288(%9), %%xmm11\n\t"


                  //Erster Dirac-Operator
                  "movapd (%7), %%xmm2\n\t"
                  "movhlps %%xmm2, %%xmm3\n\t"
                  "movlhps %%xmm2, %%xmm2\n\t"
                  "movlhps %%xmm3, %%xmm3\n\t"

		  //Collective Read
                  "movapd  0(%0), %%xmm4\n\t"
                  "movapd 16(%0), %%xmm5\n\t"
                  "movapd 32(%0), %%xmm6\n\t"
                  "movapd 48(%0), %%xmm7\n\t"
                  "prefetchnta 1024(%0)\n\t"
                  "prefetchnta 512(%7)\n\t"

                  "movapd %%xmm4,%%xmm12\n\t"
                  "movapd %%xmm5,%%xmm13\n\t"
                  "movapd %%xmm6,%%xmm14\n\t"
                  "movapd %%xmm7,%%xmm15\n\t"
	     

                  "mulpd %%xmm3,%%xmm4\n\t"
                  "mulpd %%xmm3,%%xmm5\n\t"
                  "mulpd %%xmm3,%%xmm6\n\t"
                  "mulpd %%xmm3,%%xmm7\n\t"
                  "mulpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm2,%%xmm13\n\t"
                  "mulpd %%xmm2,%%xmm14\n\t"
                  "mulpd %%xmm2,%%xmm15\n\t"

                  "movapd %%xmm4,%%xmm0\n\t"
                  "movapd %%xmm5,%%xmm1\n\t"
                  "movapd %%xmm6,%%xmm2\n\t"
                  "movapd %%xmm7,%%xmm3\n\t"

                  //Multiplikationen mit Realteil 1
                  "mulpd %%xmm11,%%xmm0\n\t"
                  "subpd %%xmm0,%%xmm14\n\t"
                  "mulpd %%xmm11,%%xmm1\n\t"
                  "addpd %%xmm1,%%xmm15\n\t"
                  "mulpd %%xmm11,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm11,%%xmm3\n\t"
                  "subpd %%xmm3,%%xmm13\n\t"


                  "movapd %%xmm4,%%xmm0\n\t"
                  "movapd %%xmm5,%%xmm1\n\t"
                  "movapd %%xmm6,%%xmm2\n\t"
                  "movapd %%xmm7,%%xmm3\n\t"
	      
                  //Multiplikationen mit Realteil 2
                  "mulpd %%xmm9,%%xmm0\n\t"
                  "subpd %%xmm0,%%xmm15\n\t"
                  "mulpd %%xmm9,%%xmm1\n\t"
                  "subpd %%xmm1,%%xmm14\n\t"
                  "mulpd %%xmm9,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm13\n\t"
                  "mulpd %%xmm9,%%xmm3\n\t"
                  "addpd %%xmm3,%%xmm12\n\t"

                  //Vertauschung Real- und Imaginärteil
                  "movhlps %%xmm4,%%xmm0\n\t"
                  "movlhps %%xmm4,%%xmm0\n\t"
                  "movhlps %%xmm5,%%xmm1\n\t"
                  "movlhps %%xmm5,%%xmm1\n\t"
                  "movhlps %%xmm6,%%xmm2\n\t"
                  "movlhps %%xmm6,%%xmm2\n\t"
                  "movhlps %%xmm7,%%xmm3\n\t"
                  "movlhps %%xmm7,%%xmm3\n\t"

                  "movapd %%xmm0,%%xmm4\n\t"
                  "movapd %%xmm1,%%xmm5\n\t"
                  "movapd %%xmm2,%%xmm6\n\t"
                  "movapd %%xmm3,%%xmm7\n\t"


                  //Multiplikationen mit Imaginärteil 1
                  "mulpd %%xmm8,%%xmm0\n\t"
                  "addpd %%xmm0,%%xmm14\n\t"
                  "mulpd %%xmm8,%%xmm1\n\t"
                  "addpd %%xmm1,%%xmm15\n\t"
                  "mulpd %%xmm8,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm8,%%xmm3\n\t"
                  "addpd %%xmm3,%%xmm13\n\t"


                  //Multiplikationen mit Imaginärteil 2
                  "mulpd %%xmm10,%%xmm4\n\t"
                  "subpd %%xmm4,%%xmm15\n\t"
                  "mulpd %%xmm10,%%xmm5\n\t"
                  "addpd %%xmm5,%%xmm14\n\t"
                  "mulpd %%xmm10,%%xmm6\n\t"
                  "addpd %%xmm6,%%xmm13\n\t"
                  "mulpd %%xmm10,%%xmm7\n\t"
                  "subpd %%xmm7,%%xmm12\n\t"

                  //Collective Write
   	          #ifdef SmallLattice
                    "movapd %%xmm12, 0(%1) \n\t"
                    "movapd %%xmm13,16(%1) \n\t"
                    "movapd %%xmm14,32(%1) \n\t"
                    "movapd %%xmm15,48(%1) \n\t"
		  #else
                    "movntpd %%xmm12, 0(%1) \n\t"
                    "movntpd %%xmm13,16(%1) \n\t"
                    "movntpd %%xmm14,32(%1) \n\t"
                    "movntpd %%xmm15,48(%1) \n\t"
		  #endif


                  //Zweiter Dirac-Operator
                  "movapd (%7), %%xmm2\n\t"
                  "movhlps %%xmm2, %%xmm3\n\t"
                  "movlhps %%xmm2, %%xmm2\n\t"
                  "movlhps %%xmm3, %%xmm3\n\t"

		  //Collective Read
                  "movapd  64(%0), %%xmm4\n\t"
                  "movapd  80(%0), %%xmm5\n\t"
                  "movapd  96(%0), %%xmm6\n\t"
                  "movapd 112(%0), %%xmm7\n\t"
                  "prefetchnta 1088(%0)\n\t"

                  "movapd %%xmm4,%%xmm12\n\t"
                  "movapd %%xmm5,%%xmm13\n\t"
                  "movapd %%xmm6,%%xmm14\n\t"
                  "movapd %%xmm7,%%xmm15\n\t"
	     

                  "mulpd %%xmm3,%%xmm4\n\t"
                  "mulpd %%xmm3,%%xmm5\n\t"
                  "mulpd %%xmm3,%%xmm6\n\t"
                  "mulpd %%xmm3,%%xmm7\n\t"
                  "mulpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm2,%%xmm13\n\t"
                  "mulpd %%xmm2,%%xmm14\n\t"
                  "mulpd %%xmm2,%%xmm15\n\t"

                  "movapd %%xmm4,%%xmm0\n\t"
                  "movapd %%xmm5,%%xmm1\n\t"
                  "movapd %%xmm6,%%xmm2\n\t"
                  "movapd %%xmm7,%%xmm3\n\t"

                  //Multiplikationen mit Realteil 1
                  "mulpd %%xmm11,%%xmm0\n\t"
                  "subpd %%xmm0,%%xmm14\n\t"
                  "mulpd %%xmm11,%%xmm1\n\t"
                  "addpd %%xmm1,%%xmm15\n\t"
                  "mulpd %%xmm11,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm11,%%xmm3\n\t"
                  "subpd %%xmm3,%%xmm13\n\t"


                  "movapd %%xmm4,%%xmm0\n\t"
                  "movapd %%xmm5,%%xmm1\n\t"
                  "movapd %%xmm6,%%xmm2\n\t"
                  "movapd %%xmm7,%%xmm3\n\t"
	      
                  //Multiplikationen mit Realteil 2
                  "mulpd %%xmm9,%%xmm0\n\t"
                  "subpd %%xmm0,%%xmm15\n\t"
                  "mulpd %%xmm9,%%xmm1\n\t"
                  "subpd %%xmm1,%%xmm14\n\t"
                  "mulpd %%xmm9,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm13\n\t"
                  "mulpd %%xmm9,%%xmm3\n\t"
                  "addpd %%xmm3,%%xmm12\n\t"

                  //Vertauschung Real- und Imaginärteil
                  "movhlps %%xmm4,%%xmm0\n\t"
                  "movlhps %%xmm4,%%xmm0\n\t"
                  "movhlps %%xmm5,%%xmm1\n\t"
                  "movlhps %%xmm5,%%xmm1\n\t"
                  "movhlps %%xmm6,%%xmm2\n\t"
                  "movlhps %%xmm6,%%xmm2\n\t"
                  "movhlps %%xmm7,%%xmm3\n\t"
                  "movlhps %%xmm7,%%xmm3\n\t"

                  "movapd %%xmm0,%%xmm4\n\t"
                  "movapd %%xmm1,%%xmm5\n\t"
                  "movapd %%xmm2,%%xmm6\n\t"
                  "movapd %%xmm3,%%xmm7\n\t"


                  //Multiplikationen mit Imaginärteil 1
                  "mulpd %%xmm8,%%xmm0\n\t"
                  "addpd %%xmm0,%%xmm14\n\t"
                  "mulpd %%xmm8,%%xmm1\n\t"
                  "addpd %%xmm1,%%xmm15\n\t"
                  "mulpd %%xmm8,%%xmm2\n\t"
                  "addpd %%xmm2,%%xmm12\n\t"
                  "mulpd %%xmm8,%%xmm3\n\t"
                  "addpd %%xmm3,%%xmm13\n\t"


                  //Multiplikationen mit Imaginärteil 2
                  "mulpd %%xmm10,%%xmm4\n\t"
                  "subpd %%xmm4,%%xmm15\n\t"
                  "mulpd %%xmm10,%%xmm5\n\t"
                  "addpd %%xmm5,%%xmm14\n\t"
                  "mulpd %%xmm10,%%xmm6\n\t"
                  "addpd %%xmm6,%%xmm13\n\t"
                  "mulpd %%xmm10,%%xmm7\n\t"
                  "subpd %%xmm7,%%xmm12\n\t"

                  //Collective Write
   	          #ifdef SmallLattice
                    "movapd %%xmm12, 64(%1) \n\t"
                    "movapd %%xmm13, 80(%1) \n\t"
                    "movapd %%xmm14, 96(%1) \n\t"
                    "movapd %%xmm15,112(%1) \n\t"
 		  #else
                    "movntpd %%xmm12, 64(%1) \n\t"
                    "movntpd %%xmm13, 80(%1) \n\t"
                    "movntpd %%xmm14, 96(%1) \n\t"
                    "movntpd %%xmm15,112(%1) \n\t"
	 	  #endif
		

                "sub $32,%9\n\t"
                "add $128, %0\n\t"
                "add $128, %1\n\t"
                "add $16, %7\n\t"
                "sub $32,%5\n\t"
              "jge SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace_schleife4\n\t"
              "mov %6, %5\n\t"
              "add 64(%10), %0\n\t"
              "add 64(%10), %1\n\t"
              "sub $32,%4\n\t"
            "jge SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace_schleife3\n\t"
            "mov 88(%10), %4\n\t"
            "add 56(%10), %0\n\t"
            "add 56(%10), %1\n\t"
            "sub $32,%3\n\t"
          "jge SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace_schleife2\n\t"
          "mov 80(%10), %3\n\t"
          "add 48(%10), %0\n\t"
          "add 48(%10), %1\n\t"
          "sub $32,%2\n\t"
        "jge SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace_schleife1"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
          : "r" (para_Comp),  "r" (auxData), "r" (sinP) );


  //Always false
  if (para_int[9]!=32*(oneDimL0-1)) { 
    int* xxx1 = (int*)sinP;
    int* xxx2 = (int*)auxData;
    printf("%d\n", *xxx1);
    printf("%d\n", *xxx2);
    printf("%d\n",*Dummy1);
  }

  delete[] para_Comp;
}


/**** SSE_ComplexVectorAdditionSPECIAL2 ***
* Performs a complex vector subtraction including a real multiplication with a constant
* acting on the result.
* y = beta*( x - y )
* Parameters:
* N                   : Vector length (must be dividable by 8)
* beta                : Real constant
* x                   : Pointer to complex vector x
* y                   : Pointer to complex vector y
*
**/
void SSE_ComplexVectorAdditionSPECIAL2(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double beta, Complex* x, Complex* y) {
  if ((L0<=0) ||(L1<=0) ||(L2<=0) ||(L3<=0)) return;
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;

  Complex* bet = new Complex[1];
  bet[0].x = beta;
  bet[0].y = beta;

  long int* para_N = new long int[10];
  para_N[0] = (long int) x;
  para_N[1] = (long int) y;
  para_N[2] = L0;
  para_N[3] = L1;
  para_N[4] = L2;
  para_N[5] = L3;
  para_N[6] = (long int) xtrS1*128*(L3+xtrS3)*(L2+xtrS2);
  para_N[7] = (long int) xtrS2*128*(L3+xtrS3);
  para_N[8] = (long int) xtrS3*128;
  para_N[9] = (long int) (bet);


__asm__ volatile("mov %10, %0\n\t"
        "mov 0(%0), %1\n\t"
        "mov 8(%0), %2\n\t"
        "mov 16(%0), %3\n\t"
        "mov 24(%0), %4\n\t"
        "mov 32(%0), %5\n\t"
        "mov 40(%0), %6\n\t"
        "mov 72(%0), %7\n\t"
        "movupd (%7), %%xmm2\n\t"

        "SSE_ComplexVectorAdditionSPECIAL2_schleife1: \n\t"
          "SSE_ComplexVectorAdditionSPECIAL2_schleife2: \n\t"
            "SSE_ComplexVectorAdditionSPECIAL2_schleife3: \n\t"
              "SSE_ComplexVectorAdditionSPECIAL2_schleife4: \n\t"
                "movapd   (%1), %%xmm8\n\t"
                "movapd 16(%1), %%xmm9\n\t"
                "movapd 32(%1), %%xmm10\n\t"
                "movapd 48(%1), %%xmm11\n\t"
                "movapd 64(%1), %%xmm12\n\t"
                "movapd 80(%1), %%xmm13\n\t"
                "movapd 96(%1), %%xmm14\n\t"
                "movapd 112(%1), %%xmm15\n\t"
                "prefetchnta 1024(%1)\n\t"
                "prefetchnta 1088(%1)\n\t"
	  

                "subpd  (%2), %%xmm8\n\t"
                "subpd  16(%2), %%xmm9\n\t"
                "subpd  32(%2), %%xmm10\n\t"
                "subpd  48(%2), %%xmm11\n\t"
                "subpd  64(%2), %%xmm12\n\t"
                "subpd  80(%2), %%xmm13\n\t"
                "subpd  96(%2), %%xmm14\n\t"
                "subpd  112(%2), %%xmm15\n\t"
                "prefetchnta 1024(%2)\n\t"
                "prefetchnta 1088(%2)\n\t"

                "mulpd %%xmm2, %%xmm8\n\t"
                "mulpd %%xmm2, %%xmm9\n\t"
                "mulpd %%xmm2, %%xmm10\n\t"
                "mulpd %%xmm2, %%xmm11\n\t"
                "mulpd %%xmm2, %%xmm12\n\t"
                "mulpd %%xmm2, %%xmm13\n\t"
                "mulpd %%xmm2, %%xmm14\n\t"
                "mulpd %%xmm2, %%xmm15\n\t"


                //Collective Write
	        #ifdef SmallLattice
                  "movapd %%xmm8,   (%2) \n\t"
                  "movapd %%xmm9,  16(%2) \n\t"
                  "movapd %%xmm10, 32(%2) \n\t"
                  "movapd %%xmm11, 48(%2) \n\t"
                  "movapd %%xmm12, 64(%2) \n\t"
                  "movapd %%xmm13, 80(%2) \n\t"
                  "movapd %%xmm14, 96(%2) \n\t"
                  "movapd %%xmm15,112(%2) \n\t"
    	        #else
                  "movntpd %%xmm8,   (%2) \n\t"
                  "movntpd %%xmm9,  16(%2) \n\t"
                  "movntpd %%xmm10, 32(%2) \n\t"
                  "movntpd %%xmm11, 48(%2) \n\t"
                  "movntpd %%xmm12, 64(%2) \n\t"
                  "movntpd %%xmm13, 80(%2) \n\t"
                  "movntpd %%xmm14, 96(%2) \n\t"
                  "movntpd %%xmm15,112(%2) \n\t"
 	        #endif
//   "prefetchnta 1024(%2)\n\t"

                "add $128, %1\n\t"
                "add $128, %2\n\t"

                "dec %6\n\t"
              "jg SSE_ComplexVectorAdditionSPECIAL2_schleife4\n\t"
	      "mov 40(%0), %6\n\t"
              "add 64(%0), %1\n\t"
              "add 64(%0), %2\n\t"
	      
              "dec %5\n\t"
            "jg SSE_ComplexVectorAdditionSPECIAL2_schleife3\n\t"
	    "mov 32(%0), %5\n\t"
            "add 56(%0), %1\n\t"
            "add 56(%0), %2\n\t"
	      
            "dec %4\n\t"
          "jg SSE_ComplexVectorAdditionSPECIAL2_schleife2\n\t"
  	  "mov 24(%0), %4\n\t"
          "add 48(%0), %1\n\t"
          "add 48(%0), %2\n\t"
	      
          "dec %3\n\t"
        "jg SSE_ComplexVectorAdditionSPECIAL2_schleife1\n\t"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
          : "r" (para_N));

  //Always false
  if (para_N[2]!=(L0)) printf("%d\n",*Dummy1);

  delete[] bet;
  delete[] para_N;
}


/**** SSE_ComplexVectorAdditionSPECIAL1 ***
* Performs a complex vector addition including a real multiplication with a constant
* acting on the vector y.
* y = beta*y + x
* Parameters:
* N                   : Vector length (must be dividable by 8)
* beta                : Real constant
* x                   : Pointer to complex vector x
* y                   : Pointer to complex vector y
*
**/
void SSE_ComplexVectorAdditionSPECIAL1(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double beta, Complex* x, Complex* y) {
  if ((L0<=0) ||(L1<=0) ||(L2<=0) ||(L3<=0)) return;
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;

  Complex* bet = new Complex[1];
  bet[0].x = beta;
  bet[0].y = beta;
  
  long int* para_N = new long int[10];
  para_N[0] = (long int) x;
  para_N[1] = (long int) y;
  para_N[2] = L0;
  para_N[3] = L1;
  para_N[4] = L2;
  para_N[5] = L3;
  para_N[6] = (long int) xtrS1*128*(L3+xtrS3)*(L2+xtrS2);
  para_N[7] = (long int) xtrS2*128*(L3+xtrS3);
  para_N[8] = (long int) xtrS3*128;
  para_N[9] = (long int) (bet);


__asm__ volatile("mov %10, %0\n\t"
        "mov 0(%0), %1\n\t"
        "mov 8(%0), %2\n\t"
        "mov 16(%0), %3\n\t"
        "mov 24(%0), %4\n\t"
        "mov 32(%0), %5\n\t"
        "mov 40(%0), %6\n\t"
        "mov 72(%0), %7\n\t"
        "movupd (%7), %%xmm2\n\t"

        "SSE_ComplexVectorAdditionSPECIAL1_schleife1: \n\t"
          "SSE_ComplexVectorAdditionSPECIAL1_schleife2: \n\t"
            "SSE_ComplexVectorAdditionSPECIAL1_schleife3: \n\t"
              "SSE_ComplexVectorAdditionSPECIAL1_schleife4: \n\t"
                "movapd   (%2), %%xmm8\n\t"
                "movapd 16(%2), %%xmm9\n\t"
                "movapd 32(%2), %%xmm10\n\t"
                "movapd 48(%2), %%xmm11\n\t"
                "movapd 64(%2), %%xmm12\n\t"
                "movapd 80(%2), %%xmm13\n\t"
                "movapd 96(%2), %%xmm14\n\t"
                "movapd 112(%2), %%xmm15\n\t"
                "prefetchnta 1024(%1)\n\t"
                "prefetchnta 1024(%2)\n\t"
                "prefetchnta 1088(%1)\n\t"
                "prefetchnta 1088(%2)\n\t"
	  

                "mulpd %%xmm2, %%xmm8\n\t"
                "mulpd %%xmm2, %%xmm9\n\t"
                "mulpd %%xmm2, %%xmm10\n\t"
                "mulpd %%xmm2, %%xmm11\n\t"
                "mulpd %%xmm2, %%xmm12\n\t"
                "mulpd %%xmm2, %%xmm13\n\t"
                "mulpd %%xmm2, %%xmm14\n\t"
                "mulpd %%xmm2, %%xmm15\n\t"

                "addpd  (%1), %%xmm8\n\t"
                "addpd  16(%1), %%xmm9\n\t"
                "addpd  32(%1), %%xmm10\n\t"
                "addpd  48(%1), %%xmm11\n\t"
                "addpd  64(%1), %%xmm12\n\t"
                "addpd  80(%1), %%xmm13\n\t"
                "addpd  96(%1), %%xmm14\n\t"
                "addpd  112(%1), %%xmm15\n\t"

                //Collective Write
	        #ifdef SmallLattice
                  "movapd %%xmm8,   (%2) \n\t"
                  "movapd %%xmm9,  16(%2) \n\t"
                  "movapd %%xmm10, 32(%2) \n\t"
                  "movapd %%xmm11, 48(%2) \n\t"
                  "movapd %%xmm12, 64(%2) \n\t"
                  "movapd %%xmm13, 80(%2) \n\t"
                  "movapd %%xmm14, 96(%2) \n\t"
                  "movapd %%xmm15,112(%2) \n\t"
	        #else
                  "movntpd %%xmm8,   (%2) \n\t"
                  "movntpd %%xmm9,  16(%2) \n\t"
                  "movntpd %%xmm10, 32(%2) \n\t"
                  "movntpd %%xmm11, 48(%2) \n\t"
                  "movntpd %%xmm12, 64(%2) \n\t"
                  "movntpd %%xmm13, 80(%2) \n\t"
                  "movntpd %%xmm14, 96(%2) \n\t"
                  "movntpd %%xmm15,112(%2) \n\t"
   	        #endif

                "add $128, %1\n\t"
                "add $128, %2\n\t"

                "dec %6\n\t"
              "jg SSE_ComplexVectorAdditionSPECIAL1_schleife4\n\t"
	      "mov 40(%0), %6\n\t"
              "add 64(%0), %1\n\t"
              "add 64(%0), %2\n\t"
	      
              "dec %5\n\t"
            "jg SSE_ComplexVectorAdditionSPECIAL1_schleife3\n\t"
	    "mov 32(%0), %5\n\t"
            "add 56(%0), %1\n\t"
            "add 56(%0), %2\n\t"
	      
            "dec %4\n\t"
          "jg SSE_ComplexVectorAdditionSPECIAL1_schleife2\n\t"
  	  "mov 24(%0), %4\n\t"
          "add 48(%0), %1\n\t"
          "add 48(%0), %2\n\t"
	      
          "dec %3\n\t"
        "jg SSE_ComplexVectorAdditionSPECIAL1_schleife1\n\t"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
          : "r" (para_N));


  //Always false
  if (para_N[2]!=(L0)) printf("%d\n",*Dummy1);

  delete[] bet;
  delete[] para_N;
}


/**** SSE_ComplexVectorAdditionWithSquaredNorm ***
* Performs a complex vector addition including a complex multiplication with a constant
* and calculates the squared norm of the result.
* y = y + alpha*x
* Parameters:
* N                   : Vector length (must be dividable by 4)
* alpha               : Complex constant
* x                   : Pointer to complex vector x
* y                   : Pointer to complex vector y
* SqrNorm             : Squared norm of result
*
**/
void SSE_ComplexVectorAdditionWithSquaredNorm(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, Complex& alpha, Complex* x, Complex* y, double& SqrNorm) {
  SqrNorm = 0;
  if ((L0<=0) ||(L1<=0) ||(L2<=0) ||(L3<=0)) return;
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;
  
  Complex* al = new Complex[2];
  al[0].x = alpha.x;
  al[0].y = alpha.x;
  al[1].x = -alpha.y;
  al[1].y = alpha.y;
  Complex* res = new Complex[1];
  res->x = 0;
  res->y = 0;
  
  long int* para_N = new long int[12];
  para_N[0] = (long int) x;
  para_N[1] = (long int) y;
  para_N[2] = L0;
  para_N[3] = L1;
  para_N[4] = L2;
  para_N[5] = L3*2;
  para_N[6] = (long int) xtrS1*128*(L3+xtrS3)*(L2+xtrS2);
  para_N[7] = (long int) xtrS2*128*(L3+xtrS3);
  para_N[8] = (long int) xtrS3*128;
  para_N[9] = (long int) (al);
  para_N[10] = (long int) (res);
  
/*  printf("para = %d\n",(long int) para_N);
  printf("x = %d\n",(long int) x);
  printf("y = %d\n",(long int) y);*/

__asm__ volatile("mov %10, %0\n\t"
        "mov 16(%0), %3\n\t"
	"mov 24(%0), %4\n\t"
	"mov 32(%0), %5\n\t"
	"mov 40(%0), %6\n\t"
        "mov 0(%0), %1\n\t"
        "mov 8(%0), %2\n\t"
	"mov 72(%0), %7\n\t"
	"mov 80(%0), %8\n\t"
	
        "movupd   (%7), %%xmm2\n\t"
        "movupd 16(%7), %%xmm3\n\t"
        "movupd (%8), %%xmm0\n\t"
        "movupd %%xmm0, %%xmm1\n\t"

        "SSE_ComplexVectorAdditionWithSquaredNorm_schleife1: \n\t"
          "SSE_ComplexVectorAdditionWithSquaredNorm_schleife2: \n\t"
            "SSE_ComplexVectorAdditionWithSquaredNorm_schleife3: \n\t"
              "SSE_ComplexVectorAdditionWithSquaredNorm_schleife4: \n\t"
                "movapd   (%1), %%xmm8\n\t"
                "movapd 16(%1), %%xmm9\n\t"
                "movapd 32(%1), %%xmm10\n\t"
                "movapd 48(%1), %%xmm11\n\t"
	  
                "movapd   (%2), %%xmm12\n\t"
                "movapd 16(%2), %%xmm13\n\t"
                "movapd 32(%2), %%xmm14\n\t"
                "movapd 48(%2), %%xmm15\n\t"
                "prefetchnta 1024(%1)\n\t"
                "prefetchnta 1024(%2)\n\t"
	  
                "movhlps %%xmm8,  %%xmm4\n\t"
                "movlhps %%xmm8,  %%xmm4\n\t"
                "movhlps %%xmm9,  %%xmm5\n\t"
                "movlhps %%xmm9,  %%xmm5\n\t"
                "movhlps %%xmm10, %%xmm6\n\t"
                "movlhps %%xmm10, %%xmm6\n\t"
                "movhlps %%xmm11, %%xmm7\n\t"
                "movlhps %%xmm11, %%xmm7\n\t"
	  
	  
                "mulpd %%xmm2, %%xmm8\n\t"
                "addpd %%xmm8, %%xmm12\n\t"
                "mulpd %%xmm2, %%xmm9\n\t"
                "addpd %%xmm9, %%xmm13\n\t"
                "mulpd %%xmm2, %%xmm10\n\t"
                "addpd %%xmm10, %%xmm14\n\t"
                "mulpd %%xmm2, %%xmm11\n\t"
                "addpd %%xmm11, %%xmm15\n\t"
	  
                "mulpd %%xmm3, %%xmm4\n\t"
                "addpd %%xmm4, %%xmm12\n\t"
                "mulpd %%xmm3, %%xmm5\n\t"
                "addpd %%xmm5, %%xmm13\n\t"
                "mulpd %%xmm3, %%xmm6\n\t"
                "addpd %%xmm6, %%xmm14\n\t"
                "mulpd %%xmm3, %%xmm7\n\t"
                "addpd %%xmm7, %%xmm15\n\t"
	  
                //Collective Write
                #ifdef SmallLattice
                  "movapd %%xmm12,   (%2) \n\t"
                  "movapd %%xmm13, 16(%2) \n\t"
                  "movapd %%xmm14, 32(%2) \n\t"
                  "movapd %%xmm15, 48(%2) \n\t"
   	        #else
                  "movntpd %%xmm12,   (%2) \n\t"
                  "movntpd %%xmm13, 16(%2) \n\t"
                  "movntpd %%xmm14, 32(%2) \n\t"
                  "movntpd %%xmm15, 48(%2) \n\t"
	        #endif

                //Skalarprodukt
                "mulpd %%xmm12, %%xmm12\n\t"
                "addpd %%xmm12, %%xmm0\n\t"
                "mulpd %%xmm13, %%xmm13\n\t"
                "addpd %%xmm13, %%xmm1\n\t"
                "mulpd %%xmm14, %%xmm14\n\t"
                "addpd %%xmm14, %%xmm0\n\t"
                "mulpd %%xmm15, %%xmm15\n\t"
                "addpd %%xmm15, %%xmm1\n\t"
	  
                "add $64, %1\n\t"
                "add $64, %2\n\t"
                "dec %6\n\t"
              "jg SSE_ComplexVectorAdditionWithSquaredNorm_schleife4\n\t"
	      "mov 40(%0), %6\n\t"
              "add 64(%0), %1\n\t"
              "add 64(%0), %2\n\t"
	      
              "dec %5\n\t"
            "jg SSE_ComplexVectorAdditionWithSquaredNorm_schleife3\n\t"
	    "mov 32(%0), %5\n\t"
            "add 56(%0), %1\n\t"
            "add 56(%0), %2\n\t"
	      
            "dec %4\n\t"
          "jg SSE_ComplexVectorAdditionWithSquaredNorm_schleife2\n\t"
  	  "mov 24(%0), %4\n\t"
          "add 48(%0), %1\n\t"
          "add 48(%0), %2\n\t"
	      
          "dec %3\n\t"
        "jg SSE_ComplexVectorAdditionWithSquaredNorm_schleife1\n\t"
        "addpd %%xmm1, %%xmm0\n\t"
        "movhlps %%xmm0, %%xmm1\n\t"
        "addpd %%xmm1, %%xmm0\n\t"
        "movupd %%xmm0, (%8)"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
          : "r" (para_N));


/*printf("reg0 %d\n",(long int) Dummy1);
printf("reg1 %d\n",(long int) Dummy2);
printf("reg2 %d\n",(long int) Dummy3);
printf("reg3 %d\n",(long int) Dummy4);
printf("reg4 %d\n",(long int) Dummy5);
printf("reg5 %d\n",(long int) Dummy6);
printf("reg6 %d\n",(long int) Dummy7);
printf("reg7 %d\n",(long int) Dummy8);
printf("reg8 %d\n",(long int) Dummy9);
printf("reg9 %d\n",(long int) Dummy10);
exit(0);*/


  //Always false
  if (para_N[2]!=(L0)) printf("%d\n",*Dummy1);

  SqrNorm = res->x;

  delete[] al;
  delete[] para_N;
  delete[] res;
}


/**** SSE_ComplexVectorAddition ***
* Performs a complex vector addition including a complex multiplication with a constant.
* y = y + alpha*x
* Parameters:
* N                   : Vector length (must be dividable by 4)
* alpha               : Complex constant
* x                   : Pointer to complex vector x
* y                   : Pointer to complex vector y
*
**/
void SSE_ComplexVectorAddition(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, Complex& alpha, Complex* x, Complex* y) {
  if ((L0<=0) ||(L1<=0) ||(L2<=0) ||(L3<=0)) return;
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;
  Complex* al = new Complex[2];
  al[0].x = alpha.x;
  al[0].y = alpha.x;
  al[1].x = -alpha.y;
  al[1].y = alpha.y;


  long int* para_N = new long int[10];
  para_N[0] = (long int) x;
  para_N[1] = (long int) y;
  para_N[2] = L0;
  para_N[3] = L1;
  para_N[4] = L2;
  para_N[5] = L3*2;
  para_N[6] = (long int) (xtrS1*128*(L3+xtrS3)*(L2+xtrS2));
  para_N[7] = (long int) (xtrS2*128*(L3+xtrS3));
  para_N[8] = (long int) (xtrS3*128);
  para_N[9] = (long int) (al);
  

__asm__ volatile("mov %10, %0\n\t"
        "mov 0(%0), %1\n\t"
        "mov 8(%0), %2\n\t"
        "mov 16(%0), %3\n\t"
        "mov 24(%0), %4\n\t"
        "mov 32(%0), %5\n\t"
        "mov 40(%0), %6\n\t"
        "mov 72(%0), %7\n\t"
        "movupd   (%7), %%xmm2\n\t"
        "movupd 16(%7), %%xmm3\n\t"

        "SSE_ComplexVectorAddition_schleife1: \n\t"
          "SSE_ComplexVectorAddition_schleife2: \n\t"
            "SSE_ComplexVectorAddition_schleife3: \n\t"
              "SSE_ComplexVectorAddition_schleife4: \n\t"
                "movapd   (%1), %%xmm8\n\t"
                "movapd 16(%1), %%xmm9\n\t"
                "movapd 32(%1), %%xmm10\n\t"
                "movapd 48(%1), %%xmm11\n\t"
     	  
                "movapd   (%2), %%xmm12\n\t"
                "movapd 16(%2), %%xmm13\n\t"
                "movapd 32(%2), %%xmm14\n\t"
                "movapd 48(%2), %%xmm15\n\t"
                "prefetchnta 1024(%1)\n\t"
                "prefetchnta 1024(%2)\n\t"
	  
                "movhlps %%xmm8,  %%xmm4\n\t"
                "movlhps %%xmm8,  %%xmm4\n\t"
                "movhlps %%xmm9,  %%xmm5\n\t"
                "movlhps %%xmm9,  %%xmm5\n\t"
                "movhlps %%xmm10, %%xmm6\n\t"
                "movlhps %%xmm10, %%xmm6\n\t"
                "movhlps %%xmm11, %%xmm7\n\t"
                "movlhps %%xmm11, %%xmm7\n\t"
	  
                "mulpd %%xmm2, %%xmm8\n\t"
                "addpd %%xmm8, %%xmm12\n\t"
                "mulpd %%xmm2, %%xmm9\n\t"
                "addpd %%xmm9, %%xmm13\n\t"
                "mulpd %%xmm2, %%xmm10\n\t"
                "addpd %%xmm10, %%xmm14\n\t"
                "mulpd %%xmm2, %%xmm11\n\t"
                "addpd %%xmm11, %%xmm15\n\t"
	  
                "mulpd %%xmm3, %%xmm4\n\t"
                "addpd %%xmm4, %%xmm12\n\t"
                "mulpd %%xmm3, %%xmm5\n\t"
                "addpd %%xmm5, %%xmm13\n\t"
                "mulpd %%xmm3, %%xmm6\n\t"
                "addpd %%xmm6, %%xmm14\n\t"
                "mulpd %%xmm3, %%xmm7\n\t"
                "addpd %%xmm7, %%xmm15\n\t"
	  
                //Collective Write
	        #ifdef SmallLattice
                  "movapd %%xmm12,   (%2) \n\t"
                  "movapd %%xmm13, 16(%2) \n\t"
                  "movapd %%xmm14, 32(%2) \n\t"
                  "movapd %%xmm15, 48(%2) \n\t"
    	        #else
                  "movntpd %%xmm12,   (%2) \n\t"
                  "movntpd %%xmm13, 16(%2) \n\t"
                  "movntpd %%xmm14, 32(%2) \n\t"
                  "movntpd %%xmm15, 48(%2) \n\t"
  	        #endif
	  
                "add $64, %1\n\t"
                "add $64, %2\n\t"

                "dec %6\n\t"
              "jg SSE_ComplexVectorAddition_schleife4\n\t"
	      "mov 40(%0), %6\n\t"
              "add 64(%0), %1\n\t"
              "add 64(%0), %2\n\t"
	      
              "dec %5\n\t"
            "jg SSE_ComplexVectorAddition_schleife3\n\t"
	    "mov 32(%0), %5\n\t"
            "add 56(%0), %1\n\t"
            "add 56(%0), %2\n\t"
	      
            "dec %4\n\t"
          "jg SSE_ComplexVectorAddition_schleife2\n\t"
  	  "mov 24(%0), %4\n\t"
          "add 48(%0), %1\n\t"
          "add 48(%0), %2\n\t"
	      
          "dec %3\n\t"
        "jg SSE_ComplexVectorAddition_schleife1\n\t"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
          : "r" (para_N));

  //Always false
  if (para_N[2]!=(L0)) printf("%d\n",*Dummy1);

  delete[] al;
  delete[] para_N;
}


/**** SSE_ComplexScalarProduct ***
* Calcs the scalar product of to vectors.
* Parameters:
* N                   : Vector length (must be dividable by 4)
* v1                  : Pointer to complex vector 1
* v2                  : Pointer to complex vector 2
* r                   : Complex results
*
**/
void SSE_ComplexScalarProduct(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, Complex* v1, Complex* v2, Complex& r) {
  r = Complex(0,0);
  if ((L0<=0) ||(L1<=0) ||(L2<=0) ||(L3<=0)) return;
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;

  Complex* res = new Complex[1];
  res->x = 0;
  res->y = 0;

  long int* para_N = new long int[12];
  para_N[0] = (long int) v1;
  para_N[1] = (long int) v2;
  para_N[2] = L0;
  para_N[3] = L1;
  para_N[4] = L2;
  para_N[5] = L3*2;
  para_N[6] = (long int) xtrS1*128*(L3+xtrS3)*(L2+xtrS2);
  para_N[7] = (long int) xtrS2*128*(L3+xtrS3);
  para_N[8] = (long int) xtrS3*128;
  para_N[9] = (long int) (res);

__asm__ volatile("mov %10, %0\n\t"
        "mov 0(%0), %1\n\t"
        "mov 8(%0), %2\n\t"
        "mov 16(%0), %3\n\t"
        "mov 24(%0), %4\n\t"
        "mov 32(%0), %5\n\t"
        "mov 40(%0), %6\n\t"
        "mov 72(%0), %7\n\t"
        "movupd (%7), %%xmm7\n\t"
        "movapd %%xmm7, %%xmm12\n\t"
        "movapd %%xmm7, %%xmm13\n\t"

        "SSE_ComplexScalarProduct_schleife1: \n\t"
          "SSE_ComplexScalarProduct_schleife2: \n\t"
            "SSE_ComplexScalarProduct_schleife3: \n\t"
              "SSE_ComplexScalarProduct_schleife4: \n\t"
                "movapd   (%1), %%xmm0\n\t"
                "movapd 16(%1), %%xmm1\n\t"
                "movapd 32(%1), %%xmm2\n\t"
                "movapd 48(%1), %%xmm3\n\t"
	  
                "movapd   (%2), %%xmm4\n\t"
                "movapd 16(%2), %%xmm5\n\t"
                "movapd 32(%2), %%xmm6\n\t"
                "movapd 48(%2), %%xmm7\n\t"
                "prefetchnta 1024(%1)\n\t"
                "prefetchnta 1024(%2)\n\t"
	  
                "movhlps  %%xmm4, %%xmm8\n\t"
                "movlhps  %%xmm4, %%xmm8\n\t"
                "movhlps  %%xmm5, %%xmm9\n\t"
                "movlhps  %%xmm5, %%xmm9\n\t"
                "movhlps  %%xmm6, %%xmm10\n\t"
                "movlhps  %%xmm6, %%xmm10\n\t"
                "movhlps  %%xmm7, %%xmm11\n\t"
                "movlhps  %%xmm7, %%xmm11\n\t"
	  
                "mulpd %%xmm0, %%xmm4\n\t"
                "mulpd %%xmm1, %%xmm5\n\t"
                "addpd %%xmm4, %%xmm5\n\t"	  
                "mulpd %%xmm2, %%xmm6\n\t"
                "mulpd %%xmm3, %%xmm7\n\t"
                "addpd %%xmm6, %%xmm7\n\t"	  
                "addpd %%xmm5, %%xmm7\n\t"	  
	  
                "mulpd %%xmm0, %%xmm8\n\t"
                "mulpd %%xmm1, %%xmm9\n\t"
                "addpd %%xmm8, %%xmm9\n\t"	  
                "mulpd %%xmm2, %%xmm10\n\t"
                "mulpd %%xmm3, %%xmm11\n\t"
                "addpd %%xmm10, %%xmm11\n\t"	  
                "addpd %%xmm9, %%xmm11\n\t"	  
	  
                "addpd %%xmm7, %%xmm12\n\t"	  
                "addpd %%xmm11, %%xmm13\n\t"	  

                "add $64, %1\n\t"
                "add $64, %2\n\t"

                "dec %6\n\t"
              "jg SSE_ComplexScalarProduct_schleife4\n\t"
	      "mov 40(%0), %6\n\t"
              "add 64(%0), %1\n\t"
              "add 64(%0), %2\n\t"
	      
              "dec %5\n\t"
            "jg SSE_ComplexScalarProduct_schleife3\n\t"
	    "mov 32(%0), %5\n\t"
            "add 56(%0), %1\n\t"
            "add 56(%0), %2\n\t"
	      
            "dec %4\n\t"
          "jg SSE_ComplexScalarProduct_schleife2\n\t"
  	  "mov 24(%0), %4\n\t"
          "add 48(%0), %1\n\t"
          "add 48(%0), %2\n\t"
	      
          "dec %3\n\t"
        "jg SSE_ComplexScalarProduct_schleife1\n\t"

        "movhlps %%xmm12, %%xmm0\n\t"
        "movhlps %%xmm13, %%xmm1\n\t"
        "addpd %%xmm0, %%xmm12\n\t"
        "subpd %%xmm1, %%xmm13\n\t"
        "shufpd $0, %%xmm13, %%xmm12\n\t"
        "movupd %%xmm12, (%7)\n\t"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
          : "r" (para_N));



  //Always false
  if (para_N[2]!=(L0)) printf("%d\n",*Dummy1);

  r.x = res[0].x;
  r.y = res[0].y;
  delete[] res;
  delete[] para_N;
}


/**** SSE_ComplexSqrNorm ***
* Calcs the squared norm of a complex vector.
* Parameters:
* N                   : Vector length (must be dividable by 8)
* v                   : Pointer to complex vector
* r                   : Complex results
*
**/
void SSE_ComplexSquareNorm(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, Complex* v, double& r) {
  r=0;
  if ((L0<=0) ||(L1<=0) ||(L2<=0) ||(L3<=0)) return;
  int* Dummy1;
  int* Dummy2;
  int* Dummy3;
  int* Dummy4;
  int* Dummy5;
  int* Dummy6;
  int* Dummy7;
  int* Dummy8;
  int* Dummy9;
  int* Dummy10;

  Complex* res = new Complex[1];
  res->x = 0;
  res->y = 0;

  long int* para_N = new long int[10];
  para_N[0] = (long int) v;
  para_N[1] = (long int) 0;
  para_N[2] = L0;
  para_N[3] = L1;
  para_N[4] = L2;
  para_N[5] = L3;
  para_N[6] = (long int) xtrS1*128*(L3+xtrS3)*(L2+xtrS2);
  para_N[7] = (long int) xtrS2*128*(L3+xtrS3);
  para_N[8] = (long int) xtrS3*128;
  para_N[9] = (long int) (res);



__asm__ volatile("mov %10, %0\n\t"
        "mov 0(%0), %1\n\t"
//        "mov 8(%0), %2\n\t"
        "mov 16(%0), %3\n\t"
        "mov 24(%0), %4\n\t"
        "mov 32(%0), %5\n\t"
        "mov 40(%0), %6\n\t"
        "mov 72(%0), %7\n\t"
        "movupd (%7), %%xmm10\n\t"

        "SSE_ComplexSquareNorm_schleife1: \n\t"
          "SSE_ComplexSquareNorm_schleife2: \n\t"
            "SSE_ComplexSquareNorm_schleife3: \n\t"
              "SSE_ComplexSquareNorm_schleife4: \n\t"
                "movapd    (%1), %%xmm0\n\t"
                "movapd  16(%1), %%xmm1\n\t"
                "movapd  32(%1), %%xmm2\n\t"
                "movapd  48(%1), %%xmm3\n\t"
                "movapd  64(%1), %%xmm4\n\t"
                "movapd  80(%1), %%xmm5\n\t"
                "movapd  96(%1), %%xmm6\n\t"
                "movapd 112(%1), %%xmm7\n\t"
                "prefetchnta 1024(%1)\n\t"
                "prefetchnta 1088(%1)\n\t"

                "mulpd %%xmm0, %%xmm0\n\t"
                "mulpd %%xmm1, %%xmm1\n\t"
                "addpd %%xmm0, %%xmm1\n\t"
	  
                "mulpd %%xmm2, %%xmm2\n\t"
                "mulpd %%xmm3, %%xmm3\n\t"
                "addpd %%xmm2, %%xmm3\n\t"
                "addpd %%xmm1, %%xmm3\n\t"
	  
                "mulpd %%xmm4, %%xmm4\n\t"
                "mulpd %%xmm5, %%xmm5\n\t"
                "addpd %%xmm4, %%xmm5\n\t"
                "addpd %%xmm3, %%xmm5\n\t"
	  
                "mulpd %%xmm6, %%xmm6\n\t"
                "mulpd %%xmm7, %%xmm7\n\t"
                "addpd %%xmm5, %%xmm10\n\t"
                "addpd %%xmm6, %%xmm10\n\t"
                "addpd %%xmm7, %%xmm10\n\t"
	  
                "add $128, %1\n\t"
	  
                "dec %6\n\t"
              "jg SSE_ComplexSquareNorm_schleife4\n\t"
	      "mov 40(%0), %6\n\t"
              "add 64(%0), %1\n\t"
//              "add 64(%0), %2\n\t"
	      
              "dec %5\n\t"
            "jg SSE_ComplexSquareNorm_schleife3\n\t"
	    "mov 32(%0), %5\n\t"
            "add 56(%0), %1\n\t"
//            "add 56(%0), %2\n\t"
	      
            "dec %4\n\t"
          "jg SSE_ComplexSquareNorm_schleife2\n\t"
  	  "mov 24(%0), %4\n\t"
          "add 48(%0), %1\n\t"
//          "add 48(%0), %2\n\t"
	      
          "dec %3\n\t"
        "jg SSE_ComplexSquareNorm_schleife1\n\t"
        "movhlps %%xmm10, %%xmm0\n\t"
        "addpd %%xmm10, %%xmm0\n\t"
        "movupd %%xmm0, (%7)"
          : "=r" (Dummy1), "=r" (Dummy2), "=r" (Dummy3), "=r" (Dummy4), "=r" (Dummy5), "=r" (Dummy6), "=r" (Dummy7), "=r" (Dummy8), "=r" (Dummy9), "=r" (Dummy10)
          : "r" (para_N));

  //Always false
  if (para_N[2]!=(L0)) {
    printf("%d\n",*Dummy1);
    res->print();
  }

  r = res->x;
  delete[] res;
  delete[] para_N;
}


/**** SSE_ZCopy ***
*
**/
void SSE_ZCopy(int count, Complex* source, int j1, Complex* dest, int j2) {
  int I;
  for (I=0; I<count; I++) {
    dest[I].x = source[I].x;
    dest[I].y = source[I].y;
  }
}


