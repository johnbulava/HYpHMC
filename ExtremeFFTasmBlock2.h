          // Forth STEP + SAVE
          "addpd %%xmm0, %%xmm8\n\t"
          "mulpd 208(%6), %%xmm0\n\t"
          "movapd %%xmm8,  256(%1,%10)\n\t"  // <- Write
          "subpd %%xmm8, %%xmm0\n\t"
          "movapd %%xmm0,  256(%12,%10)\n\t"  // <- Write
		  
          "addpd %%xmm1, %%xmm9\n\t"
          "mulpd 208(%6), %%xmm1\n\t"
          "movapd %%xmm9,  (%1,%10)\n\t"  // <- Write
          "subpd %%xmm9, %%xmm1\n\t"
          "movapd %%xmm1,  (%12,%10)\n\t"  // <- Write
		  
          "addpd %%xmm2, %%xmm10\n\t"
          "mulpd 208(%6), %%xmm2\n\t"
          "movapd %%xmm10,  256(%1,%9)\n\t"  // <- Write
          "subpd %%xmm10, %%xmm2\n\t"
          "movapd %%xmm2,  256(%12,%9)\n\t"  // <- Write

          "addpd %%xmm3, %%xmm11\n\t"
          "mulpd 208(%6), %%xmm3\n\t"
          "movapd %%xmm11,  (%1,%9)\n\t"  // <- Write
          "subpd %%xmm11, %%xmm3\n\t"
          "movapd %%xmm3,  (%12,%9)\n\t"

          "addpd %%xmm4, %%xmm12\n\t"
          "mulpd 208(%6), %%xmm4\n\t"
          "movapd %%xmm12, 256(%1,%8)\n\t"  // <- Write
          "subpd %%xmm12, %%xmm4\n\t"
          "movapd %%xmm4,  256(%12,%8)\n\t"  // <- Write
		  
          "addpd %%xmm5, %%xmm13\n\t"
          "mulpd 208(%6), %%xmm5\n\t"
          "movapd %%xmm13, (%1,%8)\n\t"  // <- Write
          "subpd %%xmm13, %%xmm5\n\t"
          "movapd %%xmm5,  (%12,%8)\n\t"  // <- Write

          "addpd %%xmm6, %%xmm14\n\t"
          "mulpd 208(%6), %%xmm6\n\t"
          "movapd %%xmm14, 256(%1)\n\t"  // <- Write
          "subpd %%xmm14, %%xmm6\n\t"
          "movapd %%xmm6, 256(%12)\n\t"  // <- Write

          "addpd %%xmm7, %%xmm15\n\t"
          "mulpd 208(%6), %%xmm7\n\t"
          "movapd %%xmm15, (%1)\n\t"  // <- Write
          "subpd %%xmm15, %%xmm7\n\t"
          "movapd %%xmm7, (%12)\n\t"  // <- Write
