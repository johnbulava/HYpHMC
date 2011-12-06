                  "addpd %%xmm0, %%xmm8\n\t"
                  "mulpd 208(%6), %%xmm0\n\t"
                  "movapd %%xmm8,  22528(%6,%9)\n\t"  // <- Write
                  "subpd %%xmm8, %%xmm0\n\t"
                  "movapd %%xmm0,  26624(%6,%9)\n\t"  // <- Write
		  
                  "addpd %%xmm1, %%xmm9\n\t"
                  "mulpd 208(%6), %%xmm1\n\t"
                  "movapd %%xmm9, 22016(%6,%9) \n\t"  // <- Write
                  "subpd %%xmm9, %%xmm1\n\t"
                  "movapd %%xmm1,  26112(%6,%9)\n\t"  // <- Write
		  
                  "addpd %%xmm2, %%xmm10\n\t"
                  "mulpd 208(%6), %%xmm2\n\t"
                  "movapd %%xmm10,  21504(%6,%9)\n\t"  // <- Write
                  "subpd %%xmm10, %%xmm2\n\t"
                  "movapd %%xmm2,  25600(%6,%9)\n\t"  // <- Write

                  "addpd %%xmm3, %%xmm11\n\t"
                  "mulpd 208(%6), %%xmm3\n\t"
                  "movapd %%xmm11, 20992(%6,%9) \n\t"  // <- Write
                  "subpd %%xmm11, %%xmm3\n\t"
                  "movapd %%xmm3,  25088(%6,%9)\n\t" // <- Writ

                  "addpd %%xmm4, %%xmm12\n\t"
                  "mulpd 208(%6), %%xmm4\n\t"
                  "movapd %%xmm12, 20480(%6,%9)\n\t"  // <- Write
                  "subpd %%xmm12, %%xmm4\n\t"
                  "movapd %%xmm4,  24576(%6,%9)\n\t"  // <- Write
		  
                  "addpd %%xmm5, %%xmm13\n\t"
                  "mulpd 208(%6), %%xmm5\n\t"
                  "movapd %%xmm13, 19968(%6,%9)\n\t"  // <- Write
                  "subpd %%xmm13, %%xmm5\n\t"
                  "movapd %%xmm5,  24064(%6,%9)\n\t"  // <- Write

                  "addpd %%xmm6, %%xmm14\n\t"
                  "mulpd 208(%6), %%xmm6\n\t"
                  "movapd %%xmm14, 19456(%6,%9)\n\t"  // <- Write
                  "subpd %%xmm14, %%xmm6\n\t"
                  "movapd %%xmm6, 23552(%6,%9)\n\t"  // <- Write

                  "addpd %%xmm7, %%xmm15\n\t"
                  "mulpd 208(%6), %%xmm7\n\t"
                  "movapd %%xmm15, 18944(%6,%9)\n\t"  // <- Write
                  "subpd %%xmm15, %%xmm7\n\t"
                  "movapd %%xmm7, 23040(%6,%9)\n\t"  // <- Write
