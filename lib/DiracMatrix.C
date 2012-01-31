#include "DiracMatrix.h"

void DiracMatrix::iniD(int OneDimLSizeL0, int OneDimLSizeL1, int OneDimLSizeL2, int OneDimLSizeL3, int NCopies) {
  if (LogLevel>2) printf("DiracMatrix: Creating Dirac operator with Lattice size %dx%dx%dx%d, %d nested copies...", OneDimLSizeL0,OneDimLSizeL1,OneDimLSizeL2,OneDimLSizeL3,NCopies);
  NestedCopies = NCopies;
  OneDimLatticeSizeL0 = OneDimLSizeL0;
  OneDimLatticeSizeL1 = OneDimLSizeL1;
  OneDimLatticeSizeL2 = OneDimLSizeL2;
  OneDimLatticeSizeL3 = OneDimLSizeL3;

  constructOperator();
  
  if (LogLevel>2) printf("sucessfully.\n");
}



DiracMatrix::DiracMatrix() : ComplexMatrix(4) {
  iniD(1,1,1,1, 1);
}


DiracMatrix::DiracMatrix(const DiracMatrix& n) : ComplexMatrix(n) {
  NestedCopies = n.NestedCopies;
  OneDimLatticeSizeL0 = n.OneDimLatticeSizeL0;
  OneDimLatticeSizeL1 = n.OneDimLatticeSizeL1;
  OneDimLatticeSizeL2 = n.OneDimLatticeSizeL2;
  OneDimLatticeSizeL3 = n.OneDimLatticeSizeL3;
}


DiracMatrix::DiracMatrix(int N) : ComplexMatrix(N) {
}


/*void DiracMatrix::constructOperator() {
  int I0,I1,I2,I3,s0,s1,s2,s3,z0,z1,z2,z3;
  vector4D p;
  Complex nup, num;
  ComplexVector v[4];
  ComplexMatrix dummy(4);
  double expNorm;
  int posS, posZ, C;
  double scalar;

  v[0].resize(4);
  v[1].resize(4);
  v[2].resize(4);
  v[3].resize(4);
  
  expNorm = OneDimLatticeSizeL0*OneDimLatticeSizeL1*OneDimLatticeSizeL2*OneDimLatticeSizeL3;

  for (I0=0; I0<OneDimLatticeSizeL0; I0++) {
    p[0] = 2*I0*pi/OneDimLatticeSizeL0;

//p[0]+= pi/OneDimLatticeSizeL0;

    for (I1=0; I1<OneDimLatticeSizeL1; I1++) {
      p[1] = 2*I1*pi/OneDimLatticeSizeL1;
//p[1]+= pi/OneDimLatticeSizeL1;
      for (I2=0; I2<OneDimLatticeSizeL2; I2++) {
        p[2] = 2*I2*pi/OneDimLatticeSizeL2;
//p[2]+= pi/OneDimLatticeSizeL2;
        for (I3=0; I3<OneDimLatticeSizeL3; I3++) {
          p[3] = 2*I3*pi/OneDimLatticeSizeL3;
p[3]+= pi/OneDimLatticeSizeL3;

          nup = analyticalEigenvalue(p);
          num = adj(nup);
          analyticalEigenvectors(p,v);

          ComplexMatrix proj = nup*ComplexMatrix(v[0]);
          proj = (nup*ComplexMatrix(v[1])) + proj;
          proj = (num*ComplexMatrix(v[2])) + proj;
          proj = (num*ComplexMatrix(v[3])) + proj;



            posZ = 0;
            for (z0=0; z0<OneDimLatticeSizeL0; z0++) {
              for (z1=0; z1<OneDimLatticeSizeL1; z1++) {
                for (z2=0; z2<OneDimLatticeSizeL2; z2++) {
                  for (z3=0; z3<OneDimLatticeSizeL3; z3++) {
                    posS = 0;
                    for (s0=0; s0<OneDimLatticeSizeL0; s0++) {
                      for (s1=0; s1<OneDimLatticeSizeL1; s1++) {
                        for (s2=0; s2<OneDimLatticeSizeL2; s2++) {
                          for (s3=0; s3<OneDimLatticeSizeL3; s3++) {
                        
                            scalar  = p[0]*(z0-s0);
                            scalar += p[1]*(z1-s1);
                            scalar += p[2]*(z2-s2);
                            scalar += p[3]*(z3-s3);
          
                            dummy = (exp(scalar*ComplexI)/expNorm) * proj;

                            for (C = 0; C<NestedCopies; C++) {
                                addMatrix(dummy,posZ,posS);
                                posS += 4;
                                posZ += 4; 
                            }
                            posZ -= NestedCopies*4;
                          } 
                        }
                      }
                    }
                    posZ += NestedCopies*4;

                  }
                }
              }
            }
        }
      }
    }
  }  
}*/


void DiracMatrix::constructOperator() {
  int I0,I1,I2,I3,s0,s1,s2,s3,z0,z1,z2,z3,phaseInt;
  vector4D p;
  Complex nup, num;
  ComplexVector v[4];
  ComplexMatrix dummy(4);
  double expNorm;
  int posS, posZ, C;
  int scalarInt,scalarInt0,scalarInt1,scalarInt2,scalarInt3;

  v[0].resize(4);
  v[1].resize(4);
  v[2].resize(4);
  v[3].resize(4);
  
  expNorm = OneDimLatticeSizeL0*OneDimLatticeSizeL1*OneDimLatticeSizeL2*OneDimLatticeSizeL3;

  for (I0=0; I0<OneDimLatticeSizeL0; I0++) {
    p[0] = 2*I0*pi/OneDimLatticeSizeL0;
    for (I1=0; I1<OneDimLatticeSizeL1; I1++) {
      p[1] = 2*I1*pi/OneDimLatticeSizeL1;
      for (I2=0; I2<OneDimLatticeSizeL2; I2++) {
        p[2] = 2*I2*pi/OneDimLatticeSizeL2;
        for (I3=0; I3<OneDimLatticeSizeL3; I3++) {
          p[3] = 2*I3*pi/OneDimLatticeSizeL3;

	  nup = analyticalEigenvalue(p);
	  num = adj(nup);
	  analyticalEigenvectors(p,v);

	  ComplexMatrix proj = nup*ComplexMatrix(v[0]);
	  proj = (nup*ComplexMatrix(v[1])) + proj;
 	  proj = (num*ComplexMatrix(v[2])) + proj;
	  proj = (num*ComplexMatrix(v[3])) + proj;


          for (phaseInt=0; phaseInt<OneDimLatticeSizeL0*OneDimLatticeSizeL1*OneDimLatticeSizeL2*OneDimLatticeSizeL3; phaseInt++) {
            bool ready = false;

            posZ = 0;
            for (z0=0; z0<OneDimLatticeSizeL0; z0++) {
              for (z1=0; z1<OneDimLatticeSizeL1; z1++) {
                for (z2=0; z2<OneDimLatticeSizeL2; z2++) {
                  for (z3=0; z3<OneDimLatticeSizeL3; z3++) {
		    posS = 0;
                    for (s0=0; s0<OneDimLatticeSizeL0; s0++) {
  		      scalarInt0 = I0*(z0-s0);
                      for (s1=0; s1<OneDimLatticeSizeL1; s1++) {
  		        scalarInt1 = I1*(z1-s1);
                        for (s2=0; s2<OneDimLatticeSizeL2; s2++) {
  		          scalarInt2 = I2*(z2-s2);
                          for (s3=0; s3<OneDimLatticeSizeL3; s3++) {
    		            scalarInt3 = I3*(z3-s3);
			
		            scalarInt = scalarInt0*OneDimLatticeSizeL1*OneDimLatticeSizeL2*OneDimLatticeSizeL3
			              + scalarInt1*OneDimLatticeSizeL0*OneDimLatticeSizeL2*OneDimLatticeSizeL3
			              + scalarInt2*OneDimLatticeSizeL0*OneDimLatticeSizeL1*OneDimLatticeSizeL3
			              + scalarInt3*OneDimLatticeSizeL0*OneDimLatticeSizeL1*OneDimLatticeSizeL2;
			    
			    if ((scalarInt-phaseInt) % (OneDimLatticeSizeL0*OneDimLatticeSizeL1*OneDimLatticeSizeL2*OneDimLatticeSizeL3) == 0) {
			      if (!ready) {
                                dummy = (exp(phaseInt*2*pi/(OneDimLatticeSizeL0*OneDimLatticeSizeL1*OneDimLatticeSizeL2*OneDimLatticeSizeL3)*ComplexI)/expNorm) * proj;
                                ready = true;			      
			      }
 			      for (C = 0; C<NestedCopies; C++) {
  			        addMatrix(dummy,posZ,posS);
			        posS += 4;
			        posZ += 4; 
			      }
			      posZ -= NestedCopies*4;
			    } else {
			      posS += NestedCopies*4;
			    }
			  }
			}
		      }
		    }
  		    posZ += NestedCopies*4;
		  }
		}
	      }
	    }
	    
	  }
	  
	}
      }
    }
  }  
}


void DiracMatrix::analyticalEigenvectors(vector4D p, ComplexVector v[4]) {
  ComplexVector xi(2);
  Quat pThetaBAR(sin(p[0]),-sin(p[1]),-sin(p[2]),-sin(p[3]));
  ComplexMatrix pThetaBARMat(pThetaBAR);
  ComplexVector dummy(2);
  double pThetaBARNorm = pThetaBAR.getNorm();
  double root2 = sqrt(2);
  int I,I2;
  int count = 0;
    
  if (pThetaBARNorm>100*eps) {
    for (I2=+1; I2>-2; I2-=2) {   
      for (I=0; I<2; I++) {
        xi.vectorElements[0].x = (1-I);
        xi.vectorElements[1].x = I;
        dummy = pThetaBARMat * xi;
      
        v[count].vectorElements[0] = xi.vectorElements[0] / root2;
        v[count].vectorElements[1] = xi.vectorElements[1] / root2;
        v[count].vectorElements[2] = I2*dummy.vectorElements[0] / (root2*pThetaBARNorm);
        v[count].vectorElements[3] = I2*dummy.vectorElements[1] / (root2*pThetaBARNorm);
        count++;
      } 
    }
  } else {
    for (I2=+1; I2>-2; I2-=2) {   
      for (I=0; I<2; I++) {
        xi.vectorElements[0].x = (1-I);
        xi.vectorElements[1].x = I;
      
        v[count].vectorElements[0] = xi.vectorElements[0] / root2;
        v[count].vectorElements[1] = xi.vectorElements[1] / root2;
        v[count].vectorElements[2] = I2*xi.vectorElements[0] / root2;
        v[count].vectorElements[3] = I2*xi.vectorElements[1] / root2;
        count++;
      } 
    }
  }
}
