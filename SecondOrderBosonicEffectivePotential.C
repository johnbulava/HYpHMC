#include "SecondOrderBosonicEffectivePotential.h"

SecondOrderBosonicEffectivePotential::SecondOrderBosonicEffectivePotential(int l0, int l1, int l2, int l3, double v, double yt, double yb, double m0sqr, double la0, double la6, double la8, double la10, bool m0iD, bool lamInDet, int pertOrd) { 
  if (LogLevel>1) printf("Initializing SecondOrderBosonicEffectivePotential with Lat=(%dx%dx%dx%d), vev=%f, yt=%f, yb=%f, m0Sqr=%f, lambdas=(%f,%f,%f,%f), m0InDet=%d, lambdasInDet=%d, perturbativeOrder=%d\n", l0, l1, l2, l3, v, yt, yb, m0sqr, la0, la6, la8, la10, m0iD, lamInDet, pertOrd);
  L0 = l0;
  L1 = l1;
  L2 = l2;
  L3 = l3;
  perturbativeOrder = pertOrd;
  Vol = L0*L1*L2*L3;
  m0InDet = m0iD;
  lambdasInDet = lamInDet;
  lam0 = NaN;
  lam6 = NaN;
  lam8 = NaN;
  lam10 = NaN; 
  vev = NaN;
  m0Sqr = NaN;  
  YukT = NaN;
  YukB = NaN;

  pHatSqr = new double[Vol];
  logBosDetDn = new double[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  FermionBubble = new double[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  pSpaceGPropDn = new double*[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  pSpaceHPropDn = new double*[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  pSpaceFPropDn = new ComplexMatrix***[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  xSpaceGPropDn = new double*[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  xSpaceHPropDn = new double*[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  xSpaceFPropDn = new ComplexMatrix***[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];  
  for (int I=0; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) {
    xSpaceGPropDn[I] = new double[Vol];
    xSpaceHPropDn[I] = new double[Vol];
    xSpaceFPropDn[I] = new ComplexMatrix**[4];
    pSpaceGPropDn[I] = new double[Vol];
    pSpaceHPropDn[I] = new double[Vol];
    pSpaceFPropDn[I] = new ComplexMatrix**[4];  
    if (perturbativeOrder>=3) {
      for (int I2=0; I2<4; I2++) {
        xSpaceFPropDn[I][I2] = new ComplexMatrix*[Vol];
        pSpaceFPropDn[I][I2] = new ComplexMatrix*[Vol];
        for (int I3=0; I3<Vol; I3++) {
          xSpaceFPropDn[I][I2][I3] = new ComplexMatrix(8);
          pSpaceFPropDn[I][I2][I3] = new ComplexMatrix(8);
	}
      }
    }
  }
  
  CouplingTermCoeffcients = new long int***[20];
  CouplingTermCoeffcientCount = new long int*[20];
  for (int I=0; I<20; I++) {
    CouplingTermCoeffcients[I] = new long int**[20];
    CouplingTermCoeffcientCount[I] = new long int[20];
    for (int I2=0; I2<20; I2++) {
      CouplingTermCoeffcients[I][I2] = new long int*[1000];
      CouplingTermCoeffcientCount[I][I2] = 0;
      for (int I3=0; I3<1000; I3++) {
        CouplingTermCoeffcients[I][I2][I3] = new long int[8];
      }
    }
  }

  xSpacePropSumDBcount = 0;
  xSpacePropSumDBid = new int*[1000000];
  xSpacePropSumDBres = new double[1000000];
  for (int I=0; I<1000000; I++) {
    xSpacePropSumDBid[I] = new int[2*(SecondOrderBosonicEffectivePotential_MAXDERIVE+1)];
  }

  FFTdummyArray = createSuperAlignedComplex(Vol);
  FFTplan = NULL;

  int* n = new int[4];
  n[0] = L0;
  n[1] = L1;
  n[2] = L2;
  n[3] = L3;
  
  int rank = 4;
  int howmany = 1;
  int* inembed = new int[4];
  inembed[0] = L0;
  inembed[1] = L1;
  inembed[2] = L2;
  inembed[3] = L3;  
  int istride = 1;
  int idist = 1;
  
  int* onembed = new int[4];
  onembed[0] = L0;
  onembed[1] = L1;
  onembed[2] = L2;
  onembed[3] = L3;    
  int ostride = 1;
  int odist = 1;

  FFTplan  = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)FFTdummyArray, inembed,
                                istride, idist,
	   			(fftw_complex*)FFTdummyArray, onembed, ostride, odist,
	  		        FFTW_FORWARD, FFTW_MEASURE);

  delete[] n;
  delete[] inembed;
  delete[] onembed;
  
  calculatePHatSqr();
  setLambdas(la0, la6, la8, la10); 
  setVeV(v);
  setM0Sqr(m0sqr);
  setYukawa(yt, yb);

  for (int N1=0; N1<=10; N1+=2) {
    for (int N2=0; N2<=10; N2+=2) {
      calcCoefficientsOfCouplingTermCombination(N1, N2);
    }
  }  
  
  //Check automatically derivatives
  double U[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  double nU[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  for (int I=0; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) {
    U[I]  = calcEffectivePotDn(I);
    nU[I] = numericalDerivativeOfEffPot_dv(vev, I, I-1, 1E-6);
    
    if (abs((U[I]-nU[I]) / U[I])>3E-5) {    
      printf("ERROR in SecondOrderBosonicEffectivePotential: Derivative d^%d/dv^%d incorrect!!!\n", I,I);
      printf("Analytical: %1.10f, Numerical: %1.10f\n", U[I], nU[I]);      
      exit(0);
    }
  }
}


SecondOrderBosonicEffectivePotential::~SecondOrderBosonicEffectivePotential() { 
  delete[] pHatSqr;
  delete[] logBosDetDn;
  delete[] FermionBubble;
  for (int I=0; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) {
    delete[] xSpaceGPropDn[I];
    delete[] xSpaceHPropDn[I];
    delete[] pSpaceGPropDn[I];
    delete[] pSpaceHPropDn[I];
    
    if (perturbativeOrder>=3) {
      for (int I2=0; I2<4; I2++) {
        for (int I3=0; I3<Vol; I3++) {
          delete xSpaceFPropDn[I][I2][I3];
          delete pSpaceFPropDn[I][I2][I3];
        }
        delete[] xSpaceFPropDn[I][I2];
        delete[] pSpaceFPropDn[I][I2];
      }  
    }
    delete[] xSpaceFPropDn[I];
    delete[] pSpaceFPropDn[I];
  }
  delete[] xSpaceGPropDn;
  delete[] xSpaceHPropDn;
  delete[] xSpaceFPropDn;
  delete[] pSpaceGPropDn;
  delete[] pSpaceHPropDn;
  delete[] pSpaceFPropDn;
  
  for (int I=0; I<20; I++) {
    for (int I2=0; I2<20; I2++) {
      for (int I3=0; I3<1000; I3++) {
        delete[] CouplingTermCoeffcients[I][I2][I3];
      }
      delete[] CouplingTermCoeffcients[I][I2];
    }
    delete[] CouplingTermCoeffcients[I];
    delete[] CouplingTermCoeffcientCount[I];
  }
  delete[] CouplingTermCoeffcients;
  delete[] CouplingTermCoeffcientCount;

  xSpacePropSumDBcount = 0;
  for (int I=0; I<1000000; I++) {
    delete[] xSpacePropSumDBid[I];
  }
  delete[] xSpacePropSumDBid;
  delete[] xSpacePropSumDBres;


  destroySuperAlignedComplex(FFTdummyArray);
  fftw_destroy_plan(FFTplan);
  FFTplan = NULL;
}


void SecondOrderBosonicEffectivePotential::calculatePHatSqr() {
  LatticeMomentumBins* latBin = new LatticeMomentumBins(L0, L1, L2, L3);
  for (int I=0; I<Vol; I++) {
    pHatSqr[I] = latBin->getLatMomSqrFromIndex(I);
  }
  delete latBin;
}


void SecondOrderBosonicEffectivePotential::calculateFermionBubble() {
  for (int I=0; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) {
    FermionBubble[I] = 0;
  }
  Complex d0(0,0);
  Complex d1(0,0);
  Complex d2(0,0);
  for (int i0=0; i0<L0; i0++) {
    for (int i1=0; i1<L1; i1++) {
      for (int i2=0; i2<L2; i2++) {
        for (int i3=0; i3<L3; i3++) {
          int q  = i3 + i2*L3 + i1*L3*L2 + i0*L3*L2*L1;
	  int qn = ((L3-i3)%L3) + ((L2-i2)%L2)*L3 + ((L1-i1)%L1)*L3*L2 + ((L0-i0)%L0)*L3*L2*L1;
	    
          d0 = d0 -xSpaceHPropDn[0][q]*((*xSpaceFPropDn[0][0][q])*(*xSpaceFPropDn[0][0][qn])).tres();	  
	  for (int j=1; j<4; j++) {
            d0 = d0 -xSpaceGPropDn[0][q]*((*xSpaceFPropDn[0][j][q])*(*xSpaceFPropDn[0][j][qn])).tres();
	  }
	  
          d1 = d1 -xSpaceHPropDn[1][q]*((*xSpaceFPropDn[0][0][q])*(*xSpaceFPropDn[0][0][qn])).tres();
          d1 = d1 -xSpaceHPropDn[0][q]*((*xSpaceFPropDn[1][0][q])*(*xSpaceFPropDn[0][0][qn])).tres();
          d1 = d1 -xSpaceHPropDn[0][q]*((*xSpaceFPropDn[0][0][q])*(*xSpaceFPropDn[1][0][qn])).tres();
	  for (int j=1; j<4; j++) {
            d1 = d1 -xSpaceGPropDn[1][q]*((*xSpaceFPropDn[0][j][q])*(*xSpaceFPropDn[0][j][qn])).tres();
            d1 = d1 -xSpaceGPropDn[0][q]*((*xSpaceFPropDn[1][j][q])*(*xSpaceFPropDn[0][j][qn])).tres();
            d1 = d1 -xSpaceGPropDn[0][q]*((*xSpaceFPropDn[0][j][q])*(*xSpaceFPropDn[1][j][qn])).tres();
	  }
	  
          d2 = d2 -xSpaceHPropDn[2][q]*((*xSpaceFPropDn[0][0][q])*(*xSpaceFPropDn[0][0][qn])).tres();
          d2 = d2 -xSpaceHPropDn[1][q]*((*xSpaceFPropDn[1][0][q])*(*xSpaceFPropDn[0][0][qn])).tres();
          d2 = d2 -xSpaceHPropDn[1][q]*((*xSpaceFPropDn[0][0][q])*(*xSpaceFPropDn[1][0][qn])).tres();
          d2 = d2 -xSpaceHPropDn[1][q]*((*xSpaceFPropDn[1][0][q])*(*xSpaceFPropDn[0][0][qn])).tres();
          d2 = d2 -xSpaceHPropDn[0][q]*((*xSpaceFPropDn[2][0][q])*(*xSpaceFPropDn[0][0][qn])).tres();
          d2 = d2 -xSpaceHPropDn[0][q]*((*xSpaceFPropDn[1][0][q])*(*xSpaceFPropDn[1][0][qn])).tres();
          d2 = d2 -xSpaceHPropDn[1][q]*((*xSpaceFPropDn[0][0][q])*(*xSpaceFPropDn[1][0][qn])).tres();
          d2 = d2 -xSpaceHPropDn[0][q]*((*xSpaceFPropDn[1][0][q])*(*xSpaceFPropDn[1][0][qn])).tres();
          d2 = d2 -xSpaceHPropDn[0][q]*((*xSpaceFPropDn[0][0][q])*(*xSpaceFPropDn[2][0][qn])).tres();
	  for (int j=1; j<4; j++) {
            d2 = d2 -xSpaceGPropDn[2][q]*((*xSpaceFPropDn[0][j][q])*(*xSpaceFPropDn[0][j][qn])).tres();
            d2 = d2 -xSpaceGPropDn[1][q]*((*xSpaceFPropDn[1][j][q])*(*xSpaceFPropDn[0][j][qn])).tres();
            d2 = d2 -xSpaceGPropDn[1][q]*((*xSpaceFPropDn[0][j][q])*(*xSpaceFPropDn[1][j][qn])).tres();
            d2 = d2 -xSpaceGPropDn[1][q]*((*xSpaceFPropDn[1][j][q])*(*xSpaceFPropDn[0][j][qn])).tres();
            d2 = d2 -xSpaceGPropDn[0][q]*((*xSpaceFPropDn[2][j][q])*(*xSpaceFPropDn[0][j][qn])).tres();
            d2 = d2 -xSpaceGPropDn[0][q]*((*xSpaceFPropDn[1][j][q])*(*xSpaceFPropDn[1][j][qn])).tres();
            d2 = d2 -xSpaceGPropDn[1][q]*((*xSpaceFPropDn[0][j][q])*(*xSpaceFPropDn[1][j][qn])).tres();
            d2 = d2 -xSpaceGPropDn[0][q]*((*xSpaceFPropDn[1][j][q])*(*xSpaceFPropDn[1][j][qn])).tres();
            d2 = d2 -xSpaceGPropDn[0][q]*((*xSpaceFPropDn[0][j][q])*(*xSpaceFPropDn[2][j][qn])).tres();
	  }
        }
      }
    }
  }
  FermionBubble[0] = d0.x;
  FermionBubble[1] = d1.x;
  FermionBubble[2] = d2.x;  
}


//Only 4x4 Propagator, ie only one fermion popagator
void SecondOrderBosonicEffectivePotential::calculateFermionProp(double rho, double r) {
  if (LogLevel>2) printf("Calculating Fermion-Propagators for second order bosonic effective potential...\n");

  if ((YukT<=0) || (isNaN(YukT)) || (YukB<=0) || (isNaN(YukB)) ||(isNaN(vev))) {
    for (int I=0; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) {
      for (int I2=0; I2<4; I2++) {
        for (int I3=0; I3<Vol; I3++) {
          pSpaceFPropDn[I][I2][I3]->setZero();
          xSpaceFPropDn[I][I2][I3]->setZero();	  
	}
      }
    }  
  } else {
    NeubergerMatrix* diracOp = new NeubergerMatrix(rho, r, 1, 1, 1, 1, 2);
    double fac = 0.5 / rho;
  
    ComplexMatrix ProjPlusSmall = getProjectorMatrix(+1);
    ComplexMatrix ProjMinusSmall = getProjectorMatrix(-1);    
    ComplexMatrix ThetaSmall[4];
    ComplexMatrix BMat[4];
    ComplexMatrix diagY(2);
    diagY.setZero();
    diagY.matrix[0][0].x = YukT;
    diagY.matrix[1][1].x = YukB;
    for (int I=0; I<4; I++) {
      ThetaSmall[I] = getThetaMatrix(I);
      BMat[I] = ComplexMatrix(diagY*ThetaSmall[I], ProjMinusSmall);
      ThetaSmall[I].dagger();
      BMat[I] = BMat[I] + ComplexMatrix(ThetaSmall[I]*diagY, ProjPlusSmall);      
      ThetaSmall[I].dagger();
    }
    
    ComplexVector vK[4];
    vK[0].resize(4);
    vK[1].resize(4);
    vK[2].resize(4);
    vK[3].resize(4);
  
    int count = 0;
    vector4D k;
    for (int i0=0; i0<L0; i0++) {
      k[0] = 2*pi*i0 / L0;
      for (int i1=0; i1<L1; i1++) {
       k[1] = 2*pi*i1 / L1;
        for (int i2=0; i2<L2; i2++) {
          k[2] = 2*pi*i2 / L2;
          for (int i3=0; i3<L3; i3++) {
            k[3] = 2*pi*i3 / L3;

            Complex ew = diracOp->analyticalEigenvalue(k);
            diracOp->analyticalEigenvectors(k, vK);
	    
            ComplexMatrix mK1(vK[0]);  
            ComplexMatrix mK2(vK[1]);  
            ComplexMatrix mK3(vK[2]);  
            ComplexMatrix mK4(vK[3]);  
	    
  	    Complex z = ComplexUnity - (fac*ew);
	    Complex nenT = ew + YukT*vev*z;
	    Complex dT[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
	    dT[0] = z / (nenT);
	    dT[1] = (-YukT)*z*z / (nenT*nenT);
	    dT[2] = (2*YukT*YukT)*z*z*z / (nenT*nenT*nenT);
	    Complex nenB = ew + YukB*vev*z;
	    Complex dB[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];	    
	    dB[0] = z / (nenB);
	    dB[1] = (-YukB)*z*z / (nenB*nenB);
	    dB[2] = (2*YukB*YukB)*z*z*z / (nenB*nenB*nenB);
	    
	    for (int I=0; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) {
              ComplexMatrix dummyMatTop = dT[I]*mK1 + dT[I]*mK2 + adj(dT[I])*mK3 + adj(dT[I])*mK4;
              ComplexMatrix dummyMatBot = dB[I]*mK1 + dB[I]*mK2 + adj(dB[I])*mK3 + adj(dB[I])*mK4;
	      
              ComplexMatrix dummyMat(8);
              dummyMat.setZero();
              dummyMat.insertMatrix(dummyMatTop, 0, 0);
              dummyMat.insertMatrix(dummyMatBot, 4, 4);
	      
              for (int I2=0; I2<4; I2++) {
  	        (*(pSpaceFPropDn[I][I2][count])) = BMat[I2] * dummyMat;
              }
	    }

	    count++;	    
   	  }
        }
      }
    }
    delete diracOp;    
  
    //x-Space Fermion-Propagator
    double vFac = 1.0 / Vol;
    for (int I=0; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) {
      for (int I2=0; I2<4; I2++) {
        for (int I3=0; I3<8; I3++) {
          for (int I4=0; I4<8; I4++) {
            for (int q=0; q<Vol; q++) {
              FFTdummyArray[q] = pSpaceFPropDn[I][I2][q]->matrix[I3][I4];
	    }

            //Do NOT ignore zero-momentum
            fftw_execute(FFTplan);   
            for (int q=0; q<Vol; q++) {
	      xSpaceFPropDn[I][I2][q]->matrix[I3][I4] = vFac*FFTdummyArray[q];
            }
          }
        }
      }
    }  
  }
}


void SecondOrderBosonicEffectivePotential::calculateProp() {
  if (LogLevel>2) printf("Calculating Propagators for second order bosonic effective potential...\n");
  double detM0Sqr = m0Sqr;
  double detLam0 = lam0;
  double detLam6 = lam6;
  double detLam8 = lam8;
  double detLam10 = lam10;
  
  if (!m0InDet) detM0Sqr = 0;
  if (!lambdasInDet) {
    detLam0 = 0;
    detLam6 = 0;
    detLam8 = 0;
    detLam10 = 0;
  }

  xSpacePropSumDBcount = 0;
  for (int I=0; I<Vol; I++) {
    double zh   = detLam0  * (12*vev*vev)
                + detLam6  * (30*vev*vev*vev*vev)
    	        + detLam8  * (56*vev*vev*vev*vev*vev*vev)
	        + detLam10 * (90*vev*vev*vev*vev*vev*vev*vev*vev);
    double zg   = detLam0  * ( 4*vev*vev)
                + detLam6  * ( 6*vev*vev*vev*vev)
	        + detLam8  * ( 8*vev*vev*vev*vev*vev*vev)
	        + detLam10 * (10*vev*vev*vev*vev*vev*vev*vev*vev);
    double dzh  = detLam0  * (24*vev)
                + detLam6  * (120*vev*vev*vev)
                + detLam8  * (336*vev*vev*vev*vev*vev)
	        + detLam10 * (720*vev*vev*vev*vev*vev*vev*vev);
    double dzg  = detLam0  * ( 8*vev)
                + detLam6  * (24*vev*vev*vev)
		+ detLam8  * (48*vev*vev*vev*vev*vev)
		+ detLam10 * (80*vev*vev*vev*vev*vev*vev*vev);
    double ddzh = detLam0  * (24)
                + detLam6  * (360*vev*vev)
		+ detLam8  * (1680*vev*vev*vev*vev)
		+ detLam10 * (5040*vev*vev*vev*vev*vev*vev);
    double ddzg = detLam0  * (8)
                + detLam6  * (72*vev*vev)
	        + detLam8  * (240*vev*vev*vev*vev)
	        + detLam10 * (560*vev*vev*vev*vev*vev*vev);

    pSpaceHPropDn[0][I] = 1.0 / (pHatSqr[I] + detM0Sqr + zh);
    pSpaceGPropDn[0][I] = 1.0 / (pHatSqr[I] + detM0Sqr + zg);   
    
    pSpaceHPropDn[1][I] = -dzh / sqr(pHatSqr[I] + detM0Sqr + zh);
    pSpaceGPropDn[1][I] = -dzg / sqr(pHatSqr[I] + detM0Sqr + zg);   
    
    pSpaceHPropDn[2][I] = (-ddzh / sqr(pHatSqr[I] + detM0Sqr + zh)) + 2*dzh*dzh/((pHatSqr[I] + detM0Sqr + zh)*sqr(pHatSqr[I] + detM0Sqr + zh));
    pSpaceGPropDn[2][I] = (-ddzg / sqr(pHatSqr[I] + detM0Sqr + zg)) + 2*dzg*dzg/((pHatSqr[I] + detM0Sqr + zg)*sqr(pHatSqr[I] + detM0Sqr + zg));   
  }
  
  //x-Space Higgs-Propagator
  double vFac = 1.0 / Vol;
  for (int D=0; D<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; D++) {
    for (int I=0; I<Vol; I++) {
      FFTdummyArray[I].x = pSpaceHPropDn[D][I];
      FFTdummyArray[I].y = 0;
    }
    //Ignore zero-momentum
    FFTdummyArray[0].x = 0;
    FFTdummyArray[0].y = 0;    
    fftw_execute(FFTplan);   
    for (int I=0; I<Vol; I++) {
      xSpaceHPropDn[D][I] = vFac*FFTdummyArray[I].x;
    }
  }

  //x-Space Goldstone-Propagator
  for (int D=0; D<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; D++) {
    for (int I=0; I<Vol; I++) {
      FFTdummyArray[I].x = pSpaceGPropDn[D][I];
      FFTdummyArray[I].y = 0;
    }
    //Ignore zero-momentum
    FFTdummyArray[0].x = 0;
    FFTdummyArray[0].y = 0;    
    fftw_execute(FFTplan);   
    for (int I=0; I<Vol; I++) {
      xSpaceGPropDn[D][I] = vFac*FFTdummyArray[I].x;
    }
  }
}


void SecondOrderBosonicEffectivePotential::calculateLogBosDet() {
  for (int D=0; D<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; D++) logBosDetDn[D] = 0;
  //Ignore zero-momentum --> start at 1
  for (int I=1; I<Vol; I++) {
    logBosDetDn[0] += -0.5*log(pSpaceHPropDn[0][I]);
    logBosDetDn[0] += -1.5*log(pSpaceGPropDn[0][I]);

    logBosDetDn[1] += -0.5*pSpaceHPropDn[1][I] / pSpaceHPropDn[0][I];
    logBosDetDn[1] += -1.5*pSpaceGPropDn[1][I] / pSpaceGPropDn[0][I];    

    logBosDetDn[2] += -0.5*(pSpaceHPropDn[2][I]*pSpaceHPropDn[0][I] - pSpaceHPropDn[1][I]*pSpaceHPropDn[1][I]) / sqr(pSpaceHPropDn[0][I]);
    logBosDetDn[2] += -1.5*(pSpaceGPropDn[2][I]*pSpaceGPropDn[0][I] - pSpaceGPropDn[1][I]*pSpaceGPropDn[1][I]) / sqr(pSpaceGPropDn[0][I]);
  }
  for (int D=0; D<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; D++) logBosDetDn[D] /= Vol;  
}


void SecondOrderBosonicEffectivePotential::calcCoefficientsOfCouplingTermRecursion(int n, int N, int* facLine, int &resultCount, long int* resultStore) {
  //Fermion contribution
  if (N==12) {
    resultCount = 1;
    resultStore[0] = 1*100 + ((long int)100000000)*100*1;
    resultStore[1] = 2;     //due to top+bottom
    return;
  }

  if (n==N) {
    int typeCounter[5];
    for (int I=0; I<5; I++) typeCounter[I]=0;
    for (int I=0; I<N; I++) typeCounter[facLine[I]]++;

    int ID = typeCounter[0] + 100*typeCounter[1] + 10000*typeCounter[2] + 1000000*typeCounter[3] + 100000000*typeCounter[4];

    bool ignore = false;
    if ((N-typeCounter[0])==2) ignore = true; //Ignore due to exact integrability
    if ((N==2) && (!m0InDet)) ignore = false;
    if ((N>=4) && (!lambdasInDet)) ignore = false;    
    if ((N-typeCounter[0])==0) ignore = true; //Ignore due to exact integrability, ie v^N    
    
    if (!ignore) {
      bool entryFound = false;
      for (int I=0; I<resultCount; I++) if (resultStore[2*I]==ID) {
        entryFound = true;
        resultStore[2*I+1] += 1;
        break;
      }
      if (!entryFound) {
        resultStore[2*resultCount+0] = ID;
        resultStore[2*resultCount+1] = 1;
        resultCount++;
      }
    }
  
    return;
  }
  facLine[n+0] = 0;
  facLine[n+1] = 0;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);
  
  facLine[n+0] = 0;
  facLine[n+1] = 1;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);
  
  facLine[n+0] = 1;
  facLine[n+1] = 0;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);

  facLine[n+0] = 1;
  facLine[n+1] = 1;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);

  facLine[n+0] = 2;
  facLine[n+1] = 2;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);

  facLine[n+0] = 3;
  facLine[n+1] = 3;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);

  facLine[n+0] = 4;
  facLine[n+1] = 4;
  calcCoefficientsOfCouplingTermRecursion(n+2, N, facLine, resultCount, resultStore);
}


void SecondOrderBosonicEffectivePotential::calcCoefficientsOfCouplingTermCombination(int N1, int N2) {
  int* facLine = new int[N1+N2];
  long int* contractFacsH = new long int[2*(N1+N2+1)];
  long int* contractFacsG1 = new long int[2*(N1+N2+1)];
  long int* contractFacsG2 = new long int[2*(N1+N2+1)];
  long int* contractFacsG3 = new long int[2*(N1+N2+1)];
  long int* resultStore1 = new long int[200000];
  long int* resultStore2 = new long int[200000];
  int resultCount1 = 0;
  int resultCount2 = 0;
  CouplingTermCoeffcientCount[N1][N2] = 0;
  calcCoefficientsOfCouplingTermRecursion(0, N1, facLine, resultCount1, resultStore1);
  calcCoefficientsOfCouplingTermRecursion(0, N2, facLine, resultCount2, resultStore2);

  int doubleSum = 1;
  if (N1==0) {
    resultCount1 = 1;
    resultStore1[0] = 0;
    resultStore1[1] = 1;
    doubleSum = 0;
  }
  if (N2==0) {
    resultCount2 = 1;
    resultStore2[0] = 0;
    resultStore2[1] = 1;
    doubleSum = 0;
  }

  for (int I1=0; I1<resultCount1; I1++) {
    for (int I2=0; I2<resultCount2; I2++) {
      long int expCode = resultStore1[2*I1+0];
      int exp1V = expCode % 100;
      expCode /= 100;
      int exp1H = expCode % 100;
      expCode /= 100;
      int exp1G1 = expCode % 100;
      expCode /= 100;
      int exp1G2 = expCode % 100;
      expCode /= 100;
      int exp1G3 = expCode % 100;
      expCode /= 100;
      int exp1F = expCode % 100;

      expCode = resultStore2[2*I2+0];
      int exp2V = expCode % 100;
      expCode /= 100;
      int exp2H = expCode % 100;
      expCode /= 100;
      int exp2G1 = expCode % 100;
      expCode /= 100;
      int exp2G2 = expCode % 100;
      expCode /= 100;
      int exp2G3 = expCode % 100;
      expCode /= 100;
      int exp2F = expCode % 100;
          
      int fac = resultStore1[2*I1+1] * resultStore2[2*I2+1];
      int iHmax = NumberOfContractionsForRealScalarFieldAtTwo(exp1H, exp2H, contractFacsH);            
      int iG1max = NumberOfContractionsForRealScalarFieldAtTwo(exp1G1, exp2G1, contractFacsG1);
      int iG2max = NumberOfContractionsForRealScalarFieldAtTwo(exp1G2, exp2G2, contractFacsG2);
      int iG3max = NumberOfContractionsForRealScalarFieldAtTwo(exp1G3, exp2G3, contractFacsG3);
      
      for (int iH=0; iH<iHmax; iH++) if (contractFacsH[iH]>0) {
        for (int iG1=0; iG1<iG1max; iG1++) if (contractFacsG1[iG1]>0) {
          for (int iG2=0; iG2<iG2max; iG2++) if (contractFacsG2[iG2]>0) {
            for (int iG3=0; iG3<iG3max; iG3++) if (contractFacsG3[iG3]>0) {
              long int termFac = fac * contractFacsH[iH] * contractFacsG1[iG1] * contractFacsG2[iG2] * contractFacsG3[iG3];
	      long int powV    = exp1V+exp2V;
	      long int powGH0  = iHmax-iH-1;
	      long int powGHXY = iH;

	      long int powGG0  = (iG1max-iG1-1) + (iG2max-iG2-1) + (iG3max-iG3-1);
	      long int powGGXY = iG1+iG2+iG3;
	      
	      long int powF = exp1F+exp2F;	      
	    
	      //Make Entry
	      bool slotFound = false;
	      for (int j=0; j<CouplingTermCoeffcientCount[N1][N2]; j++) {
	        if ((CouplingTermCoeffcients[N1][N2][j][1]==doubleSum) && (CouplingTermCoeffcients[N1][N2][j][2]==powV)
		&& (CouplingTermCoeffcients[N1][N2][j][3]==powGH0) && (CouplingTermCoeffcients[N1][N2][j][4]==powGHXY)
		&& (CouplingTermCoeffcients[N1][N2][j][5]==powGG0) && (CouplingTermCoeffcients[N1][N2][j][6]==powGGXY)
		&& (CouplingTermCoeffcients[N1][N2][j][7]==powF)) {
		  slotFound = true;
		  CouplingTermCoeffcients[N1][N2][j][0] += termFac;
		}
	      }
	      if (!slotFound) {
	        CouplingTermCoeffcients[N1][N2][CouplingTermCoeffcientCount[N1][N2]][0] = termFac;
	        CouplingTermCoeffcients[N1][N2][CouplingTermCoeffcientCount[N1][N2]][1] = doubleSum;
	        CouplingTermCoeffcients[N1][N2][CouplingTermCoeffcientCount[N1][N2]][2] = powV;
	        CouplingTermCoeffcients[N1][N2][CouplingTermCoeffcientCount[N1][N2]][3] = powGH0;
	        CouplingTermCoeffcients[N1][N2][CouplingTermCoeffcientCount[N1][N2]][4] = powGHXY;
	        CouplingTermCoeffcients[N1][N2][CouplingTermCoeffcientCount[N1][N2]][5] = powGG0;
	        CouplingTermCoeffcients[N1][N2][CouplingTermCoeffcientCount[N1][N2]][6] = powGGXY;
	        CouplingTermCoeffcients[N1][N2][CouplingTermCoeffcientCount[N1][N2]][7] = powF;		
	        CouplingTermCoeffcientCount[N1][N2]++;				
	      }
	    }
	  }
	}
      } 
    }
  }
    
  delete[] facLine;
  delete[] resultStore1;
  delete[] resultStore2;
  delete[] contractFacsH;
  delete[] contractFacsG1;
  delete[] contractFacsG2;
  delete[] contractFacsG3;
}


void SecondOrderBosonicEffectivePotential::printCoefficientsOfCouplingTermCombination(int N1, int N2) {
  printf("\n\n");
  printf("Printing coupling coefficients for second order for N1: %d, N2: %d\n",N1,N2);
  for (int I=0; I<CouplingTermCoeffcientCount[N1][N2]; I++) {
    long int fac = CouplingTermCoeffcients[N1][N2][I][0];
    long int doubleSum = CouplingTermCoeffcients[N1][N2][I][1];
    long int powV = CouplingTermCoeffcients[N1][N2][I][2];
    long int powGH0 = CouplingTermCoeffcients[N1][N2][I][3];
    long int powGHXY = CouplingTermCoeffcients[N1][N2][I][4];
    long int powGG0 = CouplingTermCoeffcients[N1][N2][I][5];
    long int powGGXY = CouplingTermCoeffcients[N1][N2][I][6];
    long int powF = CouplingTermCoeffcients[N1][N2][I][7];
  
    printf("  (%ld): v^%ld GH0^%ld GHxy^%ld GG0^%ld GGxy^%ld  GF0^%ld *  %ld\n", doubleSum, powV, powGH0, powGHXY, powGG0, powGGXY, powF, fac);
  }
}


void SecondOrderBosonicEffectivePotential::evaluateCouplingTermDn(int N1, int N2, double* resDn, bool ignoreDisconnected) {
  long int* tableCount = new long int[SecondOrderBosonicEffectivePotential_MAXDERIVE+1];
  long int*** table = new long int**[SecondOrderBosonicEffectivePotential_MAXDERIVE+1];
  int colMax = 2+6*(1+SecondOrderBosonicEffectivePotential_MAXDERIVE);
  
  //Initialize derivative table
  for (int I=0; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) {
    tableCount[I] = 0;
    int lineMax = CouplingTermCoeffcientCount[N1][N2];
    for (int I2=0; I2<I; I2++) lineMax *= (6+I2);
    table[I] = new long int*[lineMax];
    for (int I2=0; I2<lineMax; I2++) {
      table[I][I2] = new long int[colMax];
      for (int I3=0; I3<colMax; I3++) {
        table[I][I2][I3] = 0;
      }
    }
  }
  
  //Load base polynom
  for (int I=0; I<CouplingTermCoeffcientCount[N1][N2]; I++) { 
    bool ignoreTerm = false;
    if (ignoreDisconnected) {
      if ((CouplingTermCoeffcients[N1][N2][I][4] == 0) && (CouplingTermCoeffcients[N1][N2][I][6] == 0)) ignoreTerm = true;
    }

    if (!ignoreTerm) {
      for (int j=0; j<8; j++) {
        table[0][tableCount[0]][j] = CouplingTermCoeffcients[N1][N2][I][j];
      }
      tableCount[0]++;
    }
  }
  
  //Calculate derivative table
  for (int I=1; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) {
    for (int I2=0; I2<tableCount[I-1]; I2++) {
      for (int I3=2; I3<colMax; I3++) {
        if (table[I-1][I2][I3] > 0) {
          for (int I4=0; I4<colMax; I4++) {
 	    table[I][tableCount[I]][I4] = table[I-1][I2][I4];
	  }
	  table[I][tableCount[I]][0] *= table[I-1][I2][I3];
	  table[I][tableCount[I]][I3] -= 1;
	  table[I][tableCount[I]][I3+6] += 1;
	  tableCount[I]++;
	}
      }
    }
  }
  
  //Evaluate derivative table
  double vderive[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  for (int I=0; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) vderive[I] = 0;
  vderive[0] = vev;
  vderive[1] = 1;
  for (int n=0; n<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; n++) {
    resDn[n] = 0;
    for (int I=0; I<tableCount[n]; I++) {
      double contrib = table[n][I][0];
      bool doubleSum = (bool) table[n][I][1];

      //Powers of w
      for (int I2=0; I2<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I2++) {
        for (int I3=0; I3<table[n][I][2+6*I2]; I3++) {
          contrib *= vderive[I2];
        }
      }

      //Powers of GH0
      for (int I2=0; I2<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I2++) {
        for (int I3=0; I3<table[n][I][3+6*I2]; I3++) {
          contrib *= xSpaceHPropDn[I2][0];
        }
      }

      //Powers of GG0
      for (int I2=0; I2<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I2++) {
        for (int I3=0; I3<table[n][I][5+6*I2]; I3++) {
          contrib *= xSpaceGPropDn[I2][0];
        }
      }
      
      //Powers of GF0
      for (int I2=0; I2<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I2++) {
        for (int I3=0; I3<table[n][I][7+6*I2]; I3++) {
          printf("Fermion-Bubble not implemented!!!\n");
	  exit(0);
        }
      }

      
      //Search DB for result of space-time sum
      int sumFoundInDBindex = -1;
      for (int j=0; j<xSpacePropSumDBcount; j++) {
        bool match = true;
        for (int I2=0; I2<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I2++) {
	  if ((xSpacePropSumDBid[j][2*I2+0]!=table[n][I][4+6*I2]) || (xSpacePropSumDBid[j][2*I2+1]!=table[n][I][6+6*I2])) {
	    match = false;
	    break;
	  }
	}
	if (match) {
	  sumFoundInDBindex = j;
	  break;
	}
      }
      if (sumFoundInDBindex<0) {
        //Calculate sum:
	double res = 0;
	for (int q=0; q<Vol; q++) {
	  double prod = 1;
          for (int I2=0; I2<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I2++) {
            for (int I3=0; I3<table[n][I][4+6*I2]; I3++) {
  	      prod *= xSpaceHPropDn[I2][q];
	    }
            for (int I3=0; I3<table[n][I][6+6*I2]; I3++) {
  	      prod *= xSpaceGPropDn[I2][q];
	    }
	  }
	  res += prod;
	}
	
        for (int I2=0; I2<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I2++) {
  	  xSpacePropSumDBid[xSpacePropSumDBcount][2*I2+0] = table[n][I][4+6*I2];
  	  xSpacePropSumDBid[xSpacePropSumDBcount][2*I2+1] = table[n][I][6+6*I2];
	}	  
        xSpacePropSumDBres[xSpacePropSumDBcount] = res;
	sumFoundInDBindex = xSpacePropSumDBcount;
	xSpacePropSumDBcount++;
      } 
      
      if (doubleSum) contrib *= xSpacePropSumDBres[sumFoundInDBindex];
      resDn[n] += contrib;
    } 
  }

  //Desinitialize derivative table
  for (int I=0; I<1+SecondOrderBosonicEffectivePotential_MAXDERIVE; I++) {
    int lineMax = CouplingTermCoeffcientCount[N1][N2];
    for (int I2=0; I2<I; I2++) lineMax *= (6+I2);
    for (int I2=0; I2<lineMax; I2++) {
      delete[] table[I][I2];
    }
    delete[] table[I];
  }
  delete[] table;
}


void SecondOrderBosonicEffectivePotential::setLambdas(double la0, double la6, double la8, double la10) {
  if ((isNaN(lam0)) || (lam0 != la0)) {
    if ((isNaN(lam6)) || (lam6 != la6)) {
      if ((isNaN(lam8)) || (lam8 != la8)) {
        if ((isNaN(lam10)) || (lam10 != la10)) {
          lam0 = la0;
          lam6 = la6;
          lam8 = la8;
          lam10 = la10; 
          calculateProp();
          calculateLogBosDet();
	}
      }
    }
  }
}


void SecondOrderBosonicEffectivePotential::setVeV(double v) {
  if ((isNaN(vev)) || (vev != v)) {
    vev = v;
    calculateProp();
    calculateLogBosDet();
    if (perturbativeOrder>=3) {
      calculateFermionProp(1,0.5);  
      calculateFermionBubble();
    }
  }
}


void SecondOrderBosonicEffectivePotential::setM0Sqr(double m0sqr) {
  if ((isNaN(m0Sqr)) || (m0Sqr != m0sqr)) {
    m0Sqr = m0sqr; 
    if (m0InDet) {
      calculateProp();
      calculateLogBosDet();
    }
  }
}
 
void SecondOrderBosonicEffectivePotential::setYukawa(double yt, double yb) {
  if ((isNaN(YukT)) || (YukT != yt) || (isNaN(YukB)) || (YukB != yb)) {
    YukT = yt;
    YukB = yb;
    if (perturbativeOrder>=3) {
      calculateFermionProp(1,0.5);  
      calculateFermionBubble();
    }
  }
}
 
 
double SecondOrderBosonicEffectivePotential::calcEffectivePotDn(int n) {
  if (n<0) return NaN;
  if (n>SecondOrderBosonicEffectivePotential_MAXDERIVE) return NaN;
  
  double resDn[1+SecondOrderBosonicEffectivePotential_MAXDERIVE];
  double lamM[11];
  for (int I=0; I<11; I++) lamM[I] = 0;
  lamM[2] = 0.5*m0Sqr;
  lamM[4] = lam0;
  lamM[6] = lam6;
  lamM[8] = lam8;
  lamM[10] = lam10;
  double res = 0;

  if (perturbativeOrder>=0) {
    res += logBosDetDn[n];
  }
  
  if (perturbativeOrder>=1) {
    for (int N=2; N<=10; N+=2) {
      evaluateCouplingTermDn(N, 0, resDn, false);
      res += lamM[N] * resDn[n];    
    }  
  }

  if (perturbativeOrder>=2) {
    for (int N1=2; N1<=10; N1+=2) {
      for (int N2=2; N2<=10; N2+=2) {
        evaluateCouplingTermDn(N1, N2, resDn, true);
        res += -0.5 * lamM[N1] * lamM[N2] * resDn[n];    //Factor 0.5 due to 1/n! in expansion of exp, minus due to exp(-S)
      }
    }  
  }
  
  if (perturbativeOrder>=3) {
    res += -0.5*FermionBubble[n];   
  }
   
  //Tree-Level terms
  if (n==0) {
    res += 0.5*m0Sqr*vev*vev;
    res += lam0*vev*vev*vev*vev;
    res += lam6*vev*vev*vev*vev*vev*vev;
    res += lam8*vev*vev*vev*vev*vev*vev*vev*vev;
    res += lam10*vev*vev*vev*vev*vev*vev*vev*vev*vev*vev;
  }
  if (n==1) {
    res +=  m0Sqr*vev;
    res +=  4*lam0*vev*vev*vev;
    res +=  6*lam6*vev*vev*vev*vev*vev;
    res +=  8*lam8*vev*vev*vev*vev*vev*vev*vev;
    res += 10*lam10*vev*vev*vev*vev*vev*vev*vev*vev*vev;
  }
  if (n==2) {
    res += m0Sqr;
    res += 12*lam0*vev*vev;
    res += 30*lam6*vev*vev*vev*vev;
    res += 56*lam8*vev*vev*vev*vev*vev*vev;
    res += 90*lam10*vev*vev*vev*vev*vev*vev*vev*vev;
  }

  return res;
}


double SecondOrderBosonicEffectivePotential::numericalDerivativeOfEffPot_dv(double v, int order, int baseOrder, double h) {
  if (baseOrder<0) baseOrder = 0;
  if (order<0) return NaN;
  if (order<baseOrder) return NaN;
  if (order>baseOrder) return (numericalDerivativeOfEffPot_dv(v+0.5*h, order-1, baseOrder, h) - numericalDerivativeOfEffPot_dv(v-0.5*h, order-1, baseOrder, h)) / h;
  double vevMerker = vev;
  setVeV(v);
  double res = calcEffectivePotDn(baseOrder);
  setVeV(vevMerker);  
  return res;
}


void SecondOrderBosonicEffectivePotential::plotxSpacePropagators(int deriveNr) {
  char* fileName = new char[1000];
  snprintf(fileName, 1000, "xSpaceHiggsPropagatorDerivative%d.dat", deriveNr);
  FILE* file = fopen(fileName,"w");
  
  for (int I=0; I<Vol; I++) {
    int dummy = I;

    int dx3 = dummy % L3;
    if (dx3>L3/2) dx3 -= L3/2;
    dummy /= L3;

    int dx2 = dummy % L2;
    if (dx2>L2/2) dx2 -= L2/2;
    dummy /= L2;

    int dx1 = dummy % L1;
    if (dx1>L1/2) dx1 -= L1/2;
    dummy /= L1;

    int dx0 = dummy % L0;
    if (dx0>L0/2) dx0 -= L0/2;
    dummy /= L0;
    
    double deltaX = sqrt(dx0*dx0 + dx1*dx1 + dx2*dx2 + dx3*dx3);
    fprintf(file, "%1.15f %1.15f\n", deltaX, xSpaceHPropDn[deriveNr][I]);
  }
  fclose(file);

  snprintf(fileName, 1000, "xSpaceGoldstonePropagatorDerivative%d.dat", deriveNr);
  file = fopen(fileName,"w");
  
  for (int I=0; I<Vol; I++) {
    int dummy = I;

    int dx3 = dummy % L3;
    if (dx3>L3/2) dx3 -= L3/2;
    dummy /= L3;

    int dx2 = dummy % L2;
    if (dx2>L2/2) dx2 -= L2/2;
    dummy /= L2;

    int dx1 = dummy % L1;
    if (dx1>L1/2) dx1 -= L1/2;
    dummy /= L1;

    int dx0 = dummy % L0;
    if (dx0>L0/2) dx0 -= L0/2;
    dummy /= L0;
    
    double deltaX = sqrt(dx0*dx0 + dx1*dx1 + dx2*dx2 + dx3*dx3);
    fprintf(file, "%1.15f %1.15f\n", deltaX, xSpaceGPropDn[deriveNr][I]);
  }
  fclose(file);
  
  delete[] fileName;
}


void SecondOrderBosonicEffectivePotential::plotpSpacePropagators(int deriveNr) {
  char* fileName = new char[1000];
  snprintf(fileName, 1000, "pSpaceHiggsPropagatorDerivative%d.dat", deriveNr);
  FILE* file = fopen(fileName,"w");
  
  for (int I=0; I<Vol; I++) {
    fprintf(file, "%1.15f %1.15f\n", pHatSqr[I], pSpaceHPropDn[deriveNr][I]);
  }
  fclose(file);

  snprintf(fileName, 1000, "pSpaceGoldstonePropagatorDerivative%d.dat", deriveNr);
  file = fopen(fileName,"w");
  
  for (int I=0; I<Vol; I++) {
    fprintf(file, "%1.15f %1.15f\n", pHatSqr[I], pSpaceGPropDn[deriveNr][I]);
  }
  fclose(file);
  
  delete[] fileName;
}

