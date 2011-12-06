AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase::AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader, char* oName, char* nick) : AnalyzerObservable(fOps, aIOcon, SDreader, oName, nick) { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 4;
  }
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase::~AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase() {
}


bool AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  int LargestL = fermiOps->get1DSizeLargest(); 
  int timeDirection = fermiOps->getTimeDirection();
  double TOL = 1E-8;

  for (int I=0; I<getAnalyzerResultsCount(); I++) {
    analyzerResults[I] = 0;
  }

  double* phiField = phiFieldConf->getPhiFieldCopy();
  vector4D normalizedAvgPhiFieldVector;
  
  normalizedAvgPhiFieldVector[0] = 1.0;
  normalizedAvgPhiFieldVector[1] = 0;
  normalizedAvgPhiFieldVector[2] = 0;
  normalizedAvgPhiFieldVector[3] = 0;    
  
  if (fixGauge && randomGauge) {
    printf("FixGauge and RandomgGauge activated simultaneously in Observable %s\n", getObsName());
    exit(0);
  } else if (fixGauge && !randomGauge) {
    phiFieldConf->alignHiggsFieldDirection(phiField);
  } else if (!fixGauge && randomGauge) {
    phiFieldConf->randomGaugeRotation(phiField);
  } else {
    normalizedAvgPhiFieldVector[0] = phiFieldConf->getAvgPhiFieldVectorComponent(0) / phiFieldConf->getAvgPhiFieldVectorLength();  
    normalizedAvgPhiFieldVector[1] = phiFieldConf->getAvgPhiFieldVectorComponent(1) / phiFieldConf->getAvgPhiFieldVectorLength();  
    normalizedAvgPhiFieldVector[2] = phiFieldConf->getAvgPhiFieldVectorComponent(2) / phiFieldConf->getAvgPhiFieldVectorLength();  
    normalizedAvgPhiFieldVector[3] = phiFieldConf->getAvgPhiFieldVectorComponent(3) / phiFieldConf->getAvgPhiFieldVectorLength();  
  }
  if (!multiplyWithPhiMatBSelection) {  
    normalizedAvgPhiFieldVector[0] = 1.0;
    normalizedAvgPhiFieldVector[1] = 0; 
    normalizedAvgPhiFieldVector[2] = 0;
    normalizedAvgPhiFieldVector[3] = 0;    
  }

  bool daggered = false;
  if (projectorSelection>0) daggered = true;
  ComplexMatrix* normalizedAvgMatB1 = createPhiMatrix(normalizedAvgPhiFieldVector, (!daggered));
  ComplexMatrix* normalizedAvgMatB2 = createPhiMatrix(normalizedAvgPhiFieldVector, daggered);

  fermiOps->setPreconditioner(true, 1.1, 0);
  fermiOps->setQPreconditioner(true, 0.25, 0.25);    

  auxVec0 = auxVectors[0];
  auxVec1 = auxVectors[1];
  auxVec2 = auxVectors[2];
  auxVec3 = auxVectors[3];

  Complex***** TwoFPsiPsiBarMatrix = new Complex****[LargestL];
  for (int t1=0; t1<LargestL; t1++) {
    TwoFPsiPsiBarMatrix[t1] = new Complex***[LargestL];
    for (int t2=0; t2<LargestL; t2++) {
      TwoFPsiPsiBarMatrix[t1][t2] = new Complex**[8];
      for (int fInd1=0; fInd1<8; fInd1++) {
        TwoFPsiPsiBarMatrix[t1][t2][fInd1] = new Complex*[8];
        for (int fInd2=0; fInd2<8; fInd2++) {
  	  TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2] = new Complex[4];
	  for (int I=0; I<4; I++) {
            TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][I].x = NaN;
            TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][I].y = NaN;	  
	  }
        }	  
      }
    }
  }
  
  bool success = calcTwoFPsiPsiBarMatrixForPhiField(phiField, TwoFPsiPsiBarMatrix, TOL, LargestL, timeDirection, projectorSelection);


  for (int t1=0; t1<LargestL; t1++) {
    for (int t2=0; t2<LargestL; t2++) {
      int deltaT = t1-t2;
      if (deltaT<0) deltaT = -deltaT;
      int deltaT2 = LargestL-deltaT; 
      
      Complex scalar(0,0); 
      for (int s1=0; s1<8; s1++) {
        for (int s2=0; s2<8; s2++) {
          for (int s3=0; s3<8; s3++) {
            for (int s4=0; s4<8; s4++) {
              scalar = scalar + (normalizedAvgMatB1->matrix[s1][s2] *  normalizedAvgMatB2->matrix[s3][s4]) * 
       	            (  TwoFPsiPsiBarMatrix[t1][t1][s2][s1][0] * TwoFPsiPsiBarMatrix[t2][t2][s4][s3][1]  
		     - TwoFPsiPsiBarMatrix[t1][t2][s2][s3][2] * TwoFPsiPsiBarMatrix[t2][t1][s4][s1][3] ); 
	    }
	  }
	}
      }

      analyzerResults[2*deltaT+0] += scalar.x;
      analyzerResults[2*deltaT+1] += scalar.y;
      analyzerResults[2*deltaT2+0] += scalar.x;
      analyzerResults[2*deltaT2+1] += scalar.y;
    }
  }
  
  for (int I=0; I<1+LargestL; I++) {
    double normFac = 2*LargestL;
    if ((I==0) || (I==LargestL)) normFac = LargestL;
    normFac *= L0*L1*L2*L3;
    normFac /= LargestL;    
    normFac *= L0*L1*L2*L3;
    normFac /= LargestL;    

    analyzerResults[2*I+0] /= normFac;
    analyzerResults[2*I+1] /= normFac;
  }
  
  for (int I=0; I<LargestL; I++) {
    Complex scalar(0,0);   
    double normFac = 1.0;
    normFac *= L0*L1*L2*L3;
    normFac /= LargestL;    
    
    for (int s1=0; s1<8; s1++) {
      for (int s2=0; s2<8; s2++) {
        scalar = scalar + normalizedAvgMatB1->matrix[s1][s2] * TwoFPsiPsiBarMatrix[I][I][s2][s1][0];  
      }
    }
    analyzerResults[2*(LargestL+1)+2*I+0] = scalar.x / normFac;  
    analyzerResults[2*(LargestL+1)+2*I+1] = scalar.y / normFac;      
  }
  
  for (int I=0; I<LargestL; I++) {
    Complex scalar(0,0);   
    double normFac = 1.0;
    normFac *= L0*L1*L2*L3;
    normFac /= LargestL;    
    
    for (int s1=0; s1<8; s1++) {
      for (int s2=0; s2<8; s2++) {
        scalar = scalar + normalizedAvgMatB2->matrix[s1][s2] * TwoFPsiPsiBarMatrix[I][I][s2][s1][1];
      }
    }
    analyzerResults[4*LargestL+2+2*I+0] = scalar.x / normFac;  
    analyzerResults[4*LargestL+2+2*I+1] = scalar.y / normFac;  
  }

  for (int t1=0; t1<LargestL; t1++) {
    for (int t2=0; t2<LargestL; t2++) {
      for (int fInd1=0; fInd1<8; fInd1++) {
        for (int fInd2=0; fInd2<8; fInd2++) {
          delete[] TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2];
	}
        delete[] TwoFPsiPsiBarMatrix[t1][t2][fInd1];	
      }
      delete[] TwoFPsiPsiBarMatrix[t1][t2];	
    }
    delete[] TwoFPsiPsiBarMatrix[t1];
  }
  delete[] TwoFPsiPsiBarMatrix;
  delete normalizedAvgMatB1;
  delete normalizedAvgMatB2;

  return success;
}


int AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase::getNeededAuxVectorCount() {
  return 4;
}


int AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase::getAnalyzerResultsCount() {
  int LargestL = fermiOps->get1DSizeLargest(); 
  return 2*(1 + LargestL) + 2*2*LargestL;
}


void AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase::sampleFermioncScanVector(Complex* v, int t, int FermionIndex, int timeDirection) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  fermiOps->zeroFermionVector(v);
  
  int ind[4];
  for (ind[0]=0; ind[0]<L0; ind[0]++) {
    for (ind[1]=0; ind[1]<L1; ind[1]++) {
      for (ind[2]=0; ind[2]<L2; ind[2]++) {
        for (ind[3]=0; ind[3]<L3; ind[3]++) {
	  if (ind[timeDirection] == t) {
	    int index = FermionIndex;
	    index += 8 * ind[3];
	    index += 8 * (L3+xtraSize3) * ind[2];
	    index += 8 * (L3+xtraSize3) * (L2+xtraSize2) * ind[1];
	    index += 8 * (L3+xtraSize3) * (L2+xtraSize2) * (L1+xtraSize1) * ind[0];
            v[index].x = 1.0;
	  }
	}
      }
    }
  }
}


/*bool AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase::calcTwoFPsiPsiBarMatrixForPhiField(double* phiField, Complex***** TwoFPsiPsiBarMatrix, double TOL, int LargestL, int timeDirection, int projSign) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  bool success = true;  
  int neededIter;
  bool b;
  bool projPlusSign = true;
  if (projSign==-1) projPlusSign  = false;
  Complex scalar(0,0);     
  double rho, r;
  fermiOps->getDiracParameters(rho, r);


  for (int t1=0; t1<LargestL; t1++) {
    for (int fInd1=0; fInd1<8; fInd1++) {
      if (LogLevel>2) printf("2-Fermion Psi-PsiBar-Matrix for t1 = %d, fInd1 = %d\n", t1, fInd1);          
      sampleFermioncScanVector(auxVec0, t1, fInd1, timeDirection);
    
      bool projSurvived1 = false;  
      if ((((fInd1/2)%2)==0) && (projSign==1)) projSurvived1 = true;
      if ((((fInd1/2)%2)==1) && (projSign==-1)) projSurvived1 = true;	  

      if (projSurvived1) {
        //Inversion	
        fermiOps->executeGamma0(auxVec0, auxVec1);
        b = fermiOps->executeFermionMatrixDaggerInverseMatrixMultiplication(auxVec1, auxVec2, phiField, TOL, false, -1, neededIter);     
        if (!b) success = false;
	
	//Entry 0
        for (int fInd2=0; fInd2<8; fInd2++) {
          sampleFermioncScanVector(auxVec3, t1, fInd2, timeDirection);
	  fermiOps->executeProjectorHatMultiplication(auxVec3, auxVec3, projPlusSign, true); //OHNE Faktor 0.5
          fermiOps->executeGamma0(auxVec3, auxVec3);
          SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, auxVec2, auxVec3, scalar);
	  TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][0] = 0.5 * scalar * (-0.5/rho);
	}
	
	//Entry 2
        for (int t2=0; t2<LargestL; t2++) {
          for (int fInd2=0; fInd2<8; fInd2++) {
            bool projSurvived2 = false;  
            if ((((fInd2/2)%2)==0) && (projSign==1)) projSurvived2 = true;
            if ((((fInd2/2)%2)==1) && (projSign==-1)) projSurvived2 = true;	  
	    if (projSurvived2) {
              sampleFermioncScanVector(auxVec3, t2, fInd2, timeDirection);
              SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, auxVec2, auxVec3, scalar);
   	      TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][2] = scalar * (-0.5/rho);
	    } else {
   	      TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][2].x = 0;	    
   	      TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][2].y = 0;	    
	    }
	  }
	}
      } else {
        for (int fInd2=0; fInd2<8; fInd2++) {
	  TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][0].x = 0;
	  TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][0].y = 0;	
	}      
        for (int t2=0; t2<LargestL; t2++) {
          for (int fInd2=0; fInd2<8; fInd2++) {
            TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][2].x = 0;
            TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][2].y = 0;
	  }
	}
      }
      
      //Inversion
      fermiOps->executeProjectorHatMultiplication(auxVec0, auxVec1, projPlusSign, true); //OHNE Faktor 0.5
      b = fermiOps->executeFermionMatrixDaggerInverseMatrixMultiplication(auxVec1, auxVec2, phiField, TOL, false, -1, neededIter);     
      if (!b) success = false;
      
      //Entry 1
      for (int fInd2=0; fInd2<8; fInd2++) {
        bool projSurvived2 = false;  
        if ((((fInd2/2)%2)==0) && (projSign==1)) projSurvived2 = true;
        if ((((fInd2/2)%2)==1) && (projSign==-1)) projSurvived2 = true;	  
	if (projSurvived2) {
          sampleFermioncScanVector(auxVec3, t1, fInd2, timeDirection);
          SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, auxVec2, auxVec3, scalar);
          TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][1] = 0.5 * scalar * (-0.5/rho);
	} else {
          TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][1].x = 0;
          TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][1].y = 0;
	}
      }
      
      //Entry 3
      for (int t2=0; t2<LargestL; t2++) {
        for (int fInd2=0; fInd2<8; fInd2++) {
          sampleFermioncScanVector(auxVec3, t2, fInd2, timeDirection);
          fermiOps->executeProjectorHatMultiplication(auxVec3, auxVec3, projPlusSign, true); //OHNE Faktor 0.5
          fermiOps->executeGamma0(auxVec3, auxVec3);
          SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, auxVec2, auxVec3, scalar);	  
          TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][3] = 0.5 * scalar * (-0.5/rho);
	}
      }
    }   
  }
  return success;
}*/


bool AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase::calcTwoFPsiPsiBarMatrixForPhiField(double* phiField, Complex***** TwoFPsiPsiBarMatrix, double TOL, int LargestL, int timeDirection, int projSign) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  bool success = true;  
  int neededIter;
  bool b;
  bool projPlusSign = true;
  if (projSign==-1) projPlusSign  = false;
  Complex scalar(0,0);     
  double rho, r;
  fermiOps->getDiracParameters(rho, r);


  for (int t1=0; t1<LargestL; t1++) {
    for (int fInd1=0; fInd1<8; fInd1++) {
      if (LogLevel>2) printf("2-Fermion Psi-PsiBar-Matrix for t1 = %d, fInd1 = %d\n", t1, fInd1);          
      sampleFermioncScanVector(auxVec0, t1, fInd1, timeDirection);
    
      //Inversion	
      fermiOps->executeProjectorHatMultiplication(auxVec0, auxVec1, !projPlusSign, true); //OHNE Faktor 0.5
      b = fermiOps->executeFermionMatrixDaggerInverseMatrixMultiplication(auxVec1, auxVec2, phiField, TOL, false, -1, neededIter);     
      if (!b) success = false;
	
     //Entry 0
      for (int fInd2=0; fInd2<8; fInd2++) {
        int fInd2ChiralSign = +1;  
        if (((fInd2/2)%2)==1) fInd2ChiralSign = -1;
        if (fInd2ChiralSign*projSign == -1) {
          sampleFermioncScanVector(auxVec3, t1, fInd2, timeDirection);
          SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, auxVec2, auxVec3, scalar);
	  TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][0] = 0.5 * scalar * (-0.5/rho);
	} else {
	  TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][0].x = 0;
	  TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][0].y = 0;	
	}
      }
	
      //Entry 2
      for (int t2=0; t2<LargestL; t2++) {
        for (int fInd2=0; fInd2<8; fInd2++) {
          int fInd2ChiralSign = +1;  
          if (((fInd2/2)%2)==1) fInd2ChiralSign = -1;
          if (fInd2ChiralSign*projSign == +1) {
            sampleFermioncScanVector(auxVec3, t2, fInd2, timeDirection);
            SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, auxVec2, auxVec3, scalar);
            TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][2] = scalar * (-0.5/rho);
	  } else {
   	    TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][2].x = 0;	    
   	    TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][2].y = 0;	    
	  }
	}
      }

      
      //Inversion
      fermiOps->executeProjectorHatMultiplication(auxVec0, auxVec1, projPlusSign, true); //OHNE Faktor 0.5
      b = fermiOps->executeFermionMatrixDaggerInverseMatrixMultiplication(auxVec1, auxVec2, phiField, TOL, false, -1, neededIter);     
      if (!b) success = false;
      
      //Entry 1
      for (int fInd2=0; fInd2<8; fInd2++) {
        int fInd2ChiralSign = +1;  
        if (((fInd2/2)%2)==1) fInd2ChiralSign = -1;
        if (fInd2ChiralSign*projSign == +1) {
          sampleFermioncScanVector(auxVec3, t1, fInd2, timeDirection);
          SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, auxVec2, auxVec3, scalar);
          TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][1] = 0.5 * scalar * (-0.5/rho);
	} else {
          TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][1].x = 0;
          TwoFPsiPsiBarMatrix[t1][t1][fInd1][fInd2][1].y = 0;
	}
      }
      
      //Entry 3
      for (int t2=0; t2<LargestL; t2++) {
        for (int fInd2=0; fInd2<8; fInd2++) {
          int fInd2ChiralSign = +1;  
          if (((fInd2/2)%2)==1) fInd2ChiralSign = -1;
          if (fInd2ChiralSign*projSign == -1) {
            sampleFermioncScanVector(auxVec3, t2, fInd2, timeDirection);
            SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, auxVec2, auxVec3, scalar);	  
            TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][3] = 0.5 * scalar * (-0.5/rho);
	  } else {
            TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][3].x = 0;
            TwoFPsiPsiBarMatrix[t1][t2][fInd1][fInd2][3].y = 0;
	  }
	}
      }
    }   
  }
  return success;
}
