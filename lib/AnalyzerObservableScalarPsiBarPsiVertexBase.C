#include "AnalyzerObservableScalarPsiBarPsiVertexBase.h"

AnalyzerObservableScalarPsiBarPsiVertexBase::AnalyzerObservableScalarPsiBarPsiVertexBase(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader, char* oName, char* nick) : AnalyzerObservable(fOps, aIOcon, SDreader, oName, nick) { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 1;
  }
  tresSumStartIndex = 0;
  tresSumEndIndex = 7;
  PsiPsiBarMatrixInd1Start = 0;
  PsiPsiBarMatrixInd2Start = 0;
  PsiPsiBarMatrixInd1End = 7;
  PsiPsiBarMatrixInd2End = 7;
  
  ini(getAnalyzerResultsCount());
}


AnalyzerObservableScalarPsiBarPsiVertexBase::~AnalyzerObservableScalarPsiBarPsiVertexBase() {
}


bool AnalyzerObservableScalarPsiBarPsiVertexBase::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  int L0 = fermiOps->get1DSizeL1();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double volume = L0*L1*L2*L3;
  int LargestL = fermiOps->get1DSizeLargest(); 
  int timeDirection = fermiOps->getTimeDirection();
  double TOL = 1E-8;

  for (int I=0; I<getAnalyzerResultsCount(); I++) {
    analyzerResults[I] = 0;
  }

  double* phiField = phiFieldConf->getPhiFieldCopy();

  if (fixGauge) {
    phiFieldConf->alignHiggsFieldDirection(phiField);
  } 
  Complex* phiMomentumBuffer = phiFieldConf->performFourierTransform(phiField, true, 1);  
  double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);

  LeftVector = auxVectors[0];
  RightVector = auxVectors[1];
  SolutionVector = auxVectors[2];

  Complex**** PsiPsiBarMatrix = new Complex***[LargestL];
  for (int t1=0; t1<LargestL; t1++) {
    PsiPsiBarMatrix[t1] = new Complex**[LargestL];
    for (int t2=0; t2<LargestL; t2++) {
      PsiPsiBarMatrix[t1][t2] = new Complex*[8];
      for (int fInd1=0; fInd1<8; fInd1++) {
        PsiPsiBarMatrix[t1][t2][fInd1] = new Complex[8];
        for (int fInd2=0; fInd2<8; fInd2++) {
          PsiPsiBarMatrix[t1][t2][fInd1][fInd2].x = NaN;
          PsiPsiBarMatrix[t1][t2][fInd1][fInd2].y = NaN;	  
        }	  
      }
    }
  }

  fermiOps->setPreconditioner(true, 1.1, 0);
  fermiOps->setQPreconditioner(true, 0.25, 0.25);    
  
  bool success = calcPsiPsiBarMatrixForPhiField(phiField, PsiPsiBarMatrix, TOL, projectorSelection, PsiPsiBarMatrixInd1Start, PsiPsiBarMatrixInd1End, PsiPsiBarMatrixInd2Start, PsiPsiBarMatrixInd2End);

  int ind[4];
  ind[0] = 0;
  ind[1] = 0;
  ind[2] = 0;
  ind[3] = 0;  
  for (int t1=0; t1<LargestL; t1++) {
    for (int t2=0; t2<LargestL; t2++) {
      int deltaT = t2-t1;
      if (deltaT<0) deltaT += LargestL;
      ind[timeDirection] = deltaT;
      int pInd = ind[3];
      pInd += ind[2]*L3;
      pInd += ind[1]*L3*L2;
      pInd += ind[0]*L3*L2*L1;
      
      Complex scalar(0,0); 
      for (int I=tresSumStartIndex; I<=tresSumEndIndex; I++) {
        scalar = scalar + PsiPsiBarMatrix[t1][t2][I][I];
      }

      Complex scalarLeg(0,0); 
      for (int I1=PsiPsiBarMatrixInd1Start; I1<=PsiPsiBarMatrixInd1End; I1++) {
        for (int I2=PsiPsiBarMatrixInd2Start; I2<=PsiPsiBarMatrixInd2End; I2++) {
          scalarLeg = scalarLeg + PsiPsiBarMatrix[t1][t1][I1][I2] * PsiPsiBarMatrix[t2][t2][I2][I1];
	}
      }

      Complex hValP = rescale * (1.0/sqrt(volume)) * phiMomentumBuffer[pInd];

      analyzerResults[6*(t1*LargestL+t2)+0] = scalar.x;
      analyzerResults[6*(t1*LargestL+t2)+1] = scalar.y;
      analyzerResults[6*(t1*LargestL+t2)+2] = scalarLeg.x;
      analyzerResults[6*(t1*LargestL+t2)+3] = scalarLeg.y;
      analyzerResults[6*(t1*LargestL+t2)+4] = hValP.x;      
      analyzerResults[6*(t1*LargestL+t2)+5] = hValP.y;      
    }
  }

  for (int t1=0; t1<LargestL; t1++) {
    for (int t2=0; t2<LargestL; t2++) {
      for (int fInd1=0; fInd1<8; fInd1++) {
        delete[] PsiPsiBarMatrix[t1][t2][fInd1];
      }
      delete[] PsiPsiBarMatrix[t1][t2];	
    }
    delete[] PsiPsiBarMatrix[t1];
  }
  delete[] PsiPsiBarMatrix;

  return success;
}


int AnalyzerObservableScalarPsiBarPsiVertexBase::getNeededAuxVectorCount() {
  return 3;
}


int AnalyzerObservableScalarPsiBarPsiVertexBase::getAnalyzerResultsCount() {
  int LargestL = fermiOps->get1DSizeLargest(); 
  return 6*LargestL*LargestL;
}


void AnalyzerObservableScalarPsiBarPsiVertexBase::sampleFermioncScanVector(Complex* v, int pft, int FermionIndex) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  int LargestL = fermiOps->get1DSizeLargest(); 
  double volume = L0*L1*L2*L3;
  double vFac = 1.0 / sqrt(volume);
  int timeDirection = fermiOps->getTimeDirection();
  
  fermiOps->zeroFermionVector(v);
  
  int ind[4];
  for (ind[0]=0; ind[0]<L0; ind[0]++) {
    for (ind[1]=0; ind[1]<L1; ind[1]++) {
      for (ind[2]=0; ind[2]<L2; ind[2]++) {
        for (ind[3]=0; ind[3]<L3; ind[3]++) {
          int index = FermionIndex;
	  index += 8 * ind[3];
	  index += 8 * (L3+xtraSize3) * ind[2];
	  index += 8 * (L3+xtraSize3) * (L2+xtraSize2) * ind[1];
	  index += 8 * (L3+xtraSize3) * (L2+xtraSize2) * (L1+xtraSize1) * ind[0];
	  double phase = (2*pi*pft*ind[timeDirection])/LargestL;
	  
          v[index] = vFac * exp(phase*ComplexI);
	}
      }
    }
  }
}


bool AnalyzerObservableScalarPsiBarPsiVertexBase::calcPsiPsiBarMatrixForPhiField(double* phiField, Complex**** PsiPsiBarMatrix, double TOL, int projMode, int fInd1Start, int fInd1End, int fInd2Start, int fInd2End) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  int LargestL = fermiOps->get1DSizeLargest(); 
  bool Quse;
  double mu, beta;
  fermiOps->getQPreconditionerParameter(Quse, mu, beta);
  bool success = true;  
  double rho, r;
  fermiOps->getDiracParameters(rho, r);
  double fac = -2.0*rho;
  if (projMode != 0) fac *= 0.5;

  for (int t1=0; t1<LargestL; t1++) {
    for (int fInd1=fInd1Start; fInd1<=fInd1End; fInd1++) {
      if (LogLevel>2) printf("Psi-PsiBar-Matrix for pft = %d, fInd1 = %d\n", t1, fInd1);          
      sampleFermioncScanVector(LeftVector, t1, fInd1);
      
      //Berechne Projector \hat P_\pm angewendet auf LeftVector  ==> Faktor 0.5 fehlt
      if (projMode != 0) {
        fermiOps->executeGamma5(LeftVector, SolutionVector);
	fermiOps->executeDiracDaggerMatrixMultiplication(SolutionVector, RightVector, false);
	Complex alpha1(-1.0/rho, 0);
	SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, alpha1, RightVector, SolutionVector);
	Complex alpha2(projMode, 0);
	SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, alpha2, SolutionVector, LeftVector);	
      }
      
      int neededIter = 0;
      bool b = fermiOps->executeFermionMatrixDaggerInverseMatrixMultiplication(LeftVector, SolutionVector, phiField, TOL, false, -1, neededIter);
      
      if (!b) success = false;

      for (int t2=0; t2<LargestL; t2++) {
        for (int fInd2=fInd2Start; fInd2<=fInd2End; fInd2++) {
          Complex scalar(0,0);     
	  bool calc = false;
	  if (projMode==0) calc = true;
	  if ((((fInd2/2)%2)==0) && (projMode==1)) calc = true;
	  if ((((fInd2/2)%2)==1) && (projMode==-1)) calc = true;	  
	  if (calc) {
            sampleFermioncScanVector(RightVector, t2, fInd2);	  
            SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, SolutionVector, RightVector, scalar);
	  }
          PsiPsiBarMatrix[t1][t2][fInd1][fInd2] = fac*scalar;  // Hier ist Faktor 0.5 von erstem Projektor korrigiert!
        }   
      }
    }   
  }
  return success;
}
