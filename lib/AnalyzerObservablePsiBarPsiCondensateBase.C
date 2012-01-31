#include "AnalyzerObservablePsiBarPsiCondensateBase.h"

AnalyzerObservablePsiBarPsiCondensateBase::AnalyzerObservablePsiBarPsiCondensateBase(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader, char* oName, char* nick) : AnalyzerObservable(fOps, aIOcon, SDreader, oName, nick) { 
  if (fermiOps->get1DSizeLargest() >= 32) {
    analyzeEveryXXXconf = 1;
  }
  fixGauge = false;
  randomGauge = false;
  projectorSelection = 1;
  tresSumStartIndex = 0;
  tresSumEndIndex = 7;
  PsiPsiBarMatrixInd1Start = 0;
  PsiPsiBarMatrixInd2Start = 0;
  PsiPsiBarMatrixInd1End = 7;
  PsiPsiBarMatrixInd2End = 7;
  numberOfMeasurements = 5;
  stochasticalSource = true;
  LeftVector = NULL;
  RightVector = NULL;
  SolutionVector = NULL;
  SourceVector = NULL;

  ini(getAnalyzerResultsCount());
}


AnalyzerObservablePsiBarPsiCondensateBase::~AnalyzerObservablePsiBarPsiCondensateBase() {
}


bool AnalyzerObservablePsiBarPsiCondensateBase::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double volume = L0*L1*L2*L3;
  double TOL = 1E-8;

  for (int I=0; I<getAnalyzerResultsCount(); I++) {
    analyzerResults[I] = 0;
  }

  double* phiField = phiFieldConf->getPhiFieldCopy();

  if (fixGauge && randomGauge) {
    printf("FixGauge and RandomGauge activated simultaneously in Observable %s!!!\n",getObsName());
    exit(0);
  } else if (fixGauge && !randomGauge) {
    phiFieldConf->alignHiggsFieldDirection(phiField);
  } else if (!fixGauge && randomGauge) {
    phiFieldConf->randomGaugeRotation(phiField);
  } 

  LeftVector = auxVectors[0];
  RightVector = auxVectors[1];
  SolutionVector = auxVectors[2];
  SourceVector = auxVectors[3];

  Complex** PsiPsiBarMatrix = new Complex*[8];
  for (int fInd1=0; fInd1<8; fInd1++) {
    PsiPsiBarMatrix[fInd1] = new Complex[8];
    for (int fInd2=0; fInd2<8; fInd2++) {
      PsiPsiBarMatrix[fInd1][fInd2].x = NaN;
      PsiPsiBarMatrix[fInd1][fInd2].y = NaN;	  
    }
  }

  fermiOps->setPreconditioner(true, 1.1, 0);
  fermiOps->setQPreconditioner(true, 0.25, 0.25);    
  
  bool success = true;

  for (int iter=0; iter<numberOfMeasurements; iter++) {
    sampleSource();

    bool b = calcPsiPsiBarMatrixForPhiField(phiField, PsiPsiBarMatrix, TOL, projectorSelection, PsiPsiBarMatrixInd1Start, PsiPsiBarMatrixInd1End, PsiPsiBarMatrixInd2Start, PsiPsiBarMatrixInd2End);
    if (!b) success = false;

    Complex scalar(0,0); 
    for (int I=tresSumStartIndex; I<=tresSumEndIndex; I++) {
      scalar = scalar + PsiPsiBarMatrix[I][I];
    }

    analyzerResults[2*iter+0] = scalar.x;
    analyzerResults[2*iter+1] = scalar.y;
  }

  if (stochasticalSource) {
    for (int I=0; I<getAnalyzerResultsCount(); I++) {
      analyzerResults[I] /= 2.0 * volume;   //Factor 2 because both, real and imaginary part of scan-vector, are Gauss distributed with variance 1 each.
    }
  }

  for (int fInd1=0; fInd1<8; fInd1++) {
    delete[] PsiPsiBarMatrix[fInd1];
  }
  delete[] PsiPsiBarMatrix;	

  return success;
}


int AnalyzerObservablePsiBarPsiCondensateBase::getNeededAuxVectorCount() {
  return 4;
}


int AnalyzerObservablePsiBarPsiCondensateBase::getAnalyzerResultsCount() {
  return 2*numberOfMeasurements;
}


void AnalyzerObservablePsiBarPsiCondensateBase::sampleSource() {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  fermiOps->zeroFermionVector(SourceVector);

  if (stochasticalSource) {
    if (LogLevel>2) printf("Generating stochastical source.\n");
    fermiOps->fillGaussRandomVector(SourceVector, -1);
  } else {
    if (LogLevel>2) printf("Generating point source.\n");
    int i0 = (int)(AdvancedZufall(AdvancedSeed) * L0);
    int i1 = (int)(AdvancedZufall(AdvancedSeed) * L1);
    int i2 = (int)(AdvancedZufall(AdvancedSeed) * L2);
    int i3 = (int)(AdvancedZufall(AdvancedSeed) * L3);
    if (i0>=L0) i0 = L0-1;
    if (i1>=L1) i1 = L1-1;
    if (i2>=L2) i2 = L2-1;
    if (i3>=L3) i3 = L3-1;

    int index = 0;
    index += 8 * i3;
    index += 8 * (L3) * i2;
    index += 8 * (L3) * (L2) * i1;
    index += 8 * (L3) * (L2) * (L1) * i0;

    SourceVector[index].x = 1.0;
    SourceVector[index].y = 0.0;
  }

  fermiOps->transformToXtraSizeArray(SourceVector,SourceVector);
}


void AnalyzerObservablePsiBarPsiCondensateBase::sampleFermioncScanVector(Complex* v, int FermionIndex) {
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
          int index = 0;
	  index += 8 * ind[3];
	  index += 8 * (L3+xtraSize3) * ind[2];
	  index += 8 * (L3+xtraSize3) * (L2+xtraSize2) * ind[1];
	  index += 8 * (L3+xtraSize3) * (L2+xtraSize2) * (L1+xtraSize1) * ind[0];

          v[index+FermionIndex].x = SourceVector[index].x;
          v[index+FermionIndex].y = SourceVector[index].y;
	}
      }
    }
  }
}


/****
* This routine computes <LeftVector | \hat P_\pm (\D+B(1-(1/\rho)\D))^-1 P_\pm | RightVector>
* where LeftVector and RightVector are sampled according to the function sampleFermioncScanVector.
* All factors are already included.
*/
bool AnalyzerObservablePsiBarPsiCondensateBase::calcPsiPsiBarMatrixForPhiField(double* phiField, Complex** PsiPsiBarMatrix, double TOL, int projMode, int fInd1Start, int fInd1End, int fInd2Start, int fInd2End) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  bool Quse;
  double mu, beta;
  fermiOps->getQPreconditionerParameter(Quse, mu, beta);
  bool success = true;  
  double rho, r;
  fermiOps->getDiracParameters(rho, r);
  double fac = -1.0 / (2.0*rho);
  if (projMode != 0) fac *= 0.5;  //Missing factor 0.5 in projector \hat P_\pm compensated with this factor 0.5. 

  for (int fInd1=fInd1Start; fInd1<=fInd1End; fInd1++) {
    if (LogLevel>2) printf("Psi-PsiBar-Matrix for fInd1 = %d\n", fInd1);          
    sampleFermioncScanVector(LeftVector, fInd1);
      
    //Berechne daggered Projector \hat P_\pm angewendet auf LeftVector  ==> Faktor 0.5 fehlt
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

    for (int fInd2=fInd2Start; fInd2<=fInd2End; fInd2++) {
      Complex scalar(0,0);     
      bool calc = false;
      if (projMode==0) calc = true;
      if ((((fInd2/2)%2)==0) && (projMode==1)) calc = true;
      if ((((fInd2/2)%2)==1) && (projMode==-1)) calc = true;	  
      if (calc) {
        sampleFermioncScanVector(RightVector, fInd2);	  
        SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, SolutionVector, RightVector, scalar);
      }
      PsiPsiBarMatrix[fInd1][fInd2] = fac*scalar;
    }
  }
  return success;
}
