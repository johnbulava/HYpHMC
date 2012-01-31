#include "AnalyzerObservablePsiBarPhiPsiCondensateBase.h"

AnalyzerObservablePsiBarPhiPsiCondensateBase::AnalyzerObservablePsiBarPhiPsiCondensateBase(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader, char* oName, char* nick) : AnalyzerObservable(fOps, aIOcon, SDreader, oName, nick) { 
  fixGauge = false;
  randomGauge = false;
  projectorSelection = -1;
  multiplyWithPhiMatBSelection = -1;

  ini(getAnalyzerResultsCount());
}


AnalyzerObservablePsiBarPhiPsiCondensateBase::~AnalyzerObservablePsiBarPhiPsiCondensateBase() {
}


bool AnalyzerObservablePsiBarPhiPsiCondensateBase::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
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
  HelperVector = auxVectors[3];

  fermiOps->setPreconditioner(true, 1.1, 0);
  fermiOps->setQPreconditioner(true, 0.25, 0.25);    
 
  Complex PsiBarPhiPsiCond(0,0);  

  bool success = true;
  for (int row=0; row<8; row++) {
    for (int col=0; col<8; col++) {
      if (abs(row-col) % 4 == 0) {
        sampleFermionRightVector(RightVector, col);

        bool daggeredPhi = (multiplyWithPhiMatBSelection > 0);
        sampleFermionLeftVector(LeftVector, row, phiField, daggeredPhi, RightVector, col);

        Complex result(NaN, NaN);
        bool b = calcMatrixElementProjectorHatMinverseProjector(result, phiField, LeftVector, RightVector, TOL, projectorSelection);
        success = success & b;

        PsiBarPhiPsiCond = PsiBarPhiPsiCond + result;
      }
    }
  }

  analyzerResults[0] = PsiBarPhiPsiCond.x / volume;
  analyzerResults[1] = PsiBarPhiPsiCond.y / volume;

  return success;
}


int AnalyzerObservablePsiBarPhiPsiCondensateBase::getNeededAuxVectorCount() {
  return 4;
}


int AnalyzerObservablePsiBarPhiPsiCondensateBase::getAnalyzerResultsCount() {
  return 2;
}


void AnalyzerObservablePsiBarPhiPsiCondensateBase::sampleFermionRightVector(Complex* r, int FermionIndex) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double fac = 1.0 / sqrt(2.0);

  fermiOps->zeroFermionVector(r);

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

          AdvancedGaussZufall(AdvancedSeed, r[index].x, r[index].y);
          r[index].x *= fac;
          r[index].y *= fac;
        }
      }
    }
  }
}


void AnalyzerObservablePsiBarPhiPsiCondensateBase::sampleFermionLeftVector(Complex* l, int FermionIndex_l, double* phiField, bool daggeredPhi, Complex* r, int FermionIndex_r) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  fermiOps->zeroFermionVector(l);

  int qm1 = FermionIndex_l / 4;
  int qm2 = FermionIndex_r / 4;
  if (abs(FermionIndex_l-FermionIndex_r) % 4 != 0) return;
  
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

          int phiIndex = 0;
          phiIndex += 4 * ind[3];
          phiIndex += 4 * L3 * ind[2];
          phiIndex += 4 * L3*L2  * ind[1];
          phiIndex += 4 * L3*L2*L1 * ind[0];
          
          double f = 1;
          if (daggeredPhi) f = -1;
          Quat q(phiField[phiIndex+0],f*phiField[phiIndex+1],f*phiField[phiIndex+2],f*phiField[phiIndex+3]);

          ComplexMatrix qm(q);

          l[index + FermionIndex_l] = adj(qm.matrix[qm2][qm1]) * r[index + FermionIndex_r];
        }
      }
    }
  }
}


/****
* This routine computes <LeftVector | \hat P_\pm (\D+B(1-(1/\rho)\D))^-1 P_\pm | RightVector>
* All factors are already included.
*/
bool AnalyzerObservablePsiBarPhiPsiCondensateBase::calcMatrixElementProjectorHatMinverseProjector(Complex& result, double* phiField, Complex* LeftVector, Complex* RightVector, double TOL, int projMode) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  bool success = true;  
  double rho, r;
  fermiOps->getDiracParameters(rho, r);
  double fac = -1.0 / (2.0*rho);
  if (projMode != 0) fac *= 0.25;  //Missing factor 0.5 in projector \hat P_\pm and P_\pm compensated with this factor 0.25. 

  //Berechne daggered Projector \hat P_\pm angewendet auf LeftVector  ==> Faktor 0.5 fehlt
  if (projMode != 0) {
    fermiOps->executeGamma5(LeftVector, SolutionVector);
    fermiOps->executeDiracDaggerMatrixMultiplication(SolutionVector, HelperVector, false);
    Complex alpha1(-1.0/rho, 0);
    SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, alpha1, HelperVector, SolutionVector);
    Complex alpha2(projMode, 0);
    SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, alpha2, SolutionVector, LeftVector);	
  }
    
  //Apply inverse M^dagger  
  int neededIter = 0;
  bool b = fermiOps->executeFermionMatrixDaggerInverseMatrixMultiplication(LeftVector, SolutionVector, phiField, TOL, false, -1, neededIter);
  if (!b) success = false;

  //Apply P_\pm^dagger  ==> Faktor 0.5 fehlt
  if (projMode != 0) {
    bool projPlus = true;
    if (projMode<0) projPlus = false;
    fermiOps->executeProjectorMultiplication(SolutionVector, SolutionVector, projPlus); //OHNE Faktor 0.5
  }

  Complex scalar(0,0);     
	  
  SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, SolutionVector, RightVector, scalar);
  result = fac*scalar;

  return success;
}
