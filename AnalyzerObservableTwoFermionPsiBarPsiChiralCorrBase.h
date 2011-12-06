#ifndef AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase_included
#define AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "Global.h"
#include "AnalyzerIOControl.h"
#include "FermionMatrixOperations.h"
#include "StateDescriptorReader.h"
#include "AnalyzerObservable.h"
#include "AnalyzerPhiFieldConfiguration.h"


class AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase : public AnalyzerObservable {
private:
  void sampleFermioncScanVector(Complex* v, int t, int FermionIndex, int timeDirection);
  bool calcTwoFPsiPsiBarMatrixForPhiField(double* phiField, Complex***** PsiPsiBarMatrix, double TOL, int LargestL, int timeDirection, int projMode);
  Complex* auxVec0;
  Complex* auxVec1;
  Complex* auxVec2;
  Complex* auxVec3;
  

protected:
  bool fixGauge;
  bool randomGauge;
  int projectorSelection;
  bool multiplyWithPhiMatBSelection;

public:
  AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader, char* oName, char* nick); 
  ~AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableTwoFermionPsiBarPsiChiralCorrBase.C"

#endif
