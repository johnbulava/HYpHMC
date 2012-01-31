#ifndef AnalyzerObservablePsiBarPhiPsiCondensateBase_included
#define AnalyzerObservablePsiBarPhiPsiCondensateBase_included

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


class AnalyzerObservablePsiBarPhiPsiCondensateBase : public AnalyzerObservable {
private:
  void sampleFermionRightVector(Complex* r, int FermionIndex);
  void sampleFermionLeftVector(Complex* l, int FermionIndex_l, double* phiField, bool daggeredPhi, Complex* r, int FermionIndex_r);
  bool calcMatrixElementProjectorHatMinverseProjector(Complex& result, double* phiField, Complex* LeftVector, Complex* RightVector, double TOL, int projMode);

  Complex* LeftVector;
  Complex* RightVector;
  Complex* SolutionVector;
  Complex* HelperVector;


protected:
  bool fixGauge;
  bool randomGauge;
  int projectorSelection;
  int multiplyWithPhiMatBSelection;

public:
  AnalyzerObservablePsiBarPhiPsiCondensateBase(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader, char* oName, char* nick); 
  ~AnalyzerObservablePsiBarPhiPsiCondensateBase();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#endif
