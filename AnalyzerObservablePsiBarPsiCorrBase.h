#ifndef AnalyzerObservablePsiBarPsiCorrBase_included
#define AnalyzerObservablePsiBarPsiCorrBase_included

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


class AnalyzerObservablePsiBarPsiCorrBase : public AnalyzerObservable {
private:
  void sampleFermioncScanVector(Complex* v, int t, int FermionIndex, int timeDirection);
  bool calcPsiPsiBarMatrixForPhiField(double* phiField, Complex**** PsiPsiBarMatrix, double TOL, int LargestL, int timeDirection, int projMode, int fInd1Start, int fInd1End, int fInd2Start, int fInd2End);
  Complex* LeftVector;
  Complex* RightVector;
  Complex* SolutionVector;


protected:
  bool fixGauge;
  bool randomGauge;
  int projectorSelection;
  int multiplyWithPhiMatBSelection;
  int tresSumStartIndex;
  int tresSumEndIndex;
  int PsiPsiBarMatrixInd1Start;
  int PsiPsiBarMatrixInd2Start;
  int PsiPsiBarMatrixInd1End;
  int PsiPsiBarMatrixInd2End;

public:
  AnalyzerObservablePsiBarPsiCorrBase(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader, char* oName, char* nick); 
  ~AnalyzerObservablePsiBarPsiCorrBase();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#include "AnalyzerObservablePsiBarPsiCorrBase.C"

#endif
