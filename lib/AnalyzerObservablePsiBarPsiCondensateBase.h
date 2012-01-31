#ifndef AnalyzerObservablePsiBarPsiCondensateBase_included
#define AnalyzerObservablePsiBarPsiCondensateBase_included

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


class AnalyzerObservablePsiBarPsiCondensateBase : public AnalyzerObservable {
private:
  void sampleFermioncScanVector(Complex* v, int FermionIndex);
  void sampleSource();
  bool calcPsiPsiBarMatrixForPhiField(double* phiField, Complex** PsiPsiBarMatrix, double TOL, int projMode, int fInd1Start, int fInd1End, int fInd2Start, int fInd2End);
  Complex* LeftVector;
  Complex* RightVector;
  Complex* SolutionVector;
  Complex* SourceVector;


protected:
  bool fixGauge;
  bool randomGauge;
  int projectorSelection;
  int tresSumStartIndex;
  int tresSumEndIndex;
  int PsiPsiBarMatrixInd1Start;
  int PsiPsiBarMatrixInd2Start;
  int PsiPsiBarMatrixInd1End;
  int PsiPsiBarMatrixInd2End;
  int numberOfMeasurements; 
  bool stochasticalSource;

public:
  AnalyzerObservablePsiBarPsiCondensateBase(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader, char* oName, char* nick); 
  ~AnalyzerObservablePsiBarPsiCondensateBase();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#endif
