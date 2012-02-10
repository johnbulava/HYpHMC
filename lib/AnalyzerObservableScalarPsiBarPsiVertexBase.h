#ifndef AnalyzerObservableScalarPsiBarPsiVertexBase_included
#define AnalyzerObservableScalarPsiBarPsiVertexBase_included

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


class AnalyzerObservableScalarPsiBarPsiVertexBase : public AnalyzerObservable {
private:
  void sampleFermioncScanVector(Complex* v, int pft, int FermionIndex);
  bool calcPsiPsiBarMatrixForPhiField(double* phiField, Complex**** PsiPsiBarMatrix, double TOL, int projMode, int fInd1Start, int fInd1End, int fInd2Start, int fInd2End);
  Complex* LeftVector;
  Complex* RightVector;
  Complex* SolutionVector;


protected:
  bool fixGauge;
  int projectorSelection;
  int tresSumStartIndex;
  int tresSumEndIndex;
  int PsiPsiBarMatrixInd1Start;
  int PsiPsiBarMatrixInd2Start;
  int PsiPsiBarMatrixInd1End;
  int PsiPsiBarMatrixInd2End;

public:
  AnalyzerObservableScalarPsiBarPsiVertexBase(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader, const char* oName, const char* nick); 
  ~AnalyzerObservableScalarPsiBarPsiVertexBase();
  
  bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
  int getNeededAuxVectorCount();
  int getAnalyzerResultsCount();  
};


#endif
