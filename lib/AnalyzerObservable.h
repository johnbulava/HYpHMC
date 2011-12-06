#ifndef AnalyzerObservable_included
#define AnalyzerObservable_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "Global.h"
#include "AnalyzerIOControl.h"
#include "FermionMatrixOperations.h"
#include "StateDescriptorReader.h"
#include "AnalyzerPhiFieldConfiguration.h"
#include "SimulationParameterSet.h"


class AnalyzerObservable {
private:
  AnalyzerIOControl* ioControl;
  char* nickName;
  char* obsName;
  int totalPhiConfCount;


protected:
  FermionMatrixOperations* fermiOps;
  StateDescriptorReader* SDReader;
  SimulationParameterSet* SimParaSet;
  double* analyzerResults;
  int analyzeEveryXXXconf;

  
  void ini(int resCount);
  
  virtual int getAnalyzerResultsCount() = 0;
 
public:
  AnalyzerObservable(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, char* oName, char* nick); 
  virtual ~AnalyzerObservable();
  
  void saveAnalyzerResults(int ID);
  bool isConfAnalyzed(int ID);
  bool shallConfBeAnalyzed(int ID);
  int getNextConfIDforAnalysis();  //-1 if none
  bool isNick(char* nick);
  char* getObsName();
  char* getNickName();
  
  virtual bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) = 0;
  virtual int getNeededAuxVectorCount() = 0;
};


#include "AnalyzerObservable.C"

#endif
