#ifndef EvaluateObservablePsiBarPsiPropagatorBase_included
#define EvaluateObservablePsiBarPsiPropagatorBase_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "MassCorrelationMatrixAnalyzer.h"
#include "AutoCorrelation.h"
#include "NeubergerMatrix.h"


class EvaluateObservablePsiBarPsiPropagatorBase : public EvaluateObservable {
private:  
  LAPsystemPlot* createPlot1(int dataShift, double xRange);
  double* PropagatorValues;
  double* PropagatorErrors;
  double* FreePropagatorValues;
  double* pSqr;
  double* pHatSqr;
  
  double redChiSqr_Lin;
  double fittedMass_Lin;
  double fittedZ_Lin;
  double fittedMass_Lin_Error;
  double fittedZ_Lin_Error;

  double redChiSqr_FreeAna;
  double fittedMass_FreeAna;
  double fittedZ_FreeAna;
  double fittedMass_FreeAna_Error;
  double fittedZ_FreeAna_Error;
  
  Complex freeInvPropAnalytical(vector4D pvec, double mass, double Zpsi);
  
  
protected:
  double fitThreshold_Lin;
  double kickOutPHatThreshold_Lin;
  double fitThreshold_FreeAna;
  double kickOutPHatThreshold_FreeAna;
  double RescaleFactor;
  

  
public:    
  EvaluateObservablePsiBarPsiPropagatorBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, char* oName, char* nick, double relStart, double relEnd); 
  ~EvaluateObservablePsiBarPsiPropagatorBase();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      

  double chiSqrForFreeAnaFit(double mass, double Zpsi);

  
  double getPropagatorZFactor();
  double getPropagatorZFactorError();
  double getPropagatorMass();
  double getPropagatorMassError();  
};


#endif
