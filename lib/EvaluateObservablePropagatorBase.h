#ifndef EvaluateObservablePropagatorBase_included
#define EvaluateObservablePropagatorBase_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>

#include "Global.h"
#include "EvaluateObservable.h"
#include "LatticeMomentumBins.h"
#include "AutoCorrelation.h"
#include "HighPrecisionComplex.h"


class EvaluateObservablePropagatorBase : public EvaluateObservable {
private:  
  LatticeMomentumBins* latticeBins;
  AutoCorrelation* autoCorr;  
  double* pSqr;
  double* avgProp;
  double* sigmaProp;
  double* autoCorrelationTime;
  double vevContNot;
  bool kappaZeroMode;
  double PropagatorEuclideanZFactor;
  double PropagatorEuclideanZFactorError;
  double PropagatorEuclideanMass;
  double PropagatorEuclideanMassError;
  double PropagatorMinkowskiZFactor;
  double PropagatorMinkowskiZFactorError;
  double PropagatorMinkowskiMass;
  double PropagatorMinkowskiMassError;  
  double PropagatorFitConst0;
  double PropagatorFitConst0Error;
  double PropagatorFitConst1;
  double PropagatorFitConst1Error;
  double PropagatorFitConst2;
  double PropagatorFitConst2Error;
  
  

  double PropagatorEuclideanZFactorReduced;
  double PropagatorEuclideanZFactorErrorReduced;
  double PropagatorEuclideanMassReduced;
  double PropagatorEuclideanMassErrorReduced;
  double PropagatorMinkowskiZFactorReduced;
  double PropagatorMinkowskiZFactorErrorReduced;
  double PropagatorMinkowskiMassReduced;
  double PropagatorMinkowskiMassErrorReduced;  
  double PropagatorFitConst0Reduced;
  double PropagatorFitConst0ErrorReduced;
  double PropagatorFitConst1Reduced;
  double PropagatorFitConst1ErrorReduced;
  double PropagatorFitConst2Reduced;
  double PropagatorFitConst2ErrorReduced;


  double PropagatorFitConst0ReducedPole;
  double PropagatorFitConst0ErrorReducedPole;
  double PropagatorFitConst1ReducedPole;
  double PropagatorFitConst1ErrorReducedPole;
  double PropagatorFitConst2ReducedPole;
  double PropagatorFitConst2ErrorReducedPole;
  double PropagatorFitConst3ReducedPole;
  double PropagatorFitConst3ErrorReducedPole;
  

  double PropagatorZFactor;
  double PropagatorZFactorError;  
  double PropagatorZ0Factor;
  double PropagatorZ0FactorError;  
  double PropagatorMProp;
  double PropagatorMPropError;  
  double PropagatorMProp0;
  double PropagatorMProp0Error;  

  
  
  double PropagatorPoleMass;
  double PropagatorPoleMassError;  
  double PropagatorPoleDecayWidth;
  double PropagatorPoleDecayWidthError;  
  double PropagatorPoleValue;
  double PropagatorPoleValueError;  
  double PropagatorFitReducedChiSquare;
  double PropagatorFitReducedChiSquarePole;
  double PropagatorFitReducedChiSquareReduced;
  
  
  int L0;
  int L1;
  int L2;
  int L3;
  
  void calcPropAtP0(int ignoreStart, int ignoreEnd, double &avgProp0);
  void calcProp(int ignoreStart, int ignoreEnd, double &EucZ, double &Eucmass, double &MinkZ, double &Minkmass, double &fC0, double &fC1, double &fC2,  double &EucZReduced, double &EucmassReduced, double &MinkZReduced, double &MinkmassReduced, double &fC0Reduced, double &fC1Reduced, double &fC2Reduced, double &Z, double &Z0, double &mprop, double &mprop0, Complex& pole, Complex& poleVal, double &fC0P, double &fC1P, double &fC2P, double &fC3P, double &redChiSqr, double &redChiSqrReduced, double &redChiSqrPole);
  LAPsystemPlot* createPlot1(double maxP);
  LAPsystemPlot* createPlot2(double maxP, double subtractMSqr, bool considerSubTerm);
  LAPsystemPlot* createPlot3();

  void calcAvgAndSigmaFromSelectedDataBuffer(int N, int ind, double** dataBuffer, double& avg, double& sig);

protected:
  bool considerZeroMomentum;
  double selfEnergySubLamFac;
  bool doArcTanhFit;
  double gamma;
  bool HiggsPropFourParameterFit;
  
public:    
  EvaluateObservablePropagatorBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, const char* oName, const char* nick, double relStart, double relEnd); 
  ~EvaluateObservablePropagatorBase();

  bool evaluate();
  void generateLatexAndPlotsAndXML();  
  void defineObsDependencies();      
  int getAnalyzerResultsCount();      

  double getPropagatorZFactor();
  double getPropagatorZFactorError();
  double getPropagatorMass();
  double getPropagatorMassError();
  double getPropagatorEuclideanZFactor();
  double getPropagatorEuclideanZFactorError();
  double getPropagatorEuclideanMass();
  double getPropagatorEuclideanMassError();  
  double getPropagatorMinkowskiZFactor();
  double getPropagatorMinkowskiZFactorError();
  double getPropagatorMinkowskiMass();
  double getPropagatorMinkowskiMassError();
};


#endif
