#ifndef AnalyzerHiggs_included
#define AnalyzerHiggs_included

#include <stdlib.h>

#include "Global.h"
#include "Complex.h"
#include "Quat.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include FFTWIncludeFile
#include "AutoCorrelation.h"
#include "ControlLogger.h"
#include "MassCorrelationMatrixAnalyzer.h"

#define AnalyzerHiggsDataMAX 1000000
#define AnalyzerHiggsMomentumSqrSlotSize 1E-9

class AnalyzerHiggs {
private:
  void ini(int l0, int l1, int l2, int l3, ControlLogger* log, bool FLAG_GFit);
  void desini();
  void measureMags(vector4D* phiField);
  int findMomentumSqrSlot(double momSqr);
  void plotSlottedGoldstonePropagator(int ignoreStart, int ignoreEnd);
  void plotSlottedHiggsPropagator(int ignoreStart, int ignoreEnd);
  void plotHiggsTimeSliceCorrelator(int ignoreStart, int ignoreEnd);
  
  ControlLogger* MassControlLog;
  int L0,L1,L2,L3;
  int timeDirection;
  int spaceDirection1;
  int spaceDirection2;
  int spaceDirection3;
  int LargestL;
  Complex* phiFieldBuffer;
  Complex* phiMomentumBuffer;
  fftw_plan phiFieldFourierForwardPlan1Components;
  fftw_plan phiFieldFourierForwardPlan4Components;
  double* sinPSqr;
  AutoCorrelation** HiggsTimeSliceCorrelator;
  AutoCorrelation** GoldstoneTimeSliceCorrelator;
  AutoCorrelation* HiggsVEV;
  double* HiggsVEVdata;
  double** HiggsTimeSliceData;
  double** GoldstoneTimeSliceData;
  double** HiggsTimeSliceCorrelatorData;
  double** GoldstoneTimeSliceCorrelatorData;
  double* weightData;
  double* HiggsPropagator;
  double* GoldstonePropagator;
  double* MomentumSqrSlotLocations;
  int MomentumSqrSlotCount;
  double** HiggsPropagatorSlottedData;
  double** GoldstonePropagatorSlottedData;
  AutoCorrelation** HiggsPropagatorSlottedDataCorrelation;
  AutoCorrelation** GoldstonePropagatorSlottedDataCorrelation;    
  int totalN;
  double measurePhiNorm;
  double measureStaggeredPhiNorm;
  bool FLAG_GnuplotFit;
  MassCorrelationMatrixAnalyzer* GoldstoneTwoParticleMassAnalyzer;
  MassCorrelationMatrixAnalyzer* HiggsGoldstoneMassAnalyzer;
  MassCorrelationMatrixAnalyzer* HiggsMassAnalyzer;
    
public:
  AnalyzerHiggs(int l0, int l1, int l2, int l3, ControlLogger* log, bool FLAG_GFit);
  ~AnalyzerHiggs();

  double LatticeResult_VEV;
  double LatticeResult_VEVsigma;
  double LatticeResult_HiggsPropagatorMass;
  double LatticeResult_HiggsPropagatorMassSigma;
  double LatticeResult_GoldstoneZFactor;
  double LatticeResult_GoldstoneZFactorSigma;
  double LatticeResult_PhysicalHiggsMass;
  double LatticeResult_PhysicalHiggsMassSigma;
  double* LatticeResult_PhysicalEffectiveHiggsMasses;
  double* LatticeResult_PhysicalEffectiveHiggsMassesSigmas;

  
  int getTotalN();
  void analyzeHiggsField(vector4D* phiField, double weight);

  double getLastMagnetization();
  double getLastStaggeredMagnetization();
  void plotGoldstonePropagator();
  void plotSlottedGoldstonePropagator();
  void plotHiggsPropagator();
  void plotSlottedHiggsPropagator();
  void plotHiggsTimeSliceCorrelator();
  void plotGoldstoneTimeSliceCorrelator();
  void calcHiggsTimeSliceCorrelator();
  void calcGoldstoneTimeSliceCorrelator();
  void calcHiggsVEV();
  void plotGoldstone2ParticleMasses();
  void plotHiggsGoldstoneMasses();
  void plotHiggsMasses();
};

#endif
