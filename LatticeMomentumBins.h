#ifndef LatticeMomentumBins_included
#define LatticeMomentumBins_included

#include <stdlib.h>

#include "Global.h"
#include "Complex.h"

#define LatticeMomentumBinsDataMAX 1000000
#define LatticeMomentumBinsMomentumSqrSlotSize 1E-9

class LatticeMomentumBins {
private:
  void ini();
  void desini();
  int findMomentumSqrSlot(double momSqr);
  void calcAverageAndSigmaVectors();
  
  int L0,L1,L2,L3;
  double* sinPSqr;
  double* MomentumSqrSlotLocations;
  int* LatticeMomentumToSlotPointer;
  int* LatticeNegCoorPointer;
  int* MomentumSlotMultiplicity;
  int MomentumSqrSlotCount;
  
  double* dataSum;
  double* dataSqrSum;
  int* dataCount;
  double* dataAvg;
  double* dataSigma;
  double* dataAvgSigma;  
  int independentDataSetCount;
  bool initialized;

    
public:
  LatticeMomentumBins(int l0, int l1, int l2, int l3);
  ~LatticeMomentumBins();

  void clearData();
  void addDataVector(double* data);
  void addDataVectorFromfourierTrafoSPECIAL(Complex* data);  
  void addDataVectorFromOmegafourierTrafoSPECIAL(Complex* data);  
  double* getAverageVector();
  double* getSigmaVector();
  double* getAvgSigmaVector();
  void getAverageVectorInflated(double* avg);
  void getSigmaVectorInflated(double* sig);
  void getAvgSigmaVectorInflated(double* avgsig);
  void saveData(char* fileName);
  void loadData(char* fileName);
  double getLatMomSqrFromIndex(int index);
  double getLatMomSqrFromSlotNr(int slotNr);
  int getMomentumSqrSlotCount();
  int getMomentumSlotFromIndex(int index);
  int getMomentumSlotMultiplicity(int slotNr);
  bool containsData();
};

#endif
