#ifndef LatticeSiteBins_included
#define LatticeSiteBins_included

#include <stdlib.h>

#include "Global.h"
#include "Complex.h"

#define LatticeSiteBinsDataMAX 1000000
#define LatticeSiteBinsMomentumSqrSlotSize 1E-9

class LatticeSiteBins {
private:
  void ini(int l0, int l1, int l2, int l3);
  void desini();
  
  int L0,L1,L2,L3;
  double* sinPSqr;
  double* dataSum;
  double* dataSqrSum;
  int dataCount;
  int independentDataSetCount;
    
public:
  LatticeSiteBins(int l0, int l1, int l2, int l3);
  ~LatticeSiteBins();

  void clearData();
  void addDataVectorFromOmegafourierTrafo(Complex* data);  
  void getAverageVector(double* avg);
  void getSigmaVector(double* sig);
  void getAvgSigmaVector(double* avgsig);
  void saveData(char* fileName);
  void loadData(char* fileName);
  double getLatMomSqrFromIndex(int index);
  bool containsData();
};

#endif
