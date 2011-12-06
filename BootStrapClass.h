#ifndef BootStrapClass_included
#define BootStrapClass_included

#include "Global.h"
#include "Tools.C"

class BootStrapClass {
private:
double* DataPool;
long double* DataWeightPool;
int DataPoolSize;
int DataCount;
int DataPoolSizeIncrement;
char* name;
char* fileName;

public:
  BootStrapClass();
  ~BootStrapClass();
  void addToDataPool(double x, long double w);
  void resetPool();
  void setName(char* fName); 
  void saveData();
  void loadData(int setSize, int Jstart, int Jend);
  void printData();
  double getAverage();
  int getDataCount();
  void doBootStrapWithWeights(int iter, double& avg, double& sigma, double& SusAvg, double& SusSigma, double& BinderAvg, double& BinderSigma);
};

#include "BootStrapClass.C"

#endif
