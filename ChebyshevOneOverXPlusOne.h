#ifndef ChebyshevOneOverXPlusOne_included
#define ChebyshevOneOverXPlusOne_included

#include "Complex.h"
#include "Global.h"
#include <stdio.h>
#include "Tools.h"

#ifdef useBLAS
  #include CBLASIncludeFile
#endif

class ChebyshevOneOverXPlusOne {
protected:
  void ini(int setN, int setn, double setEps);
  void desini();
  void calcCoeff();
  
  
public:
  ChebyshevOneOverXPlusOne();
  ChebyshevOneOverXPlusOne(int setN, int setn, double setEps);
  ~ChebyshevOneOverXPlusOne();
  int N;
  int n;
  double epsilon;
  double* coeff;
  void recalcCoeff(int setN, int setn, double setEps);
  void printToFile();
  double calcControlPlot();
};


#endif
