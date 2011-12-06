#ifndef ChebyshevOneOverXPlusOneRelError_included
#define ChebyshevOneOverXPlusOneRelError_included

#include "Complex.h"
#include "Global.h"
#ifdef useBLAS
  #include CBLASIncludeFile
#endif

class ChebyshevOneOverXPlusOneRelError {
protected:
  void ini(int setN, double setEps);
  void desini();
  void calcCoeff();
  double epsilon;
  
  
public:
  ChebyshevOneOverXPlusOneRelError();
  ChebyshevOneOverXPlusOneRelError(int setn, double setEps);
  ~ChebyshevOneOverXPlusOneRelError();
  int n;
  double* coeff;
  void recalcCoeff(int setn, double setEps);
  void printToFile();
  double calcControlPlot();
};


#endif
