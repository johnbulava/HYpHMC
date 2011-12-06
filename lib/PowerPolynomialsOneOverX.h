#ifndef PowerPolynomialsOneOverX_included
#define PowerPolynomialsOneOverX_included

#include "Complex.h"
#include "Global.h"
#ifdef useBLAS
  #include CBLASIncludeFile
#endif
#include "Tools.h"

class PowerPolynomialsOneOverX {
protected:
  void ini(int setN, double setEps, double setAl, double setBet);
  void desini();
  void calcPCoeff();
  void calcAppCoeff();
  int N;
  double epsilon;
  double alpha;
  double beta;
  long double** Pcoeff;
  
  
public:
  PowerPolynomialsOneOverX();
  PowerPolynomialsOneOverX(int setN, double setEps, double setAl, double setBet);
  ~PowerPolynomialsOneOverX();

  long double* AppCoeff;

  void recalcCoeff(int setN, double setEps, double setAl, double setBet);
  long double PScalar(long double* p1, long double* p2);
  long double PScalarWithOneOverX(long double* p);
  
  void printToFile();
  double calcControlPlot();
};

#endif
