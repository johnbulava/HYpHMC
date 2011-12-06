#ifndef GeneralChebyshevApproximation_included
#define GeneralChebyshevApproximation_included

#include <math.h>
#include "Global.h"
#include "ComplexMatrix.h"

class GeneralChebyshevApproximation {
private:
  int coeffCount;
  double alpha;  //Recurrence - Parameter for Chebyshev - Polynomials
  double beta;   //Recurrence - Parameter for Chebyshev - Polynomials
  double* gamma;  //Chebyshev-Coefficients
  double relAccuracy;
  
  void clearData();
  double evaluatePolynomial(double x, bool withRenorm);
  bool calcApproximationStandard(double (*func)(double x), double minX, double maxX, double relAcc, int scanPoints, int& maxCoeffCount);
  bool calcApproximationRemez(double (*func)(double x), double minX, double maxX, double relAcc, int scanPoints, int& maxCoeffCount);
  bool calcApproximationKrylovBased(double (*func)(double x), double minX, double maxX, double relAcc, int scanPoints, int& maxCoeffCount);

public:
  GeneralChebyshevApproximation(); 
  ~GeneralChebyshevApproximation();
  
  void calcApproximation(double (*func)(double x), double minX, double maxX, double relAcc, int scanPoints);
  
  int getCoeffCount();
  double getAlpha();
  double getBeta();
  double* getGammas();
  double getRelAccuracy();
  void plotPolynomial(char*fileName, double (*func)(double x), double minX, double maxX, int scanPoints);
  double evaluatePolynomial(double x);
};

#endif
