#ifndef HighPrecisionComplexPolynom_included
#define HighPrecisionComplexPolynom_included


#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <time.h>
#include <sys/stat.h> 
#include <unistd.h>
#include <cmath>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/shm.h>
#include <signal.h>
#include <cln/number.h>
#include <cln/float.h>
#include <cln/complex.h>
#include <cln/real.h>
#include <cln/io.h>
#include <cln/integer_io.h>
#include <cln/float_io.h>
#include <cln/complex_io.h>
#include <cln/univpoly.h>
#include <cln/univpoly_integer.h>
#include <cln/univpoly_complex.h>
#include <cln/univpoly_real.h>

#include "Global.h"
#include "Complex.h"


class HighPrecisionComplexPolynom {
private:
  cln::float_format_t clnDIGIT;  
  cln::cl_F  ONE;
  cln::cl_F  TWO;
  cln::cl_F  ZERO;
  cln::cl_F  HALF;
  
  void ini(int len, int digit);
  void ini(const HighPrecisionComplexPolynom& pol);  
  void desini();
  cln::cl_N Lasolv(cln::cl_N* Poly, int Maxpow, cln::cl_N root, const int itemax=100);
  void Polyrootc(cln::cl_N* Poly, int Maxpow, cln::cl_N* Root);
  cln::cl_N EvalPoly(cln::cl_N* Poly, int Maxpow, cln::cl_N Valu);
  
public:
  HighPrecisionComplexPolynom(int len, int digit); 
  HighPrecisionComplexPolynom(const HighPrecisionComplexPolynom& pol);  
  ~HighPrecisionComplexPolynom();
  
  Complex evaluateAt(Complex z);
  cln::cl_N evaluateAt(cln::cl_N z);
  
  HighPrecisionComplexPolynom calcDerivative();
  void print();

  int DIGIT; 
  int length;
  cln::cl_N* coeff;
  
  int getLength();
  void setCoeff(int p, Complex val);
  void addCoeff(int p, Complex val);
  void setCoeff(int p, cln::cl_N val);
  void addCoeff(int p,cln::cl_N val);
  cln::cl_N getCoeff(int p);
  int getOrder();

  Complex* getRoots();
  cln::cl_N* getPrecRoots();
  
  HighPrecisionComplexPolynom operator + (HighPrecisionComplexPolynom);
  HighPrecisionComplexPolynom operator - (HighPrecisionComplexPolynom);
  HighPrecisionComplexPolynom operator * (HighPrecisionComplexPolynom);
  HighPrecisionComplexPolynom& operator = (HighPrecisionComplexPolynom);  
};

#endif
