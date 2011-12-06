#ifndef ComplexPolynom_included
#define ComplexPolynom_included

#include <math.h>
#include "Complex.h"

class ComplexPolynom {
private:
  void ini(int len);
  void ini(const ComplexPolynom& pol);  
  void desini();
  
public:
  ComplexPolynom(int len); 
  ComplexPolynom(const ComplexPolynom& pol);  
  ~ComplexPolynom();
  
  Complex evaluateAt(Complex z);
  ComplexPolynom calcDerivative();

  int length;
  Complex* coeff;
  
  int getLength();
  void setCoeff(int p, Complex val);
  void addCoeff(int p, Complex val);
  Complex getCoeff(int p);
  void print();
  
  ComplexPolynom operator + (ComplexPolynom);
  ComplexPolynom operator - (ComplexPolynom);
  ComplexPolynom operator * (ComplexPolynom);
  ComplexPolynom& operator = (ComplexPolynom);  
};

#endif
