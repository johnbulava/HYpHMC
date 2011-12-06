#include "ComplexPolynom.h"
#include <stdlib.h>
#include "Tools.h"


ComplexPolynom::ComplexPolynom(int len) { 
  ini(len);
}


ComplexPolynom::ComplexPolynom(const ComplexPolynom& pol) { 
  ini(pol);
}


ComplexPolynom::~ComplexPolynom() { 
  desini();
}


void ComplexPolynom::ini(int len) {
  length = len;
  coeff = new Complex[len];
  for (int I=0; I<len; I++) {
    coeff[I].x = 0;
    coeff[I].y = 0;
  }
}


void ComplexPolynom::ini(const ComplexPolynom& pol) {
  ini(pol.length);
  for (int I=0; I<pol.length; I++) {
    Complex c = pol.coeff[I];
    setCoeff(I, c);
  }
}


void ComplexPolynom::desini() {
  delete[] coeff;
  coeff = NULL;
  length = 0;
}


Complex ComplexPolynom::evaluateAt(Complex z) {
  Complex res(0,0);
  Complex dummy(1,0);
  for (int I=0; I<length; I++) {
    res = res + coeff[I]*dummy;
    dummy = dummy * z;
  }
  return res;
}


ComplexPolynom ComplexPolynom::calcDerivative() {
  int len = length-1;
  if (length < 1) len = 1;
  ComplexPolynom pol(len);
  for (int I=0; I<length-1; I++) {
    pol.setCoeff(I, (I+1.0)*coeff[I+1]);
  }
  return pol;
}


void ComplexPolynom::print() {
  for (int I=0; I<length; I++) {
    coeff[I].print();
  }
}

  
int ComplexPolynom::getLength() {
  return length;
}


void ComplexPolynom::setCoeff(int p, Complex val) {
  if ((p<0) || (p>=length)) return;
  
  coeff[p] = val;
}


void ComplexPolynom::addCoeff(int p, Complex val) {
  if ((p<0) || (p>=length)) return;
  
  coeff[p] = coeff[p] + val;
}


Complex ComplexPolynom::getCoeff(int p) {
  if ((p<0) || (p>=length)) return Complex(NaN, NaN);
  return coeff[p];
}
 
  
ComplexPolynom ComplexPolynom::operator + (ComplexPolynom pol) {
  ComplexPolynom res(length);
  for (int I=0; I<length; I++) {
    res.setCoeff(I, coeff[I] + pol.getCoeff(I));
  }
  return res;
} 


ComplexPolynom ComplexPolynom::operator - (ComplexPolynom pol) {
  ComplexPolynom res(length);
  for (int I=0; I<length; I++) {
    res.setCoeff(I, coeff[I] - pol.getCoeff(I));
  }
  return res;
}


ComplexPolynom ComplexPolynom::operator * (ComplexPolynom pol) {
  ComplexPolynom res(length+pol.getLength()-1);
  for (int I=0; I<length; I++) {
    for (int I2=0; I2<pol.getLength(); I2++) {
      res.addCoeff(I+I2, coeff[I] * pol.getCoeff(I2));
    }
  }
  return res;
}


ComplexPolynom& ComplexPolynom::operator = (ComplexPolynom pol) {
  desini();
  ini(pol);

  return *this;
}
