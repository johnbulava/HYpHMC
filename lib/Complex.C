#include "Complex.h"
#include <stdio.h>

Complex ComplexI;
Complex ComplexUnity;
Complex ComplexZero;


Complex::Complex() { }


Complex::Complex(double r, double i) {
  x = r;
  y = i;
}


Complex::~Complex() { }


void Complex::setValues(double r, double i) {
  x = r;
  y = i;
}


Complex exp(Complex x) {
  Complex res;
  double f = exp(x.x);

  res.x = f * cos(x.y);
  res.y = f * sin(x.y);

  return res;
}


void Complex::print() {
  if (x>=0) printf(" ");
  printf("%1.3f",x);
  if (y>=0) printf("+");
  printf("%1.3fi\n",y);  
}


Complex Complex::operator + (Complex B) {
  Complex res;

  res.x = x + B.x;
  res.y = y + B.y;

  return res;
}


Complex Complex::operator - (Complex B) {
  Complex res;

  res.x = x - B.x;
  res.y = y - B.y;

  return res;
}


Complex Complex::operator * (Complex B) {
  Complex res;

  res.x = x*B.x - y*B.y;
  res.y = x*B.y + y*B.x;

  return res;
}


Complex Complex::operator / (Complex B) {
  Complex res;
  double b = B.x*B.x + B.y*B.y;

  res.x = (x*B.x + y*B.y)/b;
  res.y = (y*B.x - x*B.y)/b;

  return res;
}


Complex operator + (double A, Complex B) {
  Complex res;

  res.x = A + B.x;
  res.y = B.y;

  return res;
}


Complex operator + (Complex B, double A) {
  Complex res;

  res.x = A + B.x;
  res.y = B.y;

  return res;
}


Complex operator - (double A, Complex B) {
  Complex res;

  res.x = A - B.x;
  res.y = - B.y;

  return res;
}


Complex operator - (Complex A, double B) {
  Complex res;

  res.x = A.x - B;
  res.y = A.y;

  return res;
}


Complex operator * (double A, Complex B) {
  Complex res;

  res.x = A*B.x;
  res.y = A*B.y;

  return res;
}


Complex operator * (Complex B, double A) {
  Complex res;

  res.x = A*B.x;
  res.y = A*B.y;

  return res;
}


Complex operator / (double A, Complex B) {
  Complex res;
  double b = B.x*B.x + B.y*B.y;

  res.x = (A*B.x)/b;
  res.y = (- A*B.y)/b;

  return res;
}


Complex operator / (Complex A, double B) {
  Complex res;

  res.x = A.x/B;
  res.y = A.y/B;
  
  return res;
}


Complex adj(Complex c) {
  Complex res;
  
  res.x = c.x;
  res.y = -c.y;
  
  return res;
}


Complex sqrt(Complex c) {
  Complex res(0,0);
  double norm = sqrt(c.x*c.x + c.y*c.y);
  if (norm > 0) {
    double winkel = acos(c.x/norm);
    if (c.y < 0) winkel = -winkel;
    res.x = cos(winkel/2) * sqrt(norm);
    res.y = sin(winkel/2) * sqrt(norm);
  }
  return res;
}


Complex pow(Complex c, double p) {
  Complex res(1,1);
  double norm = sqrt(c.x*c.x + c.y*c.y);
  if (norm > 0) {
    double lenFac = exp(p*log(norm));
    double winkel = acos(c.x/norm);
    if (c.y < 0) winkel = -winkel;
    res.x = cos(winkel*p) * lenFac;
    res.y = sin(winkel*p) * lenFac;
  } else {
    if (p>=0) {
      res = 0.0*res;
    } else {
      double x = 0.0;
      double NaN = 0.0/x;
      res = NaN * res;
    }
  }  
  return res;
}


double normSqr(Complex c) {
  return c.x*c.x + c.y*c.y;
}


double norm(Complex c) {
  return sqrt(c.x*c.x + c.y*c.y);
}


Complex log(Complex c) {
  Complex res(0,0);
  double norm = sqrt(c.x*c.x + c.y*c.y);
  if (norm > 0) {
    double winkel = acos(c.x/norm);
    if (c.y < 0) winkel = -winkel;
    res.x = log(norm);
    res.y = winkel;
    return res;    
  }
  double x = 0;
  double NaN = 0.0/x;
  res.x = NaN;
  res.y = NaN;
  return res;  
}


Complex arctanh(Complex c) {
  Complex dummy1(1+c.x, c.y);
  Complex dummy2(1-c.x, -c.y);
  
  return 0.5*log(dummy1/dummy2);
}


Complex arctanhAsPowerSeries(Complex c) {
  if (norm(c)<1E-13) return Complex(0,0);
  Complex res = c;
  Complex pow = c;
  Complex cSqr = c*c;
  Complex oldRes(0,0);
  int k=0;
  while (norm(res-oldRes)>1E-13) {
    oldRes = res;
    k++;
    double fac = 1.0 / (2.0*k+1.0);
    pow = pow * cSqr;
    res = res + fac * pow;
  }
  return res;
}
