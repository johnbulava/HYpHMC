#ifndef Complex_included
#define Complex_included

#include <math.h>

class Complex {
private:


public:
  double x;
  double y;
  Complex(); 
  Complex(double r, double i);
  ~Complex();
  
  void setValues(double r, double i);

  void print();
  
  Complex operator + (Complex);
  Complex operator - (Complex);  
  Complex operator * (Complex);
  Complex operator / (Complex);
};

/*
Complex ComplexI(0,1);
Complex ComplexUnity(1,0);
Complex ComplexZero(0,0);
*/

extern Complex ComplexI;
extern Complex ComplexUnity;
extern Complex ComplexZero;

Complex pow(Complex c, double p);
double normSqr(Complex c);
double norm(Complex c);
Complex log(Complex c);
Complex arctanh(Complex c);
Complex arctanhAsPowerSeries(Complex c);
Complex adj(Complex c);
Complex operator + (double A, Complex B);
Complex operator + (Complex B, double A);
Complex operator - (double A, Complex B);
Complex operator - (Complex A, double B);
Complex operator * (double A, Complex B);
Complex operator * (Complex B, double A);
Complex operator / (double A, Complex B);
Complex operator / (Complex A, double B);
Complex sqrt(Complex c);
Complex exp(Complex x);

#endif
