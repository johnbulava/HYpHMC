#include "HighPrecisionComplex.h"

HighPrecisionComplex::HighPrecisionComplex(int digit, double vx, double vy) { 
  ini(digit, vx, vy);
}


HighPrecisionComplex::HighPrecisionComplex(int digit, double vx) { 
  ini(digit, vx);
}


HighPrecisionComplex::HighPrecisionComplex(int digit) { 
  ini(digit);
}


HighPrecisionComplex::HighPrecisionComplex(int digit, cln::cl_N clnC) {
  ini(digit, clnC);
}


HighPrecisionComplex::HighPrecisionComplex(int digit, Complex c) {
  ini(digit, c.x, c.y);  
}


HighPrecisionComplex::HighPrecisionComplex(const HighPrecisionComplex& c) { 
  ini(c);
}


HighPrecisionComplex::~HighPrecisionComplex() { 
  desini();
}


void HighPrecisionComplex::setZero() {
  z = complex(ZERO, ZERO);
}


void HighPrecisionComplex::setOne() {
  z = complex(ONE, ZERO);
}


void HighPrecisionComplex::setTwo() {
  z = complex(TWO, ZERO);
}


void HighPrecisionComplex::setHalf() {
  z = complex(HALF, ZERO);
}


void HighPrecisionComplex::ini(int digit, double vx, double vy) {
  DIGIT = digit;
  clnDIGIT = cln::float_format(DIGIT);  
  
  char* xxxStr = new char[1000];
  snprintf(xxxStr,1000,"1.0e+0_%d",DIGIT);
  if (LogLevel>4) printf("Initializing ONE with: %s\n",xxxStr);
  ONE = xxxStr;

  snprintf(xxxStr,1000,"2.0e+0_%d",DIGIT);
  if (LogLevel>4) printf("Initializing TWO with: %s\n",xxxStr);
  TWO = xxxStr;

  snprintf(xxxStr,1000,"0.0e+0_%d",DIGIT);
  if (LogLevel>4) printf("Initializing ZERO with: %s\n",xxxStr);
  ZERO = xxxStr;

  snprintf(xxxStr,1000,"0.5e+0_%d",DIGIT);
  if (LogLevel>4) printf("Initializing HALF with: %s\n",xxxStr);
  HALF = xxxStr;
  delete[] xxxStr;
  
  z = complex(cl_float(vx,clnDIGIT), cl_float(vy,clnDIGIT));
}


void HighPrecisionComplex::ini(const HighPrecisionComplex& c) {
  ini(c.DIGIT, 0, 0);
  z = c.z;
}


void HighPrecisionComplex::ini(int digit, double vx) {
  ini(digit, 0, 0);  
  z = complex(cl_float(vx,clnDIGIT), ZERO);
}



void HighPrecisionComplex::ini(int digit, cln::cl_N clnC) {
  ini(digit, 0, 0);
  z = clnC;
}


void HighPrecisionComplex::ini(int digit) {
  ini(digit, 0, 0);
  z = complex(ZERO, ZERO);;
}



void HighPrecisionComplex::desini() {
}


void HighPrecisionComplex::print() {
  Complex res(double_approx(realpart(z)), double_approx(imagpart(z)));
  res.print();
} 


Complex HighPrecisionComplex::getComplex() {
  Complex res(double_approx(realpart(z)), double_approx(imagpart(z)));
  return res;
}


HighPrecisionComplex HighPrecisionComplex::operator + (HighPrecisionComplex c) {
  HighPrecisionComplex res(DIGIT, 0, 0);
  res.z = z + c.z;
  return res;
} 


HighPrecisionComplex HighPrecisionComplex::operator - (HighPrecisionComplex c) {
  HighPrecisionComplex res(DIGIT, 0, 0);
  res.z = z - c.z;
  return res;
}


HighPrecisionComplex HighPrecisionComplex::operator * (HighPrecisionComplex c) {
  HighPrecisionComplex res(DIGIT, 0, 0);
  res.z = z * c.z;
  return res;
}


HighPrecisionComplex HighPrecisionComplex::operator / (HighPrecisionComplex c) {
  HighPrecisionComplex res(DIGIT, 0, 0);
  res.z = z / c.z;
  return res;
}


HighPrecisionComplex& HighPrecisionComplex::operator = (HighPrecisionComplex c) {
  desini();
  ini(c);

  return *this;
}


HighPrecisionComplex exp(HighPrecisionComplex c) {
  HighPrecisionComplex res(c.DIGIT, exp(c.z));
  return res;
}


HighPrecisionComplex log(HighPrecisionComplex c) {
  HighPrecisionComplex res(c.DIGIT, log(c.z));
  return res;
}


HighPrecisionComplex sqrt(HighPrecisionComplex c) {
  HighPrecisionComplex res(c.DIGIT, sqrt(c.z));
  return res;
}
