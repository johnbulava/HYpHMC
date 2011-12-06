#ifndef HighPrecisionComplex_included
#define HighPrecisionComplex_included


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


class HighPrecisionComplex {
private:
  cln::float_format_t clnDIGIT;  
  cln::cl_F  ONE;
  cln::cl_F  TWO;
  cln::cl_F  ZERO;
  cln::cl_F  HALF;
  
  void ini(int digit, double vx, double vy);
  void ini(int digit, double vx); 
  void ini(const HighPrecisionComplex& c);  
  void ini(int digit, cln::cl_N clnC);
  void ini(int digit);
  void desini();
  
public:
  HighPrecisionComplex(int digit, double vx, double vy); 
  HighPrecisionComplex(int digit, Complex c); 
  HighPrecisionComplex(int digit, cln::cl_N clnC); 
  HighPrecisionComplex(int digit); 
  HighPrecisionComplex(int digit, double vx); 
  
  HighPrecisionComplex(const HighPrecisionComplex& c);  
  ~HighPrecisionComplex();
  
  void setZero();
  void setOne();
  void setTwo();
  void setHalf();
  
  void print();
  
  Complex getComplex();

  int DIGIT; 
  cln::cl_N z;
  
  
  
  HighPrecisionComplex operator + (HighPrecisionComplex);
  HighPrecisionComplex operator - (HighPrecisionComplex);
  HighPrecisionComplex operator * (HighPrecisionComplex);
  HighPrecisionComplex operator / (HighPrecisionComplex);
  
  
  
  HighPrecisionComplex& operator = (HighPrecisionComplex);  
};

#endif
