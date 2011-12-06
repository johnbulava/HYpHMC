#ifndef ComplexVector_included
#define ComplexVector_included

#include "Complex.h"
#include "Global.h"
#ifdef useBLAS
  #include CBLASIncludeFile
#endif

class ComplexVector {
protected:
  void ini(int size, Complex* vecData);
  void ini(const ComplexVector& v);
  void desini();
  bool memSelfDecalred;
 
public:
  ComplexVector();
  ComplexVector(const ComplexVector& v);  
  ComplexVector(int size);
  ComplexVector(int size, Complex* vecData);
  ~ComplexVector();
  Complex* vectorElements;
  void resize(int size);
  int vectorSize;
  void print();
  double getNorm();
  void setZero();
  
  ComplexVector operator + (ComplexVector);
  ComplexVector operator - (ComplexVector); 
  ComplexVector& operator = (ComplexVector);  
   
};

#endif
