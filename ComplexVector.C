#include "ComplexVector.h"
#include <stdio.h>

void ComplexVector::ini(const ComplexVector& v) {
  ini(v.vectorSize, (Complex*)NULL);
  vectorSize = v.vectorSize;
  int I;
  for (I=0; I<vectorSize; I++) {
    vectorElements[I] = v.vectorElements[I];
  }
}


void ComplexVector::desini() {
  vectorSize = 0;
  if (memSelfDecalred) {
    delete [] vectorElements;
  }
  memSelfDecalred = false;
}
 

ComplexVector::ComplexVector(const ComplexVector& v) {
  ini(v);
}

 
ComplexVector::ComplexVector(int size) {
  ini(size,(Complex*)NULL);
}


ComplexVector::ComplexVector() {
  ini(1,(Complex*)NULL);
}


ComplexVector::~ComplexVector() {
  desini();
}


ComplexVector::ComplexVector(int size, Complex* vecData) {
  ini(size,vecData);
}


double ComplexVector::getNorm() {
  int I;
  double s = 0;
  
  for (I=0; I<vectorSize; I++) {
    s += vectorElements[I].x*vectorElements[I].x + vectorElements[I].y*vectorElements[I].y;
  }
  return sqrt(s);
}


ComplexVector& ComplexVector::operator = (ComplexVector v) {
  desini();
  ini(v);
  return *this;
}


void ComplexVector::ini(int size, Complex* vecData) {
  if (size <= 0) size = 0;
//  printf("ComplexVector: Creating Vector with size %d...",size);
  vectorSize = size;  
  if (vecData == NULL) {
    vectorElements = new Complex[vectorSize];
    memSelfDecalred = true;
    int I;
    for (I=0; I<vectorSize; I++) {
      vectorElements[I].x = 0;
      vectorElements[I].y = 0;    
    }
  } else {
    vectorElements = vecData;
    memSelfDecalred = false;
  }
//  printf("sucessfully.\n");
}


void ComplexVector::resize(int size) {
  desini();
  ini(size);
}


void ComplexVector::setZero() {
  int I;
  for (I=0; I<vectorSize; I++) {
    vectorElements[I].x = 0.0;
    vectorElements[I].y = 0.0;
  }
}


void ComplexVector::print() {
  int I;
  
  for (I=0; I<vectorSize; I++) {
    double r = vectorElements[I].x;
    double i = vectorElements[I].y;
    if (r>=0) printf(" ");
    printf("%1.3f",r);
    if (i>=0) printf("+");
    printf("%1.3fi\n",i);
  }
}


ComplexVector ComplexVector::operator + (ComplexVector param) {
  if (param.vectorSize!=vectorSize) {
    printf("ComplexVector::operator +: ERROR: Vector sizes do not match!\n");
    return ComplexVector(0);
  }
  int I;
  ComplexVector result(vectorSize);
  for (I=0; I<vectorSize; I++) {
    result.vectorElements[I] = vectorElements[I] + param.vectorElements[I];
  }
  return result;
}


ComplexVector ComplexVector::operator - (ComplexVector param) {
  if (param.vectorSize!=vectorSize) {
    printf("ComplexVector::operator -: ERROR: Vector sizes do not match!\n");
    return ComplexVector(0);
  }
  int I;
  ComplexVector result(vectorSize);
  for (I=0; I<vectorSize; I++) {
    result.vectorElements[I] = vectorElements[I] - param.vectorElements[I];
  }
  return result;
}



ComplexVector operator * (double A, ComplexVector B) {
  ComplexVector res(B.vectorSize);
  int I;
  
  for (I=0; I<B.vectorSize; I++) {
    res.vectorElements[I] = A*B.vectorElements[I];  
  }

  return res;
}


ComplexVector operator * (Complex A, ComplexVector B) {
  ComplexVector res(B.vectorSize);
  int I;
  
  for (I=0; I<B.vectorSize; I++) {
    res.vectorElements[I] = A*B.vectorElements[I];  
  }

  return res;
}
