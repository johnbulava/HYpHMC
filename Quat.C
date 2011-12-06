#include "Quat.h"
#include <stdio.h>

Quat::Quat() { }

Quat::Quat(double valx0, double valx1, double valx2, double valx3) {
  x0 = valx0;
  x1 = valx1;
  x2 = valx2;
  x3 = valx3;
}

Quat::~Quat() { }


inline void Quat::print() {
  printf("(%1.7f, %1.7f, %1.7f, %1.7f)\n",x0,x1,x2,x3);
}


inline void Quat::setValues(double valx0, double valx1, double valx2, double valx3) {
  x0 = valx0;
  x1 = valx1;
  x2 = valx2;
  x3 = valx3;
}


  
inline Quat Quat::operator + (Quat B) {
  Quat res;

  res.x0 = x0 + B.x0;
  res.x1 = x1 + B.x1;
  res.x2 = x2 + B.x2;
  res.x3 = x3 + B.x3;

  return res;
}


inline Quat Quat::operator - (Quat B) {
  Quat res;

  res.x0 = x0 - B.x0;
  res.x1 = x1 - B.x1;
  res.x2 = x2 - B.x2;
  res.x3 = x3 - B.x3;

  return res;
}


inline Quat Quat::operator * (Quat B) {
  Quat res;

  res.x0 = x0*B.x0 - x1*B.x1 - x2*B.x2 - x3*B.x3;
  res.x1 = x0*B.x1 + x1*B.x0 + x2*B.x3 - x3*B.x2;
  res.x2 = x0*B.x2 + x2*B.x0 + x3*B.x1 - x1*B.x3;
  res.x3 = x0*B.x3 + x3*B.x0 + x1*B.x2 - x2*B.x1;

  return res;
}


inline Quat Quat::operator / (Quat B) {
  Quat res;
  double det = B.x0*B.x0 + B.x1*B.x1 + B.x2*B.x2 + B.x3*B.x3;

  res.x0 = (x0*B.x0 + x1*B.x1 + x2*B.x2 + x3*B.x3) / det;
  res.x1 = (x3*B.x2 + x1*B.x0 - x0*B.x1 - x2*B.x3) / det;
  res.x2 = (x2*B.x0 + x1*B.x3 - x0*B.x2 - x3*B.x1) / det;
  res.x3 = (x3*B.x0 + x2*B.x1 - x0*B.x3 - x1*B.x2) / det;

  return res;
}


inline double Quat::operator | (Quat B) {
  double res;

  res = x0*B.x0 - x1*B.x1 - x2*B.x2 - x3*B.x3;

  return res;
}


Quat adj(Quat A) {
  Quat res;

  res.x0 = A.x0;
  res.x1 = -A.x1;
  res.x2 = -A.x2;
  res.x3 = -A.x3;

  return res;
}


inline Quat inv (Quat A) {
  Quat res;
  double det = A.x0*A.x0 + A.x1*A.x1 + A.x2*A.x2 + A.x3*A.x3;

  res.x0 = A.x0  / det;
  res.x1 = -A.x1 / det;
  res.x2 = -A.x2 / det;
  res.x3 = -A.x3 / det;

  return res;
}


inline Quat Quat::norm() {
  Quat res;
  double n = sqrt(x0*x0 + x1*x1 + x2*x2 + x3*x3);

  res.x0 = x0 / n;
  res.x1 = x1 / n;
  res.x2 = x2 / n;
  res.x3 = x3 / n;

  return res;
}


double Quat::getNorm () {
//inline double Quat::getNorm () {
  return sqrt(x0*x0 + x1*x1 + x2*x2 + x3*x3);
}


inline double Quat::getNormSquare() {
  return x0*x0 + x1*x1 + x2*x2 + x3*x3;
}



