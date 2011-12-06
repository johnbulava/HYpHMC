#ifndef Quat_included
#define Quat_included

#include <math.h>

class Quat {
private:


public:
  Quat(); 
  Quat(double valx0, double valx1, double valx2, double valx3);
  ~Quat();
  
  double x0,x1,x2,x3;
  void setValues(double valx0, double valx1, double valx2, double valx3);

  void print();

  Quat operator + (Quat);
  Quat operator - (Quat);  
  Quat operator * (Quat);
  Quat operator / (Quat);
  double operator | (Quat);
  Quat norm();
  double getNorm ();
  double getNormSquare();
};

inline Quat operator * (double A, Quat B) {
  Quat res;

  res.x0 = A*B.x0;
  res.x1 = A*B.x1;
  res.x2 = A*B.x2;
  res.x3 = A*B.x3;

  return res;
}


inline Quat exp(Quat A) {
  double w = sqrt( A.x1*A.x1 + A.x2*A.x2 + A.x3*A.x3 );
  double sinw = sin(w);
  double cosw = cos(w);
  double ex = exp(A.x0);
  double dummy;
  Quat res;

  res.x0 = ex * cosw;
  if (w > 0) {
    dummy = ex * sinw / w;
    res.x1 = dummy * A.x1;
    res.x2 = dummy * A.x2;
    res.x3 = dummy * A.x3;
    return res;
  } else {
    res.x1 = 0;
    res.x2 = 0;
    res.x3 = 0;
  }

  return res;
}

#endif
