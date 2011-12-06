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



Quat adj(Quat A) {
  Quat res;

  res.x0 = A.x0;
  res.x1 = -A.x1;
  res.x2 = -A.x2;
  res.x3 = -A.x3;

  return res;
}


