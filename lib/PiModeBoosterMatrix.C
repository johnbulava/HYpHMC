#include "PiModeBoosterMatrix.h"

void PiModeBoosterMatrix::iniB(double b, int OneDimLSize, int NCopies) {
  if (LogLevel>2) printf("PiModeBoosterMatrix: Creating pi-modes booster matrix with 1D-Lattice size %d, %d nested copies and boost=%1.3f...\n",OneDimLSize,NCopies,b);
  boost = b;

  iniD(OneDimLSize, NCopies);
  
  if (LogLevel>2) printf("sucessfully.\n");
}


PiModeBoosterMatrix::PiModeBoosterMatrix(double b, int OneDimLSize, int NCopies) : DiracMatrix(4*NCopies*OneDimLSize*OneDimLSize*OneDimLSize*OneDimLSize) {
  iniB(b, OneDimLSize, NCopies);
}


PiModeBoosterMatrix::PiModeBoosterMatrix() : DiracMatrix(4) {
  iniB(0.5, 1, 1);
}


PiModeBoosterMatrix::PiModeBoosterMatrix(const PiModeBoosterMatrix& w) : DiracMatrix(w) {
  boost = w.boost;
}


double PiModeBoosterMatrix::getBoost() {
  return boost;
}


PiModeBoosterMatrix& PiModeBoosterMatrix::operator = (PiModeBoosterMatrix w) {
  desini();
  ini(w);
  
  boost = w.boost;
  NestedCopies = w.NestedCopies;
  OneDimLatticeSize = w.OneDimLatticeSize;
  
  return *this;
}


PiModeBoosterMatrix::~PiModeBoosterMatrix() {
}


Complex PiModeBoosterMatrix::analyticalEigenvalue(vector4D p) {
  double delta = (2.0*pi/OneDimLatticeSize) / 1E6;

  int i[4];
  int I;
  for (I=0; I<4; I++) {
    i[I] = (int)(p[I] / pi);
    if (fabs(p[I]-i[I]*pi)>delta) return Complex(1.0, 0.0);
    if (i[I] < 0) i[I] = -i[I];
    i[I] = i[I] % 2;
  }
  
  if ((i[0]==0) && (i[1]==0) && (i[2]==0) && (i[3]==0)) {
    return Complex(1.0, 0.0);
  }
  return Complex(boost, 0.0);
}
