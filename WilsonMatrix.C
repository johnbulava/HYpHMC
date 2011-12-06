#include "WilsonMatrix.h"

void WilsonMatrix::iniW(double rValue, int OneDimLSizeL0, int OneDimLSizeL1, int OneDimLSizeL2, int OneDimLSizeL3, int NCopies) {
  if (LogLevel>2) printf("WilsonMatrix: Creating Wilson operator with Lattice size %dx%dx%dx%d, %d nested copies and r=%1.3f...",OneDimLSizeL0,OneDimLSizeL1,OneDimLSizeL2,OneDimLSizeL3,NCopies,rValue);
  r = rValue;

  iniD(OneDimLSizeL0, OneDimLSizeL1, OneDimLSizeL2, OneDimLSizeL3, NCopies);
  
  if (LogLevel>2) printf("sucessfully.\n");
}


WilsonMatrix::WilsonMatrix(double rValue, int OneDimLSizeL0, int OneDimLSizeL1, int OneDimLSizeL2, int OneDimLSizeL3, int NCopies) : DiracMatrix(4*NCopies*OneDimLSizeL0*OneDimLSizeL1*OneDimLSizeL2*OneDimLSizeL3) {
  iniW(rValue, OneDimLSizeL0, OneDimLSizeL1, OneDimLSizeL2, OneDimLSizeL3, NCopies);
}


WilsonMatrix::WilsonMatrix() : DiracMatrix(4) {
  iniW(0.5, 1,1,1,1, 1);
}


WilsonMatrix::WilsonMatrix(const WilsonMatrix& w) : DiracMatrix(w) {
  r = w.r;
}


double WilsonMatrix::getR() {
  return r;
}


WilsonMatrix& WilsonMatrix::operator = (WilsonMatrix w) {
  desini();
  ini(w);
  
  r = w.r;
  NestedCopies = w.NestedCopies;
  OneDimLatticeSizeL0 = w.OneDimLatticeSizeL0;
  OneDimLatticeSizeL1 = w.OneDimLatticeSizeL1;
  OneDimLatticeSizeL2 = w.OneDimLatticeSizeL2;
  OneDimLatticeSizeL3 = w.OneDimLatticeSizeL3;
  
  return *this;
}


WilsonMatrix::~WilsonMatrix() {
}


Complex WilsonMatrix::analyticalEigenvalue(vector4D p) {
  int I;
  double pTildeSquare = 0;
  double pDoubleTildeSquare = 0;
  double dummy;
  Complex result;

  for (I=0; I<4; I++) {
    dummy = sin(p[I]);
    pTildeSquare += dummy * dummy;
    dummy = 2*sin(0.5*p[I]);
    pDoubleTildeSquare += dummy * dummy;
  }
  result.x = r * pDoubleTildeSquare;
  result.y = sqrt(pTildeSquare);
  return result;
}

