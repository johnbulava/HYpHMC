#include "NeubergerMatrix.h"

void NeubergerMatrix::iniN(double rhoValue, double rValue, int OneDimLSizeL0, int OneDimLSizeL1, int OneDimLSizeL2, int OneDimLSizeL3, int NCopies) {
  if (LogLevel>2) printf("NeubergerMatrix: Creating Neuberger operator with Lattice size %dx%dx%dx%d, %d nested copies, rho=%1.3f and  r=%1.3f...", OneDimLSizeL0, OneDimLSizeL1, OneDimLSizeL2, OneDimLSizeL3, NCopies,rhoValue,rValue);
  r = rValue;
  rho = rhoValue;

  iniD(OneDimLSizeL0, OneDimLSizeL1, OneDimLSizeL2, OneDimLSizeL3, NCopies);
  
  if (LogLevel>2) printf("sucessfully.\n");
}


NeubergerMatrix::NeubergerMatrix(double rhoValue, double rValue, int OneDimLSizeL0, int OneDimLSizeL1, int OneDimLSizeL2, int OneDimLSizeL3, int NCopies) : DiracMatrix(4*NCopies*OneDimLSizeL0*OneDimLSizeL1*OneDimLSizeL2*OneDimLSizeL3) {
  iniN(rhoValue, rValue, OneDimLSizeL0, OneDimLSizeL1, OneDimLSizeL2, OneDimLSizeL3, NCopies);
}


NeubergerMatrix::NeubergerMatrix() : DiracMatrix(4) {
  iniN(1.0, 0.5, 1,1,1,1, 1);
}


NeubergerMatrix::NeubergerMatrix(const NeubergerMatrix& n) : DiracMatrix(n) {
  r = n.r;
  rho = n.rho;
}

double NeubergerMatrix::getRho() {
  return rho;
}

double NeubergerMatrix::getR() {
  return r;
}

NeubergerMatrix& NeubergerMatrix::operator = (NeubergerMatrix n) {
  desini();
  ini(n);
  
  r = n.r;
  rho = n.rho;
  NestedCopies = n.NestedCopies;
  OneDimLatticeSizeL0 = n.OneDimLatticeSizeL0;
  OneDimLatticeSizeL1 = n.OneDimLatticeSizeL1;
  OneDimLatticeSizeL2 = n.OneDimLatticeSizeL2;
  OneDimLatticeSizeL3 = n.OneDimLatticeSizeL3;
  
  return *this;
}


NeubergerMatrix::~NeubergerMatrix() {
}


Complex NeubergerMatrix::analyticalEigenvalue(vector4D p) {
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
  double nenner = sqrt(pTildeSquare + sqr(r*pDoubleTildeSquare - rho));
  result.x = rho + rho*(r*pDoubleTildeSquare - rho)/nenner;
  result.y = rho*sqrt(pTildeSquare)/nenner;
  return result;
}

