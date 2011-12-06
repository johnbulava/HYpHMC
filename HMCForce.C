#include "HMCForce.h"

HMCForce::HMCForce(FermionMatrixOperations* fOps, bool loc): Force(fOps, loc, 0) {
  if (LogLevel>2) printf("Initializing HMC-Force-Calculator with LocalMode = %d\n", local);
  iniAdditionalFields();
}


void HMCForce::iniAdditionalFields() {
  if (local) {
    outMMdaggerInverseOmega = fermiOps->createFermionVector();
  } else {
    outMMdaggerInverseOmega = NULL;
  }
  omegaMMdaggerInverseOmegaScalarProduct = NaN;  
}


void HMCForce::desiniAdditionalFields() {
  fermiOps->destroyFermionVector(outMMdaggerInverseOmega);
}


HMCForce::~HMCForce()  {
  if (LogLevel>0) printf("Desinitializing HMC-Force-Calculator...");
  desini();
  
  if (LogLevel>0) printf("sucessfully.\n");      
}


void HMCForce::sampleOmegaField(vector4D* phi, Complex* GaussVector) {
  if (LogLevel>3) printf("HMC-Force-Calculator: Sampling omega field from Gauss-Vector.\n");
  //Nutze uebergebenes Phi-Feld!!!
  fermiOps->transformToXtraSizeArray(GaussVector, omegaField);
  fermiOps->executeFermionMatrixMultiplication(omegaField, omegaField, (double*) phi, false, NULL, NULL, 1, 0); 
}
 

/**
* ACHTUNG: Kraft ist mit Faktor 2 multipliziert!!!
* D.h Dies ist die Ableitung von /omega (MMdag)^-1 \omega (ohne den Faktor 1/2)
**/
void HMCForce::calcPhiForce(vector4D* phi, double TOL) {
  fermiOps->executeMultiplicationVectorWithDerivativesOfMMdaggerInverse(omegaField, outMMdaggerInverseOmega, (double*) dSdPhi, (double*) phi, TOL);
}


void HMCForce::calcInverse(vector4D* phi, double TOL) {
  int neededIter = 0;
  fermiOps->solveFermionMatrixLGS(omegaField, outMMdaggerInverseOmega, (double*) phi, TOL, true, false, -1, neededIter);  
  double actualAcc = fermiOps->checkLGSsolutionAccuracy(outMMdaggerInverseOmega, omegaField, (double*) phi);  
  if (LogLevel>2) {
    printf("Matrix-Inversion: %d iterations needed. Actual absolute accuracy %1.15f (requested %1.15f)\n", neededIter, actualAcc, TOL);
  }
}


void HMCForce::calcOmegaMMdaggerInverseOmegaScalarProduct() {
//  int VL = fermiOps->getVectorLength();
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  Complex dummy;  

//  cblas_zdotc_sub(VL, omegaField, 1, outMMdaggerInverseOmega, 1, &dummy);
  SSE_ComplexScalarProduct(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, omegaField, outMMdaggerInverseOmega, dummy);
  omegaMMdaggerInverseOmegaScalarProduct = dummy.x;
}


Complex* HMCForce::getMMdaggerInverseOmega() {
  return outMMdaggerInverseOmega;
}


double HMCForce::getOmegaMMdaggerInverseOmegaScalarProduct() {
  return omegaMMdaggerInverseOmegaScalarProduct;
}


void HMCForce::setOmegaMMdaggerInverseOmegaScalarProduct(double s) {
  omegaMMdaggerInverseOmegaScalarProduct = s;
}
