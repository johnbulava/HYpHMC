#include "pHMCForce.h"

pHMCForce::pHMCForce(FermionMatrixOperations* fOps, bool loc, double tht, int id, int addAuxVecsCount): Force(fOps, loc, id) {
  if (LogLevel>2) printf("Initializing pHMC-Force-Calculator %d with LocalMode = %d, theta = %1.3f, and addAuxVecsCount = %d\n", id, local, tht, addAuxVecsCount);
  theta = tht;
  nu = new double[pHMCForcePolynomialsMAX];
  int I;
  for (I=0; I<pHMCForcePolynomialsMAX; I++) nu[I] = 1.0;

  MomentaMasses = NULL;
  AdditionalAuxVectors = NULL;
  omegaForceStrengthAnalysisPREC = NULL;
  omegaForceStrengthAnalysisGLOBAL = NULL;
  avgOmegaForces = NULL;
  
  AdditionalAuxVectorsCount = addAuxVecsCount;
  VectorCollectionCount = 0;
  VectorCollection = NULL;
  dSdOmega = NULL;
  storedActForcesAvail = NULL;
  storedTraStartForcesAvail = NULL;
  storedActOmegaAction = NULL;
  storedTraStartOmegaAction = NULL;
  storedActdSdOmega = NULL;
  storedTraStartdSdOmega = NULL;
  storedActdSdPhi = NULL;
  storedTraStartdSdPhi = NULL;
  chebyApproxOfInverseSQRTofPolynomial = NULL;
  fastEvalVectors = NULL;
  fastEvalVectorsCount = 0; 

  Distributed_Phi = NULL;
  Distributed_dSdPhiDummy = NULL;
  Distributed_dSdPhiSum = NULL;
  Distributed_omegaField = NULL;
  Distributed_auxVectors = NULL;
  Distributed_fastEvalVectors = NULL;
  Distributed_AdditionalAuxVectors = NULL;
  Distributed_VectorCollection = NULL;

  quasiHermiteanMode = true;
  bMatFactorizationMode = false;

  iniAdditionalFields();
}


void pHMCForce::iniAdditionalFields() {
  if (local) {
    setOmegaFieldToZero();
    omegaMomenta = fermiOps->createFermionVector();
    omegaMomentaOLD = fermiOps->createFermionVector();
    omegaOLD = fermiOps->createFermionVector();
    auxVectors = new Complex*[4];
    auxVectors[0] = fermiOps->createFermionVector();
    auxVectors[1] = fermiOps->createFermionVector();
    auxVectors[2] = fermiOps->createFermionVector();
    auxVectors[3] = fermiOps->createFermionVector();
    
    if (!(fermiOps->isMultiThreadedOpsActivated())) {
      AdditionalAuxVectors = new Complex*[AdditionalAuxVectorsCount];
      for (int I=0; I<AdditionalAuxVectorsCount; I++) {
        AdditionalAuxVectors[I] = fermiOps->createFermionVector();
      }  
    } else {
      AdditionalAuxVectors = NULL;    
    }

    if (fermiOps->isMultiThreadedOpsActivated()) {
      Distributed_Phi = fermiOps->createDistributedUniformLatticeComplexVector(2); 
      Distributed_dSdPhiDummy = fermiOps->createDistributedUniformLatticeComplexVector(4);
      Distributed_dSdPhiSum = fermiOps->createDistributedUniformLatticeComplexVector(4);      
      Distributed_omegaField = fermiOps->createDistributedFermionVector();
      Distributed_auxVectors = new DistributedMemoryObject*[4];
      Distributed_auxVectors[0] = fermiOps->createDistributedFermionVector();
      Distributed_auxVectors[1] = fermiOps->createDistributedFermionVector();
      Distributed_auxVectors[2] = fermiOps->createDistributedFermionVector();
      Distributed_auxVectors[3] = fermiOps->createDistributedFermionVector();

      Distributed_AdditionalAuxVectors = new DistributedMemoryObject*[AdditionalAuxVectorsCount]; 
      for (int I=0; I<AdditionalAuxVectorsCount; I++) {
        Distributed_AdditionalAuxVectors[I] = fermiOps->createDistributedFermionVector();
      }  
    } else {
      Distributed_AdditionalAuxVectors = NULL;    
    }

    setVectorCollection();
  } else {
    omegaMomenta = NULL;
    omegaMomentaOLD = NULL;
    omegaOLD = NULL;
    auxVectors = NULL;
    AdditionalAuxVectors = NULL;
    VectorCollectionCount = 0;
    VectorCollection = NULL;
    
    Distributed_Phi = NULL;
    Distributed_dSdPhiDummy = NULL;
    Distributed_dSdPhiSum = NULL;      
    Distributed_omegaField = NULL;
    Distributed_auxVectors = NULL;
    Distributed_fastEvalVectors = NULL;
    Distributed_AdditionalAuxVectors = NULL;
    Distributed_VectorCollection = NULL;
  }   
  
  storedActForcesAvail = new bool[pHMCForcePolynomialsMAX];
  storedTraStartForcesAvail = new bool[pHMCForcePolynomialsMAX];
  storedActOmegaAction = new double[pHMCForcePolynomialsMAX];
  storedTraStartOmegaAction = new double[pHMCForcePolynomialsMAX];
  storedActdSdOmega = new Complex*[pHMCForcePolynomialsMAX];
  storedTraStartdSdOmega = new Complex*[pHMCForcePolynomialsMAX];
  storedActdSdPhi = new Complex*[pHMCForcePolynomialsMAX];
  storedTraStartdSdPhi = new Complex*[pHMCForcePolynomialsMAX];
  int I;
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    storedActForcesAvail[I] = false;
    storedTraStartForcesAvail[I] = false;
    storedActOmegaAction[I] = NaN;
    storedTraStartOmegaAction[I] = NaN;
    storedActdSdOmega[I] = NULL;
    storedTraStartdSdOmega[I] = NULL;
    storedActdSdPhi[I] = NULL;
    storedTraStartdSdPhi[I] = NULL;
  }

  omegaMomentumAction = NaN;
  exactOmegaMMdagInverseSQRTAction = NaN;
  approxpolyRoots = new Complex*[pHMCForcePolynomialsMAX];
  approxpolyLambda = new double[pHMCForcePolynomialsMAX];
  approxpolyDegree = new int[pHMCForcePolynomialsMAX];
  chebyApproxOfInverseSQRTofPolynomial = new GeneralChebyshevApproximation*[pHMCForcePolynomialsMAX];
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    approxpolyRoots[I]  = NULL;
    approxpolyDegree[I] = 0;    
    approxpolyLambda[I] = 0;
    chebyApproxOfInverseSQRTofPolynomial[I] = NULL;
  }

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  omegaForceStrengthAnalysisPREC = new LatticeMomentumBins*[pHMCForcePolynomialsMAX];
  omegaForceStrengthAnalysisGLOBAL = new LatticeMomentumBins*[pHMCForcePolynomialsMAX];
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    omegaForceStrengthAnalysisPREC[I] = new LatticeMomentumBins(L0,L1,L2,L3);
    omegaForceStrengthAnalysisGLOBAL[I] = new LatticeMomentumBins(L0,L1,L2,L3);
  }
  MomentaMasses = new double[L0*L1*L2*L3];
  avgOmegaForces = new double[L0*L1*L2*L3];
  for (I=0; I<L0*L1*L2*L3; I++) MomentaMasses[I] = 1.0;
  omegaMassAdaptionMode = false;

  omegaForceStrengthSumPREC = new double[pHMCForcePolynomialsMAX];
  omegaForceStrengthSqrSumPREC = new double[pHMCForcePolynomialsMAX];
  omegaForceStrengthCountPREC = new int[pHMCForcePolynomialsMAX];
  omegaForceStrengthSumGLOBAL = new double[pHMCForcePolynomialsMAX];
  omegaForceStrengthSqrSumGLOBAL = new double[pHMCForcePolynomialsMAX];
  omegaForceStrengthCountGLOBAL = new int[pHMCForcePolynomialsMAX];
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    omegaForceStrengthSumPREC[I] = 0;
    omegaForceStrengthSqrSumPREC[I] = 0;
    omegaForceStrengthCountPREC[I] = 0;
    omegaForceStrengthSumGLOBAL[I] = 0;
    omegaForceStrengthSqrSumGLOBAL[I] = 0;
    omegaForceStrengthCountGLOBAL[I] = 0;
  }
}


void pHMCForce::desiniAdditionalFields() {
  desiniFastEvalVectors();
  if (local) {
    fermiOps->destroyFermionVector(auxVectors[0]);
    fermiOps->destroyFermionVector(auxVectors[1]);  
    fermiOps->destroyFermionVector(auxVectors[2]);  
    fermiOps->destroyFermionVector(auxVectors[3]);  
    delete[] auxVectors;
    fermiOps->destroyFermionVector(omegaMomenta); 
    fermiOps->destroyFermionVector(omegaMomentaOLD);     
    fermiOps->destroyFermionVector(omegaOLD);
    
    if (!(fermiOps->isMultiThreadedOpsActivated())) {
      for (int I=0; I<AdditionalAuxVectorsCount; I++) {
        fermiOps->destroyFermionVector(AdditionalAuxVectors[I]);
      }
      delete[] AdditionalAuxVectors;
    } 
    
    if (fermiOps->isMultiThreadedOpsActivated()) {
      fermiOps->destroyDistributedComplexVector(Distributed_Phi); 
      fermiOps->destroyDistributedComplexVector(Distributed_dSdPhiDummy); 
      fermiOps->destroyDistributedComplexVector(Distributed_dSdPhiSum); 
      fermiOps->destroyDistributedComplexVector(Distributed_omegaField); 
      fermiOps->destroyDistributedComplexVector(Distributed_auxVectors[0]); 
      fermiOps->destroyDistributedComplexVector(Distributed_auxVectors[1]); 
      fermiOps->destroyDistributedComplexVector(Distributed_auxVectors[2]); 
      fermiOps->destroyDistributedComplexVector(Distributed_auxVectors[3]); 
      delete[] Distributed_auxVectors;
      fermiOps->destroyDistributedComplexVector(Distributed_omegaField); 
      for (int I=0; I<AdditionalAuxVectorsCount; I++) {
        fermiOps->destroyDistributedComplexVector(Distributed_AdditionalAuxVectors[I]); 
      }  
      delete[] Distributed_AdditionalAuxVectors;
    } 
    
    VectorCollectionCount = 0;
    if (VectorCollection!=NULL) delete[] VectorCollection;
    if (Distributed_VectorCollection!= NULL) delete[] Distributed_VectorCollection;
    
    VectorCollection = NULL;
    Distributed_VectorCollection = NULL;
  }  
  activateForceStoring(0,NULL);
  delete[] storedActForcesAvail;
  delete[] storedTraStartForcesAvail;
  delete[] storedActOmegaAction;
  delete[] storedTraStartOmegaAction;
  delete[] storedActdSdOmega;
  delete[] storedTraStartdSdOmega;
  delete[] storedActdSdPhi;
  delete[] storedTraStartdSdPhi;
  int I;
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    delete omegaForceStrengthAnalysisPREC[I];
    delete omegaForceStrengthAnalysisGLOBAL[I];
  }  
  delete[] omegaForceStrengthAnalysisPREC;
  delete[] omegaForceStrengthAnalysisGLOBAL;
  delete[] MomentaMasses;
  delete[] avgOmegaForces;
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    delete[] approxpolyRoots[I];
    delete chebyApproxOfInverseSQRTofPolynomial[I];
  }
  delete[] chebyApproxOfInverseSQRTofPolynomial;  
  delete[] approxpolyRoots;
  delete[] approxpolyDegree;
  delete[] approxpolyLambda;
  delete[] nu;
  delete[] omegaForceStrengthSumPREC;
  delete[] omegaForceStrengthSqrSumPREC;
  delete[] omegaForceStrengthCountPREC;
  delete[] omegaForceStrengthSumGLOBAL;
  delete[] omegaForceStrengthSqrSumGLOBAL;
  delete[] omegaForceStrengthCountGLOBAL;
}


pHMCForce::~pHMCForce()  {
  if (LogLevel>0) printf("Desinitializing HMC-Force-Calculator...");
  desini();
  
  if (LogLevel>0) printf("sucessfully.\n");      
}


void pHMCForce::setVectorCollection() {
  if (local) {
    VectorCollectionCount = 0;
    if (VectorCollection!=NULL) delete[] VectorCollection;
    if (Distributed_VectorCollection!= NULL) delete[] Distributed_VectorCollection;
    VectorCollection = NULL;
    Distributed_VectorCollection = NULL;
  
    VectorCollectionCount = AdditionalAuxVectorsCount + fastEvalVectorsCount;
    if (!(fermiOps->isMultiThreadedOpsActivated())) {
      VectorCollection = new Complex*[VectorCollectionCount];
      for (int I=0; I<AdditionalAuxVectorsCount; I++) {
        VectorCollection[I] = AdditionalAuxVectors[I];
      }
      for (int I=0; I<fastEvalVectorsCount; I++) {
        VectorCollection[AdditionalAuxVectorsCount+I] = fastEvalVectors[I];
      }
    }
    if (fermiOps->isMultiThreadedOpsActivated()) {
      Distributed_VectorCollection = new DistributedMemoryObject*[VectorCollectionCount];
      
      for (int I=0; I<AdditionalAuxVectorsCount; I++) {
        Distributed_VectorCollection[I] = Distributed_AdditionalAuxVectors[I];
      }
      for (int I=0; I<fastEvalVectorsCount; I++) {
        Distributed_VectorCollection[AdditionalAuxVectorsCount+I] = Distributed_fastEvalVectors[I];
      }
    }
    if (LogLevel>2) printf("VectorCollection set on force %d on node %d: %d entries.\n",ID,ownNodeID, VectorCollectionCount);
  } else {
    VectorCollectionCount = 0;
    if (VectorCollection!=NULL) delete[] VectorCollection;
    if (Distributed_VectorCollection!= NULL) delete[] Distributed_VectorCollection;
    
    VectorCollection = NULL; 
    Distributed_VectorCollection = NULL;
  }
}


void pHMCForce::iniFastEvalVectors(int degMax) {
  desiniFastEvalVectors();
  if (local) {
    if (LogLevel>3) printf("initializing fastEval vectors on Force %d on node %d...\n",ID, ownNodeID);
    fastEvalVectorsCount = degMax+1;
    if (!(fermiOps->isMultiThreadedOpsActivated())) {
      fastEvalVectors = new Complex*[fastEvalVectorsCount];
      for (int I=0; I<fastEvalVectorsCount; I++) {
        fastEvalVectors[I] = fermiOps->createFermionVector();
      }
    }
    if (fermiOps->isMultiThreadedOpsActivated()) {
      Distributed_fastEvalVectors = new DistributedMemoryObject*[fastEvalVectorsCount];
      for (int I=0; I<fastEvalVectorsCount; I++) {
        Distributed_fastEvalVectors[I] = fermiOps->createDistributedFermionVector();
      }
    }
    setVectorCollection();
    
    if (LogLevel>3) printf("Initialization of fastEval vectors on Force %d on node %d successful.\n",ID, ownNodeID);    
  }
}


void pHMCForce::desiniFastEvalVectors() {
  if (local) {
    if (!(fermiOps->isMultiThreadedOpsActivated())) {
      for (int I=0; I<fastEvalVectorsCount; I++) {
        fermiOps->destroyFermionVector(fastEvalVectors[I]);
      }  
      delete[] fastEvalVectors;
    }
    if (fermiOps->isMultiThreadedOpsActivated()) {
      for (int I=0; I<fastEvalVectorsCount; I++) {
        fermiOps->destroyDistributedComplexVector(Distributed_fastEvalVectors[I]);
      }
      delete[] Distributed_fastEvalVectors;   
    }
  }
  fastEvalVectorsCount = 0;
  fastEvalVectors = NULL;
  Distributed_fastEvalVectors = NULL;
  setVectorCollection();  
}


/**
*  Expects the roots to come in complex conjugate pairs, following each other
*  directly in the array.
**/
void pHMCForce::setApproxPolyRoots(int polynomSlot, Complex* roots, int deg, double n, double polLam) {
  Complex* dummy = approxpolyRoots[polynomSlot];
  approxpolyRoots[polynomSlot] = new Complex[deg];
  int I;
  if (LogLevel>2) printf("Setting sub-polynomial %d roots for Force-Calculator %d on node %d (degree=%d, polLam=%f):\n",polynomSlot,ID, ownNodeID,deg,polLam);
  double n1 = 0;
  for (I=0; I<deg/2; I++) {
    n1 += (roots[2*I].x - roots[2*I+1].x) *(roots[2*I].x - roots[2*I+1].x);
    n1 += (roots[2*I].y + roots[2*I+1].y) *(roots[2*I].y + roots[2*I+1].y);
  }
  double n2 = 0;
  for (I=0; I<deg/2; I++) {
    n2 += (roots[I].x - roots[deg-I-1].x) *(roots[I].x - roots[deg-I-1].x);
    n2 += (roots[I].y + roots[deg-I-1].y) *(roots[I].y + roots[deg-I-1].y);
  }
  
  if (n1<n2) {
    if (LogLevel>3) printf("Need to change order of monomials!\n");    
    for (I=0; I<deg; I++) {
      if (I<deg/2) {
        approxpolyRoots[polynomSlot][I] = roots[2*I];
      } else {
        approxpolyRoots[polynomSlot][I] = roots[2*deg-2*I-1];
      }
    }
  } else {
    if (LogLevel>3) printf("Keep order of monomials!\n");    
    for (I=0; I<deg; I++) {
      approxpolyRoots[polynomSlot][I] = roots[I];
    }
  }
  
  if (LogLevel>3) {
    for (I=0; I<deg; I++) {
      printf(" -> ");
      approxpolyRoots[polynomSlot][I].print();
    }
  }
  
  approxpolyDegree[polynomSlot] = deg;
  approxpolyLambda[polynomSlot] = polLam;
  if (polynomSlot==0) {
    iniFastEvalVectors(approxpolyDegree[0]);
  }
  delete[] dummy;
  setOmegaFieldToZero();
  setNu(polynomSlot, n);
  setChebyApproxOfInverseSQRTofPolynomial(polynomSlot);
}


Complex* pHMCForce::getApproxPolyRoots(int polynomSlot) {
  return approxpolyRoots[polynomSlot];
}


int pHMCForce::getApproxPolyDegree(int polynomSlot) {
  return approxpolyDegree[polynomSlot];
}


double pHMCForce::getApproxPolyLambda(int polynomSlot) {
  return approxpolyLambda[polynomSlot];
}


GeneralChebyshevApproximation* pHMCForce::getApproxChebyApproxOfInverseSQRTofPolynomial(int polynomSlot) {
  return chebyApproxOfInverseSQRTofPolynomial[polynomSlot];
}


double pHMCForce::getActOmegaAction(int polySlot) {
  return storedActOmegaAction[polySlot];
}


void pHMCForce::setActOmegaAction(int polySlot, double s) {
  storedActOmegaAction[polySlot] = s;
}


double pHMCForce::getExactOmegaMMdagInverseSQRTAction() {
  return exactOmegaMMdagInverseSQRTAction;
}


void pHMCForce::setExactOmegaMMdagInverseSQRTAction(double s) {
  exactOmegaMMdagInverseSQRTAction = s;
}


//Attention: input and output must be different, if only one monomial is applied!!!
Complex* pHMCForce::calcMonomialApplication(int polynomSlot, Complex* input, Complex* output, double* phi, int startM, int endM, bool inFourierSpace) {
  if (endM<startM) {
    int VLxtrSize = fermiOps->getVectorLengthXtrSize();
    if (output != NULL) {
      if (input != output) {
        SSE_ZCopy(VLxtrSize, input, 1, output, 1);  
      }
      return output;
    } else {
      SSE_ZCopy(VLxtrSize, input, 1, auxVectors[0], 1);  
      return auxVectors[0];
    }
  }

  int act = startM;
  int oneDimSizeL0 = fermiOps->get1DSizeL0();
  int oneDimSizeL1 = fermiOps->get1DSizeL1();
  int oneDimSizeL2 = fermiOps->get1DSizeL2();
  int oneDimSizeL3 = fermiOps->get1DSizeL3();
  int auxInd = 1;
  Complex* p1 = input;
  Complex* p2 = auxVectors[auxInd];
  while (act<endM) {
    if (quasiHermiteanMode) {    
      fermiOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(p1, p2, phi, true, inFourierSpace);
    } else {
      fermiOps->executeFermionMatrixFermionDaggerMatrixMultiplication(p1, p2, phi, 2, 1, inFourierSpace);
    }
      
    Complex alpha = (-1)*approxpolyRoots[polynomSlot][act];
    SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, alpha, p1, p2);
    act++;
    p1 = p2;
    auxInd = 1-auxInd;
    p2 = auxVectors[auxInd];
  }

  if (output != NULL) p2 = output;
  if (quasiHermiteanMode) {
    fermiOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(p1, p2, phi, true, inFourierSpace);
  } else {
    fermiOps->executeFermionMatrixFermionDaggerMatrixMultiplication(p1, p2, phi, 2, 1, inFourierSpace);
  }
      
  Complex alpha = (-1)*approxpolyRoots[polynomSlot][act];
  SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, alpha, p1, p2);

  return p2;
}


//Attention: input and output must be different, if only one monomial is applied!!!
DistributedMemoryObject* pHMCForce::calcDistributedMonomialApplication(int polynomSlot, DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* phi, int startM, int endM, bool inFourierSpace) {
  MultiThreadedOperations* threadedOps = fermiOps->getMultiThreadedOps();
  int oneDimSizeL0 = fermiOps->get1DSizeL0();
  int oneDimSizeL1 = fermiOps->get1DSizeL1();
  int oneDimSizeL2 = fermiOps->get1DSizeL2();
  int oneDimSizeL3 = fermiOps->get1DSizeL3();
  if (endM<startM) {
    if (output != NULL) {
      if (input != output) {
        threadedOps->copyFermionVector(input, output, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
      }
      return output;
    } else {
      threadedOps->copyFermionVector(input, Distributed_auxVectors[0], oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
      return Distributed_auxVectors[0];
    }
  }

  int act = startM;
  int auxInd = 1;
  DistributedMemoryObject* p1 = input;
  DistributedMemoryObject* p2 = Distributed_auxVectors[auxInd];
  
  while (act<endM) {
    if (quasiHermiteanMode) {
      fermiOps->executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(p1, p2, phi, true, inFourierSpace);
    } else {
      printf("ERROR in pHMCForce::calcDistributedMonomialApplication: Distributed version of executeFermionMatrixFermionDaggerMatrixMultiplication not implemented yet.\n");
      exit(0);
    }
      
    Complex alpha = (-1)*approxpolyRoots[polynomSlot][act];
    threadedOps->vectorAdditionOfFermionVectors(p1, p2, alpha, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);
    act++;
    p1 = p2;
    auxInd = 1-auxInd;
    p2 = Distributed_auxVectors[auxInd];
  }

  if (output != NULL) p2 = output;
  if (quasiHermiteanMode) {
    fermiOps->executeDistributedFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(p1, p2, phi, true, inFourierSpace);
  } else {
    printf("ERROR in pHMCForce::calcDistributedMonomialApplication: Distributed version of executeFermionMatrixFermionDaggerMatrixMultiplication not implemented yet.\n");
    exit(0);
  }
      
  Complex alpha = (-1)*approxpolyRoots[polynomSlot][act];
  threadedOps->vectorAdditionOfFermionVectors(p1, p2, alpha, oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3);

  return p2;
}


double pHMCForce::calcOmegaAction(int polynomSlot, double* phi) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  if (fermiOps->getYN() == 0) {
    setActOmegaAction(polynomSlot, 0.0);  
    addPerformanceProfilingItem("pHMCForce::calcOmegaAction", performanceProfilerStartCycle, 0);    
    return 0;
  }

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double dummy;
  double dummy2 = 0;

  if (fermiOps->isMultiThreadedOpsActivated()) {
    fermiOps->copyFermionVectorToDistributedFermionVector(omegaField, Distributed_omegaField);
    fermiOps->copyPhiFieldUniformlyToDistributedPhiField(phi, Distributed_Phi);
    DistributedMemoryObject* p = calcDistributedMonomialApplication(polynomSlot, Distributed_omegaField, NULL, Distributed_Phi, 0, (approxpolyDegree[polynomSlot]/2)-1, true);

    MultiThreadedOperations* threadedOps = fermiOps->getMultiThreadedOps();
    threadedOps->vectorNormOfFermionVector(p, dummy, L0, L1, L2, L3);
    if ((polynomSlot%2)==1) {
      threadedOps->vectorNormOfFermionVector(Distributed_omegaField, dummy2, L0, L1, L2, L3);
    }
  } else {
    Complex* p = calcMonomialApplication(polynomSlot, omegaField, NULL, phi, 0, (approxpolyDegree[polynomSlot]/2)-1, true);
    SSE_ComplexSquareNorm(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, p, dummy);
    if ((polynomSlot%2)==1) {
      SSE_ComplexSquareNorm(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, omegaField, dummy2);  
    }
  }

  setActOmegaAction(polynomSlot, 0.5*nu[polynomSlot]*dummy - 0.5*dummy2);
  addPerformanceProfilingItem("pHMCForce::calcOmegaAction", performanceProfilerStartCycle, 0);    
  return getActOmegaAction(polynomSlot);
}


/**
* ATTENTION: This force is the force without the factor 1/2.
*
**/
Complex* pHMCForce::calcAllForces(int polynomSlot, double* phi, double& S) {
  //Is Result already in Storage?
  if (storedActForcesAvail[polynomSlot]) {
    if (LogLevel>4) printf("Taking Forces from Storage for polynomial %d on force %d on node %d.\n", polynomSlot, ID, ownNodeID);
    S = storedActOmegaAction[polynomSlot];
    dSdOmega = storedActdSdOmega[polynomSlot];
    int VLSize = fermiOps->getVectorLength();
    SSE_ZCopy(VLSize/4, storedActdSdPhi[polynomSlot], 1, (Complex*)dSdPhi, 1);  
    return dSdOmega;
  }

  //Set Force to zero...
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  
  for (int I=0; I<L0*L1*L2*L3; I++) {
    dSdPhi[I][0] = 0;
    dSdPhi[I][1] = 0;
    dSdPhi[I][2] = 0;
    dSdPhi[I][3] = 0;
  }

  double y = fermiOps->getYN();
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  double VolumeNorm = 1.0 / (L0*L1*L2*L3);
  
  //If no Fermion-Coupling...
  if (y == 0) {
    Complex* p = auxVectors[0];
    for (int I=0; I<VLxtrSize; I++) {
      p[I].x = 0;
      p[I].y = 0;
    }
    S = 0;
    dSdOmega = p;

    return dSdOmega;
  }

  //Calculate results new...
  if (quasiHermiteanMode) {
    if (fermiOps->isMultiThreadedOpsActivated()) {
      dSdOmega = calcDistributedAllForcesQuasiHermiteanMode(polynomSlot, phi);
    } else {
      dSdOmega = calcAllForcesQuasiHermiteanMode(polynomSlot, phi);
    }
  } else {
    if (fermiOps->isMultiThreadedOpsActivated()) {
      printf("ERROR in pHMCForce::calcAllForces: not implemented for qhm=false yet\n");
      exit(0);
    } else {
      dSdOmega = calcAllForcesNotHermiteanMode(polynomSlot, phi);
    }
  }
  
  //Rescale Phi-Force
  double fac = 2*nu[polynomSlot]*y*VolumeNorm;
  for (int I=0; I<L0*L1*L2*L3; I++) {
    dSdPhi[I][0] *= fac;
    dSdPhi[I][1] *= fac;
    dSdPhi[I][2] *= fac;
    dSdPhi[I][3] *= fac;
  }
    
  //Omega-Action and Omega-Force
  //Omega - Action
  Complex dummy;
  SSE_ComplexScalarProduct(L0,L1,L2,L3,xtraSize1,xtraSize2,xtraSize3, omegaField, dSdOmega, dummy);
  double dummy2 = 0;
  if ((polynomSlot%2)==1) {
    SSE_ComplexSquareNorm(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, omegaField, dummy2);  
  }  
  S = 0.5*nu[polynomSlot]*dummy.x - 0.5*dummy2;
  
  //Omega - Force
  fac = nu[polynomSlot];    
  if ((polynomSlot%2)==1) {
    for (int I=0; I<VLxtrSize; I++) {
      dSdOmega[I].x = fac*dSdOmega[I].x - omegaField[I].x;
      dSdOmega[I].y = fac*dSdOmega[I].y - omegaField[I].y;
    }
  } else {
    for (int I=0; I<VLxtrSize; I++) {
      dSdOmega[I].x *= fac;
      dSdOmega[I].y *= fac;
    }
  }

  //Logging Omega-Force
  if ((omegaMassAdaptionMode) && (!tuneMode)){  
    double forceStrength = 0;
    for (int I2=0; I2<8; I2++) {
      forceStrength += dSdOmega[I2].x*dSdOmega[I2].x + dSdOmega[I2].y*dSdOmega[I2].y;
    }
    forceStrength = sqrt(forceStrength);
    
    if (LogLevel>4) printf("Omega Force Strength for polynomial %d on force %d on node %d measured:  %f\n",polynomSlot,ID,ownNodeID,forceStrength);
    omegaForceStrengthSumPREC[polynomSlot] += forceStrength;
    omegaForceStrengthSqrSumPREC[polynomSlot] += sqr(forceStrength);
    omegaForceStrengthCountPREC[polynomSlot]++;
    omegaForceStrengthSumGLOBAL[polynomSlot] += forceStrength;
    omegaForceStrengthSqrSumGLOBAL[polynomSlot] += sqr(forceStrength);
    omegaForceStrengthCountGLOBAL[polynomSlot]++;    
       
    omegaForceStrengthAnalysisPREC[polynomSlot]->addDataVectorFromOmegafourierTrafoSPECIAL(dSdOmega);
    omegaForceStrengthAnalysisGLOBAL[polynomSlot]->addDataVectorFromOmegafourierTrafoSPECIAL(dSdOmega);
  }

  //Put Result into Storage?
  storedActOmegaAction[polynomSlot] = S;
  if (storedActdSdOmega[polynomSlot] != NULL) {
    storedActForcesAvail[polynomSlot] = true;
    int VLSize = fermiOps->getVectorLength();
    int VLxtrSize = fermiOps->getVectorLengthXtrSize();
    SSE_ZCopy(VLxtrSize, dSdOmega, 1, storedActdSdOmega[polynomSlot], 1);  
    SSE_ZCopy(VLSize/4, (Complex*)dSdPhi, 1, storedActdSdPhi[polynomSlot], 1);  
  }
  
  return dSdOmega;
}



/**
* ATTENTION: This force is the force without the factor 1/2.
*
**/
Complex* pHMCForce::calcDistributedAllForcesQuasiHermiteanMode(int polynomSlot, double* phi) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  //Calculate results new...
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double VolumeNorm = 1.0 / (L0*L1*L2*L3);
  
  bool precRUse;
  double RPrecm, RPrecf;
  fermiOps->getRPreconditionerParameter(precRUse, RPrecm, RPrecf);

  //copy data to distributed memory structures
  MultiThreadedOperations* threadedOps = fermiOps->getMultiThreadedOps();
  threadedOps->copyComplexVectorToUniformlyDistributedMemory((Complex*) phi, Distributed_Phi, 2*L0*L1*L2*L3);
  threadedOps->zeroUniformlyDistributedComplexVector(Distributed_dSdPhiSum, 4*L0*L1*L2*L3);
  fermiOps->copyFermionVectorToDistributedFermionVector(omegaField, Distributed_omegaField);

  //Get Pointers to intermediate memory locations
  DistributedMemoryObject* p1 = NULL;
  DistributedMemoryObject* p2 = Distributed_auxVectors[0];
  DistributedMemoryObject* p3 = Distributed_auxVectors[1];
  DistributedMemoryObject* p4 = Distributed_auxVectors[2];
  DistributedMemoryObject* p5 = Distributed_auxVectors[3];
  DistributedMemoryObject* p6 = NULL;
  DistributedMemoryObject* p7 = NULL;
  DistributedMemoryObject* p8 = NULL;
  
  Complex a1(1,0);
  Complex a2;
  Complex a3(1,0);


  if (precRUse) {
    fermiOps->executeDistributedRPreconditionerMatrixMultiplication(Distributed_omegaField, p2, false, false, true);
  } else {
    threadedOps->copyFermionVector(Distributed_omegaField, p2, L0, L1, L2, L3);
  }  
  fermiOps->executeDistributedDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, p5, true);  
  fermiOps->executeDistributedPiModeRemoverOperator(p2, p2, true); 
  fermiOps->performDistributedFFT(p2, Distributed_fastEvalVectors[0], ExtremeFFT4D_Backward);     
  p8 = Distributed_omegaField;

  //Weitere Eintraege fuer erste Polynomhaelfte
  for (int I=0; I<approxpolyDegree[polynomSlot]/2; I++) {
    int baseInd = 2*I;
    p1 = Distributed_fastEvalVectors[baseInd];
    p6 = Distributed_fastEvalVectors[baseInd+1];
    p7 = Distributed_fastEvalVectors[baseInd+2];
    a2 = (-1.0) * approxpolyRoots[polynomSlot][I];
  
    fermiOps->executeDistributedYBDaggerOperatorMultiplication(p1, p3, VolumeNorm, Distributed_Phi);
    fermiOps->performDistributedFFT(p3, p2, ExtremeFFT4D_Forward);
    fermiOps->executeDistributedPiModeRemoverOperator(p2, p2, true);
        
    threadedOps->vectorAdditionOfFermionVectors(p5, p2, a1, L0, L1, L2, L3);
    
    if (precRUse) {
      fermiOps->executeDistributedRPreconditionerMatrixMultiplication(p2, p2, true, false, true);
    }

    fermiOps->executeDistributedDiracUnityMinusDiracQuasiInverseMatrixMultiplication(p2, p5, true);
    
    fermiOps->executeDistributedPiModeRemoverOperator(p2, p2, true);
    fermiOps->performDistributedFFT(p2, p6, ExtremeFFT4D_Backward);   
    fermiOps->executeDistributedYBOperatorMultiplication(p6, p3, VolumeNorm, Distributed_Phi);
    fermiOps->performDistributedFFT(p3, p2, ExtremeFFT4D_Forward);
    fermiOps->executeDistributedPiModeRemoverOperator(p2, p2, true);
    
    threadedOps->vectorAdditionOfFermionVectors(p2, p5, a1, L0, L1, L2, L3);    

    if (precRUse) {
      fermiOps->executeDistributedRPreconditionerMatrixMultiplication(p5, p5, false, false, true);
    }
    
    threadedOps->vectorAdditionOfFermionVectors(p8, p5, a2, L0, L1, L2, L3);
    p8 = p5;
    
    if (precRUse) {
      fermiOps->executeDistributedRPreconditionerMatrixMultiplication(p5, p2, false, false, true);
    } else {
      threadedOps->copyFermionVector(p5, p2, L0, L1, L2, L3);
    }
    
    fermiOps->executeDistributedDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, p4, true);  
    
    fermiOps->executeDistributedPiModeRemoverOperator(p2, p2, true); 
    fermiOps->performDistributedFFT(p2, p7, ExtremeFFT4D_Backward);     

    p1 = p4;
    p4 = p5;
    p5 = p1;    
  }


  p1 = Distributed_fastEvalVectors[approxpolyDegree[polynomSlot]];
  for (int I=approxpolyDegree[polynomSlot]/2; I<approxpolyDegree[polynomSlot]; I++) {
    int baseInd = 2*approxpolyDegree[polynomSlot]-2*I-1;
    p6 = Distributed_fastEvalVectors[baseInd];
    p7 = Distributed_fastEvalVectors[baseInd-1];
    a2 = (-1.0) * approxpolyRoots[polynomSlot][I];

    //Erste Kraft
    fermiOps->executeDistributedMultiplicationVectorWithDerivativesOfB(p1, p6, Distributed_dSdPhiDummy);    
    threadedOps->vectorAdditionOfUniformlyDistributedComplexVectors(Distributed_dSdPhiDummy, Distributed_dSdPhiSum, a3, 4*L0*L1*L2*L3);

    fermiOps->executeDistributedYBDaggerOperatorMultiplication(p1, p3, VolumeNorm, Distributed_Phi);
    fermiOps->performDistributedFFT(p3, p2, ExtremeFFT4D_Forward);
    fermiOps->executeDistributedPiModeRemoverOperator(p2, p2, true);

    threadedOps->vectorAdditionOfFermionVectors(p5, p2, a1, L0, L1, L2, L3);
    
    if (precRUse) {
      fermiOps->executeDistributedRPreconditionerMatrixMultiplication(p2, p2, true, false, true);
    }
    
    fermiOps->executeDistributedDiracUnityMinusDiracQuasiInverseMatrixMultiplication(p2, p5, true);
   
    fermiOps->executeDistributedPiModeRemoverOperator(p2, p2, true);
    fermiOps->performDistributedFFT(p2, p3, ExtremeFFT4D_Backward);   
   
    //Zweite Kraft
    fermiOps->executeDistributedMultiplicationVectorWithDerivativesOfB(p7, p3, Distributed_dSdPhiDummy);    
    threadedOps->vectorAdditionOfUniformlyDistributedComplexVectors(Distributed_dSdPhiDummy, Distributed_dSdPhiSum, a3, 4*L0*L1*L2*L3);

    fermiOps->executeDistributedYBOperatorMultiplication(p3, p3, VolumeNorm, Distributed_Phi);
    fermiOps->performDistributedFFT(p3, p2, ExtremeFFT4D_Forward);
    fermiOps->executeDistributedPiModeRemoverOperator(p2, p2, true);

    threadedOps->vectorAdditionOfFermionVectors(p2, p5, a1, L0, L1, L2, L3);
    
    if (precRUse) {
      fermiOps->executeDistributedRPreconditionerMatrixMultiplication(p5, p5, false, false, true);
    }

    threadedOps->vectorAdditionOfFermionVectors(p8, p5, a2, L0, L1, L2, L3);
    p8 = p5;
    
    if (I<approxpolyDegree[polynomSlot]-1) {
      if (precRUse) {
        fermiOps->executeDistributedRPreconditionerMatrixMultiplication(p5, p2, false, false, true);
      } else {
        threadedOps->copyFermionVector(p5, p2, L0, L1, L2, L3);
      }
    
      fermiOps->executeDistributedDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, p4, true);  
    
      fermiOps->executeDistributedPiModeRemoverOperator(p2, p2, true); 
      fermiOps->performDistributedFFT(p2, p7, ExtremeFFT4D_Backward);     
    }
        
    p1 = p4;
    p4 = p5;
    p5 = p1;
    p1 = p7;
  }

  //Retrieve results from distributed memory objects
  fermiOps->copyDistributedFermionVectorToFermionVector(p8, auxVectors[0]);   
  dSdOmega = auxVectors[0];
  Complex* v = auxVectors[1];  
  double* dsp = (double*) dSdPhi;
  threadedOps->calcSumOfUniformlyDistributedComplexVectors(Distributed_dSdPhiSum, v, 4*L0*L1*L2*L3);
  for (int I=0; I<4*L0*L1*L2*L3; I++) {
    dsp[I] = v[I].x;  
  }

  addPerformanceProfilingItem("pHMCForce::calcDistributedAllForcesQuasiHermiteanMode", performanceProfilerStartCycle, 0);
  return dSdOmega;
}


/**
* ATTENTION: This force is the force without the factor 1/2.
*
**/
Complex* pHMCForce::calcAllForcesQuasiHermiteanMode(int polynomSlot, double* phi) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  //Calculate results new...
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double y = fermiOps->getYN();
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  double VolumeNorm = 1.0 / (L0*L1*L2*L3);
  int I, I2;
  
  bool precRUse;
  double RPrecm, RPrecf;
  fermiOps->getRPreconditionerParameter(precRUse, RPrecm, RPrecf);

  double massSplit = fermiOps->getMassSplitRatio();

  Complex* p1 = NULL;
  Complex* p2 = auxVectors[0];
  Complex* p3 = auxVectors[1];
  Complex* p4 = auxVectors[2];
  Complex* p5 = auxVectors[3];
  Complex* p6 = NULL;
  Complex* p7 = NULL;
  Complex* p8 = NULL;
  
  Complex a1(1,0);
  Complex a2;


  if (precRUse) {
    fermiOps->executeRPreconditionerMatrixMultiplication(omegaField, p2, false, false, true);
  } else {
    SSE_ZCopy(VLxtrSize, omegaField, 1, p2, 1);
  }  
  fermiOps->executeDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, p5, true);  
  fermiOps->executePiModeRemoverOperator(p2, p2, true); 
  fermiOps->performFFT(p2, fastEvalVectors[0], ExtremeFFT4D_Backward);     
  p8 = omegaField;

  //Weitere Eintraege fuer erste Polynomhaelfte
  for (I=0; I<approxpolyDegree[polynomSlot]/2; I++) {
    int baseInd = 2*I;
    p1 = fastEvalVectors[baseInd];
    p6 = fastEvalVectors[baseInd+1];
    p7 = fastEvalVectors[baseInd+2];
    a2 = (-1.0) * approxpolyRoots[polynomSlot][I];
  
    if (massSplit == 1) {
      perform_yBDagger(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y*VolumeNorm, phi, p1, p3);    
    } else {
      perform_yBsplitDagger(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y*VolumeNorm, massSplit, phi, p1, p3);        
    }
    fermiOps->performFFT(p3, p2, ExtremeFFT4D_Forward);
    fermiOps->executePiModeRemoverOperator(p2, p2, true);
    
    SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, a1, p5, p2);
    
    if (precRUse) {
      fermiOps->executeRPreconditionerMatrixMultiplication(p2, p2, true, false, true);
    }

    fermiOps->executeDiracUnityMinusDiracQuasiInverseMatrixMultiplication(p2, p5, true);
    
    fermiOps->executePiModeRemoverOperator(p2, p2, true);
    fermiOps->performFFT(p2, p6, ExtremeFFT4D_Backward);   
    if (massSplit == 1) {
      perform_yB(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y*VolumeNorm, phi, p6, p3);    
    } else {
      perform_yBsplit(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y*VolumeNorm, massSplit, phi, p6, p3);        
    }
    fermiOps->performFFT(p3, p2, ExtremeFFT4D_Forward);
    fermiOps->executePiModeRemoverOperator(p2, p2, true);
    
    SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, a1, p2, p5);

    if (precRUse) {
      fermiOps->executeRPreconditionerMatrixMultiplication(p5, p5, false, false, true);
    }
    
    SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, a2, p8, p5);
    p8 = p5;
    
    if (precRUse) {
      fermiOps->executeRPreconditionerMatrixMultiplication(p5, p2, false, false, true);
    } else {
      SSE_ZCopy(VLxtrSize, p5, 1, p2, 1);
    }
    
    fermiOps->executeDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, p4, true);  
    
    fermiOps->executePiModeRemoverOperator(p2, p2, true); 
    fermiOps->performFFT(p2, p7, ExtremeFFT4D_Backward);     

    p1 = p4;
    p4 = p5;
    p5 = p1;    
  }


  p1 = fastEvalVectors[approxpolyDegree[polynomSlot]];
  for (I=approxpolyDegree[polynomSlot]/2; I<approxpolyDegree[polynomSlot]; I++) {
    int baseInd = 2*approxpolyDegree[polynomSlot]-2*I-1;
    p6 = fastEvalVectors[baseInd];
    p7 = fastEvalVectors[baseInd-1];
    a2 = (-1.0) * approxpolyRoots[polynomSlot][I];

    //Erste Kraft
    fermiOps->executeMultiplicationVectorWithDerivativesOfB(p1, p6, p3);
    int index=0;
    for (I2=0; I2<L0*L1*L2*L3; I2++) {
      dSdPhi[I2][0] += p3[index+0].x;
      dSdPhi[I2][1] += p3[index+1].x;
      dSdPhi[I2][2] += p3[index+2].x;
      dSdPhi[I2][3] += p3[index+3].x;
      index += 4;
    }

    if (massSplit == 1) {
      perform_yBDagger(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y*VolumeNorm, phi, p1, p3);    
    } else {
      perform_yBsplitDagger(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y*VolumeNorm, massSplit, phi, p1, p3);        
    }
    fermiOps->performFFT(p3, p2, ExtremeFFT4D_Forward);
    fermiOps->executePiModeRemoverOperator(p2, p2, true);

    SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, a1, p5, p2);
    
    if (precRUse) {
      fermiOps->executeRPreconditionerMatrixMultiplication(p2, p2, true, false, true);
    }
    
    fermiOps->executeDiracUnityMinusDiracQuasiInverseMatrixMultiplication(p2, p5, true);
   
    fermiOps->executePiModeRemoverOperator(p2, p2, true);
    fermiOps->performFFT(p2, p3, ExtremeFFT4D_Backward);   
   
    //Zweite Kraft
    fermiOps->executeMultiplicationVectorWithDerivativesOfB(p7, p3, p6);
    index=0;
    for (I2=0; I2<L0*L1*L2*L3; I2++) {
      dSdPhi[I2][0] += p6[index+0].x;
      dSdPhi[I2][1] += p6[index+1].x;
      dSdPhi[I2][2] += p6[index+2].x;
      dSdPhi[I2][3] += p6[index+3].x;
      index += 4;
    }

    if (massSplit == 1) {
      perform_yB(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y*VolumeNorm, phi, p3, p3);   
    } else {
      perform_yBsplit(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y*VolumeNorm, massSplit, phi, p3, p3);       
    }
    fermiOps->performFFT(p3, p2, ExtremeFFT4D_Forward);
    fermiOps->executePiModeRemoverOperator(p2, p2, true);

    SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, a1, p2, p5);
    
    if (precRUse) {
      fermiOps->executeRPreconditionerMatrixMultiplication(p5, p5, false, false, true);
    }

    SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, a2, p8, p5);
    p8 = p5;
    
    if (I<approxpolyDegree[polynomSlot]-1) {
      if (precRUse) {
        fermiOps->executeRPreconditionerMatrixMultiplication(p5, p2, false, false, true);
      } else {
        SSE_ZCopy(VLxtrSize, p5, 1, p2, 1);
      }
    
      fermiOps->executeDiracUnityMinusDiracQuasiInverseDaggerMatrixMultiplication(p2, p4, true);  
    
      fermiOps->executePiModeRemoverOperator(p2, p2, true); 
      fermiOps->performFFT(p2, p7, ExtremeFFT4D_Backward);     
    }
        
    p1 = p4;
    p4 = p5;
    p5 = p1;
    p1 = p7;
  }

  dSdOmega = p8;
  addPerformanceProfilingItem("pHMCForce::calcAllForcesQuasiHermiteanMode", performanceProfilerStartCycle, 0);
  return dSdOmega;
}


/**
* ATTENTION: This force is the force without the factor 1/2.
*
**/
Complex* pHMCForce::calcAllForcesNotHermiteanMode(int polynomSlot, double* phi) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  //Calculate results new...
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double y = fermiOps->getYN();
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  double VolumeNorm = 1.0 / (L0*L1*L2*L3);
  int I, I2;
  
  bool precUse;
  double precM, precS;
  fermiOps->getPreconditionerParameter(precUse, precM, precS);

  bool precQUse;
  double QPrecMu, QPrecBeta;
  fermiOps->getQPreconditionerParameter(precQUse, QPrecMu, QPrecBeta);
  
  double massSplit = fermiOps->getMassSplitRatio();

  double rho, r;
  fermiOps->getDiracParameters(rho, r);  
  Complex* p1 = NULL;
  Complex* p2 = auxVectors[0];
  Complex* p3 = auxVectors[1];
  Complex* p4 = auxVectors[2];
  Complex* p5 = auxVectors[3];
  Complex* p6 = NULL;
  Complex* p7 = NULL;
  
  Complex a1(-2*rho,0);
  Complex a2;

  //Erster Eintrag in fastEvalVectors
  SSE_ZCopy(VLxtrSize, omegaField, 1, p2, 1);  

  if (precQUse) {
    fermiOps->executeQPreconditionerMatrixMultiplication(p2, p3, false, true);
  } else {
    SSE_ZCopy(VLxtrSize, p2, 1, p3, 1);
  }
  fermiOps->performFFT(p3, fastEvalVectors[0], ExtremeFFT4D_Backward);

  //Weitere Eintraege fuer erste Polynomhaelfte
  for (I=0; I<approxpolyDegree[polynomSlot]/2; I++) {
    int baseInd = 2*I;
    p1 = fastEvalVectors[baseInd];
    p6 = fastEvalVectors[baseInd+1];
    p7 = fastEvalVectors[baseInd+2];
    a2 = (-1.0) * approxpolyRoots[polynomSlot][I];

    fermiOps->executeDiracDaggerMatrixMultiplication(p3, p4, true);
    if (massSplit == 1) {
      perform_yBDagger(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y, phi, p1, p5);
    } else {
      perform_yBsplitDagger(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y, massSplit, phi, p1, p5);
    }
    fermiOps->performFFT(p5, p3, ExtremeFFT4D_Forward);
    fermiOps->executeDiracDaggerMatrixMultiplication(p3, p5, true);
    
    double f0 = VolumeNorm;
    double f1 = 2*rho;
    if (precUse) {
      double g1 = f1 / VolumeNorm;
      for (I2=0; I2<VLxtrSize; I2++) {
        p4[I2].x = (p5[I2].x - f1*p3[I2].x) - g1*p4[I2].x;
        p4[I2].y = (p5[I2].y - f1*p3[I2].y) - g1*p4[I2].y;
      }
      fermiOps->executeFermionMatrixStaticInverseMultiplication(p4, p3, true, true);
    } else {
      for (I2=0; I2<VLxtrSize; I2++) {
        p3[I2].x = f0*(p5[I2].x - f1*p3[I2].x) - f1*p4[I2].x;
        p3[I2].y = f0*(p5[I2].y - f1*p3[I2].y) - f1*p4[I2].y;
      }
    }
    fermiOps->executeDiracMatrixMultiplication(p3, p4, true);
    for (I2=0; I2<VLxtrSize; I2++) {
      p3[I2].x = p4[I2].x - f1*p3[I2].x;
      p3[I2].y = p4[I2].y - f1*p3[I2].y;
    }
    fermiOps->performFFT(p3, p6, ExtremeFFT4D_Backward);
    if (massSplit == 1) {    
      perform_yB(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y, phi, p6, p5);
    } else {
      perform_yBsplit(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y, massSplit, phi, p6, p5);
    }
    fermiOps->performFFT(p5, p3, ExtremeFFT4D_Forward);
    for (I2=0; I2<VLxtrSize; I2++) {
      p3[I2].x = f0*p3[I2].x - f1*p4[I2].x;
      p3[I2].y = f0*p3[I2].y - f1*p4[I2].y;
    }
    if (precQUse) {
      fermiOps->executeQPreconditionerMatrixMultiplication(p3, p3, false, true);   
    }
    SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, a2, p2, p3);
    if (precQUse) {
      fermiOps->executeQPreconditionerMatrixMultiplication(p3, p4, false, true);   
    } else {
      SSE_ZCopy(VLxtrSize, p3, 1, p4, 1);
    }
    fermiOps->performFFT(p4, p7, ExtremeFFT4D_Backward);
    
    p1 = p2;
    p2 = p3;
    p3 = p4;
    p4 = p1;
  }


  p1 = fastEvalVectors[approxpolyDegree[polynomSlot]];
  for (I=approxpolyDegree[polynomSlot]/2; I<approxpolyDegree[polynomSlot]; I++) {
    int baseInd = 2*approxpolyDegree[polynomSlot]-2*I-1;
    p6 = fastEvalVectors[baseInd];
    p7 = fastEvalVectors[baseInd-1];
    a2 = (-1.0) * approxpolyRoots[polynomSlot][I];

    //Erste Kraft
    fermiOps->executeMultiplicationVectorWithDerivativesOfB(p1, p6, p4);
    int index=0;
    for (I2=0; I2<L0*L1*L2*L3; I2++) {
      dSdPhi[I2][0] += p4[index+0].x;
      dSdPhi[I2][1] += p4[index+1].x;
      dSdPhi[I2][2] += p4[index+2].x;
      dSdPhi[I2][3] += p4[index+3].x;
      index += 4;
    }

    fermiOps->executeDiracDaggerMatrixMultiplication(p3, p4, true);
    if (massSplit == 1) {     
      perform_yBDagger(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y, phi, p1, p5);
    } else {
      perform_yBsplitDagger(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y, massSplit, phi, p1, p5);    
    }
    fermiOps->performFFT(p5, p3, ExtremeFFT4D_Forward);
    fermiOps->executeDiracDaggerMatrixMultiplication(p3, p5, true);
    
    double f0 = VolumeNorm;
    double f1 = 2*rho;
    if (precUse) {
      double g1 = f1 / VolumeNorm;
      for (I2=0; I2<VLxtrSize; I2++) {
        p4[I2].x = (p5[I2].x - f1*p3[I2].x) - g1*p4[I2].x;
        p4[I2].y = (p5[I2].y - f1*p3[I2].y) - g1*p4[I2].y;
      }
      fermiOps->executeFermionMatrixStaticInverseMultiplication(p4, p3, true, true);
    } else {
      for (I2=0; I2<VLxtrSize; I2++) {
        p3[I2].x = f0*(p5[I2].x - f1*p3[I2].x) - f1*p4[I2].x;
        p3[I2].y = f0*(p5[I2].y - f1*p3[I2].y) - f1*p4[I2].y;
      }
    }
    fermiOps->executeDiracMatrixMultiplication(p3, p4, true);
    for (I2=0; I2<VLxtrSize; I2++) {
      p3[I2].x = p4[I2].x - f1*p3[I2].x;
      p3[I2].y = p4[I2].y - f1*p3[I2].y;
    }
    fermiOps->performFFT(p3, p5, ExtremeFFT4D_Backward);

    //Zweite Kraft
    fermiOps->executeMultiplicationVectorWithDerivativesOfB(p7, p5, p3);
    index=0;
    for (I2=0; I2<L0*L1*L2*L3; I2++) {
      dSdPhi[I2][0] += p3[index+0].x;
      dSdPhi[I2][1] += p3[index+1].x;
      dSdPhi[I2][2] += p3[index+2].x;
      dSdPhi[I2][3] += p3[index+3].x;
      index += 4;
    }

    if (massSplit == 1) {       
      perform_yB(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y, phi, p5, p5);
    } else {
      perform_yBsplit(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, y, massSplit, phi, p5, p5);
    }
    fermiOps->performFFT(p5, p3, ExtremeFFT4D_Forward);
    for (I2=0; I2<VLxtrSize; I2++) {
      p3[I2].x = f0*p3[I2].x - f1*p4[I2].x;
      p3[I2].y = f0*p3[I2].y - f1*p4[I2].y;
    }
    if (precQUse) {
      fermiOps->executeQPreconditionerMatrixMultiplication(p3, p3, false, true);   
    }
    SSE_ComplexVectorAddition(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, a2, p2, p3);
    if (precQUse) {
      fermiOps->executeQPreconditionerMatrixMultiplication(p3, p4, false, true);   
    } else {
      SSE_ZCopy(VLxtrSize, p3, 1, p4, 1);
    }
    
    p1 = p2;
    p2 = p3;
    p3 = p4;
    p4 = p1;
    
    fermiOps->performFFT(p3, p5, ExtremeFFT4D_Backward);
    p1 = p5;
  }

  dSdOmega = p2;
  addPerformanceProfilingItem("pHMCForce::calcAllForcesNotHermiteanMode", performanceProfilerStartCycle, 0);
  return dSdOmega;
}


void pHMCForce::sampleOmegaFieldsPurelyGaussian() {
  if (LogLevel>4) printf("Sampling omega fields purely Gaussian with own random generator on force %d on node %d\n",ID,ownNodeID);
  fermiOps->fillGaussRandomVector(auxVectors[3],-1); 
  fermiOps->transformToXtraSizeArray(auxVectors[3], omegaField);    
}


void pHMCForce::sampleOmegaFields(double* phiField) {
  if (LogLevel>4) printf("Sampling omega fields with own random generator on force %d on node %d\n",ID,ownNodeID);
  fermiOps->fillGaussRandomVector(auxVectors[3],-1); 
  sampleOmegaFields(phiField, auxVectors[3]);
}


void pHMCForce::sampleOmegaFields(double* phiField, Complex* data) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  if (LogLevel>4) printf("Sampling omega fields with given random numbers on force %d on node %d\n",ID,ownNodeID);
  fermiOps->transformToXtraSizeArray(data,auxVectors[3]);  
  data = auxVectors[3];
  
  if (fermiOps->isMultiThreadedOpsActivated()) {
    fermiOps->copyFermionVectorToDistributedFermionVector(data, Distributed_auxVectors[3]);
    fermiOps->copyPhiFieldUniformlyToDistributedPhiField(phiField, Distributed_Phi);
    applyDistributedInverseSQRTPolynomialCHEBYSHEV(Distributed_auxVectors[3], Distributed_omegaField, Distributed_Phi, 0);
    fermiOps->copyDistributedFermionVectorToFermionVector(Distributed_omegaField, omegaField);
  } else {
    applyInverseSQRTPolynomialCHEBYSHEV(data, omegaField, phiField, 0);
  }
  addPerformanceProfilingItem("pHMCForce::sampleOmegaFields", performanceProfilerStartCycle, 0);
}


void pHMCForce::sampleOmegaMomenta() {
  if (LogLevel>4) printf("Sampling omega momenta with own random generator on force %d on node %d\n",ID,ownNodeID);
  fermiOps->fillGaussRandomVector(omegaMomenta,-1); 
  sampleOmegaMomenta(omegaMomenta);
}


void pHMCForce::sampleOmegaMomenta(Complex* data) {
  if (LogLevel>4) printf("Sampling omega momenta with given random numbers on force %d on node %d\n",ID,ownNodeID);
  int I, I2;
  int VL = fermiOps->getVectorLength();
  double f = 1.0 / sqrt(theta);
  if (theta == 0) f = 0;

  int count = 0;
  for (I=0; I<VL/8; I++) {
    double fac = f*sqrt(MomentaMasses[I]);
    for (I2=0; I2<8; I2++) {
      omegaMomenta[count].x = fac*data[count].x;
      omegaMomenta[count].y = fac*data[count].y;  
      count++;
    }
  }
  fermiOps->transformToXtraSizeArray(omegaMomenta, omegaMomenta);
}


void pHMCForce::multiplyOmegaMomenta(double fac) {
  if (LogLevel>4) printf("Multiplying omega momenta with factor=%f on force %d on node %d\n",fac,ID,ownNodeID);
  int I;
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  for (I=0; I<VLxtrSize; I++) {
    omegaMomenta[I].x *= fac;
    omegaMomenta[I].y *= fac;  
  }
}


void pHMCForce::drawAndWasteOmegaMomentaRandomNumbers() {
  if (LogLevel>4) printf("Drawing and wasting omega momenta from own random generator on force %d on node %d\n",ID,ownNodeID);
  fermiOps->fillGaussRandomVector(auxVectors[0], -1); 
}


void pHMCForce::randomOmegaField() {
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  int I;

  for (I=0; I<VLxtrSize; I++) {
    omegaField[I].x = 2.0*(AdvancedZufall(AdvancedSeed)-0.5);
    omegaField[I].y = 2.0*(AdvancedZufall(AdvancedSeed)-0.5);    
  }  
}


void pHMCForce::saveOmegaField() {
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  SSE_ZCopy(VLxtrSize, omegaField, 1, omegaOLD, 1); 
  SSE_ZCopy(VLxtrSize, omegaMomenta, 1, omegaMomentaOLD, 1);    
}


void pHMCForce::restoreOmegaField() {
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  SSE_ZCopy(VLxtrSize, omegaOLD, 1, omegaField, 1);  
  SSE_ZCopy(VLxtrSize, omegaMomentaOLD, 1, omegaMomenta, 1);      
}


void pHMCForce::doOmegaStep(double eps) {
  double fac = eps*theta;

  if ((theta==0) || (eps==0)) return;
  
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  int count = 0;
  int I3, I2, I1, I0, i;
  int p = 0;
  int d0 = 8*(L3 + xtraSize3)*(L2 + xtraSize2)*(L1 + xtraSize1) - L1*8*(L3 + xtraSize3)*(L2 + xtraSize2);
  int d1 = 8*(L3 + xtraSize3)*(L2 + xtraSize2) - L2*8*(L3 + xtraSize3);
  int d2 = 8*(L3 + xtraSize3) - L3*8;
  int d3 = 8;
  for (I0=0; I0<L0; I0++) {
    for (I1=0; I1<L1; I1++) {
      for (I2=0; I2<L2; I2++) {
        for (I3=0; I3<L3; I3++) {
          double f = fac / MomentaMasses[count];
	  for (i=0; i<8; i++) {
            omegaField[p+i].x += f*omegaMomenta[p+i].x;
            omegaField[p+i].y += f*omegaMomenta[p+i].y;
  	  }
	  count++;
          p += d3;
        }
        p += d2;
      }
      p += d1;
    }
    p += d0;
  }  
  
  clearActStore();
}


void pHMCForce::doOmegaMomentumStep(double eps) {
  if ((theta==0) || (eps==0)) return;

  int I;
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  for (I=0; I<VLxtrSize; I++) {
    omegaMomenta[I].x -= eps*dSdOmega[I].x;
    omegaMomenta[I].y -= eps*dSdOmega[I].y;    
  }  
}


double pHMCForce::calcOmegaMomentumAction() {
  double S = 0;
  
  if (theta != 0) {
    int L0 = fermiOps->get1DSizeL0();
    int L1 = fermiOps->get1DSizeL1();
    int L2 = fermiOps->get1DSizeL2();
    int L3 = fermiOps->get1DSizeL3();

    if (omegaMassAdaptionMode) {
      int count = 0;
      int I3, I2, I1, I0, i;
      int p = 0;
      int d0 = 8*(L3 + xtraSize3)*(L2 + xtraSize2)*(L1 + xtraSize1) - L1*8*(L3 + xtraSize3)*(L2 + xtraSize2);
      int d1 = 8*(L3 + xtraSize3)*(L2 + xtraSize2) - L2*8*(L3 + xtraSize3);
      int d2 = 8*(L3 + xtraSize3) - L3*8;
      int d3 = 8;
      for (I0=0; I0<L0; I0++) {
        for (I1=0; I1<L1; I1++) {
          for (I2=0; I2<L2; I2++) {
            for (I3=0; I3<L3; I3++) {
              double f = 1.0 / MomentaMasses[count];
              for (i=0; i<8; i++) {
	        S += f*(omegaMomenta[p+i].x*omegaMomenta[p+i].x + omegaMomenta[p+i].y*omegaMomenta[p+i].y);
  	      }
  	      count++;
              p += d3;
            }
            p += d2;
          }
          p += d1;
        }
        p += d0;
      }  
    } else { 
      SSE_ComplexSquareNorm(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, omegaMomenta, S);
    }
  }
  
  setOmegaMomentumAction(0.5*theta*S);
  return getOmegaMomentumAction();
}


double pHMCForce_calcExactMMdagInverseSQRTomegaAction_Exponent;
double pHMCForce_calcExactMMdagInverseSQRTomegaAction_Function(double x) {
  return exp(-pHMCForce_calcExactMMdagInverseSQRTomegaAction_Exponent * log(x));
}


double pHMCForce::calcExactMMdagInverseSQRTomegaAction(double* phi, double alpha, int &neededIter) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double TOL = 1E-10;
  neededIter = 0;
  Complex res;
  bool b;
  pHMCForce_calcExactMMdagInverseSQRTomegaAction_Exponent = alpha;

  if (fermiOps->isMultiThreadedOpsActivated()) {
    fermiOps->copyFermionVectorToDistributedFermionVector(omegaField, Distributed_omegaField);
    fermiOps->copyPhiFieldUniformlyToDistributedPhiField(phi, Distributed_Phi);
    if (alpha==0.5) {
      b = fermiOps->applyDistributedFermionMatrixMMDaggerFunction(Distributed_omegaField, Distributed_auxVectors[0], Distributed_Phi, NULL, TOL, 0, neededIter, VectorCollectionCount, Distributed_VectorCollection, quasiHermiteanMode, true);
    } else {
      b = fermiOps->applyDistributedFermionMatrixMMDaggerFunction(Distributed_omegaField, Distributed_auxVectors[0], Distributed_Phi, &pHMCForce_calcExactMMdagInverseSQRTomegaAction_Function, TOL, 0, neededIter, VectorCollectionCount, Distributed_VectorCollection, quasiHermiteanMode, true);
    }
    MultiThreadedOperations* threadedOps = fermiOps->getMultiThreadedOps();
    threadedOps->scalarProductOfFermionVectors(Distributed_omegaField, Distributed_auxVectors[0], res, L0, L1, L2, L3);
  } else {
    if (alpha==0.5) {
      b = fermiOps->applyFermionMatrixMMDaggerFunction(omegaField, auxVectors[0], phi, NULL, TOL, 0, neededIter, VectorCollectionCount, VectorCollection, quasiHermiteanMode, true);
    } else {
      b = fermiOps->applyFermionMatrixMMDaggerFunction(omegaField, auxVectors[0], phi, &pHMCForce_calcExactMMdagInverseSQRTomegaAction_Function, TOL, 0, neededIter, VectorCollectionCount, VectorCollection, quasiHermiteanMode, true);
    }
    SSE_ComplexScalarProduct(L0, L1, L2, L3, xtraSize1, xtraSize2, xtraSize3, omegaField, auxVectors[0], res);
  }
  
  if (!b) {
    printf("ERROR in calcExactMMdagInverseSQRTomegaAction: No convergence\n");
    exit(0);
  }
  if (fabs(res.y)>TOL) {
    if (LogLevel>2) printf("WARNING in calcExactMMdagInverseSQRTomegaAction: Imaginary part of result greater than tolerance (%e): %e\n",TOL,res.y);
  }
  setExactOmegaMMdagInverseSQRTAction(0.5*res.x);
  addPerformanceProfilingItem("pHMCForce::calcExactMMdagInverseSQRTomegaAction", performanceProfilerStartCycle, 0);  
  return getExactOmegaMMdagInverseSQRTAction();
}


double pHMCForce::getOmegaMomentumAction() {
  return omegaMomentumAction;
}


void pHMCForce::setOmegaMomentumAction(double s) {
  omegaMomentumAction = s;
}

double pHMCForce::getNu(int polynomSlot) {
  return nu[polynomSlot];
}


double pHMCForce::getTheta() {
  return theta;
}


void pHMCForce::setTheta(double tht) {
  bool pri = (tht!=theta);
  theta = tht;
  if (pri) {
    if (LogLevel>2) printf("Setting theta to %1.5e (zero=%d) on force %d on node %d\n",theta,theta==0,ID,ownNodeID);
  }
}


void pHMCForce::plotApproxPoly(int polynomSlot, double eps, double lam, char* filename) {
  FILE* file = fopen(filename,"w");
  double p = 0;
  int I2;
  for (I2=0; I2<=1000; I2++) {
    int I;
    p = I2*lam/1000.0;
    Complex y = Complex(p,0) - approxpolyRoots[polynomSlot][0];
    for (I=1; I<approxpolyDegree[polynomSlot]; I++) {
       y = y * (Complex(p,0) - approxpolyRoots[polynomSlot][I]);
    }
    y = nu[polynomSlot]*y;
    fprintf(file,"%1.15f %1.15f %1.15f\n",p,y.x,y.y);
  }
  fclose(file);  
}


void pHMCForce::plotApproxPolyRoots(int polynomSlot, char* filename) {
  FILE* file = fopen(filename,"w");
  int I;
  for (I=0; I<approxpolyDegree[polynomSlot]; I++) {
    fprintf(file,"%1.15f %1.15f\n",approxpolyRoots[polynomSlot][I].x,approxpolyRoots[polynomSlot][I].y);
  }
  fclose(file);  
}


Complex* pHMCForce::getOmegaMomenta() {
  return omegaMomenta;
}


void pHMCForce::setNu(int polynomSlot, double n) {
  nu[polynomSlot] = n;
  if (LogLevel>3) printf("Setting nu for polynomial %d to %1.5e on force %d on node %d\n",polynomSlot, nu[polynomSlot],ID,ownNodeID);
}
  
  
void pHMCForce::setOmegaMassAdaptionMode(bool omMassAdapMode) {
  if ((!omegaMassAdaptionMode) && (omMassAdapMode)) {
    if (LogLevel>2) printf("Omega Mass Adaption mode on force %d on node %d activated!\n",ID,ownNodeID);
  }
  if ((omegaMassAdaptionMode) && (!omMassAdapMode)) {
    if (LogLevel>2) printf("Omega Mass Adaption mode on force %d on node %d deactivated!\n",ID,ownNodeID);
  }
  
  omegaMassAdaptionMode = omMassAdapMode;
}


void pHMCForce::calcOmegaMassAdaption() {
  if (omegaMassAdaptionMode) {
    int L0 = fermiOps->get1DSizeL0();
    int L1 = fermiOps->get1DSizeL1();
    int L2 = fermiOps->get1DSizeL2();
    int L3 = fermiOps->get1DSizeL3();
    if (omegaForceStrengthAnalysisPREC[0]->containsData()) {
      omegaForceStrengthAnalysisPREC[0]->getAverageVectorInflated(MomentaMasses);
      for (int I=0; I<L0*L1*L2*L3; I++) {
        MomentaMasses[I] = MomentaMasses[I]*MomentaMasses[I];
      }
    } else {
      for (int I=0; I<L0*L1*L2*L3; I++) MomentaMasses[I] = 1.0;
    }
    if (LogLevel>2) printf("Setting Omega MomentumMasses on force %d on node %d.\n",ID,ownNodeID);
  }
}


void pHMCForce::writeOmegaForceStrengthToDiskPREC(char *fileName) {
  int I;
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    char* name = new char[600];
    snprintf(name,600,"%sForce%d_Node%d_SubPol%d.dat",fileName,ID,ownNodeID,I);
    if (omegaForceStrengthAnalysisPREC[I]->containsData()) {
      omegaForceStrengthAnalysisPREC[I]->saveData(name);  
    }
    delete[] name;
  }
}


void pHMCForce::writeOmegaForceStrengthToDiskGLOBAL(char *fileName) {
  int I;
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    char* name = new char[600];
    snprintf(name,600,"%sForce%d_Node%d_SubPol%d.dat",fileName,ID,ownNodeID,I);
    if (omegaForceStrengthAnalysisGLOBAL[I]->containsData()) {
      omegaForceStrengthAnalysisGLOBAL[I]->saveData(name);  
    }
    delete[] name;
  }
}


void pHMCForce::readOmegaForceStrengthFromDiskPREC(char *fileName) {
  int I;
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    omegaForceStrengthSumPREC[I] = 0;
    omegaForceStrengthSqrSumPREC[I] = 0;
    omegaForceStrengthCountPREC[I] = 0;
  }
  
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    char* name = new char[600];
    snprintf(name,600,"%sForce%d_Node%d_SubPol%d.dat",fileName,ID,ownNodeID,I);

    omegaForceStrengthAnalysisPREC[I]->clearData();
    omegaForceStrengthAnalysisPREC[I]->loadData(name);  
    delete[] name;
  }
}


void pHMCForce::readOmegaForceStrengthFromDiskGLOBAL(char *fileName) {
  int I;
  
  for (I=0; I<pHMCForcePolynomialsMAX; I++) {
    char* name = new char[600];
    snprintf(name,600,"%sForce%d_Node%d_SubPol%d.dat",fileName,ID,ownNodeID,I);

    omegaForceStrengthAnalysisGLOBAL[I]->clearData(); 
    omegaForceStrengthAnalysisGLOBAL[I]->loadData(name);  
    delete[] name;
  }
}


void pHMCForce::activateForceStoring(int count, int* activateIndices) {
  for (int I=0; I<pHMCForcePolynomialsMAX; I++) {
    storedActOmegaAction[I] = NaN;
    storedTraStartOmegaAction[I] = NaN;
    storedActForcesAvail[I] = false;
    storedTraStartForcesAvail[I] = false;

    fermiOps->destroyFermionVector(storedActdSdOmega[I]);
    fermiOps->destroyFermionVector(storedTraStartdSdOmega[I]);
    fermiOps->destroyFermionVector(storedActdSdPhi[I]);
    fermiOps->destroyFermionVector(storedTraStartdSdPhi[I]);

    storedActdSdOmega[I] = NULL;
    storedTraStartdSdOmega[I] = NULL;
    storedActdSdPhi[I] = NULL;
    storedTraStartdSdPhi[I] = NULL;
  }
  
  if ((count==0) || (activateIndices == NULL)) return;
  if (!local) return;
    
  for (int I=0; I<count; I++) {
    storedActdSdOmega[activateIndices[I]] = fermiOps->createFermionVector();
    storedTraStartdSdOmega[activateIndices[I]] = fermiOps->createFermionVector();
    storedActdSdPhi[activateIndices[I]] = fermiOps->createFermionVector(2);
    storedTraStartdSdPhi[activateIndices[I]] = fermiOps->createFermionVector(2);
  }
}


void pHMCForce::clearActStore() {
  if (LogLevel>4) printf("Clear Act-Store on force %d on node %d.\n",ID,ownNodeID);
  for (int I=0; I<pHMCForcePolynomialsMAX; I++) {
    storedActOmegaAction[I] = NaN;
    storedActForcesAvail[I] = false;
  }
  omegaMomentumAction = NaN;
  exactOmegaMMdagInverseSQRTAction = NaN;
}


void pHMCForce::clearTraStartStore() {
  if (LogLevel>4) printf("Clear Tra-Start-Store on force %d on node %d.\n",ID,ownNodeID);
  for (int I=0; I<pHMCForcePolynomialsMAX; I++) {
    storedTraStartOmegaAction[I] = NaN;
    storedTraStartForcesAvail[I] = false;
  }
}


void pHMCForce::copyTraStartStoreToActStore() {
  if (LogLevel>4) printf("Copy Tra-Start-Store to Act-Store on force %d on node %d.\n",ID,ownNodeID);
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  int VLSize = fermiOps->getVectorLength();
  for (int I=0; I<pHMCForcePolynomialsMAX; I++) {
    if ((storedActdSdOmega[I] != NULL) && (storedTraStartdSdOmega[I] != NULL)) {
      SSE_ZCopy(VLxtrSize, storedTraStartdSdOmega[I], 1, storedActdSdOmega[I], 1);  
      storedActForcesAvail[I] = storedTraStartForcesAvail[I];
    }
    if ((storedActdSdPhi[I] != NULL) && (storedTraStartdSdPhi[I] != NULL)) {
      SSE_ZCopy(VLSize/4, storedTraStartdSdPhi[I], 1, storedActdSdPhi[I], 1);  
    }
    storedActOmegaAction[I] = storedTraStartOmegaAction[I];
  }
}


void pHMCForce::copyActStoreToTraStartStore() {
  if (LogLevel>4) printf("Copy Act-Store to Tra-Start-Store on force %d on node %d.\n",ID,ownNodeID);
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  int VLSize = fermiOps->getVectorLength();
  for (int I=0; I<pHMCForcePolynomialsMAX; I++) {
    if ((storedActdSdOmega[I] != NULL) && (storedTraStartdSdOmega[I] != NULL)) {
      SSE_ZCopy(VLxtrSize, storedActdSdOmega[I], 1, storedTraStartdSdOmega[I], 1);  
      storedTraStartForcesAvail[I] = storedActForcesAvail[I];
    }
    if ((storedActdSdPhi[I] != NULL) && (storedTraStartdSdPhi[I] != NULL)) {
      SSE_ZCopy(VLSize/4, storedActdSdPhi[I], 1, storedTraStartdSdPhi[I], 1);  
    }
    storedTraStartOmegaAction[I] = storedActOmegaAction[I];
  }
}


void pHMCForce::storeCurrentTopLevelForcesToTraStartStore() {
  int topLevelInd = -1;
  for (int I=0; I<pHMCForcePolynomialsMAX; I++) {
    if (storedTraStartdSdOmega[I] != NULL) {
      topLevelInd = I;
      break;
    }
  }
  if (topLevelInd == -1) return;
  
  int VLxtrSize = fermiOps->getVectorLengthXtrSize();
  int VLSize = fermiOps->getVectorLength();
  SSE_ZCopy(VLxtrSize, dSdOmega, 1, storedTraStartdSdOmega[topLevelInd], 1);  
  SSE_ZCopy(VLSize/4, (Complex*)dSdPhi, 1, storedTraStartdSdPhi[topLevelInd], 1);  
  storedTraStartForcesAvail[topLevelInd] = true;
}


void pHMCForce::resetOmegaForceStrengthPREC() {
  for (int I=0; I<pHMCForcePolynomialsMAX; I++) {
    omegaForceStrengthAnalysisPREC[I]->clearData();
  }
}


void pHMCForce::resetOmegaForceStrengthGLOBAL() {
  for (int I=0; I<pHMCForcePolynomialsMAX; I++) {
    omegaForceStrengthAnalysisGLOBAL[I]->clearData();
  }
}


void pHMCForce::calcAverageOmegaForceStrengthPREC(int polynomSlot) {
  omegaForceStrengthAnalysisPREC[polynomSlot]->getAverageVectorInflated(avgOmegaForces);
}


double* pHMCForce::getAverageOmegaForceStrengthPREC() {
  return avgOmegaForces;
}
  

void pHMCForce::setQuasiHermiteanMode(bool qHM) {
  quasiHermiteanMode = qHM;
  if (LogLevel>2) printf("Quasi-Hermitean-Mode set to %d on force %d on node %d.\n",quasiHermiteanMode,ID,ownNodeID);
}


void pHMCForce::setBMatFactorizationMode(bool bMatFac) {
  bMatFactorizationMode = bMatFac;
  if (LogLevel>2) printf("B-Matrix Factorization Mode set to %d on force %d on node %d.\n",bMatFactorizationMode,ID,ownNodeID);
}


void pHMCForce::setTuneMode(bool tuneM) {
  tuneMode = tuneM;
  if (LogLevel>2) printf("Tune-Mode set to %d on force %d on node %d.\n",tuneMode,ID,ownNodeID);  
}


int pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyDegree;
Complex* pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyRoots;
double pHMCForce_applyInverseSQRTPolynomial_HelperFunction_Nu;
double pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyLambda;
double pHMCForce_applyInverseSQRTPolynomial_HelperFunction(double x) {
  if ((x<0) || (x>pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyLambda)) {
    printf("ERROR in pHMCForce_applyInverseSQRTPolynomial_HelperFunction: Input larger than Poly-Lambda!\n");
    exit(0);
  }
  Complex z(pHMCForce_applyInverseSQRTPolynomial_HelperFunction_Nu, 0.0);
  Complex c(x, 0.0);
  
  for (int I=0; I<pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyDegree; I++) {
    z = z * (c - pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyRoots[I]);  
  }
  if (fabs(z.y) > 1E-11 * fabs(z.x)) {
    printf("ERROR in pHMCForce_applyInverseSQRTPolynomial_HelperFunction: Imaginary part not zero (%1.5e: %1.5e, %1.5e)!\n",x,z.x,z.y);
    exit(0);
  }
  if (z.x < 0) {
    printf("ERROR in pHMCForce_applyInverseSQRTPolynomial_HelperFunction: Real part not positive (%1.5e: %1.5e, %1.5e)!\n",x,z.x,z.y);
    exit(0);
  }
    
  return 1.0 / sqrt(z.x);  
}


void pHMCForce::setChebyApproxOfInverseSQRTofPolynomial(int polynomSlot) {
  if (polynomSlot != 0) return;
  delete chebyApproxOfInverseSQRTofPolynomial[polynomSlot];
  
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyDegree = approxpolyDegree[polynomSlot];
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyRoots = approxpolyRoots[polynomSlot];
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyLambda = approxpolyLambda[polynomSlot];
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_Nu = nu[polynomSlot];  
  
  double TOL = pHMCForceDirectOmegaSamplingRelAccuracy;
  chebyApproxOfInverseSQRTofPolynomial[polynomSlot] = new GeneralChebyshevApproximation();
  chebyApproxOfInverseSQRTofPolynomial[polynomSlot]->calcApproximation(&pHMCForce_applyInverseSQRTPolynomial_HelperFunction, 0.0, pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyLambda, TOL, 10000);
}


void pHMCForce::plotChebyApproxOfInversePolynomialSQRT(int polynomSlot, double eps, double lam, char* filename) {
  if (polynomSlot != 0) return;
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyDegree = approxpolyDegree[polynomSlot];
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyRoots = approxpolyRoots[polynomSlot];
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyLambda = approxpolyLambda[polynomSlot];
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_Nu = nu[polynomSlot];  

  chebyApproxOfInverseSQRTofPolynomial[polynomSlot]->plotPolynomial(filename, &pHMCForce_applyInverseSQRTPolynomial_HelperFunction, eps, lam, 10000);
}


void pHMCForce::applyDistributedInverseSQRTPolynomialCHEBYSHEV(DistributedMemoryObject* input, DistributedMemoryObject* output, DistributedMemoryObject* phi, int polynomSlot) {
  if (polynomSlot != 0) {
    printf("ERROR in pHMCForce::applyDistributedInverseSQRTPolynomialCHEBYSHEV: Only valid for polynomSlot 0\n");
    exit(0);
  }
  fermiOps->applyDistributedFermionMatrixMMDaggerChebyshevPolynomial(input, output, phi, chebyApproxOfInverseSQRTofPolynomial[polynomSlot], quasiHermiteanMode, true);
}


void pHMCForce::applyInverseSQRTPolynomialCHEBYSHEV(Complex* input, Complex* output, double* phi, int polynomSlot) {    
  if (polynomSlot != 0) {
    printf("ERROR in pHMCForce::applyInverseSQRTPolynomialCHEBYSHEV: Only valid for polynomSlot 0\n");
    exit(0);
  }
  fermiOps->applyFermionMatrixMMDaggerChebyshevPolynomial(input, output, phi, chebyApproxOfInverseSQRTofPolynomial[polynomSlot], quasiHermiteanMode, true);
}


void pHMCForce::applyInverseSQRTPolynomialKRYLOV(Complex* input, Complex* output, double* phi, int polynomSlot) {  
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyDegree = approxpolyDegree[polynomSlot];
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyRoots = approxpolyRoots[polynomSlot];
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_PolyLambda = approxpolyLambda[polynomSlot];
  pHMCForce_applyInverseSQRTPolynomial_HelperFunction_Nu = nu[polynomSlot];
  
  double TOL = pHMCForceDirectOmegaSamplingRelAccuracy;
  int neededIter = 0;  
 
  bool b = fermiOps->applyFermionMatrixMMDaggerFunction(input, output, phi, &pHMCForce_applyInverseSQRTPolynomial_HelperFunction, TOL, -1, neededIter, VectorCollectionCount, VectorCollection, quasiHermiteanMode, true);
  if (LogLevel>3) printf("ApplyInverseSQRTPolynomialKRYLOV: neededIter = %d\n",neededIter);  
  if (!b) {
    printf("ERROR in applyInverseSQRTPolynomial!\n");
    exit(0);
  }
}


void pHMCForce::applyFermionDoubleMatrixMMPolPol(Complex* input, Complex* output, double* phi, bool inFourierSpace) {
  int vectorLengthXtrSize = fermiOps->getVectorLengthXtrSize();
  Complex* p1 = input;
  Complex* p2 = auxVectors[2];
  
  if (quasiHermiteanMode) {
    fermiOps->executeFermionQuasiHermiteanMatrixFermionDaggerMatrixMultiplication(p1, p2, phi, true, inFourierSpace);
  } else {
    fermiOps->executeFermionMatrixFermionDaggerMatrixMultiplication(p1, p2, phi, 2, 1, inFourierSpace);
  }  
  p1 = p2;
  p2 = auxVectors[3];
  
  calcMonomialApplication(0, p1, p2, phi, 0, approxpolyDegree[0]-1, inFourierSpace);
  calcMonomialApplication(0, p2, output, phi, 0, approxpolyDegree[0]-1, inFourierSpace);
  
  double fac = nu[0] * nu[0];
  for (int I=0; I<vectorLengthXtrSize; I++) {
    output[I].x *= fac;
    output[I].y *= fac;
  }
}


bool pHMCForce::applyInverseFermionDoubleMatrixMMPolPol(Complex* input, Complex* output, double* phi, double TOL, int maxIter, int& neededIter) {
  int iterCount = 0;
  Complex errorOld;
  Complex errorNew;
  Complex alpha;
  Complex alphaZ;
  Complex alphaN;
  double beta;
  int vectorLengthXtrSize = fermiOps->getVectorLengthXtrSize();
  int oneDimSizeL0 = fermiOps->get1DSizeL0();
  int oneDimSizeL1 = fermiOps->get1DSizeL1();
  int oneDimSizeL2 = fermiOps->get1DSizeL2();
  int oneDimSizeL3 = fermiOps->get1DSizeL3();
  
  
  if ((VectorCollectionCount<4) || (fermiOps->isMultiThreadedOpsActivated())) {
    printf("ERROR in pHMCForce::applyInverseFermionDoubleMatrixMMPolPol: VectorCollectionCount must be 4 at least and multi-threaded mode must be switched off!!!\n");
    exit(0);
  }
  
  Complex* Inverse_b = VectorCollection[0];
  Complex* Inverse_rest = VectorCollection[1];
  Complex* Inverse_p = VectorCollection[2];
  Complex* Inverse_s = VectorCollection[3];

  fermiOps->performFFT(input, Inverse_b, ExtremeFFT4D_Forward);
  double normFac = 1.0 / (oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);
  for (int I=0; I<vectorLengthXtrSize; I++) {
    Inverse_b[I].x *= normFac;
    Inverse_b[I].y *= normFac;
  }

  //Startvector
  fermiOps->zeroFermionVector(output);
  SSE_ZCopy(vectorLengthXtrSize, Inverse_b, 1, Inverse_rest, 1);

  //Weitere Initialisierung
  SSE_ZCopy(vectorLengthXtrSize, Inverse_rest, 1, Inverse_p, 1);
  SSE_ComplexSquareNorm(oneDimSizeL0,oneDimSizeL1,oneDimSizeL2,oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, Inverse_rest, errorOld.x);
  if (LogLevel>5) errorOld.print();

  //Hauptschleife
  neededIter = 0;
  while ((sqrt(errorOld.x) > TOL) && ((neededIter<maxIter) || (maxIter<=0))) {
    neededIter++;
    iterCount++;
    
    applyFermionDoubleMatrixMMPolPol(Inverse_p, Inverse_s, phi, true);

    SSE_ComplexScalarProduct(oneDimSizeL0,oneDimSizeL1,oneDimSizeL2,oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Inverse_p, Inverse_rest, alphaZ);
    SSE_ComplexScalarProduct(oneDimSizeL0,oneDimSizeL1,oneDimSizeL2,oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, Inverse_p, Inverse_s, alphaN);
    
    alpha = alphaZ / alphaN;
    
    SSE_ComplexVectorAddition(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, alpha, Inverse_p, output);

    errorNew.y=0;

    alpha.x *= -1.0;
    alpha.y *= -1.0;
    SSE_ComplexVectorAdditionWithSquaredNorm(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, alpha, Inverse_s, Inverse_rest, errorNew.x);

    beta = errorNew.x / errorOld.x;    
    SSE_ComplexVectorAdditionSPECIAL1(oneDimSizeL0, oneDimSizeL1, oneDimSizeL2, oneDimSizeL3, xtraSize1,xtraSize2,xtraSize3, beta, Inverse_rest, Inverse_p);

    errorOld.x = errorNew.x;
    if (LogLevel>5) printf("%1.15f + %1.15fi\n",errorOld.x,errorOld.y);
  }

  fermiOps->performFFT(output, Inverse_p, ExtremeFFT4D_Backward);
  SSE_ZCopy(vectorLengthXtrSize, Inverse_p, 1, output, 1);

  if (LogLevel>3) printf("... ready after %d recursions with error %1.15f.\n",iterCount,sqrt(errorOld.x));
  if ((LogLevel>0) && (iterCount>oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3)) printf("ALARM: More iterations than matrix size for inverse calculation (%d)!!!\n",iterCount);
  return (sqrt(errorOld.x) <= TOL);
}


double pHMCForce::calcGaussianWeightFactor(double* phi, double TOL, int &Ncg, int &NmmdagApplications) {
  int oneDimSizeL0 = fermiOps->get1DSizeL0();
  int oneDimSizeL1 = fermiOps->get1DSizeL1();
  int oneDimSizeL2 = fermiOps->get1DSizeL2();
  int oneDimSizeL3 = fermiOps->get1DSizeL3();
  double res = NaN;

  if ((VectorCollectionCount<5) || (fermiOps->isMultiThreadedOpsActivated())) {
    printf("ERROR in pHMCForce::calcGaussianWeightFactor: VectorCollectionCount must be 5 at least and multi-threaded mode must be switched off!!!\n");
    exit(0);
  }

  Complex* output = VectorCollection[4];
  
  bool b = applyInverseFermionDoubleMatrixMMPolPol(omegaField, output, phi, TOL, -1, Ncg);
  
  NmmdagApplications = Ncg * (1+2*approxpolyDegree[0]);
  Complex s1(NaN, NaN);
  Complex s2(NaN, NaN);
  
  SSE_ComplexScalarProduct(oneDimSizeL0,oneDimSizeL1,oneDimSizeL2,oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, omegaField, output, s1);
  SSE_ComplexScalarProduct(oneDimSizeL0,oneDimSizeL1,oneDimSizeL2,oneDimSizeL3, xtraSize1, xtraSize2, xtraSize3, omegaField, omegaField, s2);
  res = exp(-0.5*(s1.x - s2.x));
  
  //Only for Testing...
  if (fabs(s1.y/s1.x)>1E-10) {
    printf("Error in pHMCForce::calcGaussianWeightFactor: Non-vanishing imaginary part!!!\n");
    exit(0);
  }
  
/*  Complex* output2 = VectorCollection[5];
  applyFermionDoubleMatrixMMPolPol(output, output2, phi, false);  
  fermiOps->transformFromXtraSizeArray(output, output);
  fermiOps->transformFromXtraSizeArray(output2, output2);
  double diff = 0;
  for (int I=0; I<oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3; I++) {
    diff += norm(omegaField[I]-output2[I]);
  }
  diff /= (oneDimSizeL0*oneDimSizeL1*oneDimSizeL2*oneDimSizeL3);*/
  
  if (!b) res = NaN;
  return sqrt(res);
}
