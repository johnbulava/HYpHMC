#include "HMCPropagator.h"

HMCPropagator::HMCPropagator(FermionMatrixOperations* fOps, double lam, double kap, int nf , double gam): Propagator(fOps, lam, kap, 0, 0, 0, 0, 0, 0, 0, nf, gam, isNaN(lam), 0.0) {
  if (LogLevel>2) printf("Initializing HMC-Propagator with lambda = %1.3f, kappa = %1.3f, Nf = %d, gamma = %1.3f\n", lam, kap, nf, gam);
  HMCforces = NULL;
  ini(fOps, lam, kap, 0, 0, 0, 0, 0, 0, 0, nf, gam, isNaN(lam), 0);
}


HMCPropagator::~HMCPropagator() {
  desini();
}


void HMCPropagator::iniAdditionalFields() {
  GaussVector = fermiOps->createFermionVector();
}


void HMCPropagator::desiniAdditionalFields() {
  fermiOps->destroyFermionVector(GaussVector);
  GaussVector = NULL;
  delete[] HMCforces;
}


int HMCPropagator::getFermionForceSubCategoryCount() {
  return 1;
}


int HMCPropagator::getFermionForceSubCategory(double x) {
  return 0;
}


void HMCPropagator::ThreadController(int nr, int mode, int& para_int, double& para_double, double* data) {
  if (LogLevel>3) printf("Thread-Controller called with nr = %d, mode = %d with paraInt = %d, paraDouble = %1.15f\n", nr, mode, para_int, para_double);

  if (mode == HMCPropagator_SAMPLE) {  
    HMCforces[nr]->sampleOmegaField(phiField, GaussVector);
  }

  if (mode == Propagator_FORCE) {
    double TOL = para_double;
    HMCforces[nr]->calcPhiForce(phiField, TOL); 
    HMCforces[nr]->calcOmegaMMdaggerInverseOmegaScalarProduct();
    para_double = HMCforces[nr]->getOmegaMMdaggerInverseOmegaScalarProduct();
  }

  if (mode == HMCPropagator_INVERSION) {
    double TOL = para_double;
    HMCforces[nr]->calcInverse(phiField, TOL);
    HMCforces[nr]->calcOmegaMMdaggerInverseOmegaScalarProduct();
    para_double = HMCforces[nr]->getOmegaMMdaggerInverseOmegaScalarProduct();
  }

  para_int = 211;
}


void HMCPropagator::SlaveController() {
  #ifdef UseMPI
  int VL = fermiOps->getVectorLength();
  int nr, mode, para_int;
  double para_double;
  int precUseInt;
  double precM, precS;
  bool precUse;
  
  if (LogLevel>3) printf("Entering Slave Control Loop...\n");
  while (true) {
    MPI_Status* status = new MPI_Status;
    MPI_Recv(&nr, 1 , MPI_INT , 0, 0, MPI_COMM_WORLD, status);
    MPI_Recv(&mode, 1 , MPI_INT , 0, 0, MPI_COMM_WORLD, status);
    MPI_Recv(&para_int, 1 , MPI_INT , 0, 0, MPI_COMM_WORLD, status);
    MPI_Recv(&para_double, 1 , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, status);
    MPI_Recv(&precUseInt, 1 , MPI_INT , 0, 0, MPI_COMM_WORLD, status);
    MPI_Recv(&precM, 1 , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, status);
    MPI_Recv(&precS, 1 , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, status);
    MPI_Recv(phiField, VL/2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, status);
    if (mode == HMCPropagator_SAMPLE) {
      MPI_Recv(GaussVector, 2*VL, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, status);
    }
    delete status;
    precUse = (bool) precUseInt;
    fermiOps->setPreconditioner(precUse, precM, precS);
//printf("               SLAVE: RECEVIED...\n");    
//printf("s");    
   ThreadController(nr, mode, para_int, para_double, NULL);
//printf("%d",ownNodeID);       
//printf("               SLAVE: CALUCULATION READY...\n");    
    
    
    MPI_Send(&nr, 1 , MPI_INT , 0, 0, MPI_COMM_WORLD);
    MPI_Send(&mode, 1 , MPI_INT , 0, 0, MPI_COMM_WORLD);
    MPI_Send(&para_int, 1 , MPI_INT , 0, 0, MPI_COMM_WORLD);
    MPI_Send(&para_double, 1 ,  MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    double* dSdPhi = (double*) HMCforces[nr]->getdSdPhi();
    MPI_Send(dSdPhi, VL/2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//printf("               SLAVE: DATA SENT...\n");    
    
    
    if (mode == Propagator_EXIT) break;
  }
  #endif 
}


void HMCPropagator::threadedExecute(int mode, double TOL) {
  if (forceCount <= 0) return;
  int* numbers = new int[forceCount];
  int* modes = new int[forceCount];
  int* para_int = new int[forceCount];
  double* para_double = new double[forceCount];
  #ifdef UseMPI
  MPI_Request* requests = new MPI_Request[8*forceCount];
  #endif 
  
  int I;
  for (I=0; I<forceCount; I++) {
    numbers[I] = I;
    modes[I] = mode;
    para_int[I] = 0;
    para_double[I] = TOL;
  }

  double precM;
  double precS;
  bool precUse;
  fermiOps->getPreconditionerParameter(precUse, precM, precS);

  int execForces = 0;
  while (execForces<forceCount) {
    int rem = forceCount-execForces;
    if (rem>nodeCount) rem = nodeCount;
    
    
//printf("Main node SEND...\n");    
    for (I=execForces+1; I<execForces+rem; I++) {
    #ifdef UseMPI
      int VL = fermiOps->getVectorLength();
      int precUseInt = (int) precUse;
      int remoteNode = I % nodeCount;
      if (LogLevel>3) printf("Remote call for nr = %d on node %d (mode = %d)...\n",I,remoteNode, mode);
      
      MPI_Send(&(numbers[I]), 1, MPI_INT, remoteNode, 0, MPI_COMM_WORLD /*, &(requests[8*I+0])*/);
      MPI_Send(&(modes[I]), 1, MPI_INT, remoteNode, 0, MPI_COMM_WORLD /*, &(requests[8*I+1])*/);
      MPI_Send(&(para_int[I]), 1, MPI_INT, remoteNode, 0, MPI_COMM_WORLD /*, &(requests[8*I+2])*/);
      MPI_Send(&(para_double[I]), 1, MPI_DOUBLE, remoteNode, 0, MPI_COMM_WORLD /*, &(requests[8*I+3])*/);
      MPI_Send(&precUseInt, 1, MPI_INT, remoteNode, 0, MPI_COMM_WORLD /*, &(requests[8*I+4])*/);
      MPI_Send(&precM, 1, MPI_DOUBLE, remoteNode, 0, MPI_COMM_WORLD /*, &(requests[8*I+5])*/);
      MPI_Send(&precS, 1, MPI_DOUBLE, remoteNode, 0, MPI_COMM_WORLD /*, &(requests[8*I+6])*/);
      MPI_Send(phiField, VL/2, MPI_DOUBLE_PRECISION, remoteNode, 0, MPI_COMM_WORLD /*, &(requests[8*I+7])*/);
      
      if (mode == HMCPropagator_SAMPLE) {
        sampleGaussVector();
	//NICHT Isend - muß hier warten
	MPI_Send(GaussVector, 2*VL, MPI_DOUBLE_PRECISION, remoteNode, 0, MPI_COMM_WORLD);
      }
    #endif 
    }
    
//printf("Main node SEND READY...\n");    
    
    //Einen Thread selber ausfuehren
    if (LogLevel>3) printf("Local call for nr = %d (mode = %d)...\n",execForces,mode);
    if (mode == HMCPropagator_SAMPLE) {
      sampleGaussVector();
    }
//printf("Main node OWN CALCULATION...\n");    
//printf("m");    
    
    ThreadController(numbers[execForces], modes[execForces], para_int[execForces], para_double[execForces], NULL);

//printf("Main node CALCULATION READY, RECEIVING...\n");    
//printf("M");    

    //Read results from processes
    for (I=execForces+1; I<execForces+rem; I++) {
    #ifdef UseMPI
      int VL = fermiOps->getVectorLength();
      int remoteNode = I % nodeCount;
      MPI_Status* status = new MPI_Status;

      MPI_Recv(&(numbers[I]), 1 , MPI_INT , remoteNode, 0, MPI_COMM_WORLD, status);
      MPI_Recv(&(modes[I]), 1 , MPI_INT , remoteNode, 0, MPI_COMM_WORLD, status);
      MPI_Recv(&(para_int[I]), 1 , MPI_INT , remoteNode, 0, MPI_COMM_WORLD, status);
      MPI_Recv(&(para_double[I]), 1 , MPI_DOUBLE, remoteNode, 0, MPI_COMM_WORLD, status);
      HMCforces[I]->setOmegaMMdaggerInverseOmegaScalarProduct(para_double[I]);
      
      double* dSdPhi = (double*) HMCforces[I]->getdSdPhi();
      MPI_Recv(dSdPhi, VL/2, MPI_DOUBLE, remoteNode, 0, MPI_COMM_WORLD, status);

      if (para_int[I] != 211) {
        printf("ERROR in MPI - Connection!!!\n");
	exit(0);
      }

      delete status;
      if (LogLevel>3) printf("Remote call for nr = %d on node %d (mode = %d) has been processed!\n",I,remoteNode, mode);
    #endif 
    }
//printf("Main node RECEIVING READY...\n");    
  
    execForces += rem;
  }
  
  #ifdef UseMPI
  delete[] requests;
  #endif 
  delete[] numbers;
  delete[] modes;
  delete[] para_int;
  delete[] para_double;
}

  
void HMCPropagator::setNf(int nf) {
  Nf = nf;
  desiniForceCalculators();
  delete[] HMCforces;
  forceCount = nf / 2;
  if (fermiOps->getYN()<=0) forceCount = 0;

  if (forceCount>0) {
    HMCforces = new HMCForce*[forceCount];
    forces = new Force*[forceCount];    
    int I;
    for (I=0; I<forceCount; I++) {
      bool local = false;
      if ((I%nodeCount) == ownNodeID) local = true;
      HMCforces[I] = new HMCForce(fermiOps, local);
      forces[I] = HMCforces[I];
    }
  }
}


void HMCPropagator::sampleGaussVector() {
  if (LogLevel>3) printf("Sampling Gauss Vector on node = %d.\n", ownNodeID);
  fermiOps->fillGaussRandomVector(GaussVector,-1);
}


void HMCPropagator::sampleOmegaFields() {
  if (LogLevel>3) printf("Sampling all %d omega fields directly.\n", forceCount);
  threadedExecute(HMCPropagator_SAMPLE,0);
  Sold = NaN;
  Sact = NaN;  
  SbeforeProp = NaN;
  SafterProp = NaN;
}


double HMCPropagator::calcOmegaAction(bool outMMdaggerInverseOmega_READY, double finalTOL) {
  double S = 0;

  if (fermiOps->getYN()>0) {
    if (!outMMdaggerInverseOmega_READY) {
      threadedExecute(HMCPropagator_INVERSION, finalTOL);
    } 
    S = 0;
    int I;
    for (I=0; I<forceCount; I++) {
      S = S + 0.5 * HMCforces[I]->getOmegaMMdaggerInverseOmegaScalarProduct();
    }
  } else {
    S = 0;
  }
  return S;
}


double HMCPropagator::calcTotalAction(double finalTOL) {
  double S = calcOmegaAction(false, finalTOL);
  S = S + calcPhiMomentumAction();
  S = S + calcPhiAction();
  return S;
}


bool HMCPropagator::LeapOmelyanMarkovStep(int iterations, double epsilon, double propTOL, double finalTOL, double lambda, double rho, double theta, double mu) {
  int I,I2;  
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double epsilonHalf = 0.5 * epsilon;
  double epsilonLambda = lambda * epsilon;
  double epsilonRho = rho * epsilon;
  double epsilonMu = mu * epsilon;
  double epsilonTheta = theta * epsilon;
  double epsilonHalfOneMinusTwoLambdaPlusTheta = 0.5*(1.0 - 2.0*(lambda + theta))*epsilon;
  double epsilonOneMinusTwoMuPlusRho = (1.0 - 2.0*(mu + rho))*epsilon;
  
  double twoEpsilonLambda = 2.0 * lambda * epsilon;
  double epsilonOneMinusTwoLambda = (1.0 - 2.0*lambda)*epsilon;
  double* dSdPhi = new double[4*L0*L1*L2*L3];
  double* dSdPhi_Dummy = NULL;


  if (!(Sold == Sold)) {
    //Sold nicht vorhanden...
    calcPhiDerivativesOfSPhi(dSdPhi);
    Sold = calcPhiMomentumAction() + calcPhiAction();
    if (fermiOps->getYN()>0) {
      threadedExecute(Propagator_FORCE, finalTOL);
      Sold = Sold + calcOmegaAction(true, finalTOL);
    
      for (I=0; I<forceCount; I++) {
        dSdPhi_Dummy = (double*) HMCforces[I]->getdSdPhi();
	for (I2=0; I2<4*L0*L1*L2*L3; I2++) dSdPhi[I2] += 0.5 * dSdPhi_Dummy[I2];
      }
    }
  } else {
    //Sold ist vorhanden...
    calcFullPhiDerivatives(dSdPhi, propTOL, 0);
  }
  savePhiFields();
  MiniPhiMomentumStep(dSdPhi, epsilonLambda);
  
  if (LogLevel>2) printf("Old Action is %f\n", Sold);
  SbeforeProp = Sold;

  for (I2=0; I2<iterations; I2++) {
    if ((lambda==0.5) && (theta==0) && (rho==0) && (mu==0)) {
      //Normaler Leap-Frog Schritt
//      if (LogLevel>2) printf("Leap-Frog Integrator\n");      
      MiniPhiStep(epsilon);
    } else if ((theta==0) && (rho==0) && (mu==0)) {
      //Omelyan - Schritt mit Order 2
//      if (LogLevel>2) printf("Omelyan Integrator with order 2\n");      
      MiniPhiStep(epsilonHalf);

      calcFullPhiDerivatives(dSdPhi, propTOL, 0);
      MiniPhiMomentumStep(dSdPhi,epsilonOneMinusTwoLambda);

      MiniPhiStep(epsilonHalf);

    } else {
      //Omelyan - Schritt mit Order 4
//      if (LogLevel>2) printf("Omelyan Integrator with order 4\n");      
      MiniPhiStep(epsilonRho);

      calcFullPhiDerivatives(dSdPhi, propTOL, 0);
      MiniPhiMomentumStep(dSdPhi, epsilonTheta);

      MiniPhiStep(epsilonMu);

      calcFullPhiDerivatives(dSdPhi, propTOL, 0);
      MiniPhiMomentumStep(dSdPhi, epsilonHalfOneMinusTwoLambdaPlusTheta);

      MiniPhiStep(epsilonOneMinusTwoMuPlusRho);

      calcFullPhiDerivatives(dSdPhi, propTOL, 0);
      MiniPhiMomentumStep(dSdPhi, epsilonHalfOneMinusTwoLambdaPlusTheta);

      MiniPhiStep(epsilonMu);
      

      calcFullPhiDerivatives(dSdPhi, propTOL, 0);
      MiniPhiMomentumStep(dSdPhi, epsilonTheta);

      MiniPhiStep(epsilonRho);
    }

    calcFullPhiDerivatives(dSdPhi, propTOL, 0);
    MiniPhiMomentumStep(dSdPhi, twoEpsilonLambda);
  }
  //Bring momenta on integer step numbers again
  MiniPhiMomentumStep(dSdPhi, -epsilonLambda);
  
  Sact = calcTotalAction(finalTOL);
  if (LogLevel>2) printf("New Action is %f\n", Sact);
  SafterProp = Sact;
  
  double deltaS = Sact-Sold;
  bool accepted = true;
  if (deltaS>0) {
    if (AdvancedZufall(AdvancedSeed) > exp(-deltaS)) {
      restorePhiFields();
      accepted = false;
      if (LogLevel>2) printf("Rejecting configuration!!!\n");
    }
  }
  
  if (accepted) {
    Sold = Sact;
  } else {
    Sact = Sold;
  }

  delete[] dSdPhi;
  return accepted;
}


void HMCPropagator::improvePreconditioningParameters(int iterGrob, int iterFein, double testTOL) {
  bool PrecUse;
  double PrecMold, PrecSold;
  fermiOps->getPreconditionerParameter(PrecUse, PrecMold, PrecSold);  

  if (!PrecUse) return;
  if (forceCount<=0) return;
  if (fermiOps->getYN()<=0) return;
  
  int bestIter;
  double bestM = PrecMold;
  double bestS = PrecSold;
  Complex* omegaField = HMCforces[0]->getOmega();
  Complex* outMMdaggerInverseOmega = HMCforces[0]->getMMdaggerInverseOmega();
  if (!(fermiOps->solveFermionMatrixLGS(omegaField, outMMdaggerInverseOmega, (double*) phiField, testTOL, true, false, -1, bestIter))) return;


  double avgPhi, avgStagPhi, avgNorm, sigmaNorm;
  measure(avgPhi, avgStagPhi, avgNorm, sigmaNorm);
  
  if (LogLevel>1) printf("Preconditioning-Improver at avgPhi = %1.3f, avgStagPhi = %1.3f: oldM = %1.3f, oldS = %1.3f needs %d iterations.\n", avgPhi, avgStagPhi, PrecMold, PrecSold, bestIter);
  double* m = new double[iterGrob+iterFein];
  double* s = new double[iterGrob+iterFein];
  int I;
  
  for (I=0; I<iterGrob; I++) {
    m[I] = 2*AdvancedZufall(AdvancedSeed)*avgPhi;
    s[I] = 2*AdvancedZufall(AdvancedSeed)*avgStagPhi;
  }
  
  for (I=0; I<iterGrob; I++) {
    fermiOps->setPreconditioner(true, m[I], s[I]);
    int newIter = 0;
    if (fermiOps->solveFermionMatrixLGS(omegaField, outMMdaggerInverseOmega, (double*) phiField, testTOL, true, false, bestIter, newIter)) {
      if (newIter < bestIter) {
        bestIter = newIter;      
        bestM = m[I];
        bestS = s[I];
        if (LogLevel>2) printf("Preconditioning improved: With PrecM = %1.3f, PrecS = %1.3f only %d iterations needed.\n", bestM, bestS, bestIter);
      }    
    }
  }
  
  for (I=0; I<iterFein; I++) {
    m[I] = abs(bestM + 2*(AdvancedZufall(AdvancedSeed)-0.5)*0.1*avgPhi);
    s[I] = abs(bestS + 2*(AdvancedZufall(AdvancedSeed)-0.5)*0.1*avgStagPhi);
  }

  for (I=0; I<iterFein; I++) {
    fermiOps->setPreconditioner(true, m[I], s[I]);
    int newIter = 0;
    if (fermiOps->solveFermionMatrixLGS(omegaField, outMMdaggerInverseOmega, (double*) phiField, testTOL, true, false, bestIter, newIter)) {
      if (newIter < bestIter) {
        bestIter = newIter;      
        bestM = m[I];
        bestS = s[I];
        if (LogLevel>2) printf("Preconditioning improved: With PrecM = %1.3f, PrecS = %1.3f only %d iterations needed.\n", bestM, bestS, bestIter);
      }    
    }
  }

  fermiOps->setPreconditioner(true, bestM, bestS);
  double PrecMact, PrecSact;
  fermiOps->getPreconditionerParameter(PrecUse, PrecMact, PrecSact);
  if (LogLevel>1) printf("Best Preconditioning data: PrecM = %1.3f, PrecS = %1.3f with %d iterations needed.\n", PrecMact, PrecSact, bestIter); 
  
  delete[] m;
  delete[] s;
}
