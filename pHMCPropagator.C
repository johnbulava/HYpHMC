#include "pHMCPropagator.h"
#include "Tools.h"

pHMCPropagator::pHMCPropagator(FermionMatrixOperations* fOps, double lam, double kap, double current, double c6, double c8, double c10, double lam6, double lam8, double lam10, int nf, double gam, bool sphMode, double sphZeta, double tht, int subPolCnt, double* polEps, double* polLam, int* polDeg, int precMCnt, double* precMss, int digit, double alpha, int maxPolDegPerNod, int addAuxVecsCount):Propagator(fOps, lam, kap, current, c6, c8, c10, lam6, lam8, lam10, nf, gam, sphMode, sphZeta) {
  fileNameIdentifier = NULL;
  if (LogLevel>2) printf("Initializing pHMC-Propagator with lambda = %1.3f, kappa = %1.3f, current = %1.3e, c6 = %1.3e, c8 = %1.3e, c10 = %1.3e, lam6 = %1.3e, lam8 = %1.3e, lam10 = %1.3e, Nf = %d, gamma = %1.3f, SphericalMode = %d, SphericalZeta = %1.3e, theta = %1.3f, SubPolCount = %d, PrecMassCount = %d, Poly.Digit = %d, Poly.Alpha = %1.3f, max.Poly.Deg. Per Node = %d, addAuxVecsCount=%d\n", lam, kap, current, c6, c8, c10, lam6, lam8, lam10, nf, gam, sphMode, sphZeta, tht, subPolCnt, precMCnt, digit, alpha, maxPolDegPerNod, addAuxVecsCount);
  int I,I2;
  for (I=0; I<precMCnt; I++) {
    if (LogLevel>2) printf(" ...with precMass %d = %1.3f\n",I,precMss[I]);
  }
  
  for (I=0; I<1+subPolCnt; I++) {
    if (LogLevel>2) printf(" ...with sub Polynomial %d: Poly.Eps. = %1.3f, Poly.Lam. = %1.3f, Poly.Deg. = %d\n",I,polEps[I], polLam[I], polDeg[I]);
  }

  pHMCforces = NULL;
  AdditionalAuxVectorsCount = addAuxVecsCount;
  subPolyCount = subPolCnt;
  precMassCount = precMCnt;
  polyRoots = new Complex**[1+precMassCount];
  approxPoly = new PolynomialApproximation*[1+precMassCount];
  polyEpsilon = new double[1+subPolyCount];
  polyLambda = new double[1+subPolyCount];
  polyDegree = new int[1+subPolyCount];
  polyDigit = digit;
  polyAlpha = alpha;
  maxPolyDegreePerNode = maxPolDegPerNod;
  if ((maxPolyDegreePerNode % 2) == 1) {
    printf("Error: Maximal Polynomial Degree per Node must be even!!!\n");
    exit(0);
  }
  for (I=0; I<1+subPolyCount; I++) {
    polyEpsilon[I] = polEps[I];
    polyLambda[I] = polLam[I];
    polyDegree[I] = polDeg[I];
    if ((maxPolyDegreePerNode<polyDegree[I]) && (subPolyCount>0)) {
      printf("ERROR: Sub Polynomials and Monomial-Distribution not compatible!!!\n");
      exit(0);
    }
  }
  for (I=0; I<subPolyCount; I++) {
    if (polyDegree[I+1]>polyDegree[I]) {
      printf("ERROR: Sub Polynomial degree must not increase!!!\n");
      exit(0);
    }
  }
  precMasses = new double[precMassCount];  
  for (I=0; I<precMassCount; I++) {
    precMasses[I] = precMss[I];
  }
  for (I=0; I<1+precMassCount; I++) {
    approxPoly[I] = NULL;
    polyRoots[I] = new Complex*[1+2*subPolyCount];
    for (I2=0; I2<1+2*subPolyCount; I2++) {
      polyRoots[I][I2] = NULL;
    }
  }
  globalForceTheta = tht;

  upperEWboundLogCount = 0;
  upperEWboundLog = new double*[pHMCPropagator_UpperEWboundLogMAX];
  for (int I=0; I<pHMCPropagator_UpperEWboundLogMAX; I++) {
    upperEWboundLog[I] = new double[4];
  }

  syncRandom = NaN;
  nodesReady = false;
  quasiHermiteanMode = true;  
  
  setPolyData();
  ini(fOps, lam, kap, current, c6, c8, c10, lam6, lam8, lam10, nf, gam, sphMode, sphZeta);
}


pHMCPropagator::~pHMCPropagator() {
  desini();
  
  for (int I=0; I<1+precMassCount; I++) {
    delete approxPoly[I];
    for (int I2=0; I2<1+2*subPolyCount; I2++) {
      delete[] polyRoots[I][I2];
    }
    delete[] polyRoots[I];
  }
  upperEWboundLogCount = 0;
  for (int I=0; I<pHMCPropagator_UpperEWboundLogMAX; I++) {
    delete[] upperEWboundLog[I];
  }
  delete[] upperEWboundLog;

  delete[] polyRoots;
  delete[] approxPoly;
  delete[] polyEpsilon;
  delete[] polyLambda;
  delete[] polyDegree;
  delete[] precMasses;
  delete[] fileNameIdentifier;
}

void pHMCPropagator::setFileNameIdentifier(char* identifier) {
  int n = strlen(identifier)+1;
  fileNameIdentifier = new char[n];
  strncpy(fileNameIdentifier, identifier, n);
}

void pHMCPropagator::iniAdditionalFields() {
}


void pHMCPropagator::desiniAdditionalFields() {
  delete[] pHMCforces;
}


int pHMCPropagator::getFermionForceSubCategoryCount() {
   return 1+2*subPolyCount;
}


int pHMCPropagator::getFermionForceSubCategory(double x) {
  return roundToInt(fabs(x));
}


void pHMCPropagator::ThreadController(int nr, int mode, int& para_int, double& para_double, double* data) {
  if (LogLevel>4) printf("Thread-Controller called with nr = %d, mode = %d with paraInt = %d, paraDouble = %1.15f\n", nr, mode, para_int, para_double);

  if (mode == pHMCPropagator_SETPOLYROOTS) {  
    pHMCforces[nr]->setApproxPolyRoots(roundToInt(data[0]), (Complex*) &(data[2]), para_int, para_double, data[1]);
    char* filename = new char[1000];
		printf ("Printing db dir\n");
		printf("Db dir = %s\n", DataBaseDirectory);
		printf("Done Printing it\n");
    snprintf(filename,1000,"%s/data/results/pHMC/miscellaneous/partPolynom_SubPoly%d_Force%d_Node%d.dat",DataBaseDirectory,roundToInt(data[0]),nr,ownNodeID);
		printf("filename = %s\n", filename);
    pHMCforces[nr]->plotApproxPoly(roundToInt(data[0]), 0, data[1], filename);
    snprintf(filename,1000,"%s/data/results/pHMC/miscellaneous/partPolynom_Roots_SubPoly%d_Force%d_Node%d.dat",DataBaseDirectory,roundToInt(data[0]), nr,ownNodeID);    
		printf("filename = %s\n", filename);
    pHMCforces[nr]->plotApproxPolyRoots(roundToInt(data[0]), filename);    
    snprintf(filename,1000,"%s/data/results/pHMC/miscellaneous/chebyApproxOfInversePolynomialSQRT_SubPoly%d_Force%d_Node%d.dat",DataBaseDirectory,roundToInt(data[0]), nr,ownNodeID);    
		printf("filename = %s\n", filename);
    pHMCforces[nr]->plotChebyApproxOfInversePolynomialSQRT(roundToInt(data[0]), 0, data[1], filename);    
    delete[] filename;
  }

  if (mode == pHMCPropagator_SETRANDOMSEED) {  
    AdvancedSeed = -para_int;
    if (LogLevel > 2) printf("Initializing Random Seed on node %d with %d\n",ownNodeID,AdvancedSeed);
    AdvancedZufall(AdvancedSeed); 
    double ran = AdvancedZufall(AdvancedSeed); 
    if (LogLevel > 2) printf("Random Generator on node %d initialized. Current Seed = %d. First Random Nr = %f\n",ownNodeID,AdvancedSeed,ran);    
  }
  
  if (mode == pHMCPropagator_DRAWRANDOMNUMBERSONSLAVENODES) {  
    int I;
    if (ownNodeID>0) {
      for (I=0; I<para_int; I++) {
        AdvancedZufall(AdvancedSeed); 
      }
    }
  }

  if (mode == pHMCPropagator_SETQUASIHERMITEANMODE) {  
    for (int I=0; I<forceCount; I++) {    
      pHMCforces[I]->setQuasiHermiteanMode((bool) para_int);
    }
    quasiHermiteanMode = (bool) para_int;
    if (LogLevel>2) printf("Quasi-Hermitean-Mode set to %d on Propagator on node %d.\n",quasiHermiteanMode,ownNodeID);    
  }
  
  if (mode == pHMCPropagator_SETTUNEMODE) {  
    for (int I=0; I<forceCount; I++) {    
      pHMCforces[I]->setTuneMode((bool) para_int);
    }
    tuneMode = (bool) para_int;
    if (LogLevel>2) printf("Tune-Mode set to %d on Propagator on node %d.\n",tuneMode,ownNodeID);    
  }
  
  if (mode == pHMCPropagator_SETMODELSELECTION) {  
  }

  if (mode == pHMCPropagator_GETAVERAGEPRECSUBPOL0OMEGAFORCESTRENGTH) {  
    pHMCforces[nr]->calcAverageOmegaForceStrengthPREC(0);  
  }
 
  if (mode == pHMCPropagator_RESETAVERAGEOMEGAFORCESTRENGTH) {  
    pHMCforces[nr]->resetOmegaForceStrengthPREC();
    pHMCforces[nr]->resetOmegaForceStrengthGLOBAL();
  }
 
  if (mode == pHMCPropagator_SYNCQPRECDATA) {  
    fermiOps->setQPreconditioner((bool) para_int, data[0], data[1]); 
  }

  if (mode == pHMCPropagator_SYNCRPRECDATA) {  
    fermiOps->setRPreconditioner((bool) para_int, data[0], data[1]); 
  }

  if (mode == pHMCPropagator_CLEARACTSTORE) {  
    pHMCforces[nr]->clearActStore();
  }
  
  if (mode == pHMCPropagator_CLEARTRASTARTSTORE) {  
    pHMCforces[nr]->clearTraStartStore();
  }
  
  if (mode == pHMCPropagator_COPYTRASTARTSTORETOACTSTORE) {  
    pHMCforces[nr]->copyTraStartStoreToActStore();
  }
  
  if (mode == pHMCPropagator_COPYACTSTORETOTRASTARTSTORE) {  
    pHMCforces[nr]->copyActStoreToTraStartStore();
  }
  
  if (mode == pHMCPropagator_STORECURRENTTOPLEVELFORCES) {  
    pHMCforces[nr]->storeCurrentTopLevelForcesToTraStartStore();
  }
  
  if (mode == pHMCPropagator_ACTIVATEFORCESTORING) {  
    int* dataInt = new int[para_int];
    for (int I=0; I<para_int; I++) dataInt[I] = roundToInt(data[I]);
    pHMCforces[nr]->activateForceStoring(para_int, dataInt);
    delete[] dataInt;
  }
  
  if (mode == pHMCPropagator_CHECKRANDOMSYNCHRON) {  
    para_double = AdvancedZufall(AdvancedSeed); 
  }
  
  if (mode == pHMCPropagator_SETPHI) { 
    if (data != NULL) {
      int VL = fermiOps->getVectorLength();
      SSE_ZCopy(VL/4, (Complex*) data, 1, (Complex*) phiField, 1);      
    }
  }
  
  if (mode == pHMCPropagator_MULTIPLYPHIMOMENTA) { 
    multiplyPhiMomenta(para_double);
  }
  
  if (mode == pHMCPropagator_MULTIPLYOMEGAMOMENTA) {  
    pHMCforces[nr]->multiplyOmegaMomenta(para_double);
  }

  if (mode == Propagator_SETFOURIERACCELERATIONTYPE) { 
    phiForceFourierType = para_int;
    phiForceFourierPara = para_double;
  }    
  
  if (mode == Propagator_SETFOURIERACCELERATIONMASSES) { 
    if (data != NULL) {
      int VL = fermiOps->getVectorLength();
      SSE_ZCopy(VL/16, (Complex*) data, 1, (Complex*) MomentaMasses, 1);      
    }
    momentumMassesDetermined = true;
  }  
  
  if (mode == pHMCPropagator_RESTOREPROP) { 
    restorePhiFields();    
  }
  
  if (mode == pHMCPropagator_RESTOREFORCE) { 
    pHMCforces[nr]->restoreOmegaField();    
  }
  
  if (mode == pHMCPropagator_SAMPLEPHIMOMENTA) { 
    samplePhiMomentumField();
  }
  
  if (mode == pHMCPropagator_SAMPLEOMEGAMOMENTA) { 
    if (para_int == nr) {
      pHMCforces[nr]->sampleOmegaMomenta();    
    } else {
      bool loc = pHMCforces[para_int]->isLocal();
      if (!loc) {
        if (nr<nodeCount) {
          pHMCforces[nr]->drawAndWasteOmegaMomentaRandomNumbers();
	}
      }
    }
  }
  
  if (mode == pHMCPropagator_SAMPLEOMEGAFIELDS) { 
    if (para_int == nr) {
      pHMCforces[nr]->sampleOmegaFields((double*) phiField); 
    } else {
      bool loc = pHMCforces[para_int]->isLocal();
      if (!loc) {
        if (nr<nodeCount) {
          pHMCforces[nr]->drawAndWasteOmegaMomentaRandomNumbers();
	}
      }
    }
  }
  
  if (mode == pHMCPropagator_PROPAGATEPHI) { 
    MiniPhiStep(para_double);  
  }
 
  if (mode == pHMCPropagator_PROPAGATEOMEGA) {     
    pHMCforces[nr]->doOmegaStep(para_double);
  }
  
  if (mode == pHMCPropagator_PROPAGATEOMEGAMOMENTA) { 
    pHMCforces[nr]->doOmegaMomentumStep(para_double);
  }
  
  if (mode == pHMCPropagator_CALCOMEGAACTION) { 
    para_double = pHMCforces[nr]->calcOmegaAction(roundToInt(para_double), (double*) phiField);
  }
  
  if (mode == pHMCPropagator_CALCOMEGAMOMENTUMACTION) { 
    para_double = pHMCforces[nr]->calcOmegaMomentumAction();
  }

  if (mode == pHMCPropagator_CALCEXACTMMDAGINVERSESQRTOMEGAACTION) { 
    int neededIter = 0;
    para_double = pHMCforces[nr]->calcExactMMdagInverseSQRTomegaAction((double*) phiField, para_double, neededIter);
  }

  if (mode == pHMCPropagator_RESETEXACTMMDAGINVERSESQRTOMEGAACTION) { 
    pHMCforces[nr]->setExactOmegaMMdagInverseSQRTAction(NaN);
    para_double = pHMCforces[nr]->getExactOmegaMMdagInverseSQRTAction();
  }

  if (mode == pHMCPropagator_SAVEPROP) { 
    savePhiFields();
  }
  
  if (mode == pHMCPropagator_SAVEFORCE) { 
    pHMCforces[nr]->saveOmegaField();
  }
  
  if (mode == pHMCPropagator_WRITEOMEGAFORCESTRENGTHTODISK) { 
    char* fileName = buildOmegaForceStrengthSaveFileName((para_int==0));
    if (para_int==0) {
      pHMCforces[nr]->writeOmegaForceStrengthToDiskPREC(fileName);
    } else {
      pHMCforces[nr]->writeOmegaForceStrengthToDiskGLOBAL(fileName);    
    }
    delete[] fileName;
  }

  if (mode == pHMCPropagator_READOMEGAFORCESTRENGTHFROMDISK) { 
    char* fileName = buildOmegaForceStrengthSaveFileName((para_int==0));
    if (para_int==0) {
      pHMCforces[nr]->readOmegaForceStrengthFromDiskPREC(fileName);
    } else { 
      pHMCforces[nr]->readOmegaForceStrengthFromDiskGLOBAL(fileName);    
    }
    delete[] fileName;
  }

  if (mode == pHMCPropagator_SETOMEGAMASSADAPTIONMODE) { 
    pHMCforces[nr]->setOmegaMassAdaptionMode((bool)para_int);
  }

  if (mode == pHMCPropagator_CALCOMEGAMASSADAPTION) { 
    pHMCforces[nr]->calcOmegaMassAdaption();
  }
 
  if (mode == Propagator_FORCE) {
    double S;
    pHMCforces[nr]->calcAllForces(para_int, (double*) phiField, S);
    para_double = S;
  }
  
  if (mode == pHMCPropagator_SETNU) { 
    pHMCforces[nr]->setNu(para_int, para_double);
  }  

  if (mode == pHMCPropagator_SETTHETA) { 
    pHMCforces[nr]->setTheta(para_double);
  }  

  if (mode == pHMCPropagator_WRITEOMEGAFIELDTODISK) { 
    char* fileName = buildForceOmegaFieldSaveFileName(nr);
    pHMCforces[nr]->writeOmegaFieldToDisk(fileName);
    delete[] fileName;
  } 

  if (mode == pHMCPropagator_READOMEGAFIELDFROMDISK) { 
    char* fileName = buildForceOmegaFieldSaveFileName(nr);
    para_int = pHMCforces[nr]->readOmegaFieldFromDisk(fileName);
    delete[] fileName;
  } 
  
  para_int = 211;
}


void pHMCPropagator::SlaveController() {
  #ifdef UseMPI
  int VL = fermiOps->getVectorLength();
  int nr, mode, para_int;
  double para_double;
  int precUseInt;
  double precM, precS;
  bool precUse;
  double* data = NULL;
  
  if (LogLevel>2) printf("Entering Slave Control Loop (on node %d)...\n",ownNodeID);
  while (true) {
    MPI_Status* status = new MPI_Status;
    MPI_Recv(&nr, 1 , MPI_INT , 0, 1, MPI_COMM_WORLD, status);
    MPI_Recv(&mode, 1 , MPI_INT , 0, 2, MPI_COMM_WORLD, status);
    MPI_Recv(&para_int, 1 , MPI_INT , 0, 3, MPI_COMM_WORLD, status);
    MPI_Recv(&para_double, 1 , MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, status);
    MPI_Recv(&precUseInt, 1 , MPI_INT , 0, 5, MPI_COMM_WORLD, status);
    MPI_Recv(&precM, 1 , MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, status);
    MPI_Recv(&precS, 1 , MPI_DOUBLE, 0, 7, MPI_COMM_WORLD, status);

    if (mode == pHMCPropagator_SETPHI) {
      MPI_Recv(phiField, VL/2, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, status);
    }
  
    if (mode == Propagator_SETFOURIERACCELERATIONMASSES) {
      MPI_Recv(MomentaMasses, VL/8, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, status);
    }
      
    if (mode == pHMCPropagator_SETPOLYROOTS) {
      data = new double[2+2*para_int];
      MPI_Recv(data, 2+2*para_int, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, status);
    }

    if (mode == pHMCPropagator_ACTIVATEFORCESTORING) {
      data = new double[para_int];
      MPI_Recv(data, para_int, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, status);
    }
    
    if (mode == pHMCPropagator_SYNCQPRECDATA) {
      data = new double[2];
      MPI_Recv(data, 2, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, status);
    }

    if (mode == pHMCPropagator_SYNCRPRECDATA) {
      data = new double[2];
      MPI_Recv(data, 2, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, status);
    }  

    delete status;
    precUse = (bool) precUseInt;
    fermiOps->setPreconditioner(precUse, precM, precS);

    ThreadController(nr, mode, para_int, para_double, data);
    
    MPI_Send(&nr, 1 , MPI_INT , 0, 11, MPI_COMM_WORLD);
    MPI_Send(&mode, 1 , MPI_INT , 0, 12, MPI_COMM_WORLD);
    MPI_Send(&para_int, 1 , MPI_INT , 0, 13, MPI_COMM_WORLD);
    MPI_Send(&para_double, 1 ,  MPI_DOUBLE, 0, 14, MPI_COMM_WORLD);

    if (mode == Propagator_FORCE) {
      double* dSdPhi = (double*) pHMCforces[nr]->getdSdPhi();
      MPI_Send(dSdPhi, VL/2, MPI_DOUBLE, 0, 15, MPI_COMM_WORLD);
    }

    if (mode == pHMCPropagator_GETAVERAGEPRECSUBPOL0OMEGAFORCESTRENGTH) {
      double* avgOmForces = pHMCforces[nr]->getAverageOmegaForceStrengthPREC();
      MPI_Send(avgOmForces, VL/8, MPI_DOUBLE, 0, 15, MPI_COMM_WORLD);
    }

    if (data != NULL) { 
      delete[] data;
    }
    data = NULL;

    if (mode == Propagator_EXIT) break;
  }
  #endif 
}


void pHMCPropagator::threadedExecutePROP(int mode, double eps) {
  int* numbers = new int[nodeCount];
  int* modes = new int[nodeCount];
  int* para_int = new int[nodeCount];
  double* para_double = new double[nodeCount];
  double** data = new double*[nodeCount];
  double** tempData = new double*[nodeCount];
  int* dataCount = new int[nodeCount];
  int I;
 

  double precM;
  double precS;
  bool precUse;
  fermiOps->getPreconditionerParameter(precUse, precM, precS);
  for (I=0; I<nodeCount; I++) {
    numbers[I] = -1;
    modes[I] = mode;
    para_int[I] = 0;
    para_double[I] = 0;
    dataCount[I] = 0;
    data[I] = NULL;
    tempData[I] = NULL;
  }

  if (mode == pHMCPropagator_SETQUASIHERMITEANMODE) {
    for (I=0; I<nodeCount; I++) {
      para_int[I] = (int) eps;
    }    
  }
  
  if (mode == pHMCPropagator_SETTUNEMODE) {
    for (I=0; I<nodeCount; I++) {
      para_int[I] = (int) eps;
    }    
  }

  if (mode == pHMCPropagator_SETMODELSELECTION) {
    for (I=0; I<nodeCount; I++) {
      para_int[I] = (int) eps;
    }    
  }

  if (mode == pHMCPropagator_SYNCQPRECDATA) {
    bool precQUse;
    double QPrecMu, QPrecBeta;
    fermiOps->getQPreconditionerParameter(precQUse, QPrecMu, QPrecBeta);
  
    for (I=0; I<nodeCount; I++) {
      para_int[I] = (int) precQUse;
      tempData[I] = new double[2];
      tempData[I][0] = QPrecMu;
      tempData[I][1] = QPrecBeta;
      dataCount[I] = 2;
      data[I] = tempData[I];
    }    
  }
  
  if (mode == pHMCPropagator_SYNCRPRECDATA) {
    bool precRUse;
    double RPrecM, RPrecF;
    fermiOps->getRPreconditionerParameter(precRUse, RPrecM, RPrecF);
  
    for (I=0; I<nodeCount; I++) {
      para_int[I] = (int) precRUse;
      tempData[I] = new double[2];
      tempData[I][0] = RPrecM;
      tempData[I][1] = RPrecF;
      dataCount[I] = 2;
      data[I] = tempData[I];
    }    
  }
  
  if (mode == pHMCPropagator_SETRANDOMSEED) {
    for (I=0; I<nodeCount; I++) {
      para_int[I] = AdvancedSeed;
    }    
  }
  if (mode == pHMCPropagator_DRAWRANDOMNUMBERSONSLAVENODES) {
    for (I=1; I<nodeCount; I++) {
      para_int[I] = (int) eps;
    }    
    para_int[0] = 0;
  }
  if (mode == pHMCPropagator_SETPHI) {
    for (I=1; I<nodeCount; I++) {
      para_int[I] = fermiOps->getVectorLength()/2;
      dataCount[I] = para_int[I];
      data[I] = (double*) phiField;
    }    
  }
  if (mode == Propagator_SETFOURIERACCELERATIONTYPE) {
    for (I=0; I<nodeCount; I++) {
      para_int[I] = phiForceFourierType;
      para_double[I] = phiForceFourierPara;
    }    
  }
  if (mode == Propagator_SETFOURIERACCELERATIONMASSES) {
    for (I=1; I<nodeCount; I++) {
      para_int[I] = fermiOps->getVectorLength()/8;
      dataCount[I] = para_int[I];
      data[I] = (double*) MomentaMasses;
    }    
  }
  if (mode == pHMCPropagator_PROPAGATEPHI) {
    for (I=0; I<nodeCount; I++) {
      para_double[I] = eps;
    }    
  }  
  if (mode == pHMCPropagator_MULTIPLYPHIMOMENTA) {
    for (I=0; I<nodeCount; I++) {
      para_double[I] = eps;
    }    
  }
  
  for (I=1; I<nodeCount; I++) {
    #ifdef UseMPI
      int precUseInt = (int) precUse;
      int remoteNode = I;
      if (LogLevel>4) printf("Remote call for node %d (mode = %d)...\n",remoteNode, mode);
      
      MPI_Send(&(numbers[I]), 1, MPI_INT, remoteNode, 1, MPI_COMM_WORLD);
      MPI_Send(&(modes[I]), 1, MPI_INT, remoteNode, 2, MPI_COMM_WORLD);
      MPI_Send(&(para_int[I]), 1, MPI_INT, remoteNode, 3, MPI_COMM_WORLD);
      MPI_Send(&(para_double[I]), 1, MPI_DOUBLE, remoteNode, 4, MPI_COMM_WORLD);
      MPI_Send(&precUseInt, 1, MPI_INT, remoteNode, 5, MPI_COMM_WORLD);
      MPI_Send(&precM, 1, MPI_DOUBLE, remoteNode, 6, MPI_COMM_WORLD);
      MPI_Send(&precS, 1, MPI_DOUBLE, remoteNode, 7, MPI_COMM_WORLD);
      if (dataCount[I] > 0) { 
        MPI_Send(data[I], dataCount[I], MPI_DOUBLE_PRECISION, remoteNode, 8, MPI_COMM_WORLD);
      }
    #endif 
  }

  //Einen Thread selber ausfuehren
  if (LogLevel>4) printf("Local call for node = %d (mode = %d)...\n",0,mode);

  ThreadController(numbers[0], modes[0], para_int[0], para_double[0], data[0]);

  //Read results from processes
  for (I=1; I<nodeCount; I++) {
    #ifdef UseMPI
      int remoteNode = I;
      MPI_Status* status = new MPI_Status;

      MPI_Recv(&(numbers[I]), 1 , MPI_INT , remoteNode, 11, MPI_COMM_WORLD, status);
      MPI_Recv(&(modes[I]), 1 , MPI_INT , remoteNode, 12, MPI_COMM_WORLD, status);
      MPI_Recv(&(para_int[I]), 1 , MPI_INT , remoteNode, 13, MPI_COMM_WORLD, status);
      MPI_Recv(&(para_double[I]), 1 , MPI_DOUBLE, remoteNode, 14, MPI_COMM_WORLD, status);

      if (para_int[I] != 211) {
        printf("ERROR in MPI - Connection!!!\n");
	exit(0);
      }

      delete status;
      if (LogLevel>4) printf("Remote call for nr = %d on node %d (mode = %d) has been processed!\n",I,remoteNode, mode);
    #endif 
  }

  if (mode == pHMCPropagator_CHECKRANDOMSYNCHRON) {
    for (I=1; I<nodeCount; I++) {
      if (para_double[I] != para_double[0]) {
        printf("ERROR: Random numbers out of sync (%f not equal to %f)!!!\n",para_double[I],para_double[0]);
	exit(0);
      }
    }
    syncRandom = para_double[0];
    if (LogLevel>4) printf("Synchronized Random Number: %f\n",syncRandom);
  }
  
  delete[] numbers;
  delete[] modes;
  delete[] para_int;
  delete[] para_double;
  delete[] data;
  for (I=0; I<nodeCount; I++) {
    if (tempData[I] != NULL) delete[] tempData[I];
  }
  delete[] tempData;
  delete[] dataCount;
}


void pHMCPropagator::threadedExecuteFORCE(int mode, double eps) {
  int* numbers = new int[forceCount];
  int* modes = new int[forceCount];
  int* para_int = new int[forceCount];
  double* para_double = new double[forceCount];
  double** data = new double*[forceCount];
  int* dataCount = new int[forceCount];
  double** tempData = new double*[forceCount];
  int I,I2;

  double precM;
  double precS;
  bool precUse;
  fermiOps->getPreconditionerParameter(precUse, precM, precS);
  
  for (I=0; I<forceCount; I++) {
    numbers[I] = I;
    modes[I] = mode;
    para_int[I] = 0;
    para_double[I] = 0;
    data[I] = NULL;
    dataCount[I] = 0;
    tempData[I] = NULL;
  }

  if (mode == pHMCPropagator_SETPOLYROOTS) {
    int polNr = roundToInt(eps);
    for (I=0; I<forceCount; I++) {
      para_double[I] = pHMCforces[I]->getNu(polNr);    
      para_int[I] = pHMCforces[I]->getApproxPolyDegree(polNr);
      dataCount[I] = 2+2*pHMCforces[I]->getApproxPolyDegree(polNr);
      data[I] = (double*) pHMCforces[I]->getApproxPolyRoots(polNr);
      
      tempData[I] = new double[dataCount[I]];
      tempData[I][0] = polNr;
      tempData[I][1] = pHMCforces[I]->getApproxPolyLambda(polNr);      
      for (I2=2; I2<dataCount[I]; I2++) tempData[I][I2] = data[I][I2-2];
      data[I] = tempData[I];
    }    
  }

  if (mode == pHMCPropagator_GETAVERAGEPRECSUBPOL0OMEGAFORCESTRENGTH) {
  }  

  if (mode == pHMCPropagator_RESETAVERAGEOMEGAFORCESTRENGTH) {
  }  
  
  if (mode == pHMCPropagator_CLEARACTSTORE) {
    for (I=0; I<forceCount; I++) {
      pHMCforces[I]->clearActStore();
    }
  }  
  
  if (mode == pHMCPropagator_CLEARTRASTARTSTORE) {
    for (I=0; I<forceCount; I++) {
      pHMCforces[I]->clearTraStartStore();
    }
  }  
  
  if (mode == pHMCPropagator_COPYTRASTARTSTORETOACTSTORE) {
    for (I=0; I<forceCount; I++) {
      pHMCforces[I]->copyTraStartStoreToActStore();
    }
  }  
  
  if (mode == pHMCPropagator_COPYACTSTORETOTRASTARTSTORE) {
    for (I=0; I<forceCount; I++) {
      pHMCforces[I]->copyActStoreToTraStartStore();
    }
  }  
  
  if (mode == pHMCPropagator_STORECURRENTTOPLEVELFORCES) {
  }  
  
  if (mode == pHMCPropagator_ACTIVATEFORCESTORING) {
    bool acti = (bool) roundToInt(eps);
    for (I=0; I<forceCount; I++) {
      if (acti) {
        para_int[I] = 1+subPolyCount;
        tempData[I] = new double[1+subPolyCount];
        data[I] = tempData[I];
        dataCount[I] = 1+subPolyCount;

        for (int I2=0; I2<subPolyCount; I2++) {
          tempData[I][I2] = 1+2*I2;
        }
        tempData[I][subPolyCount] = 2*subPolyCount;
      } else {
        para_int[I] = 0;
        data[I] = NULL;
        dataCount[I] = 0;
      }
    }    
  }  
  
  if (mode == pHMCPropagator_SETNU) {
    int polNr = roundToInt(eps);
    for (I=0; I<forceCount; I++) {
      para_int[I] = polNr;
      para_double[I] = pHMCforces[I]->getNu(polNr);
    }    
  }  
  
  if (mode == Propagator_FORCE) {
    int polNr = roundToInt(eps);
    for (I=0; I<forceCount; I++) {
      para_int[I] = polNr;
    }    
  }  
  
  if (mode == pHMCPropagator_CALCOMEGAACTION) {
    for (I=0; I<forceCount; I++) {
      para_double[I] = roundToInt(fabs(eps));
      para_int[I] = (int)false;
      if (eps<-1E-9) para_int[I] = (int)true;
    }    
  }  
  
  if (mode == pHMCPropagator_SETTHETA) {
    for (I=0; I<forceCount; I++) {
      para_double[I] = globalForceTheta;
    }    
  }    
  
  if (mode == pHMCPropagator_PROPAGATEOMEGA) {
    for (I=0; I<forceCount; I++) {
      para_double[I] = eps;
    }      
  }
  if (mode == pHMCPropagator_PROPAGATEOMEGAMOMENTA) {
    for (I=0; I<forceCount; I++) {
      para_double[I] = eps;
    }     
  }
  if (mode == pHMCPropagator_SAMPLEOMEGAMOMENTA) {
    for (I=0; I<forceCount; I++) {
      para_int[I] = (int) (eps+0.01);
    }    
  }       
  if (mode == pHMCPropagator_SAMPLEOMEGAFIELDS) {
    for (I=0; I<forceCount; I++) {
      para_int[I] = (int) (eps+0.01);
    }    
  }         
  if (mode == pHMCPropagator_MULTIPLYOMEGAMOMENTA) {
    for (I=0; I<forceCount; I++) {
      para_double[I] = eps;
    }    
  }    
  
  if (mode == pHMCPropagator_WRITEOMEGAFORCESTRENGTHTODISK) {
    for (I=0; I<forceCount; I++) {
      para_int[I] = roundToInt(eps);
    }    
  }    
  if (mode == pHMCPropagator_READOMEGAFORCESTRENGTHFROMDISK) {
    for (I=0; I<forceCount; I++) {
      para_int[I] = roundToInt(eps);
    }    
  }    
  if (mode == pHMCPropagator_SETOMEGAMASSADAPTIONMODE) {
    for (I=0; I<forceCount; I++) {
      para_int[I] = roundToInt(eps);
    }    
  }    
  if (mode == pHMCPropagator_CALCOMEGAMASSADAPTION) {
    for (I=0; I<forceCount; I++) {
      para_int[I] = roundToInt(eps);
    }    
  }    
  if (mode == pHMCPropagator_CALCEXACTMMDAGINVERSESQRTOMEGAACTION) {
    for (I=0; I<forceCount; I++) {
      para_double[I] = polyAlpha;
    }    
  }    

  int execForces = 0;
  while (execForces<forceCount) {
    int rem = forceCount-execForces;
    if (rem>nodeCount) rem = nodeCount;
    
    for (I=execForces+1; I<execForces+rem; I++) {
    #ifdef UseMPI
      int precUseInt = (int) precUse;
      int remoteNode = I % nodeCount;
      if (LogLevel>4) printf("Remote call for nr = %d on node %d (mode = %d)...\n",I,remoteNode, mode);
      
      MPI_Send(&(numbers[I]), 1, MPI_INT, remoteNode, 1, MPI_COMM_WORLD);
      MPI_Send(&(modes[I]), 1, MPI_INT, remoteNode, 2, MPI_COMM_WORLD);
      MPI_Send(&(para_int[I]), 1, MPI_INT, remoteNode, 3, MPI_COMM_WORLD);
      MPI_Send(&(para_double[I]), 1, MPI_DOUBLE, remoteNode, 4, MPI_COMM_WORLD);
      MPI_Send(&precUseInt, 1, MPI_INT, remoteNode, 5, MPI_COMM_WORLD);
      MPI_Send(&precM, 1, MPI_DOUBLE, remoteNode, 6, MPI_COMM_WORLD);
      MPI_Send(&precS, 1, MPI_DOUBLE, remoteNode, 7, MPI_COMM_WORLD);
      if (dataCount[I] > 0) { 
        MPI_Send(data[I], dataCount[I], MPI_DOUBLE_PRECISION, remoteNode, 8, MPI_COMM_WORLD);
      }
    #endif 
    }

    //Einen Thread selber ausfuehren
    if (LogLevel>4) printf("Local call for nr = %d (mode = %d)...\n",execForces,mode);

    ThreadController(numbers[execForces], modes[execForces], para_int[execForces], para_double[execForces], data[execForces]);

    //Read results from processes
    for (I=execForces+1; I<execForces+rem; I++) {
    #ifdef UseMPI
      int remoteNode = I % nodeCount;
      MPI_Status* status = new MPI_Status;

      MPI_Recv(&(numbers[I]), 1 , MPI_INT , remoteNode, 11, MPI_COMM_WORLD, status);
      MPI_Recv(&(modes[I]), 1 , MPI_INT , remoteNode, 12, MPI_COMM_WORLD, status);
      MPI_Recv(&(para_int[I]), 1 , MPI_INT , remoteNode, 13, MPI_COMM_WORLD, status);
      MPI_Recv(&(para_double[I]), 1 , MPI_DOUBLE, remoteNode, 14, MPI_COMM_WORLD, status);

      if (para_int[I] != 211) {
        printf("ERROR in MPI - Connection (return value: %d)!!!\n", para_int[I]);
	exit(0);
      }

      if (modes[I] == Propagator_FORCE) {
        int VL = fermiOps->getVectorLength();
        double* dSdPhi = (double*) pHMCforces[I]->getdSdPhi();
        MPI_Recv(dSdPhi, VL/2, MPI_DOUBLE, remoteNode, 15, MPI_COMM_WORLD, status);      
        int polNr = roundToInt(fabs(eps));
        pHMCforces[I]->setActOmegaAction(polNr, para_double[I]);
      }
      if (mode == pHMCPropagator_GETAVERAGEPRECSUBPOL0OMEGAFORCESTRENGTH) {
        int VL = fermiOps->getVectorLength();
        double* avgForce = pHMCforces[I]->getAverageOmegaForceStrengthPREC();
        MPI_Recv(avgForce, VL/8, MPI_DOUBLE, remoteNode, 15, MPI_COMM_WORLD, status);      
      }  

      if (modes[I] == pHMCPropagator_CALCOMEGAACTION) {
        int polNr = roundToInt(fabs(eps));
        pHMCforces[I]->setActOmegaAction(polNr, para_double[I]);
      }
      if (modes[I] == pHMCPropagator_CALCOMEGAMOMENTUMACTION) {
        pHMCforces[I]->setOmegaMomentumAction(para_double[I]);
      }
      if (modes[I] == pHMCPropagator_CALCEXACTMMDAGINVERSESQRTOMEGAACTION) {
        pHMCforces[I]->setExactOmegaMMdagInverseSQRTAction(para_double[I]);
      }
      if (modes[I] == pHMCPropagator_RESETEXACTMMDAGINVERSESQRTOMEGAACTION) {
        pHMCforces[I]->setExactOmegaMMdagInverseSQRTAction(para_double[I]);
      }

      delete status;
      if (LogLevel>4) printf("Remote call for nr = %d on node %d (mode = %d) has been processed!\n",I,remoteNode, mode);
    #endif 
    }

    execForces += rem;
  }
  
  delete[] numbers;
  delete[] modes;
  delete[] para_int;
  delete[] para_double;
  delete[] data;
  delete[] dataCount;
  for (I=0; I<forceCount; I++) {
    if (tempData[I] != NULL) delete[] tempData[I];
  }
  delete[] tempData;
}


void pHMCPropagator::threadedExecute(int mode, double eps) {
  if (ownNodeID != 0) return;
  
  if ((mode<0) || (mode>=50)) {
    threadedExecutePROP(mode, eps);
  }

  if (forceCount <= 0) return;
  
  if ((mode>=0) && (mode<50)) {
    threadedExecuteFORCE(mode, eps);
  }
}


void pHMCPropagator::negateAllMomenta() {
  threadedExecute(pHMCPropagator_MULTIPLYPHIMOMENTA,-1.0);
  threadedExecute(pHMCPropagator_MULTIPLYOMEGAMOMENTA,-1.0);
}


void pHMCPropagator::setTheta(double tht) {
  globalForceTheta = tht;
  threadedExecute(pHMCPropagator_SETTHETA,0);  
}


void pHMCPropagator::setPolyData() {
  int I,I2,I3;
  
  for (I=0; I<1+precMassCount; I++) {  
    double M1 = 0.0;
    double M2 = 0.0;
    if (I>0) M1 = precMasses[I-1];
    if (I<precMassCount) M2 = precMasses[I];
    approxPoly[I] = new PolynomialApproximation(subPolyCount, I, polyDegree, polyDigit, polyAlpha, M1, M2, polyEpsilon, polyLambda);
    approxPoly[I]->plotApproxPolynomials();
    for (I2=0; I2<1+2*subPolyCount; I2++) {
      polyRoots[I][I2] = approxPoly[I]->getApproxPolyRoots(I2);
    }
    if (LogLevel>3) {
      printf("Setting Polynomial roots on Propagator...\n");
      for (I2=0; I2<1+2*subPolyCount; I2++) {
        printf("Roots of Master Polynomial PrecMass%d for subpolynomial %d: \n",I,I2);
        for (I3=0; I3<polyDegree[I2/2]; I3++) {
          printf(" -> roots nr %d: ",I3);
          polyRoots[I][I2][I3].print();
	}
      }
    }
  }
}


void pHMCPropagator::distributeMonomialsToForces() {
  if ((fermiOps->getYN()<=0) || (polyDegree[0]<=0) || (maxPolyDegreePerNode <= 0)) {
    return;
  }
  
  double polLam = polyLambda[0];
  int polyPartCount = polyDegree[0] / maxPolyDegreePerNode;
  if ((polyDegree[0] % maxPolyDegreePerNode)>0) {
    polyPartCount++;
  }
  if ((subPolyCount>0) && (polyPartCount>1)) {
    printf("ERROR: Sub Polynomials and Monomial-Distribution not compatible!!!\n");
    exit(0);
  }
  if (polyPartCount>1) {
    Complex** polyParts = new Complex*[polyPartCount];
    int* polyPartLengths = new int[polyPartCount];
    double* polyPartNorms = new double[polyPartCount];
    int I,I2;
    for (I2=0; I2<1+precMassCount; I2++) {
      int pos = 0;
      int count = 0;
      while (pos<polyDegree[0]) {
        polyPartLengths[count] = maxPolyDegreePerNode;
        if (polyDegree[0]-pos<maxPolyDegreePerNode) {
          polyPartLengths[count] = polyDegree[0]-pos;    
        }
        polyParts[count] = new Complex[polyPartLengths[count]];
        polyPartNorms[count] = approxPoly[I2]->getPartPolyNormalization(0,polyPartLengths[count]);
        for (I=0; I<polyPartLengths[count]; I++) {
          polyParts[count][I] = polyRoots[I2][0][pos];
          pos++;
        }
        count++;
      }

      for (I=0; I<Nf*polyPartCount; I++) {
        pHMCforces[I2*Nf*polyPartCount + I]->setApproxPolyRoots(0,polyParts[I%polyPartCount],polyPartLengths[I%polyPartCount],polyPartNorms[I%polyPartCount], polLam);
      }
    }
    delete[] polyPartLengths;
    delete[] polyPartNorms;
    for (I=0; I<polyPartCount; I++) {
      delete[] polyParts[I];
    }
    delete[] polyParts;  
  } else {
    int count = 0;
    int I,I2,I3;    
    for (I=0; I<Nf; I++) {
      for (I2=0; I2<1+precMassCount; I2++) {
        for (I3=0; I3<1+2*subPolyCount; I3++) {
	  double n = approxPoly[I2]->getPartPolyNormalization(I3,polyDegree[I3/2]);
          pHMCforces[count]->setApproxPolyRoots(I3,polyRoots[I2][I3],polyDegree[I3/2],n, polLam);
	}
     
        count++;
      }
    }
  }
}

  
void pHMCPropagator::setNf(int nf) {
  Nf = nf;
  desiniForceCalculators();
  delete[] pHMCforces;
  if ((fermiOps->getYN()<=0) || (polyDegree[0]<=0) || (maxPolyDegreePerNode <= 0)) {
    forceCount = 0;
    return;
  }
  if (Nf<0) {
    Nf = 0;
  }

  forceCount = polyDegree[0] / maxPolyDegreePerNode;
  if ((polyDegree[0] % maxPolyDegreePerNode)>0) {
    forceCount++;
  }
  if ((subPolyCount>0) && (forceCount>1)) {
    printf("ERROR: Sub Polynomials and Monomial-Distribution not compatible!!!\n");
    exit(0);
  }
  forceCount *= Nf * (1+precMassCount);
  if (forceCount<nodeCount) {
    printf("More Nodes than Forces ==> WASTE !!!\n");
    exit(0);
  }

  pHMCforces = new pHMCForce*[forceCount];
  forces = new Force*[forceCount];    
  int I;
  for (I=0; I<forceCount; I++) {
    bool local = false;
    if ((I%nodeCount) == ownNodeID) local = true;    
    pHMCforces[I] = new pHMCForce(fermiOps, local, globalForceTheta,I,AdditionalAuxVectorsCount);    
    forces[I] = pHMCforces[I];
  }  
  distributeMonomialsToForces();
	printf ("Exciting setNf");
}


double pHMCPropagator::calcOmegaAction() {
  double S = 0;

  if (fermiOps->getYN()>0) {
    S = 0;
    for (int I2=0; I2<subPolyCount; I2++) {
      int polynomSlot = 1+2*I2;
      bool Savail = true;
      for (int I=0; I<forceCount; I++) {
        double dummy = pHMCforces[I]->getActOmegaAction(polynomSlot);
        if (isNaN(dummy)) Savail = false;
      }
      if (!Savail) {
        double para = polynomSlot+0.1;
        threadedExecute(pHMCPropagator_CALCOMEGAACTION, para);
      } else {
        if (LogLevel>4) printf("Taking Action from Data-Pool for polynom %d\n",polynomSlot);        
      }
      for (int I=0; I<forceCount; I++) {
        S = S + pHMCforces[I]->getActOmegaAction(polynomSlot);
      }
    }
    
    bool Savail = true;
    for (int I=0; I<forceCount; I++) {
      double dummy = pHMCforces[I]->getActOmegaAction(2*subPolyCount);
      if (isNaN(dummy)) Savail = false;
    }
    if (!Savail) {
      double para = 2*subPolyCount+0.1;
      threadedExecute(pHMCPropagator_CALCOMEGAACTION, para);
    } else {
      if (LogLevel>4) printf("Taking Action from Data-Pool for polynom %d\n",2*subPolyCount);        
    }
    for (int I=0; I<forceCount; I++) {
      S = S + pHMCforces[I]->getActOmegaAction(2*subPolyCount);
    }
  } else {
    S = 0;
  }
  return S;
}


double pHMCPropagator::calcOmegaMomentumAction() {
  double S = 0;

  if (fermiOps->getYN()>0) {
    threadedExecute(pHMCPropagator_CALCOMEGAMOMENTUMACTION, 0);
    S = 0;
    int I;
    for (I=0; I<forceCount; I++) {
      S = S + pHMCforces[I]->getOmegaMomentumAction();
    }
  } else {
    S = 0;
  }
  return S;
}


double pHMCPropagator::calcTotalAction(double finalTOL) {
  double S1 = calcOmegaAction();
  double S2 = calcOmegaMomentumAction();  
  double S3 = calcPhiMomentumAction();
  double S4 = calcPhiAction();
  if (LogLevel>3) printf("Action composed of OM: %f + OMp: %f + PHIp: %f + PHI: %f\n",S1,S2,S3,S4); 
  double S = S1+S2+S3+S4;
  if (!(S==S)) {
    printf("ERROR: Action is NaN!!!\n");
    if (!tuneMode) {
      exit(0);
    } else {
      S = 1E100;
    }
  }
  return S;
}


void pHMCPropagator::improveRPreconditioningParameters(int thermStepID, int ParameterAdaptionMode, double PolLambda, double upperEWsafetyFac) {
  bool RprecUse;
  double RprecM, RprecF;
  fermiOps->getRPreconditionerParameter(RprecUse,RprecM,RprecF);

  if (!RprecUse) return;
  if (forceCount<=0) return;
  if (fermiOps->getYN()<=0) return;

  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);

  if (ParameterAdaptionMode==2) {
    if (LogLevel>1) printf("R - Preconditioning-Improver: Full parameter adaption...\n");
    double avgPhi, avgStagPhi, avgNorm, sigmaNorm;
    measure(avgPhi, avgStagPhi, avgNorm, sigmaNorm);
    
    if (upperEWsafetyFac >= 0) {
      RprecM = avgPhi;
      if (RprecM < 0.35) RprecM = 0.35;
    } else {
      RprecM = 10.0 * ((-upperEWsafetyFac) - ((int)(-upperEWsafetyFac)));
    }
    synchronizedChangeOfRPreconsitionerData(RprecUse,RprecM,RprecF);
  }

  if (ParameterAdaptionMode==0) {
    if (LogLevel>1) printf("R - Preconditioning-Improver: Only measuring largest EW...\n");
  }
  if (ParameterAdaptionMode==1) {
    if (LogLevel>1) printf("R - Preconditioning-Improver: Only global factor...\n");
  }

  double TOL = 1E-2;
  Complex* EWhigh = fermiOps->calcFermionMatrixARPACKEigenValues(4, 1, (double*)phiField, TOL, true, NULL, quasiHermiteanMode, true);
  upperEWboundLog[upperEWboundLogCount][0] = thermStepID;
  upperEWboundLog[upperEWboundLogCount][1] = RprecM;
  upperEWboundLog[upperEWboundLogCount][2] = RprecF;
  upperEWboundLog[upperEWboundLogCount][3] = EWhigh[0].x;
  upperEWboundLogCount++;
  delete[] EWhigh;
  
  if (ParameterAdaptionMode>0) {
    double largestEW = upperEWboundLog[upperEWboundLogCount-1][3]/sqr(upperEWboundLog[upperEWboundLogCount-1][2]);
    int eligibleCount = 0;

    for (int I=5; I<upperEWboundLogCount; I++) {
      if ((fabs(upperEWboundLog[I][0]-thermStepID)<150) && (fabs((upperEWboundLog[I][1]-RprecM)/RprecM)<0.15) && (upperEWboundLog[I][0] > 0)) {
        eligibleCount++;
        double upEW = upperEWboundLog[I][3]/sqr(upperEWboundLog[I][2]);
        if (upEW > largestEW) largestEW = upEW;
      }
    }
    if (eligibleCount<50) {
      eligibleCount = 0;
      for (int I=5; I<upperEWboundLogCount; I++) {
        if ((fabs(upperEWboundLog[I][0]-thermStepID)<150) && (upperEWboundLog[I][0] > 0)) {
          eligibleCount++;
          double upEW = upperEWboundLog[I][3]/sqr(upperEWboundLog[I][2]);
          if (upEW > largestEW) largestEW = upEW;
        }
      }
    }
    if (eligibleCount<30) {
      eligibleCount = 0;
      for (int I=5; I<upperEWboundLogCount; I++) {
        if (fabs(upperEWboundLog[I][0]-thermStepID)<150) {
          eligibleCount++;
          double upEW = upperEWboundLog[I][3]/sqr(upperEWboundLog[I][2]);
          if (upEW > largestEW) largestEW = upEW;
        }
      }
    }

    if (upperEWsafetyFac >= 0) {
      RprecF = 1.0 / sqrt(upperEWsafetyFac*largestEW/PolLambda);
      if (LogLevel>2) printf("Upper EW-bound determined based on %d data sets.\n", eligibleCount);
    } else {
      double f = ((int)-upperEWsafetyFac);
      RprecF = 1.0 / sqrt(f/PolLambda);
      if (LogLevel>2) printf("RprecF is set to explicitly given value.\n");    
    }
    synchronizedChangeOfRPreconsitionerData(RprecUse,RprecM,RprecF);
  }
  addPerformanceProfilingItem("pHMCPropagator::improveRPreconditioningParameters", performanceProfilerStartCycle, 0);
}


void pHMCPropagator::improvePreconditioningParameters(int iterGrob, int iterFein, double testTOL) {
  bool PrecUse;
  double PrecMold, PrecSold;
  fermiOps->getPreconditionerParameter(PrecUse, PrecMold, PrecSold);  
  int randomCount = 0;

  if (!PrecUse) return;
  if (forceCount<=0) return;
  if (fermiOps->getYN()<=0) return;
  
  double bestCond;
  double bestM = PrecMold;
  double bestS = PrecSold;
  double bestEigMin, bestEigMax;
  double eigMin, eigMax;
  fermiOps->exactFermionMatrixConditionNumber((double*) phiField, bestEigMin, bestEigMax, bestCond, true, 1, quasiHermiteanMode);


  double avgPhi, avgStagPhi, avgNorm, sigmaNorm;
  measure(avgPhi, avgStagPhi, avgNorm, sigmaNorm);
  
  if (LogLevel>1) printf("Preconditioning-Improver at avgPhi = %1.3f, avgStagPhi = %1.3f, avgNorm = %1.3f: oldM = %1.3f, oldS = %1.3f with Cond. = %1.5f.\n", avgPhi, avgStagPhi, avgNorm, PrecMold, PrecSold, bestCond);
  double* m = new double[iterGrob+iterFein];
  double* s = new double[iterGrob+iterFein];
  int I;
  
  for (I=0; I<iterGrob; I++) {
/*    m[I] = 0.001 + 2*AdvancedZufall(AdvancedSeed)*avgPhi;
    s[I] = 2*AdvancedZufall(AdvancedSeed)*avgStagPhi;
    randomCount += 2;*/
    if (I>0) {
      m[I] = 0.001 + 1.5*AdvancedZufall(AdvancedSeed)*avgNorm;
      s[I] = 0;
      randomCount += 1;
    } else {
      m[I] = avgNorm;
      s[I] = 0;    
    }
  }
  
  for (I=0; I<iterGrob; I++) {
    fermiOps->setPreconditioner(true, m[I], s[I]);
    double newCond = 0;
    
    fermiOps->exactFermionMatrixConditionNumber((double*) phiField, eigMin, eigMax, newCond, true, 1, quasiHermiteanMode);    
    if (newCond>bestCond) {
      bestCond = newCond;      
      bestEigMin = eigMin;
      bestEigMax = eigMax;
      bestM = m[I];
      bestS = s[I];
      if (LogLevel>2) printf("Preconditioning improved: With PrecM = %1.3f, PrecS = %1.3f better Cond. = %1.5f achieved.\n", bestM, bestS, bestCond);
    }
  }
  
  for (I=0; I<iterFein; I++) {
/*    m[I] = 0.001 + fabs(bestM + 2*(AdvancedZufall(AdvancedSeed)-0.5)*0.1*avgPhi);
    s[I] = fabs(bestS + 2*(AdvancedZufall(AdvancedSeed)-0.5)*0.1*avgStagPhi);
    randomCount += 2;*/    
    m[I] = 0.001 + fabs(bestM + 2*(AdvancedZufall(AdvancedSeed)-0.5)*0.1*avgNorm);
    s[I] = 0;
    randomCount += 1;
  }

  for (I=0; I<iterFein; I++) {
    fermiOps->setPreconditioner(true, m[I], s[I]);
    double newCond = 0;

    fermiOps->exactFermionMatrixConditionNumber((double*) phiField, eigMin, eigMax, newCond, true, 1, quasiHermiteanMode);    
    if (newCond>bestCond) {
      bestCond = newCond;      
      bestEigMin = eigMin;
      bestEigMax = eigMax;
      bestM = m[I];
      bestS = s[I];
      if (LogLevel>2) printf("Preconditioning improved: With PrecM = %1.3f, PrecS = %1.3f better Cond. = %1.5f achieved.\n", bestM, bestS, bestCond);
    }
  }

  fermiOps->setPreconditioner(true, bestM, bestS);
  double PrecMact, PrecSact;
  fermiOps->getPreconditionerParameter(PrecUse, PrecMact, PrecSact);
  if (LogLevel>1) printf("Best Preconditioning data: PrecM = %1.3f, PrecS = %1.3f with Cond. = %1.5f (%1.5f / %1.5f)\n", PrecMact, PrecSact, bestCond, bestEigMin, bestEigMax); 

  threadedExecute(pHMCPropagator_DRAWRANDOMNUMBERSONSLAVENODES, randomCount);

  delete[] m;
  delete[] s;

  PreconditionerWasChanged();
}


void pHMCPropagator::saveALLfields() {
  threadedExecute(pHMCPropagator_SAVEFORCE,0);
  threadedExecute(pHMCPropagator_SAVEPROP,0);
}


void pHMCPropagator::restoreALLfields(bool clearStore) {
  threadedExecute(pHMCPropagator_RESTOREFORCE,0);
  threadedExecute(pHMCPropagator_RESTOREPROP,0);
  
  if (clearStore) clearForceStorage();
}


void pHMCPropagator::sampleALLMomenta() {
  threadedExecute(pHMCPropagator_SAMPLEPHIMOMENTA,0);
  
  if (globalForceTheta != 0) {
    int I;
    for (I=0; I<forceCount; I++) {  
      threadedExecute(pHMCPropagator_SAMPLEOMEGAMOMENTA,I);
    }
  }
}


void pHMCPropagator::sampleOmegaFields() {
  int I;
  for (I=0; I<forceCount; I++) {  
    threadedExecute(pHMCPropagator_SAMPLEOMEGAFIELDS,I);
  }
  clearForceStorage();
}


double pHMCPropagator::getSyncSingleRandom() {
  threadedExecute(pHMCPropagator_CHECKRANDOMSYNCHRON,0);
  return syncRandom;
}


void pHMCPropagator::getNodesReady() {
  int I;
  threadedExecute(pHMCPropagator_SETRANDOMSEED,0);
	printf("getNodesReady: randomseed set\n");
  printf("NPol = %d\n", 1+2*subPolyCount);
	for (I=0; I<1+2*subPolyCount; I++) {
		printf("i = %d\n",I); 
    threadedExecute(pHMCPropagator_SETPOLYROOTS,I);
  }
  threadedExecute(pHMCPropagator_SETPHI,0);
	printf("getNodesReady: phi set\n");
  threadedExecute(pHMCPropagator_SYNCQPRECDATA, 0);
	printf("getNodesReady: sync done\n");
  threadedExecute(pHMCPropagator_SYNCRPRECDATA, 0);
	printf("getNodesReady: done\n");
  nodesReady = true;
}


bool pHMCPropagator::LeapOmelyanMarkovStep(int iterations, double epsilon, double propTOL, double finalTOL, double lambda, double rho, double theta, double mu) {
  int I2;  
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


  threadedExecute(pHMCPropagator_SETPHI,0);
  Sold = calcTotalAction(0);
  calcFullPhiDerivatives(dSdPhi, 0, 0); //calculates also omegaForces

  saveALLfields();
  MiniPhiMomentumStep(dSdPhi, epsilonLambda);
  threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonLambda);
  if (LogLevel>2) printf("Old Action is %f\n", Sold);
  SbeforeProp = Sold;

  for (I2=0; I2<iterations; I2++) {
    if ((lambda==0.5) && (theta==0) && (rho==0) && (mu==0)) {
      //Normaler Leap-Frog Schritt

      MiniPhiStep(epsilon);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilon);            
    } else if ((theta==0) && (rho==0) && (mu==0)) {
      //Omelyan - Schritt mit Order 2

      MiniPhiStep(epsilonHalf);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonHalf);            

      threadedExecute(pHMCPropagator_SETPHI,0);
      calcFullPhiDerivatives(dSdPhi, 0, 0);
      MiniPhiMomentumStep(dSdPhi,epsilonOneMinusTwoLambda);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonOneMinusTwoLambda);

      MiniPhiStep(epsilonHalf);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonHalf);                  

    } else {
      //Omelyan - Schritt mit Order 4

      MiniPhiStep(epsilonRho);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonRho);                  

      threadedExecute(pHMCPropagator_SETPHI,0);
      calcFullPhiDerivatives(dSdPhi, 0, 0);
      MiniPhiMomentumStep(dSdPhi, epsilonTheta);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonTheta);

      MiniPhiStep(epsilonMu);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonMu);                  

      threadedExecute(pHMCPropagator_SETPHI,0);
      calcFullPhiDerivatives(dSdPhi, 0, 0);
      MiniPhiMomentumStep(dSdPhi, epsilonHalfOneMinusTwoLambdaPlusTheta);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonHalfOneMinusTwoLambdaPlusTheta);

      MiniPhiStep(epsilonOneMinusTwoMuPlusRho);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonOneMinusTwoMuPlusRho);                  

      threadedExecute(pHMCPropagator_SETPHI,0);
      calcFullPhiDerivatives(dSdPhi, 0, 0);
      MiniPhiMomentumStep(dSdPhi, epsilonHalfOneMinusTwoLambdaPlusTheta);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonHalfOneMinusTwoLambdaPlusTheta);

      MiniPhiStep(epsilonMu);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonMu);                  
      
      threadedExecute(pHMCPropagator_SETPHI,0);
      calcFullPhiDerivatives(dSdPhi, 0, 0);
      MiniPhiMomentumStep(dSdPhi, epsilonTheta);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonTheta);

      MiniPhiStep(epsilonRho);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonRho);                  
    }

    threadedExecute(pHMCPropagator_SETPHI,0);
    calcFullPhiDerivatives(dSdPhi, 0, 0);
    MiniPhiMomentumStep(dSdPhi, twoEpsilonLambda);
    threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, twoEpsilonLambda);
  }
  //Bring momenta on integer step numbers again
  MiniPhiMomentumStep(dSdPhi, -epsilonLambda);
  
  threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, -epsilonLambda);
  
  Sact = calcTotalAction(0);
  if (LogLevel>2) printf("New Action is %f\n", Sact);
  SafterProp = Sact;
  
  double deltaS = Sact-Sold;
  bool accepted = true;
  if (deltaS>0) {
    double ran = getSyncSingleRandom();  
    if (ran > exp(-deltaS)) {
      restoreALLfields(false);
      accepted = false;
      if (LogLevel>2) printf("Rejecting configuration!!!\n");
    }
  }
  
  if (accepted) {
    Sold = Sact;
    if (phiForceFourierType) {
      calcPhiChangeFourierComponents();
    }
    threadedExecute(pHMCPropagator_COPYACTSTORETOTRASTARTSTORE, 0);    
  } else {
    Sact = Sold;
    threadedExecute(pHMCPropagator_COPYTRASTARTSTORETOACTSTORE, 0);    
  }

  delete[] dSdPhi;
  return accepted;
}


void pHMCPropagator::LeapOmelyanMultiTimeScalePropagation(double* dSdPhi, int level, int* iterations, double epsilon, int* subPolNr, double* propTOL, double finalTOL, double* lambda, double* rho, double* theta, double* mu) {
  if (level == 0) {
    LeapOmelyanPhiPropagation(iterations[0], epsilon, lambda[0], rho[0], theta[0], mu[0]);
    return;
  }

  int I2;  
  double epsilonHalf = 0.5 * epsilon;
  double epsilonLambda = lambda[level] * epsilon;
  double epsilonRho = rho[level] * epsilon;
  double epsilonMu = mu[level] * epsilon;
  double epsilonTheta = theta[level] * epsilon;
  double epsilonHalfOneMinusTwoLambdaPlusTheta = 0.5*(1.0 - 2.0*(lambda[level] + theta[level]))*epsilon;
  double epsilonOneMinusTwoMuPlusRho = (1.0 - 2.0*(mu[level] + rho[level]))*epsilon;  
  double twoEpsilonLambda = 2.0 * lambda[level] * epsilon;
  double epsilonOneMinusTwoLambda = (1.0 - 2.0*lambda[level])*epsilon;
  int nestedIterations = iterations[level-1];

  threadedExecute(pHMCPropagator_SETPHI,0);
  calcFullPhiDerivatives(dSdPhi, subPolNr[level], 2); //calculates also omegaForces
  MiniPhiMomentumStep(dSdPhi, epsilonLambda);
  threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonLambda);

  for (I2=0; I2<iterations[level]; I2++) {
    if ((lambda[level]==0.5) && (theta[level]==0) && (rho[level]==0) && (mu[level]==0)) {
      //Normaler Leap-Frog Schritt

      LeapOmelyanMultiTimeScalePropagation(dSdPhi, level-1, iterations, epsilon/nestedIterations, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);
      if (level==1) threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilon);            
    } else if ((theta[level]==0) && (rho[level]==0) && (mu[level]==0)) {
      //Omelyan - Schritt mit Order 2

      LeapOmelyanMultiTimeScalePropagation(dSdPhi, level-1, iterations, epsilonHalf/nestedIterations, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);
      if (level==1) threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonHalf);            

      threadedExecute(pHMCPropagator_SETPHI,0);
      calcFullPhiDerivatives(dSdPhi, subPolNr[level], 2);
      MiniPhiMomentumStep(dSdPhi,epsilonOneMinusTwoLambda);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonOneMinusTwoLambda);

      LeapOmelyanMultiTimeScalePropagation(dSdPhi, level-1, iterations, epsilonHalf/nestedIterations, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);
      if (level==1) threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonHalf);                  

    } else {
      //Omelyan - Schritt mit Order 4

      LeapOmelyanMultiTimeScalePropagation(dSdPhi, level-1, iterations, epsilonRho/nestedIterations, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);
      if (level==1) threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonRho);                  

      threadedExecute(pHMCPropagator_SETPHI,0);
      calcFullPhiDerivatives(dSdPhi, subPolNr[level], 2);
      MiniPhiMomentumStep(dSdPhi, epsilonTheta);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonTheta);

      LeapOmelyanMultiTimeScalePropagation(dSdPhi, level-1, iterations, epsilonMu/nestedIterations, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);
      if (level==1) threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonMu);                  

      threadedExecute(pHMCPropagator_SETPHI,0);
      calcFullPhiDerivatives(dSdPhi, subPolNr[level], 2);
      MiniPhiMomentumStep(dSdPhi, epsilonHalfOneMinusTwoLambdaPlusTheta);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonHalfOneMinusTwoLambdaPlusTheta);

      LeapOmelyanMultiTimeScalePropagation(dSdPhi, level-1, iterations, epsilonOneMinusTwoMuPlusRho/nestedIterations, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);
      if (level==1) threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonOneMinusTwoMuPlusRho);                  

      threadedExecute(pHMCPropagator_SETPHI,0);
      calcFullPhiDerivatives(dSdPhi, subPolNr[level], 2);
      MiniPhiMomentumStep(dSdPhi, epsilonHalfOneMinusTwoLambdaPlusTheta);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonHalfOneMinusTwoLambdaPlusTheta);

      LeapOmelyanMultiTimeScalePropagation(dSdPhi, level-1, iterations, epsilonMu/nestedIterations, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);
      if (level==1) threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonMu);                  
      
      threadedExecute(pHMCPropagator_SETPHI,0);
      calcFullPhiDerivatives(dSdPhi, subPolNr[level], 2);
      MiniPhiMomentumStep(dSdPhi, epsilonTheta);
      threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, epsilonTheta);

      LeapOmelyanMultiTimeScalePropagation(dSdPhi, level-1, iterations, epsilonRho/nestedIterations, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);
      if (level==1) threadedExecute(pHMCPropagator_PROPAGATEOMEGA, epsilonRho);                  
    }

    threadedExecute(pHMCPropagator_SETPHI,0);
    calcFullPhiDerivatives(dSdPhi, subPolNr[level], 2);
    MiniPhiMomentumStep(dSdPhi, twoEpsilonLambda);
    threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, twoEpsilonLambda);
  }
  //Bring momenta on integer step numbers again
  MiniPhiMomentumStep(dSdPhi, -epsilonLambda);  
  threadedExecute(pHMCPropagator_PROPAGATEOMEGAMOMENTA, -epsilonLambda);
}


void pHMCPropagator::LeapOmelyanMultiTimeScalePropagation(int level, int* iterations, double epsilon, int* subPolNr, double* propTOL, double finalTOL, double* lambda, double* rho, double* theta, double* mu) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double* dSdPhi = new double[4*L0*L1*L2*L3];

  threadedExecute(pHMCPropagator_SETPHI,0);
  saveALLfields();
  Sold = calcTotalAction(0);
  SbeforeProp = Sold;
  if (LogLevel>2) printf("Old Action is %f\n", Sold);

  LeapOmelyanMultiTimeScalePropagation(dSdPhi, level, iterations, epsilon, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);
  
  Sact = calcTotalAction(0);
  SafterProp = Sact;
  if (LogLevel>2) printf("New Action is %f\n", Sact);

  delete[] dSdPhi;
}


bool pHMCPropagator::LeapOmelyanMultiTimeScaleMarkovStep(int level, int* iterations, double epsilon, int* subPolNr, double* propTOL, double finalTOL, double* lambda, double* rho, double* theta, double* mu) {

  if (level==0) return LeapOmelyanMarkovStep(iterations[0], epsilon, propTOL[0], finalTOL, lambda[0], rho[0], theta[0], mu[0]);

  LeapOmelyanMultiTimeScalePropagation(level, iterations, epsilon, subPolNr, propTOL, finalTOL, lambda, rho, theta, mu);
  
  double deltaS = Sact-Sold;
  bool accepted = true;
  if (deltaS>0) {
    double ran = getSyncSingleRandom();  
    if (ran > exp(-deltaS)) {
      restoreALLfields(false);
      accepted = false;
      if (LogLevel>2) printf("Rejecting configuration!!!\n");
    }
  }
  
  if (accepted) {
    Sold = Sact;
    if (phiForceFourierType) {
      calcPhiChangeFourierComponents();
    }
    threadedExecute(pHMCPropagator_COPYACTSTORETOTRASTARTSTORE, 0);
  } else {
    Sact = Sold;
    threadedExecute(pHMCPropagator_COPYTRASTARTSTORETOACTSTORE, 0);
  }

  return accepted;
}


char* pHMCPropagator::buildForceOmegaFieldSaveFileName(int nr) {
  char* fileName = new char[1000];
  
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double y = fermiOps->getYN();
  double rho,r;
  fermiOps->getDiracParameters(rho, r);

  if (fileNameIdentifier == NULL) {
   snprintf(fileName,1000,"%s/data/results/pHMC/omegaFields/omegaFieldL%dx%dx%dx%dNf%dKap%1.5fLam%1.5fY%1.5fRho%1.3fR%1.3fPolDeg%dPolAl%1.3f_%d.dat",
      DataBaseDirectory,L0, L1, L2, L3, Nf, kappa, lambda, y, rho, r, polyDegree[0], polyAlpha, nr);
  }
  else {
   snprintf(fileName,1000,"%s/data/results/pHMC/omegaFields/omegaField%s_%d.dat",
      DataBaseDirectory,fileNameIdentifier, nr);
  }
   
   return fileName;
}  


char* pHMCPropagator::buildOmegaForceStrengthSaveFileName(bool PREC) {
  char* fileName = new char[1000];
  
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double y = fermiOps->getYN();
  double rho,r;
  fermiOps->getDiracParameters(rho, r);
 
  if (PREC) {
    if (fileNameIdentifier == NULL) {
      snprintf(fileName,1000,"%s/data/results/pHMC/FACC/omegaForceStrengthPRECL%dx%dx%dx%dNf%dKap%1.5fLam%1.5fY%1.5fRho%1.3fR%1.3fPolDeg%dPolAl%1.3f",
	 DataBaseDirectory,L0, L1, L2, L3, Nf, kappa, lambda, y, rho, r, polyDegree[0], polyAlpha);
    }
    else {
      snprintf(fileName,1000,"%s/data/results/pHMC/FACC/omegaForceStrengthPREC%s",
	 DataBaseDirectory, fileNameIdentifier);
    }
  } else {
    if (fileNameIdentifier == NULL) {
      snprintf(fileName,1000,"%s/data/results/pHMC/FACC/omegaForceStrengthGLOBALL%dx%dx%dx%dNf%dKap%1.5fLam%1.5fY%1.5fRho%1.3fR%1.3fPolDeg%dPolAl%1.3f",
	 DataBaseDirectory,L0, L1, L2, L3, Nf, kappa, lambda, y, rho, r, polyDegree[0], polyAlpha);
    }
    else {
      snprintf(fileName,1000,"%s/data/results/pHMC/FACC/omegaForceStrengthGLOBAL%s",
	 DataBaseDirectory, fileNameIdentifier);
    }
  }
   
  return fileName;
}  


void pHMCPropagator::writeAllOmegaFieldsToDisk() {
  if (LogLevel>2) printf("Writing all Omega Fields to disk now...\n");
  threadedExecuteFORCE(pHMCPropagator_WRITEOMEGAFIELDTODISK, 0);
}


void pHMCPropagator::readAllOmegaFieldsFromDisk() {
  if (LogLevel>2) printf("Reading all Omega Fields from disk now...\n");
  threadedExecuteFORCE(pHMCPropagator_READOMEGAFIELDFROMDISK, 0);
}


void pHMCPropagator::writeOmegaForceStrengthsToDiskPREC() {
  if (LogLevel>4) printf("Writing PREC Omega Forces to Disk now...\n");
  threadedExecuteFORCE(pHMCPropagator_WRITEOMEGAFORCESTRENGTHTODISK, 0);
}


void pHMCPropagator::writeOmegaForceStrengthsToDiskGLOBAL() {
  if (LogLevel>4) printf("Writing GLOBAL Omega Forces to Disk now...\n");
  threadedExecuteFORCE(pHMCPropagator_WRITEOMEGAFORCESTRENGTHTODISK, 1);
}


void pHMCPropagator::readOmegaForceStrengthsFromDiskPREC() {
  if (LogLevel>3) printf("Reading PREC Omega Forces from Disk now...\n");
  threadedExecuteFORCE(pHMCPropagator_READOMEGAFORCESTRENGTHFROMDISK, 0);
}


void pHMCPropagator::readOmegaForceStrengthsFromDiskGLOBAL() {
  if (LogLevel>3) printf("Reading GLOBAL Omega Forces from Disk now...\n");
  threadedExecuteFORCE(pHMCPropagator_READOMEGAFORCESTRENGTHFROMDISK, 1);
}


void pHMCPropagator::setOmegaMassAdaptionMode(int omMassAdapMode) {
  if (LogLevel>4) printf("Setting Omega Mass Adaption Mode...\n");
  clearForceStorage();
  threadedExecuteFORCE(pHMCPropagator_SETOMEGAMASSADAPTIONMODE, omMassAdapMode);
}


void pHMCPropagator::calcOmegaMassAdaption() {
  if (LogLevel>3) printf("Calculating Omega Masses...\n");
  threadedExecuteFORCE(pHMCPropagator_CALCOMEGAMASSADAPTION, 0);
}
  

void pHMCPropagator::analyzeAllPolynomialForces() {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
  double* dSdPhi = new double[4*L0*L1*L2*L3];

  threadedExecute(pHMCPropagator_SETPHI,0);
  calcFullPhiDerivatives(dSdPhi, 0, 3); //calculates also omegaForces
  for (int I=1; I<1+2*subPolyCount; I++) {
    calcFullPhiDerivatives(I, 2); //calculates also omegaForces  
  }
  delete[] dSdPhi;
}


void pHMCPropagator::smoothStartPhiConfiguration() {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  for (int I=0; I<L0*L1*L2*L3; I++) {
//    double w = (getSyncSingleRandom()-0.5);
  
/*    phiField[I][0] = sin(w)*sqrt(Nf);
    phiField[I][1] = cos(w)*sqrt(Nf);
    phiField[I][2] = 0;
    phiField[I][3] = 0;*/
    phiField[I][0] = sqrt(Nf)*getSyncSingleRandom();
    phiField[I][1] = sqrt(Nf)*getSyncSingleRandom();
    phiField[I][2] = sqrt(Nf)*getSyncSingleRandom();
    phiField[I][3] = sqrt(Nf)*getSyncSingleRandom();
  }
}


void pHMCPropagator::activateForceStoring(bool act) {
  threadedExecute(pHMCPropagator_ACTIVATEFORCESTORING,(double) act);
}


void pHMCPropagator::phiFieldWasChanged() {
  threadedExecute(pHMCPropagator_CLEARACTSTORE, 0);
}


void pHMCPropagator::PreconditionerWasChanged() {
  clearForceStorage();
}


void pHMCPropagator::clearForceStorage() {
  threadedExecute(pHMCPropagator_CLEARACTSTORE, 0);
  threadedExecute(pHMCPropagator_CLEARTRASTARTSTORE, 0);
}


void pHMCPropagator::changeOfFACCtype() {
  clearForceStorage();
}


void pHMCPropagator::synchronizedChangeOfQPreconsitionerData(bool useqPrecon, double mu, double beta) {
  fermiOps->setQPreconditioner(useqPrecon, mu, beta);
  threadedExecute(pHMCPropagator_SYNCQPRECDATA, 0);
  clearForceStorage();
}


void pHMCPropagator::synchronizedChangeOfRPreconsitionerData(bool userPrecon, double m, double f) {
  fermiOps->setRPreconditioner(userPrecon, m, f);
  threadedExecute(pHMCPropagator_SYNCRPRECDATA, 0);
  clearForceStorage();
}


void pHMCPropagator::resetExactMMdagInverseSQRTOmegaAction() {
  threadedExecute(pHMCPropagator_RESETEXACTMMDAGINVERSESQRTOMEGAACTION, 0);
}


void pHMCPropagator::calcExactMMdagInverseSQRTOmegaAction() {
  threadedExecute(pHMCPropagator_CALCEXACTMMDAGINVERSESQRTOMEGAACTION, 0);
}


double pHMCPropagator::getExactReweighingFactorFromMMdagInverseSQRTOmegaAction() {
  double weight = 1.0;

  if (LogLevel>3) printf("Determining exact reweighing factor for MMdagInverseSQRT-Omega-Action...\n ");
  if (fermiOps->getYN()>0) {
    bool Savail = true;
    for (int I=0; I<forceCount; I++) {
      double dummy = pHMCforces[I]->getActOmegaAction(0);
      if (isNaN(dummy)) Savail = false;
    }
    if (!Savail) {
      double para = 0+0.1;
      threadedExecute(pHMCPropagator_CALCOMEGAACTION, para);
    } else {
      if (LogLevel>4) printf("Taking Action from Data-Pool for polynom %d\n",0);
    }
    double polS = 0;
    for (int I=0; I<forceCount; I++) {
      polS = polS + pHMCforces[I]->getActOmegaAction(0);
    }

    Savail = true;
    for (int I=0; I<forceCount; I++) {
      double dummy = pHMCforces[I]->getExactOmegaMMdagInverseSQRTAction();
      if (isNaN(dummy)) Savail = false;
    }
    if (!Savail) {
      threadedExecute(pHMCPropagator_CALCEXACTMMDAGINVERSESQRTOMEGAACTION, 0);
    } else {
      if (LogLevel>4) printf("Taking over exact MMdag-Inverse-SQRT omega-Action\n");
    }
    double exactS = 0;
    for (int I=0; I<forceCount; I++) {
      exactS = exactS + pHMCforces[I]->getExactOmegaMMdagInverseSQRTAction();
    }

    if (LogLevel>3) printf("  ==> pol. S = %1.8f, exact S = %1.8f ",polS,exactS);

    double deltaS = exactS-polS;
    weight = exp(-deltaS);
  } else {
    weight = 1.0;
  }
  if (LogLevel>3) printf("  ==> weight = %1.5f\n",weight);

  return weight;
}
  
  
void pHMCPropagator::synchronizedChangeOfQuasiHermiteanMode(bool qHM) {
  threadedExecute(pHMCPropagator_SETQUASIHERMITEANMODE, (int) qHM);
}


void pHMCPropagator::synchronizedChangeOfModelSelection(int modelSel) {
  threadedExecute(pHMCPropagator_SETMODELSELECTION, modelSel);
}


void pHMCPropagator::synchronizedChangeOfTuneMode(bool tM) {
  threadedExecute(pHMCPropagator_SETTUNEMODE, (int) tM);
}


void pHMCPropagator::saveUpperEWboundLogToDisk(char* fileName) {
  if (LogLevel>3) printf("Saving upper EW-bound log to file %s.\n", fileName);
  FILE* file = fopen(fileName, "w");

  for (int I=0; I<upperEWboundLogCount; I++) {
    fprintf(file, "%1.0f %1.15f %1.15f %1.15f\n", upperEWboundLog[I][0], upperEWboundLog[I][1], upperEWboundLog[I][2], upperEWboundLog[I][3]);
  }
  fclose(file);
}


void pHMCPropagator::loadUpperEWboundLogFromDisk(char* fileName) {
  if (LogLevel>2) printf("Loading upper EW-bound log from file %s...\n", fileName);
  upperEWboundLogCount = 0;
  FILE* file = fopen(fileName, "r");
  if (file==NULL) {  
    if (LogLevel>2) printf("File not found.\n");
    return;
  }

  while (fscanf(file, "%lf %lf %lf %lf\n", &(upperEWboundLog[upperEWboundLogCount][0]), &(upperEWboundLog[upperEWboundLogCount][1]), &(upperEWboundLog[upperEWboundLogCount][2]), &(upperEWboundLog[upperEWboundLogCount][3]))==4) {
    upperEWboundLogCount++;
  }
  fclose(file);
  if (LogLevel>2) printf("%d data sets loaded.\n", upperEWboundLogCount);
}


pHMCForce* pHMCPropagator::getForce(int forceNr) {
  if ((forceNr<0) || (forceNr>=forceCount)) return NULL;
  return pHMCforces[forceNr];
}
