#include "Tools.h"


//Variablen
ComplexMatrix ID4x4(4);
ComplexMatrix Gamma5(4);
bool numaAvail = false;
int numaNodesCount = 1;
int numaCoresCount = 1;
bool numaSupportInitialized = false;
int numaNodeAllocationMode = 1;
PerformanceProfiler* performanceProfiler = NULL;
char* GeneralUniqueFileNameExtension = NULL;
double findZeroOfBosonic1LoopInvPropagatorFit_Helper_m0 = 0;
double findZeroOfBosonic1LoopInvPropagatorFit_Helper_Z = 0;
int findZeroOfBosonic1LoopInvPropagatorFit_Helper_N = 0;
double* findZeroOfBosonic1LoopInvPropagatorFit_Helper_coeff = NULL;
int AdvancedSeed; 
double eps;
double NaN;
double epsilon4[4][4][4][4];
Quat sigmaPlus[4];            //Der Einfachheit halber sigma0 = sigma4
Quat sigmaMinus[4];           //Der Einfachheit halber sigma0 = sigma4
Complex* ASMParameterData; 
int LuescherZetaFunctionHelper2_NvecSqr;
int LuescherDerivativeOfZetaFunction_dZdqSqrHelper2_NvecSqr;
double EffectiveMassGradientMinimizationSolverHelper_relTime0;
double EffectiveMassGradientMinimizationSolverHelper_relTime1;
double EffectiveMassGradientMinimizationSolverHelper_QuotValue;

int getCurentThreadID() {
  pthread_t thread = pthread_self();
  void** p = (void**) thread;
  long int id = (long int) (p[18]);
  return (int) id;
}

void setCurrentThreadAffinityMask(long int mask) {
  int id = getCurentThreadID();
  if (LogLevel>1) printf("Set Affinity-Mask for thread (%d) to mask %ld\n", id, mask);
  int res = sched_setaffinity((pid_t) id, sizeof(mask),(cpu_set_t*) &mask);
  if (res !=0) {
    printf("ERROR in setCurrentThreadAffinityMask!\n");
    exit(0);
  }
}

long int getAffinityMaskFromNodeID(int nodeID) {  
  if (!numaSupportInitialized) {
    printf("ERROR in getAffinityMaskFromNodeID: Numa-Support not initialized!!!\n");
    exit(0);
  }
  int coresPerNode = numaCoresCount / numaNodesCount;
  long int mask = 0;
  for (int I2=0; I2<coresPerNode; I2++) {
    mask = mask | getAffinityMaskFromNthCoreOnNodeID(nodeID, I2);
  }
  return mask;
}


long int getAffinityMaskFromCoreID(int coreID) {
  if (!numaSupportInitialized) {
    printf("ERROR in getAffinityMaskFromCoreID: Numa-Support not initialized!!!\n");
    exit(0);
  }
  long int mask = 1;
  for (int I=0; I<(coreID%numaCoresCount); I++) {
    mask *= 2;
  }
  return mask;
}


long int getAffinityMaskFromNthCoreOnNodeID(int nodeID, int nthCore) {  
  if (!numaSupportInitialized) {
    printf("ERROR in getAffinityMaskFromNthCoreOnNodeID: Numa-Support not initialized!!!\n");
    exit(0);
  }
  int coresPerNode = numaCoresCount / numaNodesCount;

  int coreID = 0;
  if (numaNodeAllocationMode == 1) {
    coreID = (nodeID%numaNodesCount)*coresPerNode + (nthCore%coresPerNode);	 // (0,1) (2,3) (4,5) (6,7) etc...
  } else {
    coreID = (nodeID%numaNodesCount) + (nthCore%coresPerNode) * numaNodesCount;  // (0,4) (1,5) (2,6) (3,7) etc...
  }
  
  long int mask = getAffinityMaskFromCoreID(coreID);
  return mask;
}

void setCurrentThreadAffinityToCoreID(int coreID) {
  long int mask = getAffinityMaskFromCoreID(coreID);
  setCurrentThreadAffinityMask(mask);
}


void setCurrentThreadAffinityToNthCoreOnNodeID(int nodeID, int nthCore) {
  long int mask = getAffinityMaskFromNthCoreOnNodeID(nodeID, nthCore);
  setCurrentThreadAffinityMask(mask);
}


void setCurrentThreadAffinityToNodeID(int nodeID) {
  long int mask = getAffinityMaskFromNodeID(nodeID);
  setCurrentThreadAffinityMask(mask);
}


void initializeNUMASupport() {
  if (numaSupportInitialized) return;
  numaCoresCount = sysconf(_SC_NPROCESSORS_CONF);
  if (LogLevel>1) printf("Number of CPU-cores available: %d\n", numaCoresCount);
  numaAvail = (numa_available() >= 0);
  numaNodesCount = 1;
  numaNodeAllocationMode = 1;
  if (numaAvail) {
    numaNodesCount = 1+numa_max_node();
    if (numaNodesCount>=4) numaNodeAllocationMode = 2;    
    if (LogLevel>1) { 
      printf("NUMA is available with %d nodes with ", numaNodesCount);
      for (int I=0; I<numaNodesCount-1; I++) {
        printf("%1.2f GB, ",numa_node_size(I, NULL)/(1024.0*1024.0*1024.0));	
      }
      printf("%1.2f GB of node memory.\n",numa_node_size(numaNodesCount-1, NULL)/(1024.0*1024.0*1024.0));	
      printf("==> %d cores per node\n", numaCoresCount / numaNodesCount);
      printf("==> Numa-Node-Allocation mode = %d\n", numaNodeAllocationMode);      
    }
  } else {  
    if (LogLevel>1) printf("NUMA is NOT available!!!\n");  
  }
  numaSupportInitialized = true;
  if (LogLevel>1) { 
    printf("==> CPU-Topology: ");
    for (int I=0; I<numaNodesCount; I++) {
      printf("(%ld", getAffinityMaskFromNthCoreOnNodeID(I, 0));
      for (int I2=1; I2<numaCoresCount/numaNodesCount; I2++) {
        printf(", %ld", getAffinityMaskFromNthCoreOnNodeID(I, I2));      
      }
      printf(") ");    
    }
    printf("\n");  
  }  
}


//ALIGN muss Zweierpotenz sein
Complex* createSuperAlignedComplex(int size, int ALIGN) {
  char* p = new char[16*size + 2*ALIGN];
  long int i = (long int) (p);
  i = i & (ALIGN - 1);
  Complex* p2 = (Complex*)&(p[2*ALIGN-i]);
  char** p3 = (char**)(p2 - sizeof(char*));
  *p3 = (char*)p;
  return p2;
}


Complex* createSuperAlignedComplex(int size) {
  return createSuperAlignedComplex(size,128);
}


void destroySuperAlignedComplex(Complex* &p) {
  if (p == NULL) return;
  char** p2 = (char**)(p - sizeof(char*));
  delete[] (*p2);
  p = NULL;
}


void destroySuperAlignedIntP(int** &p) {
  if (p == NULL) return;
  char** p2 = (char**)(p - sizeof(char*));
  delete[] (*p2);
  p = NULL;
}


void destroySuperAlignedLongInt(long int* &p) {
  Complex* pC = (Complex*) p;
  destroySuperAlignedComplex(pC);
  p = NULL;
}


void destroySuperAlignedInt(int* &p) {
  Complex* pC = (Complex*) p;
  destroySuperAlignedComplex(pC);
  p = NULL;
}


void initializePerformanceProfiler(char* fileName) {
  if (performanceProfiler != NULL) return;
  performanceProfiler = new PerformanceProfiler(fileName);
}

void writePerformanceProfilingDataToDisk() {
  if (performanceProfiler == NULL) return;
  performanceProfiler->writePerformanceItemsToDisk();
}

void printBits(long int x) {
  int I;
  int bits[64];
  for (I=0; I<64; I++) {
    bits[I] = x & 1;
    x = x / 2;
  }
  for (I=0; I<64; I++) {
    printf("%d",bits[63-I]);
  }  
  printf("\n");
}


long int getCPUCycleCounter() {
  long int val = 0;
  unsigned int __a,__d; 
  asm volatile("rdtsc" : "=a" (__a), "=d" (__d)); 
  val = ((unsigned long)__a) | (((unsigned long)__d)<<32); 
  return val;
} 

void calcEPS() {
  int I;
  double h=1.0;
  double d;

  eps = 1.0;
  for (I=0; I<20; I++) {
    d = 1.0;
    d += h;
    if (d == 1.0) { return; }

    eps = h;
    h /= 10;
  }
}


double zeitwert() {
  timeval tv;
  gettimeofday(&tv, NULL);
  double secs = tv.tv_sec + (tv.tv_usec/1E6);

  return secs;
}


void delay(double timeInSecs) {
  double startTime = zeitwert();
  while (zeitwert()-startTime < timeInSecs) {
  
  }
}


double cpuTime() {
  return (1.0*clock()) / CLOCKS_PER_SEC;
}


void determineCPUCyclesPerSecond() {
  double start = zeitwert();
  long int CPUcycleStart = getCPUCycleCounter();
  bool rep = true;
  double end = 0;
  while (rep) {
    for (int I=0; I<10000000; I++) {
    
    }
    end = zeitwert();
    if (end-start > 2) rep = false;
  }
  long int CPUcycleEnd = getCPUCycleCounter();
  CPUCyclesPerSecond = (long int) ((CPUcycleEnd-CPUcycleStart) / (end-start));
}



void randomize() {
  srand(time(NULL));
}


double zufall() {
  double d = rand();
  return d / RAND_MAX;
}

double round(double x) {
  int f = (int) floor(x);
  int c = f+1;
  if ((x-f)>(c-x)) {
    return c;
  } else {
    return f;
  }
}


int roundToInt(double x) {
  int f = (int) floor(x);
  int c = f+1;
  if ((x-f)>(c-x)) {
    return c;
  } else {
    return f;
  }
}

/*
inline double abs(double x) {
  if (x>=0) return x;
  return -x;
}
*/

bool isInteger(double x) {
  if (x==0) return true;
  int r = roundToInt(x);
  if ((abs(r-x)/r)<1E-14) return true;
  return false;
}


double makeNaN() {
  double x = 0.0;
  return 0.0/x;
}


bool isNaN(double x) {
  if ((x <= 0.0) || (x >= 0.0)) { return 0; }
  return 1;
}

void print(double d) {
  printf("%1.15f\n",d);
}


void print(vector4D v) {
  printf("%1.15f\n",v[0]);
  printf("%1.15f\n",v[1]);
  printf("%1.15f\n",v[2]);
  printf("%1.15f\n",v[3]);
}


void print(int i) {
  printf("%d\n",i);
}


int calcPermutationSignum(int *per, int N) {
  int I;
  bool change = true;
  int dummy;
  int signum = 1;

  while (change) {
    change = false;
    for (I=0; I<N-1; I++) {
      if (per[I] == per[I+1]) { return 0; }
      if (per[I] > per[I+1]) {
        signum *= -1;
        dummy = per[I];
        per[I] = per[I+1];
        per[I+1] = dummy;
        change = true;
      }
    }
  }
  return signum;
}


void iniEpsilon4() {
  int *per = new int[4];
  int I1,I2,I3,I4;

  for (I1=0; I1<4; I1++) {
    for (I2=0; I2<4; I2++) {
      for (I3=0; I3<4; I3++) {
        for (I4=0; I4<4; I4++) {
          per[0] = I1;
          per[1] = I2;
          per[2] = I3;
          per[3] = I4;

          epsilon4[I1][I2][I3][I4] = calcPermutationSignum(per,4);
        }
      }
    }
  }
  delete[] per;
}


void iniTools(int randGenNum) {
  if (LogLevel>0) printf("Tools initializing...\n");

  randomize();
  AdvancedSeed = rand()*randGenNum;
  if (AdvancedSeed>0) AdvancedSeed = -AdvancedSeed;
  if (LogLevel>0) printf("Start advanced random number generator with seed: %d (fac was %d).\n",AdvancedSeed,randGenNum);
  AdvancedZufall(AdvancedSeed);
  if (LogLevel>2) {
    printf("First 5 random numbers...\n");
    double z1 = AdvancedZufall(AdvancedSeed);
    double z2 = AdvancedZufall(AdvancedSeed);
    double z3 = AdvancedZufall(AdvancedSeed);
    double z4 = AdvancedZufall(AdvancedSeed);
    double z5 = AdvancedZufall(AdvancedSeed);
    printf("%1.5f\n",z1);
    printf("%1.5f\n",z2);
    printf("%1.5f\n",z3);
    printf("%1.5f\n",z4);
    printf("%1.5f\n",z5);
  }

  NaN = makeNaN();
  calcEPS();
  if (LogLevel>0) printf("EPS = %1.2f x 10^-15\n",1E+15*eps);

  ASMParameterData = createSuperAlignedComplex(100000);

  ComplexI.setValues(0,1);
  ComplexUnity.setValues(1,0);
  ComplexZero.setValues(0,0);

  sigmaPlus[0].setValues(1,0,0,0);
  sigmaPlus[1].setValues(0,-1,0,0);
  sigmaPlus[2].setValues(0,0,-1,0);
  sigmaPlus[3].setValues(0,0,0,-1);

  sigmaMinus[0].setValues(1,0,0,0);
  sigmaMinus[1].setValues(0,1,0,0);
  sigmaMinus[2].setValues(0,0,1,0);
  sigmaMinus[3].setValues(0,0,0,1);

  iniEpsilon4();
  
  ID4x4.matrix[0][0] = Complex(1,0);
  ID4x4.matrix[1][1] = Complex(1,0);
  ID4x4.matrix[2][2] = Complex(1,0);
  ID4x4.matrix[3][3] = Complex(1,0);
  
  Gamma5.matrix[0][0] = Complex(1,0);
  Gamma5.matrix[1][1] = Complex(1,0);
  Gamma5.matrix[2][2] = Complex(-1,0);
  Gamma5.matrix[3][3] = Complex(-1,0);
  
  initializeNUMASupport();
  determineCPUCyclesPerSecond();
  if (LogLevel>0) printf("CPU-Cycles per second: %1.2f GCycles/s\n", CPUCyclesPerSecond/1E9);  
  determineCPUCyclesPerSecond();
  if (LogLevel>0) printf("CPU-Cycles per second: %1.2f GCycles/s (rechecked)\n", CPUCyclesPerSecond/1E9);  
  
  GeneralUniqueFileNameExtension = new char[100];
  snprintf(GeneralUniqueFileNameExtension,100,"%d", (int)(1000000*AdvancedZufall(AdvancedSeed)));
  if (LogLevel>0) printf("General Unique Filename Extension: %s\n", GeneralUniqueFileNameExtension);
  
  if (LogLevel>0) printf("...sucessfully.\n");
}


void desiniTools() {
  destroySuperAlignedComplex(ASMParameterData);
  delete performanceProfiler;
  delete[] GeneralUniqueFileNameExtension;
}


double calcLogDetScaledAbsNorm(ComplexMatrix& mat, int removeSmallestModesCount, double scaleFac, bool execCalcEigenV) {
  if (execCalcEigenV) {
    mat.calcEigenvalues();
  }
  int I,I2;
  
  for (I2=0; I2<removeSmallestModesCount; I2++) {
    int IMerker = 0;
    Complex nu = mat.eigenvalues[0];    
    double nuNorm = sqrt(sqr(nu.x) + sqr(nu.y))/scaleFac;    
    double nuMerker = nuNorm;    
    for (I=0; I<mat.matrixSize-I2; I++) {
      nu = mat.eigenvalues[I];
      nuNorm = sqrt(sqr(nu.x) + sqr(nu.y))/scaleFac;    
      
      if (nuNorm<nuMerker) {
        IMerker = I;
	nuMerker = nuNorm;
      }
    }
    Complex interchange = mat.eigenvalues[mat.matrixSize-I2-1];
    mat.eigenvalues[mat.matrixSize-I2-1] = mat.eigenvalues[IMerker];
    mat.eigenvalues[IMerker] = interchange;
    if (LogLevel>2) {
      printf("Sorting out eigenvalue: ");
      mat.eigenvalues[mat.matrixSize-I2-1].print();
    }
  }
  double logDeterminant = 0;
  for (I=0; I<mat.matrixSize-removeSmallestModesCount; I++) {
    Complex nu = mat.eigenvalues[I];
    double nuNorm = sqrt(sqr(nu.x) + sqr(nu.y))/scaleFac;
    
    logDeterminant += log(nuNorm);    
//    if (abs((double)(logDeterminant/log(10)))>1500) {
//      printf("ERROR: Determinant: (%lf) exceeds maximal scale!!!\n",logDeterminant/log(10));
//      exit(0);
//    }
//    if (abs((double)(logDeterminant/log(10)))<-1500) {
//      printf("ERROR: Determinant: (%lf) exceeds minimal scale!!!\n",logDeterminant/log(10));
//      exit(0);
//    }
  } 
  
  return logDeterminant;
}


void printEigenvalues(char* OPname, ComplexMatrix& op, Complex fac) {
  printf("Trying to print eigenvalues of %s-operator to disk...",OPname);
  bool b = op.calcEigenvalues();
  if (b) {
    char* fileName = new char[500];
    snprintf(fileName,500,"%sEigen.dat",OPname);
    FILE* file;
    file = fopen(fileName,"w");
    int I;
    for (I=0; I<op.matrixSize; I++) {
      fprintf(file,"%1.15f %1.15f \n", (fac*op.eigenvalues[I]).x, (fac*op.eigenvalues[I]).y);
    }
    fclose(file);    
    delete[] fileName;
    printf("...sucessfully.\n");  
  } else {
    printf("...ERROR!!!\n");  
    exit(0);
  }
}


void derive(double (*func)(double* x), int N, double EPS, double* x, double* d, int gradientMask) {
  int i;
  double dummy;
  double res;
  int m = 1;
  for (i=0; i<N; i++) {
    if (((gradientMask & m) > 0) || (gradientMask<0)) {
      dummy = x[i];
      x[i] += EPS/2;
      res = (*func) (x);
      x[i] -= EPS;
      res -= (*func) (x);
      x[i] = dummy;
      res /= EPS;
      d[i] = res;
    } else {
      d[i] = 0;
    }
    m *= 2;
  }
}

//Normalizes a given vector to 1.
int normalizeVec(double* x, int N) {
  int i;
  double l = 0;
  for (i=0; i<N; i++) {
    l += x[i] * x[i];
  }
  l = sqrt(l);
  if (l == 0) {
    return 1;
  }
  for (i=0; i<N; i++) {
    x[i] /= l;
  }
  return 0;
}


//Verifies whether a point lies within a given boundary box given by b1, b2.
//Additionally sum<bSum is checked.
int insideBox(int N, double* x, double* b1, double* b2, double* bSUM) {
  int I;
  if (b1 != NULL) {
    for (I=0; I<N; I++) {
      if (x[I] < b1[I]) {
        return false;
      }
    }
  }
  if (b2 != NULL) {
    for (I=0; I<N; I++) {
      if (x[I] > b2[I]) {
        return false;
      }
    }
  }
  if (bSUM != NULL) {
    double sum = 0.0;
    for (I=0; I<N; I++) {
      sum += x[I];
    }
    if (sum > bSUM[0]) {
      return false;
    }
  }
  return true;
}


//Implements the minimum-search. Output: pos. Finds minimum of given function, using given stepsizes, stepsize for differentiation and bounds.
bool GradientMinimization(double (*func)(double* x), int N, double StartStepSize, double MinStepSize, double DiffEPS, double* pos, double* bounds1, double* bounds2, double* boundSUM, int gradientMask, int iterMax) {
  double* derived = new double[N];
  double f1;
  double f2;
  int I;
  int loopCount;
  double resetStepSize = StartStepSize;
  double stepAv;
  double step;
  double dummy;
  int flag;
  int iterCount = 0;

  while (true) {
    if (iterCount>iterMax) {
      return true;
    }
    iterCount++;
    step = resetStepSize;
    derive(func, N, DiffEPS, pos, derived, gradientMask);
    if (normalizeVec(derived, N) == 1) {
      dummy = DiffEPS;
      flag = false;
      if (LogLevel > 2) {
        fprintf(stderr,"Derivation equals zero while minimizing --> RECALC\n\n");
      }
      for (I=0; I<10; I++) {
        dummy *= 5.0;
        derive(func, N, dummy, pos, derived, gradientMask);
        if (normalizeVec(derived, N) == 0) {
	  flag = true;
          break;
        }
      }
      if (flag == false) {
        if (LogLevel>2) fprintf(stderr,"Derivation equals zero while minimizing --> Recalc unsuccessful!!!\n\n");
        return false;
      }
    }
    f1 = (*func)(pos);
    for (I=0; I<N; I++) {
      pos[I] -= step * derived[I];
    }
    f2 = (*func)(pos);
    stepAv = 0;
    loopCount = 0;
    while ((f1>f2) && (insideBox(N,pos,bounds1,bounds2,boundSUM) == true)) {
      loopCount++;
      stepAv += step;
      step *= 2.0;
      for (I=0; I<N; I++) {
        pos[I] -= step * derived[I];
      }
      f1 = f2;
      f2 = (*func)(pos);
    }
    for (I=0; I<N; I++) {
      pos[I] += step * derived[I];
    }
    if (loopCount == 0) {
      if (resetStepSize == MinStepSize) {
        delete[] derived;
        return true;
      }
      resetStepSize /= 1000;
    } else {
      resetStepSize = (stepAv / loopCount) / 10;
    }
    if (resetStepSize < MinStepSize) {
      resetStepSize = MinStepSize;
    }
  }
}

void perform_yB(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double y, double* phi, Complex* input, Complex* output) {
  int i0,i1,i2,i3;
  Complex* working = new Complex[8];
  int count = 0;
  int phiCount = 0;
  int xtrAdd1 = xtrS1*8*(L3+xtrS3)*(L2+xtrS2);
  int xtrAdd2 = xtrS2*8*(L3+xtrS3);
  int xtrAdd3 = xtrS3*8;
  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int I2=0; I2<8; I2++) {
            working[I2].x = input[count+I2].x;
            working[I2].y = input[count+I2].y;
          }	
	
          mulWithPhiMatB(&(phi[phiCount]), working, (Complex*)&(output[count].x), y);

          count += 8;
          phiCount += 4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  delete[] working;
}


void perform_yBsplit(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double y, double split, double* phi, Complex* input, Complex* output) {
  int i0,i1,i2,i3;
  Complex* working = new Complex[8];
  int count = 0;
  int phiCount = 0;
  int xtrAdd1 = xtrS1*8*(L3+xtrS3)*(L2+xtrS2);
  int xtrAdd2 = xtrS2*8*(L3+xtrS3);
  int xtrAdd3 = xtrS3*8;
  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int I2=0; I2<8; I2++) {
            working[I2].x = input[count+I2].x;
            working[I2].y = input[count+I2].y;
          }	
	
          mulWithPhiMatB(&(phi[phiCount]), working, (Complex*)&(output[count].x), split, y);

          count += 8;
          phiCount += 4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  delete[] working;
}


void performf_YBD_2rhoD_2rhoD_AndScalarProducts(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double twoRho, double explicitMass, Complex* x, Complex* Dx, double* phi, Complex* output, Complex* Vrest, Complex* Vp, Complex& resVrests, Complex& resVps) {
  int I2;
  int i0,i1,i2,i3;
  Complex* working = new Complex[8];
  Complex* interimSource = new Complex[8];
  int count = 0;
  int phiCount = 0;
  resVrests = Complex(0,0);
  resVps = Complex(0,0);
  int xtrAdd1 = xtrS1*8*(L3+xtrS3)*(L2+xtrS2);
  int xtrAdd2 = xtrS2*8*(L3+xtrS3);
  int xtrAdd3 = xtrS3*8;
  double twoRhoMass = twoRho * explicitMass;
  
  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (I2=0; I2<8; I2++) {
            working[I2].x = Dx[count+I2].x - twoRho * x[count+I2].x;
            working[I2].y = Dx[count+I2].y - twoRho * x[count+I2].y;
            interimSource[I2].x = twoRhoMass * x[count+I2].x;
            interimSource[I2].y = twoRhoMass * x[count+I2].y;
          }
    
          mulWithPhiMatB(&(phi[phiCount]), working, (Complex*)&(output[count].x), yN);
    
          for (I2=0; I2<8; I2++) {
            output[count+I2].x -= twoRho * Dx[count+I2].x + interimSource[I2].x;
            output[count+I2].y -= twoRho * Dx[count+I2].y + interimSource[I2].y;
          }
    
          for (I2=0; I2<8; I2++) {
            resVps.x += Vp[count+I2].x*output[count+I2].x + Vp[count+I2].y*output[count+I2].y;
            resVps.y += Vp[count+I2].x*output[count+I2].y - Vp[count+I2].y*output[count+I2].x;
            resVrests.x += Vp[count+I2].x*Vrest[count+I2].x + Vp[count+I2].y*Vrest[count+I2].y;
            resVrests.y += Vp[count+I2].x*Vrest[count+I2].y - Vp[count+I2].y*Vrest[count+I2].x;
          }    
  
          count += 8;
          phiCount += 4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  delete[] working;
  delete[] interimSource;
}


void performf_YBsplitD_2rhoD_2rhoD_AndScalarProducts(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double split, double twoRho, double explicitMass, Complex* x, Complex* Dx, double* phi, Complex* output, Complex* Vrest, Complex* Vp, Complex& resVrests, Complex& resVps) {
  int I2;
  int i0,i1,i2,i3;
  Complex* working = new Complex[8];
  Complex* interimSource = new Complex[8];
  int count = 0;
  int phiCount = 0;
  resVrests = Complex(0,0);
  resVps = Complex(0,0);
  int xtrAdd1 = xtrS1*8*(L3+xtrS3)*(L2+xtrS2);
  int xtrAdd2 = xtrS2*8*(L3+xtrS3);
  int xtrAdd3 = xtrS3*8;
  double twoRhoMass = twoRho * explicitMass;
  
  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (I2=0; I2<8; I2++) {
            working[I2].x = Dx[count+I2].x - twoRho * x[count+I2].x;
            working[I2].y = Dx[count+I2].y - twoRho * x[count+I2].y;
            interimSource[I2].x = twoRhoMass * x[count+I2].x;
            interimSource[I2].y = twoRhoMass * x[count+I2].y;
          }
    
          mulWithPhiMatB(&(phi[phiCount]), working, (Complex*)&(output[count].x), split, yN);
    
          for (I2=0; I2<8; I2++) {
            output[count+I2].x -= twoRho * Dx[count+I2].x + interimSource[I2].x;
            output[count+I2].y -= twoRho * Dx[count+I2].y + interimSource[I2].y;
          }
    
          for (I2=0; I2<8; I2++) {
            resVps.x += Vp[count+I2].x*output[count+I2].x + Vp[count+I2].y*output[count+I2].y;
            resVps.y += Vp[count+I2].x*output[count+I2].y - Vp[count+I2].y*output[count+I2].x;
            resVrests.x += Vp[count+I2].x*Vrest[count+I2].x + Vp[count+I2].y*Vrest[count+I2].y;
            resVrests.y += Vp[count+I2].x*Vrest[count+I2].y - Vp[count+I2].y*Vrest[count+I2].x;
          }    
  
          count += 8;
          phiCount += 4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  delete[] working;
  delete[] interimSource;
}


void performf_YBD_2rhoD_2rhoD(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double twoRho,
double explicitMass, Complex* x, Complex* Dx, double* phi, Complex* output) {
  int I2;
  int i0,i1,i2,i3;
  Complex* working = new Complex[8];
  Complex* interimSource = new Complex[8];
  int count = 0;
  int phiCount = 0;
  int xtrAdd1 = xtrS1*8*(L3+xtrS3)*(L2+xtrS2);
  int xtrAdd2 = xtrS2*8*(L3+xtrS3);
  int xtrAdd3 = xtrS3*8;
  double twoRhoMass = twoRho * explicitMass;
  
  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {
  
          for (I2=0; I2<8; I2++) {
            working[I2].x = Dx[count+I2].x - twoRho * x[count+I2].x;
            working[I2].y = Dx[count+I2].y - twoRho * x[count+I2].y;
            interimSource[I2].x = twoRhoMass * x[count+I2].x;
            interimSource[I2].y = twoRhoMass * x[count+I2].y;
          }
    
          mulWithPhiMatB(&(phi[phiCount]), working, (Complex*)&(output[count].x), yN);
    
          for (I2=0; I2<8; I2++) {
            output[count+I2].x -= twoRho * Dx[count+I2].x + interimSource[I2].x;
            output[count+I2].y -= twoRho * Dx[count+I2].y + interimSource[I2].y;
          }
  
          count += 8;
          phiCount += 4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  delete[] working;
  delete[] interimSource;
}


void performf_YBsplitD_2rhoD_2rhoD(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double split, double twoRho, double explicitMass, Complex* x, Complex* Dx, double* phi, Complex* output) {
  int I2;
  int i0,i1,i2,i3;
  Complex* working = new Complex[8];
  Complex* interimSource = new Complex[8];
  int count = 0;
  int phiCount = 0;
  int xtrAdd1 = xtrS1*8*(L3+xtrS3)*(L2+xtrS2);
  int xtrAdd2 = xtrS2*8*(L3+xtrS3);
  int xtrAdd3 = xtrS3*8;
  double twoRhoMass = twoRho * explicitMass;
  
  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {
  
          for (I2=0; I2<8; I2++) {
            working[I2].x = Dx[count+I2].x - twoRho * x[count+I2].x;
            working[I2].y = Dx[count+I2].y - twoRho * x[count+I2].y;
            interimSource[I2].x = twoRhoMass * x[count+I2].x;
            interimSource[I2].y = twoRhoMass * x[count+I2].y;
          }
    
          mulWithPhiMatB(&(phi[phiCount]), working, (Complex*)&(output[count].x), split, yN);
    
          for (I2=0; I2<8; I2++) {
            output[count+I2].x -= twoRho * Dx[count+I2].x + interimSource[I2].x;
            output[count+I2].y -= twoRho * Dx[count+I2].y + interimSource[I2].y;
          }
  
          count += 8;
          phiCount += 4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  delete[] working;
  delete[] interimSource;
}

void perform_yBDagger(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double y, double* phi, Complex* input, Complex* output)  {
  int i0,i1,i2,i3;
  Complex* working = new Complex[8];
  int count = 0;
  int phiCount = 0;
  int xtrAdd1 = xtrS1*8*(L3+xtrS3)*(L2+xtrS2);
  int xtrAdd2 = xtrS2*8*(L3+xtrS3);
  int xtrAdd3 = xtrS3*8;
  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int I2=0; I2<8; I2++) {
            working[I2].x = input[count+I2].x;
            working[I2].y = input[count+I2].y;
          }	
	
          mulWithPhiDaggeredMatB(&(phi[phiCount]), working, (Complex*)&(output[count].x), y);

          count += 8;
          phiCount += 4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  delete[] working;
}


void perform_yBsplitDagger(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double y, double split, double* phi, Complex* input, Complex* output)  {
  int i0,i1,i2,i3;
  Complex* working = new Complex[8];
  int count = 0;
  int phiCount = 0;
  int xtrAdd1 = xtrS1*8*(L3+xtrS3)*(L2+xtrS2);
  int xtrAdd2 = xtrS2*8*(L3+xtrS3);
  int xtrAdd3 = xtrS3*8;
  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {

          for (int I2=0; I2<8; I2++) {
            working[I2].x = input[count+I2].x;
            working[I2].y = input[count+I2].y;
          }	
	
          mulWithPhiDaggeredMatB(&(phi[phiCount]), working, (Complex*)&(output[count].x), split, y);

          count += 8;
          phiCount += 4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
  delete[] working;
}


void performf_YB_2rho_AndCopyToOutput(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double twoRho, Complex* x, double* phi, Complex* interim, Complex* output) {
  int i0,i1,i2,i3;
  int count = 0;
  int phiCount = 0;
  int xtrAdd1 = xtrS1*8*(L3+xtrS3)*(L2+xtrS2);
  int xtrAdd2 = xtrS2*8*(L3+xtrS3);
  int xtrAdd3 = xtrS3*8;
  double f = 1.0 / twoRho;
  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {
          mulWithPhiDaggeredMatB(&(phi[phiCount]), (Complex*)&(x[count]), (Complex*)&(interim[count]), yN);
    
          for (int I2=0; I2<8; I2++) {
            output[count+I2].x = interim[count+I2].x;
            output[count+I2].y = interim[count+I2].y;
            interim[count+I2].x = (f * interim[count+I2].x) - x[count+I2].x;
            interim[count+I2].y = (f * interim[count+I2].y) - x[count+I2].y;
          }

          count += 8;
          phiCount += 4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
}


void performf_YBsplit_2rho_AndCopyToOutput(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double split, double twoRho, double explicitMass, Complex* x, double* phi, Complex* interim, Complex* output) {
  int i0,i1,i2,i3;
  int count = 0;
  int phiCount = 0;
  int xtrAdd1 = xtrS1*8*(L3+xtrS3)*(L2+xtrS2);
  int xtrAdd2 = xtrS2*8*(L3+xtrS3);
  int xtrAdd3 = xtrS3*8;
  double f = 1.0 / twoRho;
  
  for (i0=0; i0<L0; i0++) {
    for (i1=0; i1<L1; i1++) {
      for (i2=0; i2<L2; i2++) {
        for (i3=0; i3<L3; i3++) {
          mulWithPhiDaggeredMatB(&(phi[phiCount]), (Complex*)&(x[count]), (Complex*)&(interim[count]), split, yN);
    
          for (int I2=0; I2<8; I2++) {
            output[count+I2].x = interim[count+I2].x + explicitMass*x[count+I2].x;
            output[count+I2].y = interim[count+I2].y + explicitMass*x[count+I2].y;
            interim[count+I2].x = (f * interim[count+I2].x) - x[count+I2].x;
            interim[count+I2].y = (f * interim[count+I2].y) - x[count+I2].y;
          }

          count += 8;
          phiCount += 4;
	}
	count += xtrAdd3;
      }
      count += xtrAdd2;
    }
    count += xtrAdd1;
  }
}


ComplexMatrix* createPhiMatB(vector4D phi, bool daggered, double split ) {
  ComplexMatrix* mat = new ComplexMatrix(8);
  mat->setZero();

  ComplexVector in(8);
  ComplexVector out(8);

  for (int i2=0; i2<8; i2++) {
    in.setZero();
    in.vectorElements[i2].x = 1.0;
    if (daggered) {
      mulWithPhiDaggeredMatB((double*) phi, in.vectorElements, out.vectorElements, split, 1.0);
    } else {
      mulWithPhiMatB((double*) phi, in.vectorElements, out.vectorElements, split, 1.0);
    }
    for (int i1=0; i1<8; i1++) {  
      mat->matrix[i1][i2] = out.vectorElements[i1];
    }
  }
  return mat;
}


ComplexMatrix* createPhiMatrix(vector4D phi, bool daggered) {
  double f = 1;
  if (daggered) f = -1;
  Quat q(phi[0],f*phi[1],f*phi[2],f*phi[3]);
  ComplexMatrix qm(q);

  ComplexMatrix* mat = new ComplexMatrix(8);
  mat->setZero();

  for (int I=0; I<4; I++) mat->matrix[I+0][I+0] = qm.matrix[0][0];
  for (int I=0; I<4; I++) mat->matrix[I+0][I+4] = qm.matrix[0][1];
  for (int I=0; I<4; I++) mat->matrix[I+4][I+0] = qm.matrix[1][0];
  for (int I=0; I<4; I++) mat->matrix[I+4][I+4] = qm.matrix[1][1];

  return mat;
}


bool performGnuplotFit(char* function, char* fitCommand, int varNr, int optiNr, double* res, double* err, double& redChiSqr) {
  char* gnuFileName = new char[2000];
  char* gnuCommand = new char[2000];
  char* gnuFitFileName = new char[2000];  
  char* rmgnuFile = new char[2000];  
  char* rmgnuFitFile = new char[2000];  
  snprintf(gnuFileName,2000,"%s/data/GnuplotFit%s.gnu",DataBaseDirectory,GeneralUniqueFileNameExtension);
  snprintf(gnuCommand,2000,"gnuplot %s/data/GnuplotFit%s.gnu",DataBaseDirectory,GeneralUniqueFileNameExtension);
  snprintf(gnuFitFileName,2000,"fit%s.log",GeneralUniqueFileNameExtension);
  snprintf(rmgnuFile,2000,"rm %s/data/GnuplotFit%s.gnu",DataBaseDirectory,GeneralUniqueFileNameExtension);
  snprintf(rmgnuFitFile,2000,"rm fit%s.log",GeneralUniqueFileNameExtension);
  
  FILE* file = fopen(gnuFileName,"w");
  int I;
  for (I=0; I<varNr; I++) {
    fprintf(file, "A%d = %1.15f\n", I+1, res[I]);
  }
  fprintf(file, "set fit logfile '%s'\n",gnuFitFileName);  
  fprintf(file, "%s\n%s\n",function,fitCommand);
  fclose(file);
  
  int SysError = system(gnuCommand);
  if (SysError<0) {
    printf("System Error: %d\n",SysError);
    exit(0);
  }

  char* s = new char[1000];  
  char* dummyStr = new char[1000];
  double dummy, dummy2;
  int resCount = 0;
  int count = 0;
  bool finalRead = false;
  redChiSqr = NaN;
  
  file = fopen(gnuFitFileName,"r");
  while (true) {  
  
    fgets(s, 1000, file);
    int len = strlen(s);
  
    if (sscanf(s,"rms of residuals  %s = %s : %lf", dummyStr, dummyStr, &dummy) == 3) {
      redChiSqr = dummy * dummy; 
    }
    if (sscanf(s,"Final set of %s", dummyStr) == 1) {
      finalRead = true;
    }

    if ((finalRead) && (len>16)) {
      for (I=0; I<=len; I++) {
        s[I] = s[I+16];    
      }
    }

    if ((finalRead) && (sscanf(s,"= %lf", &dummy) == 1)) {    
      if (len>35) {
        for (I=0; I<=len-16; I++) {
          s[I] = s[I+19];    
        }
      }
      if (sscanf(s,"+/- %lf", &dummy2) == 1) {    
        res[resCount] = dummy;
        err[resCount] = dummy2;
        resCount++;
        if (resCount>=optiNr) break;
      } else {
        res[resCount] = dummy;
        err[resCount] = NaN;
        resCount++;
        if (resCount>=optiNr) break;
      }
    }
    
    count++;
    if (count>1000000) break;
  }

  fclose(file);
  if (!DebugMode) system(rmgnuFile);    
  system(rmgnuFitFile);    

  delete[] gnuFileName;
  delete[] gnuCommand;
  delete[] gnuFitFileName;
  delete[] rmgnuFile;
  delete[] rmgnuFitFile;
  delete[] s;
  delete[] dummyStr;
  return (resCount>=optiNr);
}


bool performGnuplotFit(char* functionBody, double* x, double* y, double* yErr, int pNr, int varNr, double* res, double* err, double& redChiSqr) {
  char* gnuDataFileName = new char[2000];
  char* rmgnuDataFileName = new char[2000];
  snprintf(gnuDataFileName,2000,"%s/data/GnuplotFit%s.dat",DataBaseDirectory,GeneralUniqueFileNameExtension);
  snprintf(rmgnuDataFileName,2000,"rm %s/data/GnuplotFit%s.dat",DataBaseDirectory,GeneralUniqueFileNameExtension);

  FILE* file = fopen(gnuDataFileName,"w");

  for (int I=0; I<pNr; I++) {
    fprintf(file,"%1.15f %1.15f %1.15f\n", x[I], y[I], yErr[I]);
  }
  fclose(file);
  
  char* function = new char[1000];
  char* command = new char[1000];
  char* c2 = new char[1000];
  bool b = true;
  
  snprintf(function,1000, "f(x) = %s", functionBody);
  snprintf(command,1000, "fit f(x) '%s' using 1:2:3 via", gnuDataFileName);
  for (int I=0; I<varNr; I++) {
    if (I==0) {
      snprintf(c2,1000, "%s A%d", command, I+1);
    } else {
      snprintf(c2,1000, "%s, A%d", command, I+1);
    }
    snprintf(command,1000, "%s", c2);
  }
  b = b && performGnuplotFit(function, command,  varNr, varNr, res, err, redChiSqr);
  
  delete[] function;
  delete[] command;
  delete[] c2;
  if (!DebugMode) system(rmgnuDataFileName);    
  
  delete[] gnuDataFileName;
  delete[] rmgnuDataFileName;
  return b;
}


bool performGnuplotFitWithErrorEstimateFromResampling(char* functionBody, double* x, double* y, double* yErr, int pNr, int varNr, double* res, double* err, double& redChiSqr, int iterations) {
  double* resampledY = new double[pNr];
  double* dummyFitConst = new double[varNr];
  double* dummyFitConstErr = new double[varNr];
  double dummyChiSqr = 0;
  double* sigmaHelper1 = new double[varNr];
  double* sigmaHelper2 = new double[varNr];
  bool ok = true;

  ok = ok & performGnuplotFit(functionBody, x, y, yErr, pNr, varNr, res, err, redChiSqr);

  for (int I=0; I<varNr; I++) {
    sigmaHelper1[I] = 0;
    sigmaHelper2[I] = 0;    
  }
  
  for (int iter=0;iter<iterations; iter++) {  
    for (int I=0; I<pNr; I++) {
      double g1 = 0;
      double g2 = 0;
      AdvancedGaussZufall(AdvancedSeed, g1, g2);
      resampledY[I] = y[I] + g1*yErr[I];
    }
    for (int I=0; I<varNr; I++) {
      dummyFitConst[I] = res[I];
      dummyFitConstErr[I] = err[I];
    }
    ok = ok & performGnuplotFit(functionBody, x, resampledY, yErr, pNr, varNr, dummyFitConst, dummyFitConstErr, dummyChiSqr);
    for (int I=0; I<varNr; I++) {
      sigmaHelper1[I] += dummyFitConst[I];
      sigmaHelper2[I] += sqr(dummyFitConst[I]);    
    }
  }

  for (int I=0; I<varNr; I++) {
    err[I] = sqrt(sigmaHelper2[I]/iterations - sqr(sigmaHelper1[I]/iterations));
  }
  
  delete[] resampledY;
  delete[] dummyFitConst;
  delete[] dummyFitConstErr;
  delete[] sigmaHelper1;
  delete[] sigmaHelper2;  
  return ok;
}


void executeNeubergerDiracMultiplicationInFourierSpace(int L0, int L1, int L2, int L3, Complex* input, Complex* output, Complex* sinP, Complex* auxData) {
  int I0, I1, I2, I3, I, dob;
  int count = 0;
  int countAux = 0;
  double p0,p1,p2,p3;
  Complex* work = new Complex[4];

  for (I0=2*(L0-1); I0>=0; I0-=2) {
    p0 = sinP[0*256+I0].x;
    for (I1=2*(L1-1); I1>=0; I1-=2) {
      p1 = sinP[1*256+I1].x;
      for (I2=2*(L2-1); I2>=0; I2-=2) {
        p2 = sinP[2*256+I2].x;
        for (I3=2*(L3-1); I3>=0; I3-=2) {
          p3 = sinP[3*256+I3].x;
	  
	  for (dob=0; dob<2; dob++) {
            work[0].x = auxData[countAux].x * input[count+0].x + auxData[countAux].y*(p3*input[count+2].x - p0*input[count+2].y 
	                                                       + p1*input[count+3].x + p2*input[count+3].y);
            work[0].y = auxData[countAux].x * input[count+0].y + auxData[countAux].y*(p3*input[count+2].y + p0*input[count+2].x 
	                                                       + p1*input[count+3].y - p2*input[count+3].x);
							       
            work[1].x = auxData[countAux].x * input[count+1].x + auxData[countAux].y*(p1*input[count+2].x - p2*input[count+2].y 
	                                                       - p3*input[count+3].x - p0*input[count+3].y);
            work[1].y = auxData[countAux].x * input[count+1].y + auxData[countAux].y*(p1*input[count+2].y + p2*input[count+2].x 
	                                                       - p3*input[count+3].y + p0*input[count+3].x);

            work[2].x = auxData[countAux].x * input[count+2].x + auxData[countAux].y*(-p3*input[count+0].x - p0*input[count+0].y 
	                                                       - p1*input[count+1].x - p2*input[count+1].y);
            work[2].y = auxData[countAux].x * input[count+2].y + auxData[countAux].y*(-p3*input[count+0].y + p0*input[count+0].x 
	                                                       - p1*input[count+1].y + p2*input[count+1].x);

            work[3].x = auxData[countAux].x * input[count+3].x + auxData[countAux].y*(-p1*input[count+0].x + p2*input[count+0].y 
	                                                       + p3*input[count+1].x - p0*input[count+1].y);
            work[3].y = auxData[countAux].x * input[count+3].y + auxData[countAux].y*(-p1*input[count+0].y - p2*input[count+0].x 
	                                                       + p3*input[count+1].y + p0*input[count+1].x);
				  
  	    for (I=0; I<4; I++) {
  	      output[count+I].x = work[I].x;
  	      output[count+I].y = work[I].y;
	    }
	  
	    count += 4;
	  }
          countAux++;
        }
      }
    }
  }
  delete[] work;
}


int powINT(int base, int ex) {
  int res = 1;
  int I;
  for (I=1; I<=ex; I++) {
    res *= base;
  }
  return res;
}


int bitInverter(int bits, int length) {
  int I;
  int res = 0;
  for (I=0; I<length; I++) {
    res *=2;
    res += (bits % 2);
    bits /= 2;
  }
  return res;
}


void findPrimeFactors(int number, int* &primeFactors, int& count) {
  int x = number;
  int factor = 2;
  count = 0;
  if (number == 0) count = 0;
  if (number < 0) {
    count++;
    x = -number;
  }
  
  while (x>1) {
    while ((x % factor) == 0) {
      x = x / factor;
      count++;    
    }
    factor++;  
  }

  primeFactors = new int[count];
  count = 0;
  x = number;
  factor = 2;
  if (number < 0) {
    primeFactors[0] = -1;
    count++;
    x = -number;
  }
  
  while (x>1) {
    while ((x % factor) == 0) {
      x = x / factor;
      primeFactors[count] = factor;
      count++;    
    }
    factor++;  
  }
}


int primeNumberBasedInverter(int number, int* primeNumbers, int count) {
  if (count == 0) return 0;
  int max = primeNumbers[0];
  for (int I=1; I<count; I++) max *= primeNumbers[I];
  if (number<0) number = -number;
  number = number % max;
  if (number == 0) return 0;  
  int res = 0;
  int fac = max;
  for (int I=0; I<count; I++) {
    int x = number % primeNumbers[I];
    number /= primeNumbers[I];
    fac /= primeNumbers[I];
    res += x*fac;    
  }
  return res;
}


double calcFiniteVolumeEffectiveAction(int Nf, int L0, int L1, int L2, int L3, double kappaTilde, double lambdaTilde, double m, double s) {
  if (m==s) return 1E100;
  
  double S = 0;
  S += -8.0*kappaTilde * (m*m - s*s);
  S += m*m + s*s;
  S += lambdaTilde * (m*m*m*m + s*s*s*s + 6*m*m*s*s - 2*(m*m+s*s));
  S += -4.0* log(abs(m*m-s*s));
  S +=  (-8.0/(L0*L1*L2*L3)) * log(abs(m/(m*m-s*s)));
  S += (-56.0/(L0*L1*L2*L3)) * log(abs(1/(m*m-s*s)));

  S *= Nf*L0*L1*L2*L3;
   
  return S;
}


double findFiniteVolumeEffectiveActionGroundState(int Nf, int L0, int L1, int L2, int L3, double kappa, double lambda, double& minM, double& minS) {
  int I,I2;
  int iter = 1000;
  double maxMS = 5;
  double Smin = 1E100;
  minM = NaN;
  minS = NaN;
  
  for (I=0; I<=iter; I++) {
    for (I2=0; I2<=iter; I2++) {
      double m = (maxMS*I) / iter;
      double s = (maxMS*I2) / iter;
      
      double S = calcFiniteVolumeEffectiveAction(Nf, L0, L1, L2, L3, kappa, Nf * lambda, m, s);
      if (S<Smin) {
        Smin = S;
	minM = sqrt(1.0 * Nf) * m;
	minS = sqrt(1.0 * Nf) * s;
      }
    }
  }
  return Smin;
}


void makeFiniteVolumeKappaFunction(int Nf, int L0, int L1, int L2, int L3, double lambda, double kappaMin, double kappaMax, int steps, double* &k, double* &m, double* &s) {
  double kappa;
  double minM;
  double minS;
  
  k = new double[steps];
  m = new double[steps];
  s = new double[steps];

  int I;
  for (I=0; I<steps; I++) {
    kappa = ((kappaMax-kappaMin)/steps)*I + kappaMin;
    findFiniteVolumeEffectiveActionGroundState(Nf, L0, L1, L2, L3, kappa, lambda, minM, minS);
    k[I] = kappa;
  }
}


void copyFile(char* sourceFileName, char* destFileName) {
  if (LogLevel>3) printf("Copying file %s to %s\n",sourceFileName,destFileName);
  char* command = new char[2000];
  snprintf(command,1200,"cp %s %s",sourceFileName,destFileName);
  system(command);

  delete[] command;
}


double EffectiveMassGradientMinimizationSolverHelper(double* data) {
  double arg0 = data[0] * (EffectiveMassGradientMinimizationSolverHelper_relTime0);
  double arg1 = data[0] * (EffectiveMassGradientMinimizationSolverHelper_relTime1);
  double val = (exp(arg0) + exp(-arg0)) / (exp(arg1) + exp(-arg1));
    
  double chiTotal = sqr(val - EffectiveMassGradientMinimizationSolverHelper_QuotValue);

  return chiTotal;
}


double EffectiveMassSolver(double t0, double t1, double v0, double v1, double timeExtent) {  
  EffectiveMassGradientMinimizationSolverHelper_relTime0 = t0 - 0.5*timeExtent;
  EffectiveMassGradientMinimizationSolverHelper_relTime1 = t1 - 0.5*timeExtent;
  EffectiveMassGradientMinimizationSolverHelper_QuotValue = v0 / v1;
  
  double* fitData = new double[1];
  fitData[0] = 1.0;
  double effMass;
  double lowBound = 0;
  if (GradientMinimization(&EffectiveMassGradientMinimizationSolverHelper, 1, 1E-3, 1E-6, 1E-7, fitData, &lowBound, NULL, NULL, -1, 10000)) {
    effMass = fitData[0];
  } else {
    effMass = NaN;
  }
  delete[] fitData;
  
  return effMass;
}
      

char* getHostName() {
  int nr = (int)(100000*zufall());
  char* fileName = new char[600];
  char* hName = new char[600];
  snprintf(fileName,600,"hName%d%d.txt",nr,ownNodeID);
  char* command = new char[600];
  snprintf(command,600,"hostname >> %s",fileName);
  system(command);
  FILE* file = fopen(fileName,"r");
  fscanf(file,"%s",hName);
  fclose(file);
  snprintf(command,600,"rm %s",fileName); 
  system(command);
  delete[] command;
  delete[] fileName;
  return hName;
}


char* getTuningDBFileName(int L0, int L1, int L2, int L3, bool xFFT) {
  char* tuningDBFileName = new char[1000];
  if (xFFT) {
    snprintf(tuningDBFileName, 1000, "%s/data/tuningDB/tuningDB_xFFT_%dx%dx%dx%d.dat", DataBaseDirectory,L0,L1,L2,L3);
  } else {
    snprintf(tuningDBFileName, 1000, "%s/data/tuningDB/tuningDB_FFTW_%dx%dx%dx%d.dat", DataBaseDirectory,L0,L1,L2,L3);
  }
  return tuningDBFileName;
}


char* getTuningPerformanceProfileFileName(int L0, int L1, int L2, int L3, bool xFFT) {
  char* tuningPPFileName = new char[1000];
  if (xFFT) {
    snprintf(tuningPPFileName, 1000, "%s/data/tuningDB/PerformanceProfile_xFFT_%dx%dx%dx%d.dat", DataBaseDirectory,L0,L1,L2,L3);
  } else {
    snprintf(tuningPPFileName, 1000, "%s/data/tuningDB/PerformanceProfile_FFTW_%dx%dx%dx%d.dat", DataBaseDirectory,L0,L1,L2,L3);
  }
  return tuningPPFileName;
}


bool readOptimalFermionVectorEmbeddingAndFFTPlanFromTuningDB(int L0, int L1, int L2, int L3, int threadCountPerNode, int ParaOpMode, int xFFT, int useP, int useQ, int useR, int QHM, char* &fftPlanDescriptor) {
  char* hostName = getHostName();

  if (LogLevel>0) printf("Optimizing Embedding and FFT-Plan for lattice %dx%dx%dx%d, threadCountPerNode=%d, POM=%d, xFFT=%d, P=%d, Q=%d, R=%d, QHMode=%d and host=%s...\n", L0,L1,L2,L3, threadCountPerNode, ParaOpMode, xFFT, useP, useQ, useR, QHM, hostName);
  initializeNUMASupport();
  
  char* tuningDBFileName = getTuningDBFileName(L0,L1,L2,L3, (bool)xFFT);
  if (LogLevel>0) printf("Tuning-DB is %s\n", tuningDBFileName);
  TuningDataBase* tuningDB = new TuningDataBase(tuningDBFileName);
  delete[] tuningDBFileName;

  int xtrSize1 = 0;
  int xtrSize2 = 0;
  int xtrSize3 = 0;
  
  bool b = tuningDB->queryFastestEmbedding(hostName, numaCoresCount, L0, L1, L2, L3, ParaOpMode, threadCountPerNode , xFFT, useP, useQ, useR, QHM, xtrSize1, xtrSize2, xtrSize3);
  if (b) {
    xtraSize1 = xtrSize1;
    xtraSize2 = xtrSize2;
    xtraSize3 = xtrSize3;
     
    if (LogLevel>0) printf("Optimal Ebedding set to  %d x (%d+%d) x (%d+%d) x (%d+%d)\n",L0,L1,xtraSize1,L2,xtraSize2,L3,xtraSize3);  

    b = b & tuningDB->queryFastestFFTPlanDescriptor(hostName, numaCoresCount, L0, L1, L2, L3, ParaOpMode, threadCountPerNode , xFFT, useP, useQ, useR, QHM, xtraSize1,xtraSize2,xtraSize3, fftPlanDescriptor);
    if (b) {
      if (LogLevel>0) printf("Best FFT-Plan is: %s\n",fftPlanDescriptor);  
    } else {
      if (LogLevel>0) printf("No FFT-Plan found!\n");  
    }  
  } else {
    if (LogLevel>0) printf("No entry found in Benchmark-Table. Keeping original embedding (%d,%d,%d) and FFT-plan!!!\n",xtraSize1,xtraSize2,xtraSize3);
  }
  delete tuningDB;
  delete[] hostName;
  
  return b;
}


double getAngle(double x, double y) {
  double d = sqrt(x*x + y*y);
  if (d<1E-10) return 0;
  double w = asin(x/d);
  if (y<0) w = pi - w;
  return w;
}


void getFileNameList(char* searchDir, char* searchString, char** &fileNames, int& fileCount) {
  fileCount = 0;
  fileNames = NULL;

  for (int loop=0; loop<2; loop++) {
    DIR *dp = NULL;
    struct dirent *ep;
  
    if (searchDir == NULL) {
      dp = opendir("."); 
    } else {
      dp = opendir(searchDir);
    }

    if (dp == NULL) return;
  
    while ((ep = readdir (dp)) != NULL) {
      bool match = true;

      if (searchString != NULL) {
        if (strlen(ep->d_name) < strlen(searchString)) {
	  match = false;
	} else {
          for (int I=0; I<((int)strlen(searchString)); I++) {
            if (searchString[I] != ep->d_name[I]) { 
  	      match = false;
	      break;
  	    }
	  }
        }
      }

      if (match) {
        if (loop==1) {
	  fileNames[fileCount] = new char[strlen(ep->d_name)+2+strlen(searchDir)];
	  snprintf(fileNames[fileCount], strlen(ep->d_name)+2+strlen(searchDir), "%s/%s", searchDir, ep->d_name);
	}
        fileCount++;
      }
    }
    
    if (loop==0) {
      fileNames = new char*[fileCount];
      fileCount = 0;
    }
    
    closedir(dp);
  }
}


void deleteFileNameList(char** &fileNames, int& fileCount) {
  for (int I=0; I<fileCount; I++) {
    delete[] fileNames[I];    
  }
  delete[] fileNames;
  fileNames = NULL;
  fileCount = 0;
}


int getLargestL(StateDescriptorReader* SDReader) {
  int L0 = SDReader->getL0();
  int L1 = SDReader->getL1();
  int L2 = SDReader->getL2();
  int L3 = SDReader->getL3();

  int LargestL = L0;
  if (L1>LargestL) LargestL = L1;
  if (L2>LargestL) LargestL = L2;
  if (L3>LargestL) LargestL = L3;
  
  return LargestL;
}


char* cloneString(char* s) {
  char* res = new char[1+strlen(s)];
  snprintf(res,1+strlen(s),"%s",s);
  return res;
}


ComplexMatrix getProjectorMatrix(int sign) {
  if ((sign!=1) && (sign!=-1)) {
    printf("ERROR in getProjectorMatrix: invalid parameter!!!\n");
    exit(0);  
  }
  
  ComplexMatrix mat(4);
  mat.setZero();
  mat.matrix[0][0].x = 0.5*(1+sign);
  mat.matrix[1][1].x = 0.5*(1+sign);
  mat.matrix[2][2].x = 0.5*(1-sign);
  mat.matrix[3][3].x = 0.5*(1-sign);
  
  return mat;
}


ComplexMatrix getThetaMatrix(int index) {
  if ((index<0) || (index>3)) {
    printf("ERROR in getThetaMatrix: invalid parameter!!!\n");
    exit(0);    
  }
  ComplexMatrix mat(2);
  mat.setZero();
  
  if (index==0) {
    mat.matrix[0][0].x = 1;
    mat.matrix[1][1].x = 1;    
  }
  if (index==1) {
    mat.matrix[0][1].y = -1;
    mat.matrix[1][0].y = -1;      
  }
  if (index==2) {
    mat.matrix[0][1].x = -1;
    mat.matrix[1][0].x = 1;      
  }
  if (index==3) {
    mat.matrix[0][0].y = -1;
    mat.matrix[1][1].y = 1;      
  }
  
  return mat;
}


/**
*  WARNING: Does not work if eigenvalues are degenerate.
**/
bool potentiateHermiteanComplexMatrix(ComplexMatrix mat, ComplexMatrix& res, double p) {
  bool ok = mat.calcEigenvaluesAndEigenvectors();
  res.resize(mat.matrixSize);
  res.setZero();
  for (int I=0; I<mat.matrixSize; I++) {
    ComplexMatrix proj(*(mat.rightEigenVectors[I]));
    Complex ew = mat.eigenvalues[I];
    Complex fac = pow(ew, p);
    if ((isNaN(fac.x)) || (isNaN(ew.y))) ok = false;
    res = res + fac * proj;
  }
  for (int I=0; I<mat.matrixSize; I++) {
    for (int I2=0; I2<mat.matrixSize; I2++) {
      if (I!=I2) {
        Complex sc(0,0);
        for (int I3=0; I3<mat.matrixSize; I3++) {
          sc = sc + adj((*(mat.rightEigenVectors[I])).vectorElements[I3]) * ((*(mat.rightEigenVectors[I2])).vectorElements[I3]);
        }
        if (sqrt((sc.x*sc.x + sc.y*sc.y)/mat.matrixSize) > 1E-10) ok = false;
      }
    }
  }
  
  return ok;
}


bool checkMat1IsASquareRootOfMat2(ComplexMatrix mat1, ComplexMatrix mat2) {
  double norm2 = 0;
  for (int I=0; I<mat2.matrixSize*mat2.matrixSize; I++) {
    norm2 += mat2.matrixElements[I].x*mat2.matrixElements[I].x;
    norm2 += mat2.matrixElements[I].y*mat2.matrixElements[I].y;
  }
  norm2 = sqrt(norm2);
  
  ComplexMatrix diff = mat2 - (mat1*mat1);

  double normd = 0;
  for (int I=0; I<mat2.matrixSize*mat2.matrixSize; I++) {
    normd += diff.matrixElements[I].x*diff.matrixElements[I].x;
    normd += diff.matrixElements[I].y*diff.matrixElements[I].y;
  }
  normd = sqrt(normd);

  if (normd>1E-10 * norm2) return false;
  return true;
}


bool checkMat1IsAnInverseSquareRootOfMat2(ComplexMatrix mat1, ComplexMatrix mat2) {
  double norm2 = mat2.matrixSize;
  
  ComplexMatrix diff = mat2 * (mat1*mat1);
  for (int I=0; I<mat2.matrixSize; I++) {
    diff.matrix[I][I].x -= 1.0;
  }

  double normd = 0;
  for (int I=0; I<mat2.matrixSize*mat2.matrixSize; I++) {
    normd += diff.matrixElements[I].x*diff.matrixElements[I].x;
    normd += diff.matrixElements[I].y*diff.matrixElements[I].y;
  }
  normd = sqrt(normd);

  if (normd>1E-10 * norm2) return false;
  return true;
}


double integrate(double (*func)(double p0, double p1, double p2, double p3, double para), double Parameter, vector4D start, vector4D end, double accuracy) {
  bool containsSing = false;
  if ((start[0]==0) && (start[1]==0) && (start[2]==0) && (start[3]==0)) {
    containsSing = true;
  }
  vector4D smallLengths;
  int ind[4];
  double smallBoxVol = 1;
  for (int I=0; I<4; I++) {
    smallLengths[I] = (end[I] - start[I]) / 3.0;
    smallBoxVol *= smallLengths[I];
  }
  
  if (containsSing) {
    vector4D p;
    for (int I=0; I<4; I++) {
      p[I] = 0.5*(start[I]+end[I]);
    }
    double estimate[3];
    estimate[0] = 81*smallBoxVol * ((*func)(p[0], p[1], p[2], p[3], Parameter));
  
    for (int es=1; es<3; es++) {
      estimate[es] = 0;
      double vFac = sqr(sqr((3.0/(1+es))));
      for (ind[0]=0; ind[0]<1+es; ind[0]++) {
        for (ind[1]=0; ind[1]<1+es; ind[1]++) {
          for (ind[2]=0; ind[2]<1+es; ind[2]++) {
            for (ind[3]=0; ind[3]<1+es; ind[3]++) {

              for (int I=0; I<4; I++) {
                p[I] = start[I] + (ind[I]+0.5) * (3.0/(1+es))*smallLengths[I];
	      }
              estimate[es] += vFac * smallBoxVol * ((*func)(p[0], p[1], p[2], p[3], Parameter));
	    }
          }
	}
      }
    }
    if ((abs((estimate[2]-estimate[0])) < accuracy) && (abs((estimate[2]-estimate[1])) < accuracy)) {
      return estimate[2];
    }

  } else {
    double estimate1 = 0;
    vector4D p;
    vector4D fac;
    vector4D dummyLengths;
    
    for (int I=0; I<4; I++) {
      dummyLengths[I] = 1.5* smallLengths[I];
    }    
    for (ind[0]=0; ind[0]<3; ind[0]++) {
      for (ind[1]=0; ind[1]<3; ind[1]++) {
        for (ind[2]=0; ind[2]<3; ind[2]++) {
          for (ind[3]=0; ind[3]<3; ind[3]++) {
	  
	    double fff = 1;
	    for (int I=0; I<4; I++) {
              p[I]=start[I] + ind[I]*dummyLengths[I];
              fac[I] = 1;
              if (ind[I]==1) fac[I] = 4;	    
	      fff *= fac[I];
	    }
	    
	    estimate1 += 81.0*smallBoxVol*fff * ((*func)(p[0], p[1], p[2], p[3], Parameter));
	  }
	}
      }
    }
    estimate1 /= 6*6*6*6;
    

    double estimate2 = 0;
    for (int I=0; I<4; I++) {
      dummyLengths[I] = 0.5 * dummyLengths[I];
    }    
    double dummyBoxVol = 81.0*smallBoxVol/16.0;
    for (ind[0]=0; ind[0]<5; ind[0]++) {
      for (ind[1]=0; ind[1]<5; ind[1]++) {
        for (ind[2]=0; ind[2]<5; ind[2]++) {
          for (ind[3]=0; ind[3]<5; ind[3]++) {
	  
	    double fff = 1;
	    for (int I=0; I<4; I++) {
              p[I]=start[I] + ind[I]*dummyLengths[I];
              fac[I] = 1;
              if ((ind[I]%2) == 1) fac[I] = 4;	    
	      if (ind[I]==2) fac[I] = 2;	    
	      fff *= fac[I];
	    }
	    
	    estimate2 += dummyBoxVol*fff * ((*func)(p[0], p[1], p[2], p[3], Parameter));
	  }
	}
      }
    }
    estimate2 /= 6*6*6*6;
    if (abs((estimate2-estimate1) / estimate2) < accuracy) {    
      return estimate2;
    }
  }
    
  vector4D startNEW;
  vector4D endNEW;
  double res = 0;
  for (ind[0]=0; ind[0]<3; ind[0]++) {
    for (ind[1]=0; ind[1]<3; ind[1]++) {
      for (ind[2]=0; ind[2]<3; ind[2]++) {
        for (ind[3]=0; ind[3]<3; ind[3]++) {

          for (int I=0; I<4; I++) {
	    startNEW[I] = start[I] + ind[I]*smallLengths[I];
	    endNEW[I] = start[I] + (ind[I]+1)*smallLengths[I];
	  }
	  
	  res += integrate(func, Parameter, startNEW, endNEW, accuracy);
	}
      }
    }
  }
  return res;
}


double integrate(double (*func)(double p, double para), double Parameter, double start, double end, double accuracy) {
  bool containsSing = false;
  if (start==0) {
    containsSing = true;
  }
  double smallLength = (end - start) / 3.0;
  double smallBoxVol = smallLength;

  if (containsSing) {
    double p = 0.5*(start+end);
    double estimate[3];
    estimate[0] = 3*smallBoxVol * ((*func)(p, Parameter));
  
    for (int es=1; es<3; es++) {
      estimate[es] = 0;
      double vFac = (3.0/(1+es));
      for (int ind=0; ind<1+es; ind++) {
        p = start + (ind+0.5) * (3.0/(1+es))*smallLength;
        estimate[es] += vFac * smallBoxVol * ((*func)(p, Parameter));
      }
    }
    if ((abs((estimate[2]-estimate[0])) < accuracy) && (abs((estimate[2]-estimate[1])) < accuracy)) {
      return estimate[2];
    }

  } else {
    double estimate1 = 0;
    double p;
    double fac;
    double dummyLength;
    
    dummyLength = 1.5 * smallLength;
    for (int ind=0; ind<3; ind++) {
      double fff = 1;
      p=start + ind*dummyLength;
      fac = 1;
      if (ind==1) fac = 4;	    
      fff *= fac;
	    
      estimate1 += 3.0*smallBoxVol*fff * ((*func)(p, Parameter));
    }
    estimate1 /= 6;
    

    double estimate2 = 0;
    dummyLength = 0.5 * dummyLength;
    double dummyBoxVol = 3.0*smallBoxVol/2.0;
    for (int ind=0; ind<5; ind++) {
      double fff = 1;
      p = start + ind*dummyLength;
      fac = 1;
      if ((ind%2) == 1) fac = 4;	    
      if (ind==2) fac = 2;	    
      fff *= fac;
	    
      estimate2 += dummyBoxVol*fff * ((*func)(p, Parameter));
    }
    estimate2 /= 6;
    if (abs((estimate2-estimate1) / estimate2) < accuracy) {    
      return estimate2;
    }
  }
    
  double startNEW;
  double endNEW;
  double res = 0;
  for (int ind=0; ind<3; ind++) {
    startNEW = start + ind*smallLength;
    endNEW = start + (ind+1)*smallLength;
	  
    res += integrate(func, Parameter, startNEW, endNEW, accuracy);
  }
  return res;
}


double StandardErrorFunctionHelper(double p, double para) {
  return exp(-p*p);
}


/*
*  = 1- sqrt(1/2pi)*int_{-Nsigma}^{Nsigma} exp(-x*x)
*/
double StandardErrorFunction(double Nsigma, double accuracy) {
  return 2*sqrt(1.0/pi) * integrate(&StandardErrorFunctionHelper, 0, 0, Nsigma, accuracy);  
}


double inverseStandardErrorFunction(double p, double accuracy, double accuracy2) {
  if (p>=1) return NaN;
  if (p<=0) return 0;
  
  //Find upper bound
  double lowBoundx = 0;
  double lowBoundProb = 0.0;
  double upBoundx = 1.0;
  double upBoundProb = StandardErrorFunction(upBoundx, accuracy);
  while (upBoundProb<=p) {
    upBoundx += 1.0;
    upBoundProb = StandardErrorFunction(upBoundx, accuracy);
  }
  
  //Locate x
  double x = NaN;
  while (abs(lowBoundx - upBoundx)>accuracy2) {
    x = 0.5 * (lowBoundx + upBoundx);
    double val = StandardErrorFunction(x, accuracy); 
    if (abs(val-p)<accuracy) return x;
    if (val>p) {
      upBoundx = x;
      upBoundProb = val;
    } else {
      lowBoundx = x;
      lowBoundProb = val;
    }
  }
  return x;  
}


double LuescherZetaFunctionHelper1(double p, double para) {
  double f1 = exp(-1.5 * log(4*pi*p));
  double f2 = exp(p*para) - 1;
  return 8*pi*pi*pi*f1*f2;
}


double LuescherZetaFunctionHelper2(double p, double para) {
  double f1 = exp(-1.5 * log(4*pi*p));
  double f2 = exp(p*para);
  double f3 = exp(-pi*pi*LuescherZetaFunctionHelper2_NvecSqr / p);
  return 1E4 * 8*pi*pi*pi*f1*f2*f3;
}


double calcLuescherZetaFunction(double qSqr, double accuracy) {
  double res = -2*exp(1.5*log(pi));

  res += integrate(&LuescherZetaFunctionHelper1, qSqr, 0.0, 1.0, accuracy);

  int summandsMAX = 1000000;
  double* summands1 = new double[summandsMAX];
  double* summands2 = new double[summandsMAX];
  for (int I=0; I<summandsMAX; I++) {
    summands1[I] = NaN;
    summands2[I] = NaN;    
  }
  
  double sum1 = NaN;
  double lastSum1 = NaN;
  int box1 = 5;
  while ((isNaN(lastSum1)) || (abs(sum1-lastSum1) > abs(accuracy * sum1))) {
    lastSum1 = sum1;
    sum1 = 0;
    for (int n1=-box1; n1<=box1; n1++) {
      for (int n2=-box1; n2<=box1; n2++) {
        for (int n3=-box1; n3<=box1; n3++) {
	  int nSqr = n1*n1 + n2*n2 + n3*n3;
	  if (nSqr >= summandsMAX) {
	    printf("ERROR in calcLuescherZetaFunction: Sum1 did not converge\n");
	    exit(0);
	  }
	  if (isNaN(summands1[nSqr])) {
  	    summands1[nSqr] = exp(qSqr - nSqr) / (nSqr - qSqr);
	  }
	  sum1 += summands1[nSqr];
	}
      }
    }
    box1++;
  }
  res += sum1;

  double sum2 = NaN;
  double lastSum2 = NaN;
  int box2 = 5;
  while ((isNaN(lastSum2)) || (abs(sum2-lastSum2) > abs(accuracy * sum2))) {
    lastSum2 = sum2;
    sum2 = 0;
    for (int n1=-box2; n1<=box2; n1++) {
      for (int n2=-box2; n2<=box2; n2++) {
        for (int n3=-box2; n3<=box2; n3++) {
	  int nSqr = n1*n1 + n2*n2 + n3*n3;
	  if (nSqr >= summandsMAX) {
	    printf("ERROR in calcLuescherZetaFunction: Sum2 did not converge\n");
	    exit(0);
	  }
	  if (nSqr > 0) {
  	    if (isNaN(summands2[nSqr])) {
 	      LuescherZetaFunctionHelper2_NvecSqr = nSqr;
  	      summands2[nSqr] = 1E-4 * integrate(&LuescherZetaFunctionHelper2, qSqr, 0.0, 1.0, accuracy);
	    }
	    sum2 += summands2[nSqr];
	  }
	}
      }
    }
    box2++;
  }
  res += sum2;

  delete[] summands1;
  delete[] summands2;
  
  return res / sqrt(4*pi);
}



double LuescherDerivativeOfZetaFunction_dZdqSqrHelper2(double p, double para) {
  double f1 = 1.0 / sqrt(p);
  double f2 = exp(p*para);
  double f3 = exp(-pi*pi*LuescherDerivativeOfZetaFunction_dZdqSqrHelper2_NvecSqr / p);
  return 1E4 * 0.5*pi*f1*f2*f3;
}


double calcLuescherDerivativeOfZetaFunction_dZdqSqr(double qSqr, double accuracy) {
  double res = 0;

  int summandsMAX = 1000000;
  double* summands1 = new double[summandsMAX];
  double* summands2 = new double[summandsMAX];
  for (int I=0; I<summandsMAX; I++) {
    summands1[I] = NaN;
    summands2[I] = NaN;    
  }
  
  double sum1 = NaN;
  double lastSum1 = NaN;
  int box1 = 5;
  while ((isNaN(lastSum1)) || (abs(sum1-lastSum1) > abs(accuracy * sum1))) {
    lastSum1 = sum1;
    sum1 = 0;
    for (int n1=-box1; n1<=box1; n1++) {
      for (int n2=-box1; n2<=box1; n2++) {
        for (int n3=-box1; n3<=box1; n3++) {
	  int nSqr = n1*n1 + n2*n2 + n3*n3;
	  if (nSqr >= summandsMAX) {
	    printf("ERROR in calcLuescherZetaFunction: Sum1 did not converge\n");
	    exit(0);
	  }
	  if (isNaN(summands1[nSqr])) {
  	    summands1[nSqr] = (exp(qSqr - nSqr) * (nSqr - qSqr + 1)) / sqr(nSqr - qSqr);
	  }
	  sum1 += summands1[nSqr];
	}
      }
    }
    box1++;
  }
  res += sum1 / sqrt(4*pi);

  double sum2 = NaN;
  double lastSum2 = NaN;
  int box2 = 5;
  while ((isNaN(lastSum2)) || (abs(sum2-lastSum2) > abs(accuracy * sum2))) {
    lastSum2 = sum2;
    sum2 = 0;
    for (int n1=-box2; n1<=box2; n1++) {
      for (int n2=-box2; n2<=box2; n2++) {
        for (int n3=-box2; n3<=box2; n3++) {
	  int nSqr = n1*n1 + n2*n2 + n3*n3;
	  if (nSqr >= summandsMAX) {
	    printf("ERROR in calcLuescherZetaFunction: Sum2 did not converge\n");
	    exit(0);
	  }
	  if (isNaN(summands2[nSqr])) {
 	    LuescherDerivativeOfZetaFunction_dZdqSqrHelper2_NvecSqr = nSqr;
  	    summands2[nSqr] = 1E-4 * integrate(&LuescherDerivativeOfZetaFunction_dZdqSqrHelper2, qSqr, 0.0, 1.0, accuracy);
	  }
	  sum2 += summands2[nSqr];
	}
      }
    }
    box2++;
  }
  res += sum2;

  delete[] summands1;
  delete[] summands2;
  
  return res;
}


double calcLuescherPhiFunction(double q, double accuracy) {
  double zeta = calcLuescherZetaFunction(q*q,  accuracy);
  return -atan(q*exp(1.5*log(pi)) / zeta);
}


double calcLuescherPhiDerivative(double q, double accuracy) {
  double zeta = calcLuescherZetaFunction(q*q,  accuracy);
  double dzetadqq = calcLuescherDerivativeOfZetaFunction_dZdqSqr(q*q,  accuracy);

  return sqrt(pi*pi*pi) * (2*q*q*dzetadqq - zeta) / (zeta*zeta + pi*pi*pi*q*q);
}


long int faculty(long int x) {
  if (x<=1) return 1;
  long int res = 1;
  for (int I=1; I<=x; I++) res *= I;
  return res;
}


long int NoverP(long int n, long int p) {
  if (n<0) return 0;
  if (p<0) return 0;
  if (p>n) return 0;
  return faculty(n) / (faculty(p) * faculty(n-p));
}


Complex calcGoldstone1LoopInvPropagator(Complex p, double m0, double mg, double mh, double Z, double coeff) {  
  Complex Delta = Complex(mg*mg-mh*mh, 0);
  double p0val = 1 + 0.5*log(mg*mg/(mh*mh)) * (1+2*mh*mh/Delta.x);

  if (norm(p)==0) return Complex((1.0/Z) * (m0*m0 + coeff*p0val), 0);
  
  Complex q =  (Delta+p*p)*(Delta+p*p) + 4*mh*mh*p*p;
  Complex qsqrt = sqrt(q);
  Complex X = log(mh*mh/(mg*mg)) * Delta;
  Complex Y = log( ((qsqrt+p*p)*(qsqrt+p*p) - Delta*Delta)  /  ((qsqrt-p*p)*(qsqrt-p*p) - Delta*Delta) );

  Complex res = p*p + m0*m0 + 0.5*coeff*(X + qsqrt*Y)/(p*p);
  
  return (1.0/Z) * res;
}


Complex calcGoldstone1LoopInvPropagatorWithP0PartSubtracted(Complex p, double m0, double mg, double mh, double Z, double coeff) {  
  return Complex((1.0/Z) * (m0*m0),0) + calcGoldstone1LoopInvPropagator(p, m0, mg, mh, Z, coeff) - calcGoldstone1LoopInvPropagator(ComplexZero, m0, mg, mh, Z, coeff);
}


/*Complex calcGoldstone1LoopInvPropagatorHighPrecision(Complex p, double m0, double mg, double mh, double Z, double coeff) {  
  HighPrecisionComplex pHP(1000, p);
  HighPrecisionComplex m0HP(1000, m0);
  HighPrecisionComplex mgHP(1000, mg);
  HighPrecisionComplex mhHP(1000, mh);
  HighPrecisionComplex ZHP(1000, Z);
  HighPrecisionComplex coeffHP(1000, coeff);
  HighPrecisionComplex half(1000);
  half.setHalf();
  HighPrecisionComplex one(1000);
  one.setOne();
  HighPrecisionComplex two(1000);
  two.setTwo();
    
  HighPrecisionComplex DeltaHP = mgHP*mgHP-mhHP*mhHP;  
  HighPrecisionComplex p0valHP = one + half*log(mgHP*mgHP/(mhHP*mhHP)) * (one+two*mhHP*mhHP/DeltaHP);

  if (norm(p)==0) return ((m0HP*m0HP + coeffHP*p0valHP)/ZHP).getComplex();
  
  HighPrecisionComplex qHP =  (DeltaHP+pHP*pHP)*(DeltaHP+pHP*pHP) + two*two*mhHP*mhHP*pHP*pHP;
  HighPrecisionComplex qsqrtHP = sqrt(qHP);
  HighPrecisionComplex XHP = log(mhHP*mhHP/(mgHP*mgHP)) * DeltaHP;
  HighPrecisionComplex YHP = log( ((qsqrtHP+pHP*pHP)*(qsqrtHP+pHP*pHP) - DeltaHP*DeltaHP)  /  ((qsqrtHP-pHP*pHP)*(qsqrtHP-pHP*pHP) - DeltaHP*DeltaHP) );

  HighPrecisionComplex resHP = pHP*pHP + m0HP*m0HP + half*coeffHP*(XHP + qsqrtHP*YHP)/(pHP*pHP);
  resHP = resHP / ZHP;
  
  return resHP.getComplex();
}


Complex calcGoldstone1LoopInvPropagatorWithP0PartSubtractedHighPrecision(Complex p, double m0, double mg, double mh, double Z, double coeff) {  
  return Complex((1.0/Z) * (m0*m0),0) + calcGoldstone1LoopInvPropagatorHighPrecision(p, m0, mg, mh, Z, coeff) - calcGoldstone1LoopInvPropagatorHighPrecision(ComplexZero, m0, mg, mh, Z, coeff);
}*/


Complex calcBosonic1LoopContribution(Complex p, double m0) {  
  if (norm(p)==0) return Complex(0,0);
  Complex arg = sqrt(p*p / (p*p + 4*m0*m0));
  Complex argInv = ComplexUnity / arg;
    
  return (argInv * arctanh(arg)) - Complex(1,0);
}


Complex calcBosonic1LoopContributionOnSecondSheet(Complex p, double m0) {
  if (norm(p)==0) return Complex(0,0);
  double fac = sqrt(norm(p*p / (p*p + 4*m0*m0)));
  Complex res = fac*calcBosonic1LoopContribution(p, m0);  
  if ((p.y>=2*abs(m0)) && (p.x>=0)) res.y -= pi;
  if ((p.y<=-2*abs(m0)) && (p.x<=0)) res.y -= pi;
  
  return (1.0/fac) * res;
}


Complex calcBosonic1LoopInvPropagatorFromRenPT(Complex p, Complex HiggsPole, double mG, double lamRen, double vren, int n) {
  Complex CH0 = calcBosonic1LoopContributionOnSecondSheet(HiggsPole, HiggsPole.y);
  Complex CG0 = calcBosonic1LoopContributionOnSecondSheet(HiggsPole, mG);
  
  double piInvSqr = 1/(pi*pi);
  
  return p*p - HiggsPole*HiggsPole + 36*piInvSqr*sqr(lamRen*vren)*(calcBosonic1LoopContribution(p, HiggsPole.y)  - CH0) 
  + 4*piInvSqr*(n-1)*sqr(lamRen*vren)*(calcBosonic1LoopContribution(p, mG)  - CG0);
}


Complex calcBosonic1LoopInvPropagatorOnSecondSheetFromRenPT(Complex p, Complex HiggsPole, double mG, double lamRen, double vren, int n) {
  Complex CH0 = calcBosonic1LoopContributionOnSecondSheet(HiggsPole, HiggsPole.y);
  Complex CG0 = calcBosonic1LoopContributionOnSecondSheet(HiggsPole, mG);
  double piInvSqr = 1/(pi*pi);
  
  return p*p - HiggsPole*HiggsPole + 36*piInvSqr*sqr(lamRen*vren)*(calcBosonic1LoopContributionOnSecondSheet(p, HiggsPole.y)  - CH0) 
  + 4*piInvSqr*(n-1)*sqr(lamRen*vren)*(calcBosonic1LoopContributionOnSecondSheet(p, mG)  - CG0);
}


Complex findPoleOfBosonic1LoopPropagatorFromRenPT(double mH, double mG, double lamRen, double vren, int n) {
  Complex zero(0,0);
  double GammaOld = NaN;
  double Gamma = 0;
  Complex pole(Gamma/2,mH);
  
  if (n>1) {
    while (true) {
      pole.x = Gamma/2;
      double im = calcBosonic1LoopInvPropagatorFromRenPT(zero, pole, mG, lamRen, vren, n).y;    
            
      if ((!isNaN(GammaOld)) && (abs((Gamma-GammaOld)/GammaOld)<1E-10)) break;
      GammaOld = Gamma;   
      Gamma += im / mH; 
    }
  }
  return pole;
}


double calcSpectralFunctionOfBosonic1LoopPropagatorFromRenPT(double E, Complex HiggsPole, double mG, double lamRen, double vren, int n) {
  if ((n==1) && (E<2*HiggsPole.y)) return 0;
  Complex inv1 = calcBosonic1LoopInvPropagatorFromRenPT(Complex(-0.1E-4, E), HiggsPole, mG, lamRen, vren, n);
  Complex inv2 = calcBosonic1LoopInvPropagatorFromRenPT(Complex(+0.1E-4, E), HiggsPole, mG, lamRen, vren, n);
  
  return ((ComplexUnity/inv1) - (ComplexUnity/inv2)).y;
}


Complex calcBosonic1LoopInvPropagator(Complex p, double m0, int N, double* coeff) {
  Complex res = p*p + m0*m0;
  for (int I=0; I<N; I++) {
    res = res + coeff[2*I+0]*calcBosonic1LoopContribution(p, coeff[2*I+1]);
  }
  return res;
}


Complex calcBosonic1LoopInvPropagatorOnSecondSheet(Complex p, double m0, int N, double* coeff) {
  Complex res = p*p + m0*m0;
  for (int I=0; I<N; I++) {
    res = res + coeff[2*I+0]*calcBosonic1LoopContributionOnSecondSheet(p, coeff[2*I+1]);
  }
  return res;
}


Complex calcBosonic1LoopInvPropagatorFit(Complex p, double m0, double Z, int N, double* coeff) {
  return (1.0/Z) * calcBosonic1LoopInvPropagator(p, m0, N, coeff);
}


Complex calcBosonic1LoopInvPropagatorFitOnSecondSheet(Complex p, double m0, double Z, int N, double* coeff) {
  return (1.0/Z) * calcBosonic1LoopInvPropagatorOnSecondSheet(p, m0, N, coeff);
}


Complex calcBosonic1LoopCorrelatorFit(int t, int L, double m0, double Z, int N, double* coeff) {
  Complex res(0,0);
  double fac = 2*pi/L;
  for (int Pind=0; Pind<L; Pind++) {
    Complex p(2*sin(0.5*fac*Pind),0);
    res = res + exp(t*Pind*fac*ComplexI) / calcBosonic1LoopInvPropagatorFit(p, m0, Z, N, coeff);
  }
  res = (1.0/L) * res;
  return res;
}



double findZeroOfBosonic1LoopInvPropagatorFit_Helper(double* x) {
  Complex p(x[0], x[1]);
  Complex res = calcBosonic1LoopInvPropagatorFit(p, findZeroOfBosonic1LoopInvPropagatorFit_Helper_m0, findZeroOfBosonic1LoopInvPropagatorFit_Helper_Z, findZeroOfBosonic1LoopInvPropagatorFit_Helper_N, findZeroOfBosonic1LoopInvPropagatorFit_Helper_coeff);
  return norm(res);
}
double findZeroOfBosonic1LoopInvPropagatorFitOnSecondSheet_Helper(double* x) {
  Complex p(x[0], x[1]);
  Complex res = calcBosonic1LoopInvPropagatorFitOnSecondSheet(p, findZeroOfBosonic1LoopInvPropagatorFit_Helper_m0, findZeroOfBosonic1LoopInvPropagatorFit_Helper_Z, findZeroOfBosonic1LoopInvPropagatorFit_Helper_N, findZeroOfBosonic1LoopInvPropagatorFit_Helper_coeff);
  return norm(res);
}


bool findZeroOfBosonic1LoopInvPropagatorFit(Complex& polePos, Complex& valAtPole, double m0, double Z, int N, double* coeff, double startSearchM) {
  polePos.x = NaN;
  polePos.y = NaN;
  valAtPole.x = NaN;
  valAtPole.y = NaN;
  
  findZeroOfBosonic1LoopInvPropagatorFit_Helper_m0 = m0;
  findZeroOfBosonic1LoopInvPropagatorFit_Helper_Z = Z;
  findZeroOfBosonic1LoopInvPropagatorFit_Helper_N = N;
  findZeroOfBosonic1LoopInvPropagatorFit_Helper_coeff = coeff;
  
  double bounds1[2];
  double bounds2[2];
  double pos[2];  
  bounds1[0] = 0;
  bounds1[1] = 0;
  bounds2[0] = 1;
  bounds2[1] = 1;
  pos[0] = 0.001;
  pos[1] = startSearchM;
  if ((isNaN(startSearchM)) || (startSearchM<0)) pos[1] = m0;
  
  bool succ = GradientMinimization(&findZeroOfBosonic1LoopInvPropagatorFit_Helper, 2, 3E-4, 1E-7, 1E-10, &(pos[0]), &(bounds1[0]), &(bounds2[0]), NULL, 3, 1000);

  if (succ) {
    polePos.x = pos[0];
    polePos.y = pos[1];
    valAtPole = calcBosonic1LoopInvPropagatorFit(polePos, m0, Z, N, coeff);
  }
  
  return succ;
}


bool findZeroOfBosonic1LoopInvPropagatorOnSecondSheetFit(Complex& polePos, Complex& valAtPole, double m0, double Z, int N, double* coeff, double startSearchM) {
  polePos.x = NaN;
  polePos.y = NaN;
  valAtPole.x = NaN;
  valAtPole.y = NaN;
  
  findZeroOfBosonic1LoopInvPropagatorFit_Helper_m0 = m0;
  findZeroOfBosonic1LoopInvPropagatorFit_Helper_Z = Z;
  findZeroOfBosonic1LoopInvPropagatorFit_Helper_N = N;
  findZeroOfBosonic1LoopInvPropagatorFit_Helper_coeff = coeff;
  
  double bounds1[2];
  double bounds2[2];
  double pos[2];  
  bounds1[0] = 0;
  bounds1[1] = 0;
  bounds2[0] = 1;
  bounds2[1] = 1;
  pos[0] = 0.001;
  pos[1] = startSearchM;
  if ((isNaN(startSearchM)) || (startSearchM<0)) pos[1] = m0;
  
  bool succ = GradientMinimization(&findZeroOfBosonic1LoopInvPropagatorFitOnSecondSheet_Helper, 2, 3E-4, 1E-7, 1E-10, &(pos[0]), &(bounds1[0]), &(bounds2[0]), NULL, 3, 1000);

  if (succ) {
    polePos.x = pos[0];
    polePos.y = pos[1];
    valAtPole = calcBosonic1LoopInvPropagatorFitOnSecondSheet(polePos, m0, Z, N, coeff);
  }
  
  return succ;
}


long int NumberOfContractionsForRealScalarField(int p) {
  if (p<=0) return 0;
  if ((p%2)==1) return 0;
  long int res = 1;
  for (; p>0; p-=2) res *= (p-1);
  return res;
}


int NumberOfContractionsForRealScalarFieldAtTwo(int p1, int p2, long int* fac) {
  if (p1<0) return 0;
  if (p2<0) return 0;
  if (((p1+p2)%2)==1) return 0;
  int Nc = (p1+p2) / 2;
  for (int I=0; I<Nc+1; I++) {
    fac[I] = 0;
    if ((((p1-I)%2)==0) && (((p2-I)%2)==0)) {
      long int f1 = NumberOfContractionsForRealScalarField(p1-I);
      long int f2 = NumberOfContractionsForRealScalarField(p2-I);
      if ((p1-I)==0) f1 = 1;
      if ((p2-I)==0) f2 = 1;
      fac[I] = f1 * f2;
      fac[I] *= faculty(I) * NoverP(p1,I) * NoverP(p2,I);
    }
  }
  return Nc+1;
}


void calcAverageAndStandardDeviation(int N, double* data, double& avg, double& sig) {
  avg = 0;
  sig = 0;
  for (int I=0; I<N; I++) {
    avg += data[I];
    sig += sqr(data[I]);
  }
  avg /= N;
  sig = sqrt(sig/N - sqr(avg));
}


/*
* data will be overwritten
*/
void calcAverageAndStandardDeviationWithDataSelection(int N, double* data, double probFac, double& avg, double& sig) {
  if (N<10) {
    calcAverageAndStandardDeviation(N, data, avg, sig);
    return;
  }  
  bool change = true;
  while ((change) && (N>=10)) {
    change = false;
    double cutSigma = sqrt(2.0) * inverseStandardErrorFunction(1.0 - 1.0/(N*probFac), 1E-10, 1E-5);
    calcAverageAndStandardDeviation(N, data, avg, sig);
    bool sortChange = true;
    while (sortChange) {
      sortChange = false;
      for (int I=0; I<N-1; I++) {
        if (abs(data[I]-avg) > abs(data[I+1]-avg)) {
	  double dummy = data[I];
	  data[I] = data[I+1];
	  data[I+1] = dummy;
	  sortChange = true;
	}
      }
    }
//FILE* file = fopen("calcAverageWithDataSelection.dat","w");    
//fprintf(file, "%f %f %f \n",avg, sig, sig*cutSigma);
//fclose(file);
    for (int I=0; I<N; I++) {
      if (abs(data[I]-avg)/sig > cutSigma) {
        change = true;
        N = I;
	break;
      }
    }
  }
}
