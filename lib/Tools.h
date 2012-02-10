#ifndef Tools_included
#define Tools_included

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <dirent.h>
#include <pthread.h>
#include <numa.h>

#include "Global.h"
#include "Complex.h"
#include "Quat.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "StateDescriptorReader.h"
#include "TuningDataBase.h"
#include "PerformanceProfiler.h"



//Variablen
extern double eps;
extern double NaN;
extern double epsilon4[4][4][4][4];
extern Quat sigmaPlus[4];            //Der Einfachheit halber sigma0 = sigma4
extern Quat sigmaMinus[4];           //Der Einfachheit halber sigma0 = sigma4
extern int AdvancedSeed;
extern ComplexMatrix ID4x4;
extern ComplexMatrix Gamma5;
extern Complex* ASMParameterData;
extern bool numaAvail;
extern int numaNodesCount;
extern int numaCoresCount;
extern bool numaSupportInitialized;
extern int numaNodeAllocationMode;
extern PerformanceProfiler* performanceProfiler;
extern char* GeneralUniqueFileNameExtension;
extern Complex* ASMParameterData; 

int getCurentThreadID(); 
void setCurrentThreadAffinityMask(long int mask); 
long int getAffinityMaskFromCoreID(int coreID); 
long int getAffinityMaskFromNthCoreOnNodeID(int nodeID, int nthCore); 
long int getAffinityMaskFromNodeID(int nodeID); 
void setCurrentThreadAffinityToCoreID(int coreID); 
void setCurrentThreadAffinityToNthCoreOnNodeID(int nodeID, int nthCore); 
void setCurrentThreadAffinityToNodeID(int nodeID); 
void initializeNUMASupport(); 
Complex* createSuperAlignedComplex(int size, int ALIGN); 
Complex* createSuperAlignedComplex(int size); 
void destroySuperAlignedComplex(Complex* &p); 
void destroySuperAlignedIntP(int** &p); 
void destroySuperAlignedLongInt(long int* &p); 
void destroySuperAlignedInt(int* &p); 
void initializePerformanceProfiler(char* fileName); 


inline long int getPerformanceProfilingStartCycle(int node) {
  if (performanceProfiler == NULL) return 0;
  return performanceProfiler->getTimerStartCPUcycle(node);
}

inline void addPerformanceProfilingItem(char* routineName, long int startCycle, int node) {
  if (performanceProfiler == NULL) return;
  performanceProfiler->addPerformanceItem(routineName, startCycle, node);
}


void writePerformanceProfilingDataToDisk(); 


void printBits(long int x); 

#ifdef useSSE
  #include "SSEroutines.h"
#endif


long int getCPUCycleCounter();  


void calcEPS(); 


double zeitwert(); 


void delay(double timeInSecs); 


double cpuTime(); 


void determineCPUCyclesPerSecond(); 



void randomize(); 


double zufall(); 

/**
* Long period (>2 10^18) random number generator of L'Ecuyer with Bays-Durham shuffle
* and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
* the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
* idum between successive deviates in a sequence. RNMX should approximate the largest floating
* value that is less than 1.
*
**/
inline double AdvancedZufall(int &idum) {
  const int IM1=2147483563, IM2=2147483399;
  const int IA1=40014, IA2=40692, IQ1=53668, IQ2=52774;
  const int IR1=12211, IR2=3791, NTAB=32, IMM1=IM1-1;
  const int NDIV=1+IMM1/NTAB;
  const double EPS=3.0e-16, RNMX=(1.0-EPS), AM=1.0/double(IM1);
  static int idum2=123456789,iy=0;
  static int iv[NTAB];
  int j,k;
  double temp;
  
  if (idum <= 0) {
    idum=(idum==0 ? 1 : -idum);
    idum2=idum;
    for (j=NTAB+7; j>=0; j--) {
      k=idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum += IM1;
      if (j < NTAB) iv[j] = idum;
    }
    iy = iv[0];
  }
  k = idum/IQ1;
  idum = IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum += IM1;
  k = idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy/NDIV;
  iy = iv[j]-idum2;
  iv[j] = idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}


/**
* Returns two normally distributed deviates with zero mean and unit variance,
* using AdvancedZufall(idum) as the source of uniform deviates.
**/
inline void AdvancedGaussZufall(int &idum, double& z1, double& z2) {
  double Rsqr;
  
  do {
    z1 = 2.0 * AdvancedZufall(idum) - 1.0;
    z2 = 2.0 * AdvancedZufall(idum) - 1.0;
    Rsqr = z1*z1 + z2*z2;
  } while ((Rsqr>=1.0) || (Rsqr==0.0));
  
  double fac = sqrt(-2.0 * log(Rsqr)/Rsqr);
  z1 *= fac;
  z2 *= fac;
}


/**
* Returns one sinus-distributed deviate with mean pi/2 in intervall [0,pi]
* using AdvancedZufall(idum) as the source of uniform deviates.
**/  
inline double AdvancedSinusZufall(int &idum) {
  double r = AdvancedZufall(idum);
  return acos(r);
}


/**
* Returns a SU(2)-Haar measure distributed SU(2) matrix as Quaternion
* using AdvancedZufall(idum) as the source of uniform deviates.
*
*  U(phi_1, theta, phi_2)  =
*
*        exp( i phi_3 sigma_2 ) exp( i theta sigma_1 ) exp( i phi_2 sigma_3 )
*
*
*   The normalized invariant Haar measure on the group
*   manifold is
*
*     d(mu)  =  (1/16 pi^2) sin( theta ) d(phi_1) d(theta) d(phi_2)
*
*   with coordinate ranges
*
*     0 <= phi_1 < 4 pi,  0 <= theta <=  pi,  0 <= phi_2 < 2 pi
*
**/
inline Quat AdvancedSU2QuatZufall(int &idum) {
  Quat q1(0,-AdvancedSinusZufall(idum),0,0);
  Quat q2(0,0,-4*pi*AdvancedZufall(idum),0);
  Quat q3(0,0,0,-2*pi*AdvancedZufall(idum));

  Quat res = exp(q2) * exp(q1) * exp(q3);

  return res;
}


/**
* Returns a SU(2)-Haar measure distributed SU(2) matrix
* using AdvancedZufall(idum) as the source of uniform deviates.
**/
inline ComplexMatrix AdvancedSU2Zufall(int &idum) {
  ComplexMatrix res(AdvancedSU2QuatZufall(idum));
  return res;
}


// // double round(double x);  //why needed? 


int roundToInt(double x); 

// inline double fabs(double x) {   //why needed?
//   if (x>=0) return x;
//   return -x;
// }

bool isInteger(double x); 


double makeNaN(); 


bool isNaN(double x); 

/*inline bool isInf(double x) {
  if (d == numeric_limits<double>::infinity() ) return 1;
  return 0;
}*/


inline double sqr(double x) {
  return x*x;
}

inline long double sqrl(long double x) {
  return x*x;
}


void print(double d); 


void print(vector4D v); 


void print(int i); 


int calcPermutationSignum(int *per, int N); 


void iniEpsilon4(); 


void iniTools(int randGenNum); 


void desiniTools(); 


double calcLogDetScaledAbsNorm(ComplexMatrix& mat, int removeSmallestModesCount, double scaleFac, bool execCalcEigenV); 


void printEigenvalues(char* OPname, ComplexMatrix& op, Complex fac); 


inline void calcGamma5VectorProduct(Complex* left, Complex* right, Complex& res) {
  res.x = left[0].x*right[0].x + left[0].y*right[0].y;
  res.y = left[0].x*right[0].y - left[0].y*right[0].x;

  res.x += left[1].x*right[1].x + left[1].y*right[1].y;
  res.y += left[1].x*right[1].y - left[1].y*right[1].x;

  res.x -= left[2].x*right[2].x + left[2].y*right[2].y;
  res.y -= left[2].x*right[2].y - left[2].y*right[2].x;

  res.x -= left[3].x*right[3].x + left[3].y*right[3].y;
  res.y -= left[3].x*right[3].y - left[3].y*right[3].x;
}


/*********************************************************************************
Routines to implement a minimum-search of a given 1-dimensional function.
The gradient is used.
*/
//Returns the gradient of the given function at x using stepsize EPS. Output: d
void derive(double (*func)(double* x), int N, double EPS, double* x, double* d, int gradientMask); 

//Normalizes a given vector to 1.
int normalizeVec(double* x, int N); 


//Verifies whether a point lies within a given boundary box given by b1, b2.
//Additionally sum<bSum is checked.
int insideBox(int N, double* x, double* b1, double* b2, double* bSUM); 


//Implements the minimum-search. Output: pos. Finds minimum of given function, using given stepsizes, stepsize for differentiation and bounds.
bool GradientMinimization(double (*func)(double* x), int N, double StartStepSize, double MinStepSize, double DiffEPS, double* pos, double* bounds1, double* bounds2, double* boundSUM, int gradientMask, int iterMax); 


inline void mulWithPhiMatB(double* phi, Complex* in, Complex* out, double fac2) {
  out[0].x = fac2*(phi[0]*in[0].x - phi[3]*in[0].y
           + phi[2]*in[4].x - phi[1]*in[4].y);
  out[0].y = fac2*(phi[0]*in[0].y + phi[3]*in[0].x
           + phi[2]*in[4].y + phi[1]*in[4].x);
  out[1].x = fac2*(phi[0]*in[1].x - phi[3]*in[1].y
           + phi[2]*in[5].x - phi[1]*in[5].y);
  out[1].y = fac2*(phi[0]*in[1].y + phi[3]*in[1].x
           + phi[2]*in[5].y + phi[1]*in[5].x);
  out[2].x = fac2*(phi[0]*in[2].x + phi[3]*in[2].y
           - phi[2]*in[6].x + phi[1]*in[6].y);
  out[2].y = fac2*(phi[0]*in[2].y - phi[3]*in[2].x
           - phi[2]*in[6].y - phi[1]*in[6].x);
  out[3].x = fac2*(phi[0]*in[3].x + phi[3]*in[3].y
           - phi[2]*in[7].x + phi[1]*in[7].y);
  out[3].y = fac2*(phi[0]*in[3].y - phi[3]*in[3].x
           - phi[2]*in[7].y - phi[1]*in[7].x);
	   
  out[4].x = fac2*(-phi[2]*in[0].x - phi[1]*in[0].y
           + phi[0]*in[4].x + phi[3]*in[4].y);
  out[4].y = fac2*(-phi[2]*in[0].y + phi[1]*in[0].x
           + phi[0]*in[4].y - phi[3]*in[4].x);
  out[5].x = fac2*(-phi[2]*in[1].x - phi[1]*in[1].y
           + phi[0]*in[5].x + phi[3]*in[5].y);
  out[5].y = fac2*(-phi[2]*in[1].y + phi[1]*in[1].x
           + phi[0]*in[5].y - phi[3]*in[5].x);
  out[6].x = fac2*(phi[2]*in[2].x + phi[1]*in[2].y
           + phi[0]*in[6].x - phi[3]*in[6].y);
  out[6].y = fac2*(phi[2]*in[2].y - phi[1]*in[2].x
           + phi[0]*in[6].y + phi[3]*in[6].x);
  out[7].x = fac2*(phi[2]*in[3].x + phi[1]*in[3].y
           + phi[0]*in[7].x - phi[3]*in[7].y);
  out[7].y = fac2*(phi[2]*in[3].y - phi[1]*in[3].x
           + phi[0]*in[7].y + phi[3]*in[7].x);
}


inline void mulWithPhiMatB(double* phi, Complex* in, Complex* out, double split, double fac) {
  out[0].x = fac*(phi[0]*in[0].x - phi[3]*in[0].y
           + phi[2]*in[4].x - phi[1]*in[4].y);
  out[0].y = fac*(phi[0]*in[0].y + phi[3]*in[0].x
           + phi[2]*in[4].y + phi[1]*in[4].x);
  out[1].x = fac*(phi[0]*in[1].x - phi[3]*in[1].y
           + phi[2]*in[5].x - phi[1]*in[5].y);
  out[1].y = fac*(phi[0]*in[1].y + phi[3]*in[1].x
           + phi[2]*in[5].y + phi[1]*in[5].x);
  out[2].x = fac*(phi[0]*in[2].x + phi[3]*in[2].y
           - split * (phi[2]*in[6].x - phi[1]*in[6].y));
  out[2].y = fac*(phi[0]*in[2].y - phi[3]*in[2].x
           - split * (phi[2]*in[6].y + phi[1]*in[6].x));
  out[3].x = fac*(phi[0]*in[3].x + phi[3]*in[3].y
           - split * (phi[2]*in[7].x - phi[1]*in[7].y));
  out[3].y = fac*(phi[0]*in[3].y - phi[3]*in[3].x
           - split * (phi[2]*in[7].y + phi[1]*in[7].x));
	   
  out[4].x = fac*split*(-phi[2]*in[0].x - phi[1]*in[0].y
           + phi[0]*in[4].x + phi[3]*in[4].y);
  out[4].y = fac*split*(-phi[2]*in[0].y + phi[1]*in[0].x
           + phi[0]*in[4].y - phi[3]*in[4].x);
  out[5].x = fac*split*(-phi[2]*in[1].x - phi[1]*in[1].y
           + phi[0]*in[5].x + phi[3]*in[5].y);
  out[5].y = fac*split*(-phi[2]*in[1].y + phi[1]*in[1].x
           + phi[0]*in[5].y - phi[3]*in[5].x);
  out[6].x = fac*(phi[2]*in[2].x + phi[1]*in[2].y
           + split * (phi[0]*in[6].x - phi[3]*in[6].y));
  out[6].y = fac*(phi[2]*in[2].y - phi[1]*in[2].x
           + split * (phi[0]*in[6].y + phi[3]*in[6].x));
  out[7].x = fac*(phi[2]*in[3].x + phi[1]*in[3].y
           + split * (phi[0]*in[7].x - phi[3]*in[7].y));
  out[7].y = fac*(phi[2]*in[3].y - phi[1]*in[3].x
           + split * (phi[0]*in[7].y + phi[3]*in[7].x));
}


void perform_yB(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double y, double* phi, Complex* input, Complex* output); 


void perform_yBsplit(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double y, double split, double* phi, Complex* input, Complex* output); 


void performf_YBD_2rhoD_2rhoD_AndScalarProducts(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double twoRho, double explicitMass, Complex* x, Complex* Dx, double* phi, Complex* output, Complex* Vrest, Complex* Vp, Complex& resVrests, Complex& resVps); 


void performf_YBsplitD_2rhoD_2rhoD_AndScalarProducts(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double split, double twoRho, double explicitMass, Complex* x, Complex* Dx, double* phi, Complex* output, Complex* Vrest, Complex* Vp, Complex& resVrests, Complex& resVps); 


void performf_YBD_2rhoD_2rhoD(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double twoRho,
double explicitMass, Complex* x, Complex* Dx, double* phi, Complex* output); 


void performf_YBsplitD_2rhoD_2rhoD(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double split, double twoRho, double explicitMass, Complex* x, Complex* Dx, double* phi, Complex* output); 


inline void mulWithPhiDaggeredMatB(double* phi, Complex* in, Complex* out, double fac2) {
  out[0].x = fac2*(phi[0]*in[0].x + phi[3]*in[0].y
           - phi[2]*in[4].x + phi[1]*in[4].y);
  out[0].y = fac2*(phi[0]*in[0].y - phi[3]*in[0].x
           - phi[2]*in[4].y - phi[1]*in[4].x);
  out[1].x = fac2*(phi[0]*in[1].x + phi[3]*in[1].y
           - phi[2]*in[5].x + phi[1]*in[5].y);
  out[1].y = fac2*(phi[0]*in[1].y - phi[3]*in[1].x
           - phi[2]*in[5].y - phi[1]*in[5].x);
  out[2].x = fac2*(phi[0]*in[2].x - phi[3]*in[2].y
           + phi[2]*in[6].x - phi[1]*in[6].y);
  out[2].y = fac2*(phi[0]*in[2].y + phi[3]*in[2].x
           + phi[2]*in[6].y + phi[1]*in[6].x);
  out[3].x = fac2*(phi[0]*in[3].x - phi[3]*in[3].y
           + phi[2]*in[7].x - phi[1]*in[7].y);
  out[3].y = fac2*(phi[0]*in[3].y + phi[3]*in[3].x
           + phi[2]*in[7].y + phi[1]*in[7].x);
	   
  out[4].x = fac2*(phi[2]*in[0].x + phi[1]*in[0].y
           + phi[0]*in[4].x - phi[3]*in[4].y);
  out[4].y = fac2*(phi[2]*in[0].y - phi[1]*in[0].x
           + phi[0]*in[4].y + phi[3]*in[4].x);
  out[5].x = fac2*(phi[2]*in[1].x + phi[1]*in[1].y
           + phi[0]*in[5].x - phi[3]*in[5].y);
  out[5].y = fac2*(phi[2]*in[1].y - phi[1]*in[1].x
           + phi[0]*in[5].y + phi[3]*in[5].x);
  out[6].x = fac2*(-phi[2]*in[2].x - phi[1]*in[2].y
           + phi[0]*in[6].x + phi[3]*in[6].y);
  out[6].y = fac2*(-phi[2]*in[2].y + phi[1]*in[2].x
           + phi[0]*in[6].y - phi[3]*in[6].x);
  out[7].x = fac2*(-phi[2]*in[3].x - phi[1]*in[3].y
           + phi[0]*in[7].x + phi[3]*in[7].y);
  out[7].y = fac2*(-phi[2]*in[3].y + phi[1]*in[3].x
           + phi[0]*in[7].y - phi[3]*in[7].x);
}


inline void mulWithPhiDaggeredMatB(double* phi, Complex* in, Complex* out, double split, double fac) {
  out[0].x = fac*(phi[0]*in[0].x + phi[3]*in[0].y
           - split * (phi[2]*in[4].x - phi[1]*in[4].y));
  out[0].y = fac*(phi[0]*in[0].y - phi[3]*in[0].x
           - split * (phi[2]*in[4].y + phi[1]*in[4].x));
  out[1].x = fac*(phi[0]*in[1].x + phi[3]*in[1].y
           - split * (phi[2]*in[5].x - phi[1]*in[5].y));
  out[1].y = fac*(phi[0]*in[1].y - phi[3]*in[1].x
           - split * (phi[2]*in[5].y + phi[1]*in[5].x));
  out[2].x = fac*(phi[0]*in[2].x - phi[3]*in[2].y
           + phi[2]*in[6].x - phi[1]*in[6].y);
  out[2].y = fac*(phi[0]*in[2].y + phi[3]*in[2].x
           + phi[2]*in[6].y + phi[1]*in[6].x);
  out[3].x = fac*(phi[0]*in[3].x - phi[3]*in[3].y
           + phi[2]*in[7].x - phi[1]*in[7].y);
  out[3].y = fac*(phi[0]*in[3].y + phi[3]*in[3].x
           + phi[2]*in[7].y + phi[1]*in[7].x);
	   
  out[4].x = fac*(phi[2]*in[0].x + phi[1]*in[0].y
           + split * (phi[0]*in[4].x - phi[3]*in[4].y));
  out[4].y = fac*(phi[2]*in[0].y - phi[1]*in[0].x
           + split * (phi[0]*in[4].y + phi[3]*in[4].x));
  out[5].x = fac*(phi[2]*in[1].x + phi[1]*in[1].y
           + split * (phi[0]*in[5].x - phi[3]*in[5].y));
  out[5].y = fac*(phi[2]*in[1].y - phi[1]*in[1].x
           + split * (phi[0]*in[5].y + phi[3]*in[5].x));
  out[6].x = fac*split*(-phi[2]*in[2].x - phi[1]*in[2].y
           + phi[0]*in[6].x + phi[3]*in[6].y);
  out[6].y = fac*split*(-phi[2]*in[2].y + phi[1]*in[2].x
           + phi[0]*in[6].y - phi[3]*in[6].x);
  out[7].x = fac*split*(-phi[2]*in[3].x - phi[1]*in[3].y
           + phi[0]*in[7].x + phi[3]*in[7].y);
  out[7].y = fac*split*(-phi[2]*in[3].y + phi[1]*in[3].x
           + phi[0]*in[7].y - phi[3]*in[7].x);
}


void perform_yBDagger(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double y, double* phi, Complex* input, Complex* output);  


void perform_yBsplitDagger(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double y, double split, double* phi, Complex* input, Complex* output); 


void performf_YB_2rho_AndCopyToOutput(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double twoRho, Complex* x, double* phi, Complex* interim, Complex* output); 


void performf_YBsplit_2rho_AndCopyToOutput(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double split, double twoRho, double explicitMass, Complex* x, double* phi, Complex* interim, Complex* output); 


ComplexMatrix* createPhiMatB(vector4D phi, bool daggered, double split ); 


ComplexMatrix* createPhiMatrix(vector4D phi, bool daggered); 


bool performGnuplotFit(char* function, char* fitCommand, int varNr, int optiNr, double* res, double* err, double& redChiSqr); 


bool performGnuplotFit(char* functionBody, double* x, double* y, double* yErr, int pNr, int varNr, double* res, double* err, double& redChiSqr); 


bool performGnuplotFitWithErrorEstimateFromResampling(char* functionBody, double* x, double* y, double* yErr, int pNr, int varNr, double* res, double* err, double& redChiSqr, int iterations); 


void executeNeubergerDiracMultiplicationInFourierSpace(int L0, int L1, int L2, int L3, Complex* input, Complex* output, Complex* sinP, Complex* auxData); 


int powINT(int base, int ex); 


int bitInverter(int bits, int length); 


void findPrimeFactors(int number, int* &primeFactors, int& count); 


int primeNumberBasedInverter(int number, int* primeNumbers, int count); 


double calcFiniteVolumeEffectiveAction(int Nf, int L0, int L1, int L2, int L3, double kappaTilde, double lambdaTilde, double m, double s); 


double findFiniteVolumeEffectiveActionGroundState(int Nf, int L0, int L1, int L2, int L3, double kappa, double lambda, double& minM, double& minS); 


void makeFiniteVolumeKappaFunction(int Nf, int L0, int L1, int L2, int L3, double lambda, double kappaMin, double kappaMax, int steps, double* &k, double* &m, double* &s); 


void copyFile(char* sourceFileName, char* destFileName); 


extern double EffectiveMassGradientMinimizationSolverHelper_relTime0;
extern double EffectiveMassGradientMinimizationSolverHelper_relTime1;
extern double EffectiveMassGradientMinimizationSolverHelper_QuotValue;

double EffectiveMassGradientMinimizationSolverHelper(double* data); 


double EffectiveMassSolver(double t0, double t1, double v0, double v1, double timeExtent); 
      

char* getHostName(); 


char* getTuningDBFileName(int L0, int L1, int L2, int L3, bool xFFT); 


char* getTuningPerformanceProfileFileName(int L0, int L1, int L2, int L3, bool xFFT); 


bool readOptimalFermionVectorEmbeddingAndFFTPlanFromTuningDB(int L0, int L1, int L2, int L3, int threadCountPerNode, int ParaOpMode, int xFFT, int useP, int useQ, int useR, int QHM, char* &fftPlanDescriptor); 


double getAngle(double x, double y); 


void getFileNameList(char* searchDir, char* searchString, char** &fileNames, int& fileCount); 


void deleteFileNameList(char** &fileNames, int& fileCount); 


int getLargestL(StateDescriptorReader* SDReader); 


char* cloneString(char* s); 


ComplexMatrix getProjectorMatrix(int sign); 


ComplexMatrix getThetaMatrix(int index); 


/**
*  WARNING: Does not work if eigenvalues are degenerate.
**/
bool potentiateHermiteanComplexMatrix(ComplexMatrix mat, ComplexMatrix& res, double p); 


bool checkMat1IsASquareRootOfMat2(ComplexMatrix mat1, ComplexMatrix mat2); 


bool checkMat1IsAnInverseSquareRootOfMat2(ComplexMatrix mat1, ComplexMatrix mat2); 


double integrate(double (*func)(double p0, double p1, double p2, double p3, double para), double Parameter, vector4D start, vector4D end, double accuracy); 


double integrate(double (*func)(double p, double para), double Parameter, double start, double end, double accuracy); 


double StandardErrorFunctionHelper(double p, double para); 


/*
*  = 1- sqrt(1/2pi)*int_{-Nsigma}^{Nsigma} exp(-x*x)
*/
double StandardErrorFunction(double Nsigma, double accuracy); 


double inverseStandardErrorFunction(double p, double accuracy, double accuracy2); 


double LuescherZetaFunctionHelper1(double p, double para); 


extern int LuescherZetaFunctionHelper2_NvecSqr;
extern int LuescherDerivativeOfZetaFunction_dZdqSqrHelper2_NvecSqr;

double LuescherZetaFunctionHelper2(double p, double para); 


double calcLuescherZetaFunction(double qSqr, double accuracy); 



double LuescherDerivativeOfZetaFunction_dZdqSqrHelper2(double p, double para); 


double calcLuescherDerivativeOfZetaFunction_dZdqSqr(double qSqr, double accuracy); 


double calcLuescherPhiFunction(double q, double accuracy); 


double calcLuescherPhiDerivative(double q, double accuracy); 


long int faculty(long int x); 


long int NoverP(long int n, long int p); 


Complex calcGoldstone1LoopInvPropagator(Complex p, double m0, double mg, double mh, double Z, double coeff); 


Complex calcGoldstone1LoopInvPropagatorWithP0PartSubtracted(Complex p, double m0, double mg, double mh, double Z, double coeff); 


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


Complex calcBosonic1LoopContribution(Complex p, double m0); 


Complex calcBosonic1LoopContributionOnSecondSheet(Complex p, double m0); 


Complex calcBosonic1LoopInvPropagatorFromRenPT(Complex p, Complex HiggsPole, double mG, double lamRen, double vren, int n); 


Complex calcBosonic1LoopInvPropagatorOnSecondSheetFromRenPT(Complex p, Complex HiggsPole, double mG, double lamRen, double vren, int n); 


Complex findPoleOfBosonic1LoopPropagatorFromRenPT(double mH, double mG, double lamRen, double vren, int n); 


double calcSpectralFunctionOfBosonic1LoopPropagatorFromRenPT(double E, Complex HiggsPole, double mG, double lamRen, double vren, int n); 


Complex calcBosonic1LoopInvPropagator(Complex p, double m0, int N, double* coeff); 


Complex calcBosonic1LoopInvPropagatorOnSecondSheet(Complex p, double m0, int N, double* coeff); 


Complex calcBosonic1LoopInvPropagatorFit(Complex p, double m0, double Z, int N, double* coeff); 


Complex calcBosonic1LoopInvPropagatorFitOnSecondSheet(Complex p, double m0, double Z, int N, double* coeff); 


Complex calcBosonic1LoopCorrelatorFit(int t, int L, double m0, double Z, int N, double* coeff); 


extern double findZeroOfBosonic1LoopInvPropagatorFit_Helper_m0;
extern double findZeroOfBosonic1LoopInvPropagatorFit_Helper_Z;
extern int findZeroOfBosonic1LoopInvPropagatorFit_Helper_N;
extern double* findZeroOfBosonic1LoopInvPropagatorFit_Helper_coeff;
double findZeroOfBosonic1LoopInvPropagatorFit_Helper(double* x);

double findZeroOfBosonic1LoopInvPropagatorFitOnSecondSheet_Helper(double* x); 


bool findZeroOfBosonic1LoopInvPropagatorFit(Complex& polePos, Complex& valAtPole, double m0, double Z, int N, double* coeff, double startSearchM); 


bool findZeroOfBosonic1LoopInvPropagatorOnSecondSheetFit(Complex& polePos, Complex& valAtPole, double m0, double Z, int N, double* coeff, double startSearchM); 


long int NumberOfContractionsForRealScalarField(int p); 


int NumberOfContractionsForRealScalarFieldAtTwo(int p1, int p2, long int* fac); 


void calcAverageAndStandardDeviation(int N, double* data, double& avg, double& sig); 


/*
* data will be overwritten
*/
void calcAverageAndStandardDeviationWithDataSelection(int N, double* data, double probFac, double& avg, double& sig); 

#endif
