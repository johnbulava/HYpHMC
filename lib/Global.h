#ifndef Global_included
#define Global_included

#define RunsAtZeuthen
#define SmallLattice_NO

#define Physical_VEV_GeV 246

#ifdef RunsAtZeuthen
 #define useBLAS
 #define useLAPACK
 #define useARPACK  
 #define CBLASIncludeFile <gsl/gsl_cblas.h>
 #define FFTWIncludeFile <fftw3.h>
 #define DataBaseDirectory "dataBase"
#else
 #define useBLAS
 #define useLAPACK
 #define useARPACK 
 #define CBLASIncludeFile <cblas.h>
 #define FFTWIncludeFile <fftw3.h>
 #define DataBaseDirectory "dataBase"
#endif 


#include "Complex.h"

#ifdef useLAPACK
  extern "C" {
   #define COMPLEX void
    void zgeev_(const char* jobvl, const char* jobvr, const long int* n, COMPLEX* a, const long int* lda, COMPLEX* w, const COMPLEX* vl, const long int* ldvl, const COMPLEX* vr, const long int* ldvr, COMPLEX* work, const long int* lwork, double* rwork, long int* info);
    void zgetrf_(const long int* m, const long int* n, COMPLEX* a, const long int* lda, long int* ipiv, long int* info);
    void zgetri_(const long int* n, COMPLEX* a, const long int* lda, const long int* ipiv, COMPLEX* work, const long int* lwork, long int* info);    
    void dsteqr_(const char* compz, const long int* n, double* d, double* e, const double* z, const long int*  ldz, double* work, long int* info);
   #undef COMPLEX
  }
#endif


#ifdef useARPACK
  extern "C" {  
   void znaupd_(int *ido, char *bmat, int *n, char *which,
                       int *nev, double *tol, Complex *resid,
                       int *ncv, Complex *V, int *ldv,
                       int *iparam, int *ipntr, Complex *workd,
                       Complex *workl, int *lworkl,
                       double *rwork, int *info);

  void zneupd_(int *rvec, char *HowMny, int *select,
                       Complex *d, Complex *Z, int *ldz,
                       Complex *sigma, Complex *workev,
                       char *bmat, int *n, char *which, int *nev,
                       double *tol, Complex *resid, int *ncv,
                       Complex *V, int *ldv, int *iparam,
                       int *ipntr, Complex *workd,
                       Complex *workl, int *lworkl,
                       double *rwork, int *info);
  }
#endif




//Constants
const double pi = 4.e0*atan(1.e0);
extern int LogLevel;
extern bool DebugMode; 
extern int nodeCount;
extern int ownNodeID;
extern int xtraSize3; //14;1; 4;  (16^4, 8^3x16, 16^3x32)
extern int xtraSize2; //4; 0; 2;
extern int xtraSize1; //2; 4; 2;
extern int xtraCACHESize3; //0;//xtraSize3;
extern int xtraCACHESize2; //2;//xtraSize2;
extern int xtraCACHESize1; //2;//xtraSize1;
const int L2CacheSizeInBytes = 1024*1024;
const int L1CacheSizeInBytes = 64*1024;
const int L1CacheSizePerWayInBytes = 32*1024;
const int L1Ways = 2;
const int L1CacheLineSizeInBytes = 64;
extern long int CPUCyclesPerSecond;

typedef double vector4D[4];


#endif
