#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>

#include "Complex.h"



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


void eigenVals() {
/*
c     %------------------------------------------------------%
c     | Storage Declarations:                                |
c     |                                                      |
c     | The maximum dimensions for all arrays are            |
c     | set here to accommodate a problem size of            |
c     | N .le. MAXN                                          |
c     |                                                      |
c     | NEV is the number of eigenvalues requested.          |
c     |     See specifications for ARPACK usage below.       |
c     |                                                      |
c     | NCV is the largest number of basis vectors that will |
c     |     be used in the Implicitly Restarted Arnoldi      |
c     |     Process.  Work per major iteration is            |
c     |     proportional to N*NCV*NCV.                       |
c     |                                                      |
c     | You must set:                                        |
c     |                                                      |
c     | MAXN:   Maximum dimension of the A allowed.          |
c     | MAXNEV: Maximum NEV allowed.                         |
c     | MAXNCV: Maximum NCV allowed.                         |
c     %------------------------------------------------------%
*/

  int maxn = 1024;
  int maxnev=12;
  int maxncv=100;
  int ldv=maxn;
  int I;
  
/*
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
*/

  int iparam[11];
  int ipntr[14];
  
  int* select = new int[maxncv];   
  Complex* ax = new Complex[maxn];
  Complex* d = new Complex[maxncv];
  Complex* v = new Complex[ldv*maxncv];
  Complex* workd = new Complex[3*maxn];
  Complex* workev = new Complex[2*maxncv];
  Complex* resid = new Complex[maxn];
  Complex* workl = new Complex[3*maxncv*maxncv+5*maxncv];
  
  double* rwork = new double[maxncv];
  double* rd = new double[maxncv*3];

/*
c     %---------------%
c     | Local Scalars |
c     %---------------%
*/

  char bmat[1];
  char which[2];
  int  ido, n, nx, nev, ncv, lworkl, info, ierr,
       j, ishfts, maxitr, mode1, nconv;
       
  Complex sigma;
  double tol;
  int rvec;
  

/*
c     %-------------------------------------------------%
c     | The following sets dimensions for this problem. |
c     %-------------------------------------------------%
*/

      nx    = 16;
      n     = nx*nx;

/*
c     %-----------------------------------------------%
c     |                                               | 
c     | Specifications for ARPACK usage are set       | 
c     | below:                                        |
c     |                                               |
c     |    1) NEV = 4  asks for 4 eigenvalues to be   |  
c     |       computed.                               | 
c     |                                               |
c     |    2) NCV = 20 sets the length of the Arnoldi |
c     |       factorization                           |
c     |                                               |
c     |    3) This is a standard problem              |
c     |         (indicated by bmat  = 'I')            |
c     |                                               |
c     |    4) Ask for the NEV eigenvalues of          |
c     |       largest magnitude                       |
c     |         (indicated by which = 'LM')           |
c     |       See documentation in ZNAUPD  for the     |
c     |       other options SM, LR, SR, LI, SI.       | 
c     |                                               |
c     | Note: NEV and NCV must satisfy the following  |
c     | conditions:                                   |
c     |              NEV <= MAXNEV                    |
c     |          NEV + 2 <= NCV <= MAXNCV             |
c     |                                               |
c     %-----------------------------------------------%
*/

      nev   = 4;
      ncv   = 60;
      bmat[0]  = 'I';
      which[0] = 'L';
      which[1] = 'R';
     

/*
c     %-----------------------------------------------------%
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling ZNAUPD                     |
c     |                                                     |
c     | TOL  determines the stopping criterion.             |
c     |                                                     |
c     |      Expect                                         |
c     |           fabs(lambdaC - lambdaT) < TOL*fabs(lambdaC) |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |           (machine precision) is used.              |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from ZNAUPD . (see usage below)                 |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to ZNAUPD .                                | 
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     | 
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID).  | 
c     |                                                     |
c     | The work array WORKL is used in ZNAUPD  as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
*/

      lworkl  = 3*ncv*ncv+5*ncv; 
      tol    = 1E-10;
      ido    = 0;
      info   = 0;

/*
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting IPARAM(1) = 1).             |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | ZNAUPD .                                           |
c     %---------------------------------------------------%
*/

  ishfts = 1;
  maxitr = 300;
  mode1 = 1;

  iparam[0] = ishfts;                
  iparam[2] = maxitr;          
  iparam[6] = mode1;

  
/*   
c        %---------------------------------------------%
c        | Repeatedly call the routine ZNAUPD  and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
*/


  for (I=0; I<n; I++) {
    resid[I].x = 0;
    resid[I].y = 0;  
  }


  while(true) {  
    znaupd_(&ido, &(bmat[0]), &n, &(which[0]), &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);


printf("run 1: info: %d ido: %d\n",info,ido);


    if (ido == -1) {
      int xInd = ipntr[0]-1;
      int yInd = ipntr[1]-1;
      
      for (I=0; I<n; I++) {
        workd[I+yInd].x=0;
        workd[I+yInd].y=0;
      }
    } 
    if (ido == 1) {
      int xInd = ipntr[0]-1;
      int yInd = ipntr[1]-1;
      int zInd = ipntr[2]-1;
     
     printf("indices: %d %d %d\n",xInd, yInd, zInd);
     
          
      for (I=0; I<n; I++) {
        workd[I+zInd].x=workd[I+xInd].x;
        workd[I+zInd].y=workd[I+xInd].y;

//        resid[I].print();
//        workd[I+xInd].print();
      
        workd[I+yInd]=Complex(1.2,-1.3)*workd[I+xInd];
      }
       printf("\n");
    
/*       for (I=0; I<n; I++) {
        workd[I+yInd].print();
       }
       printf("\n");
       
       for (I=0; I<n; I++) {
        workd[I+zInd].print();
       }
       printf("\n");*/
}       
   
    if ((ido == 3) || (ido == 4)) {
      printf("Not supplied mode!!!\n");
      exit(0);
    }
    
    if (ido == 99) {
      break;
    }    
  }

printf("fac size: %d\n",iparam[4]);

/* 
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
*/


  if (info<0) {
    printf("Error info code: %d\n",info);
  } else { 
/*
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using ZNEUPD .                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        |                                           |
c        | The routine ZNEUPD  now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1.)                                   |
c        |                                           |
c        %-------------------------------------------%
*/          
    rvec = false;
    char howMany = 'A';

    zneupd_(&rvec, &howMany, select, d, v, &ldv, &sigma,
            workev, &(bmat[0]), &n, &(which[0]), &nev, &tol, resid, &ncv,
            v, &ldv, iparam, ipntr, workd, workl, &lworkl,
           rwork, &ierr);

/*
c        %-----------------------------------------------%
c        | Eigenvalues are returned in the one           |
c        | dimensional array D and the corresponding     |
c        | eigenvectors are returned in the first        |
c        | NCONV (=IPARAM(5)) columns of the two         |
c        | dimensional array V if requested.  Otherwise, |
c        | an orthogonal basis for the invariant         |
c        | subspace corresponding to the eigenvalues in  |
c        | D is returned in V.                           |
c        %-----------------------------------------------%
*/
  
  printf("Eigenvalues found: \n");
  for (I=0; I<nev; I++) {
    d[I].print();
  }
  
    if (ierr!=0) {
      printf("Error in zneupd: %d\n",ierr);
    } else {


/*
c        %-------------------------------------------%
c        | Print additional convergence information. |
c        %-------------------------------------------%
*/

      if (info==1) {
        printf("Maximum number of iterations reached.\n");
      } 
      if (info==3) {
        printf("No shifts could be applied during implicit Arnoldi update, try increasing NCV.\n");
      }
    }
  }
}


int main(int argc,char **argv) {
  eigenVals();
}
