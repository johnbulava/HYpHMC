#include "PolynomialApproximation.h"

PolynomialApproximation::PolynomialApproximation(int subPolCnt, int id, int* deg, int digit, double alpha, double M1, double M2, double* eps, double* lam) {
  ini(subPolCnt, id, deg, digit, alpha, M1, M2, eps, lam);
  calcApproxPolys();
}


// Precision of cln::cl_F corresponding to MAXPOW.
//
// A good guess is:                          DIGIT = 70+2.8*MAXPOW
// See hep-lat/0302025
// but one has to check this by two runs with increasing precision.
void PolynomialApproximation::ini(int subPolCnt, int id, int* deg, int digit, double alpha, double M1, double M2, double* eps, double* lam) {
  if (LogLevel>1) printf("Initializing PolynomialApproximation with maxPow = %d, digit = %d, alpha = %f, Mass1 = %f, and Mass2 = %f\n", deg[0],digit,alpha, M1, M2);
  ID = id;
  DIGIT = digit;
  subPolyCount = subPolCnt;
  int I;
  polyDegree = new int[1+subPolyCount];
  for (I=0; I<1+subPolyCount; I++) {
    polyDegree[I] = deg[I];
    if ((polyDegree[I] % 2) == 1) {
      printf("ERROR: Only even polynomial degrees allowed!!!\n");
      exit(0);
    }
  }
  if (DIGIT<=0) DIGIT = (int) (70+2.8*polyDegree[0]);
  clnDIGIT = cln::float_format(DIGIT);  
  
  char* xxxStr = new char[1000];

  snprintf(xxxStr,1000,"1.0e+0_%d",DIGIT);
  if (LogLevel>3) printf("Initializing ONE with: %s\n",xxxStr);
  ONE = xxxStr;

  snprintf(xxxStr,1000,"2.0e+0_%d",DIGIT);
  if (LogLevel>3) printf("Initializing TWO with: %s\n",xxxStr);
  TWO = xxxStr;

  snprintf(xxxStr,1000,"0.0e+0_%d",DIGIT);
  if (LogLevel>3) printf("Initializing ZERO with: %s\n",xxxStr);
  ZERO = xxxStr;

  snprintf(xxxStr,1000,"0.5e+0_%d",DIGIT);
  if (LogLevel>3) printf("Initializing HALF with: %s\n",xxxStr);
  HALF = xxxStr;

  snprintf(xxxStr,1000,"100.0e+0_%d",DIGIT);
  if (LogLevel>3) printf("Initializing HUND with: %s\n",xxxStr);
  HUND = xxxStr;
  
  snprintf(xxxStr,1000,"%1.10fe+0_%d",alpha,DIGIT);
  if (LogLevel>3) printf("Initializing ALPHA with: %s\n",xxxStr);
  ALPHA = xxxStr;

double WeightExp = 0.0;
  snprintf(xxxStr,1000,"%1.10fe+0_%d",WeightExp, DIGIT);
  if (LogLevel>3) printf("Initializing WEIGHTEXP with: %s\n",xxxStr);
  WEIGHTEXP = xxxStr;

  snprintf(xxxStr,1000,"%1.10fe+0_%d",M1,DIGIT);
  if (LogLevel>3) printf("Initializing Mass1 with: %s",xxxStr);
  Mass1 = xxxStr;
  Mass1Zero = (M1==0);
  if (LogLevel>3) printf("  ==> Mass1Zero = %d\n",Mass1Zero);
  
  snprintf(xxxStr,1000,"%1.10fe+0_%d",M2,DIGIT);
  if (LogLevel>3) printf("Initializing Mass2 with: %s",xxxStr);
  Mass2 = xxxStr;
  Mass2Zero = (M2==0);
  if (LogLevel>3) printf("  ==> Mass2Zero = %d\n",Mass2Zero);
  
  EPSILON = new cln::cl_F[1+subPolyCount];
  LAMBDA  = new cln::cl_F[1+subPolyCount];
  for (I=0; I<1+subPolyCount; I++) {
    snprintf(xxxStr,1000,"%1.10fe+0_%d",eps[I],DIGIT);
    if (LogLevel>3) printf("Initializing EPSILON for subpolynomial %d (degree %d) with: %s\n",I,polyDegree[I],xxxStr);
    EPSILON[I] = xxxStr;

    snprintf(xxxStr,1000,"%1.10fe+0_%d",lam[I],DIGIT);
    if (LogLevel>3) printf("Initializing LAMBDA for subpolynomial %d (degree %d) with: %s\n",I,polyDegree[I],xxxStr);
    LAMBDA[I] = xxxStr;
  }

  delete[] xxxStr;
  

  Recb = new cln::cl_F*[1+subPolyCount];
  Recg = new cln::cl_F*[1+subPolyCount];
  Orth = new cln::cl_F*[1+subPolyCount];
  Coed = new cln::cl_F*[1+subPolyCount];
  for (I=0; I<1+subPolyCount; I++) {
    Recb[I] = new cln::cl_F[1+polyDegree[I]-1];
    Recg[I] = new cln::cl_F[1+polyDegree[I]-2];
    Orth[I] = new cln::cl_F[1];
    Coed[I] = new cln::cl_F[1+polyDegree[I]];
  }

  RootRepTotalNorm  = new cln::cl_F[1+2*subPolyCount];
  RootRepSingleNorm = new cln::cl_F[1+2*subPolyCount]; 
  Coef = new cln::cl_F*[1+2*subPolyCount];
  Roots = new cln::cl_N*[1+2*subPolyCount];
  RootRoots = new cln::cl_N*[1+2*subPolyCount];
  for (I=0; I<1+2*subPolyCount; I++) {
    RootRepTotalNorm[I] = ZERO;
    RootRepSingleNorm[I] = ZERO;
    Coef[I] = new cln::cl_F[1+polyDegree[I/2]];
    Roots[I] = new cln::cl_N[polyDegree[I/2]];
    RootRoots[I] = new cln::cl_N[2*polyDegree[I/2]];
  }  
}



/******************************************************************************/
//
//  Calculate the basic integral s_nu for least-square optimization.
//
//  Input parameters:
//  order of the polynomial:                               Maxpow
//
//  lower bound of interval of approximation:              Epsilon
//  upper bound of interval of approximation:              Lambda
//
//  The result is put in                                   Sint
void PolynomialApproximation::BaseIntS(int Maxpow, cln::cl_F Epsilon, cln::cl_F Lambda, cln::cl_F* Sint) {
  int  ord;
  cln::cl_F power, small = As(cln::cl_F)(expt(cln::cl_float(0.1,clnDIGIT),DIGIT));

  if (Mass1Zero && Mass2Zero) {
    for(ord = 0; ord < 2*Maxpow+1; ord++) { 
      power = As(cln::cl_F)(2*ALPHA+2*WEIGHTEXP+(ord+1));

      if(abs(power) < small) {
        Sint[ord] = As(cln::cl_F)(log(Lambda/Epsilon));
      } else {
        Sint[ord] = As(cln::cl_F)(expt(Lambda,power)-expt(Epsilon,power))/power; 
      }
    } 
  } else {
    printf("Prec Masses not implemented yet!!!\n");
    exit(0);
  }
}



/******************************************************************************/
//
//  Calculate the basic integral t_nu for least-square optimization.
//
//  Input parameters:
//  order of the polynomial:                               Maxpow
//
//  lower bound of interval of approximation:              Epsilon
//  upper bound of interval of approximation:              Lambda
//
//  The result is put in                                   Tint
void PolynomialApproximation::BaseIntT(int Maxpow, cln::cl_F Epsilon, cln::cl_F Lambda, cln::cl_F* Tint) {
  int  ord;
  cln::cl_F power, small = As(cln::cl_F)(expt(cln::cl_float(0.1,clnDIGIT),DIGIT));;

  if (Mass1Zero && Mass2Zero) {
    for(ord = 0; ord < Maxpow+1; ord++) { 
      power = As(cln::cl_F)(ALPHA+2*WEIGHTEXP+(ord+1));

      if(abs(power) < small) {
        Tint[ord] = As(cln::cl_F)(log(Lambda/Epsilon));
      } else {
        Tint[ord] = As(cln::cl_F)(expt(Lambda,power)-expt(Epsilon,power))/power; 
      }
    } 
  } else {
    printf("Prec Masses not implemented yet!!!\n");
    exit(0);
  }
}
 

 
/******************************************************************************/
//
//  Evaluate the approximate polynomial up to the order        Maxpow
//  at the variable value                                      Xval
//
//  The recurrence coefficients are in                         Recb
//  and in                                                     Recg
//  Constant term of the first orthogonal polynomial:          Orth
//  Expansion coefficients in orthogonal polynomials:          Coed
cln::cl_F PolynomialApproximation::Recurev(int Maxpow, cln::cl_F Xval, cln::cl_F* Recb, cln::cl_F* Recg, cln::cl_F* Orth, cln::cl_F* Coed) {
   int  ord;
   cln::cl_F res, orth, orthb, orthg;


// Check input

    if(Maxpow < 0)
     { cout <<endl <<"wrong order in Recurev:" <<endl;
       cout <<endl <<Maxpow <<endl; }

// Start iteration

   orth = As(cln::cl_F)(ONE);
   res = Coed[0]*orth;

    if(Maxpow > 0)
     { orthb = orth;
       orth = Xval+Orth[0];
       res = res+Coed[1]*orth; }

// Iteration for recurrence

    for(ord = 2; ord < Maxpow+1; ord++)
     { orthg = orthb;
       orthb = orth;
       orth = (Xval+Recb[ord-1])*orthb+Recg[ord-2]*orthg;
       res = res+Coed[ord]*orth; }

   return res;
}
 
 
 
/******************************************************************************/
//
// Write elements of a list of numbers in the array       Wlist
// of length                                              Leng
// as assignment statements  
// to the elements of an array                            Arrayname
//
// The ofstream of the output file is                     Ostream
//  
// At the beginning write the string                      Text
void PolynomialApproximation::WriteAssign(ostream& Ostream, char* Text, cln::cl_F* Wlist, int Leng, char* Arrayname) {
   int   ord;
   Ostream <<endl <<Text <<endl <<endl;

   Ostream.setf(ios::scientific);
   Ostream.precision(16);

   for(ord = 0; ord < Leng; ord++) {  
     Ostream <<"   " <<Arrayname <<"[" <<ord <<"] = " <<double_approx(Wlist[ord]) <<endl;
   }
}



void PolynomialApproximation::WriteRoots(char * filename, cln::cl_N *list, const int length) {
  ofstream ofs(filename, ios::out);
  ofs.setf(ios::scientific);
  ofs.precision(16);
  ofs << "Nr.           Re            Im" << endl;
  for(int i = 0; i < length; i++) {
    ofs << i << " " <<
      double_approx(realpart(list[i])) << " " <<
      double_approx(imagpart(list[i])) << endl;
  }
  ofs.close();
}



void PolynomialApproximation::WriteCoefs(cln::cl_N *list, const int length) {
  ofstream ofs("coefs.dat", ios::out);
  ofs.setf(ios::scientific);
  ofs.precision(16);
  for(int i = 0; i < length; i++) {
    ofs << i << " " <<
      double_approx(realpart(list[i])) << endl;
  }
  ofs.close();
}



/******************************************************************************/
//
//  Calculating the polynomial approximation minimizing the integral
//  of quadratic deviation from x^(-Alpha) in an interval.
//
//  The weight function in least-squares optimization
//  is the inverse of the function.
//
//  Recursion relations for orthogonal polynomials are used
//  whenever possible.
//
//  Input parameters:
//  order of the optimized polynomial:                     Maxpow  
//  lower bound of interval of approximation:              Epsilon
//  upper bound of interval of approximation:              Lambda
//
//                                                         Recb
//                                                         Recg
//                                                         Orth
//                                                         Coed
//  name of the file for writing results:                  Filename
void PolynomialApproximation::Quadropt(int Maxpow, cln::cl_F Epsilon, cln::cl_F Lambda, cln::cl_F* Recb, cln::cl_F* Recg, cln::cl_F* Orth, cln::cl_F* Coed, char* Filename) {
   int  ord, or1, or2, leng;
   cln::cl_F fmu, quad, scal, p1e, p1l;

   cln::cl_F sint[1+2*Maxpow],
        norq[1+Maxpow],
        bint[1+Maxpow],
        tint[1+Maxpow],

        rmup[1+2*Maxpow-2],
        rmuv[1+2*Maxpow-1],
        rmum[1+2*Maxpow-2],

        bmup[1+Maxpow-2],
        bmuv[1+Maxpow-1],
        bmum[1+Maxpow-2];


// Define necessary basic integrals of the weight function

   BaseIntS(Maxpow, Epsilon, Lambda, sint);

// Start iteration for determining the orthogonal polynomials
  if(LogLevel>2) {
    cout  <<"Quadropt: Calculating approximation polynomial for Alpha=" <<ALPHA <<", Epsilon=" << Epsilon
    <<", and maxpow="<<Maxpow<<endl;
  }

  if(LogLevel>3) {
    printf("Quadropt: Calculating recurrence constants for orthogonal polynomials:\n");
  }

   norq[0] = sint[0];

    if(Maxpow > 0)
     { Recb[0] = -sint[1]/sint[0];
       Orth[0] = -sint[1]/sint[0];
       norq[1] = sint[2]-expt(sint[1],2)/sint[0]; }

    if(Maxpow > 1)
     {  for(ord = 0; ord < 2*Maxpow-1; ord++)
        rmum[ord] = sint[ord];

        for(ord = 0; ord < 2*Maxpow; ord++)
        rmuv[ord] = sint[ord+1]+Orth[0]*sint[ord];

       fmu = Orth[0];
       Recb[1] = -fmu-rmuv[2]/rmuv[1];
       Recg[0] = -rmuv[1]/rmum[0]; }

// Iteration for orthogonal polynomials

     for(or1 = 1; or1 < Maxpow; or1++) { 
        Recb[or1] = -fmu-rmuv[or1+1]/rmuv[or1];
        Recg[or1-1] = -rmuv[or1]/rmum[or1-1];
        fmu = fmu+Recb[or1];

         for(or2 = 0; or2 < 2*Maxpow-or1; or2++)
         rmup[or2] = rmuv[or2+1]+Recb[or1]*rmuv[or2]+Recg[or1-1]*rmum[or2];

        norq[or1+1] = rmup[or1+1];

         for(or2 = 0; or2 < 2*Maxpow-or1-1; or2++)
         rmum[or2] = rmuv[or2];

         for(or2 = 0; or2 < 2*Maxpow-or1; or2++)
         rmuv[or2] = rmup[or2];
      }

// Calculate orthogonal expansion coefficients of optimized polynomial
// Calculate basic integrals

   BaseIntT(Maxpow, Epsilon, Lambda, tint);

// Start iteration

  if(LogLevel>3) {
    printf("Quadropt: Calculating expansion coefficients of optimized polynomials:\n");
  }

   bint[0] = tint[0];
   Coed[0] = bint[0]/norq[0];

    if(Maxpow > 0)
     { bint[1] = tint[1]+Orth[0]*tint[0];
       Coed[1] = bint[1]/norq[1]; }

    if(Maxpow > 1)
     {  for(ord = 0; ord < Maxpow-1; ord++)
        bmum[ord] = tint[ord];

        for(ord = 0; ord < Maxpow; ord++)
        bmuv[ord] = tint[ord+1]+Orth[0]*tint[ord]; }

// Perform iteration

    for(or1 = 1; or1 < Maxpow; or1++) {
        for(or2 = 0; or2 < Maxpow-or1; or2++)
        bmup[or2] = bmuv[or2+1]+Recb[or1]*bmuv[or2]+Recg[or1-1]*bmum[or2];

       bint[or1+1] = bmup[0];
       Coed[or1+1] = bint[or1+1]/norq[or1+1];

        for(or2 = 0; or2 < Maxpow-or1-1; or2++)
        bmum[or2] = bmuv[or2];

        for(or2 = 0; or2 < Maxpow-or1; or2++)
        bmuv[or2] = bmup[or2];

    }

// Print out numerical value at the minimum found

  if(LogLevel>3) {
    cout  <<"...minimum calculated with Digits [" << clnDIGIT <<"]" <<endl;
  }

   quad = As(cln::cl_F)(ONE);
   scal = Lambda-Epsilon;

    for(ord = 0; ord < Maxpow+1; ord++)
    quad = quad-Coed[ord]*bint[ord]/scal;

    if(LogLevel>2) { 
      cout <<"...quadratic deviation: " <<cln::double_approx(quad) <<endl;
      cout <<"...on interval of weight  " <<cln::double_approx(scal) <<endl; 
    }

// Check interval ends

  p1e = Recurev(Maxpow,Epsilon,Recb,Recg,Orth,Coed);
  p1e = p1e*As(cln::cl_F)(expt(Epsilon,ALPHA))-As(cln::cl_F)(ONE);

  p1l = Recurev(Maxpow,Lambda,Recb,Recg,Orth,Coed);
  p1l = p1l*As(cln::cl_F)(expt(Lambda,ALPHA))-As(cln::cl_F)(ONE);

  if(LogLevel>2) { 
    cout <<"...relative deviation of polynomial at left interval end:  ";
    cout <<cln::double_approx(p1e) <<endl;

    cout <<"...relative deviation of polynomial at right interval end:  ";
    cout <<cln::double_approx(p1l) <<endl; 
  }
  

// Write coefficients of the polynomial to file

   ofstream Ostream(Filename, ios::out);

   Ostream <<endl <<" Maxpow, Alpha, Epsilon, Lambda:" <<endl;
   Ostream <<endl <<"   " <<Maxpow <<"   " <<ALPHA 
            <<"   " <<Epsilon <<"   " <<Lambda <<endl;
   Ostream <<endl <<" calculated with Digits [" << clnDIGIT<<"]" <<endl;

   Ostream <<endl <<" Quadratic deviation: " <<quad <<endl;
   Ostream <<endl <<" on interval of weight  " <<scal <<endl;

   Ostream <<endl <<" Relative deviation of polynomial at interval ends:"
                <<endl;
   Ostream <<endl <<p1e <<endl;
   Ostream <<endl <<p1l <<endl <<endl;

   leng = 1;
   WriteAssign(Ostream," Result for starting recurrence coefficient:",
               Orth,leng,"Orth1");

   leng = Maxpow+1;
   WriteAssign(Ostream," Result for expansion coefficients:",
               Coed,leng,"Coed1");

   leng = Maxpow;
   WriteAssign(Ostream," Result for first coefficients for recurrence:",
               Recb,leng,"Recb1");

   leng = Maxpow-1;
   WriteAssign(Ostream," Result for second coefficients for recurrence:",
               Recg,leng,"Recg1");

   if (LogLevel>2) {
     printf("...polynomial is ready.\n");
   }
}
 
 
 
/******************************************************************************/
//
//  Restore polynomial from its expansion in orthogonal polynomials.
//
//  The input coefficients will be left unchanged.
//
//  Input parameters:
//
//  Maximal order of the polynomial:                          Maxpow
//
//  The recurrence coefficients are in                        Recb
//  and in                                                    Recg
//  Constant term of the first orthogonal polynomial:         Orth
//  Expansion coefficients in orthogonal polynomials:         Coed
//
//  the obtained coefficients of the polynomial:              Coef
void PolynomialApproximation::RestorePolyCoef(int Maxpow, cln::cl_F* Recb, cln::cl_F* Recg, cln::cl_F* Orth, cln::cl_F* Coed, cln::cl_F* Coef) {
   int ord, orv;

   cln::cl_F poly[1+Maxpow];
   cln::cl_F opl[1+Maxpow], oplm[1+Maxpow], oplp[1+Maxpow];


  if (LogLevel>2) {
    printf("Restoring polynomial from orthogonal representation...");
  }

// Check input

    if(Maxpow < 0) { 
      printf("wrong order in RestorePolyCoef: %d\n",Maxpow);
    }

// Start iteration

    for(ord = 0; ord < Maxpow+1; ord++)
     { poly[ord] = As(cln::cl_F)(ZERO);
       opl[ord]  = As(cln::cl_F)(ZERO);
       oplm[ord] = As(cln::cl_F)(ZERO);
       oplp[ord] = As(cln::cl_F)(ZERO); }

   oplm[0] = As(cln::cl_F)(ONE);

   opl[0] = Orth[0];
   opl[1] = As(cln::cl_F)(ONE);

   poly[0] = Coed[0]*oplm[0]+Coed[1]*opl[0];
   poly[1] = Coed[0]*oplm[1]+Coed[1]*opl[1];

// Loop over order of orthogonal polynomial

  for(ord = 1; ord < Maxpow; ord++) {

    for(orv = 0; orv < Maxpow+1; orv++) { 
      oplp[orv] = Recb[ord]*opl[orv]+Recg[ord-1]*oplm[orv];
      if(orv > 0) oplp[orv] = oplp[orv]+opl[orv-1]; 
    }

    for(orv = 0; orv < Maxpow+1; orv++) {
      poly[orv] = poly[orv]+Coed[ord+1]*oplp[orv];
    }

    for(orv = 0; orv < Maxpow+1; orv++) {
      oplm[orv] = opl[orv];
    }

    for(orv = 0; orv < Maxpow+1; orv++) {
      opl[orv] = oplp[orv]; 
    }
  }

// Extract coefficients

    for(ord = 0; ord < Maxpow+1; ord++)
    Coef[ord] = As(cln::cl_F)(poly[ord]);

  if (LogLevel>2) {
    printf("READY\n\n");
  }
}
 
 
 
/******************************************************************************/
//
//  Evaluate a complex polynomial
//
//  The polynomial is                                    Poly
//  The order is                                         Maxpow
//
//  The value of the variable:                           Valu
cln::cl_N PolynomialApproximation::EvalPoly(cln::cl_N* Poly, int Maxpow, cln::cl_N Valu) {
   int  pow;
   cln::cl_N xpow, sum;


   sum = As(cln::cl_N)(complex(ZERO,ZERO));
   xpow = As(cln::cl_N)(complex(ONE,ZERO));

    for(pow = 0; pow < Maxpow+1; pow++)
     { sum = sum+xpow*Poly[pow];
       xpow = xpow*Valu; }

   return sum;
}



// this routine evaluated the product representation of a
// polynomial with roots roots, order Maxpow at value x and
// with (overall) normalisation factor norma
cln::cl_N PolynomialApproximation::EvalPolyProd(cln::cl_N* roots, const int Maxpow, const cln::cl_N x) {
  cln::cl_N norma = complex(ONE, ZERO);

  cln::cl_N prod = As(cln::cl_N)(norma);
  for(int i = 0; i < Maxpow; i++) {
    prod = prod * (x-roots[i]);
  }
  
  return(prod);
}



/******************************************************************************/
//
//  Find a root of a complex polynomial by Laguerre iteration.
//
//  The polynomial is                                    Poly
//  The order is                                         Maxpow
//
//  The precision:                                       Digit
cln::cl_N PolynomialApproximation::Lasolv(cln::cl_N* Poly, int Maxpow, cln::cl_N root, int itemax) {
   int  pow, ite;
   
   root = complex(ZERO,ZERO);
   
   cln::cl_F angl, small = As(cln::cl_F)(expt(cln::cl_float(0.1,clnDIGIT),DIGIT/2));

   cln::cl_N dif1[Maxpow], dif2[Maxpow-1];
   cln::cl_N val0, val, val1, val2, denp, denm, las1, las2, sqrv;
   //   cln::cl_N root;
    for(pow = 0; pow < Maxpow; pow++)
    dif1[pow] = (pow+1)*Poly[pow+1];

    for(pow = 0; pow < Maxpow-1; pow++)
    dif2[pow] = (pow+1)*dif1[pow+1];

// The maximal allowed number of iterations is set here;
// this can be chosen larger, but 100 usually suffices

//   root = As(cln::cl_N)(complex(ZERO,ZERO));
   val0 = EvalPoly(Poly,Maxpow,root);

// Iteration

    for(ite = 0; ite < itemax; ite++)
     { 
       val = val0;
       val1 = EvalPoly(dif1,Maxpow-1,root);
       val2 = EvalPoly(dif2,Maxpow-2,root);

       sqrv = (Maxpow-1)*((Maxpow-1)*val1*val1-Maxpow*val0*val2);
       angl = HALF*cln::cl_float(phase(sqrv),clnDIGIT);
       sqrv = sqrt(abs(sqrv))*complex(cos(angl),sin(angl));
       denp = val1+sqrv;
       denm = val1-sqrv;

        if(denp == complex(ZERO,ZERO))
        root = root-Maxpow*val0/denm;

        else
         {  if(denm == complex(ZERO,ZERO))
            root = root-Maxpow*val0/denp;

            else
             { las1 = -Maxpow*val0/denp;
               las2 = -Maxpow*val0/denm;

                if(realpart(las1*conjugate(las1)) <
                   realpart(las2*conjugate(las2)))
                root = root+las1;

                else
                root = root+las2; } }

//  Look whether the root is good enough

       val0 = EvalPoly(Poly,Maxpow,root);

        if(abs(val0) == ZERO ||
          (abs(val0) < small) && abs(val0/val) > 0.7)
         {
            if(LogLevel>4) 
             { cout << endl << "Laguerre iterations: " << ite << endl;
               cout << endl << "root = " << root << endl;
               cout << endl << "value at root: " << val0 << endl; }

           break; } }

    if(ite >= itemax) {
      cout <<endl << "Laguerre iteration did not converge" <<endl;
      exit(5);
    }

   return root;

}
 
 
 
/******************************************************************************/
//
//  Find the complex roots of a complex polynomial       Poly
//  The order is                                         Maxpow
//  The result will be in the array                      Root
//
void PolynomialApproximation::Polyrootc(cln::cl_N* Poly, int Maxpow, cln::cl_N* Root) {
  
  int  ord, pow, fnd, pov, maxp;

  cln::cl_N poly[1+Maxpow], polc[1+Maxpow], coef[1+Maxpow], coen[1+Maxpow];


  // Put coefficients in an array

  for(pow = 0; pow < Maxpow+1; pow++) { 
    poly[pow] = As(cln::cl_N)(Poly[pow]);
    coef[pow] = As(cln::cl_N)(Poly[pow]); 
  }

  for(pow = 0; pow < Maxpow+1; pow++) {
    polc[pow] = As(cln::cl_N)(complex(ZERO,ZERO));
  }

  polc[0] = As(cln::cl_N)(complex(ONE,ZERO));
  fnd = -1;

  // Loop for finding all roots

  for(ord = 0; ord < Maxpow; ord++) {
    fnd++;
    pov = Maxpow-fnd;

    if(fnd < Maxpow) {
      if(LogLevel>4) {
	cout <<endl <<" root number: " <<fnd+1 <<endl;
      }

      if((ord%2 == 1) && (Maxpow%2 == 0) && (false)) {
	Root[fnd] = As(cln::cl_N)(conjugate(Root[fnd-1]));
 	cln::cl_N val0 = EvalPoly(poly, Maxpow, Root[fnd]);
	
 	if(LogLevel>3) {
 	  cout << endl << "root = " << Root[fnd] << endl;
 	  cout << endl << "value at root: " << val0 << endl; 
 	}
      }
      else {
	Root[fnd] = Lasolv(poly,pov, complex(ONE,ONE), 150);
      }

      for(pow = Maxpow; pow > 0; pow--) {
	polc[pow] = polc[pow-1]-Root[fnd]*polc[pow];
      }

      polc[0] = -Root[fnd]*polc[pow];

      // Divide the polynomial by the root

      maxp = Maxpow-fnd-1;
      coen[maxp] = coef[maxp+1];

      for(pow = maxp-1; pow > -1; pow--) {
	coen[pow] = coef[pow+1]+Root[fnd]*coen[pow+1];
      }

      for(pow = 0; pow < maxp+1; pow++) {
	coef[pow] = coen[pow];
	poly[pow] = coef[pow]; 
      } 
    }

    else {
      break;
    }
  }

// Compare input with product of root factors

  for(pow = 0; pow < Maxpow+1; pow++) {
    polc[pow] = Poly[pow]-poly[0]*polc[pow];
  }

  if(LogLevel>4) {
    cout <<endl <<"control polynomial should be close to zero:" <<endl;

    for(pow = 0; pow < Maxpow+1; pow++) {
      cout <<endl <<"  x^{" <<pow <<"}" <<endl;
      cout <<endl <<polc[pow] <<endl; 
    } 
  }
}



/******************************************************************************/
//
// Find the optimal order of the roots in             Root
// for calculating the polynomial with order          Maxpow
// by multiplication of the root factors.
//
// The overall factor in the polynomial is:           Coef0
// The function to be approximated is                 x^(-Alpha)
// Optimization is done in the interval               [Epsilon,Lambda]
//
// Check values in                                    Ncheck
// equidistant points of the interval.
//
// This calculation is done with precision            digit
void PolynomialApproximation::Optimord(int Maxpow, cln::cl_N* Root, cln::cl_F Coef0, cln::cl_F Epsilon, cln::cl_F Lambda, int Ncheck) {
   int nmax = Ncheck, digval = 40;
   cln::float_format_t digit = cln::float_format(digval);

   int nl, rr, rt, rv, rc, facc[Maxpow];

// Check input

    if(Ncheck < 3)
     { cout <<endl <<"number of checkpoints in Optimord should be at least 3"
            <<endl;
       nmax = 3; }

    if(Ncheck > 47)
     { cout <<endl <<"with this value of Ncheck it would take somewhat longer"
            <<endl;
       cout <<endl <<"the number of checkpoints will be set to 47" <<endl;
       nmax = 47; }

   cln::cl_F large = As(cln::cl_F)(expt(cln::cl_float(10.0,digit),digval));
   cln::cl_F small = As(cln::cl_F)(expt(cln::cl_float( 0.1,digit),digval));

   cln::cl_F xx[nmax], cc[nmax], values[Maxpow][nmax], mx[Maxpow];
   cln::cl_F max, min, mn;

   cln::cl_N root[Maxpow], ff[Maxpow][nmax], pp[nmax];


// Define the points where comparison will be made

    for(nl = 0; nl < nmax; nl++)
    xx[nl] = cln::cl_float(Epsilon*(nmax-nl-1)/(nmax-1)+Lambda*nl/(nmax-1),digit);

// Special case: Epsilon=0

    if(xx[0] == 0.0)
    xx[0] = xx[1]/100;

// The values to be compared

    for(nl = 0; nl < nmax; nl++)
    cc[nl] = cln::cl_float(Coef0*exp(ln(xx[nl])*ALPHA),digit);

// Loop over the root factors of the polynomial

    for(rr = 0; rr < Maxpow; rr++)
    root[rr] = complex(cln::cl_float(realpart(Root[rr]),digit),
                       cln::cl_float(imagpart(Root[rr]),digit));

    for(rr = 0; rr < Maxpow; rr++)
     { facc[rr] = 1;

        for(nl = 0; nl < nmax; nl++)
        ff[rr][nl] = xx[nl]-root[rr]; }

// Find optimal order by minimizing maximal ratio in partial products

    for(nl = 0; nl < nmax; nl++)
    pp[nl] = complex(cln::cl_float(1,digit),cln::cl_float(0,digit));

    if(LogLevel>3)
     { cout <<endl <<endl <<"  The logarithm of the overall constant:" <<endl;
       cout <<endl <<cln::cl_float(ln(abs(Coef0)),digit) <<endl;
       cout <<endl <<"  Maximal logarithm of ratios in partial products:"
            <<endl; }

// Loop over the roots

   rc = 0;

    for(rt = 0; rt < Maxpow; rt++)
     {
        for(rr = 0; rr < Maxpow; rr++)
         if(facc[rr] == 1)
          { max = small;
            min = large;

             for(nl = 0; nl < nmax; nl++)
              { values[rr][nl] =
                 cln::cl_float(ln(abs(cc[nl]*pp[nl]*ff[rr][nl])),digit);

                if(max < values[rr][nl]) max = values[rr][nl];
                if(min > values[rr][nl]) min = values[rr][nl]; }

            mx[rr] = max-min; }

       mn = large;
       rv = -1;

        for(rr = 0; rr < Maxpow; rr++)
         if(facc[rr] == 1)
          if(mn > mx[rr])
          { mn = mx[rr];
            rv = rr; }

        if(LogLevel>3)
        cout <<endl <<rt+1 <<":  " <<mn <<endl;

       Root[rc] = root[rv];
       rc++;
       facc[rv] = 0;

        for(nl = 0; nl < nmax; nl++)
        pp[nl] = pp[nl]*ff[rv][nl]; }

    if(LogLevel>3)
     { cout <<endl <<endl <<"  Product check of ratios:" <<endl;

        for(nl = 0; nl < nmax; nl++)
         { cout <<endl <<"  x = " <<xx[nl] <<" :" <<endl;
           cout <<endl <<cc[nl]*pp[nl] <<endl; } }

}
 
 
 
/******************************************************************************/
//
// Find the optimal order of the roots in             Root
// for calculating the polynomial with order          Maxpow
// by multiplication of the root factors.
//
// Complex conjugate pairs of roots are kept together:
// positive imaginary part first in pair.
//
// The overall factor in the polynomial is:           Coef0
// The function to be approximated is                 x^(-Alpha)
// Optimization is done in the interval               [Epsilon,Lambda]
//
// Check values in                                    Ncheck
// equidistant points of the interval.
//
// This calculation is done with precision            digit
void PolynomialApproximation::OptimordPair(int Maxpow, cln::cl_N* Root, cln::cl_F Coef0, cln::cl_F Epsilon, cln::cl_F Lambda, int Ncheck){
  int nmax = Ncheck, digval = 40;
  cln::float_format_t digit = cln::float_format(digval);

  int nl, rv, rn, r1, r2, rc, rt, rr, maxp;
  int facp[Maxpow], facr[Maxpow], facl[Maxpow];

  // Check input

  if(Ncheck < 3)
    { cout <<endl <<"number of checkpoints in Optimord should be at least 3"
	   <<endl;
      nmax = 3; }

  if(Ncheck > 137)
    { cout <<endl <<"with this value of Ncheck it would take somewhat longer"
	   <<endl;
      cout <<endl <<"the number of checkpoints will be set to 137" <<endl;
      nmax = 137; }

  cln::cl_F large = As(cln::cl_F)(expt(cln::cl_float(10.0,digit),digval));
  cln::cl_F small = As(cln::cl_F)(expt(cln::cl_float( 0.1,digit),digval));

  cln::cl_F xx[nmax], cc[nmax], values[Maxpow][nmax], mx[Maxpow];
  cln::cl_F max, min, mn;

  cln::cl_N rot[Maxpow], root[Maxpow], ff[Maxpow][nmax], pp[nmax], prd;


  // Define the points where comparison will be made

  for(nl = 0; nl < nmax; nl++)
    xx[nl] = cln::cl_float(Epsilon*(nmax-nl-1)/(nmax-1)+Lambda*nl/(nmax-1),digit);

  // Special case: Epsilon=0

  if(xx[0] == 0.0)
    xx[0] = xx[1]/100;

  // The values to be compared

  for(nl = 0; nl < nmax; nl++)
    cc[nl] = cln::cl_float(Coef0*exp(ln(xx[nl])*ALPHA),digit);

  for(rr = 0; rr < Maxpow; rr++)
    rot[rr] = complex(cln::cl_float(realpart(Root[rr]),digit),
                      cln::cl_float(imagpart(Root[rr]),digit));

  // Pair up indices of complex cojugate roots: 
  // positive imaginary part first in pair
  // Store indices of pairs and real roots

  for(rr = 0; rr < Maxpow; rr++)
    { facr[rr] = 1;
      facl[rr] = 0; }

  rn = 0;

  for(rr = 0; rr < Maxpow; rr++)
    if(abs(imagpart(rot[rr])) < 100*small)
      { facp[rn] = -rr;
        facr[rr] = 0;
        rn++; }

  // Check whether the complex roots are in pairs

  if((Maxpow-rn)%2 != 0)
    { cout <<endl <<"The number of complex roots is not even: " 
	   <<Maxpow-rn <<endl;
      return; }

  for(r1 = 0; r1 < Maxpow; r1++)
    if(facr[r1] == 1)
      {  for(r2 = r1+1; r2 < Maxpow; r2++)
          if(facr[r2] == 1)
	    {
              if(abs(conjugate(rot[r2])-rot[r1]) < 100*small)
		{
                  if(imagpart(rot[r1]) > 0)
		    { facp[rn] = r1;
		      facr[r1] = 0;
		      rn++;
		      facp[rn] = r2;
		      facr[r2] = 0;
		      rn++; }
                  else
		    { facp[rn] = r2;
		      facr[r2] = 0;
		      rn++;
		      facp[rn] = r1;
		      facr[r1] = 0;
		      rn++; }

		  break; } } }

  // Check whether the complex roots are in complex conjugate pairs

  if((Maxpow-rn) != 0)
    { cout <<endl <<"The roots are not in complex conjugate pairs: " 
	   <<Maxpow-rn <<endl;
      return; }

  // Reorder roots: first single real then complex pairs

  rc = 0;

  for(rn = 0; rn < Maxpow; rn++)
    {  if(facp[rn] < 0)
	{ root[rc] = rot[-facp[rn]];
	  facl[rc] = -1;
	  rc++; }
      else
	{ root[rc] = rot[facp[rn]];
	  facl[rc] = (imagpart(root[rc]) > 0);
	  rc++; } }

  // Calculate root factors

  for(rr = 0; rr < Maxpow; rr++)
    for(nl = 0; nl < nmax; nl++)
      ff[rr][nl] = xx[nl]-root[rr];

  for(nl = 0; nl < nmax; nl++)
    pp[nl] = complex(cln::cl_float(1,digit),cln::cl_float(0,digit));

  if(LogLevel>3)
    { cout <<endl <<endl <<"  The logarithm of the overall constant:" <<endl;
      cout <<endl <<cln::cl_float(ln(abs(Coef0)),digit) <<endl;
      cout <<endl <<"  Maximal logarithm of ratios in partial products:"
	   <<endl; }

  // Find optimal order by minimizing maximal ratio in partial products
  // Loop over complex conjugate root pairs and real roots

  maxp = 0;

  for(rr = 0; rr < Maxpow; rr++)
    { facr[rr] = 1;
      if(facl[rr] != 0) maxp++; }

  rc = 0;

  for(rt = 0; rt < maxp; rt++)
    {
      for(rr = 0; rr < Maxpow; rr++)
	if(facl[rr] != 0 && facr[rr] == 1)
          { max = small;
            min = large;

	    for(nl = 0; nl < nmax; nl++)
              { prd = ff[rr][nl];
                if(facl[rr] == 1) prd = prd*ff[rr+1][nl];

                values[rr][nl] =
		  cln::cl_float(ln(abs(cc[nl]*pp[nl]*prd)),digit);

                if(max < values[rr][nl]) max = values[rr][nl];
                if(min > values[rr][nl]) min = values[rr][nl]; }

            mx[rr] = max-min; }

      mn = large;
      rv = -1;

      for(rr = 0; rr < Maxpow; rr++)
	if(facl[rr] != 0 && facr[rr] == 1)
          if(mn > mx[rr])
	    { mn = mx[rr];
	      rv = rr; }

      if(LogLevel>3)
        cout <<endl <<rt+1 <<":  " <<mn <<endl;

      Root[rc] = root[rv];
      rc++;
      facr[rv] = 0;

      for(nl = 0; nl < nmax; nl++)
        pp[nl] = pp[nl]*ff[rv][nl]; 

      if(facl[rv] == 1)
	{ rv++;
	  Root[rc] = root[rv];
	  rc++;
	  facr[rv] = 0;

	  for(nl = 0; nl < nmax; nl++)
            pp[nl] = pp[nl]*ff[rv][nl]; } }

  if(LogLevel>3)
    { cout <<endl <<endl <<"  Product check of ratios:" <<endl;

      for(nl = 0; nl < nmax; nl++)
	{ cout <<endl <<"  x = " <<xx[nl] <<" :" <<endl;
	  cout <<endl <<cc[nl]*pp[nl] <<endl; } }

}



// Comparison function for naive ordering, See hep-lat/9805026
bool PolynomialApproximation::CompSortNaiv(const cln::cl_N & x, const cln::cl_N & y) {

  if((imagpart(x) <= ZERO) && (imagpart(y) <= ZERO)) {
    if(realpart(x) < realpart(y)) {
      return(true);
    }
  }
  if((imagpart(x) <= ZERO) && (plusp(imagpart(y)))) {
    return(true);
  }
  if((plusp(imagpart(x))) && (plusp(imagpart(y)))) {
    if(realpart(x) > realpart(y)) {
      return(true);
    }
  }
  return(false);
}



// This routine returns for idx in [0,2^digits-1] the corresponding
// bit reversal index 
int PolynomialApproximation::bitReversalRepresentation(const int idx, const int digits) {
  int res[digits];
  int num = idx;
  for(int i = 0; i < digits; i++) {
    res[i] = 0;
  }
  for(int i = 0; i < digits; i++) {
    int k = digits - i;
    int p = int(pow(2., (digits-(i+1))) );
    res[k] = num/p;
    num = num - res[k]*p;
  }
  num = 0;
  for(int i = 0; i < digits; i++) {
    num = num + res[i]*int(pow(2.,digits-(i+1)));
  }
  return(num);
}



void PolynomialApproximation::quicksort(const int n, cln::cl_N arr[], int idx[]) {
  cln::cl_N v, td;
  int i, j, l, r, ti, tos, stack[32];
  
  l = 0; r = n-1; tos = -1;
  for (;;) {
    while (r > l) {
      v = arr[r]; i = l; j = r-1;
      for (;;){
	while (CompSortNaiv(arr[i], v)) i ++;
	/* j > l prevents underflow */
	while (!CompSortNaiv(arr[j], v) && j > l) j --;
	if (i >= j) break;
	td = arr[i]; arr[i] = arr[j]; arr[j] = td;
	ti = idx[i]; idx[i] = idx[j]; idx[j] = ti;
      }
      td = arr[i]; arr[i] = arr[r]; arr[r] = td;
      ti = idx[i]; idx[i] = idx[r]; idx[r] = ti;
      if (i-l > r-i){
	stack[++tos] = l; stack[++tos] = i-1; l = i+1;
      }
      else{
	stack[++tos] = i+1; stack[++tos] = r; r = i-1;
      }
      if(tos > 31) {
	cerr << "Error in quicksort! Aborting...!" << endl;
	exit(31);
      }
    }
    if (tos == -1) break;
    r = stack[tos--]; l = stack[tos--];
  }
}



// This bring the complex list Roots of length Maxpow to
// naive order.
// See hep-lat/9805026
void PolynomialApproximation::NaiveOrder(cln::cl_N * Roots, const int Maxpow) {
  if (LogLevel>3) {
    printf("Bringing roots to naive order...\n");
  }

  int idx[Maxpow];

  for(int i =0; i < Maxpow; i++) {
    idx[i] = i;
  }
  quicksort(Maxpow, Roots, idx);
  if (LogLevel>3) {
    printf("Roots in naive order are:\n");
    for(int i = 0; i < Maxpow; i++) {
      cout << " -> " << i << ": "<< double_approx(realpart(Roots[i])) << " " << double_approx(imagpart(Roots[i])) << endl;
    }
  }

  if(LogLevel>3) {
    printf("READY\n");
  }
}


void PolynomialApproximation::enforceHermiticity(cln::cl_N *Roots, const int Maxpow) {
  if(LogLevel>3) {
    printf("Enforcing hermiticity of Polynomial...");
  }
  int I;
  for (I=0; I<Maxpow/2; I++) {
    Roots[2*I] = complex(HALF*realpart(Roots[2*I]+Roots[2*I+1]),HALF*imagpart(Roots[2*I]-Roots[2*I+1]));
    Roots[2*I+1] = complex(realpart(Roots[2*I]),imagpart(ZERO - Roots[2*I]));
  }
  if(LogLevel>3) {
    printf("READY\n");
  }  
}


void PolynomialApproximation::BitReversalOrderOfRootPairs(cln::cl_N *Roots, const int Maxpow) {
  if(LogLevel>3) {
    printf("Bringing root-pairs to Bit-reversed order...\n");
  }

  int digits = 2;
  int power = 2;
  while(power < (Maxpow/2)) {
    power=power*2;
    digits++;
  }
  
  cln::cl_N paddedRoots[2*power];
  cln::cl_N reversedRoots[2*power];

  for(int i = 0; i < Maxpow; i++) {
    paddedRoots[i] = Roots[i];
  }
  int dummyCount = 0;
  for(int i = Maxpow; i < 2*power; i++) {
    paddedRoots[i] = complex(-HUND, ZERO);
    dummyCount++;
  }
  int dummyCountHalf = dummyCount / 2;
  
  NaiveOrder(paddedRoots, 2*power);
  
  for(int i = 0; i < power; i++) {
    reversedRoots[2*i] = paddedRoots[dummyCountHalf+bitReversalRepresentation(i, digits)];
    reversedRoots[2*i+1] = complex(-HUND, ZERO);
    paddedRoots[dummyCountHalf+bitReversalRepresentation(i, digits)] = complex(-HUND, ZERO);
  }
  
  for(int i = 0; i < power; i++) {
    double difX;
    double difY;
    double dif;
    double bestDif = 1E10;
    int bestInd = 0;
    int I;
    for (I=0; I<2*power; I++) {
      difX = double_approx(realpart(reversedRoots[2*i] - paddedRoots[I]));
      difY = double_approx(imagpart(reversedRoots[2*i] + paddedRoots[I]));      
      dif = difX*difX + difY*difY;
      if (dif < bestDif) {
        bestDif = dif;
	bestInd = I;
      }    
    }
    reversedRoots[2*i+1] = paddedRoots[bestInd];
    paddedRoots[bestInd] = complex(-HUND, ZERO);
  }
  
  for(int i = 0, j=0; i < 2*power; i++) {
    if((realpart(reversedRoots[i]) != -HUND)) {
      Roots[j] = reversedRoots[i];
//    cout << j << " " << double_approx(realpart(Roots[j])) << " " <<	double_approx(imagpart(Roots[j])) <<  endl; 
      j++;
    }
  }


  if(LogLevel>3) {
    printf("...READY\n");
  }
}



// This bring the complex list Roots of length Maxpow to
// bit reversal order.
// See hep-lat/9805026
void PolynomialApproximation::BitReversalOrder(cln::cl_N *Roots, const int Maxpow) {
  if(LogLevel>3) {
    printf("Bringing roots to Bit-reversed order...\n");
  }

  int digits = 2;
  int power = 2;
  while(power < Maxpow) {
    power=power*2;
    digits++;
  }

  //  cout << "digits = " << digits << " power " << power << " " << pow(2.,digits-1) << endl;
  
  cln::cl_N paddedRoots[power];
  cln::cl_N reversedRoots[power];

  for(int i = 0; i < Maxpow; i++) {
    paddedRoots[i] = Roots[i];
  }
  for(int i = Maxpow; i < power; i++) {
    paddedRoots[i] = complex(-HUND, ZERO);
  }
  NaiveOrder(paddedRoots, power);
  for(int i = 0; i < power; i++) {
    reversedRoots[i] = paddedRoots[bitReversalRepresentation(i, digits)];
    // cout << i << " " << bitReversalRepresentation(i, digits) << endl;
  }
  
  for(int i = 0, j=0; i < power; i++) {
    if((realpart(reversedRoots[i]) != -HUND)) {
      Roots[j] = reversedRoots[i];
//    cout << j << " " << double_approx(realpart(Roots[j])) << " " <<	double_approx(imagpart(Roots[j])) <<  endl; 
      j++;
    }
  }
  if(LogLevel>3) {
    printf("...READY\n");
  }
}



/******************************************************************************/
//
//  Calculating the roots of a polynomial approximation minimizing
//  the integral of relative quadratic deviation from x^(-Alpha)
//  in an interval.
//
//  The coefficients of the polynomials are assumed to be known.

//  Input parameters:
//  order of the polynomial:                                Maxpow
//  lower bound of interval of approximation:               Epsilon
//  upper bound of interval of approximation:               Lambda
//
//  The name of a array containing the coefficients:        Coef
//                                                          Roots
//                                                          RootRoots
//                                                          RootRepTotalNorm
//                                                          RootRepSingleNorm
//
//  name of the file for writing results:                   Filename1
//  name of the file for writing results:                   Filename2
//  start file for Start=yes, otherwise append              Start
void PolynomialApproximation::ApproxiRootr(int Maxpow, cln::cl_F Epsilon, cln::cl_F Lambda, cln::cl_F* Coef, cln::cl_N* Roots, cln::cl_N* RootRoots, cln::cl_F& RootRepTotalNorm, cln::cl_F& RootRepSingleNorm, char* Filename1,char* Filename2) {
  int  ord, leng;
  
  if (LogLevel>2) {
    printf("Determining roots of approximation polynom...");
  }

  cln::cl_N Poly[1+Maxpow];
  cln::cl_F pi2 = cln::cl_float(realpart(acos(ZERO)))/HALF/HALF;
  cln::cl_F rr, ang, coef;
  
  // Check input
  if(Maxpow < 0) { 
    printf("wrong order in ApproxiRootr: %d\n",Maxpow); 
  }
  
  // Define polynomial with coefficient array
  for(ord = 0; ord < Maxpow+1; ord++) {
    Poly[ord] = complex(Coef[ord],ZERO);
  }
  // Find roots of the polynomial

  Polyrootc(Poly,Maxpow,Roots);

  BitReversalOrderOfRootPairs(Roots, Maxpow);
  enforceHermiticity(Roots, Maxpow);
      
  // Compute monomials
  for(int i = 0; i < Maxpow; i++) {
    RootRoots[i] = sqrt(Roots[i]);
    if(!plusp(imagpart(RootRoots[i]))) {
      RootRoots[i] = -RootRoots[i];
    }
    RootRoots[2*Maxpow-i-1] = conjugate(RootRoots[i]);
  }

  cln::cl_N a = EvalPoly(Poly, Maxpow, Epsilon);
  cln::cl_N b = EvalPolyProd(Roots, Maxpow, Epsilon);
  RootRepTotalNorm = As(cln::cl_F)(realpart(a/b));
  RootRepSingleNorm = As(cln::cl_F)(realpart(exp(log(RootRepTotalNorm)/cln::cl_R(Maxpow))));  
  
  if(LogLevel>3) {
    cout << endl <<"Polynomial: "<< a << endl <<endl<< "Product:"<< b << endl;
    cout << endl << "Normierungsfaktor = "<< double_approx(RootRepTotalNorm) << " ==> local = " << double_approx(RootRepSingleNorm) << endl;
  }

  leng = Maxpow;
  WriteRoots(Filename1, Roots, Maxpow);
  WriteRoots(Filename2, RootRoots, 2*Maxpow);
  WriteCoefs(Coef, Maxpow+1);
  if (LogLevel>2) {
    printf("Roots READY\n\n");  
  }
}



/******************************************************************************/
/*template <class T> cln::cl_F PolynomialApproximation::func(T &x) {
  return(ONE/sqrt(x));
}



// This routine produces coefficients for Chebycheff pol.
// of a given order and interval [epsilon, lambda]
void PolynomialApproximation::ChebyCoeff(const int order, const cln::cl_F &epsilon, const cln::cl_F &lambda, cln::cl_F * Coeff, cln::cl_F * c) {

  cln::cl_F bma = HALF*(lambda-epsilon);
  cln::cl_F bpa = HALF*(lambda+epsilon);
  cln::cl_F y;
  cln::cl_F ftable[5000];

  for(int i = 0; i < order+1; i++) {
    y = cos(cln::pi(bma)*As(cln::cl_F)((cln::cl_R(i)+HALF)/cln::cl_R(order+1)));
    ftable[i] = func(y*bma+bpa);
  }

  cln::cl_F fac = As(cln::cl_F)(TWO/cln::cl_R(order+1));
  for(int i = 0; i < order+1; i++) {
    cln::cl_F sumit = ZERO;
    for(int j = 0; j < order+1; j++) {
      sumit = sumit + ftable[j]*cos(cln::pi(bma)*(cln::cl_R(i)*(cln::cl_R(j)+HALF)/cln::cl_R(order+1)));
    }
    Coeff[i] = fac*sumit;
  }

  cln::cl_univpoly_real_ring PR = find_univpoly_ring(cln::cl_R_ring);
  cln::cl_UP_R b = PR->create(order);
  for(int i = 0; i < order+1; i++) {
    b.set_coeff(i, ZERO);
  }
  b.set_coeff(0, -HALF*Coeff[0]);
  b.finalize();

  for(int i=0; i < order+1; i++) {
    cln::cl_UP_I C = cln::tschebychev(i);
    for(int j=0; j < i+1; j++) {
      cln::cl_R c = coeff(b, j);
      b.set_coeff(j, c + Coeff[i]*cln::cl_R(coeff(C, j)));
    }
    b.finalize();
  }
  for(int i = 0; i < order+1; i++) {
    c[i] = Coeff[i];
    Coeff[i] = As(cln::cl_F)(coeff(b, i));
    cout << i << " : " << c[i] << endl;
  }
}

// here we use clenshaw to evaluate the polynomial

cln::cl_N PolynomialApproximation::EvalCheby(const int order, cln::cl_N * Coeff, cln::cl_N &x, 
	       const cln::cl_N & epsilon, const cln::cl_N lambda) {
  cln::cl_N d=complex(ZERO,ZERO), dd=complex(ZERO,ZERO), sv, z, res;
  int j;

  z = (TWO*x - epsilon - lambda)/(lambda-epsilon);

  for(j=order; j>=1; j--) {
    sv = d;
    d = TWO*z*d - dd + Coeff[j];
    dd = sv;
  }

  res = z*d - dd + HALF*Coeff[0];

  return(res);
}*/


void PolynomialApproximation::plotApproxPolynomials() {
  if (LogLevel>2) {
    printf("Plotting Approximation sub-polynomials to disk...");
  }
  int I;
  for (I=0; I<1+2*subPolyCount; I++) {
    char* fileName = new char[600];
    snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/masterApproxPoly%dSubPoly%d.dat",DataBaseDirectory,ID,I);
  
    ofstream out(fileName, ios::out);
    out.precision(20);
    cln::cl_F x = EPSILON[I/2];
    int MAXPOW = polyDegree[I/2];
    while (x < LAMBDA[I/2]) {
      cln::cl_N y1 = RootRepTotalNorm[I]*EvalPolyProd(Roots[I], MAXPOW, x);
      cln::cl_N y2 = EvalPoly(Coef[I], MAXPOW, x);
      if (I%2 == 0) {
        cln::cl_N y3 = Recurev(MAXPOW,x,Recb[I/2],Recg[I/2],Orth[I/2],Coed[I/2]);
        out << double_approx(realpart(x)) <<" " << double_approx(realpart(y1)) << " " << double_approx(realpart(y2)) << " " << double_approx(realpart(y3))<< endl;      
      } else {
        out << double_approx(realpart(x)) <<" " << double_approx(realpart(y1)) << " " << double_approx(realpart(y2)) << endl;
      }
      x = (ONE+HALF*HALF*HALF*HALF)*x;
    }
    out.close();
    delete[] fileName;
  }
  if (LogLevel>2) {
    printf("READY\n");
  }
}


void PolynomialApproximation::calcApproxPolys() {
  int I;
  for (I=0; I<1+subPolyCount; I++) {
    if (LogLevel>2) printf("Calculating subpolynomial %d of polynomial %d in orthogonal representation.\n",I,ID);
    char* fileName = new char[600];
    snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/orthogonalDataPol%dSubPoly%dApprox.dat",DataBaseDirectory,ID,I);
    Quadropt(polyDegree[I], EPSILON[I], LAMBDA[I], Recb[I], Recg[I], Orth[I], Coed[I], fileName); 
    delete[] fileName;

    RestorePolyCoef(polyDegree[I],Recb[I],Recg[I],Orth[I],Coed[I],Coef[2*I]);
  }
  for (I=0; I<subPolyCount; I++) {
    int deg1 = polyDegree[I+0];
    int deg2 = polyDegree[I+1];
    int I2;
    for (I2=0; I2<=deg1; I2++) {
      Coef[1+2*I][I2] = Coef[0+2*I][I2];
    }
    for (I2=0; I2<=deg2; I2++) {
      Coef[1+2*I][I2] = Coef[1+2*I][I2] - Coef[2+2*I][I2];
    }
    Coef[1+2*I][0] = Coef[1+2*I][0] + ONE;
  }

  for (I=0; I<1+2*subPolyCount; I++) {
    if (LogLevel>2) printf("Bringing to monomial Representation of subpolynomial %d of polynomial %d.\n",I,ID);
    char* fileName1 = new char[600];
    char* fileName2 = new char[600];
    snprintf(fileName1,600,"%s/data/results/pHMC/miscellaneous/polyApproxPol%dSubPoly%dRoots.dat",DataBaseDirectory,ID,I);
    snprintf(fileName2,600,"%s/data/results/pHMC/miscellaneous/polyApproxPol%dSubPoly%dRootsOfRoots.dat",DataBaseDirectory,ID,I);

    ApproxiRootr(polyDegree[I/2],EPSILON[I/2],LAMBDA[I/2],Coef[I],Roots[I],RootRoots[I],RootRepTotalNorm[I],RootRepSingleNorm[I],fileName1,fileName2);
    
    delete[] fileName1;
    delete[] fileName2;
  }
}


Complex* PolynomialApproximation::getApproxPolyRoots(int polynomSlot) {
  if ((polynomSlot<0) || (polynomSlot>=1+2*subPolyCount)) {
    printf("ERROR: PolynomSlot %d does not exist!!!",polynomSlot);
    exit(0);
  }
  int deg = polyDegree[polynomSlot/2];
  Complex* doublePrecRoots = new Complex[deg];
  int I;
  for (I=0; I<deg; I++) {
    doublePrecRoots[I].x = double_approx(realpart(Roots[polynomSlot][I]));
    doublePrecRoots[I].y = double_approx(imagpart(Roots[polynomSlot][I])); 
  }
  
  return doublePrecRoots;
}


int PolynomialApproximation::getApproxPolyDegree(int polynomSlot) {
  if ((polynomSlot<0) || (polynomSlot>=1+2*subPolyCount)) {
    printf("ERROR: PolynomSlot %d does not exist!!!",polynomSlot);
    exit(0);
  }
  int deg = polyDegree[polynomSlot/2];
  return deg;
}


double PolynomialApproximation::getPartPolyNormalization(int polynomSlot, int partPolyDeg) {
  if ((polynomSlot<0) || (polynomSlot>=1+2*subPolyCount)) {
    printf("ERROR: PolynomSlot %d does not exist!!!",polynomSlot);
    exit(0);
  }
  int deg = polyDegree[polynomSlot/2];
  cln::cl_F RootRepMultiNorm = As(cln::cl_F)(realpart(exp(cln::cl_R(partPolyDeg)*log(RootRepTotalNorm[polynomSlot])/cln::cl_R(deg))));  
  double res = double_approx(RootRepMultiNorm);
  return res;
}



double PolynomialApproximation::evaluatePolynomial(double x, int polynomSlot) {
  if ((polynomSlot<0) || (polynomSlot>=1+2*subPolyCount)) {
    printf("ERROR in evaluatePolynomial: PolynomSlot %d does not exist!!!",polynomSlot);
    exit(0);
  }
  
  cln::cl_F clnx = cl_float(x,clnDIGIT);
  cln::cl_N clny = EvalPoly(Coef[polynomSlot], polyDegree[polynomSlot/2], x);

  return double_approx(realpart(clny));
}
