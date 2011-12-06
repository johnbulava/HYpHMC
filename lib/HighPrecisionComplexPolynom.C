#include "HighPrecisionComplexPolynom.h"

HighPrecisionComplexPolynom::HighPrecisionComplexPolynom(int len, int digit) { 
  ini(len, digit);
}


HighPrecisionComplexPolynom::HighPrecisionComplexPolynom(const HighPrecisionComplexPolynom& pol) { 
  ini(pol);
}


HighPrecisionComplexPolynom::~HighPrecisionComplexPolynom() { 
  desini();
}


void HighPrecisionComplexPolynom::ini(int len, int digit) {
  DIGIT = digit;
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
  delete[] xxxStr;
  
  length = len;
  coeff = new cln::cl_N[len];
  for (int I=0; I<len; I++) {
    coeff[I] = complex(cl_float(0,clnDIGIT), cl_float(0,clnDIGIT));
  }
}


void HighPrecisionComplexPolynom::ini(const HighPrecisionComplexPolynom& pol) {
  ini(pol.length, pol.DIGIT);
  for (int I=0; I<pol.length; I++) {
    cln::cl_N c = pol.coeff[I];
    setCoeff(I, c);
  }
}


void HighPrecisionComplexPolynom::desini() {
  delete[] coeff;
  coeff = NULL;
  length = 0;
}


cln::cl_N HighPrecisionComplexPolynom::EvalPoly(cln::cl_N* Poly, int Maxpow, cln::cl_N Valu) {
   int  pow;
   cln::cl_N xpow, sum;

   sum = As(cln::cl_N)(complex(ZERO,ZERO));
   xpow = As(cln::cl_N)(complex(ONE,ZERO));

    for(pow = 0; pow < Maxpow+1; pow++)
     { sum = sum+xpow*Poly[pow];
       xpow = xpow*Valu; }

   return sum;
}


cln::cl_N HighPrecisionComplexPolynom::evaluateAt(cln::cl_N z) {
  return EvalPoly(coeff, getOrder(), z);
}


Complex HighPrecisionComplexPolynom::evaluateAt(Complex z) {
  cln::cl_N precZ = complex(cl_float(z.x,clnDIGIT), cl_float(z.y,clnDIGIT));
  cln::cl_N precRes = evaluateAt(precZ);
  Complex res(double_approx(realpart(precRes)), double_approx(imagpart(precRes)));

  return res;
}


HighPrecisionComplexPolynom HighPrecisionComplexPolynom::calcDerivative() {
  int len = length-1;
  if (length < 1) len = 1;
  HighPrecisionComplexPolynom pol(len, DIGIT);

  char* xxxStr = new char[1000];
  for (int I=0; I<length-1; I++) {
    snprintf(xxxStr,1000,"%d.0e+0_%d",I+1, DIGIT);    
    cln::cl_F f = xxxStr;
    cln::cl_N fac = complex(f, ZERO);
    pol.setCoeff(I, fac*coeff[I+1]);
  }
  delete[] xxxStr;
  return pol;
}

  
int HighPrecisionComplexPolynom::getLength() {
  return length;
}


void HighPrecisionComplexPolynom::setCoeff(int p, Complex val) {
  if ((p<0) || (p>=length)) return;
  
  coeff[p] = complex(cl_float(val.x,clnDIGIT), cl_float(val.y,clnDIGIT));
}


void HighPrecisionComplexPolynom::addCoeff(int p, Complex val) {
  if ((p<0) || (p>=length)) return;
  
  coeff[p] = coeff[p] + complex(cl_float(val.x,clnDIGIT), cl_float(val.y,clnDIGIT));
}

 
void HighPrecisionComplexPolynom::setCoeff(int p, cln::cl_N val) {
  if ((p<0) || (p>=length)) return;
  
  coeff[p] = val;
}


void HighPrecisionComplexPolynom::addCoeff(int p,cln::cl_N val) {
  if ((p<0) || (p>=length)) return;
  
  coeff[p] = coeff[p] + val;
}


cln::cl_N HighPrecisionComplexPolynom::getCoeff(int p) {
  if ((p<0) || (p>=length)) return complex(cl_float(0,clnDIGIT), cl_float(0,clnDIGIT));
  return coeff[p];
}
 
void HighPrecisionComplexPolynom::print() {
  for (int I=0; I<length; I++) {
    Complex res(double_approx(realpart(coeff[I])), double_approx(imagpart(coeff[I])));
    res.print();
  }
} 


HighPrecisionComplexPolynom HighPrecisionComplexPolynom::operator + (HighPrecisionComplexPolynom pol) {
  HighPrecisionComplexPolynom res(length, DIGIT);
  for (int I=0; I<length; I++) {
    res.setCoeff(I, coeff[I] + pol.getCoeff(I));
  }
  return res;
} 


HighPrecisionComplexPolynom HighPrecisionComplexPolynom::operator - (HighPrecisionComplexPolynom pol) {
  HighPrecisionComplexPolynom res(length, DIGIT);
  for (int I=0; I<length; I++) {
    res.setCoeff(I, coeff[I] - pol.getCoeff(I));
  }
  return res;
}


HighPrecisionComplexPolynom HighPrecisionComplexPolynom::operator * (HighPrecisionComplexPolynom pol) {
  HighPrecisionComplexPolynom res(length+pol.getLength()-1, DIGIT);
  for (int I=0; I<length; I++) {
    for (int I2=0; I2<pol.getLength(); I2++) {
      res.addCoeff(I+I2, coeff[I] * pol.getCoeff(I2));
    }
  }
  return res;
}


HighPrecisionComplexPolynom& HighPrecisionComplexPolynom::operator = (HighPrecisionComplexPolynom pol) {
  desini();
  ini(pol);

  return *this;
}


int HighPrecisionComplexPolynom::getOrder() {
  for (int I=length-1; I>=0; I--) {
    if (coeff[I] != ZERO) return I;    
  }
  return 0;
}


Complex* HighPrecisionComplexPolynom::getRoots() {
  if (getOrder()==0) return NULL;
  cln::cl_N* precRoots = getPrecRoots();
  Complex* roots = new Complex[getOrder()];
  for (int I=0; I<getOrder(); I++) {
    roots[I] = Complex(double_approx(realpart(precRoots[I])), double_approx(imagpart(precRoots[I])));
  }
  
  delete[] precRoots;
  return roots;
}


cln::cl_N* HighPrecisionComplexPolynom::getPrecRoots() {
  if (getOrder()==0) return NULL;
  cln::cl_N* roots = new cln::cl_N[getOrder()];
  Polyrootc(coeff, getOrder(), roots);
 
  return roots;
}

/******************************************************************************/
//
//  Find a root of a complex polynomial by Laguerre iteration.
//
//  The polynomial is                                    Poly
//  The order is                                         Maxpow
//
//  The precision:                                       Digit
cln::cl_N HighPrecisionComplexPolynom::Lasolv(cln::cl_N* Poly, int Maxpow, cln::cl_N root, int itemax) {
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

        if(abs(val0) == ZERO || (abs(val0) < small) && abs(val0/val) > 0.7) {
            if(LogLevel>4) { 
	       printf("Laguerre iterations: %d\n", ite);
               printf("root = %f +i* %f\n", double_approx(realpart(root)), double_approx(imagpart(root)));
               printf("value at root: %f +i* %f\n", double_approx(realpart(val0)), double_approx(imagpart(val0))); 
	    }

           break; 
	} 
    }

    if(ite >= itemax) {
      printf("Laguerre iteration did not converge\n");
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
void HighPrecisionComplexPolynom::Polyrootc(cln::cl_N* Poly, int Maxpow, cln::cl_N* Root) {
  
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
	printf(" root number: %d\n",fnd+1);
      }

      if((ord%2 == 1) && (Maxpow%2 == 0) && (false)) {
	Root[fnd] = As(cln::cl_N)(conjugate(Root[fnd-1]));
 	cln::cl_N val0 = EvalPoly(poly, Maxpow, Root[fnd]);
	
 	if(LogLevel>3) {
 	  printf("root = %f +i*%f\n", double_approx(realpart(Root[fnd])), double_approx(imagpart(Root[fnd])));
 	  printf("value at root: %f +i*%f\n", double_approx(realpart(val0)), double_approx(imagpart(val0))); 
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
    printf("control polynomial should be close to zero:\n");

    for(pow = 0; pow < Maxpow+1; pow++) {
      printf("  x^{%d}\n",pow);
      printf("%1.15f +i*%1.15f\n",double_approx(realpart(polc[pow])), double_approx(imagpart(polc[pow])));
    } 
  }
}
