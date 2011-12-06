void AutoCorrelation::ini(int obsCount, int GWmax) {
  observableCount = obsCount;
  GammaWMax = GWmax;
  Gammas = new double**[obsCount];
  GammaSigmas = new double**[obsCount];
  int I,I2;
  for (I=0; I<obsCount; I++) {
    Gammas[I] = new double*[obsCount];
    GammaSigmas[I] = new double*[obsCount];
    for (I2=0; I2<obsCount; I2++) {
      Gammas[I][I2] = new double[GammaWMax];
      GammaSigmas[I][I2] = new double[GammaWMax];
    }
  }
  CombinedGamma = new double[GammaWMax];
  CombinedGammaSigma = new double[GammaWMax];
  WstepArray = new double[GammaWMax];
  CombinedCFunction  = new double[GammaWMax];
  CombinedReducedCFunction  = new double[GammaWMax];
  CombinedCFunctionErrors = new double[GammaWMax];
  CombinedReducedCFunctionErrors = new double[GammaWMax];
  averages = new double[observableCount];
  totalN = 0;
  AutoCorrelationTime = 0;
}


void AutoCorrelation::desini() {
  int I,I2;
  for (I=0; I<observableCount; I++) {
    for (I2=0; I2<observableCount; I2++) {
      delete[] Gammas[I][I2];
      delete[] GammaSigmas[I][I2];      
    }
    delete[] Gammas[I];
    delete[] GammaSigmas[I];    
  }  
  delete[] Gammas;
  delete[] GammaSigmas;
  delete[] CombinedGamma;
  delete[] CombinedGammaSigma;
  delete[] WstepArray;
  delete[] CombinedCFunction;
  delete[] CombinedReducedCFunction;
  delete[] CombinedCFunctionErrors;
  delete[] CombinedReducedCFunctionErrors;
  delete[] averages;
}


AutoCorrelation::AutoCorrelation() {
  ini(4,100);
}


AutoCorrelation::AutoCorrelation(int obsCount, int GWmax) {
  ini(obsCount, GWmax);
}


AutoCorrelation::~AutoCorrelation() {
  desini();
}


double AutoCorrelation::getAverage(int nr) {
  return averages[nr];
}


int AutoCorrelation::getGammaWMax() {
  return GammaWMax;
}


int AutoCorrelation::getTotalN() {
  return totalN;
}


double* AutoCorrelation::getCombinedGamma() {
  return CombinedGamma;
}

double* AutoCorrelation::getCombinedCFunction() {
  return CombinedCFunction;
}

double* AutoCorrelation::getCombinedReducedCFunction() {
  return CombinedReducedCFunction;
}


double* AutoCorrelation::getCombinedCFunctionErrors() {
  return CombinedCFunctionErrors;
}


double* AutoCorrelation::getCombinedReducedCFunctionErrors() {
  return CombinedReducedCFunctionErrors;
}


void AutoCorrelation::loadData(int RunCount, int* RunLengths, double* detData, double* measureData) {
  int I,I2, ob, ob1, ob2, t, R, count, index;
  totalN = 0;
  for (I=0; I<RunCount; I++) totalN += RunLengths[I];
  if (totalN<=0) {
    // Keine Daten !!!
    totalN = 0;
    GammaWMax = 0;
    return;
  }
  for (I=0; I<observableCount; I++) {
    averages[I] = 0;
  }
  double* dummy1 = new double[observableCount];
  double* dummy2 = new double[observableCount];
  for (I2=0; I2<totalN; I2++) {
    dummy1[0] = detData[I2];
    averages[0] += detData[I2];
    for (I=1; I<observableCount; I++) {
      dummy1[I] = measureData[I2] * dummy1[I-1];
      averages[I] += dummy1[I];
    }
  }
  double avgWeight = averages[0] / totalN;
  for (I=0; I<observableCount; I++) {
    averages[I] /= (totalN * avgWeight);
  }
  
  for (t=GammaWMax-1; t>=0; t--) {
    for (ob1=0; ob1<observableCount; ob1++) {
      for (ob2=ob1; ob2<observableCount; ob2++) {
        Gammas[ob1][ob2][t] = 0;
        GammaSigmas[ob1][ob2][t] = 0;
      }
    }
	
    count = 0;
    index = 0;
    for (R=0; R<RunCount; R++) {
      for (I=0; I<RunLengths[R]-t; I++) {
        dummy1[0] = detData[I+index] / avgWeight;
        dummy2[0] = detData[I+index+t] / avgWeight;
        for (ob=1; ob<observableCount; ob++) {
	  dummy1[ob] = measureData[I+index] * dummy1[ob-1];
	  dummy2[ob] = measureData[I+index+t] * dummy2[ob-1];
	}

        for (ob1=0; ob1<observableCount; ob1++) {
          for (ob2=ob1; ob2<observableCount; ob2++) {
            Gammas[ob1][ob2][t] += (dummy1[ob1]-averages[ob1]) * (dummy2[ob2]-averages[ob2]);
            GammaSigmas[ob1][ob2][t] += sqr((dummy1[ob1]-averages[ob1]) * (dummy2[ob2]-averages[ob2]));
	  }
	}
        count++;
      }
      index += RunLengths[R];
    }
    for (ob1=0; ob1<observableCount; ob1++) {
      for (ob2=ob1; ob2<observableCount; ob2++) {
        if (count>0) {
          Gammas[ob1][ob2][t] /= count;
          GammaSigmas[ob1][ob2][t] = sqrt((GammaSigmas[ob1][ob2][t]/count) - sqr(Gammas[ob1][ob2][t]))/sqrt(count);
          Gammas[ob2][ob1][t] = Gammas[ob1][ob2][t];
          GammaSigmas[ob2][ob1][t] = GammaSigmas[ob1][ob2][t];
	} else {
	  GammaWMax = t;
	}
      }
    }
  }
  if (totalN<GammaWMax-10) GammaWMax = (int) (totalN-10);
  if (GammaWMax<=0) GammaWMax = 0;
  
  delete[] dummy1;
  delete[] dummy2;
}


void AutoCorrelation::calcCombinedGammaFunction(ComplexVector& derivatives) {
  int t, I, I2;
  for (t=0; t<GammaWMax; t++) {
    CombinedGamma[t] = 0;
    CombinedGammaSigma[t] = 0;
    for (I=0; I<observableCount; I++) {
      for (I2=0; I2<observableCount; I2++) {
        CombinedGamma[t] += Gammas[I][I2][t] * derivatives.vectorElements[I].x * derivatives.vectorElements[I2].x;
        CombinedGammaSigma[t] += sqr(GammaSigmas[I][I2][t]) * derivatives.vectorElements[I].x * derivatives.vectorElements[I2].x;	
      }
    }
    CombinedGammaSigma[t] = sqrt(CombinedGammaSigma[t]); 
  }
}


void AutoCorrelation::calcCombinedCFunction(ComplexVector& derivatives) {
  calcCombinedGammaFunction(derivatives);
  int W, I;    
  for (W=0; W<GammaWMax; W++) {
    CombinedCFunction[W] = CombinedGamma[0];
    CombinedCFunctionErrors[W] = sqr(CombinedGammaSigma[0]);
    CombinedReducedCFunction[0] = 1; 
    CombinedReducedCFunctionErrors[0] = 0;
    for (I=1; I<=W; I++) {    
      CombinedCFunction[W] += 2*CombinedGamma[I];
      CombinedCFunctionErrors[W] += 2*CombinedGammaSigma[I];
      CombinedReducedCFunction[W] = CombinedCFunction[W] / CombinedCFunction[0];      
      CombinedReducedCFunctionErrors[W] = CombinedCFunctionErrors[W] / CombinedCFunction[0];      
    }
  }
}


int findInsetIndexOfPlateau(double* data, int count, double slopeFac) {
  if (count <= 0) return -1;
  if (data == NULL) return -1;
  
  double smallestSlope = 0;
  if  (count>1) {
    smallestSlope = data[1] - data[0];
  }
  int p = 0;
  for (int I=0; I<count-1; I++) {
    if (data[I+1] <= data[I]) break;
    double s = 0;
    s = data[I+1] - data[I];
    if (s>slopeFac*smallestSlope) break;
    if (s<smallestSlope) smallestSlope = s;
    p++;  
  }  
  
  return p;
}


double AutoCorrelation::estimateCombinedError(ComplexVector& derivatives) {
  if (GammaWMax==0) {
    return NaN;
  }
  calcCombinedCFunction(derivatives);
  
  double C = 0;
  int pInd = findInsetIndexOfPlateau(CombinedCFunction, GammaWMax, 5.0);
  if (pInd>=0) {
    C = CombinedCFunction[pInd];
  }
  
  return sqrt(C / totalN);
}
  
  
double AutoCorrelation::estimateAutoCorrelationTime() {
  ComplexVector derivatives(5);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;

  return estimateAutoCorrelationTime(derivatives);
}

  
double AutoCorrelation::estimateAutoCorrelationTime(ComplexVector& derivatives) {
  double C = sqr(estimateCombinedError(derivatives))*totalN;
  double Gamma0 = CombinedCFunction[0];
  double aTime = 1/log(1+2*Gamma0/(C-Gamma0));
  return aTime;
}
  

double* AutoCorrelation::getWstepArray() {
  int I;
  for (I=0; I<GammaWMax; I++) {
    WstepArray[I] = I;
  }
  return WstepArray;
}
