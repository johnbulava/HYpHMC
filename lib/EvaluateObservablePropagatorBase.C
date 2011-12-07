#include "EvaluateObservablePropagatorBase.h"

EvaluateObservablePropagatorBase::EvaluateObservablePropagatorBase(AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, EvaluateObservable* obsWeight, EvaluateObservable* obsDetSign, char* oName, char* nick, double relStart, double relEnd) : EvaluateObservable(aIOcon, sdr, oName, nick, relStart, relEnd) { 
  L0 = SDReader->getL0();  
  L1 = SDReader->getL1();  
  L2 = SDReader->getL2();  
  L3 = SDReader->getL3();  
  latticeBins = new LatticeMomentumBins(L0,L1,L2,L3);  
  autoCorr = new AutoCorrelation(2, 100);
  avgProp = new double[latticeBins->getMomentumSqrSlotCount()];
  sigmaProp = new double[latticeBins->getMomentumSqrSlotCount()];
  autoCorrelationTime = new double[latticeBins->getMomentumSqrSlotCount()];
  pSqr = new double[latticeBins->getMomentumSqrSlotCount()]; 
  kappaZeroMode = false;
  if (abs(SDReader->getKappa()) < 1E-8) kappaZeroMode = true;
  
  PropagatorEuclideanZFactor = NaN;
  PropagatorEuclideanZFactorError = NaN;
  PropagatorEuclideanMass = NaN;
  PropagatorEuclideanMassError = NaN;
  PropagatorMinkowskiZFactor = NaN;
  PropagatorMinkowskiZFactorError = NaN;
  PropagatorMinkowskiMass = NaN;
  PropagatorMinkowskiMassError = NaN;  
  PropagatorPoleMass = NaN;
  PropagatorPoleMassError = NaN;  
  PropagatorPoleDecayWidth = NaN;
  PropagatorPoleDecayWidthError = NaN;  
  PropagatorPoleValue = NaN;
  PropagatorPoleValueError = NaN;
  PropagatorFitConst0 = NaN;
  PropagatorFitConst0Error = NaN;
  PropagatorFitConst1 = NaN;
  PropagatorFitConst1Error = NaN;
  PropagatorFitConst2 = NaN;
  PropagatorFitConst2Error = NaN;
  
  PropagatorEuclideanZFactorReduced = NaN;
  PropagatorEuclideanZFactorErrorReduced = NaN;
  PropagatorEuclideanMassReduced = NaN;
  PropagatorEuclideanMassErrorReduced = NaN;
  PropagatorMinkowskiZFactorReduced = NaN;
  PropagatorMinkowskiZFactorErrorReduced = NaN;
  PropagatorMinkowskiMassReduced = NaN;
  PropagatorMinkowskiMassErrorReduced = NaN;  
  PropagatorFitConst0Reduced = NaN;
  PropagatorFitConst0ErrorReduced = NaN;
  PropagatorFitConst1Reduced = NaN;
  PropagatorFitConst1ErrorReduced = NaN;
  PropagatorFitConst2Reduced = NaN;
  PropagatorFitConst2ErrorReduced = NaN;

  PropagatorFitConst0ReducedPole = NaN;
  PropagatorFitConst0ErrorReducedPole = NaN;
  PropagatorFitConst1ReducedPole = NaN;
  PropagatorFitConst1ErrorReducedPole = NaN;
  PropagatorFitConst2ReducedPole = NaN;
  PropagatorFitConst2ErrorReducedPole = NaN;
  PropagatorFitConst3ReducedPole = NaN;
  PropagatorFitConst3ErrorReducedPole = NaN;

  PropagatorZFactor = NaN;
  PropagatorZFactorError = NaN;  
  PropagatorZ0Factor = NaN;
  PropagatorZ0FactorError = NaN;  
  PropagatorMProp = NaN;
  PropagatorMPropError = NaN;  
  PropagatorMProp0 = NaN;
  PropagatorMProp0Error = NaN;  
  
  
  selfEnergySubLamFac = NaN;
  doArcTanhFit = false;
  gamma = NaN;
  vevContNot = NaN;
  HiggsPropFourParameterFit = false;
  
  considerZeroMomentum = true;
  
  PropagatorFitReducedChiSquare = NaN;
  PropagatorFitReducedChiSquareReduced = NaN;
  PropagatorFitReducedChiSquarePole = NaN;
    
  ini(getAnalyzerResultsCount(), obsWeight, obsDetSign);
}


EvaluateObservablePropagatorBase::~EvaluateObservablePropagatorBase() {
  delete latticeBins;
  delete autoCorr;
  delete[] avgProp;  
  delete[] sigmaProp;
  delete[] autoCorrelationTime;
  delete[] pSqr;
}


void EvaluateObservablePropagatorBase::defineObsDependencies() { 
  addDependOnObsByName("Magnetizations");
  if (!isObs("GoldstonePropagator")) {
    addDependOnObsByName("GoldstonePropagator");
  }
}


int EvaluateObservablePropagatorBase::getAnalyzerResultsCount() {
  return latticeBins->getMomentumSqrSlotCount();
}


void EvaluateObservablePropagatorBase::calcPropAtP0(int ignoreStart, int ignoreEnd, double &avgProp0) {
  ComplexVector derivatives(2);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;

  double* data = new double[dataAvailCount];
  double* dataDet = new double[dataAvailCount];  
  
  int dataStartOffsetForMag = getAnalyzerResultsCount();  

  int count = 0;
  for (int I2=0; I2<dataAvailCount; I2++) {
    if ((I2<ignoreStart) || (I2>ignoreEnd)) {
      data[count] = dataAvail[I2][dataStartOffsetForMag+0];
      dataDet[count] = dataAvailWeightAndSign[I2];	
      count++;
    }
  }
  autoCorr->loadData(1, &(count), dataDet, data);
  double vevMode = autoCorr->getAverage(1);
  double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);
  rescale *= sqrt(L0*L1*L2*L3);
  vevMode *= rescale;
  double vevModeSqr = vevMode * vevMode;
  
  count = 0;
  for (int I2=0; I2<dataAvailCount; I2++) {
    if ((I2<ignoreStart) || (I2>ignoreEnd)) {
      data[count] = dataAvail[I2][0];
      dataDet[count] = dataAvailWeightAndSign[I2];	
      count++;
    }
  }
  autoCorr->loadData(1, &(count), dataDet, data);
  avgProp0 = autoCorr->getAverage(1);
  avgProp0 -= vevModeSqr;

  delete[] data;
  delete[] dataDet;
}


Complex calcGoldstone1LoopInvPropagatorHighPrecision(Complex p, double m0, double mg, double mh, double Z, double coeff) {  
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
}


void EvaluateObservablePropagatorBase::calcProp(int ignoreStart, int ignoreEnd, double &EucZ, double &Eucmass, double &MinkZ, double &Minkmass, double &fC0, double &fC1, double &fC2,  double &EucZReduced, double &EucmassReduced, double &MinkZReduced, double &MinkmassReduced, double &fC0Reduced, double &fC1Reduced, double &fC2Reduced, double &Z, double &Z0, double &mprop, double &mprop0, Complex& pole, Complex& poleVal, double &fC0P, double &fC1P, double &fC2P, double &fC3P, double &redChiSqr, double &redChiSqrReduced, double &redChiSqrPole) {
  ComplexVector derivatives(2);
  derivatives.setZero();
  derivatives.vectorElements[1].x = 1;

  double* data = new double[dataAvailCount];
  double* dataDet = new double[dataAvailCount];  
  double* y = new double[latticeBins->getMomentumSqrSlotCount()];
  double* yerr = new double[latticeBins->getMomentumSqrSlotCount()];
  
  int dataStartOffsetForMag = getAnalyzerResultsCount();  

  bool isGoldStone = false;
  if (isObs("GoldstonePropagator")) {
    isGoldStone = true;
  }

  int count = 0;
  for (int I2=0; I2<dataAvailCount; I2++) {
    if ((I2<ignoreStart) || (I2>ignoreEnd)) {
      data[count] = dataAvail[I2][dataStartOffsetForMag+0];
      dataDet[count] = dataAvailWeightAndSign[I2];	
      count++;
    }
  }
  autoCorr->loadData(1, &(count), dataDet, data);
  double vevMode = autoCorr->getAverage(1);
  double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);
  vevContNot = vevMode * rescale;
  rescale *= sqrt(L0*L1*L2*L3);
  vevMode *= rescale;
  double vevModeSqr = vevMode * vevMode;
  
  for (int I=0; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    int count = 0;
    for (int I2=0; I2<dataAvailCount; I2++) {
      if ((I2<ignoreStart) || (I2>ignoreEnd)) {
        data[count] = dataAvail[I2][I];
        dataDet[count] = dataAvailWeightAndSign[I2];	
	count++;
      }
    }
    autoCorr->loadData(1, &(count), dataDet, data);
    pSqr[I] = latticeBins->getLatMomSqrFromSlotNr(I);
    avgProp[I] = autoCorr->getAverage(1);
    autoCorrelationTime[I] = autoCorr->estimateAutoCorrelationTime();
    
    if ((I==0) && (!isGoldStone)) {
      avgProp[I] -= vevModeSqr;
    } else {
      sigmaProp[I] = autoCorr->estimateCombinedError(derivatives);    
    }
    
    y[I] = 1.0 / avgProp[I];
    yerr[I] = abs(sigmaProp[I] / sqr(avgProp[I]));
  }
  delete[] data;
  delete[] dataDet;
  
 
  char* functionBody = new char[1000];
  double* fitRes = new double[6];
  double* fitErr = new double[6];
  
  int startP = 1;
//  if (considerZeroMomentum) startP = 0;
  snprintf(functionBody,1000,"(x+A2*A2)/A1");
  fitRes[0] = 1.0;
  fitRes[1] = 1.0;
  fitRes[2] = 0.0;
  fitRes[3] = 0.0;
  fitRes[4] = 0.0;
  fitRes[5] = 0.0;
  bool b = performGnuplotFit(functionBody, &(pSqr[startP]), &(y[startP]), &(yerr[startP]), latticeBins->getMomentumSqrSlotCount()-startP, 2, fitRes, fitErr, redChiSqr);
  
  snprintf(functionBody,1000,"(x+A2*A2+A3*x/(x+A4))/A1");
  fitRes[2] = 1.0;
  fitRes[3] = 1.0;
  b = b & performGnuplotFit(functionBody, &(pSqr[startP]), &(y[startP]), &(yerr[startP]), latticeBins->getMomentumSqrSlotCount()-startP, 4, fitRes, fitErr, redChiSqr);

  EucZ = fitRes[0];
  if (abs(fitRes[2]) > 1E-10) EucZ = fitRes[0]*fitRes[3] / (fitRes[2] + fitRes[3]);
  Eucmass = sqrt(EucZ) * abs(fitRes[1]) / sqrt(fitRes[0]);
  fC0 = fitRes[0];
  fC1 = fitRes[2];
  fC2 = fitRes[3];

  //Find Minkoski - Mass and Z-Factor
  double bound1=-2;
  double bound2=0;
  while (bound2-bound1>1E-10) {
    double middle = 0.5*(bound1+bound2);
    double val = middle + sqr(fitRes[1]) + fitRes[2]*middle/(middle+fitRes[3]);
    if (val>0) {
      bound2 = middle;
    } else {
      bound1 = middle;
    }
  }
  Minkmass = sqrt(-0.5*(bound1+bound2));
  MinkZ = fitRes[0] / (1 + fitRes[2]*fitRes[3]/sqr(fitRes[3]-Minkmass*Minkmass));

  if ((!b) || (isNaN(fitRes[0])) || (isNaN(fitRes[1])) || (isNaN(fitRes[2])) || (isNaN(fitRes[3]))) {
    EucZ = NaN;
    Eucmass = NaN;
    MinkZ = NaN;
    Minkmass = NaN;    
    fC0 = NaN;
    fC1 = NaN;
    fC2 = NaN;
  }
  
  
  //Calculate Reduced data set
  int reducedDataCount = 0;
  int loopStartI = 0;
  if (isObs("GoldstonePropagator")) {
    loopStartI = 1;
  }
  for (int I=loopStartI; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    if (pSqr[I]<=gamma) reducedDataCount++;
  }
  double* reducedDataPSqr = new double[reducedDataCount];
  double* reducedDataY = new double[reducedDataCount];
  double* reducedDataYerr = new double[reducedDataCount];
  reducedDataCount = 0;
  for (int I=loopStartI; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    if (pSqr[I]<=gamma) {
      reducedDataPSqr[reducedDataCount] = pSqr[I];
      reducedDataY[reducedDataCount] = y[I];
      reducedDataYerr[reducedDataCount] = yerr[I];      
      reducedDataCount++;
    }
  }
  
  
  snprintf(functionBody,1000,"(x+A2*A2)/A1");
  fitRes[0] = 1.0;
  fitRes[1] = 1.0;
  fitRes[2] = 0.0;
  fitRes[3] = 0.0;
  fitRes[4] = 0.0;
  fitRes[5] = 0.0;
  //Nur Linearer Fit
  b = performGnuplotFit(functionBody, reducedDataPSqr, reducedDataY, reducedDataYerr, reducedDataCount, 2, fitRes, fitErr, redChiSqrReduced);
  
//  snprintf(functionBody,1000,"(x+A2*A2+A3*x/(x+A4))/A1");
//  fitRes[2] = 1.0;
//  fitRes[3] = 1.0;
//  b = b & performGnuplotFit(functionBody, reducedDataPSqr, reducedDataY, reducedDataYerr, reducedDataCount, 4, fitRes, fitErr, redChiSqrReduced);

  EucZReduced = fitRes[0];
  if (abs(fitRes[2]) > 1E-10) EucZReduced = fitRes[0]*fitRes[3] / (fitRes[2] + fitRes[3]);
  EucmassReduced = sqrt(EucZReduced) * abs(fitRes[1]) / sqrt(fitRes[0]);
  fC0Reduced = fitRes[0];
  fC1Reduced = fitRes[2];
  fC2Reduced = fitRes[3];

  //Find Minkoski - Mass and Z-Factor
  bound1=-2;
  bound2=0;
  while (bound2-bound1>1E-10) {
    double middle = 0.5*(bound1+bound2);
    double val = middle + sqr(fitRes[1]) + fitRes[2]*middle/(middle+fitRes[3]);
    if (val>0) {
      bound2 = middle;
    } else {
      bound1 = middle;
    }
  }
  MinkmassReduced = sqrt(-0.5*(bound1+bound2));
  MinkZReduced = fitRes[0] / (1 + fitRes[2]*fitRes[3]/sqr(fitRes[3]-MinkmassReduced*MinkmassReduced));

  if ((!b) || (isNaN(fitRes[0])) || (isNaN(fitRes[1])) || (isNaN(fitRes[2])) || (isNaN(fitRes[3]))) {
    EucZReduced = NaN;
    EucmassReduced = NaN;
    MinkZReduced = NaN;
    MinkmassReduced = NaN;    
    fC0Reduced = NaN;
    fC1Reduced = NaN;
    fC2Reduced = NaN;
  }

  //Do double ArcTanh-Continuation
  if (doArcTanhFit) {
    if (isObs("GoldstonePropagator")) {
      //Goldstone-Propagator Fit
      pole = Complex(NaN,NaN);
      poleVal = Complex(NaN,NaN);
      double eps = 1E-10;
      
      double fitFac = 1;
      for (int fitCounter=0; fitCounter<5; fitCounter++) {
        fitRes[0] = 1;
        fitRes[1] = 1;
        fitRes[2] = 1E-6;
        fitRes[3] = 0.5 / fitFac;            
	fitFac *= 2;
      
        if (kappaZeroMode) {
          snprintf(functionBody,1000,"0.1*(x+A2*A2) / A1");	
	} else {
          snprintf(functionBody,1000,"(x+A2*A2) / A1");
	}
        b = performGnuplotFit(functionBody, reducedDataPSqr, reducedDataY, reducedDataYerr, reducedDataCount, 2, fitRes, fitErr, redChiSqrPole);
      
        if (kappaZeroMode) {
          snprintf(functionBody,1000,"0.1*(x+A2*A2 + A3*(( log(A4*A4/(A2*A2+%1.15f))*(A2*A2+%1.15f-A4*A4)  \
                                                        + sqrt(((A2*A2+%1.15f-A4*A4) + x)**2 + 4*A4*A4*x) \
	    				  	        * log( ((sqrt(((A2*A2+%1.15f-A4*A4) + x)**2 + 4*A4*A4*x) + x)**2 - (A2*A2+%1.15f-A4*A4)**2) \
						             / ((sqrt(((A2*A2+%1.15f-A4*A4) + x)**2 + 4*A4*A4*x) - x)**2 - (A2*A2+%1.15f-A4*A4)**2) ) ) / x \
						        - (2 + log((A2*A2+%1.15f)/(A4*A4)) * (1+2*A4*A4/(A2*A2+%1.15f-A4*A4))) )   ) / A1",
						        eps,eps,eps,eps,eps,eps,eps,eps,eps);
        } else {
          snprintf(functionBody,1000,"(x+A2*A2 + A3*(( log(A4*A4/(A2*A2+%1.15f))*(A2*A2+%1.15f-A4*A4)  \
                                                    + sqrt(((A2*A2+%1.15f-A4*A4) + x)**2 + 4*A4*A4*x) \
	    				  	    * log( ((sqrt(((A2*A2+%1.15f-A4*A4) + x)**2 + 4*A4*A4*x) + x)**2 - (A2*A2+%1.15f-A4*A4)**2) \
						         / ((sqrt(((A2*A2+%1.15f-A4*A4) + x)**2 + 4*A4*A4*x) - x)**2 - (A2*A2+%1.15f-A4*A4)**2) ) ) / x \
						    - (2 + log((A2*A2+%1.15f)/(A4*A4)) * (1+2*A4*A4/(A2*A2+%1.15f-A4*A4))) )   ) / A1",
						    eps,eps,eps,eps,eps,eps,eps,eps,eps);
	
	
	}
      
//        b = b & performGnuplotFit(functionBody, reducedDataPSqr, reducedDataY, reducedDataYerr, reducedDataCount, 3, fitRes, fitErr, redChiSqrPole);
        b = b & performGnuplotFit(functionBody, reducedDataPSqr, reducedDataY, reducedDataYerr, reducedDataCount, 4, fitRes, fitErr, redChiSqrPole);

        if (kappaZeroMode) {
          fitRes[0] *= 10;
	}

        if (b) break;
      }

      fC0P = fitRes[0];
      fC1P = fitRes[1];
      fC2P = fitRes[2];
      fC3P = fitRes[3];
      
      double h = 1E-2;
      Z0 = calcGoldstone1LoopInvPropagatorWithP0PartSubtractedHighPrecision(Complex(0,0), fitRes[1], fitRes[1], fitRes[3], fitRes[0], fitRes[2]).x;
      Z0 -= calcGoldstone1LoopInvPropagatorWithP0PartSubtractedHighPrecision(Complex(0,h), fitRes[1], fitRes[1], fitRes[3], fitRes[0], fitRes[2]).x;
      Z0 /= h*h;
      Z0 = 1.0 / Z0;
      
      mprop0 = sqrt(Z0) * abs(fitRes[1]) / sqrt(fitRes[0]);

      bound1=0;
      bound2=2;
      while (bound2-bound1>1E-10) {
        double middle = 0.5*(bound1+bound2);
        double val = calcGoldstone1LoopInvPropagatorWithP0PartSubtractedHighPrecision(Complex(0,middle), fitRes[1], fitRes[1], fitRes[3], fitRes[0], fitRes[2]).x;
        if (val>0) {
          bound1 = middle;
        } else {
          bound2 = middle;
        }
      }
      mprop = 0.5*(bound1+bound2);

      Z = calcGoldstone1LoopInvPropagatorWithP0PartSubtractedHighPrecision(Complex(0,mprop), fitRes[1], fitRes[1], fitRes[3], fitRes[0], fitRes[2]).x;
      Z -= calcGoldstone1LoopInvPropagatorWithP0PartSubtractedHighPrecision(Complex(0,mprop+h), fitRes[1], fitRes[1], fitRes[3], fitRes[0], fitRes[2]).x;
      Z /= (sqr(mprop+h) - sqr(mprop));
      Z = 1.0 / Z;
    } else {  
      //Higgs-Propagator Fit
      double mG = NaN;
      double eps = 1E-10;
      mG = ((EvaluateObservablePropagatorBase*)(getDependOnObsByName("GoldstonePropagator")))->getPropagatorMass();

      double fitFac = 1;
      for (int fitCounter=0; fitCounter<5; fitCounter++) {
        fitRes[0] = 1;
        fitRes[1] = 1/fitFac;
        fitRes[2] = 1E-6/fitFac;
        fitRes[3] = mG;
	
	fitFac *= 2;
      
        if (kappaZeroMode) {
          snprintf(functionBody,1000,"0.1*(x+A2*A2) / A1");
	} else {
          snprintf(functionBody,1000,"(x+A2*A2) / A1");	
	}
        b = performGnuplotFit(functionBody, reducedDataPSqr, reducedDataY, reducedDataYerr, reducedDataCount, 2, fitRes, fitErr, redChiSqrPole);
  
        if (!HiggsPropFourParameterFit) {
          if (kappaZeroMode) {
            snprintf(functionBody,1000,"0.1*(x+A2*A2 + A3*(288*(sqrt((x+4*A2*A2+%1.15f)/(x+%1.15f))*atanh(sqrt((x+%1.15f)/(x+4*A2*A2+%1.15f))) -1) \
                                                          + 96*(sqrt((x+4*%1.15f*%1.15f+%1.15f)/(x+%1.15f))*atanh(sqrt((x+%1.15f)/(x+4*%1.15f*%1.15f+%1.15f))) -1)) )/A1",
                                                          eps,eps,eps,eps,mG,mG,eps,eps,eps,mG,mG,eps);
          } else {
            snprintf(functionBody,1000,"(x+A2*A2 + A3*(288*(sqrt((x+4*A2*A2+%1.15f)/(x+%1.15f))*atanh(sqrt((x+%1.15f)/(x+4*A2*A2+%1.15f))) -1) \
                                                      + 96*(sqrt((x+4*%1.15f*%1.15f+%1.15f)/(x+%1.15f))*atanh(sqrt((x+%1.15f)/(x+4*%1.15f*%1.15f+%1.15f))) -1)) )/A1",
                                                      eps,eps,eps,eps,mG,mG,eps,eps,eps,mG,mG,eps);
	  }

          b = b & performGnuplotFit(functionBody, reducedDataPSqr, reducedDataY, reducedDataYerr, reducedDataCount, 3, fitRes, fitErr, redChiSqrPole);

          if (kappaZeroMode) {
            fitRes[0] *= 10;
	  }
	} else {
          if (kappaZeroMode) {
            snprintf(functionBody,1000,"0.1*(x+A2*A2 + A3*(288*(sqrt((x+4*A2*A2+%1.15f)/(x+%1.15f))*atanh(sqrt((x+%1.15f)/(x+4*A2*A2+%1.15f))) -1) \
                                                          + 96*(sqrt((x+4*A4*A4+%1.15f)/(x+%1.15f))*atanh(sqrt((x+%1.15f)/(x+4*A4*A4+%1.15f))) -1)) )/A1",
                                                          eps,eps,eps,eps,eps,eps,eps,eps);
          } else { 
            snprintf(functionBody,1000,"(x+A2*A2 + A3*(288*(sqrt((x+4*A2*A2+%1.15f)/(x+%1.15f))*atanh(sqrt((x+%1.15f)/(x+4*A2*A2+%1.15f))) -1) \
                                                      + 96*(sqrt((x+4*A4*A4+%1.15f)/(x+%1.15f))*atanh(sqrt((x+%1.15f)/(x+4*A4*A4+%1.15f))) -1)) )/A1",
                                                      eps,eps,eps,eps,eps,eps,eps,eps);
	  }
          b = b & performGnuplotFit(functionBody, reducedDataPSqr, reducedDataY, reducedDataYerr, reducedDataCount, 4, fitRes, fitErr, redChiSqrPole);

          if (kappaZeroMode) {
            fitRes[0] *= 10;
	  }
	}
	
        if (b) break;
      }

      fC0P = fitRes[0];
      fC1P = fitRes[1];
      fC2P = fitRes[2];
      fC3P = fitRes[3];
      mG = fitRes[3];
  
      fitRes[3] = fitRes[1];
      fitRes[4] = 96*fitRes[2];
      fitRes[2] = 288*fitRes[2];
      fitRes[5] = mG;
      double m0 = abs(fitRes[1]);

      b = b & findZeroOfBosonic1LoopInvPropagatorOnSecondSheetFit(pole, poleVal, m0, fitRes[0], 2, &(fitRes[2]), 1.2*m0); 
    
      double h = 1E-2;
      Z0 = calcBosonic1LoopInvPropagatorFit(Complex(0,0), m0, fitRes[0], 2, &(fitRes[2])).x;
      Z0 -= calcBosonic1LoopInvPropagatorFit(Complex(0,h), m0, fitRes[0], 2, &(fitRes[2])).x;
      Z0 /= h*h;
      Z0 = 1.0 / Z0;

      mprop0 = sqrt(Z0) * abs(fitRes[1]) / sqrt(fitRes[0]);

      bound1=0;
      bound2=2;
      while (bound2-bound1>1E-10) {
        double middle = 0.5*(bound1+bound2);
        double val = calcBosonic1LoopInvPropagatorFit(Complex(0,middle), m0, fitRes[0], 2, &(fitRes[2])).x;
        if (val>0) {
          bound1 = middle;
        } else {
          bound2 = middle;
        }
      }
      mprop = 0.5*(bound1+bound2);

      Z = calcBosonic1LoopInvPropagatorFit(Complex(0,mprop), m0, fitRes[0], 2, &(fitRes[2])).x;
      Z -= calcBosonic1LoopInvPropagatorFit(Complex(0,mprop+h), m0, fitRes[0], 2, &(fitRes[2])).x;
      Z /= (sqr(mprop+h) - sqr(mprop));
      Z = 1.0 / Z;
    }
  }

  if ((!doArcTanhFit) || (!b) || (isNaN(fitRes[0])) || (isNaN(fitRes[1])) || (isNaN(fitRes[2]))) {
    pole.x = NaN;
    pole.y = NaN;
    poleVal.x = NaN;
    poleVal.y = NaN;
    Z = NaN;
    Z0 = NaN;
    mprop = NaN;
    mprop0 = NaN;
    fC0P = NaN;
    fC1P = NaN;
    fC2P = NaN;
    fC3P = NaN;
  } 
    
  delete[] functionBody;
  delete[] fitRes;
  delete[] fitErr;    
  delete[] y;
  delete[] yerr;
  delete[] reducedDataPSqr;
  delete[] reducedDataY;
  delete[] reducedDataYerr;
}


bool EvaluateObservablePropagatorBase::evaluate() {
  int reps = 0;
  if (dataAvailCount <= 10) {
    reps = dataAvailCount;
  } else if (dataAvailCount <= 100) {
    reps = 10;
  } else {
    reps = ((int) sqrt(dataAvailCount));  
  }
  if (reps>100) reps = 100;
      
  double errorRescaleFactor =  sqrt(reps) * (1.0 - 1.0/reps);
  int blockSize = dataAvailCount / reps;
  if (blockSize<1) blockSize = 1;

  PropagatorEuclideanZFactor = 0.0;
  PropagatorEuclideanZFactorError = 0.0;
  PropagatorEuclideanMass = 0.0;
  PropagatorEuclideanMassError = 0.0;
  PropagatorMinkowskiZFactor = 0.0;
  PropagatorMinkowskiZFactorError = 0.0;
  PropagatorMinkowskiMass = 0.0;
  PropagatorMinkowskiMassError = 0.0;
  PropagatorPoleMass = 0.0;
  PropagatorPoleMassError = 0.0;  
  PropagatorPoleDecayWidth = 0.0;
  PropagatorPoleDecayWidthError = 0.0;  
  PropagatorPoleValue = 0.0;
  PropagatorPoleValueError = 0.0;

  PropagatorEuclideanZFactorReduced = 0;
  PropagatorEuclideanZFactorErrorReduced = 0;
  PropagatorEuclideanMassReduced = 0;
  PropagatorEuclideanMassErrorReduced = 0;
  PropagatorMinkowskiZFactorReduced = 0;
  PropagatorMinkowskiZFactorErrorReduced = 0;
  PropagatorMinkowskiMassReduced = 0;
  PropagatorMinkowskiMassErrorReduced = 0;  
  PropagatorFitConst0Reduced = 0;
  PropagatorFitConst0ErrorReduced = 0;
  PropagatorFitConst1Reduced = 0;
  PropagatorFitConst1ErrorReduced = 0;
  PropagatorFitConst2Reduced = 0;
  PropagatorFitConst2ErrorReduced = 0;

  PropagatorFitConst0ReducedPole = 0;
  PropagatorFitConst0ErrorReducedPole = 0;
  PropagatorFitConst1ReducedPole = 0;
  PropagatorFitConst1ErrorReducedPole = 0;
  PropagatorFitConst2ReducedPole = 0;
  PropagatorFitConst2ErrorReducedPole = 0;
  PropagatorFitConst3ReducedPole = 0;
  PropagatorFitConst3ErrorReducedPole = 0;

  PropagatorZFactor = 0;
  PropagatorZFactorError = 0;  
  PropagatorZ0Factor = 0;
  PropagatorZ0FactorError = 0;  
  PropagatorMProp = 0;
  PropagatorMPropError = 0;  
  PropagatorMProp0 = 0;
  PropagatorMProp0Error = 0;  

  
  if (!isObs("GoldstonePropagator")) {
    double sigmaP0Helper1 = 0;
    double sigmaP0Helper2 = 0;  
    for (int I=0; I<reps; I++) {
      int igStart = I * blockSize;
      int igEnd = (I+1) * blockSize;
      if (igEnd>=dataAvailCount) igEnd = dataAvailCount-1;
      double avgProp0 = 0;
      calcPropAtP0(igStart, igEnd, avgProp0);
      sigmaP0Helper1 += avgProp0;
      sigmaP0Helper2 += sqr(avgProp0);
    }
    sigmaProp[0] = errorRescaleFactor * sqrt(sigmaP0Helper2/reps - sqr(sigmaP0Helper1/reps));
  } else {
    sigmaProp[0] = NaN;
  }
 

  double** dataBuffer = new double*[reps];

  for (int I=0; I<reps; I++) {
    int igStart = I * blockSize;
    int igEnd = (I+1) * blockSize;
    if (igEnd>=dataAvailCount) igEnd = dataAvailCount-1;
    double EucZ, EucM, MinkZ, MinkM, fC0, fC1, fC2;
    double EucZR, EucMR, MinkZR, MinkMR, fC0R, fC1R, fC2R;
    double Z, Z0, mprop, mprop0, fC0P, fC1P, fC2P, fC3P;
    Complex pole, poleVal;
    calcProp(igStart, igEnd, EucZ, EucM, MinkZ, MinkM, fC0, fC1, fC2, EucZR, EucMR, MinkZR, MinkMR, fC0R, fC1R, fC2R,
             Z, Z0, mprop, mprop0, pole, poleVal, fC0P, fC1P, fC2P, fC3P, PropagatorFitReducedChiSquare, PropagatorFitReducedChiSquareReduced, PropagatorFitReducedChiSquarePole);

    dataBuffer[I] = new double[29];
    dataBuffer[I][0] = EucZ;
    dataBuffer[I][1] = EucM;
    dataBuffer[I][2] = MinkZ;
    dataBuffer[I][3] = MinkM;
    dataBuffer[I][4] = fC0;
    dataBuffer[I][5] = fC1;
    dataBuffer[I][6] = fC2;
    dataBuffer[I][7] = EucZR;
    dataBuffer[I][8] = EucMR;
    dataBuffer[I][9] = MinkZR;
    dataBuffer[I][10] = MinkMR;
    dataBuffer[I][11] = fC0R;
    dataBuffer[I][12] = fC1R;
    dataBuffer[I][13] = fC2R;
    dataBuffer[I][14] = Z;
    dataBuffer[I][15] = Z0;
    dataBuffer[I][16] = mprop;
    dataBuffer[I][17] = mprop0;
    dataBuffer[I][18] = pole.x;
    dataBuffer[I][19] = pole.y;
    dataBuffer[I][20] = poleVal.x;
    dataBuffer[I][21] = poleVal.y;
    dataBuffer[I][22] = fC0P;
    dataBuffer[I][23] = fC1P;
    dataBuffer[I][24] = fC2P;
    dataBuffer[I][25] = fC3P;
    dataBuffer[I][26] = PropagatorFitReducedChiSquare;
    dataBuffer[I][27] = PropagatorFitReducedChiSquareReduced;
    dataBuffer[I][28] = PropagatorFitReducedChiSquarePole;

    
    PropagatorEuclideanZFactor += EucZ;
    PropagatorEuclideanZFactorError += EucZ*EucZ;
    PropagatorEuclideanMass += EucM;
    PropagatorEuclideanMassError += EucM*EucM;
    PropagatorMinkowskiZFactor += MinkZ;
    PropagatorMinkowskiZFactorError += MinkZ*MinkZ;
    PropagatorMinkowskiMass += MinkM;
    PropagatorMinkowskiMassError += MinkM*MinkM;
    
    PropagatorPoleMass += abs(pole.y);
    PropagatorPoleMassError += pole.y*pole.y;  
    PropagatorPoleDecayWidth += 2*abs(pole.x);
    PropagatorPoleDecayWidthError += 4*pole.x*pole.x;  
    PropagatorPoleValue += norm(poleVal);
    PropagatorPoleValueError += sqr(norm(poleVal));
    
    PropagatorFitConst0 += fC0;
    PropagatorFitConst0Error += fC0*fC0;
    PropagatorFitConst1 += fC1;
    PropagatorFitConst1Error += fC1*fC1;
    PropagatorFitConst2 += fC2;
    PropagatorFitConst2Error += fC2*fC2;    

    PropagatorEuclideanZFactorReduced += EucZR;
    PropagatorEuclideanZFactorErrorReduced += EucZR*EucZR;
    PropagatorEuclideanMassReduced += EucMR;
    PropagatorEuclideanMassErrorReduced += EucMR*EucMR;
    PropagatorMinkowskiZFactorReduced += MinkZR;
    PropagatorMinkowskiZFactorErrorReduced += MinkZR*MinkZR;
    PropagatorMinkowskiMassReduced += MinkMR;
    PropagatorMinkowskiMassErrorReduced += MinkMR*MinkMR;  
    PropagatorFitConst0Reduced += fC0R;
    PropagatorFitConst0ErrorReduced += fC0R*fC0R;
    PropagatorFitConst1Reduced += fC1R;
    PropagatorFitConst1ErrorReduced += fC1R*fC1R;
    PropagatorFitConst2Reduced += fC2R;
    PropagatorFitConst2ErrorReduced += fC2R*fC2R;

    PropagatorFitConst0ReducedPole += fC0P;
    PropagatorFitConst0ErrorReducedPole += fC0P*fC0P;
    PropagatorFitConst1ReducedPole += fC1P;
    PropagatorFitConst1ErrorReducedPole += fC1P*fC1P;
    PropagatorFitConst2ReducedPole += fC2P;
    PropagatorFitConst2ErrorReducedPole += fC2P*fC2P;
    PropagatorFitConst3ReducedPole += fC3P;
    PropagatorFitConst3ErrorReducedPole += fC3P*fC3P;

    PropagatorZFactor += Z;
    PropagatorZFactorError += Z*Z;  
    PropagatorZ0Factor += Z0;
    PropagatorZ0FactorError += Z0*Z0;  
    PropagatorMProp += mprop;
    PropagatorMPropError += mprop*mprop;  
    PropagatorMProp0 += mprop0;
    PropagatorMProp0Error += mprop0*mprop0;  
  }
  
  PropagatorEuclideanZFactorError = errorRescaleFactor * sqrt(PropagatorEuclideanZFactorError/reps - sqr(PropagatorEuclideanZFactor/reps));
  PropagatorEuclideanMassError = errorRescaleFactor * sqrt(PropagatorEuclideanMassError/reps - sqr(PropagatorEuclideanMass/reps));
  PropagatorMinkowskiZFactorError = errorRescaleFactor * sqrt(PropagatorMinkowskiZFactorError/reps - sqr(PropagatorMinkowskiZFactor/reps));
  PropagatorMinkowskiMassError = errorRescaleFactor * sqrt(PropagatorMinkowskiMassError/reps - sqr(PropagatorMinkowskiMass/reps));

  PropagatorPoleMassError = errorRescaleFactor * sqrt(PropagatorPoleMassError/reps - sqr(PropagatorPoleMass/reps));
  PropagatorPoleDecayWidthError = errorRescaleFactor * sqrt(PropagatorPoleDecayWidthError/reps - sqr(PropagatorPoleDecayWidth/reps));
  PropagatorPoleValueError = errorRescaleFactor * sqrt(PropagatorPoleValueError/reps - sqr(PropagatorPoleValue/reps));

  PropagatorFitConst0Error = errorRescaleFactor * sqrt(PropagatorFitConst0Error/reps - sqr(PropagatorFitConst0/reps));
  PropagatorFitConst1Error = errorRescaleFactor * sqrt(PropagatorFitConst1Error/reps - sqr(PropagatorFitConst1/reps));
  PropagatorFitConst2Error = errorRescaleFactor * sqrt(PropagatorFitConst2Error/reps - sqr(PropagatorFitConst2/reps));  

  PropagatorEuclideanZFactorErrorReduced = errorRescaleFactor * sqrt(PropagatorEuclideanZFactorErrorReduced/reps - sqr(PropagatorEuclideanZFactorReduced/reps));
  PropagatorEuclideanMassErrorReduced = errorRescaleFactor * sqrt(PropagatorEuclideanMassErrorReduced/reps - sqr(PropagatorEuclideanMassReduced/reps));
  PropagatorMinkowskiZFactorErrorReduced = errorRescaleFactor * sqrt(PropagatorMinkowskiZFactorErrorReduced/reps - sqr(PropagatorMinkowskiZFactorReduced/reps));
  PropagatorMinkowskiMassErrorReduced = errorRescaleFactor * sqrt(PropagatorMinkowskiMassErrorReduced/reps - sqr(PropagatorMinkowskiMassReduced/reps));
  PropagatorFitConst0ErrorReduced = errorRescaleFactor * sqrt(PropagatorFitConst0ErrorReduced/reps - sqr(PropagatorFitConst0Reduced/reps));
  PropagatorFitConst1ErrorReduced = errorRescaleFactor * sqrt(PropagatorFitConst1ErrorReduced/reps - sqr(PropagatorFitConst1Reduced/reps));
  PropagatorFitConst2ErrorReduced = errorRescaleFactor * sqrt(PropagatorFitConst2ErrorReduced/reps - sqr(PropagatorFitConst2Reduced/reps));

  PropagatorFitConst0ErrorReducedPole = errorRescaleFactor * sqrt(PropagatorFitConst0ErrorReducedPole/reps - sqr(PropagatorFitConst0ReducedPole/reps));
  PropagatorFitConst1ErrorReducedPole = errorRescaleFactor * sqrt(PropagatorFitConst1ErrorReducedPole/reps - sqr(PropagatorFitConst1ReducedPole/reps));
  PropagatorFitConst2ErrorReducedPole = errorRescaleFactor * sqrt(PropagatorFitConst2ErrorReducedPole/reps - sqr(PropagatorFitConst2ReducedPole/reps));
  PropagatorFitConst3ErrorReducedPole = errorRescaleFactor * sqrt(PropagatorFitConst3ErrorReducedPole/reps - sqr(PropagatorFitConst3ReducedPole/reps));


  PropagatorZFactorError = errorRescaleFactor * sqrt(PropagatorZFactorError/reps - sqr(PropagatorZFactor/reps));
  PropagatorZ0FactorError = errorRescaleFactor * sqrt(PropagatorZ0FactorError/reps - sqr(PropagatorZ0Factor/reps));
  PropagatorMPropError = errorRescaleFactor * sqrt(PropagatorMPropError/reps - sqr(PropagatorMProp/reps));
  PropagatorMProp0Error = errorRescaleFactor * sqrt(PropagatorMProp0Error/reps - sqr(PropagatorMProp0/reps));

  Complex pole, poleVal;  
  
  calcProp(-1, -1, PropagatorEuclideanZFactor, PropagatorEuclideanMass, PropagatorMinkowskiZFactor, PropagatorMinkowskiMass, PropagatorFitConst0, PropagatorFitConst1, PropagatorFitConst2,
            PropagatorEuclideanZFactorReduced, PropagatorEuclideanMassReduced, PropagatorMinkowskiZFactorReduced, PropagatorMinkowskiMassReduced, PropagatorFitConst0Reduced, PropagatorFitConst1Reduced, PropagatorFitConst2Reduced,
            PropagatorZFactor, PropagatorZ0Factor, PropagatorMProp, PropagatorMProp0, pole, poleVal, PropagatorFitConst0ReducedPole, PropagatorFitConst1ReducedPole,
	    PropagatorFitConst2ReducedPole, PropagatorFitConst3ReducedPole, PropagatorFitReducedChiSquare, PropagatorFitReducedChiSquareReduced, PropagatorFitReducedChiSquarePole);

  PropagatorPoleMass = abs(pole.y);
  PropagatorPoleDecayWidth = 2*abs(pole.x);
  PropagatorPoleValue = norm(poleVal);

/*
  //Only consider successful fits...
  if (reps>=6) {
    bool change = true;
    while (change) {
      change = false;
      for (int I=0; I<reps-1; I++) {
        if (dataBuffer[I][28]>dataBuffer[I+1][28]) {
	  change = true;
	  for (int I2=0; I2<29; I2++) {
	    double dummy = dataBuffer[I][I2];
	    dataBuffer[I][I2] = dataBuffer[I+1][I2];
	    dataBuffer[I+1][I2] = dummy;
	  }
	}
      }
    }
    
    int LJpos = 4;
    double LJ = 0;
    
    for (int I=4; I<reps-1; I++) {
      double diff = abs(dataBuffer[I][28] - dataBuffer[I+1][28]);
      if (diff>LJ) {
        LJpos = I;
	LJ = diff;
      }
    }
    
    //Determine sigma of ChiSqr for data below LJ
    double avgChiSqr = 0;
    double sigmaChiSqr = 0;
    for (int I=0; I<LJpos+1; I++) {
      avgChiSqr += dataBuffer[I][28];
      sigmaChiSqr += sqr(dataBuffer[I][28]);
    }
    avgChiSqr /= (LJpos+1);
    sigmaChiSqr = sqrt(sigmaChiSqr/(LJpos+1) - sqr(avgChiSqr));
     
    //If LJ>3*sigma then consider only up to LJpos
    int considerIndex = reps;   
    if (LJ > 3*sigmaChiSqr) considerIndex = LJpos+1;   


FILE* file = NULL;

if (isObs("GoldstonePropagator")) {
  file=fopen("PropOutputGoldstone.dat","w");
} else {
  file=fopen("PropOutputHiggs.dat","w");
}
for (int I=0; I<reps; I++) {
  for (int I2=0; I2<29; I2++) {
    fprintf(file, "%f ", dataBuffer[I][I2]);
  }

  fprintf(file, "\n");

}
fprintf(file, "%f %d %f %f %d\n", LJ, LJpos, avgChiSqr,sigmaChiSqr, considerIndex);
fclose(file);

    //Average results up to considerIndex
    PropagatorFitReducedChiSquarePole = 0;
    
    for (int I=0; I<considerIndex; I++) {
      PropagatorFitReducedChiSquarePole += dataBuffer[I][28];    
    }
        
    PropagatorFitReducedChiSquarePole /= considerIndex;
    
    calcAvgAndSigmaFromSelectedDataBuffer(considerIndex, 14, dataBuffer, PropagatorZFactor, PropagatorZFactorError);
    calcAvgAndSigmaFromSelectedDataBuffer(considerIndex, 15, dataBuffer, PropagatorZ0Factor, PropagatorZ0FactorError);
    calcAvgAndSigmaFromSelectedDataBuffer(considerIndex, 17, dataBuffer, PropagatorMProp0, PropagatorMProp0Error);
    calcAvgAndSigmaFromSelectedDataBuffer(considerIndex, 16, dataBuffer, PropagatorMProp, PropagatorMPropError);
    
    PropagatorZFactorError *= errorRescaleFactor;
    PropagatorZ0FactorError *= errorRescaleFactor;
    PropagatorMPropError *= errorRescaleFactor;
    PropagatorMProp0Error *= errorRescaleFactor;
  }
  */

  return true;
}


LAPsystemPlot* EvaluateObservablePropagatorBase::createPlot1(double maxP) {
  char* name = new char[1000];
  snprintf(name,1000,"%smaxP%1.2f",getObsName(),maxP);
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  int startP = 1;
  if (considerZeroMomentum) startP = 0;
  char* plotCmd = new char[10000];
  double** plotData = new double*[latticeBins->getMomentumSqrSlotCount()-startP];
  for (int I=startP; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    plotData[I-startP] = new double[3];
    plotData[I-startP][0] = pSqr[I];
    plotData[I-startP][1] = 1.0 / avgProp[I];
    plotData[I-startP][2] = abs(sigmaProp[I] / sqr(avgProp[I]));
  }
  plot->setPlotData(latticeBins->getMomentumSqrSlotCount()-startP, 3, plotData);
  plot->setXLabel("$\\\\hat p^2$");
  plot->setYLabel("Propagator");
  plot->setCaption("Caption");
  plot->setTitle("");
  plot->setPlotTitle(name);
  plot->setYLogScale(false);
  plot->setSize(0.8, 0.8);
  plot->setXRange(0, maxP);
  plot->setYErrorBars(true);  
  plot->setPointSize(0.3);
  plot->setPointType(5);
  plot->setLineType(0);  
  plot->plotData("1:2:3");

  //Show only linear and double arctanh-fits of reduced data [that is up to threshold value gamma].
//  if ((!isNaN(PropagatorEuclideanMass)) && (!isNaN(PropagatorFitConst0)) && (!isNaN(PropagatorFitConst1)) && (!isNaN(PropagatorFitConst2)) && (!isNaN(PropagatorEuclideanZFactor)) && (!isNaN(PropagatorMinkowskiMass)) && (!isNaN(PropagatorMinkowskiZFactor))) {
//    snprintf(plotCmd,10000,"replot (x+%1.15e*%1.15e+%1.15e*x/(x+%1.15e))/%1.15e notitle", PropagatorEuclideanMass, PropagatorEuclideanMass, PropagatorFitConst1, PropagatorFitConst2, PropagatorFitConst0);
//    plot->plotDirect(plotCmd);
//  }

  if ((!isNaN(PropagatorEuclideanMassReduced)) && (!isNaN(PropagatorFitConst0Reduced)) && (!isNaN(PropagatorFitConst1Reduced)) && (!isNaN(PropagatorFitConst2Reduced)) && (!isNaN(PropagatorEuclideanZFactorReduced)) && (!isNaN(PropagatorMinkowskiMassReduced)) && (!isNaN(PropagatorMinkowskiZFactorReduced))) {
    snprintf(plotCmd,10000,"replot (x+%1.15e*%1.15e+%1.15e*x/(x+%1.15e))/%1.15e notitle", PropagatorEuclideanMassReduced, PropagatorEuclideanMassReduced, PropagatorFitConst1Reduced, PropagatorFitConst2Reduced, PropagatorFitConst0Reduced);
    plot->plotDirect(plotCmd);
  }

  if ((!isNaN(PropagatorFitConst0ReducedPole)) && (!isNaN(PropagatorFitConst1ReducedPole)) && (!isNaN(PropagatorFitConst2ReducedPole)) && (!isNaN(PropagatorFitConst3ReducedPole))) {
    if (isObs("GoldstonePropagator")) {
      double A1 = PropagatorFitConst0ReducedPole;
      double A2 = PropagatorFitConst1ReducedPole;
      double A3 = PropagatorFitConst2ReducedPole;
      double A4 = PropagatorFitConst3ReducedPole;
      double eps = 1E-10;

      snprintf(plotCmd,10000,"replot (x+%1.15f*%1.15f + %1.15f*(( log(%1.15f*%1.15f/(%1.15f*%1.15f+%1.15f))*(%1.15f*%1.15f+%1.15f-%1.15f*%1.15f)  \
                                                + sqrt(((%1.15f*%1.15f+%1.15f-%1.15f*%1.15f) + x)**2 + 4*%1.15f*%1.15f*x) \
						* log( ((sqrt(((%1.15f*%1.15f+%1.15f-%1.15f*%1.15f) + x)**2 + 4*%1.15f*%1.15f*x) + x)**2 - (%1.15f*%1.15f+%1.15f-%1.15f*%1.15f)**2) \
						     / ((sqrt(((%1.15f*%1.15f+%1.15f-%1.15f*%1.15f) + x)**2 + 4*%1.15f*%1.15f*x) - x)**2 - (%1.15f*%1.15f+%1.15f-%1.15f*%1.15f)**2) ) ) / x \
						- (2 + log((%1.15f*%1.15f+%1.15f)/(%1.15f*%1.15f)) * (1+2*%1.15f*%1.15f/(%1.15f*%1.15f+%1.15f-%1.15f*%1.15f))) )   ) / %1.15f notitle",
                                                A2,A2,A3,A4,A4,A2,A2,eps,A2,A2,eps,A4,A4,A2,A2,eps,A4,A4,A4,A4,A2,A2,eps,A4,A4,A4,A4,A2,A2,eps,A4,A4,A2,A2,eps,A4,
						A4,A4,A4,A2,A2,eps,A4,A4,A2,A2,eps,A4,A4,A4,A4,A2,A2,eps,A4,A4,A1);

      plot->plotDirect(plotCmd);
    } else {
      double eps = 1E-10;
      double A1 = PropagatorFitConst0ReducedPole;
      double A2 = PropagatorFitConst1ReducedPole;
      double A3 = PropagatorFitConst2ReducedPole;
      double mG = PropagatorFitConst3ReducedPole;

      snprintf(plotCmd,10000,"replot (x+%1.15f*%1.15f + %1.15f*(288*(sqrt((x+4*%1.15f*%1.15f+%1.15f)/(x+%1.15f))*atanh(sqrt((x+%1.15f)/(x+4*%1.15f*%1.15f+%1.15f))) -1) \
                                                + 96*(sqrt((x+4*%1.15f*%1.15f+%1.15f)/(x+%1.15f))*atanh(sqrt((x+%1.15f)/(x+4*%1.15f*%1.15f+%1.15f))) -1)) )/%1.15f notitle",
                                                A2,A2,A3,A2,A2,eps,eps,eps,A2,A2,eps,mG,mG,eps,eps,eps,mG,mG,eps,A1);
      plot->plotDirect(plotCmd);
    }
  }
    
  for (int I=0; I<latticeBins->getMomentumSqrSlotCount()-startP; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


LAPsystemPlot* EvaluateObservablePropagatorBase::createPlot2(double maxP, double subtractMSqr, bool considerSubTerm) {
  char* name = new char[1000];
  snprintf(name,1000,"%sSelfEnergymaxP%1.2f",getObsName(),maxP);
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  int startP = 1;
  if (considerZeroMomentum) startP = 0;
  char* plotCmd = new char[1000];
  double** plotData = new double*[latticeBins->getMomentumSqrSlotCount()-startP];
  for (int I=startP; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    plotData[I-startP] = new double[3];
    plotData[I-startP][0] = pSqr[I];
    plotData[I-startP][1] = (1/*PropagatorZFactor */ / avgProp[I]) - pSqr[I] - subtractMSqr;
    plotData[I-startP][2] = sqrt(sqr(sigmaProp[I] / sqr(avgProp[I])) /*+ sqr(PropagatorZFactorError/PropagatorZFactor)*/);
  }
  plot->setPlotData(latticeBins->getMomentumSqrSlotCount()-startP, 3, plotData);
  plot->setXLabel("$\\\\hat p^2$");
  plot->setYLabel("Neg. Self-Energy $-\\\\Sigma(p^2)$");
  char* caption = new char[1000];
  if (considerSubTerm) {
    snprintf(caption, 1000, "Subtracted constant: $m_0^2+12\\lambda v^2=%1.8f$", subtractMSqr);
  } else {
    snprintf(caption, 1000, "Subtracted constant: $%f$", subtractMSqr);  
  }
  plot->setCaption(caption);
  delete[] caption;
  plot->setTitle("");
  plot->setPlotTitle(name);
  plot->setYLogScale(false);
  plot->setSize(1.0, 1.0);
  plot->setXRange(0, maxP);
  plot->setYErrorBars(true);  
  plot->setPointSize(0.3);
  plot->setPointType(5);
  plot->setLineType(0);  
  plot->plotData("1:2:3");
  
  for (int I=0; I<latticeBins->getMomentumSqrSlotCount()-startP; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


LAPsystemPlot* EvaluateObservablePropagatorBase::createPlot3() {
  char* name = new char[1000];
  snprintf(name,1000,"%sAutoCorrelationTimes",getObsName());
  LAPsystemPlot* plot = LAPsystem->createNewPlot(name);

  int startP = 0;
//  if (considerZeroMomentum) startP = 0;
  char* plotCmd = new char[1000];
  double** plotData = new double*[latticeBins->getMomentumSqrSlotCount()-startP];
  for (int I=startP; I<latticeBins->getMomentumSqrSlotCount(); I++) {
    plotData[I-startP] = new double[3];
    plotData[I-startP][0] = pSqr[I];
    plotData[I-startP][1] = autoCorrelationTime[I];
    plotData[I-startP][2] = 0;
  }
  plot->setPlotData(latticeBins->getMomentumSqrSlotCount()-startP, 3, plotData);
  plot->setXLabel("$\\\\hat p^2$");
  plot->setYLabel("Auto-Correlation Times");
  plot->setCaption("");
  plot->setTitle("");
  plot->setPlotTitle(name);
  plot->setYLogScale(false);
  plot->setSize(1.0, 1.0);
//  plot->setXRange(0, maxP);
  plot->setYErrorBars(false);  
  plot->setPointSize(0.3);
  plot->setPointType(5);
  plot->setLineType(0);  
  plot->plotData("1:2");
  
  for (int I=0; I<latticeBins->getMomentumSqrSlotCount()-startP; I++) {
    delete[] plotData[I];
  }
  delete[] plotData;
  delete[] plotCmd;
  delete[] name;
  
  return plot;
}


void EvaluateObservablePropagatorBase::generateLatexAndPlotsAndXML() {
  LAPsystemPlot* plot1 = createPlot1(16);  
  LAPsystem->addPlot(plot1);

  LAPsystemPlot* plot2 = createPlot1(1);  
  LAPsystem->addPlot(plot2);
  
  if (!isNaN(selfEnergySubLamFac)) {
    double m0Sqr = SimParaSet->getM0Squared();
    double fac = selfEnergySubLamFac;
    double lambda = SimParaSet->getLambda0();
    bool considerSubTerm = true;
    
    double subConst =  m0Sqr + fac * lambda * sqr(vevContNot);
    if (isNaN(subConst)) {
      subConst = 0;
      considerSubTerm = false;
    }
  
    LAPsystemPlot* plot3 = createPlot2(16, subConst, considerSubTerm);  
    LAPsystem->addPlot(plot3);
    LAPsystemPlot* plot4 = createPlot2(1, subConst, considerSubTerm);  
    LAPsystem->addPlot(plot4);    
  }

  LAPsystemPlot* plot5 = createPlot3();  
  LAPsystem->addPlot(plot5);

  startLatexOutputSummaryTable();
  
  
  addXML_And_LatexOutputSummaryTableLine("ZFactor", "Field renormalization factor Z", "$Z$", getPropagatorZFactor(), getPropagatorZFactorError(), NULL, "%1.5f");
  addXML_And_LatexOutputSummaryTableLine("LatPropMass", "Pole mass of propagator in lattice units", "$m_{lat}$", getPropagatorMass(), getPropagatorMassError(), NULL, "%1.3f");  
  if (!isNaN(physicalScaleInGEV)) {
    addXML_And_LatexOutputSummaryTableLine("PhysPropMass", "Pole mass of propagator in GeV", "$m$",  physicalScaleInGEV*getPropagatorMass(), sqrt(sqr(physicalScaleInGEV*getPropagatorMassError()) + sqr(physicalScaleErrorInGEV*getPropagatorMass())), "GeV", "%1.1f");
  }
  addXML_And_LatexOutputSummaryTableLine("ZFactorEuc", "Euclidean field renormalization factor Z (derivative)", "$Z$", PropagatorEuclideanZFactor, PropagatorEuclideanZFactorError, NULL, "%1.5f");
  addXML_And_LatexOutputSummaryTableLine("LatPropMassEuc", "Euclidean Pole mass of propagator in lattice units", "$m_{lat}$", PropagatorEuclideanMass, PropagatorEuclideanMassError, NULL, "%1.3f");  
  if (!isNaN(physicalScaleInGEV)) {
    addXML_And_LatexOutputSummaryTableLine("PhysPropMassEuc", "Euclidean Pole mass of propagator in GeV", "$m$", physicalScaleInGEV*PropagatorEuclideanMass, sqrt(sqr(physicalScaleInGEV*PropagatorEuclideanMassError) + sqr(physicalScaleErrorInGEV*PropagatorEuclideanMass)), "GeV", "%1.1f");
  }  
  addXML_And_LatexOutputSummaryTableLine("ZFactorMink", "Minkowski field renormalization factor Z (derivative)", "$Z$", PropagatorMinkowskiZFactor, PropagatorMinkowskiZFactorError, NULL, "%1.5f");
  addXML_And_LatexOutputSummaryTableLine("LatPropMassMink", "Minkowski Pole mass of propagator in lattice units", "$m_{lat}$", PropagatorMinkowskiMass, PropagatorMinkowskiMassError, NULL, "%1.3f");  
  if (!isNaN(physicalScaleInGEV)) {
    addXML_And_LatexOutputSummaryTableLine("PhysPropMassMink", "Minkowski Pole mass of propagator in GeV", "$m$", physicalScaleInGEV*PropagatorMinkowskiMass, sqrt(sqr(physicalScaleInGEV*PropagatorMinkowskiMassError) + sqr(physicalScaleErrorInGEV*PropagatorMinkowskiMass)), "GeV", "%1.1f");
  }
  addXML_And_LatexOutputSummaryTableLine("RedChiSqrAlgFit", "Reduced $\\chi^2$ from algebraic fit", "$\\chi^2_{alg}$", PropagatorFitReducedChiSquare, 0, NULL, "%1.5f");


  addXML_And_LatexOutputSummaryTableLine("ZFactorEucRed", "Euclidean field ren. factor Z (derivative) (lin)", "$Z$",  PropagatorEuclideanZFactorReduced, PropagatorEuclideanZFactorErrorReduced, NULL, "%1.5f");
  addXML_And_LatexOutputSummaryTableLine("LatPropMassEucRed", "Euclidean Pole mass of prop. in lattice units (lin)", "$m_{lat}$", PropagatorEuclideanMassReduced, PropagatorEuclideanMassErrorReduced, NULL, "%1.3f");  
  if (!isNaN(physicalScaleInGEV)) {
    addXML_And_LatexOutputSummaryTableLine("PhysPropMassEucRed", "Euclidean Pole mass of prop. in GeV (lin)", "$m$", physicalScaleInGEV*PropagatorEuclideanMassReduced, sqrt(sqr(physicalScaleInGEV*PropagatorEuclideanMassErrorReduced) + sqr(physicalScaleErrorInGEV*PropagatorEuclideanMassReduced)), "GeV", "%1.1f");
  }  
  addXML_And_LatexOutputSummaryTableLine("ZFactorMinkRed", "Minkowski field ren. factor Z (derivative) (lin)", "$Z$", PropagatorMinkowskiZFactorReduced, PropagatorMinkowskiZFactorErrorReduced, NULL, "%1.5f");
  addXML_And_LatexOutputSummaryTableLine("LatPropMassMinkRed", "Minkowski Pole mass of prop. in lattice units (lin)", "$m_{lat}$", PropagatorMinkowskiMassReduced, PropagatorMinkowskiMassErrorReduced, NULL, "%1.3f");  
  if (!isNaN(physicalScaleInGEV)) {
    addXML_And_LatexOutputSummaryTableLine("PhysPropMassMinkRed", "Minkowski Pole mass of prop. in GeV (lin)", "$m$", physicalScaleInGEV*PropagatorMinkowskiMassReduced, sqrt(sqr(physicalScaleInGEV*PropagatorMinkowskiMassErrorReduced) + sqr(physicalScaleErrorInGEV*PropagatorMinkowskiMassReduced)), "GeV", "%1.1f");
  }
  addXML_And_LatexOutputSummaryTableLine("RedChiSqrAlgFitRed", "Reduced $\\chi^2$ from algebraic fit (lin)", "$\\chi^2_{alg}$", PropagatorFitReducedChiSquareReduced, 0, NULL, "%1.5f");

  if (doArcTanhFit) {
    addXML_And_LatexOutputSummaryTableLine("ZFactorArcFit", "Field renormalization factor Z from arctanf fit", "$Z$", PropagatorZFactor, PropagatorZFactorError, NULL, "%1.5f");
    addXML_And_LatexOutputSummaryTableLine("Z0FactorArcFit", "Field renormalization factor Z0 from arctanf fit", "$Z_0$", PropagatorZ0Factor, PropagatorZ0FactorError, NULL, "%1.5f");
    addXML_And_LatexOutputSummaryTableLine("LatPropMassArcFit", "Propagator mass in lattice units from arctanf fit", "$m_{lat}$", PropagatorMProp, PropagatorMPropError, NULL, "%1.3f");  
    if (!isNaN(physicalScaleInGEV)) {
      addXML_And_LatexOutputSummaryTableLine("PhysPropMassArcFit", "Propagator mass in GeV from arctanf fit", "$m$", physicalScaleInGEV*PropagatorMProp, sqrt(sqr(physicalScaleInGEV*PropagatorMPropError) + sqr(physicalScaleErrorInGEV*PropagatorMProp)), "GeV", "%1.1f");
    }
    addXML_And_LatexOutputSummaryTableLine("LatPropMass0ArcFit", "Propagator mass0 in lattice units from arctanf fit", "$m_{lat}$", PropagatorMProp0, PropagatorMProp0Error, NULL, "%1.3f");  
    if (!isNaN(physicalScaleInGEV)) {
      addXML_And_LatexOutputSummaryTableLine("PhysPropMass0ArcFit", "Propagator mass0 in GeV from arctanf fit", "$m$", physicalScaleInGEV*PropagatorMProp0, sqrt(sqr(physicalScaleInGEV*PropagatorMProp0Error) + sqr(physicalScaleErrorInGEV*PropagatorMProp0)), "GeV", "%1.1f");
    }

    if (!isObs("GoldstonePropagator")) {
      addXML_And_LatexOutputSummaryTableLine("LatPropMassPole", "Complex extrapolation of Pole mass in lattice units", "$m_{lat}$", PropagatorPoleMass, PropagatorPoleMassError, NULL, "%1.3f");  
      if (!isNaN(physicalScaleInGEV)) {
        addXML_And_LatexOutputSummaryTableLine("PhysPropMassPole", "Complex extrapolation of Pole mass in GeV", "$m$", physicalScaleInGEV*PropagatorPoleMass, sqrt(sqr(physicalScaleInGEV*PropagatorPoleMassError) + sqr(physicalScaleErrorInGEV*PropagatorPoleMass)), "GeV", "%1.1f");
      }
      addXML_And_LatexOutputSummaryTableLine("LatDecayWidthPole", "Complex extrapolation of decay width in lattice units", "$\\Gamma_{lat}$", PropagatorPoleDecayWidth, PropagatorPoleDecayWidthError, NULL, "%1.3f");  
      if (!isNaN(physicalScaleInGEV)) {
        addXML_And_LatexOutputSummaryTableLine("PhysDecayWidthPole", "Complex extrapolation of decay width in GeV", "$\\Gamma$", physicalScaleInGEV*PropagatorPoleDecayWidth, sqrt(sqr(physicalScaleInGEV*PropagatorPoleDecayWidthError) + sqr(physicalScaleErrorInGEV*PropagatorPoleDecayWidth)), "GeV", "%1.1f");
      } 
      addXML_And_LatexOutputSummaryTableLine("LatInvValAtPole", "Inverse value at pole", "$G^{-1}(pole)$", PropagatorPoleValue, PropagatorPoleValueError, NULL, "%1.3f");  
    }
    addXML_And_LatexOutputSummaryTableLine("RedChiSqrArcTanHFit", "Reduced $\\chi^2$ from double-arctanh fit", "$\\chi^2_{atanh}$", PropagatorFitReducedChiSquarePole, 0, NULL, "%1.5f");
  }
  
  addXML_And_LatexOutputSummaryTableLine("gamma", "Gamma-Value for fit procedure", "$\\gamma$", gamma, 0, NULL, "%1.2f");
 
  endLatexOutputSummaryTable();
}


double EvaluateObservablePropagatorBase::getPropagatorEuclideanZFactor() {
  return PropagatorEuclideanZFactor;
}


double EvaluateObservablePropagatorBase::getPropagatorEuclideanZFactorError() {
  return PropagatorEuclideanZFactorError;
}


double EvaluateObservablePropagatorBase::getPropagatorEuclideanMass() {
  return PropagatorEuclideanMass;
}


double EvaluateObservablePropagatorBase::getPropagatorEuclideanMassError() {
  return PropagatorEuclideanMassError;
}


double EvaluateObservablePropagatorBase::getPropagatorMinkowskiZFactor() {
  return PropagatorMinkowskiZFactor;
}


double EvaluateObservablePropagatorBase::getPropagatorMinkowskiZFactorError() {
  return PropagatorMinkowskiZFactorError;
}


double EvaluateObservablePropagatorBase::getPropagatorMinkowskiMass() {
  return PropagatorMinkowskiMass;
}


double EvaluateObservablePropagatorBase::getPropagatorMinkowskiMassError() {
  return PropagatorMinkowskiMassError;
}


double EvaluateObservablePropagatorBase::getPropagatorZFactor() {
  return PropagatorZFactor;
}


double EvaluateObservablePropagatorBase::getPropagatorZFactorError() {
  return PropagatorZFactorError;
}


double EvaluateObservablePropagatorBase::getPropagatorMass() {
  return PropagatorMProp;
}


double EvaluateObservablePropagatorBase::getPropagatorMassError() {
  return PropagatorMPropError;
}


void EvaluateObservablePropagatorBase::calcAvgAndSigmaFromSelectedDataBuffer(int N, int ind, double** dataBuffer, double& avg, double& sig) {
  double* data = new double[N];
  for (int I=0; I<N; I++) data[I] = dataBuffer[I][ind];
  calcAverageAndStandardDeviationWithDataSelection(N, data, 10, avg, sig);
  delete[] data;
}
