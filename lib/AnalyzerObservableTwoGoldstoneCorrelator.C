#include "AnalyzerObservableTwoGoldstoneCorrelator.h"

AnalyzerObservableTwoGoldstoneCorrelator::AnalyzerObservableTwoGoldstoneCorrelator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader) : AnalyzerObservable(fOps, aIOcon, SDreader, "TwoGoldstoneCorrelator", "twogoldstonecorr") {
	ini(getAnalyzerResultsCount());
}


AnalyzerObservableTwoGoldstoneCorrelator::~AnalyzerObservableTwoGoldstoneCorrelator() {
}


bool AnalyzerObservableTwoGoldstoneCorrelator::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
	for (int I=0; I<getAnalyzerResultsCount(); I++) {
		analyzerResults[I] = 0;
	}
	
	double* phiField = phiFieldConf->getPhiFieldCopy(); 
	double* phiFieldTransform = phiFieldConf->getPhiFieldCopy();
	double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);
	phiFieldConf->multiplyHiggsFieldWithConst(phiField, rescale);
	phiFieldConf->multiplyHiggsFieldWithConst(phiFieldTransform, rescale);
	phiFieldConf->performFourierTransform(phiFieldTransform, true, 4);	
	
	int* momentumP = new int[4];
	momentumP[0] =  0;
	momentumP[1] = +1;
	momentumP[2] = -1;
	momentumP[3] =  0;	
	
	int* momentumTP = new int[4];
	momentumTP[0] = momentumP[0];
	momentumTP[1] = momentumP[1];
	momentumTP[2] = momentumP[2];
	momentumTP[3] = momentumP[3];
	
	int L0 = fermiOps->get1DSizeL0();
	int L1 = fermiOps->get1DSizeL1();
	int L2 = fermiOps->get1DSizeL2();
	int L3 = fermiOps->get1DSizeL3();
	int timeDirection = fermiOps->getTimeDirection();
	int LargestL = fermiOps->get1DSizeLargest();
	  
	Complex* phiFieldOfTAndP = new Complex[4*L0*L1*L2*L3];
	for (int i = 0; i < 4*L0*L1*L2*L3; i++) {
		phiFieldOfTAndP[i].x = phiFieldTransform[i];
	}
	int indexP = 0;
	int indexTP = 0;
	double expArg = 0.0;
	for (int a = 0; a < 4; a++) {
		for (int t = 0; t < LargestL; t++) {
			momentumTP[timeDirection] = t;
			indexTP = transformCoordinatesToIndex(a, momentumTP[0], momentumTP[1], momentumTP[2], momentumTP[3]);
			for (momentumP[timeDirection] = 0; momentumP[timeDirection] < LargestL; momentumP[timeDirection]++) {
				indexP = transformCoordinatesToIndex(a, momentumP[0], momentumP[1], momentumP[2], momentumP[3]);
				expArg = 2.0 * pi * momentumP[timeDirection] * t / LargestL;
				phiFieldOfTAndP[indexTP] = phiFieldTransform[indexP] * exp(expArg*ComplexI); 
			}
		}
	}
	
	Complex* phiOfT_0 = new Complex[LargestL];	 
	Complex* phiOfT_1 = new Complex[LargestL];
	Complex* phiOfT_2 = new Complex[LargestL];
	Complex* phiOfT_3 = new Complex[LargestL];
	
	Complex* phiDaggerOfT_0 = new Complex[LargestL];	 
	Complex* phiDaggerOfT_1 = new Complex[LargestL];
	Complex* phiDaggerOfT_2 = new Complex[LargestL];
	Complex* phiDaggerOfT_3 = new Complex[LargestL];
	
	
	for (int i = 0; i < LargestL; i++) {
		phiOfT_0[i] = ComplexZero;
		phiOfT_1[i] = ComplexZero;
		phiOfT_2[i] = ComplexZero;
		phiOfT_3[i] = ComplexZero;
		
		phiDaggerOfT_0[i] = ComplexZero;
		phiDaggerOfT_1[i] = ComplexZero;
		phiDaggerOfT_2[i] = ComplexZero;
		phiDaggerOfT_3[i] = ComplexZero;
	}
	
	
	Complex expFuncOfP = ComplexUnity;
	Complex expFuncOfMinusP = ComplexUnity;
	expArg = 0.0;
	
	int* index = new int[4];
	int count = 0;
	for (index[0] = 0; index[0] < L0; index[0]++) {
		expArg += momentumP[0] * index[0] * 2.0 * pi/L0;
		for (index[1] = 0; index[1] < L1; index[1]++) {
			expArg += momentumP[1] * index[1] * 2.0 * pi/L1;
			for (index[2] = 0; index[2] < L2; index[2]++) {
				expArg += momentumP[2] * index[2] * 2.0 * pi/L2;
				for (index[3] = 0; index[3] < L3; index[3]++) {
					expArg += momentumP[3] * index[3] * 2.0 * pi/L3;
					expArg = expArg - momentumP[timeDirection]*index[timeDirection] * 2.0 * pi/LargestL;
					expFuncOfP = exp(expArg * ComplexI);
					expFuncOfMinusP = exp(-1.0 * expArg * ComplexI);
					
					phiOfT_0[index[timeDirection]] = phiOfT_0[index[timeDirection]] + expFuncOfP * phiField[4*count+0];
					phiOfT_1[index[timeDirection]] = phiOfT_0[index[timeDirection]] + expFuncOfP * phiField[4*count+1];
					phiOfT_2[index[timeDirection]] = phiOfT_0[index[timeDirection]] + expFuncOfP * phiField[4*count+2];
					phiOfT_3[index[timeDirection]] = phiOfT_0[index[timeDirection]] + expFuncOfP * phiField[4*count+3];
					
					phiDaggerOfT_0[index[timeDirection]] = phiDaggerOfT_0[index[timeDirection]] + expFuncOfMinusP * phiField[4*count+0];
					phiDaggerOfT_1[index[timeDirection]] = phiDaggerOfT_0[index[timeDirection]] + expFuncOfMinusP * phiField[4*count+1];
					phiDaggerOfT_2[index[timeDirection]] = phiDaggerOfT_0[index[timeDirection]] + expFuncOfMinusP * phiField[4*count+2];
					phiDaggerOfT_3[index[timeDirection]] = phiDaggerOfT_0[index[timeDirection]] + expFuncOfMinusP * phiField[4*count+3];
					count++;
				}
			}
		}
	}
		
	Complex phiSqr = ComplexUnity;
	for (int i=0; i < getAnalyzerResultsCount(); i++) {
		phiSqr = phiOfT_1[i] * phiDaggerOfT_1[i] + phiOfT_2[i] * phiDaggerOfT_2[i] + phiOfT_3[i] * phiDaggerOfT_3[i];
		analyzerResults[i] = 1.0/3.0 * phiSqr.x;
		analyzerResults[i] /= L0*L1*L2*L3;
		analyzerResults[i] *= LargestL;
	}
		
	delete[] momentumP;
	delete[] momentumTP;
	
	delete[] phiOfT_0;	
	delete[] phiOfT_1;
	delete[] phiOfT_2;
	delete[] phiOfT_3;
	
	delete[] phiDaggerOfT_0;	
	delete[] phiDaggerOfT_1;
	delete[] phiDaggerOfT_2;
	delete[] phiDaggerOfT_3;
	
	delete[] phiFieldOfTAndP;
	return true;
}

int AnalyzerObservableTwoGoldstoneCorrelator::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableTwoGoldstoneCorrelator::getAnalyzerResultsCount() {
  int LargestL = fermiOps->get1DSizeLargest();
  return LargestL;
}

void AnalyzerObservableTwoGoldstoneCorrelator::performTimeEnergyFourierTransform(bool dir, double* phiSource, double* phiDest, int L0, int L1, int L2, int L3) {
	
}


