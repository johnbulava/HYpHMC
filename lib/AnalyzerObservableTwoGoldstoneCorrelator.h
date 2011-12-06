#ifndef AnalyzerObservableTwoGoldstoneCorrelator_included
#define AnalyzerObservableTwoGoldstoneCorrelator_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "Global.h"
#include "AnalyzerIOControl.h"
#include "FermionMatrixOperations.h"
#include "StateDescriptorReader.h"
#include "AnalyzerObservable.h"
#include "AnalyzerPhiFieldConfiguration.h"
#include "LatticeMomentumBins.h"

class AnalyzerObservableTwoGoldstoneCorrelator : public AnalyzerObservable {
	private:
		void performTimeEnergyFourierTransform(bool dir, double* phiSource, double* phiDest, int L0, int L1, int L2, int L3);
		inline int AnalyzerObservableTwoGoldstoneCorrelator::transformCoordinatesToIndex(int a, int n0, int n1, int n2, int n3);
		
	public:
		AnalyzerObservableTwoGoldstoneCorrelator(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDreader); 
		~AnalyzerObservableTwoGoldstoneCorrelator();
		
		bool analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors);
		int getNeededAuxVectorCount();
		int getAnalyzerResultsCount();  
};


#include "AnalyzerObservableTwoGoldstoneCorrelator.C"

#endif
