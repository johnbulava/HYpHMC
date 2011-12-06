#ifndef StateDescriptorReader_included
#define StateDescriptorReader_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdarg.h>
#include <sys/stat.h>

class StateDescriptorReader {
	private:
		char* fileName;
		va_list args;
		bool fileNameIsValid;
		char* stateDescriptorContent;
		
		void read();
		
	public:
		StateDescriptorReader(const char* stateDescriptorPath); 
		~StateDescriptorReader();
		
		void checkEndLineCharacter(char* str);
		char* getStrParameter(int tagCount, char** tags);
		bool isValidFileName();  
		double getDoubleParameter(int tagCount, ...);
		double getDoubleParameter(double defaultValue, int tagCount, ...);
		int getIntParameter(int tagCount, ...);
		int getIntParameter(int defaultValue, int tagCount, ...);
		char* getStrParameter(int tagCount, ...);  //creates a new char*-String. Must be deleted by calling function
		
		//Helper-Functions: simply specify tags and use above functions
		char* getFileNameSuffix();
		char* getFileNameExtension();
		int getL0();
		int getL1();
		int getL2();
		int getL3();
		int getNf();
		double getKappa();
		double getLambda();
		double getRho();
		double getR();
		double getYN();  
		bool getUseP();
		bool getUseQ();
		bool getUseR();
		bool getUseQHM();
		double getRPrecondParameterM();
		double getRPrecondParameterF();
		double getPPrecondParameterM();
		double getPPrecondParameterS();
		double getMassSplit();
		double getExternalCurrent();
		
		double getAveragePhi();
		double getAverageStaggeredPhi();
		double getExplicitFermionMass();
		double getPolynomLowerBound_P0();
		double getPolynomLowerBound_P1();
		double getPolynomLowerBound_P2();
		double getPolynomLowerBound_P3();
		double getPolynomLowerBound_P4();
		double getPolynomUpperBound();
		int getPolynomDegree_P0();
		int getPolynomDegree_P1();
		int getPolynomDegree_P2();
		int getPolynomDegree_P3();
		int getPolynomDegree_P4();
		double getAlpha(); //Exponent of function to be approximated: f(x) = 1/x^alpha
		int getPolynomDigits();
		int getMaximumPolynomDegreePerNode();
		bool useRandomGauge();
		bool useExactReweighing();
		double getUpperEWBoundSafetyFactor();
		int getAdditionalAuxVectorCount();
		bool useDirectSamplingOfOmegaFields();
		double getTheta();
		bool useThetaScan();
		double getThetaMin();
		double getThetaMax();
		double getThetaAdaptionFactor();
		bool useModelSelection();
		bool useXFFT();
		double getMolecularDynamicsEpsilon();
		int getSubPolynomCount();
		int getPolynomIterations_P0();
		int getPolynomIterations_P1();
		int getPolynomIterations_P2();
		int getPolynomIterations_P3();
		int getPolynomIterations_P4();
		int getHiggsDynamicsIterations();
		int getPolynomIntegrationType_P0();
		int getPolynomIntegrationType_P1();
		int getPolynomIntegrationType_P2();
		int getPolynomIntegrationType_P3();
		int getPolynomIntegrationType_P4();
		int getHiggsDynamicsIntegrationType();
		
		
		
		int getFourierAccelerationType();
		double getFourieAccelerationParameter();
		int getPerformedPreconditioningSteps();
		int getPerformedThermalizingSteps();
		int getPerformedMeasuringSteps();		
		int getPreconditioningSteps();
		int getThermalizingSteps();
		int getMeasuringSteps();
		double getQPrecondParameterMu();
		double getQPrecondParameterBeta();
		int getSphericalHiggsIntegrationMode();
		double getZetaForHiggsIntegrationMode(); 
		double getModelParameterC6();
		double getModelParameterC8();
		double getModelParameterC10();
		double getModelParameterLambda6();
		double getModelParameterLambda8();
		double getModelParameterLambda10();		
};

#endif
