#include "StateDescriptorReader.h"
#include <cstring>

StateDescriptorReader::StateDescriptorReader(const char* stateDescriptorPath) { 
	fileName = NULL;
	fileNameIsValid = false;
			
	if(stateDescriptorPath != NULL) {
		fileName = new char[2000];	 	
		strcpy(fileName, stateDescriptorPath);				
		if (isValidFileName()) {
			fileNameIsValid = true;	 	
		}
	}	
	stateDescriptorContent = NULL;
	read();
	
}

StateDescriptorReader::~StateDescriptorReader() {
	delete[] fileName;	
	delete[] stateDescriptorContent;
	stateDescriptorContent = NULL;
	fileName = NULL;
}

void StateDescriptorReader::checkEndLineCharacter(char* str) {
	if (str != NULL && strlen(str) > 1) {
		int len = strlen(str);			
		if (str[len-1] == 0x0A && str[len-2] == 0x0D)
			str[len-2] = 0x0A;
	}	
}

bool StateDescriptorReader::isValidFileName() {
	bool res = false;
	if (fileName != NULL) {
		struct stat stat_struct;
		if (!stat(fileName, &stat_struct) && (stat_struct.st_mode & S_IFMT) == S_IFREG) {
			res = true;
		}
		else 
			printf("%s is NOT a regular file!\n", fileName);
	}
	return res;
}

double StateDescriptorReader::getDoubleParameter(int tagCount, ...) {  
	double res = 0.0;	
	if (tagCount != 0) {
		char* tmp = NULL;	
		char** tags = new char*[tagCount];
		va_start(args, tagCount);
		for (int i=0; i < tagCount; i++)  {
			tmp = va_arg(args, char*);
			tags[i] = new char[strlen(tmp)+1];
			strcpy(tags[i], tmp);
		}
		va_end(args);
		
		char* str = getStrParameter(tagCount, tags);

		int readSuccess = 0; 		
		if (str == NULL || (readSuccess = sscanf(str, "%lf", &res)) == EOF || readSuccess == 0) {
			res = 0.0;
		}
		
		delete[] str;
		for (int i = 0; i < tagCount; i++) {
			delete[] tags[i];
		}
		delete[] tags;
	}		
	return res;
}

double StateDescriptorReader::getDoubleParameter(double defaultValue, int tagCount, ...) {  
	double res = defaultValue;	
	if (tagCount != 0) {
		char* tmp = NULL;	
		char** tags = new char*[tagCount];
		va_start(args, tagCount);
		for (int i=0; i < tagCount; i++)  {
			tmp = va_arg(args, char*);
			tags[i] = new char[strlen(tmp)+1];
			strcpy(tags[i], tmp);
		}
		va_end(args);
				
		char* str = getStrParameter(tagCount, tags);
		
		int readSuccess = 0; 		
		if (str == NULL || (readSuccess = sscanf(str, "%lf", &res)) == EOF || readSuccess == 0) {
			res = defaultValue;
		}
		
		delete[] str;
		for (int i = 0; i < tagCount; i++) {
			delete[] tags[i];
		}
		delete[] tags;
	}		
	return res;
}


int StateDescriptorReader::getIntParameter(int tagCount, ...){
	int res = 0;	
	if (tagCount != 0) {
		char* tmp = NULL;	
		char** tags = new char*[tagCount];
		va_start(args, tagCount);
		for (int i=0; i < tagCount; i++)  {
			tmp = va_arg(args, char*);
			tags[i] = new char[strlen(tmp)+1];
			strcpy(tags[i], tmp);
		}
		va_end(args);
		
		char* str = getStrParameter(tagCount, tags);

		int readSuccess = 0; 		
		if (str == NULL || (readSuccess = sscanf(str, "%d", &res)) == EOF || readSuccess == 0) {
			res = 0;
		}
		
		
		delete[] str;
		for (int i = 0; i < tagCount; i++) {
			delete[] tags[i];
		}
		delete[] tags;
	}		
	return res;
}

int StateDescriptorReader::getIntParameter(int defaultValue, int tagCount, ...){
	int res = defaultValue;	
	if (tagCount != 0) {
		char* tmp = NULL;	
		char** tags = new char*[tagCount];
		va_start(args, tagCount);
		for (int i=0; i < tagCount; i++)  {
			tmp = va_arg(args, char*);
			tags[i] = new char[strlen(tmp)+1];
			strcpy(tags[i], tmp);
		}
		va_end(args);
		
		char* str = getStrParameter(tagCount, tags);

		int readSuccess = 0; 		
		if (str == NULL || (readSuccess = sscanf(str, "%d", &res)) == EOF || readSuccess == 0) {
			res = defaultValue;
		}
		
		
		delete[] str;
		for (int i = 0; i < tagCount; i++) {
			delete[] tags[i];
		}
		delete[] tags;
	}		
	return res;
}

char* StateDescriptorReader::getStrParameter(int tagCount, char** tags) {
	char* res = NULL;
	if (stateDescriptorContent != NULL) {
		bool match;		
		char * pch;
		char* match_str;
		char* stateDescriptorLines = new char[strlen(stateDescriptorContent)+1];
		strcpy(stateDescriptorLines, stateDescriptorContent);
		pch = strtok (stateDescriptorLines,"\n");
		while (pch != NULL) {		
			match = true;				
			match_str = new char[strlen(pch)+2];
			strcpy(match_str, pch);
			strcat(match_str, "\n");

			for (int i = 0; i < tagCount; i++) {				
				checkEndLineCharacter(tags[i]);
				if (strstr(match_str, tags[i]) == NULL) {
					match = false;
				}
			}
			
			delete[] match_str;
							
			if (match) {				
				res = new char[1024];	
				sprintf(res, "%s", "");				
				sscanf(pch, "%s", res);									
				break;
			}
		
			pch = strtok (NULL, "\n");
		}
		delete[] stateDescriptorLines;						
	}		
	return res;
}

void StateDescriptorReader::read() {
	if (stateDescriptorContent != NULL) {
		delete[] stateDescriptorContent;
		stateDescriptorContent = NULL;
	}
	if (fileNameIsValid) {		
		stateDescriptorContent = new char[1000000];
		sprintf(stateDescriptorContent, "%s", "");
		FILE* stateDescriptorFile = fopen(fileName, "r");
		char* str_buffer = new char[2048];
		bool isHeader = true;		
			
		while (isHeader && stateDescriptorFile != NULL && !feof(stateDescriptorFile)) {
			fgets (str_buffer, 2048, stateDescriptorFile);			
			if(strstr(str_buffer, "***") != NULL && strstr(str_buffer, "Phi") != NULL && strstr(str_buffer, "Field") != NULL) {
				isHeader = false;
			}
			else {
				checkEndLineCharacter(str_buffer);
				strcat(stateDescriptorContent, str_buffer);
			}			
		}
		fclose(stateDescriptorFile);		
		delete[] str_buffer;		
	}
		
}

char* StateDescriptorReader::getStrParameter(int tagCount, ...) {
	char* res = NULL;
	if (tagCount != 0) {
		char* tmp = NULL;	
		char** tags = new char*[tagCount];
		va_start(args, tagCount);
		for (int i=0; i < tagCount; i++)  {
			tmp = va_arg(args, char*);
			tags[i] = new char[strlen(tmp)+1];
			strcpy(tags[i], tmp);
		}
		va_end(args);
			
		res = getStrParameter(tagCount, tags);
		
		for (int i = 0; i < tagCount; i++) {
			delete[] tags[i];
		}
		delete[] tags;								
	}			
	return res;
}

char* StateDescriptorReader::getFileNameSuffix() {
  return getStrParameter(2, "Filename", "Suffix");
}

char* StateDescriptorReader::getFileNameExtension() {
  return getStrParameter(2, "Filename", "Extension");
}

int StateDescriptorReader::getL0() {
	return getIntParameter(1, ": L0\n");
}

int StateDescriptorReader::getL1() {
	return getIntParameter(1, ": L1\n");
}

int StateDescriptorReader::getL2() {
	return getIntParameter(1, ": L2\n");
}

int StateDescriptorReader::getL3() {
	return getIntParameter(1, ": L3\n");
}

int StateDescriptorReader::getNf() {
	return getIntParameter(1, ": Nf\n");
}

double StateDescriptorReader::getRho() {
	return getDoubleParameter(1, ": Rho\n");
}

double StateDescriptorReader::getR() {
	return getDoubleParameter(1, ": R\n");	
}

double StateDescriptorReader::getYN() {
	return getDoubleParameter(1, ": Y\n");
}

double StateDescriptorReader::getKappa() {
	return getDoubleParameter(1, ": Kappa\n");
}

double StateDescriptorReader::getLambda() {
	return getDoubleParameter(1, ": Lambda\n");
}

bool StateDescriptorReader::getUseP() {	
	char* fileNameSuffix = getFileNameSuffix();
	bool res = false;	
	if (strcmp(fileNameSuffix, "level8") == 0) 
		res = (bool)getIntParameter(2, "Use", "Preconditioner");
	else {
		res = (bool)getIntParameter(4, "Flag", "Use", " P ", "Preconditioner");
	}
		
	delete[] fileNameSuffix;
	return res; 
}

bool StateDescriptorReader::getUseQ() {
	char* fileNameSuffix = getFileNameSuffix();
	bool res = false;	
	if (strcmp(fileNameSuffix, "level8") == 0) 
		res = false;
	else {
		res = (bool)getIntParameter(4, "Flag", "Use", " Q ", "Preconditioner");
	}
		
	delete[] fileNameSuffix;
	return res;
}

bool StateDescriptorReader::getUseR() {
  char* fileNameSuffix = getFileNameSuffix();
	bool res = false;	
	if (strcmp(fileNameSuffix, "level8") == 0) 
		res = false;
	else {
		res = (bool)getIntParameter(4, "Flag", "Use", " R ", "Preconditioner");
	}
		
	delete[] fileNameSuffix;
	return res;
}

bool StateDescriptorReader::getUseQHM() {
	bool res = false;
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") != 0) //StateDescriptor is not level8 (so 9 or 10) 
		res = (bool)getIntParameter(4, "Flag", "Quasi", "Hermitean", "Mode");
		
	delete[] fileNameSuffix;
	return res; 
}

double StateDescriptorReader::getRPrecondParameterM() {
	double m = 0.0;
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") != 0) 
		m = getDoubleParameter(5, "Current", " R ", "Preconditioner", "parameter", " M");
		
	delete[] fileNameSuffix;
	return m; 
}

double StateDescriptorReader::getRPrecondParameterF() {
	double s = 0.0;
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") != 0) 
		s = getDoubleParameter(5, "Current", " R ", "Preconditioner", "parameter", " F");
		
	delete[] fileNameSuffix;
	return s; 
}

double StateDescriptorReader::getPPrecondParameterM() {
	double m = 0.0;
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") == 0) {
		m = getDoubleParameter(3, "Preconditioner", "parameter", " M");
	}
	else { 
		m = getDoubleParameter(5, "Current", " P ", "Preconditioner", "parameter", " M");
	}
		
	delete[] fileNameSuffix;
	return m; 
}

double StateDescriptorReader::getPPrecondParameterS() {
	double s = 0.0;
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") == 0) {
		s = getDoubleParameter(3, "Preconditioner", "parameter", " S");
	}
	else { 
		s = getDoubleParameter(5, "Current", " P ", "Preconditioner", "parameter", " S");
	}
		
	delete[] fileNameSuffix;
	return s; 
}

double StateDescriptorReader::getMassSplit() {
	double split = 1.0;
	char* fileNameSuffix = getFileNameSuffix();
	
	if (strcmp(fileNameSuffix, "level8") == 0) {
		split = 1.0;
	}
	else { 
		split = getDoubleParameter(1, "Mass-Split");
	}
		
	delete[] fileNameSuffix;
	return split;	
}

double StateDescriptorReader::getExternalCurrent() {
	double extJ = 0.0;
	extJ= getDoubleParameter(0.0, 3, " Explicit ", "Current", " J");
	return extJ;
}

double StateDescriptorReader::getAveragePhi() {
	double avgPhi = getDoubleParameter(1, " Average Phi\n");
	return avgPhi;	
}

double StateDescriptorReader::getAverageStaggeredPhi() {
	double res = getDoubleParameter(1, "Average Staggered Phi\n");
	return res;
}

double StateDescriptorReader::getExplicitFermionMass() {
	double res = 0.0;
	char* fileNameSuffix = getFileNameSuffix();
	if ((strcmp(fileNameSuffix, "level1") == 0) ||
	    (strcmp(fileNameSuffix, "level2") == 0) ||
	    (strcmp(fileNameSuffix, "level3") == 0) ||
	    (strcmp(fileNameSuffix, "level4") == 0) ||
	    (strcmp(fileNameSuffix, "level5") == 0) ||
	    (strcmp(fileNameSuffix, "level6") == 0) ||
	    (strcmp(fileNameSuffix, "level7") == 0) ||
	    (strcmp(fileNameSuffix, "level8") == 0) ||
	    (strcmp(fileNameSuffix, "level9") == 0)) {
  	      delete[] fileNameSuffix;
	      return res;
	}

	res = getDoubleParameter(2, "Explicit Fermion-Mass", " m_F");
		
	delete[] fileNameSuffix;
	return res;
}

double StateDescriptorReader::getPolynomLowerBound_P0() {
	double res = getDoubleParameter(2, "Lower Bound", "Approximation sub-Polynomial P0\n");
	return res;
}

double StateDescriptorReader::getPolynomLowerBound_P1() {
	double res = getDoubleParameter(2, "Lower Bound", "Approximation sub-Polynomial P1\n");
	return res;
}

double StateDescriptorReader::getPolynomLowerBound_P2() {
	double res = getDoubleParameter(2, "Lower Bound", "Approximation sub-Polynomial P2\n");
	return res;
}

double StateDescriptorReader::getPolynomLowerBound_P3() {
	double res = getDoubleParameter(2, "Lower Bound", "Approximation sub-Polynomial P3\n");
	return res;
}

double StateDescriptorReader::getPolynomLowerBound_P4() {
	double res = getDoubleParameter(2, "Lower Bound", "Approximation sub-Polynomial P4\n");
	return res;
}

double StateDescriptorReader::getPolynomUpperBound() {
	double res = getDoubleParameter(2, "Upper Bound", "all Approximation Polynomials\n");
	return res;
}

int StateDescriptorReader::getPolynomDegree_P0() {
	int res = getIntParameter(1, " Degree of Approximation sub-Polynomials P0\n");
	return res;
}

int StateDescriptorReader::getPolynomDegree_P1() {
	int  res = getIntParameter(1, " Degree of Approximation sub-Polynomials P1\n");
	return res;
}

int StateDescriptorReader::getPolynomDegree_P2() {
	int res = getIntParameter(1, " Degree of Approximation sub-Polynomials P2\n");
	return res;
}

int StateDescriptorReader::getPolynomDegree_P3() {
	int res = getIntParameter(1, " Degree of Approximation sub-Polynomials P3\n");
	return res;
}

int StateDescriptorReader::getPolynomDegree_P4() {
	int res = getIntParameter(1, " Degree of Approximation sub-Polynomials P4\n");
	return res;
}

double StateDescriptorReader::getAlpha() {
	double res = getDoubleParameter(2, "Exponent of function", "approximated");
	return res;
}

int StateDescriptorReader::getPolynomDigits() {
	int res = getIntParameter(1, " Digits for Calculation of Polynomial Roots");
	return res;
}

int StateDescriptorReader::getMaximumPolynomDegreePerNode() {
	int res = getIntParameter(1, " Maximum polynomial degree per node");
	return res;
}

bool StateDescriptorReader::useRandomGauge() {
	bool res = false;
	char* fileNameSuffix = getFileNameSuffix();
	
	if (strcmp(fileNameSuffix, "level8") != 0 && strcmp(fileNameSuffix, "level9") != 0) {
		res = (bool)getIntParameter(1, " Flag: Random Gauge");
	}
		
	delete[] fileNameSuffix;	
	return res;
}

bool StateDescriptorReader::useExactReweighing() {
	bool res = false;
	
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") != 0) 
		res = (bool)getIntParameter(1, " Flag: Exact Reweighing");
		
	delete[] fileNameSuffix;
	return res;
}

double StateDescriptorReader::getUpperEWBoundSafetyFactor() {
	double res = -1.0;
	
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") != 0) 
		res = getDoubleParameter(1, " Upper EW bound safety factor");
		
	delete[] fileNameSuffix;
	return res;
}

int StateDescriptorReader::getAdditionalAuxVectorCount() {
	int res = 0;	
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") != 0) 
		res = getIntParameter(1, " Number of additional auxiliary vectors");
		
	delete[] fileNameSuffix;	
	return res;
}

bool StateDescriptorReader::useDirectSamplingOfOmegaFields() {
	bool res = false;
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") != 0) 
		res = (bool)getIntParameter(1, " Flag: Direct Sampling of Omega Fields");
		
	delete[] fileNameSuffix;
	return res;
}

double StateDescriptorReader::getTheta() {
	double res = getDoubleParameter(1, " Current value of Theta");
	return res;
}

bool StateDescriptorReader::useThetaScan() {
	bool res = (bool)getIntParameter(1, " Flag - Activate Theta-Scan");
	return res;
}

double StateDescriptorReader::getThetaMin() {
	double res = getDoubleParameter(1, " Minimal Theta");
	return res;
}

double StateDescriptorReader::getThetaMax() {
	double res = getDoubleParameter(1, " Maximal Theta");
	return res;
}

double StateDescriptorReader::getThetaAdaptionFactor() {
	double res = getDoubleParameter(1, " Theta-Adaption Factor");
	return res;
}

bool StateDescriptorReader::useModelSelection() {
	bool res = false;
		
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") != 0) 
		res = (bool)getIntParameter(1, " Flag: Model Selection");
		
	delete[] fileNameSuffix;	
	return res;
}

bool StateDescriptorReader::useAntiPeriodicBoundaryConditionsInTimeDirection() {
        bool res = false;

        res = (bool)getIntParameter(0, 1, " Flag: Anti-Periodic Boundary Conditions in time direction");

        return res;
}


bool StateDescriptorReader::useXFFT() {
	bool res = (bool)getIntParameter(1, " Flag: Use of xFFT");
	return res;
}

double StateDescriptorReader::getMolecularDynamicsEpsilon() {
	double res = getDoubleParameter(1, " Epsilon for molecular dynamics");
	return res;
}

int StateDescriptorReader::getSubPolynomCount() {
	int res = getIntParameter(1, " Number of sub-Polynomials");
	return res;
}

int StateDescriptorReader::getPolynomIterations_P0() {
	int res = getIntParameter(1, " Iterations for sub-Polynomial P0");
	return res;
}

int StateDescriptorReader::getPolynomIterations_P1() {
	int res = getIntParameter(1, " Iterations for sub-Polynomial P1");
	return res;
}

int StateDescriptorReader::getPolynomIterations_P2() {
	int res = getIntParameter(1, " Iterations for sub-Polynomial P2");
	return res;
}

int StateDescriptorReader::getPolynomIterations_P3() {
	int res = getIntParameter(1, " Iterations for sub-Polynomial P3");
	return res;
}

int StateDescriptorReader::getPolynomIterations_P4() {
	int res = getIntParameter(1, " Iterations for sub-Polynomial P4");
	return res;
}

int StateDescriptorReader::getHiggsDynamicsIterations() {
	int res = getIntParameter(1, " Iterations for Higgs-Dynamics");
	return res;
}


int StateDescriptorReader::getPolynomIntegrationType_P0() {
  int res = getIntParameter(1, " Type of integrator being used for P0");
  return res;
}

int StateDescriptorReader::getPolynomIntegrationType_P1() {
  int res = getIntParameter(1, " Type of integrator being used for P1");
  return res;
}

int StateDescriptorReader::getPolynomIntegrationType_P2() {
  int res = getIntParameter(1, " Type of integrator being used for P2");
  return res;
}

int StateDescriptorReader::getPolynomIntegrationType_P3() {
  int res = getIntParameter(1, " Type of integrator being used for P3");
  return res;
}

int StateDescriptorReader::getPolynomIntegrationType_P4() {
  int res = getIntParameter(1, " Type of integrator being used for P4");
  return res;
}

int StateDescriptorReader::getHiggsDynamicsIntegrationType() {
  int res = getIntParameter(1, " Type of integrator being used for Higgs-Dynamics");
  return res;
}


int StateDescriptorReader::getFourierAccelerationType() {
	int res = getIntParameter(1, " Fourier Acceleration type");
	return res;
}

double StateDescriptorReader::getFourieAccelerationParameter() {
	double res = getDoubleParameter(1, " Fourier Acceleration parameter");
	return res;
}

int StateDescriptorReader::getPerformedPreconditioningSteps() {
	int res = getIntParameter(1, " Performed Preconditioning successful Metros");
	return res; 	
}

int StateDescriptorReader::getPerformedThermalizingSteps() {
	int res = getIntParameter(1, " Performed Thermalizing successful Metros");
	return res;
}

int StateDescriptorReader::getPerformedMeasuringSteps() {
	int res = getIntParameter(1, " Performed Measuring successful Metros Steps");
	return res;
}
		
int StateDescriptorReader::getPreconditioningSteps() {
	int res = getIntParameter(1, "Number of Metropolis steps for automatic Preconditioning");
	return res;
}

int StateDescriptorReader::getThermalizingSteps() {
	int res = getIntParameter(1, "Number of Metropolis steps for thermalizing");
	return res;
}

int StateDescriptorReader::getMeasuringSteps() {
	int res = getIntParameter(1, " Number of total data to be collected");
	return res;
}

double StateDescriptorReader::getQPrecondParameterMu() {
	double res = 0.0;
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") != 0) 
		res = getDoubleParameter(2, "Current Q", "Preconditioner parameter Mu");
		
	delete[] fileNameSuffix;
	return res;	
}

double StateDescriptorReader::getQPrecondParameterBeta() {
	double res = 0.0;
	char* fileNameSuffix = getFileNameSuffix();
	if (strcmp(fileNameSuffix, "level8") != 0) 
		res = getDoubleParameter(2, "Current Q", "Preconditioner parameter Beta");
		
	delete[] fileNameSuffix;
	return res; 	
}

int StateDescriptorReader::getSphericalHiggsIntegrationMode() {
	int res = getIntParameter(0, 1, "Spherical Higgs-Integration Mode");
	return res; 	
}

double StateDescriptorReader::getZetaForHiggsIntegrationMode() {
	double res = getDoubleParameter(0.0, 1, "Zeta-value for Higgs-Integration Mode");
	return res;
}

double StateDescriptorReader::getModelParameterC6() {
	double res = getDoubleParameter(0.0, 1, "Model-Extension Parameter c6");
	return res;
}

double StateDescriptorReader::getModelParameterC8() {
	double res = getDoubleParameter(0.0, 1, "Model-Extension Parameter c8");
	return res;
}

double StateDescriptorReader::getModelParameterC10() {
	double res = getDoubleParameter(0.0, 1, "Model-Extension Parameter c10");
	return res;
}

double StateDescriptorReader::getModelParameterLambda6() {
	double res = getDoubleParameter(0.0, 1, "Model-Extension Parameter lambda6");
	return res;
}

double StateDescriptorReader::getModelParameterLambda8() {
	double res = getDoubleParameter(0.0, 1, "Model-Extension Parameter lambda8");
	return res;
}

double StateDescriptorReader::getModelParameterLambda10() {
	double res = getDoubleParameter(0.0, 1, "Model-Extension Parameter lambda10");
	return res;
}
