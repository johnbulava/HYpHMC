#include <cstring>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <unistd.h>

#ifdef UseMPI
  #include <mpi.h>
#endif 


#include "Complex.h"
#include "Quat.h"
#include "Global.h"
#include "ComplexVector.h"
#include "ComplexMatrix.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"
#include "pHMCPropagator.h"

#define AutomaticAdaption_HighPro 0.94
#define AutomaticAdaption_LowPro 0.80
#define AutomaticAdaption_LowPro2 0.10

#define AutomaticIntegratorAdaption_MaxLeapFrogInversions 10
#define ThetaIterIntegratorTableMAX 100000
#define SubPolynomMAX 4
#define MinimalIterationCountForP0BeforeMassDet 3
#define IntegrationSchemeAdaption


//Variables
int Parameter_L0 = 0;
int Parameter_L1 = 0;
int Parameter_L2 = 0;
int Parameter_L3 = 0;
int Parameter_Nf = 0;
double Parameter_RHO = 0;
double Parameter_R = 0;
double Parameter_Lambda = 0;
double Parameter_Y = 0;
double Parameter_Kappa = 0;
double Parameter_MassSplit = 0;
double Parameter_ExplicitMass = 0;
double Parameter_ExplicitCurrent = 0;
double Parameter_Epsilon = 0;
int Parameter_SubPolyCount = 0;
int Parameter_PolyDegree[1+SubPolynomMAX];
double Parameter_PolyEpsilon[1+SubPolynomMAX];
double Parameter_PolyLambda = 0;
double Parameter_PolyAlpha = 0;
int Parameter_PolDigit = 0;
double Parameter_Theta = 0;
double Parameter_ThetaMin = 0;
double Parameter_ThetaMax = 0;
double Parameter_ThetaFactor = 0;
int Parameter_FLAG_ThetaScan = 0;
int Parameter_AdditionalAuxVectors = 0;
int Parameter_FLAG_QuasiHermiteanMode = 0;
int Parameter_FLAG_FactorizationOfMatrixB = 0;
int Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime = 0;
int Parameter_FLAG_Use_P_Preconditioner = 0;
int Parameter_FLAG_Use_Q_Preconditioner = 0;
int Parameter_FLAG_Use_R_Preconditioner = 0;
int Parameter_FLAG_RandomGauge = 0;
double Parameter_QPreconditioner_Mu = 0;
double Parameter_QPreconditioner_Beta = 0;
double Parameter_UpperEWboundSafetyFactor = 0;
int Parameter_FLAG_DirectOmegaSampling = 0;
int Parameter_MaxPolDegPerNode = 0;
int Parameter_Iterations[2+SubPolynomMAX];
int Parameter_Measurements = 0;
int Parameter_AutomaticPreconditioningMetros = 0;
int Parameter_ThermalizingMetros = 0;
int Parameter_TotalData = 0;
int Parameter_SaveConfigurationEveryXXXResult = 0;
int Parameter_SaveStateDescriptorEveryXXXResult = 0;
int Parameter_DetermineEWeveryXXXconfsPREC = 0;
int Parameter_DetermineEWeveryXXXconfsMEASURE = 0;
char* Parameter_filenameSuffix;
int Parameter_FLAG_ModelSelection = 0;
int Parameter_FLAG_xFFT = 0;
int Parameter_FACC_type = 0;
int Parameter_OmegaMassAdaptionMode = 0;
double Parameter_FACC_parameter = 0;
int Parameter_StartRandSeed = 0;
int Parameter_ReversibilityCheckFreqPrec = 0;
int Parameter_ReversibilityCheckFreqMeas = 0;
double Parameter_MaxRunTime = 0;
int Parameter_FLAG_ExactReweighing = 0;
int Parameter_IntegratorType[2+SubPolynomMAX];
int Parameter_SphericalHiggsMode = 0;
double Parameter_SphericalHiggsZeta = 0;
int Parameter_MultiThreadedOpMode = 0;
int Parameter_FFTWThreadCount = 0;
int Parameter_xFFTThreadCount = 0;
double Parameter_c6 = 0;
double Parameter_c8 = 0;
double Parameter_c10 = 0;
double Parameter_lambda6 = 0;
double Parameter_lambda8 = 0;
double Parameter_lambda10 = 0;
FermionMatrixOperations* fermiOps = NULL;
pHMCPropagator* pHMCProp = NULL;
int ThermalizingRunCount = 0;
int AutomaticPreconRunCount = 0;
int TotallyMeasuredConfigurationsCount = 0;
char* outputFileNameExtension = NULL;
double startTime;
int stateDescriptorCounter = 0;
int phiConfCounter = -1;
int phiConfIDnr = 0;
int JobNr = -1;
double** ThetaIterIntegratorTable = new double*[ThetaIterIntegratorTableMAX];
int ThetaIterIntegratorTableCount = 0;


void startTimer() {
  startTime = zeitwert();
}


double timePassed() {
  double time = (zeitwert()-startTime);
  return time;
}


bool timeOver() {
  if (timePassed()>Parameter_MaxRunTime*3600) return true;
  return false;
}


void loadParameters(int MultiProcessNr, bool ll) {
  FILE* file;
  if (LogLevel>0) printf("Trying to read parameter-file 'SimulationParameterspHMC.txt'...\n");
  file = fopen("SimulationParameterspHMC.txt","r");
  char *Comment = new char[500];
  bool error = false;
  Parameter_filenameSuffix = new char[500];
  int I;

  if (fscanf(file,"%s",Parameter_filenameSuffix)==1) {
    if (LogLevel>0) printf(" -> Filename-Suffix: %s\n",Parameter_filenameSuffix);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_L0)==1) {
    if (LogLevel>0) printf(" -> L0: %d\n",Parameter_L0);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_L1)==1) {
    if (LogLevel>0) printf(" -> L1: %d\n",Parameter_L1);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_L2)==1) {
    if (LogLevel>0) printf(" -> L2: %d\n",Parameter_L2);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_L3)==1) {
    if (LogLevel>0) printf(" -> L3: %d\n",Parameter_L3);
  } else error = true;
  fgets(Comment, 500, file);  
  
  if (fscanf(file,"%d",&Parameter_Nf)==1) {
    if (LogLevel>0) printf(" -> Nf: %d\n",Parameter_Nf);
  } else error = true;
  fgets(Comment, 500, file);
    
  if (fscanf(file,"%lf",&Parameter_RHO)==1) {
    if (LogLevel>0) printf(" -> Rho: %f\n",Parameter_RHO);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%lf",&Parameter_R)==1) {
    if (LogLevel>0) printf(" -> R: %f\n",Parameter_R);
  } else error = true;
  fgets(Comment, 500, file);
  
    
  double lambdaMin, lambdaMax;
  int lambdaTics;
  if (fscanf(file,"%lf",&lambdaMin)==1) {
    if (LogLevel>0) printf(" -> lambdaMin: %f\n",lambdaMin);
  } else error = true;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf",&lambdaMax)==1) {
    if (LogLevel>0) printf(" -> lambdaMax: %f\n",lambdaMax);
  } else error = true;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d",&lambdaTics)==1) {
    if (LogLevel>0) printf(" -> lambdaTics: %d\n",lambdaTics);
  } else error = true;
  fgets(Comment, 500, file);
  
  double yMin, yMax;
  int yTics;
  if (fscanf(file,"%lf",&yMin)==1) {
    if (LogLevel>0) printf(" -> yMin: %f\n",yMin);
  } else error = true;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf",&yMax)==1) {
    if (LogLevel>0) printf(" -> yMax: %f\n",yMax);
  } else error = true;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d",&yTics)==1) {
    if (LogLevel>0) printf(" -> yTics: %d\n",yTics);
  } else error = true;
  fgets(Comment, 500, file);
  
  double kappaMin, kappaMax;
  int kappaTics;
  if (fscanf(file,"%lf",&kappaMin)==1) {
    if (LogLevel>0) printf(" -> kappaMin: %f\n",kappaMin);
  } else error = true;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf",&kappaMax)==1) {
    if (LogLevel>0) printf(" -> kappaMax: %f\n",kappaMax);
  } else error = true;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d",&kappaTics)==1) {
    if (LogLevel>0) printf(" -> kappaMTics: %d\n",kappaTics);
  } else error = true;
  fgets(Comment, 500, file);
  
  MultiProcessNr = MultiProcessNr % (lambdaTics*yTics*kappaTics);
 
  if (kappaTics <= 1) {
    Parameter_Kappa = kappaMin;
  } else {
    Parameter_Kappa = kappaMin+(kappaMax-kappaMin)*(MultiProcessNr % kappaTics)/(kappaTics-1);
    MultiProcessNr /= kappaTics;
  }
  if (yTics <= 1) {
    Parameter_Y = yMin;
  } else {
    Parameter_Y = yMin+(yMax-yMin)*(MultiProcessNr % yTics)/(yTics-1);
    MultiProcessNr /= yTics;    
  }
  if (lambdaTics <= 1) {
    Parameter_Lambda = lambdaMin;
  } else {
    Parameter_Lambda = lambdaMin+(lambdaMax-lambdaMin)*(MultiProcessNr % lambdaTics)/(lambdaTics-1);
    MultiProcessNr /= lambdaTics;    
  }
  
  if (LogLevel>0) printf(" -> Lambda (selected): %f (infinity=%d)\n",Parameter_Lambda, isNaN(Parameter_Lambda));  
  if (LogLevel>0) printf(" -> Y (selected): %f\n",Parameter_Y);
  if (LogLevel>0) printf(" -> Kappa (selected): %f\n",Parameter_Kappa);

  if (fscanf(file,"%lf",&Parameter_MassSplit)==1) {
    if (LogLevel>0) printf(" -> Mass-Split y_b / y_t (is one: %d): %f\n", Parameter_MassSplit==1, Parameter_MassSplit);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_ExplicitMass)==1) {
    if (LogLevel>0) printf(" -> Explicit Fermion-Mass m_F (is zero: %d): %f\n", Parameter_ExplicitMass==0, Parameter_ExplicitMass);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_ExplicitCurrent)==1) {
    if (LogLevel>0) printf(" -> Explicit Current J (is zero: %d): %f\n", Parameter_ExplicitCurrent==0, Parameter_ExplicitCurrent);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_c6)==1) {
    if (LogLevel>0) printf(" -> Model-Extension Parameter c6 : %f\n",Parameter_c6);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_c8)==1) {
    if (LogLevel>0) printf(" -> Model-Extension Parameter c8 : %f\n",Parameter_c8);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_c10)==1) {
    if (LogLevel>0) printf(" -> Model-Extension Parameter c10 : %f\n",Parameter_c10);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_lambda6)==1) {
    if (LogLevel>0) printf(" -> Model-Extension Parameter lambda6 : %f\n",Parameter_lambda6);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_lambda8)==1) {
    if (LogLevel>0) printf(" -> Model-Extension Parameter lambda8 : %f\n",Parameter_lambda8);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_lambda10)==1) {
    if (LogLevel>0) printf(" -> Model-Extension Parameter lambda10 : %f\n",Parameter_lambda10);
  } else error = true;
  fgets(Comment, 500, file);

  for (I=0; I<1+SubPolynomMAX; I++) {
    if (fscanf(file,"%lf",&(Parameter_PolyEpsilon[I]))==1) {
      if (LogLevel>0) printf(" -> Epsilon of Approximation sub-Polynomial P%d: %f\n",I, Parameter_PolyEpsilon[I]);
    } else error = true;
    fgets(Comment, 500, file);
  }

  if (fscanf(file,"%lf",&Parameter_PolyLambda)==1) {
    if (LogLevel>0) printf(" -> Lambda for all Approximation Polynomials: %f\n",Parameter_PolyLambda);
  } else error = true;
  fgets(Comment, 500, file);

  for (I=0; I<1+SubPolynomMAX; I++) {
    if (fscanf(file,"%d",&(Parameter_PolyDegree[I]))==1) {
      if (LogLevel>0) printf(" -> Degree of Approximation sub-Polynomial P%d: %d\n",I, Parameter_PolyDegree[I]);
    } else error = true;
    fgets(Comment, 500, file);
  }

  if (fscanf(file,"%lf",&Parameter_PolyAlpha)==1) {
    if (LogLevel>0) printf(" -> Alpha of Approximation Polynomial: %f\n",Parameter_PolyAlpha);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_PolDigit)==1) {
    if (LogLevel>0) printf(" -> Digits for calculation of polynomial roots: %d\n",Parameter_PolDigit);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_MaxPolDegPerNode)==1) {
    if (LogLevel>0) printf(" -> Maximal Polynomial Degree on each node: %d\n",Parameter_MaxPolDegPerNode);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_RandomGauge)==1) {
    if (LogLevel>0) printf(" -> Random Gauge (0: deactivated): %d\n",Parameter_FLAG_RandomGauge);
  } else error = true;
  fgets(Comment, 500, file);
  if (Parameter_FLAG_RandomGauge!=0) {
    printf("Random Gauge not implemented yet!!!\n");
    exit(0);
  }
  
  if (fscanf(file,"%d",&Parameter_FLAG_ExactReweighing)==1) {
    if (LogLevel>0) printf(" -> FLAG: Exact reweighing: %d\n",Parameter_FLAG_ExactReweighing);
  } else error = true;
  fgets(Comment, 500, file);
   
  if (fscanf(file,"%d",&Parameter_FLAG_QuasiHermiteanMode)==1) {
    if (LogLevel>0) printf(" -> FLAG: Quasi-Hermitean Mode: %d\n",Parameter_FLAG_QuasiHermiteanMode);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_FactorizationOfMatrixB)==1) {
    if (LogLevel>0) printf(" -> FLAG: Parameter_FLAG_FactorizationOfMatrixB: %d\n",Parameter_FLAG_FactorizationOfMatrixB);
  } else {
		Parameter_FLAG_FactorizationOfMatrixB = 0; 
		if (LogLevel>0) printf(" -> FLAG: Parameter_FLAG_FactorizationOfMatrixB: %d\n",Parameter_FLAG_FactorizationOfMatrixB);
	}
		
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime)==1) {
    if (LogLevel>0) printf(" -> FLAG: Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime: %d\n",Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime);
  } else {
		Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime = 0;
		if (LogLevel>0) printf(" -> FLAG: Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime: %d\n",Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime);
	}

  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_Use_P_Preconditioner)==1) {
    if (LogLevel>0) printf(" -> FLAG: Use P-Preconditioner: %d\n",Parameter_FLAG_Use_P_Preconditioner);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_FLAG_Use_Q_Preconditioner)==1) {
    if (LogLevel>0) printf(" -> FLAG: Use Q-Preconditioner: %d\n",Parameter_FLAG_Use_Q_Preconditioner);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_Use_R_Preconditioner)==1) {
    if (LogLevel>0) printf(" -> FLAG: Use R-Preconditioner: %d\n",Parameter_FLAG_Use_R_Preconditioner);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%lf",&Parameter_QPreconditioner_Mu)==1) {
    if (LogLevel>0) printf(" -> Q-Preconditioner mu: %f\n",Parameter_QPreconditioner_Mu);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_QPreconditioner_Beta)==1) {
    if (LogLevel>0) printf(" -> Q-Preconditioner beta: %f\n",Parameter_QPreconditioner_Beta);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_UpperEWboundSafetyFactor)==1) {
    if (LogLevel>0) printf(" -> Upper EW bound safety factor: %f\n",Parameter_UpperEWboundSafetyFactor);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_AdditionalAuxVectors)==1) {
    if (LogLevel>0) printf(" -> Number of additional auxiliary vectors: %d\n",Parameter_AdditionalAuxVectors);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_DirectOmegaSampling)==1) {
    if (LogLevel>0) printf(" -> FLAG: Sample omega fields directly: %d\n",Parameter_FLAG_DirectOmegaSampling);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_Theta)==1) {
    if (LogLevel>0) printf(" -> Start value of Theta: %f",Parameter_Theta);
    if ((LogLevel>0) && (Parameter_Theta == 0)) printf(" (zero)\n");
    if ((LogLevel>0) && (Parameter_Theta != 0)) printf(" (non-zero)\n");
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_ThetaScan)==1) {
    if (LogLevel>0) printf(" -> Flag - Activate Theta-Scan: %d\n",Parameter_FLAG_ThetaScan);
  } else error = true;
  fgets(Comment, 500, file); 
 
  if (fscanf(file,"%lf",&Parameter_ThetaMin)==1) {
    if (LogLevel>0) printf(" -> Minimal Theta: %f\n",Parameter_ThetaMin);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_ThetaMax)==1) {
    if (LogLevel>0) printf(" -> Maximal Theta: %f\n",Parameter_ThetaMax);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_ThetaFactor)==1) {
    if (LogLevel>0) printf(" -> Theta-Adaption Factor: %f\n",Parameter_ThetaFactor);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_FLAG_ModelSelection)==1) {
    if (LogLevel>0) printf(" -> FLAG: Model Selection (0: Wilson-Model, 1: Neuberger - Model): %d\n",Parameter_FLAG_ModelSelection);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_xFFT)==1) {
    if (LogLevel>0) printf(" -> Flag - Use of xFFT: %d\n",Parameter_FLAG_xFFT);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_Epsilon)==1) {
    if (LogLevel>0) printf(" -> Epsilon for molecular dynamics: %f\n",Parameter_Epsilon);
  } else error = true;
  fgets(Comment, 500, file);

  for (I=0; I<1+SubPolynomMAX; I++) {
    Parameter_Iterations[I] = 0;
    if (fscanf(file,"%d",&(Parameter_Iterations[I]))==1) {
      if (LogLevel>0) printf(" -> Iterations for sub-Polynom P%d: %d\n",I,Parameter_Iterations[I]);
    } else error = true;
    fgets(Comment, 500, file);
  }

  Parameter_SubPolyCount = SubPolynomMAX;
  for (I=SubPolynomMAX; I>=0; I--) {
    if (Parameter_Iterations[I]==0) Parameter_SubPolyCount = I-1;
  }
  if (LogLevel>0) printf(" -> Number of sub-Polynomials: %d\n",Parameter_SubPolyCount);
  if (Parameter_SubPolyCount<0) {
    printf("ERROR: Number of sub-polynomials must not be negative!!!\n");
    exit(0);
  }

  if (fscanf(file,"%d",&(Parameter_Iterations[1+SubPolynomMAX]))==1) {
    if (LogLevel>0) printf(" -> Iterations for Higgs-Dynamics: %d\n",Parameter_Iterations[1+SubPolynomMAX]);
  } else error = true;
  fgets(Comment, 500, file);

  for (I=0; I<1+SubPolynomMAX; I++) {
    if (fscanf(file,"%d",&(Parameter_IntegratorType[I]))==1) {
      if (LogLevel>0) printf(" -> Type of integrator being used for sub-Polynom P%d: %d\n",I,Parameter_IntegratorType[I]);
    } else error = true;
    fgets(Comment, 500, file);
  }

  if (fscanf(file,"%d",&(Parameter_IntegratorType[1+SubPolynomMAX]))==1) {
    if (LogLevel>0) printf(" -> Type of integrator being used Higgs-Dynamics: %d\n",Parameter_IntegratorType[1+SubPolynomMAX]);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_Measurements)==1) {
    if (LogLevel>0) printf(" -> Measurement attempts before integrator adaption: %d\n",Parameter_Measurements);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_AutomaticPreconditioningMetros)==1) {
    if (LogLevel>0) printf(" -> Number of Metropolis steps for automatic Preconditioning: %d\n",Parameter_AutomaticPreconditioningMetros);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_ThermalizingMetros)==1) {
    if (LogLevel>0) printf(" -> Number of Metropolis steps for thermalizing: %d\n",Parameter_ThermalizingMetros);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_TotalData)==1) {
    if (LogLevel>0) printf(" -> Number of total data to be collected: %d\n",Parameter_TotalData);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_SaveConfigurationEveryXXXResult)==1) {
    if (LogLevel>0) printf(" -> Save Phi field every xxx results: %d\n",Parameter_SaveConfigurationEveryXXXResult);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_SaveStateDescriptorEveryXXXResult)==1) {
    if (LogLevel>0) printf(" -> Save State-Descriptor field every xxx results: %d\n",Parameter_SaveStateDescriptorEveryXXXResult);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_DetermineEWeveryXXXconfsPREC)==1) {
    if (LogLevel>0) printf(" -> Determine Eigenvalues every xxx configurations (PREC-phase): %d\n",Parameter_DetermineEWeveryXXXconfsPREC);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_DetermineEWeveryXXXconfsMEASURE)==1) {
    if (LogLevel>0) printf(" -> Determine Eigenvalues every xxx configurations (MEASURE-phase): %d\n",Parameter_DetermineEWeveryXXXconfsMEASURE);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_FACC_type)==1) {
    if (LogLevel>0) printf(" -> Fourier Acceleration Type: %d\n",Parameter_FACC_type);
  } else error = true;
  fgets(Comment, 500, file);    

  if (fscanf(file,"%lf",&Parameter_FACC_parameter)==1) {
    if (LogLevel>0) printf(" -> Fourier Acceleration parameter: %f\n",Parameter_FACC_parameter);
  } else error = true;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d",&Parameter_SphericalHiggsMode)==1) {
    if (LogLevel>0) printf(" -> Spherical Higgs-Integration Mode: %d\n",Parameter_SphericalHiggsMode);
  } else error = true;
  fgets(Comment, 500, file);    

  if (fscanf(file,"%lf",&Parameter_SphericalHiggsZeta)==1) {
    if (LogLevel>0) printf(" -> Zeta-value for Higgs-Integration Mode: %1.3e\n",Parameter_SphericalHiggsZeta);
  } else error = true;
  fgets(Comment, 500, file);
    
  if (fscanf(file,"%d",&Parameter_OmegaMassAdaptionMode)==1) {
    if (LogLevel>0) printf(" -> Omega Mass Adaption Mode: %d\n",Parameter_OmegaMassAdaptionMode);
  } else error = true;
  fgets(Comment, 500, file);    
  
  if (fscanf(file,"%d",&Parameter_ReversibilityCheckFreqPrec)==1) {
    if (LogLevel>0) printf(" -> Reversibility check frequency, PREC-phase: %d\n",Parameter_ReversibilityCheckFreqPrec);
  } else error = true;
  fgets(Comment, 500, file);    
  
  if (fscanf(file,"%d",&Parameter_ReversibilityCheckFreqMeas)==1) {
    if (LogLevel>0) printf(" -> Reversibility check frequency, MEASURE-phase: %d\n",Parameter_ReversibilityCheckFreqMeas);
  } else error = true;
  fgets(Comment, 500, file);    
  
  if (fscanf(file,"%d",&Parameter_StartRandSeed)==1) {
    if (LogLevel>0) printf(" -> Start value for RandSeed: %d\n",Parameter_StartRandSeed);
  } else error = true;
  fgets(Comment, 500, file);    

  if (fscanf(file,"%d",&Parameter_MultiThreadedOpMode)==1) {
    if (LogLevel>0) printf(" -> Multi-Threaded Operation Mode: %d\n",Parameter_MultiThreadedOpMode);
  } else error = true;
  fgets(Comment, 500, file);    

  if (fscanf(file,"%d",&Parameter_FFTWThreadCount)==1) {
    if (LogLevel>0) printf(" -> Number of threads used for each FFTW Transformation: %d\n",Parameter_FFTWThreadCount);
  } else error = true;
  fgets(Comment, 500, file);    

  if (fscanf(file,"%d",&Parameter_xFFTThreadCount)==1) {
    if (LogLevel>0) printf(" -> Number of threads used for each xFFT Transformation: %d\n",Parameter_xFFTThreadCount);
  } else error = true;
  fgets(Comment, 500, file);    

  if (fscanf(file,"%lf",&Parameter_MaxRunTime)==1) {
    if (LogLevel>0) printf(" -> Maximal RunTime in hours: %1.15f\n",Parameter_MaxRunTime);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&LogLevel)==1) {
    if ((LogLevel>0) && (ll)) printf(" -> Log - Level: %d\n",LogLevel);
  } else error = true;
  fgets(Comment, 500, file);

  fclose(file);
  if (error) {
    printf("... ERROR in paramters file!\n");
    exit(0);
  }  
  delete [] Comment;
  if ((LogLevel>0) && (ll)) printf("...sucessfully.\n");
}


void writeCurrentStateDescriptor(int running, bool forcedWrite) {  
  if (Parameter_SaveStateDescriptorEveryXXXResult == 0) {
    return;
  }
  if (!forcedWrite) {
    stateDescriptorCounter++;
    if (stateDescriptorCounter < Parameter_SaveStateDescriptorEveryXXXResult) {
      return;
    }
    stateDescriptorCounter = 0;
  }
  
  if (LogLevel>2) printf("Writing StateDescriptor to disk...\n");
  double measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm;
  pHMCProp->measure(measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm);
  
  bool usePrec;
  double precM, precS;
  fermiOps->getPreconditionerParameter(usePrec, precM, precS);
  
  bool useQPrec;
  double QprecMu, QprecBeta;
  fermiOps->getQPreconditionerParameter(useQPrec, QprecMu, QprecBeta);
  
  bool useRPrec;
  double RprecM, RprecF;
  fermiOps->getRPreconditionerParameter(useRPrec, RprecM, RprecF);
  
  int I;
  char* fileName = new char[600]; 
  snprintf(fileName,600,"%s/data/results/pHMC/states/StateDescriptor%s.dat",DataBaseDirectory,outputFileNameExtension);
  FILE* file;
  file = fopen(fileName,"w");

  fprintf(file,"%d                                                      : Running\n", running);
  fprintf(file,"%1.15f                                      : Average Phi\n", measurePhiNorm);
  fprintf(file,"%1.15f                                      : Average Staggered Phi\n", measureStaggeredPhiNorm);
  fprintf(file,"%1.15f                                      : Average Vector Length of Phi\n", avgNorm);
  fprintf(file,"%1.15f                                      : Sigma of vector Length of Phi\n", sigmaNorm);  
  fprintf(file,"%d                                                      : Performed Preconditioning successful Metros\n", AutomaticPreconRunCount);
  fprintf(file,"%d                                                      : Performed Thermalizing successful Metros\n", ThermalizingRunCount);
  fprintf(file,"%d                                                      : Performed Measuring successful Metros Steps\n", TotallyMeasuredConfigurationsCount);
  fprintf(file,"%d                                                      : Current P - Preconditioner usage\n", usePrec);  
  fprintf(file,"%1.15f                                      : Current P - Preconditioner parameter M\n", precM);
  fprintf(file,"%1.15f                                      : Current P - Preconditioner parameter S\n", precS);  
  fprintf(file,"%d                                                      : Current Q - Preconditioner usage\n", useQPrec);
  fprintf(file,"%1.15f                                      : Current Q - Preconditioner parameter Mu\n", QprecMu);
  fprintf(file,"%1.15f                                      : Current Q - Preconditioner parameter Beta\n", QprecBeta);  
  fprintf(file,"%d                                                      : Current R - Preconditioner usage\n", useRPrec);
  fprintf(file,"%1.15f                                      : Current R - Preconditioner parameter M\n", RprecM);
  fprintf(file,"%1.15f                                      : Current R - Preconditioner parameter F\n", RprecF);
  fprintf(file,"%s                                                 : Filename-Suffix\n", Parameter_filenameSuffix);
  fprintf(file,"%s          : Filename-Extension\n", outputFileNameExtension);
  fprintf(file,"%d                                                      : L0\n", Parameter_L0);
  fprintf(file,"%d                                                      : L1\n", Parameter_L1);
  fprintf(file,"%d                                                      : L2\n", Parameter_L2);
  fprintf(file,"%d                                                      : L3\n", Parameter_L3);
  fprintf(file,"%d                                                      : Nf\n", Parameter_Nf);
  fprintf(file,"%1.15f                                      : Rho\n", Parameter_RHO);
  fprintf(file,"%1.15f                                      : R\n", Parameter_R);
  fprintf(file,"%1.15f                                      : Lambda\n", Parameter_Lambda);
  fprintf(file,"%1.15f                                      : Kappa\n", Parameter_Kappa );
  fprintf(file,"%1.15f                                      : Y\n", Parameter_Y);
  fprintf(file,"%1.15f                                      : Mass-Split y_b / y_t\n", Parameter_MassSplit);
  fprintf(file,"%1.15f                                      : Explicit Fermion-Mass m_F\n", Parameter_ExplicitMass);
  fprintf(file,"%1.15f                                      : Explicit Current J\n", Parameter_ExplicitCurrent);      
  fprintf(file,"%1.15f                                      : Model-Extension Parameter c6\n", Parameter_c6);  
  fprintf(file,"%1.15f                                      : Model-Extension Parameter c8\n", Parameter_c8);  
  fprintf(file,"%1.15f                                      : Model-Extension Parameter c10\n", Parameter_c10);  
  fprintf(file,"%1.15f                                      : Model-Extension Parameter lambda6\n", Parameter_lambda6);  
  fprintf(file,"%1.15f                                      : Model-Extension Parameter lambda8\n", Parameter_lambda8);  
  fprintf(file,"%1.15f                                      : Model-Extension Parameter lambda10\n", Parameter_lambda10);  
  for (I=0; I<1+SubPolynomMAX; I++) {
    fprintf(file,"%1.15f                                      : Lower Bound of Approximation sub-Polynomial P%d\n", Parameter_PolyEpsilon[I], I);
  }
  fprintf(file,"%1.15f                                      : Upper Bound for all Approximation Polynomials\n", Parameter_PolyLambda);
  for (I=0; I<1+SubPolynomMAX; I++) {
    fprintf(file,"%d                                                      : Degree of Approximation sub-Polynomials P%d\n", Parameter_PolyDegree[I], I);
  }
  fprintf(file,"%1.15f                                      : Exponent of function to be approximated: f(x) = 1/x^alpha\n", Parameter_PolyAlpha);
  fprintf(file,"%d                                                      : Digits for Calculation of Polynomial Roots\n", Parameter_PolDigit);
  fprintf(file,"%d                                                      : Maximum polynomial degree per node\n", Parameter_MaxPolDegPerNode);  
  fprintf(file,"%d                                                      : Flag: Random Gauge\n", Parameter_FLAG_RandomGauge);
  fprintf(file,"%d                                                      : Flag: Exact Reweighing\n", Parameter_FLAG_ExactReweighing);
  fprintf(file,"%d                                                      : Flag: Quasi-Hermitean Mode\n", Parameter_FLAG_QuasiHermiteanMode);
  fprintf(file,"%d                                                      : Flag: Factorization of matrix B\n", Parameter_FLAG_FactorizationOfMatrixB);
  fprintf(file,"%d                                                      : Flag: Anti-Periodic Boundary Conditions in time direction\n", Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime);
  fprintf(file,"%d                                                      : Flag: Use P - Preconditioner\n", Parameter_FLAG_Use_P_Preconditioner);
  fprintf(file,"%d                                                      : Flag: Use Q - Preconditioner\n", Parameter_FLAG_Use_Q_Preconditioner);
  fprintf(file,"%d                                                      : Flag: Use R - Preconditioner\n", Parameter_FLAG_Use_R_Preconditioner);
  fprintf(file,"%1.15f                                      : Q-Preconditioning: Mu\n", Parameter_QPreconditioner_Mu);
  fprintf(file,"%1.15f                                      : Q-Preconditioning: Beta\n", Parameter_QPreconditioner_Beta);  
  fprintf(file,"%1.15f                                      : Upper EW bound safety factor\n", Parameter_UpperEWboundSafetyFactor);
  fprintf(file,"%d                                                      : Number of additional auxiliary vectors\n", Parameter_AdditionalAuxVectors);  
  fprintf(file,"%d                                                      : Flag: Direct Sampling of Omega Fields\n", Parameter_FLAG_DirectOmegaSampling);
  fprintf(file,"%1.15f                                      : Current value of Theta\n", Parameter_Theta);
  fprintf(file,"%d                                                      : Flag - Activate Theta-Scan\n", Parameter_FLAG_ThetaScan);  
  fprintf(file,"%1.15f                                      : Minimal Theta\n", Parameter_ThetaMin);
  fprintf(file,"%1.15f                                      : Maximal Theta\n", Parameter_ThetaMax);
  fprintf(file,"%1.15f                                      : Theta-Adaption Factor\n", Parameter_ThetaFactor);
  fprintf(file,"%d                                                      : Flag: Model Selection\n", Parameter_FLAG_ModelSelection);
  fprintf(file,"%d                                                      : Flag: Use of xFFT\n", Parameter_FLAG_xFFT);
  fprintf(file,"%1.15f                                      : Epsilon for molecular dynamics \n",Parameter_Epsilon );


  for (I=0; I<1+SubPolynomMAX; I++) {
    fprintf(file,"%d                                                      : Type of integrator being used for P%d\n", Parameter_IntegratorType[I], I);
  }
  fprintf(file,"%d                                                      : Type of integrator being used for Higgs-Dynamics\n", Parameter_IntegratorType[1+SubPolynomMAX]);
  fprintf(file,"%d                                                     : Number of sub-Polynomials\n", Parameter_SubPolyCount);
  for (I=0; I<1+SubPolynomMAX; I++) {
    fprintf(file,"%d                                                     : Iterations for sub-Polynomial P%d\n", Parameter_Iterations[I], I);
  }
  fprintf(file,"%d                                                     : Iterations for Higgs-Dynamics\n", Parameter_Iterations[1+SubPolynomMAX]);
  fprintf(file,"%d                                                     : Measurement attempts before adaption of integrator\n", Parameter_Measurements);
  fprintf(file,"%d                                                    : Number of Metropolis steps for automatic Preconditioning \n", Parameter_AutomaticPreconditioningMetros);
  fprintf(file,"%d                                                    : Number of Metropolis steps for thermalizing \n", Parameter_ThermalizingMetros);
  fprintf(file,"%d                                                  : Number of total data to be collected \n", Parameter_TotalData);
  fprintf(file,"%d                                                      : Save Phi field every xxx results\n", Parameter_SaveConfigurationEveryXXXResult);
  fprintf(file,"%d                                                      : Save State-Descriptor every xxx results\n", Parameter_SaveStateDescriptorEveryXXXResult);
  fprintf(file,"%d                                                      : Determine EW every xxx confs. (PREC-phase)\n", Parameter_DetermineEWeveryXXXconfsPREC);
  fprintf(file,"%d                                                      : Determine EW every xxx confs. (MEASURE-phase)\n", Parameter_DetermineEWeveryXXXconfsMEASURE);
  fprintf(file,"%d                                                      : Fourier Acceleration type\n", Parameter_FACC_type);
  fprintf(file,"%1.15f                                      : Fourier Acceleration parameter\n", Parameter_FACC_parameter);      
  fprintf(file,"%d                                                      : Spherical Higgs-Integration Mode\n", Parameter_SphericalHiggsMode);
  fprintf(file,"%1.15f                                      : Zeta-value for Higgs-Integration Mode\n", Parameter_SphericalHiggsZeta);    
  fprintf(file,"%d                                                      : Omega Mass Adaption Mode\n", Parameter_OmegaMassAdaptionMode);
  fprintf(file,"%d                                                      : Reversibility check frequency, PREC-phase\n", Parameter_ReversibilityCheckFreqPrec);
  fprintf(file,"%d                                                      : Reversibility check frequency, MEASURE-phase\n", Parameter_ReversibilityCheckFreqMeas);
  fprintf(file,"%d                                                      : Start randseed (0: Random Start)\n", Parameter_StartRandSeed);
  fprintf(file,"%d                                                      : Current Phi-Conf-ID\n", phiConfIDnr);  
  fprintf(file,"%d                                                      : Multi-Threaded Operation Mode\n", Parameter_MultiThreadedOpMode);
  fprintf(file,"%d                                                      : Number of threads used for each FFTW Transformation\n", Parameter_FFTWThreadCount);
  fprintf(file,"%d                                                      : Number of threads used for each xFFT Transformation\n", Parameter_xFFTThreadCount);      
  fprintf(file,"%1.15f                                      : Maximal RunTime in hours\n", Parameter_MaxRunTime);
  fprintf(file,"%d                                                      : Log - Level\n", LogLevel);
  fprintf(file,"%d                                             : Batch-Job-Nr \n", JobNr);  
  fprintf(file,"%d                                             : Current RandSeed \n", AdvancedSeed);
  
  
  fprintf(file,"\n   *** Phi - Field ***\n");
  int VL = fermiOps->getVectorLength();
  for (I=0; I<VL/8; I++) {
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f\n",pHMCProp->phiField[I][0],pHMCProp->phiField[I][1],pHMCProp->phiField[I][2],pHMCProp->phiField[I][3]);
  }  
  fclose(file);  
  delete [] fileName;
}


void writePhiConfigurationToDisk(char *confFileName, double weight) {
  if (LogLevel>2) printf("Writing Phi - configuration to disk...\n");

  std::fstream confFile;
  confFile.open(confFileName, std::ios::out);

  if (!confFile.is_open()) {
    printf("ERROR: Could not open conf-file: %s\n",confFileName);
    exit(0);
  }
  
  double* phiField = (double*) pHMCProp->phiField;
  int VL = fermiOps->getVectorLength();
  confFile.write((char*)phiField, 4*VL);
  confFile.write((char*)(&weight), 8);
  confFile.flush();
  confFile.close();
}


//Call this routine AFTER writeStateDescriptor !!!
void writePhiConfigurationToDisk(double weight) {
  if (Parameter_SaveConfigurationEveryXXXResult == 0) {
    return;
  }
  if (phiConfCounter == -1) {
    phiConfCounter = stateDescriptorCounter;
  } 
  phiConfCounter++;

  if (phiConfCounter < Parameter_SaveConfigurationEveryXXXResult) {
    return;
  }
  phiConfCounter = 0;
  
  phiConfIDnr++;
  char* fileName = new char[1200]; 
  
  snprintf(fileName,1200,"%s/data/results/pHMC/configurations/subFolder%s/PhiConf%s_%d.dat",DataBaseDirectory,outputFileNameExtension,outputFileNameExtension,phiConfIDnr);
  if (phiConfIDnr <= 1) {
    char* dirCreationCommand = new char[1200]; 
    snprintf(dirCreationCommand,1200,"mkdir %s/data/results/pHMC/configurations/subFolder%s",DataBaseDirectory,outputFileNameExtension);
    system(dirCreationCommand);
    delete[] dirCreationCommand;
  }
  
  writePhiConfigurationToDisk(fileName, weight);
  
  delete[] fileName;
}


void buildOutputFileNameExtension() {
  outputFileNameExtension = new char[300];
  
  snprintf(outputFileNameExtension,300,"L%dx%dx%dx%dNf%dKap%1.5fLam%1.5fY%1.5fRho%1.3fR%1.3fPolDeg%dPolAl%1.3f_%s",
   Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3,Parameter_Nf,Parameter_Kappa,Parameter_Lambda,Parameter_Y,Parameter_RHO,Parameter_R,Parameter_PolyDegree[0],Parameter_PolyAlpha,Parameter_filenameSuffix);	 
}  


void measureAndWrite(double weight) {
  double measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm;
  pHMCProp->measure(measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm);  
    
  char* fileName = new char[600]; 
  snprintf(fileName,600,"%s/data/results/pHMC/measure/measure%s.dat",DataBaseDirectory,outputFileNameExtension);	 
  FILE* file;

  file = fopen(fileName,"a");
  double logWeight = 0;
  if ((weight>0) && (!isNaN(weight))) logWeight = log(weight);
  fprintf(file,"%f %f %f nan nan nan nan nan nan nan\n", measurePhiNorm,measureStaggeredPhiNorm,logWeight);
  fclose(file);  

  delete[] fileName;
}


void writeConditionNumber(int nr, int phase) {
  double eigMin;
  double eigMax;
  double invCond;

  double* phiField = (double*) pHMCProp->phiField;
  fermiOps->exactFermionMatrixConditionNumber(phiField, eigMin, eigMax, invCond, true, 1, Parameter_FLAG_QuasiHermiteanMode);

  char* fileName = new char[600]; 
  if (phase == 0) {
    snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/%s%s_Prec.dat",DataBaseDirectory,"conditionNumber",outputFileNameExtension);	 
  } else if (phase == 1) {
    snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/%s%s_Therm.dat",DataBaseDirectory,"conditionNumber",outputFileNameExtension);	 
  } else {
    snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/%s%s_Meas.dat",DataBaseDirectory,"conditionNumber",outputFileNameExtension);	 
  }
  
  FILE* file;

  file = fopen(fileName,"a");
  fprintf(file,"%d %f %f %f\n", nr, eigMin, eigMax, invCond);
  fclose(file);  

  delete[] fileName;
}


void writeWeightFactor(int nr, int phase, double weight) {
  char* fileName = new char[600]; 
  if (phase == 0) {
    snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/%s%s_Prec.dat",DataBaseDirectory,"weightFactor",outputFileNameExtension);	 
  } else if (phase == 1) {
    snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/%s%s_Therm.dat",DataBaseDirectory,"weightFactor",outputFileNameExtension);	 
  } else {
    snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/%s%s_Meas.dat",DataBaseDirectory,"weightFactor",outputFileNameExtension);	 
  }
  
  FILE* file;

  file = fopen(fileName,"a");
  fprintf(file,"%d %1.15f\n", nr, weight);
  fclose(file);

  delete[] fileName;
}


bool markovStep() {
  int I;

  int LevelCount = 2 + Parameter_SubPolyCount;
  int* iterations = new int[LevelCount];
  for (I=0; I<1+Parameter_SubPolyCount; I++) {
    iterations[1+Parameter_SubPolyCount-I] = Parameter_Iterations[I];
  }
  iterations[0] = Parameter_Iterations[1+SubPolynomMAX];
  int* integrators = new int[LevelCount];
   for (I=0; I<1+Parameter_SubPolyCount; I++) {
    integrators[1+Parameter_SubPolyCount-I] = Parameter_IntegratorType[I];
  }
  integrators[0] = Parameter_IntegratorType[1+SubPolynomMAX];
  int* subPolNr = new int[LevelCount];
  for (I=2; I<LevelCount; I++) subPolNr[I] = 2*Parameter_SubPolyCount-2*I+3;
  subPolNr[0] = 0;
  subPolNr[1] = 2*Parameter_SubPolyCount;  
  double* propTOL = NULL;
  double finalTOL = 0;

  bool b = pHMCProp->multiTimeScaleMarkovStep(LevelCount, iterations, Parameter_Epsilon, integrators, subPolNr, propTOL, finalTOL);

  return b;
}


void purePropagation() {
  int I;

  int LevelCount = 2 + Parameter_SubPolyCount;
  int* iterations = new int[LevelCount];
  for (I=0; I<1+Parameter_SubPolyCount; I++) {
    iterations[1+Parameter_SubPolyCount-I] = Parameter_Iterations[I];
  }
  iterations[0] = Parameter_Iterations[1+SubPolynomMAX];
  int* integrators = new int[LevelCount];
   for (I=0; I<1+Parameter_SubPolyCount; I++) {
    integrators[1+Parameter_SubPolyCount-I] = Parameter_IntegratorType[I];
  }
  integrators[0] = Parameter_IntegratorType[1+SubPolynomMAX];
  int* subPolNr = new int[LevelCount];
  for (I=2; I<LevelCount; I++) subPolNr[I] = 2*Parameter_SubPolyCount-2*I+3;
  subPolNr[0] = 0;
  subPolNr[1] = 2*Parameter_SubPolyCount;  
  double* propTOL = NULL;
  double finalTOL = 0;

  pHMCProp->multiTimeScalePropagation(LevelCount, iterations, Parameter_Epsilon, integrators, subPolNr, propTOL, finalTOL);
}


void firstNMarkovSteps(int nSteps, int facIter, double facEps, double facSafety) {
  int iter0Merker = Parameter_Iterations[0];
  double epsMerker = Parameter_Epsilon;

  Parameter_Iterations[0] *= facIter;
  Parameter_Epsilon /= facEps;

  if (LogLevel>1) printf("Performing first %d Markov Steps with facIter = %d, facEps = %1.2f, facSafety = %1.2f.\n",nSteps, facIter, facEps,facSafety);
  for (int I=0; I<nSteps; I++) {
    if (Parameter_FLAG_QuasiHermiteanMode) pHMCProp->improveRPreconditioningParameters(-1, 1, Parameter_PolyLambda, facSafety*Parameter_UpperEWboundSafetyFactor);
    if (!Parameter_FLAG_QuasiHermiteanMode) pHMCProp->improvePreconditioningParametersFAST();
    pHMCProp->sampleALLMomenta();
    if (Parameter_FLAG_DirectOmegaSampling > 0) pHMCProp->sampleOmegaFields();

    purePropagation();
  }
  Parameter_Iterations[0] = iter0Merker;
  Parameter_Epsilon = epsMerker;
  if (LogLevel>1) printf("First %d Markov Steps Ready.\n",nSteps);
}


void tuneIntegrators(double allowedDis, bool afterMassDet) {
  if (LogLevel>2) printf("Tuning Integrators...\n");

  pHMCProp->saveALLfields();
  pHMCProp->synchronizedChangeOfTuneMode(true);
  
  markovStep();
  double deltaS0 = pHMCProp->SafterProp - pHMCProp->SbeforeProp;
  if ((pHMCProp->SafterProp==0) || (pHMCProp->SbeforeProp==0)) deltaS0 = NaN;
  pHMCProp->restoreALLfields(true);
  
  //Adaption of sub-Polynomial-Integrators
  for (int I=0; I<1+Parameter_SubPolyCount; I++) {
    int Ip1 = I+1;
    if (Ip1 == 1+Parameter_SubPolyCount) Ip1 = 1+SubPolynomMAX;    
    
    while (true) {
      if (Parameter_Iterations[I]==1) break;
      if ((I==0) && (!afterMassDet) && (Parameter_Iterations[I]<=MinimalIterationCountForP0BeforeMassDet)) break;
    
      int deltaI = Parameter_Iterations[I]/10;
      if (deltaI<=0) deltaI=1;
      Parameter_Iterations[I] -= deltaI;
      int itMerker = Parameter_Iterations[Ip1];
      double epsMerker = Parameter_Epsilon;
      if (I==0) Parameter_Epsilon = (Parameter_Epsilon*(Parameter_Iterations[I]+deltaI)) / Parameter_Iterations[I];
      Parameter_Iterations[Ip1] = 1+((Parameter_Iterations[Ip1]*(Parameter_Iterations[I]+deltaI)) / Parameter_Iterations[I]);
    
      markovStep();
      double deltaS1 = pHMCProp->SafterProp - pHMCProp->SbeforeProp;
      if ((pHMCProp->SafterProp==0) || (pHMCProp->SbeforeProp==0)) deltaS1 = NaN;
      pHMCProp->restoreALLfields(true);
      
      if (LogLevel>2) {
        printf("Setting:  theta = %f, eps=%f, iter= ", Parameter_Theta, Parameter_Epsilon);  
        for (int I2=0; I2<1+Parameter_SubPolyCount; I2++) printf("%d, ",Parameter_Iterations[I2]);
        printf("%d (Higgs) ",Parameter_Iterations[1+SubPolynomMAX]);
        printf("tried with DeltaDeltaS = %f\n",abs(deltaS0-deltaS1));
      }
      
      if ((isNaN(deltaS0-deltaS1)) || (abs(deltaS0-deltaS1) > allowedDis)) {
        Parameter_Iterations[I] += deltaI;
        Parameter_Iterations[Ip1] = itMerker;
	Parameter_Epsilon = epsMerker;
	break;
      }
    }
  }
  //Adaption of Higgs-Integrator
  while (true) {
    if (Parameter_Iterations[1+SubPolynomMAX]==1) break;
      
    int deltaI = Parameter_Iterations[1+SubPolynomMAX]/5;
    if (deltaI<=0) deltaI=1;
    Parameter_Iterations[1+SubPolynomMAX] -= deltaI;
    
    markovStep();
    double deltaS1 = pHMCProp->SafterProp - pHMCProp->SbeforeProp;
    if ((pHMCProp->SafterProp==0) || (pHMCProp->SbeforeProp==0)) deltaS1 = NaN;
    pHMCProp->restoreALLfields(true);

    if (LogLevel>2) {
      printf("Setting:  theta = %f, eps=%f, iter= ", Parameter_Theta, Parameter_Epsilon);  
      for (int I2=0; I2<1+Parameter_SubPolyCount; I2++) printf("%d, ",Parameter_Iterations[I2]);
      printf("%d (Higgs) ",Parameter_Iterations[1+SubPolynomMAX]);
      printf("tried with DeltaDeltaS = %f\n",abs(deltaS0-deltaS1));
    }
      
    if ((isNaN(deltaS0-deltaS1)) || (abs(deltaS0-deltaS1) > allowedDis)) {
      Parameter_Iterations[1+SubPolynomMAX] += deltaI;
      break;
    }
  }
  //Adaption of Theta-Setting
  if ((afterMassDet) && (Parameter_FLAG_ThetaScan>0)) {
    while (true) {
      if (Parameter_Theta*Parameter_ThetaFactor>Parameter_ThetaMax) break;

      Parameter_Theta *= Parameter_ThetaFactor;
      pHMCProp->setTheta(Parameter_Theta);  
      
      markovStep();
      double deltaS1 = pHMCProp->SafterProp - pHMCProp->SbeforeProp;
      if ((pHMCProp->SafterProp==0) || (pHMCProp->SbeforeProp==0)) deltaS1 = NaN;
      pHMCProp->restoreALLfields(true);

      if (LogLevel>2) {
        printf("Setting:  theta = %f, eps=%f, iter= ", Parameter_Theta, Parameter_Epsilon);  
        for (int I2=0; I2<1+Parameter_SubPolyCount; I2++) printf("%d, ",Parameter_Iterations[I2]);
        printf("%d (Higgs) ",Parameter_Iterations[1+SubPolynomMAX]);
        printf("tried with DeltaDeltaS = %f\n",abs(deltaS0-deltaS1));
      }
      
      if ((isNaN(deltaS0-deltaS1)) || (abs(deltaS0-deltaS1) > allowedDis)) {
        Parameter_Theta /= Parameter_ThetaFactor;
        break;
      }
    }
    pHMCProp->setTheta(Parameter_Theta);  
  }
  pHMCProp->synchronizedChangeOfTuneMode(false);
  //Final Setting
  if (LogLevel>2) {
    printf("Tuned Integrator setting is:  theta = %f, eps=%f, iter = ", Parameter_Theta,Parameter_Epsilon);  
    for (int I=0; I<1+Parameter_SubPolyCount; I++) printf("%d, ",Parameter_Iterations[I]);
    printf("%d (Higgs)\n\n",Parameter_Iterations[1+SubPolynomMAX]);
  }
}


void automaticEpsilonAdaption(double acceptRate, bool afterMassDet) {
#ifdef IntegrationSchemeAdaption
  if (acceptRate<AutomaticAdaption_LowPro) {
    int iterNeu = (int)(Parameter_Iterations[0]/0.90);
    if (iterNeu <= Parameter_Iterations[0]) iterNeu = Parameter_Iterations[0]+1;
    Parameter_Epsilon = (Parameter_Iterations[0]*Parameter_Epsilon) / iterNeu;    
    Parameter_Iterations[0] = iterNeu;
    for (int I=1; I<1+Parameter_SubPolyCount; I++) {
      Parameter_Iterations[I]++;
    }
    Parameter_Iterations[1+SubPolynomMAX]++;
    if ((afterMassDet) && (Parameter_FLAG_ThetaScan>0)) Parameter_Theta /= Parameter_ThetaFactor; 
    if ((Parameter_Theta < Parameter_ThetaMin) && (Parameter_FLAG_ThetaScan>0)) Parameter_Theta = Parameter_ThetaMin;
    if (LogLevel>2) printf("==> Decreasing epsilon to %f and increasing iterations[0] to %d.\n",Parameter_Epsilon, Parameter_Iterations[0]);
    tuneIntegrators(0.15, afterMassDet);
  } else if (acceptRate<AutomaticAdaption_LowPro2) {
    int iterNeu = (int)(2.00 * Parameter_Iterations[0]);
    if (iterNeu <= Parameter_Iterations[0]) iterNeu = Parameter_Iterations[0]+1;    
    Parameter_Epsilon = (Parameter_Iterations[0]*Parameter_Epsilon) / iterNeu;    
    Parameter_Iterations[0] = iterNeu;
    for (int I=1; I<1+Parameter_SubPolyCount; I++) {
      Parameter_Iterations[I]++;
    }
    Parameter_Iterations[1+SubPolynomMAX]++;
    if ((afterMassDet) && (Parameter_FLAG_ThetaScan>0)) Parameter_Theta = Parameter_ThetaMin;
    if (LogLevel>2) printf("==> Decreasing epsilon to %f and increasing iterations[0] to %d.\n",Parameter_Epsilon, Parameter_Iterations[0]);
    tuneIntegrators(0.20, afterMassDet);
  } else if (acceptRate>AutomaticAdaption_HighPro) {
    int polNr = (int)((2+Parameter_SubPolyCount)*pHMCProp->getSyncSingleRandom());  
    polNr--;
    if ((Parameter_FLAG_ThetaScan == 0) && (polNr<0)) polNr = 0;    
    if (polNr>Parameter_SubPolyCount) polNr = Parameter_SubPolyCount;
    if (polNr>=0) {
      if (Parameter_Iterations[polNr] == 1) {
        polNr = -1;
        for (int I=0; I<1+Parameter_SubPolyCount; I++) {
          if (Parameter_Iterations[I] > 1) {
	    polNr = I;
	    break;
	  }
	}
      }
    }
    if (polNr == -1) {
      if ((afterMassDet) && (Parameter_FLAG_ThetaScan>0)) {
        Parameter_Theta *= Parameter_ThetaFactor;
        tuneIntegrators(0.001, afterMassDet);  
      } else {
        tuneIntegrators(0.2, afterMassDet);  
      }
    } else {
      if ((polNr!=0) || (afterMassDet) || (Parameter_Iterations[polNr]>MinimalIterationCountForP0BeforeMassDet)) {
        int Ip1 = polNr+1;
        if (Ip1 == 1+Parameter_SubPolyCount) Ip1 = 1+SubPolynomMAX;    
        Parameter_Iterations[polNr]--;      
        if (polNr==0) Parameter_Epsilon = (Parameter_Epsilon*(Parameter_Iterations[polNr]+1)) / Parameter_Iterations[polNr];
        Parameter_Iterations[Ip1] = 2+((Parameter_Iterations[Ip1]*(Parameter_Iterations[polNr]+1)) / Parameter_Iterations[polNr]);

        if (LogLevel>2) {
          printf("==> Improving Integrators to: eps=%f, iter = ",Parameter_Epsilon);
          for (int I=0; I<1+Parameter_SubPolyCount; I++) printf("%d, ",Parameter_Iterations[I]);
          printf("%d (Higgs)\n\n",Parameter_Iterations[1+SubPolynomMAX]);
        }
        tuneIntegrators(0.001, afterMassDet);
      }
    }
  } else {
    tuneIntegrators(0.10, afterMassDet);  
  }
  if ((afterMassDet) && (Parameter_FLAG_ThetaScan>0)) {
    int level = roundToInt((log(Parameter_Theta/Parameter_ThetaMin) / log(Parameter_ThetaFactor)));
    Parameter_Theta = Parameter_ThetaMin*exp(level*log(Parameter_ThetaFactor));
    if (Parameter_Theta<Parameter_ThetaMin) Parameter_Theta = Parameter_ThetaMin;
    if (Parameter_Theta>Parameter_ThetaMax) Parameter_Theta = Parameter_ThetaMax;
    pHMCProp->setTheta(Parameter_Theta);   
  }
#endif
}


// ONLY call after accepted integration!!!
bool calcReversibilityNorm(double& difNorm, double& difS) {
  if (LogLevel>2) printf("Reversibility Test...\n");
  vector4D* savedPhi = (vector4D*)fermiOps->createFermionVector(2);
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  pHMCProp->getCopyOfSavedPhi((double*)savedPhi);
  pHMCProp->saveALLfields();
  difNorm = 0;
  bool acc = false;  
  double S1 = pHMCProp->SbeforeProp;
  difS = 0;
  
  pHMCProp->negateAllMomenta();
  if (markovStep()) {
    double S2 = pHMCProp->SafterProp;
    difS = S2 - S1;
    int I;
    for (I=0; I<L0*L1*L2*L3; I++) {
      int I2;
      for (I2=0; I2<4; I2++) {
        double d = pHMCProp->phiField[I][I2]-savedPhi[I][I2];
	difNorm += d*d;
      }
    }
    difNorm = difNorm / (4*L0*L1*L2*L3);
    difNorm = sqrt(difNorm);
    acc = true;
    if (LogLevel>2) printf("Reversibility: Dif. Norm: %1.3e, Dif. Action: %1.3e\n",difNorm,difS);    
  }
  pHMCProp->restoreALLfields(true);

  Complex* c1 = (Complex*) savedPhi;
  fermiOps->destroyFermionVector(c1);
  return acc;
}


void writeReversibilityNorm(int runPhase) {
  char* fileName = new char[600]; 
  if (runPhase==0) {
    snprintf(fileName,600,"%s/data/results/pHMC/revers/ReversibilityCheckPREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  }
  if (runPhase==1) {
    snprintf(fileName,600,"%s/data/results/pHMC/revers/ReversibilityCheckTHERM%s.dat",DataBaseDirectory,outputFileNameExtension);
  }
  if (runPhase==2) {
    snprintf(fileName,600,"%s/data/results/pHMC/revers/ReversibilityCheckMEAS%s.dat",DataBaseDirectory,outputFileNameExtension);
  } 
  
  double difReversNorm = 0;
  double divReversAction = 0;
  if (calcReversibilityNorm(difReversNorm,divReversAction)) {
    FILE* file = fopen(fileName,"a");
    fprintf(file,"%1.20f %1.20f\n",difReversNorm,divReversAction);
    fclose(file);
  }
  
  delete[] fileName;
}


void iniMPI(int argc,char **argv) {
#ifdef UseMPI
  int namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];  //Initialisierung von MPICH...
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nodeCount);
  MPI_Comm_rank(MPI_COMM_WORLD,&ownNodeID);
  MPI_Get_processor_name(processor_name,&namelen);
  printf("Simulation runs on %d nodes.\n", nodeCount);
  printf("Process %d runs on node %s.\n", ownNodeID, processor_name);
#else
  printf("Simulation runs on %d nodes. OwnID = %d\n", nodeCount, ownNodeID);
#endif 
}


void desini() {
  delete[] outputFileNameExtension;
  pHMCProp->killSlaves();
  delete pHMCProp;
  delete fermiOps;
  #ifdef UseMPI
  MPI_Abort(MPI_COMM_WORLD, 0);
  #endif
}


bool readCurrentStateDescriptor() {
  //Default Settings if no Descriptor-File is found
  AutomaticPreconRunCount = 0;
  ThermalizingRunCount = 0;
  TotallyMeasuredConfigurationsCount = 0;
  //Default Setting END
  
  char* fileName = new char[600]; 
  snprintf(fileName,600,"%s/data/results/pHMC/states/StateDescriptor%s.dat",DataBaseDirectory,outputFileNameExtension);
  FILE* file;
  file = fopen(fileName,"r");
  if (LogLevel>1) printf("Trying to read DescriptorStateFile: %s\n",fileName);  
  delete[] fileName;
  char *Comment = new char[500];


  if (file==NULL) {
    if (LogLevel>0) printf("StateDescriptor File does not exist.\n");
    return false;
  }


  int running = -1;
  if ((fscanf(file,"%d",&running)!=1) || (running<0) || (running>1)) {
    if (LogLevel>0) printf("Could not read from StateDescriptor File!!! ==> Exiting \n");
    fclose(file);
    desini();
    exit(0);
  }
  if (running==1) {
    if (LogLevel>0) printf("CurrentStateDescriptor indicates running!!! ==> Exiting\n");
    fclose(file);
    desini();
    exit(0);
  }
  fgets(Comment, 500, file);

  
  int readError = 0;
  double READmeasurePhiNorm,READmeasureStaggeredPhiNorm,READprecM,READprecS,READParameter_RHO,READParameter_R, 
         READParameter_Lambda,READParameter_Kappa,READParameter_Y,READParameter_Epsilon,
         READmeasurePhiVectorLength, READmeasurePhiVectorLengthSigma;
  
  double READParameter_PolyEpsilon[1+SubPolynomMAX],READParameter_PolyAlpha,
         READParameter_PolyLambda, READParameter_Theta,
	 READParameter_ThetaMin,READParameter_ThetaMax,READParameter_ThetaFactor,
         READParameter_MaxRunTime,READParameter_FACC_parameter,
	 READQprecMu,READQprecBeta,READRprecM,READRprecF,READParameter_QPreconditioner_Mu,
	 READParameter_QPreconditioner_Beta,READParameter_UpperEWboundSafetyFactor,
	 READParameter_MassSplit, READParameter_ExplicitMass, READParameter_ExplicitCurrent,
	 READParameter_SphericalHiggsZeta, READParameter_c6, READParameter_c8, READParameter_c10,
         READParameter_lambda6, READParameter_lambda8, READParameter_lambda10;	 

    
  int READusePrec,READAutomaticPreconRunCount,READThermalizingRunCount,READTotallyMeasuredConfigurationsCount,
      READParameter_L0,READParameter_L1,READParameter_L2,READParameter_L3,
      READParameter_Nf,READParameter_FLAG_ThetaScan, READParameter_IntegratorType[2+SubPolynomMAX],
      READParameter_Iterations[2+SubPolynomMAX],READParameter_Measurements,READParameter_AutomaticPreconditioningMetros,
      READParameter_ThermalizingMetros, READParameter_TotalData,READAdvancedSeed,
      READParameter_FACC_type,READParameter_StartRandSeed,READParameter_ReversibilityCheckFreqPrec,
      READParameter_ReversibilityCheckFreqMeas, READParameter_OmegaMassAdaptionMode;
  
  int READParameter_PolyDegree[1+SubPolynomMAX],READParameter_PolDigit,  
      READParameter_MaxPolDegPerNode, READParameter_FLAG_xFFT, 
      READParameter_SaveConfigurationEveryXXXResult,
      READphiConfIDnr,READJobNr,READLogLevel, READParameter_SubPolyCount,
      READuseQPrec,READuseRPrec,READParameter_FLAG_ExactReweighing,READParameter_FLAG_QuasiHermiteanMode,
      READParameter_FLAG_FactorizationOfMatrixB, READParameter_FLAG_AntiPeriodicBoundaryConditionsInTime,
      READParameter_FLAG_Use_P_Preconditioner,READParameter_FLAG_Use_Q_Preconditioner,
      READParameter_FLAG_Use_R_Preconditioner,READParameter_AdditionalAuxVectors,
      READParameter_FLAG_DirectOmegaSampling,READParameter_FLAG_ModelSelection,
      READParameter_SaveStateDescriptorEveryXXXResult,READParameter_DetermineEWeveryXXXconfsPREC,
      READParameter_DetermineEWeveryXXXconfsMEASURE, READParameter_FLAG_RandomGauge,
      READParameter_SphericalHiggsMode,READParameter_MultiThreadedOpMode,READParameter_FFTWThreadCount,READParameter_xFFTThreadCount;
  
  
  int I;
  char* READParameter_filenameSuffix = new char[300];
  char* READoutputFileNameExtension = new char[300];  
  
  if (fscanf(file,"%lf", &READmeasurePhiNorm)!=1) readError = 1;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READmeasureStaggeredPhiNorm)!=1) readError = 2;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READmeasurePhiVectorLength)!=1) readError = 3;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READmeasurePhiVectorLengthSigma)!=1) readError = 4;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READAutomaticPreconRunCount)!=1) readError = 5;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READThermalizingRunCount)!=1) readError = 6;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READTotallyMeasuredConfigurationsCount)!=1) readError = 7;
  fgets(Comment, 500, file);
  
  if (fscanf(file,"%d", &READusePrec)!=1) readError = 8;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READprecM)!=1) readError = 9;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READprecS)!=1) readError = 10;
  fgets(Comment, 500, file);
    
  if (fscanf(file,"%d", &READuseQPrec)!=1) readError = 11;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READQprecMu)!=1) readError = 12;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READQprecBeta)!=1) readError = 13;
  fgets(Comment, 500, file);
    
  if (fscanf(file,"%d", &READuseRPrec)!=1) readError = 14;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READRprecM)!=1) readError = 15;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READRprecF)!=1) readError = 16;
  fgets(Comment, 500, file);    
  
  if (fscanf(file,"%s", READParameter_filenameSuffix)!=1) readError = 17;
  fgets(Comment, 500, file);
  if (fscanf(file,"%s", READoutputFileNameExtension)!=1) readError = 18;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_L0)!=1) readError = 19;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_L1)!=1) readError = 20;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_L2)!=1) readError = 21;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_L3)!=1) readError = 22;
  fgets(Comment, 500, file);  
  if (fscanf(file,"%d", &READParameter_Nf)!=1) readError = 23;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_RHO)!=1) readError = 24;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_R)!=1) readError = 25;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_Lambda)!=1) readError = 26;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_Kappa)!=1) readError = 27;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_Y)!=1) readError = 28;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf", &READParameter_MassSplit)!=1) readError = 281;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_ExplicitMass)!=1) readError = 282;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_ExplicitCurrent)!=1) readError = 283;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_c6)!=1) readError = 284;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_c8)!=1) readError = 285;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_c10)!=1) readError = 286;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_lambda6)!=1) readError = 287;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_lambda8)!=1) readError = 288;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_lambda10)!=1) readError = 289;
  fgets(Comment, 500, file);
  
  for (I=0; I<1+SubPolynomMAX; I++) {
    if (fscanf(file,"%lf", &(READParameter_PolyEpsilon[I]))!=1) readError = 29;
    fgets(Comment, 500, file);
  }
  if (fscanf(file,"%lf", &READParameter_PolyLambda)!=1) readError = 30;
  fgets(Comment, 500, file);
  for (I=0; I<1+SubPolynomMAX; I++) {
    if (fscanf(file,"%d", &(READParameter_PolyDegree[I]))!=1) readError = 31;
    fgets(Comment, 500, file);
  }
  if (fscanf(file,"%lf", &READParameter_PolyAlpha)!=1) readError = 32;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_PolDigit)!=1) readError = 33;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_MaxPolDegPerNode)!=1) readError = 34;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_RandomGauge)!=1) readError = 341;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_ExactReweighing)!=1) readError = 35;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_QuasiHermiteanMode)!=1) readError = 36;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_FactorizationOfMatrixB)!=1) readError = 361;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_AntiPeriodicBoundaryConditionsInTime)!=1) readError = 362;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_Use_P_Preconditioner)!=1) readError = 37;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_Use_Q_Preconditioner)!=1) readError = 38;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_Use_R_Preconditioner)!=1) readError = 39;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_QPreconditioner_Mu)!=1) readError = 40;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_QPreconditioner_Beta)!=1) readError = 41;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_UpperEWboundSafetyFactor)!=1) readError = 42;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_AdditionalAuxVectors)!=1) readError = 43;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_DirectOmegaSampling)!=1) readError = 44;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_Theta)!=1) readError = 45;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_ThetaScan)!=1) readError = 46;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_ThetaMin)!=1) readError = 47;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_ThetaMax)!=1) readError = 48;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_ThetaFactor)!=1) readError = 49;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_ModelSelection)!=1) readError = 50;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_xFFT)!=1) readError = 51;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf",&READParameter_Epsilon)!=1) readError = 52;
  fgets(Comment, 500, file);
  for (I=0; I<2+SubPolynomMAX; I++) {
    if (fscanf(file,"%d", &(READParameter_IntegratorType[I]))!=1) readError = 53;
    fgets(Comment, 500, file);
  }
  if (fscanf(file,"%d", &READParameter_SubPolyCount)!=1) readError = 54;
  fgets(Comment, 500, file);
  for (I=0; I<2+SubPolynomMAX; I++) {
    if (fscanf(file,"%d", &(READParameter_Iterations[I]))!=1) readError = 55;
    fgets(Comment, 500, file);
  }
  if (fscanf(file,"%d", &READParameter_Measurements)!=1) readError = 56;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_AutomaticPreconditioningMetros)!=1) readError = 57;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_ThermalizingMetros)!=1) readError = 58;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_TotalData)!=1) readError = 59;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_SaveConfigurationEveryXXXResult)!=1) readError = 60;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_SaveStateDescriptorEveryXXXResult)!=1) readError = 61;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_DetermineEWeveryXXXconfsPREC)!=1) readError = 62;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_DetermineEWeveryXXXconfsMEASURE)!=1) readError = 63;
  fgets(Comment, 500, file);  
  if (fscanf(file,"%d", &READParameter_FACC_type)!=1) readError = 64;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf",&READParameter_FACC_parameter)!=1) readError = 65;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_SphericalHiggsMode)!=1) readError = 651;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf",&READParameter_SphericalHiggsZeta)!=1) readError = 652;
  fgets(Comment, 500, file);    
  if (fscanf(file,"%d", &READParameter_OmegaMassAdaptionMode)!=1) readError = 66;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_ReversibilityCheckFreqPrec)!=1) readError = 67;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_ReversibilityCheckFreqMeas)!=1) readError = 68;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_StartRandSeed)!=1) readError = 69;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READphiConfIDnr)!=1) readError = 70;
  fgets(Comment, 500, file);    
  if (fscanf(file,"%d", &READParameter_MultiThreadedOpMode)!=1) readError = 701;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FFTWThreadCount)!=1) readError = 702;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_xFFTThreadCount)!=1) readError = 703;
  fgets(Comment, 500, file);    
  if (fscanf(file,"%lf", &READParameter_MaxRunTime)!=1) readError = 71;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READLogLevel)!=1) readError = 72;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READJobNr)!=1) readError = 73;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READAdvancedSeed)!=1) readError = 74;
  fgets(Comment, 500, file);

  if (readError>0) {
    fclose(file);
    if (LogLevel>0) printf("Error (line %d) while reading File-Descriptor. ==> Exiting\n",readError);
    desini();
    exit(0);  
  }
  
  //Consistency Check
  int consistencyError=0;
  double ceps = 1E-13;
  
  if (strlen(Parameter_filenameSuffix) != strlen(READParameter_filenameSuffix)) {
    consistencyError=1;
  } else {
    for (I=0; I<((int)strlen(Parameter_filenameSuffix)); I++) {
      if (Parameter_filenameSuffix[I] != READParameter_filenameSuffix[I]) {
        consistencyError=2;
	break;
      }
    }
  }
  
  if (Parameter_L0!=READParameter_L0) consistencyError=3;
  if (Parameter_L1!=READParameter_L1) consistencyError=4;
  if (Parameter_L2!=READParameter_L2) consistencyError=5;
  if (Parameter_L3!=READParameter_L3) consistencyError=6;
  if (Parameter_Nf!=READParameter_Nf) consistencyError=7;
  if (abs(Parameter_RHO-READParameter_RHO)>ceps) consistencyError=8;
  if (abs(Parameter_R-READParameter_R)>ceps) consistencyError=9;
  if (abs(Parameter_Lambda-READParameter_Lambda)>ceps) consistencyError=10;
  if (abs(Parameter_Kappa-READParameter_Kappa)>ceps) consistencyError=11;
  if (abs(Parameter_Y-READParameter_Y)>ceps) consistencyError=12;
  if (abs(Parameter_MassSplit-READParameter_MassSplit)>ceps) consistencyError=121;
  if (abs(Parameter_ExplicitMass-READParameter_ExplicitMass)>ceps) consistencyError=122;
  if (abs(Parameter_ExplicitCurrent-READParameter_ExplicitCurrent)>ceps) consistencyError=123;
  if (abs(Parameter_c6-READParameter_c6)>ceps) consistencyError=124;
  if (abs(Parameter_c8-READParameter_c8)>ceps) consistencyError=125;
  if (abs(Parameter_c10-READParameter_c10)>ceps) consistencyError=126;
  if (abs(Parameter_lambda6-READParameter_lambda6)>ceps) consistencyError=127;
  if (abs(Parameter_lambda8-READParameter_lambda8)>ceps) consistencyError=128;
  if (abs(Parameter_lambda10-READParameter_lambda10)>ceps) consistencyError=129;
  
  for (I=0; I<1+SubPolynomMAX; I++) {  
    if (abs(Parameter_PolyEpsilon[I]-READParameter_PolyEpsilon[I])>ceps) consistencyError=13;
  }
  if (abs(Parameter_PolyLambda-READParameter_PolyLambda)>ceps) consistencyError=14;
  for (I=0; I<1+SubPolynomMAX; I++) {
    if (Parameter_PolyDegree[I]!=READParameter_PolyDegree[I]) consistencyError=15;
  }
  if (abs(Parameter_PolyAlpha-READParameter_PolyAlpha)>ceps) consistencyError=16;
  if (Parameter_PolDigit!=READParameter_PolDigit) consistencyError=17;
  if (Parameter_MaxPolDegPerNode!=READParameter_MaxPolDegPerNode) consistencyError=18;
  if (Parameter_FLAG_RandomGauge!=READParameter_FLAG_RandomGauge) consistencyError=181;
  if (Parameter_FLAG_ExactReweighing!=READParameter_FLAG_ExactReweighing) consistencyError=19;
  if (Parameter_FLAG_QuasiHermiteanMode!=READParameter_FLAG_QuasiHermiteanMode) consistencyError=20;
  if (Parameter_FLAG_FactorizationOfMatrixB!=READParameter_FLAG_FactorizationOfMatrixB) consistencyError=201;
  if (Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime!=READParameter_FLAG_AntiPeriodicBoundaryConditionsInTime) consistencyError=202;
  if (Parameter_FLAG_Use_P_Preconditioner!=READParameter_FLAG_Use_P_Preconditioner) consistencyError=21;
  if (Parameter_FLAG_Use_Q_Preconditioner!=READParameter_FLAG_Use_Q_Preconditioner) consistencyError=22;
  if (Parameter_FLAG_Use_R_Preconditioner!=READParameter_FLAG_Use_R_Preconditioner) consistencyError=23;
  if (abs(Parameter_QPreconditioner_Mu-READParameter_QPreconditioner_Mu)>ceps) consistencyError=24;
  if (abs(Parameter_QPreconditioner_Beta-READParameter_QPreconditioner_Beta)>ceps) consistencyError=25;
  if (abs(Parameter_UpperEWboundSafetyFactor-READParameter_UpperEWboundSafetyFactor)>ceps) consistencyError=26;
//  if (Parameter_AdditionalAuxVectors!=READParameter_AdditionalAuxVectors) consistencyError=27;
  if (Parameter_FLAG_DirectOmegaSampling!=READParameter_FLAG_DirectOmegaSampling) consistencyError=28;
  if (Parameter_FLAG_ThetaScan!=READParameter_FLAG_ThetaScan) consistencyError=29;    
  if (abs(Parameter_ThetaMin-READParameter_ThetaMin)>ceps) consistencyError=30;
  if (abs(Parameter_ThetaMax-READParameter_ThetaMax)>ceps) consistencyError=31;
  if (abs(Parameter_ThetaFactor-READParameter_ThetaFactor)>ceps) consistencyError=32;
  if (Parameter_FLAG_ModelSelection!=READParameter_FLAG_ModelSelection) consistencyError=33;
//  if (Parameter_FLAG_xFFT!=READParameter_FLAG_xFFT) consistencyError=34;
  if (Parameter_SubPolyCount!=READParameter_SubPolyCount) consistencyError=35;
  if (Parameter_Measurements!=READParameter_Measurements) consistencyError=36;
  if (Parameter_AutomaticPreconditioningMetros!=READParameter_AutomaticPreconditioningMetros) consistencyError=37;
  if (Parameter_ThermalizingMetros!=READParameter_ThermalizingMetros) consistencyError=38;
  if (Parameter_SaveConfigurationEveryXXXResult!=READParameter_SaveConfigurationEveryXXXResult) consistencyError=39;
  if (Parameter_SaveStateDescriptorEveryXXXResult!=READParameter_SaveStateDescriptorEveryXXXResult) consistencyError=40;
  if (Parameter_DetermineEWeveryXXXconfsPREC!=READParameter_DetermineEWeveryXXXconfsPREC) consistencyError=41;
  if (Parameter_DetermineEWeveryXXXconfsMEASURE!=READParameter_DetermineEWeveryXXXconfsMEASURE) consistencyError=42;
  if (Parameter_FACC_type!=READParameter_FACC_type) consistencyError=43;
  if (abs(Parameter_FACC_parameter-READParameter_FACC_parameter)>ceps) consistencyError=44;
  if (Parameter_OmegaMassAdaptionMode!=READParameter_OmegaMassAdaptionMode) consistencyError=45;
  if (Parameter_ReversibilityCheckFreqPrec!=READParameter_ReversibilityCheckFreqPrec) consistencyError=46;
  if (Parameter_ReversibilityCheckFreqMeas!=READParameter_ReversibilityCheckFreqMeas) consistencyError=47;
  if (Parameter_StartRandSeed!=READParameter_StartRandSeed) consistencyError=48;  
  if (Parameter_SphericalHiggsMode!=READParameter_SphericalHiggsMode) consistencyError=49;  
  if (abs(Parameter_SphericalHiggsZeta-READParameter_SphericalHiggsZeta)>ceps) consistencyError=50;  

 
  if (consistencyError>0) {
    fclose(file);
    if (LogLevel>0) printf("Consistency-Error (check nr. %d) with data read from File-Descriptor. ==> Exiting\n",consistencyError);
    desini();
    exit(0);  
  }


  //Take over Data
  if (LogLevel>0) printf("Taking over values from StateDescriptor File...\n");
  fermiOps->setPreconditioner(READusePrec, READprecM, READprecS);
  fermiOps->setQPreconditioner(READuseQPrec, READQprecMu, READQprecBeta);
  fermiOps->setRPreconditioner(READuseRPrec, READRprecM, READRprecF);
    
  AutomaticPreconRunCount = READAutomaticPreconRunCount;
  if (LogLevel>0) printf("  ***> AutomaticPreconRunCount=%d\n",AutomaticPreconRunCount);
  ThermalizingRunCount = READThermalizingRunCount;
  if (LogLevel>0) printf("  ***> ThermalizingRunCount=%d\n",ThermalizingRunCount);
  TotallyMeasuredConfigurationsCount = READTotallyMeasuredConfigurationsCount;
  if (LogLevel>0) printf("  ***> TotallyMeasuredConfigurationsCount=%d\n",TotallyMeasuredConfigurationsCount);
  
  snprintf(Parameter_filenameSuffix,600,"%s",READParameter_filenameSuffix);
  if (LogLevel>0) printf("  ***> Parameter_filenameSuffix=%s\n",Parameter_filenameSuffix); 
  snprintf(outputFileNameExtension,600,"%s",READoutputFileNameExtension);
  if (LogLevel>0) printf("  ***> outputFileNameExtension=%s\n",outputFileNameExtension);   
  
  Parameter_L0 = READParameter_L0;
  if (LogLevel>0) printf("  ***> Parameter_L0=%d\n",Parameter_L0);
  Parameter_L1 = READParameter_L1;
  if (LogLevel>0) printf("  ***> Parameter_L1=%d\n",Parameter_L1);
  Parameter_L2 = READParameter_L2;
  if (LogLevel>0) printf("  ***> Parameter_L2=%d\n",Parameter_L2);
  Parameter_L3 = READParameter_L3;
  if (LogLevel>0) printf("  ***> Parameter_L3=%d\n",Parameter_L3);
  Parameter_Nf = READParameter_Nf;
  if (LogLevel>0) printf("  ***> Parameter_Nf=%d\n",Parameter_Nf);
  Parameter_RHO = READParameter_RHO;
  if (LogLevel>0) printf("  ***> Parameter_RHO=%f\n",Parameter_RHO);
  Parameter_R = READParameter_R;
  if (LogLevel>0) printf("  ***> Parameter_R=%f\n",Parameter_R);
  Parameter_Lambda = READParameter_Lambda;
  if (LogLevel>0) printf("  ***> Parameter_Lambda=%f\n",Parameter_Lambda);
  Parameter_Kappa = READParameter_Kappa;
  if (LogLevel>0) printf("  ***> Parameter_Kappa=%f\n",Parameter_Kappa);
  Parameter_Y = READParameter_Y;
  if (LogLevel>0) printf("  ***> Parameter_Y=%f\n",Parameter_Y);

  Parameter_MassSplit = READParameter_MassSplit;
  if (LogLevel>0) printf("  ***> Parameter_MassSplit=%f\n",Parameter_MassSplit);
  Parameter_ExplicitMass = READParameter_ExplicitMass;
  if (LogLevel>0) printf("  ***> Parameter_ExplicitMass=%f\n",Parameter_ExplicitMass);
  Parameter_ExplicitCurrent = READParameter_ExplicitCurrent;
  if (LogLevel>0) printf("  ***> Parameter_ExplicitCurrent=%f\n",Parameter_ExplicitCurrent);
    
  Parameter_c6 = READParameter_c6;
  if (LogLevel>0) printf("  ***> Parameter_c6=%f\n",Parameter_c6);
  Parameter_c8 = READParameter_c8;
  if (LogLevel>0) printf("  ***> Parameter_c8=%f\n",Parameter_c8);
  Parameter_c10 = READParameter_c10;
  if (LogLevel>0) printf("  ***> Parameter_c10=%f\n",Parameter_c10);
  Parameter_lambda6 = READParameter_lambda6;
  if (LogLevel>0) printf("  ***> Parameter_lambda6=%f\n",Parameter_lambda6);
  Parameter_lambda8 = READParameter_lambda8;
  if (LogLevel>0) printf("  ***> Parameter_lambda8=%f\n",Parameter_lambda8);
  Parameter_lambda10 = READParameter_lambda10;
  if (LogLevel>0) printf("  ***> Parameter_lambda10=%f\n",Parameter_lambda10);      
  
  for (I=0; I<1+SubPolynomMAX; I++) {
    Parameter_PolyEpsilon[I] = READParameter_PolyEpsilon[I];
    if (LogLevel>0) printf("  ***> Parameter_PolyEpsilon[%d]=%f\n",I,Parameter_PolyEpsilon[I]);
  }
  Parameter_PolyLambda = READParameter_PolyLambda;
  if (LogLevel>0) printf("  ***> Parameter_PolyLambda=%f\n",Parameter_PolyLambda);
  for (I=0; I<1+SubPolynomMAX; I++) {
    Parameter_PolyDegree[I] = READParameter_PolyDegree[I];
    if (LogLevel>0) printf("  ***> Parameter_PolyDegree[%d]=%d\n",I,Parameter_PolyDegree[I]);
  }
  Parameter_PolyAlpha = READParameter_PolyAlpha;
  if (LogLevel>0) printf("  ***> Parameter_PolyAlpha=%f\n",Parameter_PolyAlpha);
  Parameter_PolDigit = READParameter_PolDigit;
  if (LogLevel>0) printf("  ***> Parameter_PolDigit=%d\n",Parameter_PolDigit);
  Parameter_MaxPolDegPerNode = READParameter_MaxPolDegPerNode;
  if (LogLevel>0) printf("  ***> Parameter_MaxPolDegPerNode=%d\n",Parameter_MaxPolDegPerNode);

  Parameter_FLAG_RandomGauge = READParameter_FLAG_RandomGauge;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_RandomGauge=%d\n",Parameter_FLAG_RandomGauge);  
  Parameter_FLAG_ExactReweighing = READParameter_FLAG_ExactReweighing;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_ExactReweighing=%d\n",Parameter_FLAG_ExactReweighing);
  Parameter_FLAG_QuasiHermiteanMode = READParameter_FLAG_QuasiHermiteanMode;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_QuasiHermiteanMode=%d\n",Parameter_FLAG_QuasiHermiteanMode);
  Parameter_FLAG_FactorizationOfMatrixB = READParameter_FLAG_FactorizationOfMatrixB;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_FactorizationOfMatrixB=%d\n",Parameter_FLAG_FactorizationOfMatrixB);
  Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime = READParameter_FLAG_AntiPeriodicBoundaryConditionsInTime;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime=%d\n",Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime);

  Parameter_FLAG_Use_P_Preconditioner = READParameter_FLAG_Use_P_Preconditioner;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_Use_P_Preconditioner=%d\n",Parameter_FLAG_Use_P_Preconditioner);
  Parameter_FLAG_Use_Q_Preconditioner = READParameter_FLAG_Use_Q_Preconditioner;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_Use_Q_Preconditioner=%d\n",Parameter_FLAG_Use_Q_Preconditioner);
  Parameter_FLAG_Use_R_Preconditioner = READParameter_FLAG_Use_R_Preconditioner;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_Use_R_Preconditioner=%d\n",Parameter_FLAG_Use_R_Preconditioner);
  Parameter_QPreconditioner_Mu = READParameter_QPreconditioner_Mu;
  if (LogLevel>0) printf("  ***> Parameter_QPreconditioner_Mu=%f\n",Parameter_QPreconditioner_Mu);
  Parameter_QPreconditioner_Beta = READParameter_QPreconditioner_Beta;
  if (LogLevel>0) printf("  ***> Parameter_QPreconditioner_Beta=%f\n",Parameter_QPreconditioner_Beta);
  
  Parameter_UpperEWboundSafetyFactor = READParameter_UpperEWboundSafetyFactor;
  if (LogLevel>0) printf("  ***> Parameter_UpperEWboundSafetyFactor=%f\n",Parameter_UpperEWboundSafetyFactor);
  Parameter_FLAG_DirectOmegaSampling = READParameter_FLAG_DirectOmegaSampling;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_DirectOmegaSampling=%d\n",Parameter_FLAG_DirectOmegaSampling);
  
  Parameter_Theta = READParameter_Theta;
  if (LogLevel>0) printf("  ***> Parameter_Theta=%f\n",Parameter_Theta);
  Parameter_FLAG_ThetaScan = READParameter_FLAG_ThetaScan;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_ThetaScan=%d\n",Parameter_FLAG_ThetaScan);  
  Parameter_ThetaMin = READParameter_ThetaMin;
  if (LogLevel>0) printf("  ***> Parameter_ThetaMin=%f\n",Parameter_ThetaMin);
  Parameter_ThetaMax = READParameter_ThetaMax;
  if (LogLevel>0) printf("  ***> Parameter_ThetaMax=%f\n",Parameter_ThetaMax);
  Parameter_ThetaFactor = READParameter_ThetaFactor;
  if (LogLevel>0) printf("  ***> Parameter_ThetaFactor=%f\n",Parameter_ThetaFactor);
  
  Parameter_FLAG_ModelSelection = READParameter_FLAG_ModelSelection;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_ModelSelection=%d\n",Parameter_FLAG_ModelSelection);

  Parameter_Epsilon = READParameter_Epsilon;
  if (LogLevel>0) printf("  ***> Parameter_Epsilon=%f\n",Parameter_Epsilon);
  for (I=0; I<2+SubPolynomMAX; I++) {
    Parameter_IntegratorType[I] = READParameter_IntegratorType[I];
    if (LogLevel>0) printf("  ***> Parameter_IntegratorType[%d]=%d\n",I,Parameter_IntegratorType[I]);
  }
  Parameter_SubPolyCount = READParameter_SubPolyCount;
  if (LogLevel>0) printf("  ***> Parameter_SubPolyCount=%d\n",Parameter_SubPolyCount);
  for (I=0; I<2+SubPolynomMAX; I++) {
    Parameter_Iterations[I] = READParameter_Iterations[I];
    if (LogLevel>0) printf("  ***> Parameter_Iterations[%d]=%d\n",I,Parameter_Iterations[I]);
  }
  
  Parameter_Measurements = READParameter_Measurements;
  if (LogLevel>0) printf("  ***> Parameter_Measurements=%d\n",Parameter_Measurements);
  Parameter_AutomaticPreconditioningMetros = READParameter_AutomaticPreconditioningMetros;
  if (LogLevel>0) printf("  ***> Parameter_AutomaticPreconditioningMetros=%d\n",Parameter_AutomaticPreconditioningMetros);
  Parameter_ThermalizingMetros = READParameter_ThermalizingMetros;
  if (LogLevel>0) printf("  ***> Parameter_ThermalizingMetros=%d\n",Parameter_ThermalizingMetros);
  
  Parameter_SaveConfigurationEveryXXXResult = READParameter_SaveConfigurationEveryXXXResult;
  if (LogLevel>0) printf("  ***> Parameter_SaveConfigurationEveryXXXResult=%d\n",Parameter_SaveConfigurationEveryXXXResult);
  Parameter_SaveStateDescriptorEveryXXXResult = READParameter_SaveStateDescriptorEveryXXXResult;
  if (LogLevel>0) printf("  ***> Parameter_SaveStateDescriptorEveryXXXResult=%d\n",Parameter_SaveStateDescriptorEveryXXXResult);
  
  Parameter_DetermineEWeveryXXXconfsPREC = READParameter_DetermineEWeveryXXXconfsPREC;
  if (LogLevel>0) printf("  ***> Parameter_DetermineEWeveryXXXconfsPREC=%d\n",Parameter_DetermineEWeveryXXXconfsPREC);
  Parameter_DetermineEWeveryXXXconfsMEASURE = READParameter_DetermineEWeveryXXXconfsMEASURE;
  if (LogLevel>0) printf("  ***> Parameter_DetermineEWeveryXXXconfsMEASURE=%d\n",Parameter_DetermineEWeveryXXXconfsMEASURE);

  Parameter_FACC_type = READParameter_FACC_type;
  if (LogLevel>0) printf("  ***> Parameter_FACC_type=%d\n",Parameter_FACC_type);
  Parameter_FACC_parameter = READParameter_FACC_parameter;
  if (LogLevel>0) printf("  ***> Parameter_FACC_parameter=%f\n",Parameter_FACC_parameter);  
  Parameter_OmegaMassAdaptionMode = READParameter_OmegaMassAdaptionMode;
  if (LogLevel>0) printf("  ***> Parameter_OmegaMassAdaptionMode=%d\n",Parameter_OmegaMassAdaptionMode);
  Parameter_ReversibilityCheckFreqPrec = READParameter_ReversibilityCheckFreqPrec;
  if (LogLevel>0) printf("  ***> Parameter_ReversibilityCheckFreqPrec=%d\n",Parameter_ReversibilityCheckFreqPrec);
  Parameter_ReversibilityCheckFreqMeas = READParameter_ReversibilityCheckFreqMeas;
  if (LogLevel>0) printf("  ***> Parameter_ReversibilityCheckFreqMeas=%d\n",Parameter_ReversibilityCheckFreqMeas);
  Parameter_StartRandSeed = READParameter_StartRandSeed;
  if (LogLevel>0) printf("  ***> Parameter_StartRandSeed=%d\n",Parameter_StartRandSeed);
  
  phiConfIDnr = READphiConfIDnr;
  if (LogLevel>0) printf("  ***> phiConfIDnr=%d\n",phiConfIDnr);    
  Parameter_SphericalHiggsMode = READParameter_SphericalHiggsMode;
  if (LogLevel>0) printf("  ***> Parameter_SphericalHiggsMode=%d\n",Parameter_SphericalHiggsMode);    
  Parameter_SphericalHiggsZeta = READParameter_SphericalHiggsZeta;
  if (LogLevel>0) printf("  ***> Parameter_SphericalHiggsZeta=%f\n",Parameter_SphericalHiggsZeta);  


  //Setting Rand-Seed
  AdvancedSeed = -READAdvancedSeed;
  if (LogLevel>0) printf("  ***> AdvancedSeed=%d\n",AdvancedSeed);
  AdvancedZufall(AdvancedSeed);
  if (LogLevel>0) printf("  ***> AdvancedSeed=%d after ini, and first random=%f\n",AdvancedSeed, AdvancedZufall(AdvancedSeed));


  //Taking some parameters over from ini-file
  if (LogLevel>0) printf("Taking over values from Ini-File...\n");
  if (LogLevel>0) printf("  ***> Parameter_AdditionalAuxVectors=%d\n",Parameter_AdditionalAuxVectors);
  if (LogLevel>0) printf("  ***> Parameter_FLAG_xFFT=%d\n",Parameter_FLAG_xFFT);
  if (LogLevel>0) printf("  ***> Parameter_TotalData=%d\n",Parameter_TotalData);  
  if (LogLevel>0) printf("  ***> Parameter_MaxRunTime=%f\n",Parameter_MaxRunTime);    
  if (LogLevel>0) printf("  ***> LogLevel=%d\n",LogLevel);  


  //Reading Phi-Field 
  fgets(Comment, 500, file);
  fgets(Comment, 500, file);
  if (LogLevel>0) printf("\nReading Phi - Field:\n");
  
  int VL = fermiOps->getVectorLength();
  readError = -1;
  for (I=0; I<VL/8; I++) {
    if (fscanf(file,"%lf %lf %lf %lf",&(pHMCProp->phiField[I][0]),&(pHMCProp->phiField[I][1]),&(pHMCProp->phiField[I][2]),&(pHMCProp->phiField[I][3]))!=4) {
      readError = I;
      break;
    }
  }
  fclose(file);
  delete[] Comment;
  delete[] READParameter_filenameSuffix;
  delete[] READoutputFileNameExtension;
  
  if (readError>=0) {
    if (LogLevel>0) printf("Error (index %d) while reading Phi-Field from File-Descriptor. ==> Exiting\n",readError);
    desini();
    exit(0);  
  }

  //Consistency-Check on Phi-Field
  double measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm;  
  pHMCProp->measure(measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm);
  if (LogLevel>0) printf("Configuration read with mag=%1.15f and stag. mag=%1.15f\n",measurePhiNorm,measureStaggeredPhiNorm);
  if (LogLevel>0) printf("Old result              mag=%1.15f and stag. mag=%1.15f\n",READmeasurePhiNorm,READmeasureStaggeredPhiNorm); 
  
  if ((abs(READmeasurePhiNorm-measurePhiNorm)>ceps) || (abs(READmeasureStaggeredPhiNorm-measureStaggeredPhiNorm)>ceps)) {
    if (LogLevel>0) printf("Consistency-Error for read Phi-Field from File-Descriptor. ==> Exiting\n");
    desini();
    exit(0);  
  }
  return true;
}


void iniPerformanceProfiler() {
  char* fileName = new char[2000];
  snprintf(fileName, 2000, "%s/data/results/pHMC/miscellaneous/PerformanceProfile%s.dat", DataBaseDirectory,outputFileNameExtension);
  initializePerformanceProfiler(fileName);  
  delete[] fileName;
}


void readFourierAccelerationData() {
  char* fileName1 = new char[600]; 
  char* fileName2 = new char[600]; 
  char* fileName3 = new char[600]; 
  char* fileName4 = new char[600]; 
  char* fileName5 = new char[600]; 
  char* fileName6 = new char[600]; 
  char* fileName7 = new char[600]; 
  char* fileName8 = new char[600]; 
  char* fileName9 = new char[600]; 
  char* fileName10 = new char[600]; 
  snprintf(fileName1,600,"%s/data/results/pHMC/FACC/PhiTotalForcePREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName2,600,"%s/data/results/pHMC/FACC/PhiTotalForceGLOBAL%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName3,600,"%s/data/results/pHMC/FACC/PhiHiggsForcePREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName4,600,"%s/data/results/pHMC/FACC/PhiHiggsForceGLOBAL%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName5,600,"%s/data/results/pHMC/FACC/PhiFermionForcePREC%s",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName6,600,"%s/data/results/pHMC/FACC/PhiFermionForceGLOBAL%s",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName7,600,"%s/data/results/pHMC/FACC/PhiChangeNoFACC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName8,600,"%s/data/results/pHMC/FACC/PhiChangeFACC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName9,600,"%s/data/results/pHMC/FACC/PhiTotalSphericalProjectedForcePREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName10,600,"%s/data/results/pHMC/FACC/PhiTotalSphericalProjectedForceGLOBAL%s.dat",DataBaseDirectory,outputFileNameExtension);

  
  pHMCProp->loadPhiTotalForceFourierComponentsPREC(fileName1);
  pHMCProp->loadPhiTotalForceFourierComponentsGLOBAL(fileName2);
  pHMCProp->loadPhiHiggsForceFourierComponentsPREC(fileName3);
  pHMCProp->loadPhiHiggsForceFourierComponentsGLOBAL(fileName4);
  pHMCProp->loadPhiFermionForceFourierComponentsPREC(fileName5);
  pHMCProp->loadPhiFermionForceFourierComponentsGLOBAL(fileName6);
  pHMCProp->loadPhiChangeFourierComponentsNoFACC(fileName7);   
  pHMCProp->loadPhiChangeFourierComponentsFACC(fileName8);   
  pHMCProp->loadPhiTotalSphericalProjectedForceFourierComponentsPREC(fileName9);   
  pHMCProp->loadPhiTotalSphericalProjectedForceFourierComponentsGLOBAL(fileName10);   
  
  
  delete[] fileName1;
  delete[] fileName2;
  delete[] fileName3;
  delete[] fileName4;
  delete[] fileName5;
  delete[] fileName6;
  delete[] fileName7;
  delete[] fileName8;
  delete[] fileName9;
  delete[] fileName10;
}


void writeFourierAccelerationData(bool prec) {
  char* fileName1 = new char[600]; 
  char* fileName2 = new char[600]; 
  char* fileName3 = new char[600]; 
  char* fileName4 = new char[600]; 
  char* fileName5 = new char[600]; 
  char* fileName6 = new char[600]; 
  char* fileName7 = new char[600]; 
  char* fileName8 = new char[600]; 
  char* fileName9 = new char[600]; 
  char* fileName10 = new char[600]; 
  char* fileName11 = new char[600]; 
  snprintf(fileName1,600,"%s/data/results/pHMC/FACC/PhiTotalForcePREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName2,600,"%s/data/results/pHMC/FACC/PhiTotalForceGLOBAL%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName3,600,"%s/data/results/pHMC/FACC/PhiHiggsForcePREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName4,600,"%s/data/results/pHMC/FACC/PhiHiggsForceGLOBAL%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName5,600,"%s/data/results/pHMC/FACC/PhiFermionForcePREC%s",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName6,600,"%s/data/results/pHMC/FACC/PhiFermionForceGLOBAL%s",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName7,600,"%s/data/results/pHMC/FACC/PhiChangeNoFACC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName8,600,"%s/data/results/pHMC/FACC/PhiChangeFACC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName9,600,"%s/data/results/pHMC/FACC/MomentumMasses%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName10,600,"%s/data/results/pHMC/FACC/PhiTotalSphericalProjectedForcePREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  snprintf(fileName11,600,"%s/data/results/pHMC/FACC/PhiTotalSphericalProjectedForceGLOBAL%s.dat",DataBaseDirectory,outputFileNameExtension);
    
  if (prec) {
    pHMCProp->savePhiTotalForceFourierComponentsPREC(fileName1);
    pHMCProp->savePhiHiggsForceFourierComponentsPREC(fileName3);
    pHMCProp->savePhiFermionForceFourierComponentsPREC(fileName5);
    pHMCProp->savePhiTotalSphericalProjectedForceFourierComponentsPREC(fileName10);   
  }
  pHMCProp->savePhiTotalForceFourierComponentsGLOBAL(fileName2);
  pHMCProp->savePhiHiggsForceFourierComponentsGLOBAL(fileName4);
  pHMCProp->savePhiFermionForceFourierComponentsGLOBAL(fileName6);
  pHMCProp->savePhiChangeFourierComponentsNoFACC(fileName7);   
  pHMCProp->savePhiChangeFourierComponentsFACC(fileName8);   
  pHMCProp->writeMomentumMasses(fileName9);
  pHMCProp->savePhiTotalSphericalProjectedForceFourierComponentsGLOBAL(fileName11);   
  
  delete[] fileName1;
  delete[] fileName2;
  delete[] fileName3;
  delete[] fileName4;
  delete[] fileName5;
  delete[] fileName6;
  delete[] fileName7;
  delete[] fileName8;
  delete[] fileName9;
  delete[] fileName10;
  delete[] fileName11;
}


void readUpperEWboundLogData() {
  char* fileName1 = new char[600]; 
  snprintf(fileName1,600,"%s/data/results/pHMC/miscellaneous/UpperEWboundLogPREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  pHMCProp->loadUpperEWboundLogFromDisk(fileName1);
  delete[] fileName1;
}


void writeUpperEWboundLogData() {
  char* fileName1 = new char[600]; 
  snprintf(fileName1,600,"%s/data/results/pHMC/miscellaneous/UpperEWboundLogPREC%s.dat",DataBaseDirectory,outputFileNameExtension);
  pHMCProp->saveUpperEWboundLogToDisk(fileName1);
  delete[] fileName1;
}


int howManyMatrixApplicationsForIntegrator(double* iter) {
  int howMany = 0;
  
  if (Parameter_Y>0) {
    for (int I=SubPolynomMAX; I>=0; I--) {
      double add = Parameter_PolyDegree[I];
    
      for (int I2=I; I2>=0; I2--) {
        add *= iter[I2];
        if (Parameter_IntegratorType[I2]==1) add *= 2;
        if (Parameter_IntegratorType[I2]==2) add *= 5;
      }

      howMany += roundToInt(add);
    }
  } else {
    double add = 1;
    
    for (int I2=SubPolynomMAX+1; I2>=0; I2--) if (iter[I2]>0) {
      add *= iter[I2];
      if (Parameter_IntegratorType[I2]==1) add *= 2;
      if (Parameter_IntegratorType[I2]==2) add *= 5;
    }

    howMany += roundToInt(add);
  }
  
  return howMany;
}


int howManyMatrixApplicationsForDirectOmegaSampling() {
  int howMany = 0;
  if (Parameter_FLAG_DirectOmegaSampling>0) {
    howMany = 3 * Parameter_PolyDegree[0]; 
  }
  if (Parameter_Y==0) howMany = 0;  
  return howMany;
}


int howManyMatrixApplicationsForExactReweighing() {
  int howMany = 0;  
  if (Parameter_FLAG_ExactReweighing > 0) {
    howMany = 3 * Parameter_PolyDegree[0]; 
  }
  if (Parameter_Y==0) howMany = 0;  
  return howMany;
}


bool finalIntegratorSelection() {
#ifdef IntegrationSchemeAdaption
  if (LogLevel>1) printf("  *** Integrator Selection from Data Pool ***\n");
  double** integratorTableEval = new double*[ThetaIterIntegratorTableMAX];
  int tableEntries = 0;
  double ceps = 1E-9;
  
  for (int I=0; I<ThetaIterIntegratorTableCount; I++) {
    if (abs(ThetaIterIntegratorTable[I][6+SubPolynomMAX]-1)<ceps) {
      bool found = false;
      for (int I2=0; I2<tableEntries; I2++) {
        found = true;
	for (int I3=0; I3<3+SubPolynomMAX; I3++) {
          if (abs(integratorTableEval[I2][2+I3]-ThetaIterIntegratorTable[I][2+I3])>ceps) found = false;
	}
	if (found) {
	  integratorTableEval[I2][0]++;
	  integratorTableEval[I2][1]+= ThetaIterIntegratorTable[I][7+SubPolynomMAX];
	  break;
	}
      }
      if (!found) {
        integratorTableEval[tableEntries] = new double[5+SubPolynomMAX];
	integratorTableEval[tableEntries][0] = 1;
	integratorTableEval[tableEntries][1] = ThetaIterIntegratorTable[I][7+SubPolynomMAX];
	for (int I3=0; I3<3+SubPolynomMAX; I3++) {
          integratorTableEval[tableEntries][2+I3] = ThetaIterIntegratorTable[I][2+I3];
	}
	tableEntries++;
      }
    }
  }
  
  FILE* file;
  char* fileName = new char[600]; 
  snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/IntegratorSelection%s.dat",DataBaseDirectory,outputFileNameExtension);
  file = fopen(fileName,"w");
  delete[] fileName;
  double bestSpeed = 0;
  int bestInd = -1;
  for (int I=0; I<tableEntries; I++) {
    int thetaLevel = roundToInt(log(integratorTableEval[I][4+SubPolynomMAX]/Parameter_ThetaMin) / log(Parameter_ThetaFactor));
    if ((integratorTableEval[I][4+SubPolynomMAX]==0) || (Parameter_ThetaMin==0) || (Parameter_ThetaFactor==0) || (Parameter_FLAG_ThetaScan==0)) {
      thetaLevel = 0;
    }
    if (integratorTableEval[I][0] < 20) {
      integratorTableEval[I][0] = 0;
      integratorTableEval[I][1] = 1;
    } else {
      integratorTableEval[I][0] = integratorTableEval[I][1] / integratorTableEval[I][0];
      integratorTableEval[I][1] = howManyMatrixApplicationsForDirectOmegaSampling() + howManyMatrixApplicationsForExactReweighing() + howManyMatrixApplicationsForIntegrator(&(integratorTableEval[I][2]));
      
      if (LogLevel>2) printf("Tracked Integrator setting: theta = %f (L%d), iter = ",integratorTableEval[I][4+SubPolynomMAX],thetaLevel);
      fprintf(file,"Tracked Integrator setting: theta = %f (L%d), iter = ",integratorTableEval[I][4+SubPolynomMAX],thetaLevel);
      for (int I2=0; I2<2+SubPolynomMAX; I2++) {
        if (LogLevel>2) printf("%d ", roundToInt(integratorTableEval[I][2+I2]));
        fprintf(file,"%d ", roundToInt(integratorTableEval[I][2+I2]));
      }
      double speed = integratorTableEval[I][0] / integratorTableEval[I][1];
      if (LogLevel>2) printf("(Higgs) with acceptRate=%f and %d=%d+%d+%d matApp per conf ==> Speed: %f.\n\n",integratorTableEval[I][0],roundToInt(integratorTableEval[I][1]),roundToInt(howManyMatrixApplicationsForDirectOmegaSampling()), roundToInt(howManyMatrixApplicationsForExactReweighing()), roundToInt(integratorTableEval[I][1]-howManyMatrixApplicationsForDirectOmegaSampling()- howManyMatrixApplicationsForExactReweighing()), speed);
      fprintf(file,"(Higgs) with acceptRate=%f and %d=%d+%d+%d matApp per conf ==> Speed: %f.\n\n",integratorTableEval[I][0], roundToInt(integratorTableEval[I][1]), roundToInt(howManyMatrixApplicationsForDirectOmegaSampling()), roundToInt(howManyMatrixApplicationsForExactReweighing()), roundToInt(integratorTableEval[I][1]-howManyMatrixApplicationsForDirectOmegaSampling()- howManyMatrixApplicationsForExactReweighing()), speed);
    }
    if (integratorTableEval[I][0]<AutomaticAdaption_LowPro) {
      integratorTableEval[I][0] = 0;
      integratorTableEval[I][1] = 1;
    }
    double xtrThetaFac = exp(0.75*log(Parameter_ThetaFactor)*thetaLevel);
    if (xtrThetaFac*integratorTableEval[I][0]/integratorTableEval[I][1]-ceps > bestSpeed) {
      bestSpeed = xtrThetaFac*integratorTableEval[I][0]/integratorTableEval[I][1];
      bestInd = I;
    }
  }
  
  bool bestFound = true;
  if (bestInd==-1) bestFound = false;
  
  if (bestFound) {
    Parameter_Epsilon *= Parameter_Iterations[0];
    for (int I=0; I<2+SubPolynomMAX; I++) {
      Parameter_Iterations[I] = roundToInt(integratorTableEval[bestInd][2+I]);
    }
    Parameter_Epsilon /= Parameter_Iterations[0];
    if (Parameter_FLAG_ThetaScan>0) {
      Parameter_Theta = integratorTableEval[bestInd][4+SubPolynomMAX];
      pHMCProp->setTheta(Parameter_Theta);   
    }
    double accRate = integratorTableEval[bestInd][0];
    int matApp = roundToInt(integratorTableEval[bestInd][1]);
    if (LogLevel>1) printf("Best Integrator setting: theta = %f, eps=%f, iter = ",Parameter_Theta,Parameter_Epsilon);
    fprintf(file, "Best Integrator setting: theta = %f, eps=%f, iter = ",Parameter_Theta,Parameter_Epsilon);
    for (int I=0; I<2+SubPolynomMAX; I++) {
      if (LogLevel>1) printf("%d ", Parameter_Iterations[I]);
      fprintf(file, "%d ", Parameter_Iterations[I]);
    }
    double speed = accRate/matApp;
    if (LogLevel>1) printf("(Higgs) with acceptRate=%f and %d=%d+%d+%d matApp per conf ==> Speed: %f.\n\n",accRate,matApp, roundToInt(howManyMatrixApplicationsForDirectOmegaSampling()), roundToInt(howManyMatrixApplicationsForExactReweighing()), roundToInt(matApp-howManyMatrixApplicationsForDirectOmegaSampling()- howManyMatrixApplicationsForExactReweighing()), speed);
    fprintf(file,"(Higgs) with acceptRate=%f and %d=%d+%d+%d matApp per conf ==> Speed: %f.\n\n",accRate,matApp, roundToInt(howManyMatrixApplicationsForDirectOmegaSampling()), roundToInt(howManyMatrixApplicationsForExactReweighing()), roundToInt(matApp-howManyMatrixApplicationsForDirectOmegaSampling()- howManyMatrixApplicationsForExactReweighing()), speed);
  } else {
    if (LogLevel>1) printf("  ***WARNING*** No acceptable integrator scheme found. Leaving old one.\n");
    fprintf(file,"  ***WARNING*** No acceptable integrator scheme found. Leaving old one.\n");
  }
  fclose(file);
  
  for (int I=0; I<tableEntries; I++) delete[] integratorTableEval[I];
  delete[] integratorTableEval;
  return bestFound;
#else
  if (LogLevel>1) printf("  *** Integrator Selection deactivated ***\n");
  return true;
#endif
}


void addThetaIterIntegratorTableEntry(int accepted, bool afterMassDet) {  
  ThetaIterIntegratorTable[ThetaIterIntegratorTableCount] = new double[6+2+SubPolynomMAX];
  ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][0] = AutomaticPreconRunCount;
  ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][1] = Parameter_Epsilon;
  for (int I=0; I<2+SubPolynomMAX; I++) {
    ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][2+I] = Parameter_Iterations[I];
  }
  ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][4+SubPolynomMAX] = Parameter_Theta;
  ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][5+SubPolynomMAX] = pHMCProp->SafterProp-pHMCProp->SbeforeProp;  
  ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][6+SubPolynomMAX] = afterMassDet;
  ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][7+SubPolynomMAX] = accepted;
  
  
  FILE* file;
  char* fileName = new char[600]; 
  snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/ThetaIterIntegratorTable%s.dat",DataBaseDirectory,outputFileNameExtension);
  file = fopen(fileName,"a");
  delete[] fileName;
  
  fprintf(file,"%d %f ",(int)(ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][0]), ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][1]);                                          
  for (int I=0; I<2+SubPolynomMAX; I++) {
    fprintf(file,"%d ",(int)(ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][2+I]));                                          
  }
  fprintf(file,"%f ",ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][4+SubPolynomMAX]);    
  if (ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][5+SubPolynomMAX]>=0) {
    fprintf(file," "); 
  }
  fprintf(file,"%f ",ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][5+SubPolynomMAX]);    
  fprintf(file,"%d ",(int)(ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][6+SubPolynomMAX]));    
  fprintf(file,"%d\n",(int)(ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][7+SubPolynomMAX]));    

  fclose(file);  
  
  ThetaIterIntegratorTableCount++;
}


void readThetaIterIntegratorTable() {  
  FILE* file;
  char* fileName = new char[600]; 
  snprintf(fileName,600,"%s/data/results/pHMC/miscellaneous/ThetaIterIntegratorTable%s.dat",DataBaseDirectory,outputFileNameExtension);
  file = fopen(fileName,"r");

  ThetaIterIntegratorTableCount = 0;

  if (file != NULL) {
    int i0;
    double d0;
    bool ok = true;
  
    while (ok) {  
      ThetaIterIntegratorTable[ThetaIterIntegratorTableCount] = new double[6+2+SubPolynomMAX];
      if (fscanf(file,"%d ",&i0) != 1) ok = false;
      ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][0] = i0;
      if (fscanf(file,"%lf ",&d0) != 1) ok = false;
      ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][1] = d0;
      for (int I=0; I<2+SubPolynomMAX; I++) {
        if (fscanf(file,"%d ",&i0) != 1) ok = false;
        ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][2+I] = i0;
      }
      if (fscanf(file,"%lf ",&d0) != 1) ok = false;
      ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][4+SubPolynomMAX] = d0;
      if (fscanf(file,"%lf ",&d0) != 1) ok = false;
      ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][5+SubPolynomMAX] = d0;
      if (fscanf(file,"%d ",&i0) != 1) ok = false;
      ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][6+SubPolynomMAX] = i0;
      if (fscanf(file,"%d\n",&i0) != 1) ok = false;
      ThetaIterIntegratorTable[ThetaIterIntegratorTableCount][7+SubPolynomMAX] = i0;

      if (!ok) {
        delete[] ThetaIterIntegratorTable[ThetaIterIntegratorTableCount];
      } else {
        ThetaIterIntegratorTableCount++;
      }
    }  
    fclose(file);
  }
  
  if (LogLevel>2) printf("%d entries of ThetaIterIntegratorTable read from file %s\n",ThetaIterIntegratorTableCount,fileName);
  delete[] fileName;
}


void automaticBackup() {
  if (LogLevel>2) printf("Automatic Backup...\n");
  char* sourceFileName = new char[600]; 
  snprintf(sourceFileName,600,"%s/data/results/pHMC/states/StateDescriptor%s.dat",DataBaseDirectory,outputFileNameExtension);
  char* destFileName = new char[600]; 
  snprintf(destFileName,600,"%s/data/results/pHMC/states/BACKUP/StateDescriptor%s.dat",DataBaseDirectory,outputFileNameExtension);

  copyFile(sourceFileName, destFileName);

  delete[] sourceFileName;
  delete[] destFileName;
  if (LogLevel>2) printf("...Automatic Backup complete!\n\n");
}



int main(int argc,char **argv) {
  double acceptRate, RunCount, acceptCount;
  int ResampleCount, MeasureCount;
  int ranIniValue = 1;
  int I;

  iniMPI(argc, argv);
  loadParameters(0,false);  
  if ((Parameter_FLAG_DirectOmegaSampling == 0) && (Parameter_Theta == 0)) {
    printf("ERROR: Direct omega sampling is off and theta = 0!!!\n");
    exit(0);
  }
  if (LogLevel>1) printf("Initialisiere FFTW with %d threads.\n",Parameter_FFTWThreadCount);
  fftw_init_threads();
  fftw_plan_with_nthreads(Parameter_FFTWThreadCount);  

  startTimer();
  if (LogLevel>1) printf("Timer started. Passed time: %f\n",timePassed());
  
  if (LogLevel>1) printf("Number of arguments = %d\n",argc);
  for (I=0; I<argc; I++) {
    if (LogLevel>1) printf("Argument %d: %s\n",I+1,argv[I]);  
  }
  
  JobNr = -1;
  if (argc>=2) {
    if (sscanf(argv[1],"%d",&JobNr)!=1)  JobNr = -1;
  }
  int MultiProcessNr = -1;
  if (argc>=3) {
    if (sscanf(argv[2],"%d",&MultiProcessNr)!=1)  MultiProcessNr= -1;
  }
  
  if (JobNr>0) {
    ranIniValue = JobNr * (ownNodeID+1);
    if (LogLevel>1) printf("Use Job-Nr*(ownID+1) =  %d*(%d+1) = %d for rand seed generation.\n",JobNr,ownNodeID, ranIniValue);
  } else {
    ranIniValue = 5517*(ownNodeID+1);
    if (LogLevel>1) printf("Start rand seed from system time only.\n");
    JobNr = 1;
  }
  if (MultiProcessNr>0) {
    if (LogLevel>1) printf("Use multi-start number: %d for simulation parameter determination.\n",MultiProcessNr);
  } else {
    if (LogLevel>1) printf("Single Job-start.\n");
    MultiProcessNr = 0;
  }
  if (LogLevel>0) printf("\n");
  
  iniTools(ranIniValue);
  int threadCountPerNode = Parameter_FFTWThreadCount;
  if (Parameter_FLAG_xFFT==1) threadCountPerNode = Parameter_xFFTThreadCount;
  char* fftPlanDescriptor = NULL;
  bool fftPlanAvail = readOptimalFermionVectorEmbeddingAndFFTPlanFromTuningDB(Parameter_L0, Parameter_L1, Parameter_L2, Parameter_L3, threadCountPerNode, Parameter_MultiThreadedOpMode, Parameter_FLAG_xFFT, Parameter_FLAG_Use_P_Preconditioner, Parameter_FLAG_Use_Q_Preconditioner, Parameter_FLAG_Use_R_Preconditioner, Parameter_FLAG_QuasiHermiteanMode, fftPlanDescriptor);

  loadParameters(MultiProcessNr,true);
  if (Parameter_StartRandSeed>0) {
    if (LogLevel>0) printf("\n\n *** WARNING ***  Start value for randseed set to %d\n\n",Parameter_StartRandSeed);  
    AdvancedSeed  = -Parameter_StartRandSeed;
    AdvancedZufall(AdvancedSeed);
  }
  buildOutputFileNameExtension();
  

  fermiOps = new FermionMatrixOperations(Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3, Parameter_RHO, Parameter_R, Parameter_Y);
  fermiOps->setMassSplitRatio(Parameter_MassSplit);
  fermiOps->setExplicitMass(Parameter_ExplicitMass);
  fermiOps->setAntiPeriodicBoundaryConditionsInTime(Parameter_FLAG_AntiPeriodicBoundaryConditionsInTime);

  if ((Parameter_MultiThreadedOpMode>=0) && (Parameter_MultiThreadedOpMode<=2)) {
    fermiOps->activateMultiThreadedOps(Parameter_MultiThreadedOpMode, false);  
  }
  if (Parameter_FLAG_xFFT==1) {
    fermiOps->setxFFTusage(true);
    fermiOps->setxFFT_DistributedFFT_ThreadCount(Parameter_xFFTThreadCount);
    if (fftPlanAvail) {
      bool b = fermiOps->setDistributedFFTPlan(fftPlanDescriptor);
      if (!b) {
        printf("ERROR could not set FFT-Plan: %s\n", fftPlanDescriptor);
	exit(0);
      }
    } else {
      fermiOps->tuneDistributedFFT(ExtremeFFT4D_TuneLevel_Low);
    }
    
    fermiOps->testFourierTrafo(true);
    fermiOps->testFourierTrafo(false);
    if ((Parameter_MultiThreadedOpMode>=0) && (Parameter_MultiThreadedOpMode<=2)) {  
      fermiOps->testDistributedFourierTrafo(true);
      fermiOps->testDistributedFourierTrafo(false);
    }
  } else {
    fermiOps->setxFFTusage(false);
  }
  delete[] fftPlanDescriptor;
  
  double* polLam = new double[1+Parameter_SubPolyCount];
  for (I=0; I<1+Parameter_SubPolyCount; I++) polLam[I] = Parameter_PolyLambda;
  pHMCProp = new pHMCPropagator(fermiOps, Parameter_Lambda, Parameter_Kappa, Parameter_ExplicitCurrent, Parameter_c6, Parameter_c8, Parameter_c10, Parameter_lambda6, Parameter_lambda8, Parameter_lambda10, Parameter_Nf, 1.0, Parameter_SphericalHiggsMode, Parameter_SphericalHiggsZeta, Parameter_Theta, Parameter_SubPolyCount, Parameter_PolyEpsilon, polLam, Parameter_PolyDegree, 0, NULL, Parameter_PolDigit, Parameter_PolyAlpha, Parameter_MaxPolDegPerNode,Parameter_AdditionalAuxVectors);
  pHMCProp->activateForceStoring(true);
    
 
  if (ownNodeID>0) {
    pHMCProp->SlaveController();
    if (LogLevel>1) printf("pHMC Slave finished!!!  ==>  EXITING !!!\n");
    exit(0);
  }

  
  bool descriptorFileRead = readCurrentStateDescriptor();
  pHMCProp->setTheta(Parameter_Theta);  
  pHMCProp->setPhiForceFourierType(Parameter_FACC_type,Parameter_FACC_parameter);  
  pHMCProp->setOmegaMassAdaptionMode(Parameter_OmegaMassAdaptionMode);
  pHMCProp->resetExactMMdagInverseSQRTOmegaAction();
  pHMCProp->synchronizedChangeOfQuasiHermiteanMode(Parameter_FLAG_QuasiHermiteanMode);
  pHMCProp->synchronizedChangeOfModelSelection(Parameter_FLAG_ModelSelection);
  pHMCProp->synchronizedChangeOfBMatrixFactorizationMode(Parameter_FLAG_FactorizationOfMatrixB);

  pHMCProp->getNodesReady();

  if (descriptorFileRead) pHMCProp->readAllOmegaFieldsFromDisk();
  writeCurrentStateDescriptor(1, true);
  if (Parameter_FACC_type>0) {
    readFourierAccelerationData();
  }
  readUpperEWboundLogData();
  if (Parameter_OmegaMassAdaptionMode>0) {
    pHMCProp->readOmegaForceStrengthsFromDiskPREC();
    pHMCProp->readOmegaForceStrengthsFromDiskGLOBAL();
  }
  readThetaIterIntegratorTable();
  
  if (LogLevel>2) iniPerformanceProfiler();

  if ((Parameter_AutomaticPreconditioningMetros>0) && (AutomaticPreconRunCount<Parameter_AutomaticPreconditioningMetros)) {
    double swapFrac  = 0.10;
    int FACCmassDetTime = (int)(0.25*Parameter_AutomaticPreconditioningMetros);
    if (AutomaticPreconRunCount==0) {
      fermiOps->setPreconditioner(Parameter_FLAG_Use_P_Preconditioner, sqrt(Parameter_Nf), 0);
      pHMCProp->synchronizedChangeOfQPreconsitionerData(false, Parameter_QPreconditioner_Mu, Parameter_QPreconditioner_Beta);      
      pHMCProp->synchronizedChangeOfRPreconsitionerData(Parameter_FLAG_Use_R_Preconditioner, sqrt(Parameter_Nf), 1.0);
      pHMCProp->setPhiForceFourierType(0,0);
      pHMCProp->setOmegaMassAdaptionMode(0);
      
      firstNMarkovSteps(1,1,20,5);
      firstNMarkovSteps(1,1,10,4);
      firstNMarkovSteps(1,1,5,3);
      firstNMarkovSteps(5,1,2,2);
      firstNMarkovSteps(5,2,2,1.5);
    }
    pHMCProp->synchronizedChangeOfQPreconsitionerData(Parameter_FLAG_Use_Q_Preconditioner, Parameter_QPreconditioner_Mu, Parameter_QPreconditioner_Beta);      
    if (LogLevel>1) printf("Performing %d Automatic Preconditioning Metropolis steps...\n",Parameter_AutomaticPreconditioningMetros-AutomaticPreconRunCount);
    if ((Parameter_FACC_type>0) && (AutomaticPreconRunCount>=FACCmassDetTime)) {
      pHMCProp->calcPhiMomentumMasses(Parameter_Epsilon*Parameter_Iterations[0]);    
    }
    if ((Parameter_OmegaMassAdaptionMode>0) && (AutomaticPreconRunCount>=FACCmassDetTime)) {
      pHMCProp->calcOmegaMassAdaption();    
    }
    if (AutomaticPreconRunCount<(FACCmassDetTime/2)) {
      pHMCProp->setPhiForceFourierType(0,0);
      pHMCProp->setOmegaMassAdaptionMode(0);	
    } else {
      pHMCProp->setPhiForceFourierType(Parameter_FACC_type,Parameter_FACC_parameter);	
      pHMCProp->setOmegaMassAdaptionMode(Parameter_OmegaMassAdaptionMode);
    }
    
    for (ResampleCount=0; true; ResampleCount++) {
      RunCount = 0;
      acceptCount = 0;
      
      for (MeasureCount=0; MeasureCount<Parameter_Measurements; MeasureCount++) {
        if (AutomaticPreconRunCount==(FACCmassDetTime/2)) {
          pHMCProp->setPhiForceFourierType(Parameter_FACC_type,Parameter_FACC_parameter);	
          pHMCProp->setOmegaMassAdaptionMode(Parameter_OmegaMassAdaptionMode);
	}
        pHMCProp->setTheta(Parameter_Theta);        
        pHMCProp->sampleALLMomenta();
        if (Parameter_FLAG_DirectOmegaSampling > 0) pHMCProp->sampleOmegaFields();
        if (markovStep()) {
          double exactWeight = NaN;
          if (Parameter_FLAG_ExactReweighing>0) {
            pHMCProp->calcExactMMdagInverseSQRTOmegaAction();
            exactWeight = pHMCProp->getExactReweighingFactorFromMMdagInverseSQRTOmegaAction();
            writeWeightFactor(AutomaticPreconRunCount, 0, exactWeight);
          }
          if (AutomaticPreconRunCount>=(FACCmassDetTime/2)) {
  	    pHMCProp->analyzeAllPolynomialForces();
	  }
	  if ((Parameter_DetermineEWeveryXXXconfsPREC>0) && (AutomaticPreconRunCount % Parameter_DetermineEWeveryXXXconfsPREC == 0)) {
            writeConditionNumber(AutomaticPreconRunCount, 0);
	  }
	  if ((Parameter_FACC_type>0) && (AutomaticPreconRunCount>=FACCmassDetTime) && (AutomaticPreconRunCount % FACCmassDetTime == 0)) {
	    pHMCProp->calcPhiMomentumMasses(Parameter_Epsilon*Parameter_Iterations[0]);     
	  }
	  if ((Parameter_OmegaMassAdaptionMode>0) && (AutomaticPreconRunCount>=FACCmassDetTime) && (AutomaticPreconRunCount % FACCmassDetTime == 0)) {
	    pHMCProp->calcOmegaMassAdaption(); 
	    if ((AutomaticPreconRunCount==FACCmassDetTime) && (Parameter_FLAG_ThetaScan>0)) Parameter_Theta = Parameter_ThetaMin;
	  }
	  	
          addThetaIterIntegratorTableEntry(1, (AutomaticPreconRunCount>=FACCmassDetTime));

          acceptCount++;
          AutomaticPreconRunCount++;
          writeCurrentStateDescriptor(1, false);
          if (Parameter_FACC_type>0) {
            writeFourierAccelerationData(true);	 
	  }
          writeUpperEWboundLogData();
          if (Parameter_OmegaMassAdaptionMode>0) {
	    pHMCProp->writeOmegaForceStrengthsToDiskPREC();    
	    pHMCProp->writeOmegaForceStrengthsToDiskGLOBAL();    
	  }
	  if ((Parameter_ReversibilityCheckFreqPrec>0) && (((AutomaticPreconRunCount+1) % Parameter_ReversibilityCheckFreqPrec) == 0)) {
    	    writeReversibilityNorm(0);	  
	  }	
	  if (AutomaticPreconRunCount<swapFrac*Parameter_AutomaticPreconditioningMetros) {
            if (!Parameter_FLAG_QuasiHermiteanMode) pHMCProp->improvePreconditioningParametersFAST();
	    if (Parameter_FLAG_QuasiHermiteanMode) pHMCProp->improveRPreconditioningParameters(AutomaticPreconRunCount, 2, Parameter_PolyLambda, Parameter_UpperEWboundSafetyFactor);
	  } else {
	    if ((Parameter_FLAG_QuasiHermiteanMode) && (MeasureCount<Parameter_Measurements-1)) {
	      if (AutomaticPreconRunCount<2*FACCmassDetTime) {
 	        pHMCProp->improveRPreconditioningParameters(AutomaticPreconRunCount, 1, Parameter_PolyLambda, Parameter_UpperEWboundSafetyFactor);
	      } else {
 	        pHMCProp->improveRPreconditioningParameters(AutomaticPreconRunCount, 0, Parameter_PolyLambda, Parameter_UpperEWboundSafetyFactor);
	      }
	    }
          }
        } else {
          addThetaIterIntegratorTableEntry(0, (AutomaticPreconRunCount>=FACCmassDetTime));	
	}
        RunCount++;
      }
      writePerformanceProfilingDataToDisk();
      
      if (AutomaticPreconRunCount>=swapFrac*Parameter_AutomaticPreconditioningMetros) {
        if (!Parameter_FLAG_QuasiHermiteanMode) pHMCProp->improvePreconditioningParametersFAST();
	if (Parameter_FLAG_QuasiHermiteanMode) pHMCProp->improveRPreconditioningParameters(AutomaticPreconRunCount, 2, Parameter_PolyLambda, Parameter_UpperEWboundSafetyFactor);
      }
      acceptRate = acceptCount / RunCount;
      if (LogLevel>2) printf("Average (preconditioning step: %d) update quote: %1.2f\n",AutomaticPreconRunCount,acceptRate);
      
      if (AutomaticPreconRunCount>Parameter_AutomaticPreconditioningMetros-3*Parameter_Measurements) {
        finalIntegratorSelection();
      } else {
        automaticEpsilonAdaption(acceptRate, AutomaticPreconRunCount>=FACCmassDetTime);
      }
      pHMCProp->setTheta(Parameter_Theta);        

      if (AutomaticPreconRunCount>Parameter_AutomaticPreconditioningMetros) break;
      if (timeOver()) {
        pHMCProp->writeAllOmegaFieldsToDisk();
        writeCurrentStateDescriptor(0, true);
        if (Parameter_SaveStateDescriptorEveryXXXResult>0) automaticBackup();
        if (LogLevel>1) printf("Time limit reached. ==> EXITING!\n");
        desini();
	exit(0);
      }
    }
    if (LogLevel>1) printf("...Preconditioning ready. Best Preconditioning Parameters:\n");    
    fermiOps->printPreconditionerParameter();
  }

  if (Parameter_FACC_type>0) {
    pHMCProp->calcPhiMomentumMasses(Parameter_Epsilon*Parameter_Iterations[0]);    
  }
  if (Parameter_OmegaMassAdaptionMode>0) {
    pHMCProp->calcOmegaMassAdaption();
  }    


  if (ThermalizingRunCount<Parameter_ThermalizingMetros) {
    if (LogLevel>1) printf("Performing %d thermalizing Metroplis steps...\n",Parameter_ThermalizingMetros-ThermalizingRunCount);
    pHMCProp->setTheta(Parameter_Theta);  
    fermiOps->printPreconditionerParameter();
    for (ResampleCount=0; true; ResampleCount++) {
      RunCount = 0;
      acceptCount = 0;
      for (MeasureCount=0; MeasureCount<Parameter_Measurements; MeasureCount++) {
        pHMCProp->sampleALLMomenta();
        if (Parameter_FLAG_DirectOmegaSampling > 0) pHMCProp->sampleOmegaFields();
        if (markovStep()) {
          acceptCount++;
          ThermalizingRunCount++;
          writeCurrentStateDescriptor(1, false);
          if (Parameter_FACC_type>0) {
            writeFourierAccelerationData(false);	 
	  }
          if (Parameter_OmegaMassAdaptionMode>0) {
	    pHMCProp->writeOmegaForceStrengthsToDiskGLOBAL();    
	  }
        }
        RunCount++;
      }
      writePerformanceProfilingDataToDisk();
      acceptRate = acceptCount / RunCount;
      if (LogLevel>2) printf("Average (thermalizing step: %d) update quote: %1.2f\n",ThermalizingRunCount,acceptRate);

      if (ThermalizingRunCount>Parameter_ThermalizingMetros) break;
      if (timeOver()) {
        pHMCProp->writeAllOmegaFieldsToDisk();
        writeCurrentStateDescriptor(0, true);
        if (Parameter_SaveStateDescriptorEveryXXXResult>0) automaticBackup();
        if (LogLevel>1) printf("Time limit reached. ==> EXITING!\n");
        desini();
	exit(0);
      }
    }
    if (LogLevel>1) printf("...Thermalizing ready.\n");
  }
  
  if (TotallyMeasuredConfigurationsCount<Parameter_TotalData) {
    if (LogLevel>1) printf("Performing %d Measurements...\n",Parameter_TotalData-TotallyMeasuredConfigurationsCount);
    pHMCProp->setTheta(Parameter_Theta);  
    fermiOps->printPreconditionerParameter();
    for (ResampleCount=0; true; ResampleCount++) {
      RunCount = 0;
      acceptCount = 0;
      for (MeasureCount=0; MeasureCount<Parameter_Measurements; MeasureCount++) {
        pHMCProp->sampleALLMomenta();
        if (Parameter_FLAG_DirectOmegaSampling > 0) pHMCProp->sampleOmegaFields();
        if (markovStep()) {
          double exactWeight = NaN;
          if (Parameter_FLAG_ExactReweighing>0) {
            pHMCProp->calcExactMMdagInverseSQRTOmegaAction();
            exactWeight = pHMCProp->getExactReweighingFactorFromMMdagInverseSQRTOmegaAction();
            writeWeightFactor(TotallyMeasuredConfigurationsCount, 2, exactWeight);
          }
	  if ((Parameter_ReversibilityCheckFreqMeas>0) && (((TotallyMeasuredConfigurationsCount+1) % Parameter_ReversibilityCheckFreqMeas) == 0)) {
    	    writeReversibilityNorm(2);	  
	  }	
	  if ((Parameter_DetermineEWeveryXXXconfsMEASURE>0) && (TotallyMeasuredConfigurationsCount % Parameter_DetermineEWeveryXXXconfsMEASURE == 0)) {
            writeConditionNumber(TotallyMeasuredConfigurationsCount, 2);
	  }
          measureAndWrite(exactWeight);
          acceptCount++;
  	  TotallyMeasuredConfigurationsCount++;
	  writePhiConfigurationToDisk(exactWeight);
	  writeCurrentStateDescriptor(1, false);
          if (Parameter_FACC_type>0) {
            writeFourierAccelerationData(false);	 
	  }
          if (Parameter_OmegaMassAdaptionMode>0) {
	    pHMCProp->writeOmegaForceStrengthsToDiskGLOBAL();    
	  }
        } else {
          if (LogLevel>1) printf("Proposed Configuration Rejected.\n");
          double exactWeight = NaN;
          if (Parameter_FLAG_ExactReweighing>0) {
            exactWeight = pHMCProp->getExactReweighingFactorFromMMdagInverseSQRTOmegaAction();
          }
          measureAndWrite(exactWeight);
	  writePhiConfigurationToDisk(exactWeight);
	  writeCurrentStateDescriptor(1, false);
        }
        RunCount++;
      }
      writePerformanceProfilingDataToDisk();
      acceptRate = acceptCount / RunCount;
      if (LogLevel>1) printf("-->Update quote: %1.2f\n\n",acceptRate);    

      if (TotallyMeasuredConfigurationsCount>=Parameter_TotalData) break;
      if (timeOver()) {
        pHMCProp->writeAllOmegaFieldsToDisk();
        writeCurrentStateDescriptor(0, true);
        if (Parameter_SaveStateDescriptorEveryXXXResult>0) automaticBackup();
        if (LogLevel>1) printf("Time limit reached. ==> EXITING!\n");
        desini();
	exit(0);
      }
    }
  }
  
  pHMCProp->writeAllOmegaFieldsToDisk();
  writeCurrentStateDescriptor(0, true);
  if (Parameter_SaveStateDescriptorEveryXXXResult>0) automaticBackup();
  
  if (LogLevel>1) printf("pHMC ready. Collected %d data samples!!!\n", TotallyMeasuredConfigurationsCount);

  desini();
  desiniTools();
  fftw_cleanup_threads();
  if (LogLevel>1) printf("Synchronizing caches...\n"); 
  sync();
  if (LogLevel>1) printf("pHMC terminated correctly.\n");
}
