#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>

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
#include "HMCPropagator.h"

#define AutomaticAdaption_HighPro 0.98
#define AutomaticAdaption_LowPro 0.80
#define AutomaticAdaption_LowPro2 0.10

#define AutomaticIntegratorAdaption_MaxLeapFrogInversions 40
#define AutomaticPropTolAdaption_MaxPropTol 0.1


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
double Parameter_Epsilon = 0;
double Parameter_PropTOL = 0;
double Parameter_FinalTOL = 0;
int Parameter_Iterations = 0;
int Parameter_Measurements = 0;
int Parameter_AutomaticPreconditioningMetros = 0;
int Parameter_ThermalizingMetros = 0;
int Parameter_TotalData = 0;
char* Parameter_filenamePrefix;
int Parameter_FLAG_CalcNeubergerDetWithXi = 0;
int Parameter_FLAG_CalcWilsonDetWithoutXi = 0;
int Parameter_FLAG_CalcNeubergerDetWithoutXi = 0;
int Parameter_FLAG_CalcWilsonDetWithXi = 0;
int Parameter_FLAG_ReadWriteStateDescriptor = 0;
int Parameter_FLAG_WriteConditionNumber = 0;
int Parameter_FLAG_xFFT = 0;
double Parameter_MaxRunTime = 0;
int Parameter_IntegratorType = 0;
int Automatic_Adaption_LastAction = 1;
int Automatic_Adaption_LastObject = 1;
FermionMatrixOperations* fermiOps = NULL;
HMCPropagator* HMCProp = NULL;
int ThermalizingRunCount = 0;
int AutomaticPreconRunCount = 0;
int TotallyMeasuredConfigurationsCount = 0;
char* outputFileNameExtension1 = NULL;
char* outputFileNameExtension2 = NULL;
clock_t startTime;


void startTimer() {
  startTime = clock();
}


double timePassed() {
  double time = (clock()-startTime)/CLOCKS_PER_SEC;
  return time;
}

bool timeOver() {
  if (timePassed()>Parameter_MaxRunTime*3600) return true;
  return false;
}


void loadParameters(int MultiProcessNr, bool ll) {
  FILE* file;
  if (LogLevel>0) printf("Trying to read parameter-file 'SimulationParametersHybrid.txt'...\n");
  file = fopen("SimulationParametersHybrid.txt","r");
  char *Comment = new char[500];
  bool error = false;
  Parameter_filenamePrefix = new char[500];

  if (fscanf(file,"%s",Parameter_filenamePrefix)==1) {
    if (LogLevel>0) printf(" -> Filename-Prefix: %s\n",Parameter_filenamePrefix);
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
  
  if (fscanf(file,"%lf",&Parameter_Lambda)==1) {
    if (LogLevel>0) printf(" -> Lambda: %f\n",Parameter_Lambda);
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
  
  MultiProcessNr = MultiProcessNr % (yTics*kappaTics);
  if (kappaTics <= 1) {
    Parameter_Kappa = kappaMin;
  } else {
    Parameter_Kappa = kappaMin+(kappaMax-kappaMin)*(MultiProcessNr % kappaTics)/(kappaTics-1);
  }
  if (yTics <= 1) {
    Parameter_Y = yMin;
  } else {
    Parameter_Y = yMin+(yMax-yMin)*(MultiProcessNr / kappaTics)/(yTics-1);
  }
  if (LogLevel>0) printf(" -> Y (selected): %f\n",Parameter_Y);
  if (LogLevel>0) printf(" -> Kappa (selected): %f\n",Parameter_Kappa);

  if (fscanf(file,"%d",&Parameter_FLAG_CalcNeubergerDetWithXi)==1) {
    if (LogLevel>0) printf(" -> Flag - Calculate Neuberger Determinant with Xi-fields: %d\n",Parameter_FLAG_CalcNeubergerDetWithXi);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_CalcWilsonDetWithoutXi)==1) {
    if (LogLevel>0) printf(" -> Flag - Calculate Wilson Determinant without Xi-fields: %d\n",Parameter_FLAG_CalcWilsonDetWithoutXi);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_CalcNeubergerDetWithoutXi)==1) {
    if (LogLevel>0) printf(" -> Flag - Calculate Neuberger Determinant without Xi-fields: %d\n",Parameter_FLAG_CalcNeubergerDetWithoutXi);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_CalcWilsonDetWithXi)==1) {
    if (LogLevel>0) printf(" -> Flag - Calculate Wilson Determinant with Xi-fields: %d\n",Parameter_FLAG_CalcWilsonDetWithXi);
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

  if (fscanf(file,"%lf",&Parameter_PropTOL)==1) {
    if (LogLevel>0) printf(" -> Absolute Precision of Matrix Inversion during propagation: %1.15f\n",Parameter_PropTOL);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%lf",&Parameter_FinalTOL)==1) {
    if (LogLevel>0) printf(" -> Absolute Precision of Matrix Inversion during action evaluation: %1.15f\n",Parameter_FinalTOL);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_IntegratorType)==1) {
    if (LogLevel>0) printf(" -> Type of integrator being used: %d\n",Parameter_IntegratorType);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_Iterations)==1) {
    if (LogLevel>0) printf(" -> Iterations between Metropolis steps: %d\n",Parameter_Iterations);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_Measurements)==1) {
    if (LogLevel>0) printf(" -> Measurement attempts before resampling of p and omega: %d\n",Parameter_Measurements);
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

  if (fscanf(file,"%d",&Parameter_FLAG_ReadWriteStateDescriptor)==1) {
    if (LogLevel>0) printf(" -> Read/Write State-Descriptor: %d\n",Parameter_FLAG_ReadWriteStateDescriptor);
  } else error = true;
  fgets(Comment, 500, file);

  if (fscanf(file,"%d",&Parameter_FLAG_WriteConditionNumber)==1) {
    if (LogLevel>0) printf(" -> Write Condition-Number: %d\n",Parameter_FLAG_WriteConditionNumber);
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

  //Check parameters
  if ((Parameter_Nf % 2 == 1) && (Parameter_Y>0)) {
    printf("This Hybrid Monte Carlo can only run with even Nf if fermions included!!!\n");
    exit(0);
  }

  fclose(file);
  if (error) {
    printf("... ERROR!\n");
    exit(0);
  }  
  delete [] Comment;
  if ((LogLevel>0) && (ll)) printf("...sucessfully.\n");
}


void desini() {
  delete[] outputFileNameExtension1;
  delete[] outputFileNameExtension2;
  HMCProp->killSlaves();
  delete HMCProp;
  delete fermiOps;
  #ifdef UseMPI
  MPI_Abort(MPI_COMM_WORLD, 0);
  #endif
}


void writeCurrentStateDescriptor(int running) {
  if (Parameter_FLAG_ReadWriteStateDescriptor == 0) {
    return;
  }
  
  if (LogLevel>2) printf("Writing StateDescriptor to disk...\n");
  double measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm;
  HMCProp->measure(measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm);
  
  bool usePrec;
  double precM, precS;
  fermiOps->getPreconditionerParameter(usePrec, precM, precS);
  
  char* fileName = new char[600]; 
  snprintf(fileName,600,"%s/data/results/hybrid/states/StateDescriptor%s",DataBaseDirectory,outputFileNameExtension1);
  FILE* file;
  file = fopen(fileName,"w");

  fprintf(file,"%d                                                      : Running\n", running);
  fprintf(file,"%1.15f                                      : Average Phi\n", measurePhiNorm);
  fprintf(file,"%1.15f                                      : Average Staggered Phi\n", measureStaggeredPhiNorm);
  fprintf(file,"%d                                                      : Use Preconditioner\n", usePrec);
  fprintf(file,"%1.15f                                      : Preconditioner parameter M\n", precM);
  fprintf(file,"%1.15f                                      : Preconditioner parameter S\n", precS);
  fprintf(file,"%d                                                      : Performed Preconditioning successful Metros\n", AutomaticPreconRunCount);
  fprintf(file,"%d                                                      : Performed Thermalizing successful Metros\n", ThermalizingRunCount);
  fprintf(file,"%d                                                      : Performed Measuring successful Metros Steps\n", TotallyMeasuredConfigurationsCount);
  fprintf(file,"%s                                                 : Filename-Prefix\n", Parameter_filenamePrefix);
  fprintf(file,"%s          : Filename-Extension no Job-Extension\n", outputFileNameExtension1);
  fprintf(file,"%s    : Filename-Extension Full\n", outputFileNameExtension2);
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
  fprintf(file,"%d                                                      : Flag - Calculate Neuberger Determinant with Xi-fields\n", Parameter_FLAG_CalcNeubergerDetWithXi);
  fprintf(file,"%d                                                      : Flag - Calculate Wilson Determinant without Xi-fields \n", Parameter_FLAG_CalcWilsonDetWithoutXi);
  fprintf(file,"%d                                                      : Flag - Calculate Neuberger Determinant without Xi-fields  \n", Parameter_FLAG_CalcNeubergerDetWithoutXi);
  fprintf(file,"%d                                                      : Flag - Calculate Wilson Determinant with Xi-fields  \n", Parameter_FLAG_CalcWilsonDetWithXi);
  fprintf(file,"%1.15f                                      : Epsilon for molecular dynamics \n",Parameter_Epsilon );
  fprintf(file,"%1.15f                                      : Absolute Precision of Matrix Inversion during propagation \n", Parameter_PropTOL);
  fprintf(file,"%1.15f                                      : Absolute Precision of Matrix Inversion during action evaluation \n", Parameter_FinalTOL);
  fprintf(file,"%d                                                      : Type of integrator being used \n", Parameter_IntegratorType);
  fprintf(file,"%d                                                     : Iterations between Metropolis steps \n", Parameter_Iterations);
  fprintf(file,"%d                                                     : Measurement attempts before resampling of p and omega \n", Parameter_Measurements);
  fprintf(file,"%d                                                    : Number of Metropolis steps for automatic Preconditioning \n", Parameter_AutomaticPreconditioningMetros);
  fprintf(file,"%d                                                    : Number of Metropolis steps for thermalizing \n", Parameter_ThermalizingMetros);
  fprintf(file,"%d                                                  : Number of total data to be collected \n", Parameter_TotalData);
  fprintf(file,"%d                                             : Current RandSeed \n", AdvancedSeed);
  
  
  fprintf(file,"\n   *** Phi - Field ***\n");
  int I;
  int VL = fermiOps->getVectorLength();
  for (I=0; I<VL/8; I++) {
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f\n",HMCProp->phiField[I][0],HMCProp->phiField[I][1],HMCProp->phiField[I][2],HMCProp->phiField[I][3]);
  }  
  fclose(file);  
  delete [] fileName;
}


void readCurrentStateDescriptor() {
  //Default Settings if no Descriptor-File is found
  AutomaticPreconRunCount = 0;
  ThermalizingRunCount = 0;
  TotallyMeasuredConfigurationsCount = 0;
  //Default Setting END
  
  if (Parameter_FLAG_ReadWriteStateDescriptor == 0) {
    return;
  }

  char *Comment = new char[500];
  char* fileName = new char[600];
  snprintf(fileName,600,"%s/data/results/hybrid/states/StateDescriptor%s",DataBaseDirectory,outputFileNameExtension1);
  if (LogLevel>1) printf("Trying to read DescriptorStateFile: %s\n",fileName);  
  FILE* file;
  file = fopen(fileName,"r");
  delete[] fileName;


  if (file==NULL) {
    if (LogLevel>0) printf("Could not open StateDescriptor File.\n");
    return;
  }


  int running = -1;
  if ((fscanf(file,"%d",&running)!=1) || (running<0) || (running>1)) {
    if (LogLevel>0) printf("Could not open StateDescriptor File.\n");
    return;
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
         READParameter_Lambda,READParameter_Kappa,READParameter_Y,READParameter_Epsilon,READParameter_PropTOL,
         READParameter_FinalTOL;
  
  int READusePrec,READAutomaticPreconRunCount,READThermalizingRunCount,READTotallyMeasuredConfigurationsCount,
      READParameter_L0,READParameter_L1,READParameter_L2,READParameter_L3,
      READParameter_Nf,READParameter_FLAG_CalcNeubergerDetWithXi,READParameter_FLAG_CalcWilsonDetWithoutXi,
      READParameter_FLAG_CalcNeubergerDetWithoutXi,READParameter_FLAG_CalcWilsonDetWithXi,READParameter_IntegratorType,
      READParameter_Iterations,READParameter_Measurements,READParameter_AutomaticPreconditioningMetros,READParameter_ThermalizingMetros,
      READParameter_TotalData,READAdvancedSeed;
  
  char* READParameter_filenamePrefix = new char[300];
  char* READoutputFileNameExtension1 = new char[300];  
  char* READoutputFileNameExtension2 = new char[300];  
  
  if (fscanf(file,"%lf", &READmeasurePhiNorm)!=1) readError = 1;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READmeasureStaggeredPhiNorm)!=1) readError = 2;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READusePrec)!=1) readError = 3;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READprecM)!=1) readError = 4;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READprecS)!=1) readError = 5;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READAutomaticPreconRunCount)!=1) readError = 6;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READThermalizingRunCount)!=1) readError = 7;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READTotallyMeasuredConfigurationsCount)!=1) readError = 8;
  fgets(Comment, 500, file);
  if (fscanf(file,"%s", READParameter_filenamePrefix)!=1) readError = 9;
  fgets(Comment, 500, file);
  if (fscanf(file,"%s", READoutputFileNameExtension1)!=1) readError = 10;
  fgets(Comment, 500, file);
  if (fscanf(file,"%s", READoutputFileNameExtension2)!=1) readError = 11;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_L0)!=1) readError = 120;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_L1)!=1) readError = 121;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_L2)!=1) readError = 122;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_L3)!=1) readError = 123;
  fgets(Comment, 500, file);  
  if (fscanf(file,"%d", &READParameter_Nf)!=1) readError = 13;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_RHO)!=1) readError = 14;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_R)!=1) readError = 15;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_Lambda)!=1) readError = 16;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_Kappa)!=1) readError = 17;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_Y)!=1) readError = 18;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_CalcNeubergerDetWithXi)!=1) readError = 19;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_CalcWilsonDetWithoutXi)!=1) readError = 20;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_CalcNeubergerDetWithoutXi)!=1) readError = 21;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_FLAG_CalcWilsonDetWithXi)!=1) readError = 22;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf",&READParameter_Epsilon)!=1) readError = 23;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_PropTOL)!=1) readError = 24;
  fgets(Comment, 500, file);
  if (fscanf(file,"%lf", &READParameter_FinalTOL)!=1) readError = 25;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_IntegratorType)!=1) readError = 26;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_Iterations)!=1) readError = 27;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_Measurements)!=1) readError = 28;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_AutomaticPreconditioningMetros)!=1) readError = 29;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_ThermalizingMetros)!=1) readError = 30;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READParameter_TotalData)!=1) readError = 31;
  fgets(Comment, 500, file);
  if (fscanf(file,"%d", &READAdvancedSeed)!=1) readError = 32;
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
  if (Parameter_L0!=READParameter_L0) consistencyError=100;
  if (Parameter_L1!=READParameter_L1) consistencyError=101;
  if (Parameter_L2!=READParameter_L2) consistencyError=102;
  if (Parameter_L3!=READParameter_L3) consistencyError=103;
  if (Parameter_Nf!=READParameter_Nf) consistencyError=2;
  if (fabs(Parameter_RHO-READParameter_RHO)>ceps) consistencyError=3;
  if (fabs(Parameter_R-READParameter_R)>ceps) consistencyError=4;
  if (fabs(Parameter_Lambda-READParameter_Lambda)>ceps) consistencyError=5;
  if (fabs(Parameter_Kappa-READParameter_Kappa)>ceps) consistencyError=6;
  if (fabs(Parameter_Y-READParameter_Y)>ceps) consistencyError=7;
  if (Parameter_FLAG_CalcNeubergerDetWithXi!=READParameter_FLAG_CalcNeubergerDetWithXi) consistencyError=8;
  if (Parameter_FLAG_CalcWilsonDetWithoutXi!=READParameter_FLAG_CalcWilsonDetWithoutXi) consistencyError=9;
  if (Parameter_FLAG_CalcNeubergerDetWithoutXi!=READParameter_FLAG_CalcNeubergerDetWithoutXi) consistencyError=10;
  if (Parameter_FLAG_CalcWilsonDetWithXi!=READParameter_FLAG_CalcWilsonDetWithXi) consistencyError=11;
  if (fabs(Parameter_FinalTOL-READParameter_FinalTOL)>ceps) consistencyError=12;
  if (Parameter_Measurements!=READParameter_Measurements) consistencyError=13;
  if (Parameter_AutomaticPreconditioningMetros!=READParameter_AutomaticPreconditioningMetros) consistencyError=14;
  if (Parameter_ThermalizingMetros!=READParameter_ThermalizingMetros) consistencyError=15;
  
  if (consistencyError>0) {
    fclose(file);
    if (LogLevel>0) printf("Consistency-Error (check nr. %d) with data read from File-Descriptor. ==> Exiting\n",consistencyError);
    desini();
    exit(0);  
  }

  //Take over Data
  if (LogLevel>0) printf("Taking over values from StateDescriptor File...\n");
  fermiOps->setPreconditioner(READusePrec, READprecM, READprecS);
    
  AutomaticPreconRunCount = READAutomaticPreconRunCount;
  if (LogLevel>0) printf("  ***> AutomaticPreconRunCount=%d\n",AutomaticPreconRunCount);
  ThermalizingRunCount = READThermalizingRunCount;
  if (LogLevel>0) printf("  ***> ThermalizingRunCount=%d\n",ThermalizingRunCount);
  TotallyMeasuredConfigurationsCount = READTotallyMeasuredConfigurationsCount;
  if (LogLevel>0) printf("  ***> TotallyMeasuredConfigurationsCount=%d\n",TotallyMeasuredConfigurationsCount);
  
  snprintf(Parameter_filenamePrefix,600,"%s",READParameter_filenamePrefix);
  if (LogLevel>0) printf("  ***> Parameter_filenamePrefix=%s\n",Parameter_filenamePrefix); 
  snprintf(outputFileNameExtension1,600,"%s",READoutputFileNameExtension1);
  if (LogLevel>0) printf("  ***> outputFileNameExtension1=%s\n",outputFileNameExtension1);   
  snprintf(outputFileNameExtension2,600,"%s",READoutputFileNameExtension2);
  if (LogLevel>0) printf("  ***> outputFileNameExtension2=%s\n",outputFileNameExtension2);   
  
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
  Parameter_FLAG_CalcNeubergerDetWithXi = READParameter_FLAG_CalcNeubergerDetWithXi;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_CalcNeubergerDetWithXi=%d\n",Parameter_FLAG_CalcNeubergerDetWithXi);
  Parameter_FLAG_CalcWilsonDetWithoutXi = READParameter_FLAG_CalcWilsonDetWithoutXi;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_CalcWilsonDetWithoutXi=%d\n",Parameter_FLAG_CalcWilsonDetWithoutXi);
  Parameter_FLAG_CalcNeubergerDetWithoutXi = READParameter_FLAG_CalcNeubergerDetWithoutXi;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_CalcNeubergerDetWithoutXi=%d\n",Parameter_FLAG_CalcNeubergerDetWithoutXi);
  Parameter_FLAG_CalcWilsonDetWithXi = READParameter_FLAG_CalcWilsonDetWithXi;
  if (LogLevel>0) printf("  ***> Parameter_FLAG_CalcWilsonDetWithXi=%d\n",Parameter_FLAG_CalcWilsonDetWithXi);
  Parameter_Epsilon = READParameter_Epsilon;
  if (LogLevel>0) printf("  ***> Parameter_Epsilon=%f\n",Parameter_Epsilon);
  Parameter_PropTOL = READParameter_PropTOL;
  if (LogLevel>0) printf("  ***> Parameter_PropTOL=%f\n",Parameter_PropTOL);
  Parameter_FinalTOL = READParameter_FinalTOL;
  if (LogLevel>0) printf("  ***> Parameter_FinalTOL=%1.15f\n",Parameter_FinalTOL);
  Parameter_IntegratorType = READParameter_IntegratorType;
  if (LogLevel>0) printf("  ***> Parameter_IntegratorType=%d\n",Parameter_IntegratorType);
  Parameter_Iterations = READParameter_Iterations;
  if (LogLevel>0) printf("  ***> Parameter_Iterations=%d\n",Parameter_Iterations);
  Parameter_Measurements = READParameter_Measurements;
  if (LogLevel>0) printf("  ***> Parameter_Measurements=%d\n",Parameter_Measurements);
  Parameter_AutomaticPreconditioningMetros = READParameter_AutomaticPreconditioningMetros;
  if (LogLevel>0) printf("  ***> Parameter_AutomaticPreconditioningMetros=%d\n",Parameter_AutomaticPreconditioningMetros);
  Parameter_ThermalizingMetros = READParameter_ThermalizingMetros;
  if (LogLevel>0) printf("  ***> Parameter_ThermalizingMetros=%d\n",Parameter_ThermalizingMetros);
  Parameter_TotalData = READParameter_TotalData;
  if (LogLevel>0) printf("  ***> Parameter_TotalData=%d\n",Parameter_TotalData);

  AdvancedSeed = -READAdvancedSeed;
  if (LogLevel>0) printf("  ***> AdvancedSeed=%d\n",AdvancedSeed);
  AdvancedZufall(AdvancedSeed);
  if (LogLevel>0) printf("  ***> AdvancedSeed=%d after ini, and first random=%f\n",AdvancedSeed, AdvancedZufall(AdvancedSeed));
  
  
  //Reading Phi-Field 
  fgets(Comment, 500, file);
  fgets(Comment, 500, file);
  
  int I;
  int VL = fermiOps->getVectorLength();
  readError = -1;
  for (I=0; I<VL/8; I++) {
    if (fscanf(file,"%lf %lf %lf %lf",&(HMCProp->phiField[I][0]),&(HMCProp->phiField[I][1]),&(HMCProp->phiField[I][2]),&(HMCProp->phiField[I][3]))!=4) {
      readError = I;
      break;
    }
  }
  fclose(file);
  delete[] Comment;
  
  if (readError>=0) {
    if (LogLevel>0) printf("Error (index %d) while reading Phi-Field from File-Descriptor. ==> Exiting\n",readError);
    desini();
    exit(0);  
  }

  //Consistency-Check on Phi-Field
  double measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm;  
  HMCProp->measure(measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm);
  if (LogLevel>0) printf("Configuration read with mag=%1.15f and stag. mag=%1.15f\n",measurePhiNorm,measureStaggeredPhiNorm);
  if (LogLevel>0) printf("Old result              mag=%1.15f and stag. mag=%1.15f\n",READmeasurePhiNorm,READmeasureStaggeredPhiNorm); 
  
  if ((fabs(READmeasurePhiNorm-measurePhiNorm)>ceps) || (fabs(READmeasureStaggeredPhiNorm-measureStaggeredPhiNorm)>ceps)) {
    if (LogLevel>0) printf("Consistency-Error for read Phi-Field from File-Descriptor. ==> Exiting\n");
    desini();
    exit(0);  
  }
}


void buildOutputFileNameExtension(int valForJobID) {
  outputFileNameExtension1 = new char[300];
  outputFileNameExtension2 = new char[300];
  snprintf(outputFileNameExtension1,300,"L%dx%dx%dx%dNf%dKap%1.3fLam%1.3fY%1.3fRho%1.3fR%1.3f.dat",
  Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3,Parameter_Nf,Parameter_Kappa,Parameter_Lambda,Parameter_Y,Parameter_RHO,Parameter_R);	 
  snprintf(outputFileNameExtension2,300,"L%dx%dx%dx%dNf%dKap%1.3fLam%1.3fY%1.3fRho%1.3fR%1.3fJob%d.dat",
  Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3,Parameter_Nf,Parameter_Kappa,Parameter_Lambda,Parameter_Y,Parameter_RHO,Parameter_R,valForJobID);	 
}  


void measureAndWrite() {
  double measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm;
  HMCProp->measure(measurePhiNorm, measureStaggeredPhiNorm, avgNorm, sigmaNorm);  
    
  char* fileName = new char[600]; 
  snprintf(fileName,600,"%s/data/results/hybrid/%s%s",DataBaseDirectory,Parameter_filenamePrefix,outputFileNameExtension2);	 
  FILE* file;

  file = fopen(fileName,"a");
  fprintf(file,"%f %f 0 nan nan nan nan nan nan nan\n", measurePhiNorm,measureStaggeredPhiNorm);
  fclose(file);  

  delete[] fileName;
}


void writeConditionNumber() {
  double eigMin;
  double eigMax;
  double invCond;

  double* phiField = (double*) HMCProp->phiField;

  fermiOps->exactFermionMatrixConditionNumber(phiField, eigMin, eigMax, invCond, true, 1, false);

  char* fileName = new char[600]; 
  snprintf(fileName,600,"%s/data/results/hybrid/%s%s",DataBaseDirectory,"conditionNumber",outputFileNameExtension2);	 
  FILE* file;

  file = fopen(fileName,"a");
  fprintf(file,"%f %f %f\n", eigMin, eigMax, invCond);
  fclose(file);  

  delete[] fileName;
}



void automaticEpsilonAdaption(double acceptRate) {
  int intOrder = 2;
  if (Parameter_IntegratorType==2) intOrder = 4;
  
  if (acceptRate<AutomaticAdaption_LowPro) {
    int iterNeu = (int)((Parameter_Iterations+1)/0.90) + 1;
    Parameter_Epsilon = (Parameter_Iterations*Parameter_Epsilon) / iterNeu;    
    Parameter_Iterations = iterNeu;
    if (LogLevel>2) printf("==> Decreasing epsilon to %f and increasing iterations to %d (Integrator order %d).\n",Parameter_Epsilon, Parameter_Iterations,intOrder);
  }
  if (acceptRate<AutomaticAdaption_LowPro2) {
    int iterNeu = (int)(2.00 * Parameter_Iterations) + 1;
    Parameter_Epsilon = (Parameter_Iterations*Parameter_Epsilon) / iterNeu;    
    Parameter_Iterations = iterNeu;
    if (LogLevel>2) printf("==> Decreasing epsilon to %f and increasing iterations to %d (Integrator order %d).\n",Parameter_Epsilon, Parameter_Iterations,intOrder);
  }
  if (acceptRate>AutomaticAdaption_HighPro) {
    int iterNeu = (int)(0.90 * Parameter_Iterations) - 1;
    if (iterNeu<=1) iterNeu = 2;
    Parameter_Epsilon = (Parameter_Iterations*Parameter_Epsilon) / iterNeu;    
    Parameter_Iterations = iterNeu;
    if (LogLevel>2) printf("==> Increasing epsilon to %f and decreasing iterations to %d (Integrator order %d).\n",Parameter_Epsilon, Parameter_Iterations,intOrder);
  }
  
  if ((Parameter_IntegratorType==2) && (Parameter_Iterations*5<AutomaticIntegratorAdaption_MaxLeapFrogInversions)) {
    if (LogLevel>2) printf("==> Change integration scheme from O(4) to O(2), since iterations=%d.\n",Parameter_Iterations);
    int iterNeu = (int)(5 * Parameter_Iterations);
    if (iterNeu<=1) iterNeu = 2;
    Parameter_Epsilon = (Parameter_Iterations*Parameter_Epsilon) / iterNeu;    
    Parameter_Iterations = iterNeu;
    Parameter_IntegratorType = 0;
    if (LogLevel>2) printf("==> New epsilon = %f and iterations = %d.\n",Parameter_Epsilon, Parameter_Iterations);
  } else if (Parameter_IntegratorType==0) {
    int wouldBeIter = Parameter_Iterations / 5;
    if (wouldBeIter <= 1) wouldBeIter = 2;
    if (AutomaticIntegratorAdaption_MaxLeapFrogInversions <= wouldBeIter*5) {
      if (LogLevel>2) printf("==> Change integration scheme from O(2) to O(4), since iterations=%d.\n",Parameter_Iterations);
      Parameter_Epsilon = (Parameter_Iterations*Parameter_Epsilon) / wouldBeIter;    
      Parameter_Iterations = wouldBeIter;
      Parameter_IntegratorType = 2;
      if (LogLevel>2) printf("==> New epsilon = %f and iterations = %d.\n",Parameter_Epsilon, Parameter_Iterations);
    }
  }
}


void automaticPropTolAdaption(double acceptRate) {
  if (acceptRate<AutomaticAdaption_LowPro) {
    double fac = 1.15 + (AdvancedZufall(AdvancedSeed)/10.0);
    fac *= 1.5*(fabs(log(Parameter_PropTOL)/log(AutomaticPropTolAdaption_MaxPropTol)));
    Parameter_PropTOL /= fac;
    if (Parameter_PropTOL <= Parameter_FinalTOL) Parameter_PropTOL = Parameter_FinalTOL;    
    if (LogLevel>2) printf("==> Decreasing PropTol to %1.15f (FinalTol = %1.15f).\n",Parameter_PropTOL, Parameter_FinalTOL);
  }
  if (acceptRate<AutomaticAdaption_LowPro2) {
    double fac = 1.15 + (AdvancedZufall(AdvancedSeed)/10.0);
    fac *= 10;
    Parameter_PropTOL /= fac;
    if (Parameter_PropTOL <= Parameter_FinalTOL) Parameter_PropTOL = Parameter_FinalTOL;    
    if (LogLevel>2) printf("==> Decreasing PropTol to %1.15f (FinalTol = %1.15f).\n",Parameter_PropTOL, Parameter_FinalTOL);
  }
  if (acceptRate>AutomaticAdaption_HighPro) {
    double fac = 1.15 + (AdvancedZufall(AdvancedSeed)/10.0);
    fac *= (fabs(log(Parameter_PropTOL)/log(AutomaticPropTolAdaption_MaxPropTol)));
    Parameter_PropTOL *= fac;
    if (Parameter_PropTOL >= AutomaticPropTolAdaption_MaxPropTol) Parameter_PropTOL = AutomaticPropTolAdaption_MaxPropTol;    
    if (LogLevel>2) printf("==> Increasing PropTol to %1.15f (FinalTol = %1.15f).\n",Parameter_PropTOL, Parameter_FinalTOL);
  }
}


void automaticAdaption(double acceptRate) {
  if (acceptRate<AutomaticAdaption_LowPro) {
    if (Automatic_Adaption_LastAction == 1) {
      if (Automatic_Adaption_LastObject == 1) {
        automaticEpsilonAdaption(acceptRate);
	Automatic_Adaption_LastObject = 1;
	Automatic_Adaption_LastAction = -1;
	return;
      } else {
        automaticPropTolAdaption(acceptRate);
	Automatic_Adaption_LastObject = 2;
	Automatic_Adaption_LastAction = -1;
	return;
      }
    } else {
      if (Automatic_Adaption_LastObject == 1) {
        automaticPropTolAdaption(acceptRate);
	Automatic_Adaption_LastObject = 2;
	Automatic_Adaption_LastAction = -1;
	return;
      } else {
        automaticEpsilonAdaption(acceptRate);
	Automatic_Adaption_LastObject = 1;
	Automatic_Adaption_LastAction = -1;
	return;
      }
    }
  }
  if (acceptRate>AutomaticAdaption_HighPro) {
    if (Automatic_Adaption_LastAction == 1) {
      if (Automatic_Adaption_LastObject == 1) {
        automaticPropTolAdaption(acceptRate);
	Automatic_Adaption_LastObject = 2;
	Automatic_Adaption_LastAction = 1;
	return;
      } else {
        automaticEpsilonAdaption(acceptRate);
	Automatic_Adaption_LastObject = 1;
	Automatic_Adaption_LastAction = 1;
	return;
      }
    } else {
      if (Automatic_Adaption_LastObject == 1) {
        automaticPropTolAdaption(acceptRate);
	Automatic_Adaption_LastObject = 2;
	Automatic_Adaption_LastAction = 1;
	return;
      } else {
        automaticEpsilonAdaption(acceptRate);
	Automatic_Adaption_LastObject = 1;
	Automatic_Adaption_LastAction = 1;
	return;
      }
    }
  }
}


bool propagateFields(int iter, double eps, double propTOL, double finalTOL) {
  bool b;
  if (Parameter_IntegratorType == 0) {
    b = HMCProp->LeapFrogMarkovStep(iter, eps, propTOL, finalTOL);
  } else if (Parameter_IntegratorType == 1) {
    b = HMCProp->OmelyanO2MarkovStep(iter, eps, propTOL, finalTOL);
  } else if (Parameter_IntegratorType == 2) {
    b = HMCProp->OmelyanO4MarkovStep(iter, eps, propTOL, finalTOL);
  } else {
    printf("ERROR: Invalid integrator!!!\n");
    exit(0);
  }

  return b;
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


int main(int argc,char **argv) {
  double acceptRate, RunCount, acceptCount;
  int ResampleCount, MeasureCount;
  int ranIniValue = 1;

  iniMPI(argc, argv);
  loadParameters(0,false);
  startTimer();
  if (LogLevel>1) printf("Timer started. Passed time: %f\n",timePassed());
  
  if (LogLevel>1) printf("Number of arguments = %d\n",argc);
  int I;
  for (I=0; I<argc; I++) {
    if (LogLevel>1) printf("Argument %d: %s\n",I+1,argv[I]);  
  }
  
  int JobNr = -1;
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
  loadParameters(MultiProcessNr,true);
  buildOutputFileNameExtension(JobNr);
  
  fermiOps = new FermionMatrixOperations(Parameter_L0,Parameter_L1,Parameter_L2,Parameter_L3, Parameter_RHO, Parameter_R, Parameter_Y);
  HMCProp = new HMCPropagator(fermiOps, Parameter_Lambda, Parameter_Kappa, Parameter_Nf, 1.0);

  if (Parameter_FLAG_xFFT==1) {
    fermiOps->setxFFTusage(true);
  } else {
    fermiOps->setxFFTusage(false);
  }

  if (ownNodeID>0) {
    HMCProp->SlaveController();
    if (LogLevel>1) printf("Hybrid Monte Carlo Slave finished!!!  ==>  EXITING !!!\n");
    exit(0);
  }


  readCurrentStateDescriptor();

  fermiOps->testFourierTrafo(true);
  
  if ((Parameter_AutomaticPreconditioningMetros>0) && (AutomaticPreconRunCount<Parameter_AutomaticPreconditioningMetros)) {
    double swapFrac = 0.25;
    if (LogLevel>1) printf("Performing %d Automatic Preconditioning Metroplis steps...\n",Parameter_AutomaticPreconditioningMetros-AutomaticPreconRunCount);
    if (AutomaticPreconRunCount==0) {
      fermiOps->setPreconditioner(true, 1, 0);
    }
    for (ResampleCount=0; true; ResampleCount++) {
      RunCount = 0;
      acceptCount = 0;
      
      for (MeasureCount=0; MeasureCount<Parameter_Measurements; MeasureCount++) {
        HMCProp->samplePhiMomentumField();
        HMCProp->sampleOmegaFields();
      
        if (propagateFields(Parameter_Iterations, Parameter_Epsilon, Parameter_PropTOL, Parameter_FinalTOL)) {
	  if (AutomaticPreconRunCount<swapFrac*Parameter_AutomaticPreconditioningMetros) {
            HMCProp->improvePreconditioningParameters(10, 5, Parameter_PropTOL);
	  }
          acceptCount++;
          AutomaticPreconRunCount++;
          writeCurrentStateDescriptor(1);
        }
        RunCount++;
      }
      if (AutomaticPreconRunCount>=swapFrac*Parameter_AutomaticPreconditioningMetros) {
        HMCProp->improvePreconditioningParameters(10*Parameter_Measurements, 10*Parameter_Measurements,	Parameter_PropTOL);
      }
      acceptRate = acceptCount / RunCount;
      if (LogLevel>2) printf("Average (preconditioning) update quote: %1.2f\n",acceptRate);
      automaticAdaption(acceptRate);
      if (AutomaticPreconRunCount>Parameter_AutomaticPreconditioningMetros) break;
      if (timeOver()) {
        writeCurrentStateDescriptor(0);
        if (LogLevel>1) printf("Time limit reached. ==> EXITING!\n");
        desini();
	exit(0);
      }
    }
    double PrecM, PrecS;
    bool PrecUse;
    fermiOps->getPreconditionerParameter(PrecUse, PrecM, PrecS);
    if (LogLevel>1) printf("...Preconditioning ready. Best Preconditioning Parameters PrecM = %1.3f, PrecS = %1.3f\n", PrecM, PrecS);    
  }

  if (ThermalizingRunCount<Parameter_ThermalizingMetros) {
    if (LogLevel>1) printf("Performing %d thermalizing Metroplis steps...\n",Parameter_ThermalizingMetros-ThermalizingRunCount);
    fermiOps->printPreconditionerParameter();
    for (ResampleCount=0; true; ResampleCount++) {
      RunCount = 0;
      acceptCount = 0;
      for (MeasureCount=0; MeasureCount<Parameter_Measurements; MeasureCount++) {
        HMCProp->samplePhiMomentumField();
        HMCProp->sampleOmegaFields();
      
        if (propagateFields(Parameter_Iterations, Parameter_Epsilon, Parameter_PropTOL, Parameter_FinalTOL)) {
          acceptCount++;
          ThermalizingRunCount++;
          writeCurrentStateDescriptor(1);
        }
        RunCount++;
      }
      acceptRate = acceptCount / RunCount;
      if (LogLevel>2) printf("Average (thermalizing) update quote: %1.2f\n",acceptRate);
      automaticAdaption(acceptRate);
      if (ThermalizingRunCount>Parameter_ThermalizingMetros) break;
      if (timeOver()) {
        writeCurrentStateDescriptor(0);
        if (LogLevel>1) printf("Time limit reached. ==> EXITING!\n");
        desini();
	exit(0);
      }
    }
    if (LogLevel>1) printf("...Thermalizing ready.\n");
  }
  
  if (TotallyMeasuredConfigurationsCount<Parameter_TotalData) {
    if (LogLevel>1) printf("Performing %d Measurements...\n",Parameter_TotalData-TotallyMeasuredConfigurationsCount);
    fermiOps->printPreconditionerParameter();
    for (ResampleCount=0; true; ResampleCount++) {
      RunCount = 0;
      acceptCount = 0;
      for (MeasureCount=0; MeasureCount<Parameter_Measurements; MeasureCount++) {
        HMCProp->samplePhiMomentumField();
        HMCProp->sampleOmegaFields();
        if (propagateFields(Parameter_Iterations, Parameter_Epsilon, Parameter_PropTOL, Parameter_FinalTOL)) {
	  if (Parameter_FLAG_WriteConditionNumber==1) {
            writeConditionNumber();
	  }
          measureAndWrite();
          acceptCount++;
  	  TotallyMeasuredConfigurationsCount++;
	  writeCurrentStateDescriptor(1);
        } else {
          if (LogLevel>1) printf("No measurement.\n");
        }
        RunCount++;
      }
      acceptRate = acceptCount / RunCount;
      if (LogLevel>1) printf("-->Update quote: %1.2f\n\n",acceptRate);    
      automaticAdaption(acceptRate);
      if (TotallyMeasuredConfigurationsCount>=Parameter_TotalData) break;
      if (timeOver()) {
        writeCurrentStateDescriptor(0);
        if (LogLevel>1) printf("Time limit reached. ==> EXITING!\n");
        desini();
	exit(0);
      }
    }
  }

  writeCurrentStateDescriptor(0);
  
  if (LogLevel>1) printf("Hybrid Monte Carlo ready. Collected %d data samples!!!\n", TotallyMeasuredConfigurationsCount);

  desini();
}
