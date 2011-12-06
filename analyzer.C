#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

#ifdef UseMPI
  #include <mpi.h>
#endif 


#include "Global.h"
#include "Tools.h"
#include "FermionMatrixOperations.h"
#include "AnalyzerIOControl.h"
#include "StateDescriptorReader.h"
#include "AnalyzerObservable.h"
#include "AnalyzerObservableDetSign.h"
#include "AnalyzerPhiFieldConfiguration.h"
#include "AnalyzerObservableWeight.h"
#include "AnalyzerObservablePsiBarPsiCorr.h"
#include "AnalyzerObservablePsiBarPsiCorrGauged.h"
#include "AnalyzerObservableMagnetizations.h"
#include "AnalyzerObservableGoldstonePropagator.h"
#include "AnalyzerObservableHiggsPropagator.h"
#include "AnalyzerObservableHiggsCorrelator.h"
#include "AnalyzerObservablePsiBarPsiChiralLeftHandedCorr.h"
#include "AnalyzerObservablePsiBarPsiChiralRightHandedCorr.h"
#include "AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr.h"
#include "AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged.h"
#include "AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged.h"
#include "AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged.h"
#include "AnalyzerObservablePsiBarPsiChiralRightHandedCorrGauged.h"
#include "AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged.h"
#include "AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged.h"
#include "AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged.h"
#include "AnalyzerObservableTopBarTopChiralLeftHandedCorrGauged.h"
#include "AnalyzerObservableFermionMatrixConditionNumber.h"
#include "AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning.h"
#include "AnalyzerObservableFermionMatrixSingleMConditionNumberNoPreconditioning.h"
#include "AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning.h"
#include "AnalyzerObservableFermionMatrixSingleMSpectrumPPreconditioning.h"
#include "AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning.h"
#include "AnalyzerObservableHiggsRenormalizedQuarticCoupling.h"
#include "AnalyzerObservablePhiRenormalizedQuarticCoupling.h"
#include "AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling.h"
#include "AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator.h"
#include "AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator.h"
#include "AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator.h"
#include "AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator.h"
#include "AnalyzerObservableFermionMatrixSingleMFullSpectrumNoPreconditioning.h"
#include "AnalyzerObservableFermionMatrixSingleMFullRanSpectrumNoPreconditioning.h"
#include "AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning.h"
#include "AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning.h"
#include "AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning.h"
#include "AnalyzerObservableFermionMatrixConditionNumberPPreconditioning.h"
#include "AnalyzerObservableGaussianWeightEstimate.h"
#include "AnalyzerObservableMultipleTimeScaleIntegration.h"
#include "AnalyzerObservableMultipleTimeScaleIntegration2.h"
#include "AnalyzerObservableMultipleTimeScaleIntegration3.h"
#include "AnalyzerObservableMultipleTimeScaleIntegration4.h"
#include "AnalyzerObservableMultipleTimeScaleIntegrationTest.h"



//Variables
int Parameter_Job_ID = 0;
int Parameter_SD_Selector = 0;
double Parameter_MaxRunTime = 0;
char* Parameter_SD_FileName = NULL;
char** Parameter_tags = NULL;
int Parameter_tagCount = 0;
double startTime = 0;
AnalyzerIOControl* ioControl;
FermionMatrixOperations* fermiOps;
StateDescriptorReader* SDReader;
AnalyzerObservable** analyzerObs = NULL;
Complex** auxVecs = NULL;
int auxVecCount = 0;
int analyzerObsCount = 0;


void startTimer() {
  startTime = zeitwert();
}


double timePassed() {
  double time = (zeitwert()-startTime);
  return time;
}


void randomDelay(double maxTimeInSecs) {
  double timeDelay = maxTimeInSecs * AdvancedZufall(AdvancedSeed);
  if (LogLevel>2) printf("Random Delay: %1.2f secs.\n",timeDelay);
  delay(timeDelay);
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


void loadCommandLineParameters(int argc,char **argv) {
  if (LogLevel>1) printf("Number of arguments = %d\n",argc);
  for (int I=0; I<argc; I++) {
    if (LogLevel>1) printf("Argument %d: %s\n",I+1,argv[I]);  
  }
  
  if (argc<4) {
    printf("Job_ID, MaxRunTime, and StateDescriptor-Selector required!!!\n");
    exit(0);
  }

  bool error = false;
  if (sscanf(argv[1],"%d",&Parameter_Job_ID)!=1) error = true;
  if (sscanf(argv[2],"%lf",&Parameter_MaxRunTime)!=1) error = true;  
  if (sscanf(argv[3],"%d",&Parameter_SD_Selector)!=1) error = true;
  if (error) {
    printf("Parameters could not be read!!!\n");
    exit(0);  
  }

  //Consistency-Check
  if (Parameter_MaxRunTime<=0) error = true;
  if (Parameter_SD_Selector<0) error = true;
  if (error) {
    printf("Invalid Parameter setting !!!\n");
    exit(0);  
  }

  if (LogLevel>1) printf("Parameters: \n");
  if (LogLevel>1) printf("  --> Job_ID           = %d\n", Parameter_Job_ID);
  if (LogLevel>1) printf("  --> MaxRunTime       = %1.2f\n", Parameter_MaxRunTime);
  if (LogLevel>1) printf("  --> SD_Selector      = %d\n", Parameter_SD_Selector);
  
#ifdef UseMPI
  Parameter_Job_ID = 1000*Parameter_Job_ID + ownNodeID;
  Parameter_SD_Selector = nodeCount*Parameter_SD_Selector + ownNodeID;

  if (LogLevel>1) printf("Rescaled Parameters for Parallel-Mode: \n");
  if (LogLevel>1) printf("  --> Job_ID           = %d\n", Parameter_Job_ID);
  if (LogLevel>1) printf("  --> MaxRunTime       = %1.2f\n", Parameter_MaxRunTime);
  if (LogLevel>1) printf("  --> SD_Selector      = %d\n", Parameter_SD_Selector);
#endif   
}


void loadDataFromAnalyzerToDoList() {
  FILE* file = fopen("AnalyzerToDoList.dat","r");
  char* restLine = new char[1000];
  char* dummyStr = new char[1000];  
  Parameter_SD_FileName = new char[1000]; 
  int lineCount = 0;
  Parameter_tags = new char*[1000];
  Parameter_tagCount = 0;
  while (fscanf(file,"%s",dummyStr)==1) {
    fgets(restLine, 1000, file);
    lineCount++;
  }  
  fclose(file);
  
  if (LogLevel>1) printf("AnalyzerToDoList.dat contains %d data lines.\n",lineCount);
  if (lineCount<=0) {
    printf("ERROR: AnalyzerToDoList empty!!!\n");
    exit(0);  
  }  
  
  int jumpOverLines = Parameter_SD_Selector % lineCount;
  file = fopen("AnalyzerToDoList.dat","r");  
  for (int I=0; I<jumpOverLines; I++) {
    fscanf(file,"%s",dummyStr);
    fgets(restLine, 1000, file);    
  }
  fscanf(file,"%s",dummyStr);
  fgets(restLine, 1000, file);    
  snprintf(Parameter_SD_FileName,1000,"%s/data/results/pHMC/states/%s",DataBaseDirectory,dummyStr);
  Parameter_tags[0] = new char[100];
  while (sscanf(restLine,"%s",Parameter_tags[Parameter_tagCount])==1) {
    char* pos = strstr(restLine,Parameter_tags[Parameter_tagCount]);
    for (int I=0; I<(int)strlen(Parameter_tags[Parameter_tagCount]); I++) pos[I] = ' ';
    pos += 1+strlen(Parameter_tags[Parameter_tagCount]);
    Parameter_tagCount++;
    Parameter_tags[Parameter_tagCount] = new char[100];    
  }
  delete[] Parameter_tags[Parameter_tagCount];

  delete[] restLine;
  delete[] dummyStr;
  fclose(file);
  if (LogLevel>1) {
    printf("SD-Filename = %s\n",Parameter_SD_FileName);
    printf("Observable tags...\n");
    for (int I=0; I<Parameter_tagCount; I++) {
      printf("tag%d = %s\n",I,Parameter_tags[I]);
    }  
  }
}


void iniAnalyzerObs() {
  analyzerObsCount = 44;
  analyzerObs = new AnalyzerObservable*[analyzerObsCount];
  analyzerObs[0] = new AnalyzerObservableDetSign(fermiOps, ioControl, SDReader); 
  analyzerObs[1] = new AnalyzerObservableWeight(fermiOps, ioControl, SDReader);
  analyzerObs[2] = new AnalyzerObservablePsiBarPsiCorr(fermiOps, ioControl, SDReader);
  analyzerObs[3] = new AnalyzerObservablePsiBarPsiCorrGauged(fermiOps, ioControl, SDReader);
  analyzerObs[4] = new AnalyzerObservableMagnetizations(fermiOps, ioControl, SDReader); 
  analyzerObs[5] = new AnalyzerObservableGoldstonePropagator(fermiOps, ioControl, SDReader); 
  analyzerObs[6] = new AnalyzerObservableHiggsPropagator(fermiOps, ioControl, SDReader); 
  analyzerObs[7] = new AnalyzerObservableHiggsCorrelator(fermiOps, ioControl, SDReader); 
  analyzerObs[8] = new AnalyzerObservablePsiBarPsiChiralLeftHandedCorr(fermiOps, ioControl, SDReader); 
  analyzerObs[9] = new AnalyzerObservablePsiBarPsiChiralRightHandedCorr(fermiOps, ioControl, SDReader); 
  analyzerObs[10] = new AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorr(fermiOps, ioControl, SDReader); 
  analyzerObs[11] = new AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged(fermiOps, ioControl, SDReader); 
  analyzerObs[12] = new AnalyzerObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged(fermiOps, ioControl, SDReader);   
  analyzerObs[13] = new AnalyzerObservablePsiBarPsiChiralLeftHandedCorrGauged(fermiOps, ioControl, SDReader); 
  analyzerObs[14] = new AnalyzerObservablePsiBarPsiChiralRightHandedCorrGauged(fermiOps, ioControl, SDReader); 
  analyzerObs[15] = new AnalyzerObservablePsiBarPsiChiralLeftHandedCorrRandomGauged(fermiOps, ioControl, SDReader); 
  analyzerObs[16] = new AnalyzerObservablePsiBarPsiChiralRightHandedCorrRandomGauged(fermiOps, ioControl, SDReader);  
  analyzerObs[17] = new AnalyzerObservableBottomBarBottomChiralLeftHandedCorrGauged(fermiOps, ioControl, SDReader); 
  analyzerObs[18] = new AnalyzerObservableTopBarTopChiralLeftHandedCorrGauged(fermiOps, ioControl, SDReader); 
  analyzerObs[19] = new AnalyzerObservableFermionMatrixConditionNumber(fermiOps, ioControl, SDReader); 
  analyzerObs[20] = new AnalyzerObservableFermionMatrixConditionNumberNoPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[21] = new AnalyzerObservableFermionMatrixSingleMConditionNumberNoPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[22] = new AnalyzerObservableFermionMatrixSingleMConditionNumberPPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[23] = new AnalyzerObservableFermionMatrixSingleMSpectrumPPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[24] = new AnalyzerObservableFermionMatrixSingleMSpectrumNoPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[25] = new AnalyzerObservableHiggsRenormalizedQuarticCoupling(fermiOps, ioControl, SDReader); 
  analyzerObs[26] = new AnalyzerObservablePhiRenormalizedQuarticCoupling(fermiOps, ioControl, SDReader); 
  analyzerObs[27] = new AnalyzerObservablePhiRenormalizedUnsymmetricQuarticCoupling(fermiOps, ioControl, SDReader); 
  analyzerObs[28] = new AnalyzerObservableHiggsDefinedOnTimeSliceCorrelator(fermiOps, ioControl, SDReader); 
  analyzerObs[29] = new AnalyzerObservableHiggsGoldstoneUnrotatedCorrelator(fermiOps, ioControl, SDReader); 
  analyzerObs[30] = new AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator(fermiOps, ioControl, SDReader); 
  analyzerObs[31] = new AnalyzerObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator(fermiOps, ioControl, SDReader); 
  analyzerObs[32] = new AnalyzerObservableFermionMatrixSingleMFullSpectrumNoPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[33] = new AnalyzerObservableFermionMatrixSingleMFullRanSpectrumNoPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[34] = new AnalyzerObservableFermionMatrixSingleMLowHighSpectrumNoPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[35] = new AnalyzerObservableFermionMatrixSingleMLowHighSpectrumPPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[36] = new AnalyzerObservableFermionMatrixConditionNumberQHMPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[37] = new AnalyzerObservableFermionMatrixConditionNumberPPreconditioning(fermiOps, ioControl, SDReader); 
  analyzerObs[38] = new AnalyzerObservableGaussianWeightEstimate(fermiOps, ioControl, SDReader); 
  analyzerObs[39] = new AnalyzerObservableMultipleTimeScaleIntegration(fermiOps, ioControl, SDReader); 
  analyzerObs[40] = new AnalyzerObservableMultipleTimeScaleIntegration2(fermiOps, ioControl, SDReader); 
  analyzerObs[41] = new AnalyzerObservableMultipleTimeScaleIntegration3(fermiOps, ioControl, SDReader); 
  analyzerObs[42] = new AnalyzerObservableMultipleTimeScaleIntegration4(fermiOps, ioControl, SDReader); 
  analyzerObs[43] = new AnalyzerObservableMultipleTimeScaleIntegrationTest(fermiOps, ioControl, SDReader); 
}


void iniAuxVecs() {
  auxVecCount = 0;
  for (int I=0; I<analyzerObsCount; I++) {
    int auV = analyzerObs[I]->getNeededAuxVectorCount();
    if (auV > auxVecCount) {
      auxVecCount = auV;
    }
  }
  
  auxVecs = new Complex*[auxVecCount];
  for (int I=0; I<auxVecCount; I++) {
    auxVecs[I] = fermiOps->createFermionVector();
  }
}


void desini() {
  if (LogLevel>2) printf("Desinitializing...\n");
  for (int I=0; I<analyzerObsCount; I++) {
    delete analyzerObs[I];
  }
  delete[] analyzerObs;
  for (int I=0; I<auxVecCount; I++) {
    fermiOps->destroyFermionVector(auxVecs[I]);
  }
  delete[] auxVecs;  
  delete SDReader;
  delete fermiOps;
  delete ioControl;
  delete[] Parameter_SD_FileName;
  for (int I=0; I<Parameter_tagCount; I++) {
    delete[] Parameter_tags[I];
  }
  delete[] Parameter_tags;  
  #ifdef UseMPI
  MPI_Abort(MPI_COMM_WORLD, 0);
  #endif
}


int getNextConfIDforAnalysis() {
  int lowestNextID = -1;

  for (int I=0; I<analyzerObsCount; I++) {
    bool activeObs = false;
    for (int I2=0; I2<Parameter_tagCount; I2++) {
      if (analyzerObs[I]->isNick(Parameter_tags[I2])) activeObs = true;
    }
    if (activeObs) {
      int nextID = analyzerObs[I]->getNextConfIDforAnalysis();
      if (((nextID >= 0) && (nextID < lowestNextID)) || (lowestNextID == -1)) {
        lowestNextID = nextID;
      }
    }
  }
  return lowestNextID;
}


void createFolders() {
  ioControl->createWorkFolder();
  for (int I=0; I<analyzerObsCount; I++) {
    ioControl->createObsFolder(analyzerObs[I]->getObsName());  
  }
}


int main(int argc,char **argv) {
  LogLevel = 3;
  iniMPI(argc, argv);
  fftw_init_threads();
  fftw_plan_with_nthreads(1);  

	printf("ini done\n");
  loadCommandLineParameters(argc,argv);
	printf("loaded command line params\n");
  loadDataFromAnalyzerToDoList();
	printf("loaded data from analyzer\n");

  iniTools(537*Parameter_Job_ID);
  startTimer();
  if (Parameter_Job_ID>0) randomDelay(1000);
	 
	
  SDReader = new StateDescriptorReader(Parameter_SD_FileName); 
  int threadCountPerNode = 1;
  int ParaOpMode = 0;
  char* fftPlanDescriptor = NULL;
  bool usexFFT = false;
  readOptimalFermionVectorEmbeddingAndFFTPlanFromTuningDB(SDReader->getL0(), SDReader->getL1(), SDReader->getL2(), SDReader->getL3(), threadCountPerNode, ParaOpMode, usexFFT, SDReader->getUseP(), SDReader->getUseQ(), SDReader->getUseR(), SDReader->getUseQHM(), fftPlanDescriptor);
  delete[] fftPlanDescriptor;

  char* fileNameExtension = SDReader->getFileNameExtension();
  ioControl = new AnalyzerIOControl(fileNameExtension, Parameter_Job_ID);
  delete[] fileNameExtension;
  fermiOps = new FermionMatrixOperations(SDReader->getL0(), SDReader->getL1(), SDReader->getL2(), SDReader->getL3(), SDReader->getRho(), SDReader->getR(), SDReader->getYN());  
  fermiOps->setMassSplitRatio(SDReader->getMassSplit());
  
  iniAnalyzerObs();
  iniAuxVecs();
  createFolders();
  ioControl->removeDeprecatedInProgressFiles(24.0);

  if (LogLevel>2) printf("\nEntering main Analyzer-Loop.\n");
  int confID = -1;
  while (((confID=getNextConfIDforAnalysis())>=0) && (timePassed()<3600*Parameter_MaxRunTime)) {
    ioControl->markAsInProgress(confID);
    bool keepMarked = false;
    printf("\033[31mSTARTING WITH CONF ID=%d\033[0m\n", confID);
    char* confFileName = ioControl->getPhiConfFileName(confID);
    AnalyzerPhiFieldConfiguration* phiFieldConf = new AnalyzerPhiFieldConfiguration(confFileName, fermiOps); 
    delete[] confFileName;

    if (phiFieldConf->getErrorState() == 0) {
      for (int I=0; I<analyzerObsCount; I++) {
        bool activeObs = false;
        for (int I2=0; I2<Parameter_tagCount; I2++) {
          if (analyzerObs[I]->isNick(Parameter_tags[I2])) activeObs = true;
        }
        if (activeObs) {
          if ((!analyzerObs[I]->isConfAnalyzed(confID)) && (analyzerObs[I]->shallConfBeAnalyzed(confID))) {
	    if (LogLevel>2) printf("Analyzing Observable %s...\n",analyzerObs[I]->getObsName());
            bool b = analyzerObs[I]->analyze(phiFieldConf, auxVecs);
  	    phiFieldConf->clearPhiFieldCopies();
            if (b) {
	      analyzerObs[I]->saveAnalyzerResults(confID);
	    } else {
	      keepMarked = true;
	      if (LogLevel>1) {
	        printf(" *** Could not analyze observable %s on configuration %d. Continue with next configuration.\n",analyzerObs[I]->getObsName(),confID);
	      }
	    }
 	  }
        }
      } 
    } else {
      keepMarked = true;
      if (LogLevel>1) {
        printf(" *** Could not load configuration file %d --> Skipping configuration file\n",confID);
      }
    }
  
    delete phiFieldConf;  
    if (!keepMarked) ioControl->unmarkAsInProgress(confID);
    if (Parameter_Job_ID>0) randomDelay(2);
  }
  if (LogLevel>2) printf("\nLeaving main Analyzer-Loop.\n\n");

  desini();
  desiniTools();
  fftw_cleanup_threads();
  if (LogLevel>1) printf("Analyzer terminated correctly.\n");
}
