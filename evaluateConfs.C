#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>



#include "Global.h"
#include FFTWIncludeFile
#include "Tools.C"
#include "LatexSummary.h"
#include "AnalyzerIOControl.h"
#include "StateDescriptorReader.h"
#include "EvaluateObservable.h"
#include "EvaluateObservableDetSign.h"
#include "EvaluateObservablePsiBarPsiCorr.h"
#include "EvaluateObservablePsiBarPsiCorrGauged.h"
#include "EvaluateObservableMagnetizations.h"
#include "EvaluateObservableGoldstonePropagator.h"
#include "EvaluateObservableWeight.h"
#include "SimulationParameterSet.h"
#include "EvaluateObservableHiggsPropagator.h"
#include "EvaluateObservableHiggsCorrelator.h"
#include "EvaluateObservablePsiBarPsiCorrChiralLeftHanded.h"
#include "EvaluateObservablePsiBarPsiCorrChiralRightHanded.h"
#include "EvaluateObservablePsiBarPsiCorrChiralLeftHandedGauged.h"
#include "EvaluateObservablePsiBarPsiCorrChiralRightHandedGauged.h"
#include "EvaluateObservablePsiBarPsiCorrChiralLeftHandedRandomGauged.h"
#include "EvaluateObservablePsiBarPsiCorrChiralRightHandedRandomGauged.h"
#include "EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorr.h"
#include "EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged.h"
#include "EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged.h"
#include "EvaluateObservableBottomBarBottomCorrChiralLeftHandedGauged.h"
#include "EvaluateObservableTopBarTopCorrChiralLeftHandedGauged.h"
#include "EvaluateObservableFermionMatrixConditionNumber.h"
#include "EvaluateObservableFermionMatrixConditionNumberNoPreconditioning.h"
#include "EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning.h"
#include "EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning.h"
#include "EvaluateObservableFermionMatrixSpectrumPPreconditioning.h"
#include "EvaluateObservableFermionMatrixSpectrumNoPreconditioning.h"
#include "EvaluateObservableHiggsRenormalizedQuarticCoupling.h"
#include "EvaluateObservablePhiRenormalizedQuarticCoupling.h"
#include "EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling.h"
#include "EvaluateObservableHiggsDefinedOnTimeSliceCorrelator.h"
#include "EvaluateObservableHiggsGoldstoneUnrotatedCorrelator.h"
#include "EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator.h"
#include "EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator.h"
#include "EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning.h"
#include "EvaluateObservableFermionMatrixSingleMFullRanSpectrumNoPreconditioning.h"
#include "EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning.h"
#include "EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning.h"
#include "EvaluateObservableFermionMatrixConditionNumberPPreconditioning.h"
#include "EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning.h"
#include "EvaluateObservableGaussianWeightEstimate.h"
#include "EvaluateObservableMultipleTimeScaleIntegration.h"
#include "EvaluateObservableMultipleTimeScaleIntegration2.h"
#include "EvaluateObservableMultipleTimeScaleIntegration3.h"
#include "EvaluateObservableMultipleTimeScaleIntegration4.h"
#include "EvaluateObservableMultipleTimeScaleIntegrationTest.h"



//Variables
int Parameter_SD_Selector = 0;
double Parameter_EvalRelativeStartPosition = 0;
double Parameter_EvalRelativeEndPosition = 1;
char* Parameter_SD_FileName = NULL;
char** Parameter_tags = NULL;
int Parameter_tagCount = 0;
AnalyzerIOControl* ioControl;
StateDescriptorReader* SDReader;
EvaluateObservable** evaluateObs = NULL;
int evaluateObsCount = 0;
char** xmlOutput_Tag;
char** xmlOutput_Description;
double* xmlOutput_Value;
double* xmlOutput_ValueError;  
int xmlOutput_Count;
double GoldstoneRenormFactor;
double GoldstoneRenormFactorError;
bool GoldstoneRenormFactorDetermined;
double physicalScale;
double physicalScaleError;
SimulationParameterSet* SimParaSet;
LatexSummary* latexSummary;


void loadCommandLineParameters(int argc,char **argv) {
  if (LogLevel>1) printf("Number of arguments = %d\n",argc);
  for (int I=0; I<argc; I++) {
    if (LogLevel>1) printf("Argument %d: %s\n",I+1,argv[I]);  
  }
  
  if (argc<2) {
    printf("StateDescriptor-Selector required!!!\n");
    exit(0);
  }

  bool error = false;
  if (sscanf(argv[1],"%d",&Parameter_SD_Selector)!=1) error = true;
  if (argc>=3) {
    if (sscanf(argv[2],"%lf",&Parameter_EvalRelativeStartPosition)!=1) error = true;  
  }
  if (argc>=4) {
    if (sscanf(argv[3],"%lf",&Parameter_EvalRelativeEndPosition)!=1) error = true;  
  }
  
  if (error) {
    printf("Parameters could not be read!!!\n");
    exit(0);  
  }

  //Consistency-Check
  if (Parameter_SD_Selector<0) error = true;
  if (Parameter_EvalRelativeStartPosition<0) error = true;
  if (Parameter_EvalRelativeStartPosition>1) error = true;
  if (Parameter_EvalRelativeEndPosition<0) error = true;
  if (Parameter_EvalRelativeEndPosition>1) error = true;
  if (Parameter_EvalRelativeStartPosition>=Parameter_EvalRelativeEndPosition) error = true;
  
  if (LogLevel>1) printf("Parameters: \n");
  if (LogLevel>1) printf("  --> SD_Selector                    = %d\n", Parameter_SD_Selector);
  if (LogLevel>1) printf("  --> EvalRelativeStartPosition      = %f\n", Parameter_EvalRelativeStartPosition);
  if (LogLevel>1) printf("  --> EvalRelativeEndPosition        = %f\n", Parameter_EvalRelativeEndPosition);
  
  if (error) {
    printf("Invalid Parameter setting !!!\n");
    exit(0);  
  }
}


void loadDataFromEvaluateConfsToDoList() {
  FILE* file = fopen("EvaluateConfsToDoList.dat","r");
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
  
  if (LogLevel>1) printf("EvaluateConfsToDoList.dat contains %d data lines.\n",lineCount);
  if (lineCount<=0) {
    printf("ERROR: EvaluateConfsToDoList empty!!!\n");
    exit(0);  
  }  
  
  int jumpOverLines = Parameter_SD_Selector % lineCount;
  file = fopen("EvaluateConfsToDoList.dat","r");  
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


void ini() {
  LogLevel = 2;
  fftw_init_threads();   //<-- Desini aufrufen
  fftw_plan_with_nthreads(1);  

  iniTools(5517+Parameter_SD_Selector*12345);
  
  physicalScale = NaN;
  physicalScaleError = NaN;
  GoldstoneRenormFactor = 1.0;
  GoldstoneRenormFactorError = 0.0;
  GoldstoneRenormFactorDetermined = false;
  xmlOutput_Tag = new char*[10000];
  xmlOutput_Description = new char*[10000];
  xmlOutput_Value = new double[10000];
  xmlOutput_ValueError = new double[10000]; 
  xmlOutput_Count = 0;
}


void iniEvaluateObs() {
  evaluateObsCount = 44;
  EvaluateObservable* evalWeight = NULL;
  EvaluateObservable* evalDetSign = NULL;  
  evaluateObs = new EvaluateObservable*[evaluateObsCount];
  evaluateObs[0] = new EvaluateObservableWeight(ioControl, SDReader, NULL, NULL, 0, 1); 

  for (int I=0; I<Parameter_tagCount; I++) {
    if (evaluateObs[0]->isNick(Parameter_tags[I])) evalWeight = evaluateObs[0];
  }
  
  evaluateObs[1] = new EvaluateObservableDetSign(ioControl, SDReader, evalWeight, NULL, 0, 1); 
  
  for (int I=0; I<Parameter_tagCount; I++) {
    if (evaluateObs[1]->isNick(Parameter_tags[I])) evalDetSign = evaluateObs[1];
  }
  
  evaluateObs[2] = new EvaluateObservableGoldstonePropagator(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[3] = new EvaluateObservableMagnetizations(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[4] = new EvaluateObservablePsiBarPsiCorr(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[5] = new EvaluateObservablePsiBarPsiCorrGauged(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[6] = new EvaluateObservableHiggsPropagator(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[7] = new EvaluateObservableHiggsCorrelator(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[8] = new EvaluateObservablePsiBarPsiCorrChiralLeftHanded(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[9] = new EvaluateObservablePsiBarPsiCorrChiralRightHanded(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[10] = new EvaluateObservablePsiBarPsiCorrChiralLeftHandedGauged(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[11] = new EvaluateObservablePsiBarPsiCorrChiralRightHandedGauged(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[12] = new EvaluateObservablePsiBarPsiCorrChiralLeftHandedRandomGauged(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[13] = new EvaluateObservablePsiBarPsiCorrChiralRightHandedRandomGauged(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[14] = new EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorr(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[15] = new EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrGauged(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[16] = new EvaluateObservableTwoFermionPsiBarPsiChiralLeftHandedCorrRandomGauged(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[17] = new EvaluateObservableBottomBarBottomCorrChiralLeftHandedGauged(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[18] = new EvaluateObservableTopBarTopCorrChiralLeftHandedGauged(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[19] = new EvaluateObservableFermionMatrixConditionNumber(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[20] = new EvaluateObservableFermionMatrixConditionNumberNoPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[21] = new EvaluateObservableFermionMatrixSingleMConditionNumberPPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[22] = new EvaluateObservableFermionMatrixSingleMConditionNumberNoPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[23] = new EvaluateObservableFermionMatrixSpectrumPPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[24] = new EvaluateObservableFermionMatrixSpectrumNoPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[25] = new EvaluateObservableHiggsRenormalizedQuarticCoupling(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[26] = new EvaluateObservablePhiRenormalizedQuarticCoupling(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[27] = new EvaluateObservablePhiRenormalizedUnsymmetricQuarticCoupling(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[28] = new EvaluateObservableHiggsDefinedOnTimeSliceCorrelator(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[29] = new EvaluateObservableHiggsGoldstoneUnrotatedCorrelator(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[30] = new EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK0Correlator(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition); 
  evaluateObs[31] = new EvaluateObservableHiggsGoldstoneUnrotatedTwoParticleK012Correlator(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[32] = new EvaluateObservableFermionMatrixSingleMFullSpectrumNoPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[33] = new EvaluateObservableFermionMatrixSingleMFullRanSpectrumNoPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[34] = new EvaluateObservableFermionMatrixLowHighSpectrumNoPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[35] = new EvaluateObservableFermionMatrixLowHighSpectrumPPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[36] = new EvaluateObservableFermionMatrixConditionNumberPPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[37] = new EvaluateObservableFermionMatrixConditionNumberQHMPreconditioning(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[38] = new EvaluateObservableGaussianWeightEstimate(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[39] = new EvaluateObservableMultipleTimeScaleIntegration(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[40] = new EvaluateObservableMultipleTimeScaleIntegration2(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[41] = new EvaluateObservableMultipleTimeScaleIntegration3(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[42] = new EvaluateObservableMultipleTimeScaleIntegration4(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   
  evaluateObs[43] = new EvaluateObservableMultipleTimeScaleIntegrationTest(ioControl, SDReader, evalWeight, evalDetSign, Parameter_EvalRelativeStartPosition, Parameter_EvalRelativeEndPosition);   

    
  for (int I=0; I<evaluateObsCount; I++) {
    for (int I2=0; I2<evaluateObsCount; I2++) {
      evaluateObs[I]->addOtherObsPointer(evaluateObs[I2]);
    }
  }
  for (int I=0; I<evaluateObsCount; I++) {
    evaluateObs[I]->defineObsDependencies();
  }
  bool change = true;
  while (change) {
    change = false;
    for (int I=0; I<evaluateObsCount; I++) {
      if (evaluateObs[I]->reduceDataDueToDependencies()) {
        change = true;
      }
    }
  }

  //Initialisiere LatexSummary
  char** obsNames = new char*[evaluateObsCount];
  char** LatexBodyLocalFileNames = new char*[evaluateObsCount];
  char* LatexSummaryBaseDirName = ioControl->getWorkFolderName();
  char* LatexFileName = ioControl->getLatexSummaryBaseName();
  for (int I=0; I<evaluateObsCount; I++) {
    obsNames[I] = new char[1000];
    snprintf(obsNames[I], 1000, "%s",evaluateObs[I]->getObsName());
    LatexBodyLocalFileNames[I] = ioControl->getObservableLatexBodyLocalFileName(evaluateObs[I]->getObsName());
  }

  latexSummary = new LatexSummary(LatexFileName, LatexSummaryBaseDirName, evaluateObsCount, obsNames, LatexBodyLocalFileNames);
  
  for (int I=0; I<evaluateObsCount; I++) {
    delete[] obsNames[I];
    delete[] LatexBodyLocalFileNames[I];
  }  
  delete[] obsNames;
  delete[] LatexSummaryBaseDirName;
  delete[] LatexBodyLocalFileNames;  
  delete[] LatexFileName;
}


void desini() {
  if (LogLevel>1) printf("\nDesinitializing...\n");
  for (int I=0; I<xmlOutput_Count; I++) {
    delete[] xmlOutput_Tag[I];
    delete[] xmlOutput_Description[I];
  }
  delete[] xmlOutput_Value;
  delete[] xmlOutput_ValueError;
  delete[] xmlOutput_Tag;
  delete[] xmlOutput_Description;
  xmlOutput_Count = 0;
  
  for (int I=0; I<evaluateObsCount; I++) {
    delete evaluateObs[I];
  }
  delete[] evaluateObs;
    
  delete latexSummary;    
      
  delete SDReader;
  delete ioControl;
  delete SimParaSet;
  delete[] Parameter_SD_FileName;
  for (int I=0; I<Parameter_tagCount; I++) {
    delete[] Parameter_tags[I];
  }
  delete[] Parameter_tags;  
}


void addXmlOutput(char* tag, char* des, double val, double err) {
  xmlOutput_Tag[xmlOutput_Count] = cloneString(tag);
  xmlOutput_Description[xmlOutput_Count] = cloneString(des);
  xmlOutput_Value[xmlOutput_Count] = val;
  xmlOutput_ValueError[xmlOutput_Count] = err;
  xmlOutput_Count++;
}


void writeSummaryXmlFile() {
  if (LogLevel>2) printf("Writing SummaryXmlFile to disk...");
  char* sumFileName = ioControl->getMainXmlFileName();  
  FILE* sumFile = fopen(sumFileName,"w");
  delete[] sumFileName;
  
  fprintf(sumFile,"<?xml version=\'1.0\'?>\n");
  char* ext = SDReader->getFileNameExtension();
  char* mainTag = new char[1000];
  snprintf(mainTag,1000,"Summary%s",ext);
  delete[] ext;
  fprintf(sumFile,"<%s>\n",mainTag);    
  fprintf(sumFile,"<General>\n");  
  for (int I=0; I<xmlOutput_Count; I++) {
    fprintf(sumFile,"  <%s>\n",xmlOutput_Tag[I]);
    if ((!isInteger(xmlOutput_Value[I])) || (isNaN(xmlOutput_ValueError[I])) || (xmlOutput_ValueError[I]!=0)) {
      fprintf(sumFile,"    <value>%1.15e</value>\n",xmlOutput_Value[I]);
      if ((isNaN(xmlOutput_ValueError[I])) || (xmlOutput_ValueError[I]!=0)) {
        fprintf(sumFile,"    <error>%1.15e</error>\n",xmlOutput_ValueError[I]);    
      }
    } else {
      fprintf(sumFile,"    <value>%d</value>\n",roundToInt(xmlOutput_Value[I]));    
    }
    fprintf(sumFile,"    <description>%s</description>\n",xmlOutput_Description[I]);
    fprintf(sumFile,"  </%s>\n",xmlOutput_Tag[I]);    
  }
  fprintf(sumFile,"</General>\n");  
  
  for (int I=0; I<evaluateObsCount; I++) { 
    char* obsFileName = ioControl->getObservableXmlFileName(evaluateObs[I]->getObsName());    
    FILE* obsFile = fopen(obsFileName,"r");
    delete[] obsFileName;
    if (obsFile!=NULL) {
      char* line = new char[10000];
      fgets(line, 10000, obsFile);

      while (!feof(obsFile)) {
        fgets(line, 10000, obsFile);
	if (!feof(obsFile)) fputs(line, sumFile);        
      }
      fprintf(sumFile,"\n");

      delete[] line;
      fclose(obsFile);
    }
  }
  fprintf(sumFile,"</%s>\n",mainTag);    
  delete[] mainTag;
  fclose(sumFile);
  if (LogLevel>2) printf("ready!\n");
}


void addXML_And_LatexSummaryInfoTableLine(char* xmltag, char* des, char* shortCut, double val, double error, char* unit, char* valFormat) {
  char* line = new char[2000];
  char* formatStr = new char[1000];
  
  if ((error!=0) && (unit!=NULL)) {
    snprintf(formatStr,1000,"%s & %s & ($%s \\pm %s$)\\,%s \\\\ \n",des,shortCut,valFormat,valFormat, unit);
    snprintf(line,2000,formatStr,val,error);
  }
  if ((error==0) && (unit!=NULL)) {
    snprintf(formatStr,1000,"%s & %s & %s \\,%s \\\\ \n",des,shortCut,valFormat, unit);
    snprintf(line,2000,formatStr,val,error);
  }
  if ((error!=0) && (unit==NULL)) {
    snprintf(formatStr,1000,"%s & %s & $%s \\pm %s$ \\\\ \n",des,shortCut,valFormat,valFormat);
    snprintf(line,2000,formatStr,val,error);
  }
  if ((error==0) && (unit==NULL)) {
    snprintf(formatStr,1000,"%s & %s & %s \\\\ \n",des,shortCut,valFormat);
    snprintf(line,2000,formatStr,val,error);
  }
  
  latexSummary->addDirectTextBeforeIncludes(line);
  
  addXmlOutput(xmltag, des, val, error);
  
  delete[] formatStr;
  delete[] line;
}


void composeGeneralSummaryInfo() {
  latexSummary->addDirectTextBeforeIncludes("\\section{General Information}\n\n");

  latexSummary->addDirectTextBeforeIncludes("\\begin{center}\n");
  latexSummary->addDirectTextBeforeIncludes("\\begin{tabular}{lcr}\n");
  
  addXML_And_LatexSummaryInfoTableLine("PhysicalScale", "Inverse lattice spacing in GeV","$a^{-1}$", physicalScale, physicalScaleError, "GeV", "%1.1f");
  addXML_And_LatexSummaryInfoTableLine("L0", "Lattice size L0","L0", SDReader->getL0(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("L1", "Lattice size L1","L1", SDReader->getL1(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("L2", "Lattice size L2","L2", SDReader->getL2(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("L3", "Lattice size L3","L3", SDReader->getL3(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("Nf", "Number of fermion generations","$N_f$", SDReader->getNf(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("Y0", "Bare Yukawa coupling in continuum notation","$y_0$", SimParaSet->getY0(), 0, NULL, "%1.5f");
  addXML_And_LatexSummaryInfoTableLine("Lambda0", "Bare quartic self-coupling in continuum notation","$\\lambda_0$", SimParaSet->getLambda0(), 0, NULL, "%1.5f");
  addXML_And_LatexSummaryInfoTableLine("M0Squared", "Squared bare Higgs mass in continuum notation","$m_0^2$", SimParaSet->getM0Squared(), 0, NULL, "%1.5f");
  addXML_And_LatexSummaryInfoTableLine("YN", "Bare Yukawa coupling in Nf-notation","$y_N$", SimParaSet->getYN(), 0, NULL, "%1.5f");
  addXML_And_LatexSummaryInfoTableLine("LambdaN", "Bare quartic self-coupling in Nf-notation","$\\lambda_N$", SimParaSet->getLambdaN(), 0, NULL, "%1.5f");
  addXML_And_LatexSummaryInfoTableLine("KappaN", "Bare Higgs-Hopping parameter kappa in Nf-notation","$\\kappa_N$", SimParaSet->getKappaN(), 0, NULL, "%1.6f");
  addXML_And_LatexSummaryInfoTableLine("MassSplit", "Ratio of Yukawa couplings for bottom and top","$\\frac{y_B}{y_T}$", SDReader->getMassSplit(), 0, NULL, "%1.5f");
  addXML_And_LatexSummaryInfoTableLine("ExternalCurrent", "Bare external current in lattice units","$J_0$", SDReader->getExternalCurrent(), 0, NULL, "%1.5f");


  addXML_And_LatexSummaryInfoTableLine("Rho", "Rho-Parameter of the Neuberger-Dirac operator","$\\rho$", SDReader->getRho(), 0, NULL, "%1.3f");
  addXML_And_LatexSummaryInfoTableLine("R", "R-Parameter of the Neuberger-Dirac operator","$r$", SDReader->getR(), 0, NULL, "%1.3f");
  addXML_And_LatexSummaryInfoTableLine("ExplicitFermiMass", "Explicit Fermion Mass","$m_F$", SDReader->getExplicitFermionMass(), 0, NULL, "%1.5f");
  addXML_And_LatexSummaryInfoTableLine("Alpha", "Exponent of MMdagger","$\\alpha$", SDReader->getAlpha(), 0, NULL, "%1.3f");
  addXML_And_LatexSummaryInfoTableLine("RandomGauge", "Use random gauging","", SDReader->useRandomGauge(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("ExactReweight", "Use exact reweighing","", SDReader->useExactReweighing(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("DirectSampling", "Use direct sampling of omega-fields","", SDReader->useDirectSamplingOfOmegaFields(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("Theta", "Theta-Parameter for omega-propagation","$\\theta$", SDReader->getTheta(), 0, NULL, "%1.3f");
  addXML_And_LatexSummaryInfoTableLine("TotalDataCount", "Total amount of independent measured data","", SDReader->getPerformedMeasuringSteps(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("IntStepSize", "Step size integrator for polynomial P0 ","$\\epsilon_{int}$", SDReader->getMolecularDynamicsEpsilon(), 0, NULL, "%1.5f");  
  addXML_And_LatexSummaryInfoTableLine("IntIterP0", "Iterations of integrator for polynomial P0","$N_{P_0}$", SDReader->getPolynomIterations_P0(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("IntIterHiggs", "Iterations of integrator for Higgs-dynamics","$N_H$", SDReader->getHiggsDynamicsIterations(), 0, NULL, "%1.0f");    
  addXML_And_LatexSummaryInfoTableLine("IntTraLen", "Trajectory length for integrators","$N_{P_0}\\epsilon_{P_0}$", SDReader->getPolynomIterations_P0()*SDReader->getMolecularDynamicsEpsilon(), 0, NULL, "%1.3f");  
  addXML_And_LatexSummaryInfoTableLine("PolLowerBoundP0", "Lower bound for polynomial P0","$\\epsilon_{P_0}$", SDReader->getPolynomLowerBound_P0(), 0, NULL, "%1.2e");
  addXML_And_LatexSummaryInfoTableLine("PolUpperBound", "Lower bound for all polynomials","$\\lambda_{pol}$", SDReader->getPolynomUpperBound(), 0, NULL, "%1.3f");
  addXML_And_LatexSummaryInfoTableLine("PolDegP0", "Degree of polynomial P0","", SDReader->getPolynomDegree_P0(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("UpperEWSafetyFac", "Upper EW-safety factor","", SDReader->getUpperEWBoundSafetyFactor(), 0, NULL, "%1.2f");
  addXML_And_LatexSummaryInfoTableLine("ModelSelection", "Luescher-Model selected","", SDReader->useModelSelection(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("FourierAccType", "Fourier-Acceleration Type","", SDReader->getFourierAccelerationType(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("FourierAccPara", "Fourier-Acceleration Parameter","", SDReader->getFourieAccelerationParameter(), 0, NULL, "%1.3f");

  latexSummary->addDirectTextBeforeIncludes("\\end{tabular}\n");
  latexSummary->addDirectTextBeforeIncludes("\\end{center}\n");


  //Next Table
  latexSummary->addDirectTextBeforeIncludes("\\begin{center}\n");
  latexSummary->addDirectTextBeforeIncludes("\\begin{tabular}{lcr}\n");
  
  addXML_And_LatexSummaryInfoTableLine("QuasiHermitMode", "Usage of Quasi-Hermitean-Mode","", SDReader->getUseQHM(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("UseP", "Usage of preconditioner P","", SDReader->getUseP(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("UseQ", "Usage of preconditioner Q","", SDReader->getUseQ(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("UseR", "Usage of preconditioner R","", SDReader->getUseR(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("PPrecM", "P-Preconditioner parameter M","", SDReader->getPPrecondParameterM(), 0, NULL, "%1.3f");
  addXML_And_LatexSummaryInfoTableLine("PPrecS", "P-Preconditioner parameter S","", SDReader->getPPrecondParameterS(), 0, NULL, "%1.3f");
  addXML_And_LatexSummaryInfoTableLine("PPrecM", "R-Preconditioner parameter M","", SDReader->getRPrecondParameterM(), 0, NULL, "%1.3f");
  addXML_And_LatexSummaryInfoTableLine("PPrecF", "R-Preconditioner parameter F","", SDReader->getRPrecondParameterF(), 0, NULL, "%1.3f");
  
  latexSummary->addDirectTextBeforeIncludes("\\end{tabular}\n");
  latexSummary->addDirectTextBeforeIncludes("\\end{center}\n");


  //Next Table
  latexSummary->addDirectTextBeforeIncludes("\\begin{center}\n");
  latexSummary->addDirectTextBeforeIncludes("\\begin{tabular}{lcr}\n");
  
  addXML_And_LatexSummaryInfoTableLine("SphIntMode", "Spherical Integation Mode","", SDReader->getSphericalHiggsIntegrationMode(), 0, NULL, "%1.0f");
  addXML_And_LatexSummaryInfoTableLine("SphIntZeta", "Spherical Integration parameter","", SDReader->getZetaForHiggsIntegrationMode(), 0, NULL, "%1.6f");
  addXML_And_LatexSummaryInfoTableLine("ModExtC6", "Model Extension Parameter C6","", SDReader->getModelParameterC6(), 0, NULL, "%1.6f");
  addXML_And_LatexSummaryInfoTableLine("ModExtC8", "Model Extension Parameter C8","", SDReader->getModelParameterC8(), 0, NULL, "%1.6f");
  addXML_And_LatexSummaryInfoTableLine("ModExtC10", "Model Extension Parameter C10","", SDReader->getModelParameterC10(), 0, NULL, "%1.6f");
  
  addXML_And_LatexSummaryInfoTableLine("ModExtLam6", "Model Extension Parameter Lambda6","", SDReader->getModelParameterLambda6(), 0, NULL, "%1.6f");
  addXML_And_LatexSummaryInfoTableLine("ModExtLam8", "Model Extension Parameter Lambda8","", SDReader->getModelParameterLambda8(), 0, NULL, "%1.6f");
  addXML_And_LatexSummaryInfoTableLine("ModExtLam10", "Model Extension Parameter Lambda10","", SDReader->getModelParameterLambda10(), 0, NULL, "%1.6f");
  
  latexSummary->addDirectTextBeforeIncludes("\\end{tabular}\n");
  latexSummary->addDirectTextBeforeIncludes("\\end{center}\n");
 

  latexSummary->addDirectTextBeforeIncludes("\\clearpage\n");  
}


void createAndPrepareFolders() {
  ioControl->createWorkFolder();
  for (int I=0; I<evaluateObsCount; I++) {
    ioControl->createAndPrepareObsPlotsFolder(evaluateObs[I]->getObsName());  
  }
}


void afterEvaluateControl(EvaluateObservable* obs) {
  if (obs->isObs("GoldstonePropagator")) {
    EvaluateObservableGoldstonePropagator* obsGProp = (EvaluateObservableGoldstonePropagator*) obs;
    GoldstoneRenormFactor = obsGProp->getPropagatorZFactor();
    GoldstoneRenormFactorError = obsGProp->getPropagatorZFactorError();
    if ((GoldstoneRenormFactor==1.0) && (GoldstoneRenormFactorError==0.0)) {
      GoldstoneRenormFactorDetermined = false;
    } else {
      GoldstoneRenormFactorDetermined = true;    
    }
  }
  
  if (obs->isObs("Magnetizations")) {
    EvaluateObservableMagnetizations* obsMag= (EvaluateObservableMagnetizations*) obs;
    double rescale = SimParaSet->reparametrize_HiggsField(SimulationParameterSet_ContinuumNotation);
    double vev = rescale*obsMag->getMagnetizationM();
    double vevError = rescale*obsMag->getMagnetizationMError();

    physicalScale = Physical_VEV_GeV*sqrt(GoldstoneRenormFactor) / vev;
    physicalScaleError = physicalScale*sqrt(sqr(0.5*GoldstoneRenormFactorError/GoldstoneRenormFactor) + sqr(vevError/vev));   
  }
}


int main(int argc,char **argv) {
  loadCommandLineParameters(argc,argv);

  ini();

  loadDataFromEvaluateConfsToDoList();

  SDReader = new StateDescriptorReader(Parameter_SD_FileName); 
  int threadCountPerNode = 1;
  int ParaOpMode = 0;
  char* fftPlanDescriptor = NULL;
  bool usexFFT = false;
  readOptimalFermionVectorEmbeddingAndFFTPlanFromTuningDB(SDReader->getL0(), SDReader->getL1(), SDReader->getL2(), SDReader->getL3(), threadCountPerNode, ParaOpMode, usexFFT, SDReader->getUseP(), SDReader->getUseQ(), SDReader->getUseR(), SDReader->getUseQHM(), fftPlanDescriptor);
  delete[] fftPlanDescriptor;
  char* fileNameExtension = SDReader->getFileNameExtension();
  ioControl = new AnalyzerIOControl(fileNameExtension, -1);
  delete[] fileNameExtension;
  SimParaSet = new SimulationParameterSet(SDReader, SimulationParameterSet_NfNotation);

  iniEvaluateObs();
  createAndPrepareFolders();

  if (LogLevel>1) printf("\nEntering main Evaluate-Loop.\n");
  bool* evalSuccess = new bool[evaluateObsCount];
  for (int I=0; I<evaluateObsCount; I++) {
    evalSuccess[I] = false;
    bool activeObs = false;
    for (int I2=0; I2<Parameter_tagCount; I2++) {
      if (evaluateObs[I]->isNick(Parameter_tags[I2])) activeObs = true;
    }
    if (activeObs) {
      if (LogLevel>1) printf("Evaluating Observable %s...\n",evaluateObs[I]->getObsName());
      evaluateObs[I]->setPhysicalScale(physicalScale, physicalScaleError, GoldstoneRenormFactorDetermined);
      evalSuccess[I] = evaluateObs[I]->evaluateWrapper();
      afterEvaluateControl(evaluateObs[I]);
      if (evalSuccess[I]) {
	evaluateObs[I]->generateLatexAndPlotsAndXMLWrapper();
	evaluateObs[I]->writeXmlOutput();
      } else {
        if (LogLevel>1) printf(" *** Could not evaluate observable %s. Continue with next observable.\n",evaluateObs[I]->getObsName());
      }
    } 
  }
  if (LogLevel>1) printf("\nLeaving main Evaluate-Loop.\n");
  composeGeneralSummaryInfo();
  writeSummaryXmlFile();
  
  if (LogLevel>1) {
    printf("\n\n *** Evaluation Status Summary: ***\n\n");
    for (int I=0; I<evaluateObsCount; I++) {
      if (evalSuccess[I]) {
        bool activeObs = false;
        for (int I2=0; I2<Parameter_tagCount; I2++) {
          if (evaluateObs[I]->isNick(Parameter_tags[I2])) activeObs = true;
        }
        if (activeObs) {
          evaluateObs[I]->printStatusSummaryToScreen();
        }
      }
    }
  }
  delete[] evalSuccess;

  desini();
  desiniTools();
  fftw_cleanup_threads();
  if (LogLevel>1) printf("\n   *** EvaluateConfs terminated correctly. ***\n\n");
}
