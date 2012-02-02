#include "AnalyzerObservable.h"


AnalyzerObservable::AnalyzerObservable(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* sdr, char* oName, char* nick) { 
  if (LogLevel>1) printf("AnalyzerObservable %s (%s) initializing ...\n",oName,nick);
  ioControl = aIOcon;
  fermiOps = fOps;
  SDReader = sdr;
  nickName = new char[1+strlen(nick)];
  snprintf(nickName, 1+strlen(nick), "%s", nick);
  obsName = new char[1+strlen(oName)];
  snprintf(obsName, 1+strlen(oName), "%s", oName);  
  totalPhiConfCount = ioControl->getTotalPhiConfCount();
  if (LogLevel>2) printf("Total number of configurations on disk: %d\n",totalPhiConfCount);
  analyzerResults = NULL;
  analyzeEveryXXXconf = 1;
  SimParaSet = new SimulationParameterSet(SDReader, SimulationParameterSet_NfNotation);
}


AnalyzerObservable::~AnalyzerObservable() {
  if (LogLevel>2) printf("Desinitializing Observable %s...",obsName);
  delete[] nickName;
  delete[] obsName;
  delete[] analyzerResults;  
  delete SimParaSet;
  if (LogLevel>2) printf("sucessfully\n");
}


void AnalyzerObservable::ini(int resCount) {
  if (LogLevel>2) printf("Reserving %d doubles for result-storage...\n", resCount);
  analyzerResults = new double[resCount];
  for (int I=0; I<resCount; I++) analyzerResults[I] = NaN;
}


void AnalyzerObservable::saveAnalyzerResults(int ID) {
  char* fileName = ioControl->getObservableFileName(obsName);
  FILE* file = fopen(fileName,"a");
  delete[] fileName;
  
  fprintf(file, "%d ",ID);
  int count = getAnalyzerResultsCount();
  if (count<100) {
    for (int I=0; I<count; I++) {
      fprintf(file, "%1.15e ", analyzerResults[I]);
    }
  } else {
    fprintf(file, "%d ", ioControl->getJobID());
    char* outFileName = ioControl->getObservableOutsourceFileName(obsName, ID);
    std::fstream outFile;
    outFile.open(outFileName, std::ios::out);
    outFile.write((char*)analyzerResults, 8*count);
    outFile.flush();
    outFile.close();
    delete[] outFileName;
  }
  fprintf(file, "\n");
  
  fclose(file);
}


bool AnalyzerObservable::shallConfBeAnalyzed(int ID) {
  if ((ID % analyzeEveryXXXconf) == 0) {
    return true;
  } else {
    return false;
  }
}


bool AnalyzerObservable::isConfAnalyzed(int ID) {
  char** fileNames = NULL;
  int fileCount = 0;
  char* restLine = new char[10000];
  bool isAna = false;

  ioControl->getObservableFileNameList(obsName, fileNames, fileCount);
  for (int I=0; I<fileCount; I++) {
    FILE* file = fopen(fileNames[I],"r");  
    int id = -1;
    while (fscanf(file,"%d",&id)==1) {
      fgets(restLine, 10000, file);
      
      if (id == ID) {
        isAna = true;
	I = fileCount;
	break;
      }
    }
  
    fclose(file);
  }
  deleteFileNameList(fileNames, fileCount);
  delete[] restLine;
  
  return isAna;
}


int AnalyzerObservable::getNextConfIDforAnalysis() {
  char** fileNames = NULL;
  int fileCount = 0;
  char* restLine = new char[10000];
  int statusDataMAX = 1000000;
  int* statusData = new int[statusDataMAX];

  for (int I=0; I<statusDataMAX; I++) statusData[I] = 0;
  ioControl->getObservableFileNameList(obsName, fileNames, fileCount);
  for (int I=0; I<fileCount; I++) {
    FILE* file = fopen(fileNames[I],"r");  
    int id = -1;
    while (fscanf(file,"%d",&id)==1) {
      fgets(restLine, 10000, file);
      
      if (id>=statusDataMAX) {
        printf("ERROR: AnalyzerObservable::getNextConfIDforAnalysis: ID=%d higher than %d\n",id,statusDataMAX);
	exit(0);
      }
      statusData[id] = 1;
    }
  
    fclose(file);
  }
  deleteFileNameList(fileNames, fileCount);
  delete[] restLine;
  
  int confInProgressIDcount = 0;
  int* confInProgressIDs = NULL; 
  ioControl->getConfInProgressIDs(confInProgressIDcount, confInProgressIDs);
  for (int I=0; I<confInProgressIDcount; I++) {
    int id = confInProgressIDs[I];
    if (id>=statusDataMAX) {
      printf("ERROR: AnalyzerObservable::getNextConfIDforAnalysis (II): ID=%d higher than %d\n",id,statusDataMAX);
      exit(0);
    }
    statusData[id]++;
  }
  delete[] confInProgressIDs;  
  
  int nextID = -1;
  for (int I=1; I<statusDataMAX; I++) {
    if ((statusData[I] == 0) && (shallConfBeAnalyzed(I))) {
      nextID = I;
      break;
    }
  }
  if (nextID>totalPhiConfCount) nextID = -1;
  
  delete[] statusData;
  return nextID;
}


bool AnalyzerObservable::isNick(char* nick) {
  if (strlen(nick) != strlen(nickName)) return false;
  
  for (int I=0; I<(int)strlen(nick); I++) {
    if (nick[I] != nickName[I]) return false;
  }

  return true;
}

  
char* AnalyzerObservable::getObsName() {
  return obsName;
}


char* AnalyzerObservable::getNickName() {
  return nickName;
}
