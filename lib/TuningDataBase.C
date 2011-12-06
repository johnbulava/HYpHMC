#include "TuningDataBase.h"
#include "Global.h"
#include <string.h>

TuningDataBase::TuningDataBase(char* dbFname) { 
  DBFileName = new char[2000];
  snprintf(DBFileName,2000, "%s", dbFname);
  FILE* file = fopen(DBFileName, "r");
  if (LogLevel>2) printf("TuningDataBase initialized with DBname = %s\n", DBFileName);
  if (file == NULL) {
    if (LogLevel>2) printf("  ==> TuningDataBase does not exist yet!!!\n");
  } else {
    fclose(file);
    if (LogLevel>2) printf("  ==> TuningDataBase already exists!!!\n");
  }
}


TuningDataBase::~TuningDataBase() { 
  delete[] DBFileName;
}


void TuningDataBase::freeDBEntry(DBEntryType &entry) {
  delete entry.hostName;
  entry.hostName = NULL;
  delete entry.fftPlanDescriptor;
  entry.fftPlanDescriptor = NULL;
}


bool TuningDataBase::readNextEntry(FILE* file, DBEntryType &entry) {
  entry.hostName = new char[4];
  entry.fftPlanDescriptor = new char[10000];
  
  if (file==NULL) return false;
  
  if (fscanf(file,"%s", entry.hostName) != 1) return false;
  if (fscanf(file,"%d", &(entry.procCount)) != 1) return false;
  if (fscanf(file,"%d", &(entry.L0)) != 1) return false;
  if (fscanf(file,"%d", &(entry.L1)) != 1) return false;
  if (fscanf(file,"%d", &(entry.L2)) != 1) return false;
  if (fscanf(file,"%d", &(entry.L3)) != 1) return false;
  if (fscanf(file,"%d", &(entry.paraOpMode)) != 1) return false;
  if (fscanf(file,"%d", &(entry.threadPerNodeCount)) != 1) return false;
  if (fscanf(file,"%d", &(entry.xFFT)) != 1) return false;
  if (fscanf(file,"%d", &(entry.useP)) != 1) return false;
  if (fscanf(file,"%d", &(entry.useQ)) != 1) return false;
  if (fscanf(file,"%d", &(entry.useR)) != 1) return false;
  if (fscanf(file,"%d", &(entry.QHM)) != 1) return false;
  if (fscanf(file,"%d", &(entry.xtrSize1)) != 1) return false;
  if (fscanf(file,"%d", &(entry.xtrSize2)) != 1) return false;
  if (fscanf(file,"%d", &(entry.xtrSize3)) != 1) return false;
  if (fscanf(file,"%d", &(entry.runCount)) != 1) return false;
  if (fscanf(file,"%lf", &(entry.deltaT)) != 1) return false;
  fgets(entry.fftPlanDescriptor, 10000, file);
  
  int SpaceCount = 0;
  while ((entry.fftPlanDescriptor[SpaceCount]==' ') && (SpaceCount<((int)strlen(entry.fftPlanDescriptor)))) {
    SpaceCount++;
  }
  for (int I=0; I<=((int)strlen(entry.fftPlanDescriptor)-SpaceCount); I++) {
    entry.fftPlanDescriptor[I] = entry.fftPlanDescriptor[I+SpaceCount];
  }
  
  return true;
}

    
void TuningDataBase::writeDBEntry(char* hostName, int procCount, int L0, int L1, int L2, int L3, int paraOpMode, int threadPerNodeCount, int xFFT, int useP, int useQ, int useR, int QHM, int xtrSize1, int xtrSize2, int xtrSize3, int runCount, double deltaT, char* fftPlanDescriptor) {
  FILE* file = fopen(DBFileName, "a");
  int hnlen = strlen(hostName);
  if (hnlen>4) hnlen=4;
  for (int I=0; I<hnlen; I++) fprintf(file,"%c",hostName[I]);
  fprintf(file," %d",procCount);
  fprintf(file," %d",L0);
  fprintf(file," %d",L1);
  fprintf(file," %d",L2);
  fprintf(file," %d",L3);
  fprintf(file," %d",paraOpMode);
  fprintf(file," %d",threadPerNodeCount);
  fprintf(file," %d",xFFT);
  fprintf(file," %d",useP);
  fprintf(file," %d",useQ);
  fprintf(file," %d",useR);
  fprintf(file," %d",QHM);
  fprintf(file," %d",xtrSize1);
  fprintf(file," %d",xtrSize2);
  fprintf(file," %d",xtrSize3);
  fprintf(file," %d",runCount);
  fprintf(file," %1.6e",deltaT);
  fprintf(file," %s",fftPlanDescriptor);
  
  fprintf(file,"\n");
  fclose(file);
}


bool TuningDataBase::matchHostNames(char* hn1, char* hn2) {
  bool mismatch = false;
  int hnlen = strlen(hn1);
  if (((int)strlen(hn2))<hnlen) hnlen = strlen(hn2);
  if (hnlen>4) hnlen = 4;
  for (int I=0; I<hnlen; I++) if (hn1[I]!=hn2[I]) mismatch = true;
  return !mismatch;
}


bool TuningDataBase::isDBEntryAvail(char* hostName, int procCount, int L0, int L1, int L2, int L3, int paraOpMode, int threadPerNodeCount, int xFFT, int useP, int useQ, int useR, int QHM, int xtrSize1, int xtrSize2, int xtrSize3) {
  FILE* file = fopen(DBFileName, "r");
  if (file == NULL) return false;
  bool avail = false;
  DBEntryType entry;

  while (readNextEntry(file, entry)) {
    bool mismatch = false;
    if (!matchHostNames(hostName, entry.hostName)) mismatch = true;
    if (procCount != entry.procCount) mismatch = true;
    if (L0 != entry.L0) mismatch = true;
    if (L1 != entry.L1) mismatch = true;
    if (L2 != entry.L2) mismatch = true;
    if (L3 != entry.L3) mismatch = true;
    if (paraOpMode != entry.paraOpMode) mismatch = true;
    if (threadPerNodeCount != entry.threadPerNodeCount) mismatch = true;
    if (xFFT != entry.xFFT) mismatch = true;
    if ((useP != entry.useP) && (!QHM)) mismatch = true;
    if ((useQ != entry.useQ) && (!QHM)) mismatch = true;
    if ((useR != entry.useR) && (QHM)) mismatch = true;
    if (QHM != entry.QHM) mismatch = true;
    if (xtrSize1 != entry.xtrSize1) mismatch = true;
    if (xtrSize2 != entry.xtrSize2) mismatch = true;
    if (xtrSize3 != entry.xtrSize3) mismatch = true;
  
    freeDBEntry(entry);
    if (!mismatch) {
      avail = true;
      break;
    }
  }
    
  if (!avail) freeDBEntry(entry);
  fclose(file);
  return avail;
}


bool TuningDataBase::queryFastestEmbedding(char* hostName, int procCount, int L0, int L1, int L2, int L3, int paraOpMode, int threadPerNodeCount, int xFFT, int useP, int useQ, int useR, int QHM, int& xtrSize1, int& xtrSize2, int& xtrSize3) {
  xtrSize1 = 0;
  xtrSize2 = 0;
  xtrSize3 = 0;
  double bestTiming = 0;

  FILE* file = fopen(DBFileName, "r");
  if (file == NULL) return false;
  bool avail = false;
  DBEntryType entry;
  
  while (readNextEntry(file, entry)) {
    bool mismatch = false;
    if (!matchHostNames(hostName, entry.hostName)) mismatch = true;
//    if (procCount != entry.procCount) mismatch = true;
    if (L0 != entry.L0) mismatch = true;
    if (L1 != entry.L1) mismatch = true;
    if (L2 != entry.L2) mismatch = true;
    if (L3 != entry.L3) mismatch = true;
    if (paraOpMode != entry.paraOpMode) mismatch = true;
    if (threadPerNodeCount != entry.threadPerNodeCount) mismatch = true;
    if (xFFT != entry.xFFT) mismatch = true;
    if ((useP != entry.useP) && (!QHM)) mismatch = true;
    if ((useQ != entry.useQ) && (!QHM)) mismatch = true;
    if ((useR != entry.useR) && (QHM)) mismatch = true;
    if (QHM != entry.QHM) mismatch = true;
    
    if (!mismatch) {
      double timing = entry.runCount / entry.deltaT;
      if ((!avail) || (timing > bestTiming)) {
        xtrSize1 = entry.xtrSize1;
        xtrSize2 = entry.xtrSize2;
        xtrSize3 = entry.xtrSize3;
	bestTiming = timing;
      }
    
      avail = true;
    }
  
    freeDBEntry(entry);
  }
    
  freeDBEntry(entry);

  fclose(file);
  return avail;
}


bool TuningDataBase::queryFastestFFTPlanDescriptor(char* hostName, int procCount, int L0, int L1, int L2, int L3, int paraOpMode, int threadPerNodeCount, int xFFT, int useP, int useQ, int useR, int QHM, int xtrSize1, int xtrSize2, int xtrSize3, char* &fftPlanDescriptor) {
  fftPlanDescriptor = new char[10000];
  snprintf(fftPlanDescriptor, 10000, " ");
  double bestTiming = 0;

  FILE* file = fopen(DBFileName, "r");
  if (file == NULL) return false;
  bool avail = false;
  DBEntryType entry;
  
  while (readNextEntry(file, entry)) {
    bool mismatch = false;
    if (!matchHostNames(hostName, entry.hostName)) mismatch = true;
//    if (procCount != entry.procCount) mismatch = true;
    if (L0 != entry.L0) mismatch = true;
    if (L1 != entry.L1) mismatch = true;
    if (L2 != entry.L2) mismatch = true;
    if (L3 != entry.L3) mismatch = true;
    if (paraOpMode != entry.paraOpMode) mismatch = true;
    if (threadPerNodeCount != entry.threadPerNodeCount) mismatch = true;
    if (xFFT != entry.xFFT) mismatch = true;
    if ((useP != entry.useP) && (!QHM)) mismatch = true;
    if ((useQ != entry.useQ) && (!QHM)) mismatch = true;
    if ((useR != entry.useR) && (QHM)) mismatch = true;
    if (QHM != entry.QHM) mismatch = true;
    if (xtrSize1 != entry.xtrSize1) mismatch = true;
    if (xtrSize2 != entry.xtrSize2) mismatch = true;
    if (xtrSize3 != entry.xtrSize3) mismatch = true;
    
    if (!mismatch) {
      double timing = entry.runCount / entry.deltaT;
      if ((!avail) || (timing > bestTiming)) {
        snprintf(fftPlanDescriptor, 10000, "%s", entry.fftPlanDescriptor);
	bestTiming = timing;
      }
    
      avail = true;
    }
  
    freeDBEntry(entry);
  }
    
  freeDBEntry(entry);

  fclose(file);
  return avail;
}
