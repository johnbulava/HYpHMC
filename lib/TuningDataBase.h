#ifndef TuningDataBase_included
#define TuningDataBase_included

#include <math.h>
#include <stdlib.h>
#include <stdio.h>





class TuningDataBase {
private:
  struct DBEntryType {
    char* hostName;
    int procCount;
    int L0;
    int L1;
    int L2;
    int L3;
    int paraOpMode;
    int threadPerNodeCount;
    int xFFT;
    int useP; 
    int useQ; 
    int useR;
    int QHM;
    int xtrSize1;
    int xtrSize2;
    int xtrSize3;
    int runCount;
    double deltaT;
    char* fftPlanDescriptor;
  };

  char* DBFileName;
  
  
  void freeDBEntry(DBEntryType &entry);
  bool readNextEntry(FILE* file, DBEntryType &entry);
  bool matchHostNames(char* hn1, char* hn2);

public:
  TuningDataBase(char* dbFname); 
  ~TuningDataBase();

  void writeDBEntry(char* hostName, int procCount, int L0, int L1, int L2, int L3, int paraOpMode, int threadPerNodeCount, int xFFT, int useP, int useQ, int useR, int QHM, int xtrSize1, int xtrSize2, int xtrSize3, int runCount, double deltaT, char* fftPlanDescriptor);
  bool isDBEntryAvail(char* hostName, int procCount, int L0, int L1, int L2, int L3, int paraOpMode, int threadPerNodeCount, int xFFT, int useP, int useQ, int useR, int QHM, int xtrSize1, int xtrSize2, int xtrSize3);
  bool queryFastestEmbedding(char* hostName, int procCount, int L0, int L1, int L2, int L3, int paraOpMode, int threadPerNodeCount, int xFFT, int useP, int useQ, int useR, int QHM, int& xtrSize1, int& xtrSize2, int& xtrSize3);
  bool queryFastestFFTPlanDescriptor(char* hostName, int procCount, int L0, int L1, int L2, int L3, int paraOpMode, int threadPerNodeCount, int xFFT, int useP, int useQ, int useR, int QHM, int xtrSize1, int xtrSize2, int xtrSize3, char* &fftPlanDescriptor);
};


#endif
