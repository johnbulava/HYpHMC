#include "PerformanceProfiler.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Global.h"

PerformanceProfiler::PerformanceProfiler(const char* fileN) { 
  fileName = new char[2000];
  snprintf(fileName, 2000, "%s", fileN);
  node0StartedTimerCount = 0;
  performanceItems = new PerformanceItemType*[PerformanceProfiler_ItemMAX];
  for (int I=0; I<PerformanceProfiler_ItemMAX; I++) {
    performanceItems[I] = new PerformanceItemType[PerformanceProfiler_NodeMAX];
    for (int I2=0; I2<PerformanceProfiler_NodeMAX; I2++) {
      performanceItems[I][I2].stringLength = 0;
      performanceItems[I][I2].routineName = NULL;
      performanceItems[I][I2].callCount = 0;
      performanceItems[I][I2].cycleSum = 0;
    }
  }
  performanceItemCount = 0;
  StartCycle = getCPUCycleCounter();
  coveredCyclesByNode0 = 0;
  pthread_mutex_init (&pThreadMutex, NULL);
  if (LogLevel>2) printf("PerformanceProfiler initialized with fileName = %s\n", fileName);
}


PerformanceProfiler::~PerformanceProfiler() { 
  delete[] fileName;
  for (int I=0; I<PerformanceProfiler_ItemMAX; I++) {
    if (performanceItems[I] != NULL) {
      for (int I2=0; I2<PerformanceProfiler_NodeMAX; I2++) {
        delete[] performanceItems[I][I2].routineName;
      }
      delete[] performanceItems[I];
    }
  }
  delete[] performanceItems;
}


void PerformanceProfiler::addPerformanceItem(const char* rName, long int startCycle, int node) {
  long int cycles = getCPUCycleCounter() - startCycle;
  if (cycles<0) cycles = 0;
  
  if (node==0) {
    node0StartedTimerCount--;
    if (node0StartedTimerCount == 0) {
      coveredCyclesByNode0 += cycles;
    }
  }
  if (node0StartedTimerCount<0) {
    printf("ERROR in PerformanceProfiler::addPerformanceItem: node0StartedTimerCount<0 \n");
    exit(0);
  }
  if ((node<0) || (node>=PerformanceProfiler_NodeMAX)) {
    printf("ERROR in PerformanceProfiler::addPerformanceItem: invalid node number %d\n", node);
    exit(0);  
  }
  
  int stringLength = strlen(rName);

  for (int I=0; I<performanceItemCount; I++) {
    if (performanceItems[I][0].stringLength == stringLength) {
      bool match = true;
      for (int I2=0; I2<stringLength; I2++) {
        if (rName[I2] != performanceItems[I][0].routineName[I2]) {
	  match = false;
	  break;
	}      
      }
      if (match) {      
        addCyclesToItem(performanceItems[I][node], cycles);
	return;
      }
    }
  }

  pthread_mutex_lock(&pThreadMutex);
  int pC = performanceItemCount;
  for (int I=0; I<pC; I++) {
    if (performanceItems[I][0].stringLength == stringLength) {
      bool match = true;
      for (int I2=0; I2<stringLength; I2++) {
        if (rName[I2] != performanceItems[I][0].routineName[I2]) {
	  match = false;
	  break;
	}      
      }
      if (match) {
        addCyclesToItem(performanceItems[I][node], cycles);
        pthread_mutex_unlock(&pThreadMutex);	
	return;
      }
    }
  }

  for (int I=0; I<PerformanceProfiler_NodeMAX; I++) {
    performanceItems[pC][I].stringLength = stringLength;
    performanceItems[pC][I].routineName = new char[stringLength+10];
    snprintf(performanceItems[pC][I].routineName, stringLength+10, "%s", rName);
    performanceItems[pC][I].callCount = 0;
    performanceItems[pC][I].cycleSum = 0;
  }
  addCyclesToItem(performanceItems[pC][node], cycles);
  performanceItemCount = pC+1;
  pthread_mutex_unlock(&pThreadMutex);
}


void PerformanceProfiler::writePerformanceItemsToDisk() {
  double totalCycles = getCPUCycleCounter() - StartCycle;
  FILE* file = fopen(fileName, "w");
  int alignPOS = 80;
  
  fprintf(file, "Total GCycles: %1.2f \n", totalCycles/1E9);  
  fprintf(file, "Percentage of covered cycles: %1.2f\n", 100*coveredCyclesByNode0 / totalCycles);
  fprintf(file, "CPU-GCycles per second: %1.2f\n", CPUCyclesPerSecond/1E9);
  fprintf(file, "Number of registered items: %d \n\n", performanceItemCount);
    
  for (int I=0; I<performanceItemCount; I++) {
    for (int I2=0; I2<PerformanceProfiler_NodeMAX; I2++) {
      if (performanceItems[I][I2].callCount>0) {
        fprintf(file, "%s ",performanceItems[I][I2].routineName);
	for (int I3=0; I3<(int)(alignPOS-strlen(performanceItems[I][I2].routineName)); I3++) {
          fprintf(file, " ");	
	}
      
        fprintf(file, "(Node %d):    %1.2f Percentage, Calls: %d, Total GCycles: %1.2f, Average MCycles %1.2f\n",I2,100*performanceItems[I][I2].cycleSum/totalCycles,performanceItems[I][I2].callCount,performanceItems[I][I2].cycleSum/1E9,performanceItems[I][I2].cycleSum/(performanceItems[I][I2].callCount*1E6));
      }
    }
  }
  
  fclose(file);
}
