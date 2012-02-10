#ifndef PerformanceProfiler_included
#define PerformanceProfiler_included

#include <math.h>
#include <pthread.h>
//#include "Tools.h"


#define PerformanceProfiler_NodeMAX 10
#define PerformanceProfiler_ItemMAX 300


class PerformanceProfiler {
private:
  char* fileName;
  int node0StartedTimerCount;  
  long int StartCycle;
  long int coveredCyclesByNode0;
  pthread_mutex_t pThreadMutex;
  
  struct PerformanceItemType {
    int stringLength;
    char* routineName;
    int callCount;
    long int cycleSum;
  };
  
  PerformanceItemType** performanceItems;
  int performanceItemCount;
  
  void addCyclesToItem(PerformanceItemType &item, long int cycles);
  long int getCPUCycleCounter();
  
public:
  PerformanceProfiler(const char* fileN); 
  ~PerformanceProfiler();

  long int getTimerStartCPUcycle(int node);
  void addPerformanceItem(const char* rName, long int startCycle, int node);
  void writePerformanceItemsToDisk();
};

inline long int PerformanceProfiler::getCPUCycleCounter() {
  long int val = 0;
  unsigned int __a,__d; 
  asm volatile("rdtsc" : "=a" (__a), "=d" (__d)); 
  val = ((unsigned long)__a) | (((unsigned long)__d)<<32); 
  return val;
} 


inline long int PerformanceProfiler::getTimerStartCPUcycle(int node) {
  if (node==0) {
    node0StartedTimerCount++;
  }
  return getCPUCycleCounter();
}


inline void PerformanceProfiler::addCyclesToItem(PerformanceItemType &item, long int cycles) {
  if (item.callCount == 0) {
    if (cycles<0) {
      item.callCount++;
      item.cycleSum += 0;
    } else {
      item.callCount++;
      item.cycleSum += cycles;  
    }
  } else {
    if ((cycles<0) && (cycles>((10*item.cycleSum)/item.callCount))) {
      item.callCount++;
      item.cycleSum += 0;
    } else {
      item.callCount++;
      item.cycleSum += cycles;  
    }  
  }
}


#endif
