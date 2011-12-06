#include "CacheSimulator.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Global.h"

CacheSimulator::CacheSimulator(int cacheS, int ways, int cachLineS) { 
  CacheSize = cacheS;
  CacheWays = ways;
  CacheLineSize = cachLineS;
  CacheLinesPerWay = CacheSize / (CacheWays*CacheLineSize);
  
  if (LogLevel>4) printf("CacheSimulator initializing with CacheSize: %d, CacheWays: %d, CacheLineSize: %d, CacheLinesPerWay: %d\n",CacheSize,CacheWays,CacheLineSize,CacheLinesPerWay);
  
  //Consistency-Checks
  if (CacheLinesPerWay*CacheWays*CacheLineSize != CacheSize) {
    printf("ERROR in CacheSimulator: CacheSize not dividable by CacheWays*cachLineS!!!\n");
    exit(0);
  }
  
  CachedAddresses = new long int*[CacheLinesPerWay];
  for (int I=0; I<CacheLinesPerWay; I++) {
    CachedAddresses[I] = new long int[2*CacheWays];
  }
  LastAccessed = new long int[CacheLinesPerWay];
  
  reset();
}


CacheSimulator::~CacheSimulator() { 
  delete[] LastAccessed;
  for (int I=0; I<CacheLinesPerWay; I++) {
    delete[] CachedAddresses[I];
  }
  delete[] CachedAddresses;
}


void CacheSimulator::accessAddress(long int addr) {
  if (addr<0) return;
  accesses++;
  
  long int p =  addr % (CacheLinesPerWay*CacheLineSize);
  p /= CacheLineSize;
  addr /= CacheLineSize;
  addr *= CacheLineSize;
  
  if (LastAccessed[p] != addr) lastAccessMisses++;
  LastAccessed[p] = addr;
  
  //Check whether datum already in Cache
  for (int I=0; I<CacheWays; I++) {
    if (CachedAddresses[p][2*I+1] == addr) {
      CachedAddresses[p][2*I+0] = accesses;
      return;
    }
  }
  misses++;
  
  //Kick out oldest datum
  int oldestTime = CachedAddresses[p][0];
  for (int I=0; I<CacheWays; I++) {
    if (CachedAddresses[p][2*I+0] < oldestTime) oldestTime = CachedAddresses[p][2*I+0];
  }
  for (int I=0; I<CacheWays; I++) {
    if (CachedAddresses[p][2*I+0] == oldestTime) {
      CachedAddresses[p][2*I+0] = accesses;
      CachedAddresses[p][2*I+1] = addr;
      return;
    }
  }
}


void CacheSimulator::reset() {
  accesses = 0;
  misses = 0;
  lastAccessMisses = 0;
  for (int I=0; I<CacheLinesPerWay; I++) {
    LastAccessed[I] = -1;
    for (int I2=0; I2<2*CacheWays; I2++) CachedAddresses[I][I2] = -1;
  }
}


int CacheSimulator::getMisses() {
  return misses;
}


int CacheSimulator::getLastAccessMisses() {
  return lastAccessMisses;
}


int CacheSimulator::getAccesses() {
  return accesses;
}
