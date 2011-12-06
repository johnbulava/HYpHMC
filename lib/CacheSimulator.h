#ifndef CacheSimulator_included
#define CacheSimulator_included


class CacheSimulator {
private:
  int CacheSize;      //in Bytes
  int CacheWays;
  int CacheLineSize;  //in Bytes
  int CacheLinesPerWay;
  
  long int** CachedAddresses;
  long int* LastAccessed;
  
  int accesses;
  int misses;
  int lastAccessMisses;

public:
  CacheSimulator(int cacheS, int ways, int cachLineS); 
  ~CacheSimulator();
  
  void accessAddress(long int addr);
  void reset();
  int getMisses();
  int getLastAccessMisses();
  int getAccesses();


};


#endif
