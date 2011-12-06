#ifndef DistributedMemoryObject_included
#define DistributedMemoryObject_included

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Global.h"
#include "Complex.h"
#include "Tools.h"


class DistributedMemoryObject {
private:  

public:    
  DistributedMemoryObject(int datStCnt); 
  ~DistributedMemoryObject();  

  Complex** dataSets;
  int dataSetCount;
};

#endif
