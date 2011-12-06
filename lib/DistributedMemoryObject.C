#include "DistributedMemoryObject.h"

DistributedMemoryObject::DistributedMemoryObject(int datStCnt) {
  dataSetCount = datStCnt;
  dataSets = new Complex*[dataSetCount];
  for (int I=0; I<dataSetCount; I++) {
    dataSets[I] = NULL;
  }
}


DistributedMemoryObject::~DistributedMemoryObject() {	
  for (int I=0; I<dataSetCount; I++) {
    destroySuperAlignedComplex(dataSets[I]);
  }
  delete[] dataSets;
}
