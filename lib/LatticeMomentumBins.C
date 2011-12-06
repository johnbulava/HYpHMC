#include "LatticeMomentumBins.h"
#include "Tools.h"

void LatticeMomentumBins::ini() {
  if (initialized) return;
  sinPSqr = new double[L0*L1*L2*L3];
  MomentumSqrSlotLocations = new double[LatticeMomentumBinsDataMAX];
  LatticeMomentumToSlotPointer = new int[L0*L1*L2*L3];
  MomentumSqrSlotCount = 0;
  LatticeNegCoorPointer = new int[L0*L1*L2*L3];   
   
  int I0,I1,I2,I3;
  vector4D p;
  int count = 0;
  for (I0=0; I0<L0; I0++) {
    p[0] = 2*I0*pi/L0;
    int n0 = (L0-I0) % L0;
    for (I1=0; I1<L1; I1++) {
      p[1] = 2*I1*pi/L1;
      int n1 = (L1-I1) % L1;      
      for (I2=0; I2<L2; I2++) {
        p[2] = 2*I2*pi/L2;
        int n2 = (L2-I2) % L2;
        for (I3=0; I3<L3; I3++) {
          p[3] = 2*I3*pi/L3;
          int n3 = (L3-I3) % L3;
	  
	  int nPos = n3 + n2*L3 + n1*L3*L2 + n0*L3*L2*L1;
	  LatticeNegCoorPointer[count] = nPos;
	  sinPSqr[count] = 4.0 * (sqr(sin(0.5*p[0])) + sqr(sin(0.5*p[1])) + sqr(sin(0.5*p[2])) + sqr(sin(0.5*p[3])));
	  
	  int slotNr = findMomentumSqrSlot(sinPSqr[count]);
	  LatticeMomentumToSlotPointer[count] = slotNr;
	  if (slotNr<0) {
	    MomentumSqrSlotLocations[MomentumSqrSlotCount] = sinPSqr[count];
	    LatticeMomentumToSlotPointer[count] = MomentumSqrSlotCount;
	    MomentumSqrSlotCount++;
	  }
	  count++;
	}
      }
    }
  }
  dataSum = new double[MomentumSqrSlotCount];
  dataSqrSum = new double[MomentumSqrSlotCount];
  dataCount = new int[MomentumSqrSlotCount];
  dataAvg = new double[MomentumSqrSlotCount]; 
  dataSigma = new double[MomentumSqrSlotCount]; 
  dataAvgSigma = new double[MomentumSqrSlotCount]; 
  initialized = true;  
  clearData();
}


void LatticeMomentumBins::desini() {
  if (initialized) {
    delete[] sinPSqr;
    delete[] MomentumSqrSlotLocations;
    delete[] LatticeMomentumToSlotPointer;
    delete[] dataSum;
    delete[] dataSqrSum;
    delete[] dataCount;
    delete[] dataAvg;
    delete[] dataSigma;
    delete[] dataAvgSigma;
    delete[] LatticeNegCoorPointer;
    delete[] MomentumSlotMultiplicity;
  }
}


LatticeMomentumBins::LatticeMomentumBins(int l0, int l1, int l2, int l3) {
  L0 = l0;
  L1 = l1;
  L2 = l2;
  L3 = l3;  
  
  sinPSqr = NULL;
  MomentumSqrSlotLocations = NULL;
  LatticeMomentumToSlotPointer = NULL;
  dataSum = NULL;
  dataSqrSum = NULL;
  dataCount = NULL;
  dataAvg = NULL;
  dataSigma = NULL;
  dataAvgSigma = NULL;
  LatticeNegCoorPointer = NULL;
  MomentumSlotMultiplicity = NULL;
  
  initialized = false;
}


LatticeMomentumBins::~LatticeMomentumBins() {
  desini();
}


int LatticeMomentumBins::findMomentumSqrSlot(double momSqr) {
  int I;
  for(I=0; I<MomentumSqrSlotCount; I++) {
    if (abs(MomentumSqrSlotLocations[I]-momSqr) < LatticeMomentumBinsMomentumSqrSlotSize) {
      return I;    
    }  
  }
  return -1;
}


void LatticeMomentumBins::calcAverageAndSigmaVectors() {
  for (int I=0; I<MomentumSqrSlotCount; I++) {
    if (dataCount[I]>0) {
      dataAvg[I] = dataSum[I] / dataCount[I];
      dataSigma[I] = sqrt((dataSqrSum[I]/dataCount[I] - sqr(dataAvg[I])));
      dataAvgSigma[I] = sqrt((dataSqrSum[I]/dataCount[I] - sqr(dataAvg[I])) / independentDataSetCount);
    } else {
      dataAvg[I] = 0;
      dataSigma[I] = 0;
      dataAvgSigma[I] = 0;
    }
  }  
}


void LatticeMomentumBins::clearData() {
  if (!initialized) return;
  for (int I=0; I<MomentumSqrSlotCount; I++) {
    dataSum[I] = 0;
    dataSqrSum[I] = 0;    
    dataCount[I] = 0;
    dataAvg[I] = 0;    
    dataSigma[I] = 0;    
    dataAvgSigma[I] = 0;    
  }
  independentDataSetCount = 0;
}


void LatticeMomentumBins::addDataVector(double* data) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  if (!initialized) ini();
  for (int I=0; I<L0*L1*L2*L3; I++) {
    int index = LatticeMomentumToSlotPointer[I];
    dataSum[index] += data[I];
    dataSqrSum[index] += sqr(data[I]);    
    dataCount[index]++;
  }
  independentDataSetCount++;
  addPerformanceProfilingItem("LatticeMomentumBins::addDataVector", performanceProfilerStartCycle, 0);
}


double* LatticeMomentumBins::getAverageVector() {
  if (!initialized) ini();
  calcAverageAndSigmaVectors();  
  return dataAvg;
}


double* LatticeMomentumBins::getSigmaVector() {
  if (!initialized) ini();
  calcAverageAndSigmaVectors();
  return dataSigma;
}


double* LatticeMomentumBins::getAvgSigmaVector() {
  if (!initialized) ini();
  calcAverageAndSigmaVectors();
  return dataAvgSigma;
}


void LatticeMomentumBins::getAverageVectorInflated(double* avg) {
  if (!initialized) ini();
  calcAverageAndSigmaVectors();
  for (int I=0; I<L0*L1*L2*L3; I++) {
    int index = LatticeMomentumToSlotPointer[I];
    avg[I] = dataAvg[index];
  }
}


void LatticeMomentumBins::getSigmaVectorInflated(double* sig) {
  if (!initialized) ini();
  calcAverageAndSigmaVectors();
  for (int I=0; I<L0*L1*L2*L3; I++) {
    int index = LatticeMomentumToSlotPointer[I];
    sig[I] = dataSigma[index];
  }
}


void LatticeMomentumBins::getAvgSigmaVectorInflated(double* avgsig) {
  if (!initialized) ini();
  calcAverageAndSigmaVectors();
  for (int I=0; I<L0*L1*L2*L3; I++) {
    int index = LatticeMomentumToSlotPointer[I];
    avgsig[I] = dataAvgSigma[index];
  }
}


void LatticeMomentumBins::saveData(char* fileName) {
  if (!initialized) ini();
  calcAverageAndSigmaVectors();
  if (LogLevel>4) printf("Saving LatticeMomentumBins to file: %s...",fileName);
  FILE* file = fopen(fileName, "w");
  int I;
  for (I=0; I<MomentumSqrSlotCount; I++) {
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %d %d\n",MomentumSqrSlotLocations[I], dataAvg[I], dataAvgSigma[I], dataSigma[I], dataSum[I], dataSqrSum[I], dataCount[I], independentDataSetCount);  
  }
  fclose(file);
  if (LogLevel>4) printf("ready.\n");  
}


void LatticeMomentumBins::loadData(char* fileName) {
  if (LogLevel>2) printf("Loading LatticeMomentumBins from file: %s...",fileName);
  FILE* file = fopen(fileName, "r");
  if (file==NULL) {
    if (LogLevel>2) printf("Data File (%s) not found.\n", fileName);
    return;
  }

  if (!initialized) ini();
  clearData();
  
  for (int I=0; I<MomentumSqrSlotCount; I++) {
    fscanf(file,"%lf %lf %lf %lf %lf %lf %d %d\n",&(MomentumSqrSlotLocations[I]), &(dataAvg[I]), &(dataAvgSigma[I]), &(dataSigma[I]), &(dataSum[I]), &(dataSqrSum[I]), &(dataCount[I]), &independentDataSetCount);  
  }
  fclose(file);
  if (LogLevel>2) printf("ready.\n");  
}


void LatticeMomentumBins::addDataVectorFromfourierTrafoSPECIAL(Complex* data) {
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  if (!initialized) ini();
  for (int I=0; I<L0*L1*L2*L3; I++) {
    int index = LatticeMomentumToSlotPointer[I];
    int nI = LatticeNegCoorPointer[I];
    Complex c1, c2;
    int I2;
    for (I2=0; I2<2; I2++) {
      c1.x = 0.5*(data[2*I+I2].x + data[2*nI+I2].x);
      c1.y = 0.5*(data[2*I+I2].y - data[2*nI+I2].y);
      c2.x = 0.5*(data[2*I+I2].x - data[2*nI+I2].x);
      c2.y = 0.5*(data[2*I+I2].y + data[2*nI+I2].y);
      
      double sc1 = c1.x*c1.x + c1.y*c1.y;
      double sc2 = c2.x*c2.x + c2.y*c2.y;
      
      dataSum[index] += sqrt(sc1);
      dataSum[index] += sqrt(sc2);
      dataSqrSum[index] += sc1 + sc2;
      dataCount[index] += 2;
    }
  }
  independentDataSetCount++;
  addPerformanceProfilingItem("LatticeMomentumBins::addDataVectorFromfourierTrafoSPECIAL", performanceProfilerStartCycle, 0);
}


void LatticeMomentumBins::addDataVectorFromOmegafourierTrafoSPECIAL(Complex* data) {  
  long int performanceProfilerStartCycle = getPerformanceProfilingStartCycle(0);
  if (!initialized) ini();
  double normFac = 1.0;
  int count = 0;
  int I3, I2, I1, I0, i;
  int p = 0;
  int d0 = 8*(L3 + xtraSize3)*(L2 + xtraSize2)*(L1 + xtraSize1) - L1*8*(L3 + xtraSize3)*(L2 + xtraSize2);
  int d1 = 8*(L3 + xtraSize3)*(L2 + xtraSize2) - L2*8*(L3 + xtraSize3);
  int d2 = 8*(L3 + xtraSize3) - L3*8;
  int d3 = 8;
  for (I0=0; I0<L0; I0++) {
    for (I1=0; I1<L1; I1++) {
      for (I2=0; I2<L2; I2++) {
        for (I3=0; I3<L3; I3++) {
          int index = LatticeMomentumToSlotPointer[count];
	  for (i=0; i<8; i++) {
            dataSum[index] += normFac*sqrt(sqr(data[p+i].x) + sqr(data[p+i].y));	    
            dataSqrSum[index] += sqr(normFac*sqrt(sqr(data[p+i].x) + sqr(data[p+i].y)));
            dataCount[index] += 1;
	  }
	  count++;
          p += d3;
        }
        p += d2;
      }
      p += d1;
    }
    p += d0;
  }  
  independentDataSetCount++;
  addPerformanceProfilingItem("LatticeMomentumBins::addDataVectorFromOmegafourierTrafoSPECIAL", performanceProfilerStartCycle, 0);
}


double LatticeMomentumBins::getLatMomSqrFromIndex(int index) {
  if (!initialized) ini();
  return sinPSqr[index];
}


bool LatticeMomentumBins::containsData() {
  if (!initialized) return false;
  return (independentDataSetCount>0);
}


int LatticeMomentumBins::getMomentumSqrSlotCount() {
  if (!initialized) ini();
  return MomentumSqrSlotCount;
}


double LatticeMomentumBins::getLatMomSqrFromSlotNr(int slotNr) {
  if (!initialized) ini();
  if (slotNr<0) return NaN;
  if (slotNr>=MomentumSqrSlotCount) return NaN;
  return MomentumSqrSlotLocations[slotNr];
}


int LatticeMomentumBins::getMomentumSlotFromIndex(int index) {
  if (!initialized) ini();  
  return LatticeMomentumToSlotPointer[index];
}


int LatticeMomentumBins::getMomentumSlotMultiplicity(int slotNr) {
  if (!initialized) ini();
  if (slotNr<0) return 0;
  if (slotNr>=MomentumSqrSlotCount) return 0;
  
  if (MomentumSlotMultiplicity==NULL) {
    MomentumSlotMultiplicity = new int[MomentumSqrSlotCount];
    for (int I=0; I<MomentumSqrSlotCount; I++) {
      MomentumSlotMultiplicity[I] = 0;
    }
    for (int I=0; I<L0*L1*L2*L3; I++) {
      MomentumSlotMultiplicity[LatticeMomentumToSlotPointer[I]]++;
    }  
  }
  return MomentumSlotMultiplicity[slotNr];
}
