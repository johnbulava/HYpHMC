#include "LatticeSiteBins.h"
#include <stdio.h>
#include "Tools.h"

void LatticeSiteBins::ini(int l0, int l1, int l2, int l3) {
  L0 = l0;
  L1 = l1;
  L2 = l2;
  L3 = l3;  

  sinPSqr = new double[L0*L1*L2*L3];
   
  int I0,I1,I2,I3;
  vector4D p;
  int count = 0;
  for (I0=0; I0<L0; I0++) {
    p[0] = 2*I0*pi/L0;
    for (I1=0; I1<L1; I1++) {
      p[1] = 2*I1*pi/L1;
      for (I2=0; I2<L2; I2++) {
        p[2] = 2*I2*pi/L2;
        for (I3=0; I3<L3; I3++) {
          p[3] = 2*I3*pi/L3;
	  sinPSqr[count] = 4.0 * (sqr(sin(0.5*p[0])) + sqr(sin(0.5*p[1])) + sqr(sin(0.5*p[2])) + sqr(sin(0.5*p[3])));
	  count++;
	}
      }
    }
  }
  dataSum = new double[L0*L1*L2*L3];
  dataSqrSum = new double[L0*L1*L2*L3];
  dataCount = 0;
  independentDataSetCount = 0;
  clearData();
}


void LatticeSiteBins::desini() {
  delete[] sinPSqr;
  delete[] dataSum;
  delete[] dataSqrSum;
}


LatticeSiteBins::LatticeSiteBins(int l0, int l1, int l2, int l3) {
  ini(l0,l1,l2,l3);
}


LatticeSiteBins::~LatticeSiteBins() {
  desini();
}


void LatticeSiteBins::clearData() {
  int I;
  for (I=0; I<L0*L1*L2*L3; I++) {
    dataSum[I] = 0;
    dataSqrSum[I] = 0;    
  }
  independentDataSetCount = 0;
  dataCount = 0;
}


void LatticeSiteBins::getAverageVector(double* avg) {
  int I;
  for(I=0; I<L0*L1*L2*L3; I++) {
    if (dataCount>0) {
      avg[I] = dataSum[I] / dataCount;
    } else {
      avg[I] = 0;
    }
  } 
}


void LatticeSiteBins::getSigmaVector(double* sig) {
  int I;
  for(I=0; I<L0*L1*L2*L3; I++) {
    if (dataCount>0) {
      double avg = dataSum[I] / dataCount;
      sig[I] = sqrt((dataSqrSum[I]/dataCount - sqr(avg)));
    } else {
      sig[I] = 0;
    }
  } 
}


void LatticeSiteBins::getAvgSigmaVector(double* avgsig) {
  int I;
  for(I=0; I<L0*L1*L2*L3; I++) {
    if (dataCount>0) {
      double avg = dataSum[I] / dataCount;
      avgsig[I] = sqrt((dataSqrSum[I]/dataCount - sqr(avg)) / independentDataSetCount);
    } else {
      avgsig[I] = 0;
    }
  } 
}


void LatticeSiteBins::saveData(char* fileName) {
  if (LogLevel>4) printf("Saving LatticeSiteBins to file: %s...",fileName);
  FILE* file = fopen(fileName, "w");
  int I;
  for (I=0; I<L0*L1*L2*L3; I++) {
    double avg = 0;
    double sig = 0;
    double avgsig = 0;
    if (dataCount>0) {
      avg = dataSum[I] / dataCount;
      sig = sqrt((dataSqrSum[I]/dataCount - sqr(avg)));
      avgsig = sqrt((dataSqrSum[I]/dataCount - sqr(avg)) / independentDataSetCount);
    }
  
    fprintf(file,"%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %d %d\n", sinPSqr[I], avg, sig, avgsig, dataSum[I], dataSqrSum[I], dataCount, independentDataSetCount);  
  }
  fclose(file);
  if (LogLevel>4) printf("ready.\n");  
}


void LatticeSiteBins::loadData(char* fileName) {
  if (LogLevel>2) printf("Loading LatticeSiteBins from file: %s...",fileName);
  clearData();
  
  FILE* file = fopen(fileName, "r");
  if (file==NULL) {
    if (LogLevel>2) printf("Data File (%s) not found.\n", fileName);
    return;
  }
  int I;
  double d1,d2,d3,d4;
  bool error = false;
  for (I=0; I<L0*L1*L2*L3; I++) {
    if (fscanf(file,"%lf %lf %lf %lf %lf %lf %d %d\n",&d1, &d2, &d3, &d4, &(dataSum[I]), &(dataSqrSum[I]), &dataCount, &independentDataSetCount) != 8) {
      error = true;    
    }  
  }
  fclose(file);
  if (error) {
    clearData();
    if (LogLevel>2) printf("Error. Data Cleared.\n");      
  } else {
    if (LogLevel>2) printf("ready.\n");  
  }
}


void LatticeSiteBins::addDataVectorFromOmegafourierTrafo(Complex* data) {
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
	  for (i=0; i<8; i++) {
            dataSum[count] += normFac*sqrt(sqr(data[p+i].x) + sqr(data[p+i].y));
            dataSqrSum[count] += sqr(normFac*sqrt(sqr(data[p+i].x) + sqr(data[p+i].y)));
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
  dataCount+=8;
  independentDataSetCount++;
}


double LatticeSiteBins::getLatMomSqrFromIndex(int index) {
  return sinPSqr[index];
}


bool LatticeSiteBins::containsData() {
  return (independentDataSetCount>0);
}
