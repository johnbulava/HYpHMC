#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>


#define dim 4
#define regs 2


int** vecs;
int Max;
int* actBranch;
int* bestBranch;
int bestResult;


int pow(int base, int ex) {
  int res = 1;
  int I;
  for (I=1; I<=ex; I++) {
    res *= base;
  }
  return res;
}


void iniVecs() {
  int I,I2,k;
  vecs = new int*[Max+1];
  for (I=0; I<=Max; I++) {
    vecs[I] = new int[Max+1];
    for (I2=0; I2<=Max; I2++) {
      vecs[I][I2] = 0;
    }
  }
  for (I=1; I<=Max; I++) {
    for (I2=1; I2<=Max; I2++) {
      int m = I;
      int n = I2;
      int s = 0;
      for (k=1; k<=dim; k++) {
        s += (m%2)*(n%2);
	m /= 2;
	n /= 2;
      }
      if ((s%2)==1) {
        vecs[I][I2] = -1; 
      } else {
        vecs[I][I2] = 1;
      }
    }
  }
  actBranch = new int[Max+regs];
  bestBranch = new int[Max+regs];
  for (I=0; I<Max+regs; I++) {
    actBranch[I] = 0;
    bestBranch[I] = 0;
  }
  bestResult = 1000000;
}


void printVecs() {
  int I,I2;
  for (I=1; I<=Max; I++) {
    for (I2=1; I2<=Max; I2++) {
      if (vecs[I][I2]==1) {
        printf(" +");
      } else {
        printf(" -");
      }
    }
    printf("\n");
  }
}


int diffSteps(int start, int target) {
  int I;
  int diff = 0;
  for (I=1; I<=Max; I++) {
    if (vecs[start][I]!=vecs[target][I]) diff++;
  }
  return diff;
}


int diffRegs(int* state, int* target) {
  int I;
  int d;
  int dmax = 0;
  for (I=0; I<regs; I++) {
    d = diffSteps(state[I], target[I]);
    if (d>dmax) dmax = d;
  }
  return dmax;
}


int scanBranches(int depth, int alreadyUsedOps) {
  int I,k,k2,m,I2,I3;
  int M = pow(Max,regs);
  int** ind = new int*[M];
  int indCount = 0;
  bool b;
  int min = 1000000;
  int max = 0;
  
  ind[0] = new int[regs+1];
  for (I=0; I<M; I++) {
    m = I;
    for (k=0; k<regs; k++) {
      ind[indCount][k] = 1 + (m % Max);
      m /= Max;
    }
    b = true;
    
    
    for (k=0; k<regs; k++) {
      if (vecs[ind[indCount][k]][0] == 1) b = false;
      
      for (k2=k+1; k2<regs; k2++) {
        if (ind[indCount][k] == ind[indCount][k2]) b = false;
      }
    }
    if (b) {
      ind[indCount][regs] = diffRegs(ind[indCount], &(actBranch[depth]));
      if (ind[indCount][regs]<min) min = ind[indCount][regs];
      if (ind[indCount][regs]>max) max = ind[indCount][regs];
      indCount++;
      ind[indCount] = new int[regs+1];
    }
  }

//printf("%d %d\n", min, max);
  int res;
  int best = 1000000;
  for (I3=min; I3<=max; I3++) {
    for (I=0; I<indCount; I++) {
      if ((alreadyUsedOps + ind[I][regs]<bestResult) && (ind[I][regs]==I3)){
        for (I2=0; I2<regs; I2++) {
          actBranch[depth+regs+I2] = ind[I][I2];
          vecs[ind[I][I2]][0] = 1;
        }
        res = ind[I][regs] + scanBranches(depth+regs, alreadyUsedOps+ind[I][regs]);
        for (I2=0; I2<regs; I2++) {
          vecs[ind[I][I2]][0] = 0;
        }
        if (res<best) {
          best = res;
        }
      }
    }
  }

  if (indCount==0) {
    if (alreadyUsedOps<bestResult) {
      bestResult = alreadyUsedOps;
      printf("Best found with %d\n",bestResult);
      
      for (I2=0; I2<Max+regs; I2++) {
        bestBranch[I2] = actBranch[I2];
        printf("%d ",bestBranch[I2]);
      }
      printf("\n");
    }
  }
  
  for (I=0; I<=indCount; I++) {
    delete[] ind[I];
  }
  delete[] ind;
  return best;
}


int main(int argc,char **argv) {
  Max = pow(2,dim);
  iniVecs();
  printVecs();
  scanBranches(0,0);
  
  diffSteps(1,2);
}
