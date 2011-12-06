BootStrapClass::BootStrapClass() {
  DataPoolSizeIncrement = 10000;
  DataPoolSize = DataPoolSizeIncrement;
  DataPool = new double[DataPoolSize];
  DataWeightPool = new long double[DataPoolSize];
  DataCount = 0;
  name = new char[100];
  fileName = new char[100];
  setName("BootStrapData");
}


BootStrapClass::~BootStrapClass() {
  delete[] DataPool;
  delete[] DataWeightPool;
  delete[] name;
  delete[] fileName;
}


void BootStrapClass::setName(char* n) {
  snprintf(name,100,"%s",n);
  snprintf(fileName,100,"data/BootStrap/%s.dat",n);
}


int BootStrapClass::getDataCount() {
  return DataCount;
}


void BootStrapClass::saveData() {
  FILE *file;
  int I;

  if (DataCount == 0) return;
  file = fopen(fileName,"w");
  for (I=0; I<DataCount; I++) {  
    fprintf(file,"%1.15f %Lf\n",DataPool[I],DataWeightPool[I]);
  }
  fclose(file);
}

 
void BootStrapClass::loadData(int setSize, int Jstart, int Jend) {
  FILE *file;
  int size;
  double value;
  long double value2;
  int count = 0;

  resetPool();
  
  file = fopen(fileName,"r");
  if (file != NULL) {
    double* dummy = new double;  
    long double* dummy2 = new long double;  
    while (fscanf(file,"%lf %Lf\n",dummy, dummy2)==2) { 
      value = *dummy;
      value2 = *dummy2;
      count++;
    }
    free(dummy);
    free(dummy2);
    fclose(file);
  }
  size = count;
  count = 0;
  
  file = fopen(fileName,"r");
  if (file != NULL) {
    double* dummy = new double;  
    long double* dummy2 = new long double;  
    while (fscanf(file,"%lf %Lf\n",dummy, dummy2)==2) { 
      value = *dummy;
      value2 = *dummy2;
      count++;
      if ((count>setSize*Jstart) && (count<=size-setSize*Jend)) {
        addToDataPool(value, value2);
      }
    }
    free(dummy);
    free(dummy2);
    fclose(file);
  }
  printf("Anzahl geladener Zeilen: %d\n",DataCount);
}


void BootStrapClass::printData() {
  int I;
  printf("Print Data belonging to BootStrap Pool: %s --> %d values\n",name,DataCount);
  for (I=0; I<DataCount; I++) {  
    printf("%1.15f %Lf\n",DataPool[I],DataWeightPool[I]);
  }
  printf("\n");
}
 
  
void BootStrapClass::addToDataPool(double x, long double w) {
  if (DataCount>=DataPoolSize) {
    double* buffer = new double[DataCount];
    long double* buffer2 = new long double[DataCount];
    int I;
    for (I=0; I<DataCount; I++) {
      buffer[I] = DataPool[I];
      buffer2[I] = DataWeightPool[I];
    }
    free(DataPool);
    free(DataWeightPool);
    DataPoolSize += DataPoolSizeIncrement;
    DataPoolSizeIncrement *= 2;
    DataPool = new double[DataPoolSize];
    DataWeightPool = new long double[DataPoolSize];
    for (I=0; I<DataCount; I++) {
      DataPool[I] = buffer[I];
      DataWeightPool[I] = buffer2[I];
    }
    delete[] buffer;
    delete[] buffer2;
  }
  DataPool[DataCount] = x;
  DataWeightPool[DataCount] = w;
  DataCount++;
}


void BootStrapClass::resetPool() {
  delete[] DataPool;
  delete[] DataWeightPool;
  DataPoolSizeIncrement = 10;
  DataPool = new double[DataPoolSizeIncrement];
  DataWeightPool = new long double[DataPoolSizeIncrement];
  DataPoolSize = DataPoolSizeIncrement;
  DataCount = 0;  
}


double BootStrapClass::getAverage() {
  int I;
  long double z1 = 0;
  long double z2 = 0;
  
  if (DataCount == 0) return NaN;
  for (I=0; I<DataCount; I++) {
    z1 += (DataPool[I] * DataWeightPool[I]) / DataCount;  
    z2 += DataWeightPool[I] / DataCount;  
  }
  return (double) (z1/z2);
}


void BootStrapClass::doBootStrapWithWeights(int iter, double& avg, double& sigma, double& SusAvg, double& SusSigma, double& BinderAvg, double& BinderSigma) {
  int I, I2;
  
  if ((DataCount<=0) || (iter<=0)) {
    avg = NaN;
    sigma = NaN;  
    SusAvg = NaN;
    SusSigma = NaN;
    BinderAvg = NaN;
    BinderSigma = NaN;
  }
  
  avg = 0;
  sigma = 0;
  SusAvg = 0;
  SusSigma = 0;
  BinderAvg = 0;
  BinderSigma = 0;
  
  for (I2=0; I2<iter; I2++) {
    long double z1 = 0;
    long double z2 = 0;
    long double z3 = 0;
    long double z4 = 0;
    for (I=0; I<DataCount; I++) {
      int index = (int)  (DataCount*AdvancedZufall(AdvancedSeed));
      if (index<0) index = 0;
      if (index>=DataCount) index = DataCount-1;
      
      z1 += (DataPool[index]*DataWeightPool[index]) / DataCount;
      z2 += DataWeightPool[index] / DataCount;
      z3 += (DataPool[index]*DataPool[index]*DataWeightPool[index]) / DataCount;
      z4 += (DataPool[index]*DataPool[index]*DataPool[index]*DataPool[index]*DataWeightPool[index]) / DataCount;
    }
    avg += (double) (z1/z2);
    sigma += (double) ((z1/z2)*(z1/z2));
    SusAvg += (double) ((z3-z1*z1)/z2);
    SusSigma += (double) (((z3-z1*z1)/z2)*((z3-z1*z1)/z2));
    BinderAvg += (double) (z2*z4/(z3*z3));
    BinderSigma += (double) ((z2*z4/(z3*z3)) * (z2*z4/(z3*z3)));
  }
  
  avg /= iter;
  sigma = sqrt(sigma/iter - avg*avg);
  SusAvg /= iter;
  SusSigma = sqrt(SusSigma/iter - SusAvg*SusAvg);
  BinderAvg /= iter;
  BinderSigma = sqrt(BinderSigma/iter - BinderAvg*BinderAvg);
  BinderAvg = 2 - BinderAvg;
}

