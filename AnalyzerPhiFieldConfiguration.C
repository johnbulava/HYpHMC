AnalyzerPhiFieldConfiguration::AnalyzerPhiFieldConfiguration(char* fileName, FermionMatrixOperations* fOps) {
  fermiOps = fOps;
  
  phiField = (double*) fermiOps->createFermionVector(2);
  phiFieldCopies = new double*[1000];
  phiFieldCopyCount = 0;

  weight = NaN;
  errorState = 0;
  weightAvail = false;
  magnetizationM = NaN;
  magnetizationS = NaN;
  phiFieldAvgVectorLength = NaN;
  phiFieldAvgVectorLengthVariation = NaN;
  avgPhiFieldVector[0] = NaN;
  avgPhiFieldVector[1] = NaN;
  avgPhiFieldVector[2] = NaN;
  avgPhiFieldVector[3] = NaN;
  
  loadPhiconfiguration(fileName);
  measureMagnetizations();
}


AnalyzerPhiFieldConfiguration::~AnalyzerPhiFieldConfiguration() { 
  Complex* c1 = (Complex*) phiField;
  fermiOps->destroyFermionVector(c1);
  phiField = NULL;  
  
  clearPhiFieldCopies();
  delete[] phiFieldCopies;
}


double* AnalyzerPhiFieldConfiguration::getPhiFieldCopy() {
  phiFieldCopies[phiFieldCopyCount] = (double*) fermiOps->createFermionVector(2);

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
   
  SSE_ZCopy(2*L0*L1*L2*L3, (Complex*) phiField, 1, (Complex*) (phiFieldCopies[phiFieldCopyCount]), 1);
  
  phiFieldCopyCount++;
  return phiFieldCopies[phiFieldCopyCount-1];
}


void AnalyzerPhiFieldConfiguration::clearPhiFieldCopies() {
  for (int I=0; I<phiFieldCopyCount; I++) {
    Complex* c = (Complex*) phiFieldCopies[I];
    fermiOps->destroyFermionVector(c);
    phiFieldCopies[I] = NULL;  
  }
  phiFieldCopyCount = 0;
}


int AnalyzerPhiFieldConfiguration::getErrorState() {
  return errorState;
}


bool AnalyzerPhiFieldConfiguration::isWeightAvail() {
  return weightAvail;
}


double AnalyzerPhiFieldConfiguration::getWeight() {
  return weight;
}


void AnalyzerPhiFieldConfiguration::loadPhiconfiguration(char* fileName) {
  if (LogLevel>2) printf("Loading Configuration: %s...",fileName);  

  std::fstream confFile;
  confFile.open(fileName, std::ios::in);
  if (!confFile.good()) {
    if (LogLevel>2) printf("ERROR !!!\n");
    errorState = 1;
    return;
  }

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  confFile.read((char*)phiField, 32*L0*L1*L2*L3);
  if (confFile.eof()) {
    if (LogLevel>2) printf("ERROR !!!\n");
    errorState = 2;    
    return;
  }

  confFile.read((char*)(&weight),8);
  if (!confFile.eof()) {
    if (LogLevel>2) printf(" with weight: %f ",weight);
    weightAvail = true;
  } else {
    if (LogLevel>2) printf(" without weight ");
    weight = NaN;
    weightAvail = false;
  }

  confFile.close();

  if (LogLevel>2) printf("successfully.\n");
  errorState = 0;    
} 


void AnalyzerPhiFieldConfiguration::multiplyHiggsFieldWithConst(double* pField, double fac) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  for (int I=0; I<4*L0*L1*L2*L3; I++) {
    pField[I] *= fac;  
  }
}


void AnalyzerPhiFieldConfiguration::randomizeHiggsField(double* pField) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  for (int I=0; I<4*L0*L1*L2*L3; I++) {
    pField[I] = 2.0 * (AdvancedZufall(AdvancedSeed) - 0.5);  
  }
}


void AnalyzerPhiFieldConfiguration::randomizeGaussHiggsField(double* pField) {
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();

  fermiOps->fillGaussRandomVector((Complex*) pField, 2*L0*L1*L2*L3);
}


void AnalyzerPhiFieldConfiguration::getHiggsFieldDirection(double* pField, vector4D dir) {
  dir[0] = 0;
  dir[1] = 0;
  dir[2] = 0;
  dir[3] = 0;
  
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
    
  int count = 0;
  for (int I=0; I<L0*L1*L2*L3; I++) {
    dir[0] += pField[count+0];
    dir[1] += pField[count+1];
    dir[2] += pField[count+2];
    dir[3] += pField[count+3];
    count += 4;
  }
  
  dir[0] /= (L0*L1*L2*L3);
  dir[1] /= (L0*L1*L2*L3);
  dir[2] /= (L0*L1*L2*L3);
  dir[3] /= (L0*L1*L2*L3);
  if (LogLevel>2) {
    printf("Higgs-Field-Direction: (%f,%f,%f,%f)\n",dir[0],dir[1],dir[2],dir[3]);
  }
}


Complex* AnalyzerPhiFieldConfiguration::performFourierTransform(double* pField, bool forward, int howmanyComponents) {
  if ((howmanyComponents<=0) || (howmanyComponents>4)) {
    printf("ERROR in AnalyzerPhiFieldConfiguration::performFourierTransform: Invalid number for howmanyComponents (%d)\n", howmanyComponents);  
    exit(0);  
  }

  phiFieldCopies[phiFieldCopyCount] = (double*) fermiOps->createFermionVector(howmanyComponents);
  Complex* output = (Complex*)phiFieldCopies[phiFieldCopyCount];
  phiFieldCopyCount++;

  int* n = new int[4];
  n[0] = fermiOps->get1DSizeL0();
  n[1] = fermiOps->get1DSizeL1();
  n[2] = fermiOps->get1DSizeL2();
  n[3] = fermiOps->get1DSizeL3();
  
  int rank = 4;
  int* inembed = new int[4];
  inembed[0] = fermiOps->get1DSizeL0();
  inembed[1] = fermiOps->get1DSizeL1();
  inembed[2] = fermiOps->get1DSizeL2();
  inembed[3] = fermiOps->get1DSizeL3();  
  int istride = howmanyComponents;
  int idist = 1;
  
  int* onembed = new int[4];
  onembed[0] = fermiOps->get1DSizeL0();
  onembed[1] = fermiOps->get1DSizeL1();
  onembed[2] = fermiOps->get1DSizeL2();
  onembed[3] = fermiOps->get1DSizeL3();    
  int ostride = howmanyComponents;
  int odist = 1;

  fftw_plan fftwPlan = NULL;
  if (forward) {
    fftwPlan = fftw_plan_many_dft(rank, n, howmanyComponents, (fftw_complex*) output, inembed,
                                  istride, idist,
	   		          (fftw_complex*)output, onembed, ostride, odist,
	  			  FFTW_FORWARD, FFTW_MEASURE);
  } else {
    fftwPlan = fftw_plan_many_dft(rank, n, howmanyComponents, (fftw_complex*) output, inembed,
                                  istride, idist,
	   		          (fftw_complex*)output, onembed, ostride, odist,
	  			  FFTW_BACKWARD, FFTW_MEASURE);
  }
  delete[] n;
  delete[] inembed;
  delete[] onembed;


  int countO = 0;
  int countP = 0;
  for (int I0=0; I0<fermiOps->get1DSizeL0(); I0++) {
    for (int I1=0; I1<fermiOps->get1DSizeL1(); I1++) {
      for (int I2=0; I2<fermiOps->get1DSizeL2(); I2++) {
        for (int I3=0; I3<fermiOps->get1DSizeL3(); I3++) {
          for (int i=0; i<howmanyComponents; i++) {
	    output[countO].x = pField[countP+i];
	    output[countO].y = 0;
	    countO++;
	  }
	  countP += 4;
        }
      }
    }
  }

  fftw_execute(fftwPlan);

  fftw_destroy_plan(fftwPlan);
  
  return output;
}


void AnalyzerPhiFieldConfiguration::rotateHiggsField(double* pField, int ind1, int ind2, double w) {
  double c = cos(w);
  double s = sin(w);
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
   
  int count = 0;
  for (int I=0; I<L0*L1*L2*L3; I++) {
    double x = pField[count+ind1];
    double y = pField[count+ind2];
    pField[count+ind1] = c*x+s*y;
    pField[count+ind2] = -s*x+c*y;    
    count += 4;
  }
}


void AnalyzerPhiFieldConfiguration::alignHiggsFieldDirection(double* pField) {
  vector4D dir;
  double w;
  getHiggsFieldDirection(pField, dir);
  w = getAngle(dir[1], dir[0]);
  rotateHiggsField(pField, 1, 0, -w);

  getHiggsFieldDirection(pField, dir);
  w = getAngle(dir[2], dir[0]);
  rotateHiggsField(pField, 2, 0, -w);

  getHiggsFieldDirection(pField, dir);
  w = getAngle(dir[3], dir[0]);
  rotateHiggsField(pField, 3, 0, -w);
  
  getHiggsFieldDirection(pField, dir);
}


void AnalyzerPhiFieldConfiguration::randomGaugeRotation(double* pField) {
  Quat rotQ = AdvancedSU2QuatZufall(AdvancedSeed);

  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
   
  int count = 0;
  for (int I=0; I<L0*L1*L2*L3; I++) {
    Quat phiQ(pField[count+0],pField[count+1],pField[count+2],pField[count+3]);
    Quat resQ = phiQ * rotQ;
    
    pField[count+0] = resQ.x0;
    pField[count+1] = resQ.x1;
    pField[count+2] = resQ.x2;
    pField[count+3] = resQ.x3;

    count += 4;
  }
}


void AnalyzerPhiFieldConfiguration::measureMagnetizations() {
  int I1,I2,I3,I4;
  int count = 0;
  vector4D measurePhi, measureStaggeredPhi;
  double avgNorm, sigmaNorm;
  double dummy;
  
  measurePhi[0] = 0;
  measurePhi[1] = 0;
  measurePhi[2] = 0;
  measurePhi[3] = 0;
  measureStaggeredPhi[0] = 0;
  measureStaggeredPhi[1] = 0;
  measureStaggeredPhi[2] = 0;
  measureStaggeredPhi[3] = 0;
  avgNorm = 0;
  sigmaNorm = 0;
  
  int L0 = fermiOps->get1DSizeL0();
  int L1 = fermiOps->get1DSizeL1();
  int L2 = fermiOps->get1DSizeL2();
  int L3 = fermiOps->get1DSizeL3();
   
  count = 0;
  for (I1=0; I1<L0; I1++) {
    for (I2=0; I2<L1; I2++) {
      for (I3=0; I3<L2; I3++) {
        for (I4=0; I4<L3; I4++) {
          measurePhi[0] += phiField[count+0];
          measurePhi[1] += phiField[count+1];
          measurePhi[2] += phiField[count+2];
          measurePhi[3] += phiField[count+3];
	  
	  avgNorm += dummy = sqrt( phiField[count+0]*phiField[count+0] 
	                          +phiField[count+1]*phiField[count+1]
	                          +phiField[count+2]*phiField[count+2]
	                          +phiField[count+3]*phiField[count+3]);
          sigmaNorm += dummy*dummy;
	  
	  double stagFac = 1;
	  if (((I1+I2+I3+I4) % 2) == 1) stagFac = -1;
          measureStaggeredPhi[0] += stagFac*phiField[count+0];
          measureStaggeredPhi[1] += stagFac*phiField[count+1];
          measureStaggeredPhi[2] += stagFac*phiField[count+2];
          measureStaggeredPhi[3] += stagFac*phiField[count+3];
	  count+=4;
	}
      }
    }
  }
  double norm = L0*L1*L2*L3;
  measurePhi[0] = measurePhi[0]/norm;
  measurePhi[1] = measurePhi[1]/norm;
  measurePhi[2] = measurePhi[2]/norm;
  measurePhi[3] = measurePhi[3]/norm;
  measureStaggeredPhi[0] = measureStaggeredPhi[0]/norm;
  measureStaggeredPhi[1] = measureStaggeredPhi[1]/norm;
  measureStaggeredPhi[2] = measureStaggeredPhi[2]/norm;
  measureStaggeredPhi[3] = measureStaggeredPhi[3]/norm;
  
  avgNorm /= norm;
  sigmaNorm = sqrt(sigmaNorm/norm - avgNorm*avgNorm);
  
  magnetizationM = measurePhi[0]*measurePhi[0]
                 + measurePhi[1]*measurePhi[1]
                 + measurePhi[2]*measurePhi[2]
                 + measurePhi[3]*measurePhi[3];
  magnetizationS = measureStaggeredPhi[0]*measureStaggeredPhi[0]
                          + measureStaggeredPhi[1]*measureStaggeredPhi[1]
                          + measureStaggeredPhi[2]*measureStaggeredPhi[2]
                          + measureStaggeredPhi[3]*measureStaggeredPhi[3];
		 
  magnetizationM = sqrt(magnetizationM);
  magnetizationS = sqrt(magnetizationS);
  
  phiFieldAvgVectorLength = avgNorm;
  phiFieldAvgVectorLengthVariation = sigmaNorm;
  avgPhiFieldVector[0] = measurePhi[0];
  avgPhiFieldVector[1] = measurePhi[1];
  avgPhiFieldVector[2] = measurePhi[2];
  avgPhiFieldVector[3] = measurePhi[3];
    
  if (LogLevel>4) {
    printf("m = %1.15f\n", magnetizationM);
    printf("s = %1.15f\n", magnetizationS);
  }
  if (LogLevel>4) printf("\n");
  
  if ((isNaN(magnetizationM)) || (isNaN(magnetizationS))) {
    errorState = 3;
    return;
  }
}


double AnalyzerPhiFieldConfiguration::getMagnetizationM() {
  return magnetizationM;
}


double AnalyzerPhiFieldConfiguration::getMagnetizationS() {
  return magnetizationS;
}


double AnalyzerPhiFieldConfiguration::getPhiFieldAvgVectorLength() {
  return phiFieldAvgVectorLength;
}


double AnalyzerPhiFieldConfiguration::getPhiFieldAvgVectorLengthVariation() {
  return phiFieldAvgVectorLengthVariation;
}


double AnalyzerPhiFieldConfiguration::getAvgPhiFieldVectorComponent(int comp) {
  if (comp<0) return NaN;
  if (comp>3) return NaN;
  return avgPhiFieldVector[comp];
}


double AnalyzerPhiFieldConfiguration::getAvgPhiFieldVectorLength() {
  double norm = sqr(avgPhiFieldVector[0]);
  norm += sqr(avgPhiFieldVector[1]);
  norm += sqr(avgPhiFieldVector[2]);
  norm += sqr(avgPhiFieldVector[3]);
  return sqrt(norm);
}
