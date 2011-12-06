AnalyzerObservableDetSign::AnalyzerObservableDetSign(FermionMatrixOperations* fOps, AnalyzerIOControl* aIOcon, StateDescriptorReader* SDReader) : AnalyzerObservable(fOps, aIOcon, SDReader, "DetSign", "det") { 
        ini(getAnalyzerResultsCount());
	max_iterations = 10;
	curr_eigenvalue_count = 6;
	eigenvalue_increment = 4;
	
	str = new char[2048];
	eigenvalue_prec = 1E-2;
	logLevel = LogLevel; //defined in Global.h
	usePreconditionerP = true;
	extPhiField = NULL;
}


AnalyzerObservableDetSign::~AnalyzerObservableDetSign() {
  delete[] str;
}


bool AnalyzerObservableDetSign::analyze(AnalyzerPhiFieldConfiguration* phiFieldConf, Complex** auxVectors) {
	usePreconditionerP = true;
	double m = phiFieldConf->getMagnetizationM();
	//double m = 0.3;
	fermiOps->setPreconditioner(usePreconditionerP, m, 0.0);
	fermiOps->setQPreconditioner(false, 0.25, 0.25);
	fermiOps->setRPreconditioner(false, 1.0, 1.0);
	
	double* phiField = NULL;
	phiField_2save = NULL;

	if (extPhiField != NULL)
		phiField = extPhiField;
	else
		phiField = phiFieldConf->getPhiFieldCopy();
	
	if (phiField == NULL)
		println(2, "FATAL ERROR: phiField is NULL!");
	
	Complex det(1.0, 0.0);	
	int result = 0;	
	int counter = 0;	
	while (result == 0) {				
		int eVCnt = curr_eigenvalue_count + eigenvalue_increment;
		double prec = eigenvalue_prec / (10.0 * counter + 1.0);
		result = getSignOfDeterminant (*fermiOps, phiField, eVCnt, prec, det);				
		counter++;
	}	
		
	if(analyzerResults != NULL) {
		analyzerResults[0] = (double)result;
		analyzerResults[1] = det.x;
		analyzerResults[2] = det.y;
	}
	
	if (result == 1)
		println(2, "Sign of determinant is positive.");
	
	if (result == -1) {
		println(2, "Sign of determinant is NEGATIVE!.");
		//neg_det = true;
		//phiField_2save = phiField;
	}
	return true;
}


int AnalyzerObservableDetSign::getNeededAuxVectorCount() {
  return 0;
}


int AnalyzerObservableDetSign::getAnalyzerResultsCount() {
  return 3;
}

void AnalyzerObservableDetSign::println(int level, char * format, ...) {	
	va_start (args, format);
	vsprintf (str,format, args);	
	va_end (args);
 	if (level <= logLevel) {
 		printf("%s\n", str); 			 	
	} 	
}

int AnalyzerObservableDetSign::getSignOfDeterminant (FermionMatrixOperations& fermiOps, double* phiField, int initialEigenValueCount, double tol, Complex& det) {
	/*
	PROCEDURE:
	1. Get necessary amount of eigen values (at least 2 more than the min. amount of eigen value, 
		where the largest real part is pos.).
	2. Take the complex product of an even number of eigen values, where the largest real part is 
		positive w.r.t. errors.
	3. If the above product has an imaginary part unequal to zero w.r.t. errors, go back to step 2 
		and consider at least two more eigen values.
	4. The Sign of the determinant of M is given by the sign of the real part of the above product.
	*/
	
//_______________1.)		
	bool usePrecond = true;
	double m, s;
	fermiOps.getPreconditionerParameter(usePrecond, m, s);		
	int what = LARGEST_REAL_VALUE; // depends on preconditioning of M
	double factor = -1.0;	// depends on preconditioning of M
	if (usePrecond) {
		factor = 1.0;
		what = SMALLEST_REAL_VALUE;
	}
			
	int eVCnt = initialEigenValueCount;
	if (eVCnt%2 == 1)
		eVCnt++;		
	if (eVCnt < eigenvalue_increment)
		eVCnt = eigenvalue_increment;
	
	double rel_error = 0.0;
	Complex* ew = NULL;
	bool stop = false;
	while (!stop) {	
		//println(2, "eVCnt=%d, what=%d, tol=%1.4f", eVCnt, what, tol);			
		ew = fermiOps.calcFermionMatrixARPACKEigenValues(what, eVCnt,  phiField, tol, false, NULL, false, true);
		for (int i = eVCnt-1; i >= 4; i--) { 
			rel_error = sqrt(ew[i].x*ew[i].x + ew[i].y*ew[i].y) * tol;					
			if (factor * ew[i].x > rel_error) {
				stop = true;							
			}
		}		
		if (!stop) {
			eVCnt += eigenvalue_increment;
			delete[] ew;
			ew = NULL;
		}						
	}	
	println(2 ,"Current eigen value count: %d. Eigen value margin: max = %+1.8f %+1.8f i (%1.8f) : min = %+1.8f %+1.8f i. STOP=%d", eVCnt, factor*ew[0].x, ew[0].y, rel_error, factor*ew[eVCnt-1].x, ew[eVCnt-1].y, stop);
	// ew is an array of eVCnt complex numbers, where the largest real part is definitely pos.	

//_______________2+3.)
	stop = false;
	double re_old = 1.0;
	double re = 1.0;
	double im = 0.0;
	double abs_det = 0.0;
	double det_error = 0.0;
	
	int counter = 0;
	for (int i = eVCnt-1; i >= 0 && !stop; i--) {
		counter++;
		//complex product
		re_old = re;
		re = 	re * (factor * ew[i].x) - im * ew[i].y;
		im = re_old * ew[i].y + im * (factor * ew[i].x);
			
		// error
		abs_det = sqrt(re*re + im*im);
		det_error = abs_det * tol * sqrt(counter);
		rel_error = sqrt(ew[i].x*ew[i].x + ew[i].y*ew[i].y) * tol;
		
		println(3, "Partial product is : %+1.15f %+1.15f i (error: %1.8f)", re, im, det_error);
		if (factor * ew[i].x > rel_error && fabs(im) < det_error) {
			stop = true;
		}	
	}
	det.x = re;
	det.y = im;
	
//_______________4.)
	bool printEV = true;
	int sign = 0;
	if (stop) {
		if (re + det_error < 0) {			
			sign = -1;			
		}
		else if (re >= det_error)
			sign = +1;
		else
			sign = 0;
	}				
	if (printEV) {
		println(3, "%10s", "Set of eigen values:");
		for(int i = 0; i < eVCnt; i++) {
			rel_error = sqrt(ew[i].x*ew[i].x + ew[i].y*ew[i].y) * tol;
			println(3, "%10s %+1.8f %1.8f i (error: %1.8f)", "->", factor * ew[i].x, ew[i].y, rel_error);
		}
		println(3, "%10s", "end of set");
	}
	
	println(3, "Sign is : %d", sign);
	
	if (sign != 0) {
		curr_eigenvalue_count = counter + 4;
		if (curr_eigenvalue_count < eigenvalue_increment)
			curr_eigenvalue_count = eigenvalue_increment;
	}
	else 
		curr_eigenvalue_count = eVCnt;
	
	delete[] ew;
	return sign;
}

int AnalyzerObservableDetSign::getDeterminant_old (FermionMatrixOperations& fermiOps, double* phiField, int initialEigenValueCount, double tol, Complex& det) {		
	Complex* ew = NULL;	
	double re = 1.0;
	double im = 0.0;
	double re_old = re;
	double rel_error = 0.0;
	double factor = 1.0;	
	
	bool usePrecond = false;
	double m, s;
	fermiOps.getPreconditionerParameter(usePrecond, m, s);
	
	int what = LARGEST_REAL_VALUE;
	if (usePrecond) {
		factor = 1.0;
		what = SMALLEST_REAL_VALUE;
	}
	else 
		factor = -1.0;
	
	// find no of necessary eigen values (eVCnt) dynamically
	int eVCnt = initialEigenValueCount;
	if (eVCnt%2 == 1)
		eVCnt++;
		
	int stepSize = 1;
	double eVPrec = 1E-1;
	bool stop = false;
			
	while (!stop) {		
		ew = fermiOps.calcFermionMatrixARPACKEigenValues(what, eVCnt,  phiField, eVPrec, false, NULL, false, true);
		
		for (int I=0; I<eVCnt; I++) {
  		  rel_error = sqrt(ew[I].x*ew[I].x + ew[I].y*ew[I].y) * eVPrec;		
  		  if (factor * ew[I].x > rel_error) {
			 stop = true;			
		  }
		}
		if (!stop) {
			eVCnt += 10*stepSize;
		}
		println(3 ,"Current eigen value count: %d. Eigen value margin: max = %+1.8f %+1.8f i (%1.8f) : min = %+1.8f %+1.8f i", eVCnt, factor*ew[0].x, ew[0].y, rel_error, factor*ew[eVCnt-1].x, ew[eVCnt-1].y);												
		delete[] ew;
		ew = NULL;
	}
			 		
	eVPrec = tol;	
	re = 1.0;
	im = 0.0;
	int real_axis_ev_count = 0;
	double abs_det = 1.0;
	double det_error = 0.0;
	double det_sign = 1.0;
	
	ew = fermiOps.calcFermionMatrixARPACKEigenValues(what, eVCnt,  phiField, eVPrec, false, NULL, false, true);				
	for (int i = 0; i < eVCnt; i++) {
		ew[i].x = ew[i].x * factor;
		rel_error = sqrt(ew[i].x*ew[i].x + ew[i].y*ew[i].y)*eVPrec;
		if ( fabs(ew[i].y) <= rel_error) {
			real_axis_ev_count++;
			
			re_old = re;
			re = 	re * ew[i].x - im * ew[i].y;
			im = re_old * ew[i].y + im * ew[i].x;
			
			det_sign *= ew[i].x;
			
			abs_det = sqrt(re*re + im*im);
			det_error = abs_det * eVPrec * sqrt(real_axis_ev_count);
			println(3, "Eigenvalue along real axis: %+1.8f %+1.8f i, curr partial product: (%+1.8f %1.8f i error: %1.8f)", ew[i].x, ew[i].y, re, im, det_error);			
		}
		else
			println(4, "Eigenvalue : %+1.8f %+1.8f i", ew[i].x, ew[i].y);		
	}	
	println(4, "Partial product is : %+1.15f %+1.15f i (error: %1.8f)", re, im, det_error);
	bool printEV=false;	
	int result = 0;
	if (re + det_error < 0) {
		println(3, "Determinant is negative! %+1.15f %+1.15f i (%1.8f)", re, im, det_error);
		printEV=false;			
		result = -1;
	} else if (re - det_error > 0) {
		println(4, "positive det: %+1.15f %+1.15f i (%1.8f), precision: %1.3f, eigen value count: %d", re, im, det_error, eVPrec, eVCnt);
		result = +1;
	} else {
		println(3, "UNDETERMINED DET: %+1.15f %+1.15f i (%1.8f), precision: %1.3f, eigen value count: %d", re, im, det_error, eVPrec, eVCnt);
		printEV = true;
		result = 0;
	}

	if (printEV) {
		println(3, "%10s", "Set of eigen values:");
		for(int i = 0; i < eVCnt; i++) {
			println(3, "%10s %+1.8f %1.8f i", "->", ew[i].x, ew[i].y);
		}
		println(3, "%10s", "end of set");
	}
	
	if (result != 0) {
		curr_eigenvalue_count = eVCnt - 10*stepSize;
		if (curr_eigenvalue_count <= 4)
			curr_eigenvalue_count = 4;
	}
	else 
		curr_eigenvalue_count = eVCnt;
	
	delete[] ew;
	det.x = re;
	det.y = im;
	return result;
}
