#ifndef SSEroutines_included
#define SSEroutines_included


#include "Complex.h"
#include "Tools.h"





/**** SSE_FourierTrafoForthFFTstep ***
* Performs the third FFT step. 
* ATTENTION: Restriction to oneDim = 128. Input and Output X3-Increments
* are fixed. Input and Output are exprected to have same size.
* Input and Output must not be the same!!!
*
* Parameters: The following data are expected in ASMParameterData
* OneDim                : Lattice size in one dimension
* data                  : Input-/Output-Vector
* Increments           : Increments for output base and output pi-shifted
* ASMParameterData[2576...] : Relative Lese-, Schreibe-Adressen
*
**/
void SSE_FourierTrafoForthFFTstep(); 


/**** SSE_FourierTrafoThirdFFTstep ***
* Performs the third FFT step. 
* ATTENTION: Restriction to oneDim = 128. Input and Output X3-Increments
* are fixed.
*
* Parameters: The following data are expected in ASMParameterData
* OneDim                : Lattice size in one dimension
* data                  : Input-/Output-Vector
* Incerements           : Increments for output base and output pi-shifted
*
**/
void SSE_FourierTrafoThirdFFTstep(); 


/**** SSE_FourierTrafoSecondFFTstep ***
* Performs the second FFT step. Input data are expected to be in fixed
* size array.
* ATTENTION: Restriction to oneDim = 128
*
* Parameters: The following data are expected in ASMParameterData
* OneDim                : Lattice size in one dimension
* data                  : Input-/Output-Vector
* Incerements           : Increments for output base and output pi-shifted
*
**/
void SSE_FourierTrafoSecondFFTstep(); 


/**** SSE_FourierTrafoRearrangerWithFirstFFTstep ***
* Rearranges the input-vector as prescribed by the FFT-algorithm.
* ATTENTION: Restriction to oneDim = 128
*
* Parameters: The following data are expected in ASMParameterData
* OneDim                : Lattice size in one dimension
* data                  : Input-/Output-Vector
* BitInverse            : Bit-Inversion data
* b0,b1,b2,b3           : Integer position of this block (beginning)
*
**/
void SSE_FourierTrafoRearrangerWithFirstFFTstep(); 


/**** SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace ***
* Performs the application of the inverse of MdagM for static Phi in the Fourier-Space.
*
* Parameters:
* oneDim                   : Lattice size in one dimension
* input                    : Input Vector
* output                   : Output Vector
* sinP                     : Table, containing tilde p
* auxData                  : Table, containing eigenvalues and eigenvalues divided by sqrt(p^2)...
*
**/
void SSE_ExecuteFermionMatrixStaticInverseMultiplicationInFourierSpace(int oneDimL0, int oneDimL1, int oneDimL2, int oneDimL3, Complex* input, Complex* output, Complex* sinP, Complex* auxData, long int* PiPiPiPiIndices); 


/**** SSE_Performf_YB_2rho_AndCopyToOutput ***
* Applies yNB^+/2rho -1 to x and stores result in interim and stores yN B^+ applied on x in output.
* ATTENTION: Watch out for factors 2rho!!!
*
* Parameters:
* Lto4                  : Lattice size in one dimension
* yN                    : Parameter yN
* twoRho                : 2 * Parameter Rho
* x                     : Input Vector
* phi                   : Phi - Feld
* interim               : Zwischenspeichervektor
* output                : Output Vector
*
**/
void SSE_Performf_YB_2rho_AndCopyToOutput(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double yN, double twoRho, Complex* x, double* phi, Complex* interim, Complex* output); 


/**** SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace ***
* Performs the application of the Dirac-Operator in the Fourier-Space.
*
*            work[0] = work[0] + auxData[countAux].y*(Complex(p3,p0)*input[count+2] + Complex(p1,-p2)*input[count+3]);
*            work[1] = work[1] + auxData[countAux].y*(Complex(p1,p2)*input[count+2] + Complex(-p3,p0)*input[count+3]);
*            work[2] = work[2] + auxData[countAux].y*(Complex(-p3,p0)*input[count+0] + Complex(-p1,p2)*input[count+1]);
*            work[3] = work[3] + auxData[countAux].y*(Complex(-p1,-p2)*input[count+0] + Complex(p3,p0)*input[count+1]);
*
* Parameters:
* oneDim                   : Lattice size in one dimension
* input                    : Input Vector
* output                   : Output Vector
* sinP                     : Table, containing tilde p
* auxData                  : Table, containing eigenvalues and eigenvalues divided by sqrt(p^2)...
*
**/
void SSE_ExecuteNeubergerDiracMultiplicationInFourierSpace(int oneDimL0, int oneDimL1, int oneDimL2, int oneDimL3, Complex* input, Complex* output, Complex* sinP, Complex* auxData); 


/**** SSE_ComplexVectorAdditionSPECIAL2 ***
* Performs a complex vector subtraction including a real multiplication with a constant
* acting on the result.
* y = beta*( x - y )
* Parameters:
* N                   : Vector length (must be dividable by 8)
* beta                : Real constant
* x                   : Pointer to complex vector x
* y                   : Pointer to complex vector y
*
**/
void SSE_ComplexVectorAdditionSPECIAL2(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double beta, Complex* x, Complex* y); 


/**** SSE_ComplexVectorAdditionSPECIAL1 ***
* Performs a complex vector addition including a real multiplication with a constant
* acting on the vector y.
* y = beta*y + x
* Parameters:
* N                   : Vector length (must be dividable by 8)
* beta                : Real constant
* x                   : Pointer to complex vector x
* y                   : Pointer to complex vector y
*
**/
void SSE_ComplexVectorAdditionSPECIAL1(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, double beta, Complex* x, Complex* y); 


/**** SSE_ComplexVectorAdditionWithSquaredNorm ***
* Performs a complex vector addition including a complex multiplication with a constant
* and calculates the squared norm of the result.
* y = y + alpha*x
* Parameters:
* N                   : Vector length (must be dividable by 4)
* alpha               : Complex constant
* x                   : Pointer to complex vector x
* y                   : Pointer to complex vector y
* SqrNorm             : Squared norm of result
*
**/
void SSE_ComplexVectorAdditionWithSquaredNorm(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, Complex& alpha, Complex* x, Complex* y, double& SqrNorm); 


/**** SSE_ComplexVectorAddition ***
* Performs a complex vector addition including a complex multiplication with a constant.
* y = y + alpha*x
* Parameters:
* N                   : Vector length (must be dividable by 4)
* alpha               : Complex constant
* x                   : Pointer to complex vector x
* y                   : Pointer to complex vector y
*
**/
void SSE_ComplexVectorAddition(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, Complex& alpha, Complex* x, Complex* y); 


/**** SSE_ComplexScalarProduct ***
* Calcs the scalar product of to vectors.
* Parameters:
* N                   : Vector length (must be dividable by 4)
* v1                  : Pointer to complex vector 1
* v2                  : Pointer to complex vector 2
* r                   : Complex results
*
**/
void SSE_ComplexScalarProduct(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, Complex* v1, Complex* v2, Complex& r); 


/**** SSE_ComplexSqrNorm ***
* Calcs the squared norm of a complex vector.
* Parameters:
* N                   : Vector length (must be dividable by 8)
* v                   : Pointer to complex vector
* r                   : Complex results
*
**/
void SSE_ComplexSquareNorm(int L0, int L1, int L2, int L3, int xtrS1, int xtrS2, int xtrS3, Complex* v, double& r); 


/**** SSE_ZCopy ***
*
**/
void SSE_ZCopy(int count, Complex* source, int j1, Complex* dest, int j2); 


#endif
