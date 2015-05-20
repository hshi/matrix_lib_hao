#ifndef MATRIX_HAO_BL_CPU_H
#define MATRIX_HAO_BL_CPU_H


namespace matrix_hao_lib
{

 /*************************************/
 /*Matrix Multiply C=alpha*A.B+beta*C */
 /*************************************/

 void gmm_cpu(const Matrix<float,2>& A, const Matrix<float,2>& B, Matrix<float,2>& C, 
          char TRANSA='N', char TRANSB='N', float alpha=1, float beta=0);

 void gmm_cpu(const Matrix<double,2>& A, const Matrix<double,2>& B, Matrix<double,2>& C,
          char TRANSA='N', char TRANSB='N', double alpha=1, double beta=0);

 void gmm_cpu(const Matrix<std::complex<float>,2>& A, const Matrix<std::complex<float>,2>& B, Matrix<std::complex<float>,2>& C,
          char TRANSA='N', char TRANSB='N', std::complex<float> alpha=1,std::complex<float> beta=0);

 void gmm_cpu(const Matrix<std::complex<double>,2>& A, const Matrix<std::complex<double>,2>& B, Matrix<std::complex<double>,2>& C,
          char TRANSA='N', char TRANSB='N', std::complex<double> alpha=1, std::complex<double> beta=0);


 /*****************************************************/
 /*Diagonalize Real symmetric/ Hermitian Matrix********/
 /*****************************************************/
 void eigen_cpu(Matrix<double,2>& A, Matrix<double,1>& W, char JOBZ='V', char UPLO='U');
 void eigen_cpu(Matrix<std::complex<double>,2>& A, Matrix<double,1>& W, char JOBZ='V', char UPLO='U');


 /*************************************/
 /*Construct LUDecomp with CPU ********/
 /*************************************/
 LUDecomp< std::complex<double> > LUconstruct_cpu(const Matrix<std::complex<double>,2>& x);


 /********************/
 /*Inverse of  Matrix*/
 /********************/
 Matrix<std::complex<double>,2> inverse_cpu(const LUDecomp<std::complex<double>>& x);


 /******************************************************/
 /*Solve Linear Equation of the matrix A*M=B: M=A^{-1}B*/
 /******************************************************/
 Matrix<std::complex<double>,2> solve_lineq_cpu(const LUDecomp<std::complex<double>>& x, const Matrix<std::complex<double>,2>& B
                                           ,char TRANS='N');


 /***********************************************************/
 /*QR decompostion of matrix ph, return the determinant of R*/
 /***********************************************************/
 double QRMatrix_cpu(Matrix<std::complex<double>,2>& ph);


}//end namespace matrix_hao_lib

#endif
