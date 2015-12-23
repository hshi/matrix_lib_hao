#include "matrix_all.h"

using namespace std;

namespace matrix_hao_lib
{


 /*************************************/
 /*Matrix Multiply C=alpha*A.B+beta*C */
 /*************************************/

 void gmm_cpu(const Matrix<float,2>& A, const Matrix<float,2>& B, Matrix<float,2>& C, 
          char TRANSA, char TRANSB, float alpha, float beta)
 {
     BL_INT  M, N, K, LDA, LDB, LDC;
     M=(TRANSA=='N') ? A.L1:A.L2;
     K=(TRANSA=='N') ? A.L2:A.L1;
     N=(TRANSB=='N') ? B.L2:B.L1;
     LDA=A.L1;
     LDB=B.L1;
     LDC=C.L1;
     FORTRAN_NAME(sgemm)(&TRANSA, &TRANSB, &M, &N, &K, 
                         (BL_FLOAT* )&alpha, (BL_FLOAT* )A.base_array, &LDA, 
                         (BL_FLOAT* )B.base_array, &LDB, 
                         (BL_FLOAT* )&beta,  (BL_FLOAT* )C.base_array, &LDC);
 }



 void gmm_cpu(const Matrix<double,2>& A, const Matrix<double,2>& B, Matrix<double,2>& C,
          char TRANSA, char TRANSB, double alpha, double beta)
 {
     BL_INT  M, N, K, LDA, LDB, LDC;
     M=(TRANSA=='N') ? A.L1:A.L2;
     K=(TRANSA=='N') ? A.L2:A.L1;
     N=(TRANSB=='N') ? B.L2:B.L1;
     LDA=A.L1;
     LDB=B.L1;
     LDC=C.L1;
     FORTRAN_NAME(dgemm)(&TRANSA, &TRANSB, &M, &N, &K,
                         (BL_DOUBLE* )&alpha, (BL_DOUBLE* )A.base_array, &LDA,
                         (BL_DOUBLE* )B.base_array, &LDB,
                         (BL_DOUBLE* )&beta,  (BL_DOUBLE* )C.base_array, &LDC);
 }


 void gmm_cpu(const Matrix<complex<float>,2>& A, const Matrix<complex<float>,2>& B, Matrix<complex<float>,2>& C,
          char TRANSA, char TRANSB, complex<float> alpha,complex<float> beta)
 {
     BL_INT  M, N, K, LDA, LDB, LDC;
     M=(TRANSA=='N') ? A.L1:A.L2;
     K=(TRANSA=='N') ? A.L2:A.L1;
     N=(TRANSB=='N') ? B.L2:B.L1;
     LDA=A.L1;
     LDB=B.L1;
     LDC=C.L1;
     FORTRAN_NAME(cgemm)(&TRANSA, &TRANSB, &M, &N, &K,
                         (BL_COMPLEX8* )&alpha, (BL_COMPLEX8* )A.base_array, &LDA,
                         (BL_COMPLEX8* )B.base_array, &LDB,
                         (BL_COMPLEX8* )&beta,  (BL_COMPLEX8* )C.base_array, &LDC);
 }

 void gmm_cpu(const Matrix<complex<double>,2>& A, const Matrix<complex<double>,2>& B, Matrix<complex<double>,2>& C,
          char TRANSA, char TRANSB, complex<double> alpha, complex<double> beta)
 {
     BL_INT  M, N, K, LDA, LDB, LDC;
     M=(TRANSA=='N') ? A.L1:A.L2;
     K=(TRANSA=='N') ? A.L2:A.L1;
     N=(TRANSB=='N') ? B.L2:B.L1;
     LDA=A.L1;
     LDB=B.L1;
     LDC=C.L1;
     FORTRAN_NAME(zgemm)(&TRANSA, &TRANSB, &M, &N, &K,
                         (BL_COMPLEX16* )&alpha, (BL_COMPLEX16* )A.base_array, &LDA,
                         (BL_COMPLEX16* )B.base_array, &LDB,
                         (BL_COMPLEX16* )&beta,  (BL_COMPLEX16* )C.base_array, &LDC);
 }


 /******************************************/
 /*Diagonalize Real symmetric Matrix********/
 /******************************************/
 void eigen_cpu(Matrix<double,2>& A, Matrix<double,1>& W, char JOBZ, char UPLO)
 {
     if(A.L1!=A.L2) {cout<<"Input for eigen is not square matrix!\n"; exit(1);}
     BL_INT N=A.L1; BL_INT info;
     double work_test[1]; int iwork_test[1];
     BL_INT lwork=-1; BL_INT liwork=-1;

     FORTRAN_NAME(dsyevd)(&JOBZ,&UPLO,&N,(BL_DOUBLE* )A.base_array,&N,(BL_DOUBLE* )W.base_array,
                          (BL_DOUBLE* )work_test, &lwork, (BL_INT* )iwork_test, &liwork ,&info);

     lwork=lround(work_test[0]); liwork=iwork_test[0];
     double* work= new double[lwork]; int* iwork=new int[liwork];
     FORTRAN_NAME(dsyevd)(&JOBZ,&UPLO,&N,(BL_DOUBLE* )A.base_array,&N,(BL_DOUBLE* )W.base_array,
                          (BL_DOUBLE* )work, &lwork, (BL_INT* )iwork, &liwork ,&info);

     delete[] work; delete[] iwork;

     if(info!=0) {cout<<"Dsyevd failed: info= "<< info<<"\n"; exit(1);}
 }



 /*************************************/
 /*Diagonalize Hermitian Matrix********/
 /*************************************/
 void eigen_cpu(Matrix<complex<double>,2>& A, Matrix<double,1>& W, char JOBZ, char UPLO)
 {
     if(A.L1!=A.L2) {cout<<"Input for eigen is not square matrix!\n"; exit(1);}
     BL_INT N=A.L1; BL_INT info; 
     complex<double> work_test[1]; double rwork_test[1]; int iwork_test[1]; 
     BL_INT lwork=-1; BL_INT lrwork=-1; BL_INT liwork=-1;

     FORTRAN_NAME(zheevd)(&JOBZ,&UPLO,&N,(BL_COMPLEX16* )A.base_array,&N,(BL_DOUBLE* )W.base_array,
      (BL_COMPLEX16* )work_test, &lwork, (BL_DOUBLE* )rwork_test, &lrwork, (BL_INT* )iwork_test, &liwork ,&info);

     lwork=lround(work_test[0].real()); lrwork=lround(rwork_test[0]); liwork=iwork_test[0];
     complex<double>* work= new complex<double>[lwork]; double* rwork=new double[lrwork]; int* iwork=new int[liwork];
     FORTRAN_NAME(zheevd)(&JOBZ,&UPLO,&N,(BL_COMPLEX16* )A.base_array,&N,(BL_DOUBLE* )W.base_array,
      (BL_COMPLEX16* )work, &lwork, (BL_DOUBLE* )rwork, &lrwork, (BL_INT* )iwork, &liwork ,&info);

     delete[] work; delete[] rwork; delete[] iwork;

     if(info!=0) {cout<<"Zheevd failed: info= "<< info<<"\n"; exit(1);}
 }


 /******************************************/
 /*LU Decomposition a complex square Matrix*/
 /******************************************/
 LUDecomp<complex<double>> LUconstruct_cpu(const Matrix<complex<double>,2>& x)
 {
     if(x.L1!=x.L2) {cout<<"Input for LU is not square matrix!\n"; exit(1);}
     LUDecomp<complex<double>> y; y.A=x; y.ipiv=Matrix<BL_INT,1>(x.L1);
     BL_INT N=x.L1;
     FORTRAN_NAME(zgetrf)(&N,&N,(BL_COMPLEX16* )y.A.base_array,&N,y.ipiv.base_array,&y.info);
     if(y.info<0) {cout<<"The "<<y.info<<"-th parameter is illegal!\n"; exit(1);}
     return y; 
 }
 

 
 /*********************************************************************************************************************/
 /*Inverse of  Matrix: If determinant of the matrix is outof machine precision, inverse should be fine, since it solve*
  *The linear equation, every small value is well defined                                                             */
 /*********************************************************************************************************************/
 void inverse_cpu_utilities(const Matrix<complex<double>,2>& A, const LUDecomp<complex<double>>& x)
 {
     BL_INT N=A.L1; BL_INT info;

     BL_INT lwork=-1; complex<double> work_test[1];
     FORTRAN_NAME(zgetri)(&N,(BL_COMPLEX16* )A.base_array,&N,x.ipiv.base_array,(BL_COMPLEX16* )work_test,&lwork,&info);

     lwork=lround(work_test[0].real());
     complex<double>* work= new complex<double>[lwork];
     FORTRAN_NAME(zgetri)(&N,(BL_COMPLEX16* )A.base_array,&N,x.ipiv.base_array,(BL_COMPLEX16* )work,&lwork,&info);
     if(info<0) {cout<<"The "<<info<<"-th parameter is illegal!\n"; exit(1);}
     delete[] work;
 }

 Matrix<complex<double>,2> inverse_cpu(const LUDecomp<complex<double>>& x)
 {
     Matrix<complex<double>,2> A=x.A; //We know x.A own the matrix
     inverse_cpu_utilities(A, x); 
     return A;
 }

 Matrix<complex<double>,2> inverse_cpu(LUDecomp<complex<double>>&& x)
 {
     Matrix<complex<double>,2> A=move(x.A); //We know x.A own the matrix
     inverse_cpu_utilities(A, x);
     return A;
 }

 /******************************************************/
 /*Solve Linear Equation of the matrix A*M=B: M=A^{-1}B*/
 /******************************************************/
 Matrix<complex<double>,2> solve_lineq_cpu(const LUDecomp<complex<double>>& x, const Matrix<complex<double>,2>& B, char TRANS)
 {
     if(x.A.L1!=B.L1)  {cout<<"Input size for solving linear equation is not consistent!\n"; exit(1);}
     Matrix<complex<double>,2> M; M=B;  //B might not own itself, we'd better use equal assigment instead of equal construction.
     BL_INT N=B.L1; BL_INT NRHS=B.L2; BL_INT info;
     FORTRAN_NAME(zgetrs)(&TRANS,&N,&NRHS,(BL_COMPLEX16* )x.A.base_array,&N,x.ipiv.base_array,(BL_COMPLEX16* )M.base_array,&N,&info);
     if(info!=0) 
     {
         cout<<"Solve linear equation is not suceesful: "<<info<<"-th parameter is illegal! \n"; 
         exit(1);
     }
     return M;
 }


 /******************************/
 /*QR decompostion of matrix ph*/
 /******************************/
 double QRMatrix_cpu(Matrix<complex<double>,2>& ph)
 {
     BL_INT L=ph.L1; BL_INT N=ph.L2; BL_INT info;
     BL_INT lwork=-1; complex<double> work_test[1];
     complex<double>* tau= new complex<double>[N];

     FORTRAN_NAME(zgeqrf) (&L,&N,(BL_COMPLEX16* )ph.base_array,&L,(BL_COMPLEX16* )tau,(BL_COMPLEX16* )work_test,&lwork,&info);

     lwork=lround(work_test[0].real());
     complex<double>* work= new complex<double>[lwork];
     FORTRAN_NAME(zgeqrf) (&L,&N,(BL_COMPLEX16* )ph.base_array,&L,(BL_COMPLEX16* )tau,(BL_COMPLEX16* )work,&lwork,&info);
     if(info!=0) {cout<<"QR run is not suceesful: "<<info<<"-th parameter is illegal! \n"; exit(1);}

     complex<double> det={1.0,0.0}; for (size_t i=0; i<ph.L2; i++)  det*=ph(i,i); 

     FORTRAN_NAME(zungqr) (&L,&N,&N,(BL_COMPLEX16* )ph.base_array,&L,(BL_COMPLEX16* )tau,(BL_COMPLEX16* )work,&lwork,&info);
     if(info!=0) {cout<<"Zungqr run is not suceesful: "<<info<<"-th parameter is illegal! \n"; exit(1);}

     if(det.real()<0) {det=-det; for(size_t i=0; i<ph.L1; i++) ph(i,0)=-ph(i,0);}

     delete[] tau;delete[] work;

     return det.real();
 }

 double QRMatrix_cpu(Matrix<complex<double>,2>& ph, vector<double>& detVec)
 {
     BL_INT L=ph.L1; BL_INT N=ph.L2; BL_INT info;
     BL_INT lwork=-1; complex<double> work_test[1];
     complex<double>* tau= new complex<double>[N];

     FORTRAN_NAME(zgeqrf) (&L,&N,(BL_COMPLEX16* )ph.base_array,&L,(BL_COMPLEX16* )tau,(BL_COMPLEX16* )work_test,&lwork,&info);

     lwork=lround(work_test[0].real());
     complex<double>* work= new complex<double>[lwork];
     FORTRAN_NAME(zgeqrf) (&L,&N,(BL_COMPLEX16* )ph.base_array,&L,(BL_COMPLEX16* )tau,(BL_COMPLEX16* )work,&lwork,&info);
     if(info!=0) {cout<<"QR run is not suceesful: "<<info<<"-th parameter is illegal! \n"; exit(1);}

     detVec.resize(ph.L2); complex<double> det={1.0,0.0};
     for (size_t i=0; i<ph.L2; i++)  {detVec[i]=ph(i,i).real(); det*=ph(i,i);}

     FORTRAN_NAME(zungqr) (&L,&N,&N,(BL_COMPLEX16* )ph.base_array,&L,(BL_COMPLEX16* )tau,(BL_COMPLEX16* )work,&lwork,&info);
     if(info!=0) {cout<<"Zungqr run is not suceesful: "<<info<<"-th parameter is illegal! \n"; exit(1);}

     if(det.real()<0) 
     {
        det=-det; 
        detVec[0]=-detVec[0]; 
        for(size_t i=0; i<ph.L1; i++) ph(i,0)=-ph(i,0);
     }

     delete[] tau;delete[] work;

     return det.real();
 }

} //end namespace matrix_hao_lib
