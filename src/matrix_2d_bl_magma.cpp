#ifdef USE_MAGMA
#include "matrix_all.h"

using namespace std;

namespace matrix_hao_lib
{

 static magmaFloatComplex  _cast_C(const complex<float>&  Z)  { return MAGMA_C_MAKE(real(Z), imag(Z)); }
 static magmaDoubleComplex _cast_Z(const complex<double>& Z)  { return MAGMA_Z_MAKE(real(Z), imag(Z)); }


 /*************************************/
 /*Matrix Multiply C=alpha*A.B+beta*C */
 /*************************************/
 //Magma_*gemm only support GPU interface.

 void gmm_magma(const Matrix<float,2>& A, const Matrix<float,2>& B, Matrix<float,2>& C, 
          char TRANSA, char TRANSB, float alpha, float beta)
 {
     magma_int_t M, N, K, LDA, LDB, LDC;
     magma_trans_t transA=magma_trans_const(TRANSA), transB=magma_trans_const(TRANSB);
     magmaFloat_ptr d_A, d_B, d_C;

     //Set LDA, LDB, and LDC, round up to multiple of 32 for best GPU performance
     LDA = ((A.L1+31)/32)*32; LDB = ((B.L1+31)/32)*32; LDC = ((C.L1+31)/32)*32;   

     // Allocate memory for the matrices on GPU 
     magma_smalloc(&d_A, LDA*A.L2 );
     magma_smalloc(&d_B, LDB*B.L2 );
     magma_smalloc(&d_C, LDC*C.L2 );

     // Copy data from host (CPU) to device (GPU)
     magma_ssetmatrix( A.L1, A.L2, A.base_array, A.L1, d_A, LDA );
     magma_ssetmatrix( B.L1, B.L2, B.base_array, B.L1, d_B, LDB );
     if( abs(beta)>1e-32 ) magma_ssetmatrix( C.L1, C.L2, C.base_array, C.L1, d_C, LDC );

     //Call magma_sgemm
     M=(TRANSA=='N') ? A.L1:A.L2;
     K=(TRANSA=='N') ? A.L2:A.L1;
     N=(TRANSB=='N') ? B.L2:B.L1;
     magma_sgemm(transA, transB, M, N, K, alpha, d_A, LDA, d_B, LDB, beta,d_C, LDC);

     // Copy solution from device (GPU) to host (CPU)
     magma_sgetmatrix(C.L1, C.L2, d_C, LDC, C.base_array, C.L1);

     // Free memory on GPU
     magma_free(d_A); magma_free(d_B); magma_free(d_C);
 }


 void gmm_magma(const Matrix<double,2>& A, const Matrix<double,2>& B, Matrix<double,2>& C,
          char TRANSA, char TRANSB, double alpha, double beta)
 {
     magma_int_t M, N, K, LDA, LDB, LDC;
     magma_trans_t transA=magma_trans_const(TRANSA), transB=magma_trans_const(TRANSB);
     magmaDouble_ptr d_A, d_B, d_C;

     //Set LDA, LDB, and LDC, round up to multiple of 32 for best GPU performance
     LDA = ((A.L1+31)/32)*32; LDB = ((B.L1+31)/32)*32; LDC = ((C.L1+31)/32)*32;

     // Allocate memory for the matrices on GPU 
     magma_dmalloc(&d_A, LDA*A.L2 );
     magma_dmalloc(&d_B, LDB*B.L2 );
     magma_dmalloc(&d_C, LDC*C.L2 );

     // Copy data from host (CPU) to device (GPU)
     magma_dsetmatrix( A.L1, A.L2, A.base_array, A.L1, d_A, LDA );
     magma_dsetmatrix( B.L1, B.L2, B.base_array, B.L1, d_B, LDB );
     if( abs(beta)>1e-32 ) magma_dsetmatrix( C.L1, C.L2, C.base_array, C.L1, d_C, LDC );

     //Call magma_sgemm
     M=(TRANSA=='N') ? A.L1:A.L2;
     K=(TRANSA=='N') ? A.L2:A.L1;
     N=(TRANSB=='N') ? B.L2:B.L1;
     magma_dgemm(transA, transB, M, N, K, alpha, d_A, LDA, d_B, LDB, beta,d_C, LDC);

     // Copy solution from device (GPU) to host (CPU)
     magma_dgetmatrix(C.L1, C.L2, d_C, LDC, C.base_array, C.L1);

     // Free memory on GPU
     magma_free(d_A); magma_free(d_B); magma_free(d_C);
 } 


 void gmm_magma(const Matrix<complex<float>,2>& A, const Matrix<complex<float>,2>& B, Matrix<complex<float>,2>& C,
          char TRANSA, char TRANSB, complex<float> alpha,complex<float> beta)
 {
     magma_int_t M, N, K, LDA, LDB, LDC;
     magma_trans_t transA=magma_trans_const(TRANSA), transB=magma_trans_const(TRANSB);
     magmaFloatComplex_ptr d_A, d_B, d_C;

     //Set LDA, LDB, and LDC, round up to multiple of 32 for best GPU performance
     LDA = ((A.L1+31)/32)*32; LDB = ((B.L1+31)/32)*32; LDC = ((C.L1+31)/32)*32;

     // Allocate memory for the matrices on GPU 
     magma_cmalloc(&d_A, LDA*A.L2 );
     magma_cmalloc(&d_B, LDB*B.L2 );
     magma_cmalloc(&d_C, LDC*C.L2 );

     // Copy data from host (CPU) to device (GPU)
     magma_csetmatrix( A.L1, A.L2, (magmaFloatComplex* ) A.base_array, A.L1, d_A, LDA );
     magma_csetmatrix( B.L1, B.L2, (magmaFloatComplex* ) B.base_array, B.L1, d_B, LDB );
     if( abs(beta)>1e-32 ) magma_csetmatrix( C.L1, C.L2, (magmaFloatComplex* ) C.base_array, C.L1, d_C, LDC );

     //Call magma_sgemm
     M=(TRANSA=='N') ? A.L1:A.L2;
     K=(TRANSA=='N') ? A.L2:A.L1;
     N=(TRANSB=='N') ? B.L2:B.L1;
     magma_cgemm(transA, transB, M, N, K, _cast_C(alpha), d_A, LDA, d_B, LDB, _cast_C(beta),d_C, LDC);

     // Copy solution from device (GPU) to host (CPU)
     magma_cgetmatrix(C.L1, C.L2, d_C, LDC, (magmaFloatComplex* ) C.base_array, C.L1);

     // Free memory on GPU
     magma_free(d_A); magma_free(d_B); magma_free(d_C);
 }


 void gmm_magma(const Matrix<complex<double>,2>& A, const Matrix<complex<double>,2>& B, Matrix<complex<double>,2>& C,
          char TRANSA, char TRANSB, complex<double> alpha, complex<double> beta)
 {
     magma_int_t M, N, K, LDA, LDB, LDC;
     magma_trans_t transA=magma_trans_const(TRANSA), transB=magma_trans_const(TRANSB);
     magmaDoubleComplex_ptr d_A, d_B, d_C;

     //Set LDA, LDB, and LDC, round up to multiple of 32 for best GPU performance
     LDA = ((A.L1+31)/32)*32; LDB = ((B.L1+31)/32)*32; LDC = ((C.L1+31)/32)*32;

     // Allocate memory for the matrices on GPU 
     magma_zmalloc(&d_A, LDA*A.L2 );
     magma_zmalloc(&d_B, LDB*B.L2 );
     magma_zmalloc(&d_C, LDC*C.L2 );

     // Copy data from host (CPU) to device (GPU)
     magma_zsetmatrix( A.L1, A.L2, (magmaDoubleComplex* ) A.base_array, A.L1, d_A, LDA );
     magma_zsetmatrix( B.L1, B.L2, (magmaDoubleComplex* ) B.base_array, B.L1, d_B, LDB );
     if( abs(beta)>1e-32 ) magma_zsetmatrix( C.L1, C.L2, (magmaDoubleComplex* ) C.base_array, C.L1, d_C, LDC );

     //Call magma_sgemm
     M=(TRANSA=='N') ? A.L1:A.L2;
     K=(TRANSA=='N') ? A.L2:A.L1;
     N=(TRANSB=='N') ? B.L2:B.L1;
     magma_zgemm(transA, transB, M, N, K, _cast_Z(alpha), d_A, LDA, d_B, LDB, _cast_Z(beta),d_C, LDC);

     // Copy solution from device (GPU) to host (CPU)
     magma_zgetmatrix(C.L1, C.L2, d_C, LDC, (magmaDoubleComplex* ) C.base_array, C.L1);

     // Free memory on GPU
     magma_free(d_A); magma_free(d_B); magma_free(d_C);
 }


 /******************************************/
 /*Diagonalize Real symmetric Matrix********/
 /******************************************/
 void eigen_magma(Matrix<double,2>& A, Matrix<double,1>& W, char JOBZ, char UPLO)
 {
     if(A.L1!=A.L2) {cout<<"Input for eigen is not square matrix!\n"; exit(1);}

     magma_vec_t jobz = magma_vec_const(JOBZ); magma_uplo_t uplo = magma_uplo_const(UPLO);
     magma_int_t N=A.L1; magma_int_t info;

     double work_test[1]; magma_int_t iwork_test[1];
     magma_int_t lwork=-1; magma_int_t liwork=-1;   
     magma_dsyevd( jobz, uplo, N, A.base_array, N, W.base_array, work_test, lwork, iwork_test, liwork, &info );

     lwork=lround(work_test[0]); liwork=iwork_test[0];
     double* work; magma_int_t* iwork;
     magma_dmalloc_cpu(&work, lwork); magma_imalloc_cpu(&iwork, liwork);
     magma_dsyevd( jobz, uplo, N, A.base_array, N, W.base_array, work, lwork, iwork, liwork, &info );
     magma_free_cpu(work); magma_free_cpu(iwork);

     if(info!=0) {cout<<"Dsyevd failed: info= "<< info<<"\n"; exit(1);}
 }


 /*************************************/
 /*Diagonalize Hermitian Matrix********/
 /*************************************/
 void eigen_magma(Matrix<complex<double>,2>& A, Matrix<double,1>& W, char JOBZ, char UPLO)
 {
     if(A.L1!=A.L2) {cout<<"Input for eigen is not square matrix!\n"; exit(1);}

     magma_vec_t jobz = magma_vec_const(JOBZ); magma_uplo_t uplo = magma_uplo_const(UPLO);
     magma_int_t N=A.L1; magma_int_t info;

     magmaDoubleComplex work_test[1]; double rwork_test[1]; magma_int_t iwork_test[1];
     magma_int_t lwork=-1; magma_int_t lrwork=-1; magma_int_t liwork=-1;
     magma_zheevd( jobz, uplo, N, (magmaDoubleComplex* ) A.base_array, N, W.base_array, 
                   work_test, lwork, rwork_test, lrwork, iwork_test, liwork, &info );

     lwork=lround( MAGMA_Z_REAL(work_test[0]) ); lrwork=lround(rwork_test[0]); liwork=iwork_test[0];
     magmaDoubleComplex* work; double* rwork; magma_int_t* iwork;
     magma_zmalloc_cpu(&work, lwork); magma_dmalloc_cpu(&rwork, lrwork); magma_imalloc_cpu(&iwork, liwork);
     magma_zheevd( jobz, uplo, N, (magmaDoubleComplex* ) A.base_array, N, W.base_array,
                   work, lwork, rwork, lrwork, iwork, liwork, &info );

     magma_free_cpu(work); magma_free_cpu(rwork); magma_free_cpu(iwork);
     if(info!=0) {cout<<"Zheevd failed: info= "<< info<<"\n"; exit(1);}
 }


 /******************************************/
 /*LU Decomposition a complex square Matrix*/
 /******************************************/
 LUDecomp<complex<double>> LUconstruct_magma(const Matrix<complex<double>,2>& x)
 {
     if(x.L1!=x.L2) {cout<<"Input for LU is not square matrix!\n"; exit(1);}
     LUDecomp<complex<double>> y; y.A=x; y.ipiv=Matrix<BL_INT,1>(x.L1);
     magma_int_t N=x.L1; magma_int_t info;
     magma_zgetrf(N, N, (magmaDoubleComplex* )y.A.base_array, N, (magma_int_t*)y.ipiv.base_array, &info);
     y.info=info;
     if(y.info<0) {cout<<"The "<<y.info<<"-th parameter is illegal!\n"; exit(1);}
     return y;
 }


 /*********************************************************************************************************************/
 /*Inverse of  Matrix: If determinant of the matrix is outof machine precision, inverse should be fine, since it solve*
  *The linear equation, every small value is well defined                                                             */
 /*********************************************************************************************************************/
 Matrix<complex<double>,2> inverse_magma(const LUDecomp<complex<double>>& x)
 {
     magma_int_t N=x.A.L1; magma_int_t info;

     magmaDoubleComplex_ptr d_A , dwork;
     magma_int_t lda, ldwork;
     lda = ((N+31)/32)*32;             //round up to multiple of 32 for best GPU performance
     ldwork = N*magma_get_zgetri_nb(N); // magma_get_zgetri_nb optimizes the blocksize
     magma_zmalloc( &d_A, lda*N ); magma_zmalloc( &dwork, ldwork );

     //copy matrix from CPU to GPU
     magma_zsetmatrix( N, N, (magmaDoubleComplex* )x.A.base_array, N, d_A, lda );

     //calculate the inverse matrix with zgetri
     magma_zgetri_gpu( N, d_A, lda, (magma_int_t*) x.ipiv.base_array, dwork, ldwork, &info );
     if(info<0) {cout<<"The "<<info<<"-th parameter is illegal!\n"; exit(1);}
     
     //copy matrix from GPU to CPU
     Matrix<complex<double>,2> A(N,N); 
     magma_zgetmatrix( N, N, d_A, lda, (magmaDoubleComplex* )A.base_array, N );

     magma_free(d_A); magma_free(dwork);

     return A;
 }


 /******************************************************/
 /*Solve Linear Equation of the matrix A*M=B: M=A^{-1}B*/
 /******************************************************/
 Matrix<complex<double>,2> solve_lineq_magma(const LUDecomp<complex<double>>& x, const Matrix<complex<double>,2>& B, char TRANS)
 {
     if(x.A.L1!=B.L1)  {cout<<"Input size for solving linear equation is not consistent!\n"; exit(1);}
     magma_int_t N=B.L1; magma_int_t NRHS=B.L2; magma_int_t info;

     magma_trans_t Trans = magma_trans_const(TRANS);
     magmaDoubleComplex_ptr d_A, d_B;
     magma_int_t lda, ldb;
     lda = ((N+31)/32)*32;
     ldb = ((N+31)/32)*32;

     //allocate memory on GPU
     magma_zmalloc( &d_A, lda*N );
     magma_zmalloc( &d_B, ldb*NRHS );

     //copy matrix from CPU to GPU
     magma_zsetmatrix( N, N,    (magmaDoubleComplex* )x.A.base_array, N, d_A, lda );
     magma_zsetmatrix( N, NRHS, (magmaDoubleComplex* )B.base_array,   N, d_B, ldb );

     //Solve the equation
     magma_zgetrs_gpu( Trans, N, NRHS, d_A, lda, (magma_int_t*)x.ipiv.base_array, d_B, ldb, &info );
     if(info!=0)
     {
         cout<<"Solve linear equation is not suceesful: "<<info<<"-th parameter is illegal! \n";
         exit(1);
     }


     //copy matrix from GPU to CPU
     Matrix<complex<double>,2> M(N,NRHS);
     magma_zgetmatrix( N, NRHS, d_B, ldb, (magmaDoubleComplex* ) M.base_array, N );

     //free memory
     magma_free( d_A );  magma_free( d_B );     

     return M;
 }


 /******************************/
 /*QR decompostion of matrix ph*/
 /******************************/


//This QR does not give correct det(R) in large system size. 
//Though it works well for small system size. Need to use the
//CPU interface

/*
 double QRMatrix_magma(Matrix<complex<double>,2>& ph)
 {
     //If we need to call magma two or more times, it's better 
     //to use GPU interface, this will avoid transfer between
     //CPU and GPU.
     magma_int_t L=ph.L1; magma_int_t N=ph.L2; magma_int_t info;

     //Prepare for zgeqrf_gpu
     magmaDoubleComplex_ptr d_A;
     magma_int_t LDA = ((L+31)/32)*32;
     magma_zmalloc(&d_A, LDA*N );
     magma_zsetmatrix( L, N, (magmaDoubleComplex*) ph.base_array, L, d_A, LDA );

     magmaDoubleComplex* tau;  magma_zmalloc_cpu( &tau, N );

     magmaDoubleComplex_ptr dT;
     magma_int_t nb=magma_get_zgeqrf_nb(L);
     magma_zmalloc( &dT, (2*N + ((N+31)/32)*32 )*nb );

     //QR with zgeqrf_gpu
     magma_zgeqrf_gpu(L,N,d_A,LDA,tau,dT,&info);
     if(info!=0) {cout<<"QR run is not suceesful: "<<info<<"-th parameter is illegal! \n"; exit(1);}

     //calculate det
     magma_zgetmatrix(L, N, d_A, LDA, (magmaDoubleComplex* ) ph.base_array, L);
     complex<double> det={1.0,0.0}; for (size_t i=0; i<ph.L2; i++)  det*=ph(i,i);

     //zungqr_gpu get the correct form ph 
     magma_zungqr_gpu(L,N,N,d_A,LDA,tau,dT,nb,&info);
     if(info!=0) {cout<<"magma_zungqr_gpu run is not suceesful: "<<info<<"-th parameter is illegal! \n"; exit(1);}
     magma_zgetmatrix(L, N, d_A, LDA, (magmaDoubleComplex* ) ph.base_array, L);

     //Reshape the phi to get positive det
     if(det.real()<0) {det=-det; for(size_t i=0; i<ph.L1; i++) ph(i,0)=-ph(i,0);}

     magma_free(d_A); magma_free_cpu(tau); magma_free(dT);

     return det.real();
 }
*/

 double QRMatrix_magma(Matrix<complex<double>,2>& ph)
 {
     magma_int_t L=ph.L1; magma_int_t N=ph.L2; magma_int_t info;

     magmaDoubleComplex* tau;  magma_zmalloc_cpu( &tau, N );

     magmaDoubleComplex work_test[1]; magma_int_t lwork=-1;
     magma_zgeqrf(L, N, (magmaDoubleComplex *)ph.base_array, L, tau, work_test, lwork, &info);

     lwork=lround( MAGMA_Z_REAL(work_test[0]) );
     magmaDoubleComplex* work;  magma_zmalloc_cpu( &work, lwork );
     magma_zgeqrf(L, N, (magmaDoubleComplex *)ph.base_array, L, tau, work, lwork, &info);
     if(info!=0) {cout<<"QR run is not suceesful: "<<info<<"-th parameter is illegal! \n"; exit(1);}

     complex<double> det={1.0,0.0}; for (size_t i=0; i<ph.L2; i++)  det*=ph(i,i);
     magma_zungqr2(L, N, N, (magmaDoubleComplex *)ph.base_array, L, tau, &info );
     if(info!=0) {cout<<"magma_zungqr2 run is not suceesful: "<<info<<"-th parameter is illegal! \n"; exit(1);}

     //Reshape the phi to get positive det
     if(det.real()<0) {det=-det; for(size_t i=0; i<ph.L1; i++) ph(i,0)=-ph(i,0);}

     magma_free_cpu(tau); magma_free_cpu(work);

     return det.real();
 }


} //end namespace matrix_hao_lib

#endif
