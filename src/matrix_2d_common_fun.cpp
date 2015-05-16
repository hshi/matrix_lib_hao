#include "matrix_all.h"

using namespace std;

namespace matrix_hao_lib
{

 /*************************************/
 /*Check Hermitian of the matrix*******/
 /*************************************/
 void check_Hermitian(const Matrix<complex<double>,2>& A)
 {
     if(A.L1!=A.L2) {cout<<"Input for Hermitian is not square matrix! \n"; exit(1);}
     double error=0; double norm=0;
     for(size_t i=0; i<A.L1; i++)
     {
         for(size_t j=i; j<A.L2; j++)
         {
             error+=abs(A(i,j)-conj(A(j,i)));
             norm+=abs(A(i,j));
         }
     }
     norm/=(A.L_f()*1.0);
     if(error/norm>1e-12) cout<<"Warning!!!!!Matrix is not Hermition!";
 }


 /************************/
 /*Determinant of  Matrix*/
 /************************/
 complex<double> determinant(const LUDecomp<complex<double>>& x)
 {
     if(x.info>0) return 0;
    
     complex<double> det={1,0};
     BL_INT L=x.ipiv.L1;
     for(BL_INT i=0;i<L;i++)
     {
         if(x.ipiv(i)!=(i+1)) det*=(-x.A(i,i));
         else det*=x.A(i,i);
     }
     return det;
 }

 /*****************************************/
 /*Get Log(|det|) and det/|det| of  Matrix*/
 /*****************************************/
 void lognorm_phase_determinant(const LUDecomp<complex<double>>& x, complex<double>& lognorm, complex<double>& phase)
 {
     if(x.info>0)
     {
         cout<<"WARNING!!!! lognorm_phase_determinant function has zero determinant!\n";
         lognorm=complex<double>(-1e300,0.0);
         phase=complex<double>(1.0,0.0);
         return;
     }

     lognorm=complex<double>(0.0,0.0); phase=complex<double>(1.0,0.0);
     BL_INT L=x.ipiv.L1;
     for(BL_INT i=0;i<L;i++)
     {
         lognorm+=log(abs(x.A(i,i)));
         if(x.ipiv(i)!=(i+1)) phase*=(-x.A(i,i)/abs(x.A(i,i)));
         else phase*=(x.A(i,i)/abs(x.A(i,i)));
     }
     return;
 }


 /****************************/
 /*Log Determinant of  Matrix*/
 /****************************/
 complex<double> log_determinant(const LUDecomp<complex<double>>& x)
 {
     complex<double> log_det,phase; 
     lognorm_phase_determinant(x,log_det,phase);
     log_det+=log(phase);
     return log_det;
 }


 /*******************************/
 /*Diagonal array multipy matrix*/
 /*******************************/
 Matrix<complex<double>,2> D_Multi_Matrix(const Matrix<complex<double>,1>& D,const Matrix<complex<double>,2>& ph)
 {
     if(D.L1!=ph.L1) {cout<<"D_Multi_Matrix input error: D.L1!=ph.L1! \n"; exit(1);} 
     Matrix<complex<double>,2> ph_new(ph.L1,ph.L2);
     for(size_t i=0; i<ph.L1; i++)
     {
         for(size_t j=0; j<ph.L2; j++) ph_new(i,j)=D(i)*ph(i,j);
     }
     return ph_new; 
 }
 
} //end namespace matrix_hao_lib
