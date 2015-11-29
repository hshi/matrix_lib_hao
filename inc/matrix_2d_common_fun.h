#ifndef MATRIX_2D_COMMON_FUN_H
#define MATRIX_2D_COMMON_FUN_H


namespace matrix_hao_lib
{

 /*************************************/
 /*Diagonalize Hermitian Matrix********/
 /*************************************/
 void check_Hermitian(const Matrix<std::complex<double>,2>& A);



 /*******************************************/
 /*LU Decomposition of Complex double Matrix*/
 /*******************************************/
 template <class T> class LUDecomp
 {
     public:
     Matrix<T,2> A;
     Matrix<BL_INT,1> ipiv;
     BL_INT info;
   
     LUDecomp() {}
     LUDecomp(const LUDecomp<T>& x) {A=x.A;ipiv=x.ipiv;info=x.info;}
     LUDecomp(LUDecomp<T>&& x) {A=std::move(x.A);ipiv=std::move(x.ipiv);info=x.info;}
     ~LUDecomp() {}
     LUDecomp<T>& operator = (const LUDecomp<T>& x) {A=x.A;ipiv=x.ipiv;info=x.info;return *this;}
     LUDecomp<T>& operator = (LUDecomp<T>&& x) {A=std::move(x.A);ipiv=std::move(x.ipiv);info=x.info;return *this;}
 };


 /************************/
 /*Determinant of  Matrix*/
 /************************/
 std::complex<double> determinant(const LUDecomp<std::complex<double>>& x);
 void lognorm_phase_determinant(const LUDecomp<std::complex<double>>& x, std::complex<double>& lognorm, std::complex<double>& phase);
 std::complex<double> log_determinant(const LUDecomp<std::complex<double>>& x);

 /*******************************/
 /*Diagonal array multipy matrix*/
 /*******************************/
 Matrix<std::complex<double>,2> D_Multi_Matrix(const Matrix<std::complex<double>,1>& D,const Matrix<std::complex<double>,2>& ph);

  /**********************/
 /*Pfaffian of a matrix*/
 /**********************/
 void check_skew_symmetric(const Matrix<std::complex<double>,2>& A);
 std::complex<double> Pfaffian(Matrix<std::complex<double>,2>& A);
 
}//end namespace matrix_hao_lib

#endif
