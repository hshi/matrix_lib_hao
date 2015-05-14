#ifndef MATRIX_HAO_ELEMENT_WISE
#define MATRIX_HAO_ELEMENT_WISE

#include <cmath>


namespace matrix_hao_lib
{

 //for add: (array+array)
 template <class T, int D>
 Matrix<T,D> operator + (const Matrix<T,D>& A,const Matrix<T,D>& B) {Matrix<T,D> C;C=A; C+=B; return C;} 
 template <class T, int D>
 Matrix<T,D> operator + (const Matrix<T,D>& A, Matrix<T,D>&& B)     {Matrix<T,D> C;C=std::move(B); C+=A; return C;} 
 template <class T, int D>
 Matrix<T,D> operator + (Matrix<T,D>&& A,const Matrix<T,D>& B)      {Matrix<T,D> C;C=std::move(A); C+=B; return C;} 
 template <class T, int D>
 Matrix<T,D> operator + (Matrix<T,D>&& A,Matrix<T,D>&& B)           {Matrix<T,D> C;C=std::move(A); C+=B; return C;} 


 //for minus:(array-array)
 template <class T, int D>
 Matrix<T,D> operator - (const Matrix<T,D>& A,const Matrix<T,D>& B) {Matrix<T,D> C;C=A; C-=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator - (const Matrix<T,D>& A, Matrix<T,D>&& B){Matrix<T,D> C;C=std::move(B); C.min_add_equal(A); return C;} 
 template <class T, int D>
 Matrix<T,D> operator - (Matrix<T,D>&& A,const Matrix<T,D>& B)      {Matrix<T,D> C;C=std::move(A); C-=B; return C;} 
 template <class T, int D>
 Matrix<T,D> operator - (Matrix<T,D>&& A,Matrix<T,D>&& B)           {Matrix<T,D> C;C=std::move(A); C-=B; return C;} 


 //for time: (array*array)
 template <class T, int D>
 Matrix<T,D> operator * (const Matrix<T,D>& A,const Matrix<T,D>& B) {Matrix<T,D> C;C=A; C*=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator * (const Matrix<T,D>& A, Matrix<T,D>&& B)     {Matrix<T,D> C;C=std::move(B); C*=A; return C;}
 template <class T, int D>
 Matrix<T,D> operator * (Matrix<T,D>&& A,const Matrix<T,D>& B)      {Matrix<T,D> C;C=std::move(A); C*=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator * (Matrix<T,D>&& A,Matrix<T,D>&& B)           {Matrix<T,D> C;C=std::move(A); C*=B; return C;}



 //for div:(array/array)
 template <class T, int D>
 Matrix<T,D> operator / (const Matrix<T,D>& A,const Matrix<T,D>& B) {Matrix<T,D> C;C=A; C/=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator / (const Matrix<T,D>& A, Matrix<T,D>&& B){Matrix<T,D> C;C=std::move(B); C.inv_div_equal(A); return C;}
 template <class T, int D>
 Matrix<T,D> operator / (Matrix<T,D>&& A,const Matrix<T,D>& B)      {Matrix<T,D> C;C=std::move(A); C/=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator / (Matrix<T,D>&& A,Matrix<T,D>&& B)           {Matrix<T,D> C;C=std::move(A); C/=B; return C;}


 //for add: (array+scalar)
 template <class T, int D>
 Matrix<T,D> operator + (const Matrix<T,D>& A,const T & B) {Matrix<T,D> C;C=A; C+=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator + (Matrix<T,D>&& A, const T & B)     {Matrix<T,D> C;C=std::move(A); C+=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator + (const T & B,const Matrix<T,D>& A) {Matrix<T,D> C;C=A; C+=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator + (const T & B,Matrix<T,D>&& A)      {Matrix<T,D> C;C=std::move(A); C+=B; return C;}


 //for minus: (array-scalar)
 template <class T, int D>
 Matrix<T,D> operator - (const Matrix<T,D>& A,const T & B) {Matrix<T,D> C;C=A; C-=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator - (Matrix<T,D>&& A, const T & B)     {Matrix<T,D> C;C=std::move(A); C-=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator - (const T & B,const Matrix<T,D>& A) {Matrix<T,D> C;C=A; C.min_add_equal(B);return C;}
 template <class T, int D>
 Matrix<T,D> operator - (const T & B,Matrix<T,D>&& A)      {Matrix<T,D> C;C=std::move(A);C.min_add_equal(B);return C;}



 //for time: (array*scalar)
 template <class T, int D>
 Matrix<T,D> operator * (const Matrix<T,D>& A,const T & B) {Matrix<T,D> C;C=A; C*=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator * (Matrix<T,D>&& A, const T & B)     {Matrix<T,D> C;C=std::move(A); C*=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator * (const T & B,const Matrix<T,D>& A) {Matrix<T,D> C;C=A; C*=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator * (const T & B,Matrix<T,D>&& A)      {Matrix<T,D> C;C=std::move(A); C*=B; return C;}


 //for div: (array/scalar)
 template <class T, int D>
 Matrix<T,D> operator / (const Matrix<T,D>& A,const T & B) {Matrix<T,D> C;C=A; C/=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator / (Matrix<T,D>&& A, const T & B)     {Matrix<T,D> C;C=std::move(A); C/=B; return C;}
 template <class T, int D>
 Matrix<T,D> operator / (const T & B,const Matrix<T,D>& A) {Matrix<T,D> C;C=A; C.inv_div_equal(B);return C;}
 template <class T, int D>
 Matrix<T,D> operator / (const T & B,Matrix<T,D>&& A)      {Matrix<T,D> C;C=std::move(A);C.inv_div_equal(B);return C;}

 //for diff: return the different elements in two matrix
 template <class T, int D>
 size_t diff(const Matrix<T,D>& A, const Matrix<T,D>& B, double eta)
 {
    //Use std::abs instead of abs is very important here since it is in header file
    //For normal cpp file, if we use " using namespace std; ", we are using std:abs
    //To see the difference, run the following line:
    //std::cout<<abs(0.123)<<" "<<std::abs(0.123)<<std::endl;

    size_t flag=0; double abs_eta=std::abs(eta);
    for(size_t i=0; i<A.L_f(); i++) {if( std::abs( A.base_array[i]-B.base_array[i] )> abs_eta  ) flag++;}
    return flag;
 }

} //end namespace matrix_hao_lib

#endif
