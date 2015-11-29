#include "matrix_all.h"

using namespace std;

namespace matrix_hao_lib
{

 void LUDecomp_test()
 {
     Matrix<complex<double>,2> X={3,3,{ {1.0,0.0} ,   {3.0,4.0},    {2.123,3.11}, 
                                        {3.0,-2.0},   {2.0,0.0},    {5.123,3.11}, 
                                        {2.123,-5.11},{5.123,-6.11},{3,0.0} } };
     LUDecomp<complex<double>> LU=LUconstruct_cpu(X);

     Matrix<complex<double>,2> A_exact=LU.A;

     size_t flag=0;

     LUDecomp<complex<double>> LUC(LU);
     flag+=diff(LUC.A,A_exact,1e-13);

     LUDecomp<complex<double>> LUR(std::move(LU));
     flag+=diff(LUR.A,A_exact,1e-13);

     LUDecomp<complex<double>> LUEC;LUEC=LUC;
     flag+=diff(LUEC.A,A_exact,1e-13);

     LUDecomp<complex<double>> LUER;LUER=std::move(LUR);
     flag+=diff(LUER.A,A_exact,1e-13);


     if(flag==0) cout<<"LUDecomp passed complex double test! \n";
     else cout<<"WARNING!!!!!!!!! LUDecomp failed complex double test! \n";
 }


 void determinant_test()
 {
     Matrix<complex<double>,2> X={3,3,{ {1.0,0.0} ,   {3.0,4.0},    {2.123,3.11},
                                        {3.0,-2.0},   {2.0,0.0},    {5.123,3.11},
                                        {2.123,-5.11},{5.123,-6.11},{3,0.0} } };
     complex<double> det=determinant( LUconstruct_cpu(X) );
     complex<double> det_exact={123.11968700000003,3.3324580000000115};
     if(abs(det-det_exact)<1e-12) cout<<"Determinant passed complex double test in cpu! \n";
     else cout<<"WARNING!!!!!!!!! Determinant failed complex double test in cpu! \n";
#ifdef USE_MAGMA
     det=determinant( LUconstruct_magma(X) );
     if(abs(det-det_exact)<1e-12) cout<<"Determinant passed complex double test in magma! \n";
     else cout<<"WARNING!!!!!!!!! Determinant failed complex double test in magma! \n";

#endif
     //cout<<setprecision(16);
     //cout<<det<<"\n";
 }

 void log_determinant_test()
 {
     Matrix<complex<double>,2> X={3,3,{ {1.0,0.0} ,   {3.0,4.0},    {2.123,3.11},
                                        {3.0,-2.0},   {2.0,0.0},    {5.123,3.11},
                                        {2.123,-5.11},{5.123,-6.11},{3,0.0} } };
     X*=1.e103;
     complex<double> logdet=log_determinant( LUconstruct_cpu(X) );
     complex<double> logdet_exact={716.3123168546207,0.027060209772387683};
     if(abs(logdet-logdet_exact)<1e-12) cout<<"Log_determinant passed complex double test in cpu! \n";
     else cout<<"WARNING!!!!!!!!! Log_determinant failed complex double test in cpu! \n";
#ifdef USE_MAGMA
     logdet=log_determinant( LUconstruct_magma(X) );
     if(abs(logdet-logdet_exact)<1e-12) cout<<"Log_determinant passed complex double test in magma! \n";
     else cout<<"WARNING!!!!!!!!! Log_determinant failed complex double test in magma! \n";
#endif
     //cout<<abs(logdet-logdet_exact)<<"\n";
     //cout<<setprecision(16);
     //cout<<logdet<<"\n";
 }


 void D_Multi_Matrix_test()
 {
     Matrix<complex<double>,2> A={3,2,{ {2.0,0.0} ,   {3.0,5.0},    {3.123,3.11},
                                        {3.0,-6.0},   {2.0,1.0},    {6.123,3.11},} };
     Matrix<complex<double>,1> D={3, { {1.2,0.0},{2.0,0.0},{3.0,0.0} } };
     Matrix<complex<double>,2> B=D_Multi_Matrix(D,A);
     Matrix<complex<double>,2> B_exact={3,2,{ {2.4,0.0} ,   {6.0,10.0},    {9.369,9.33},
                                              {3.6,-7.2},   {4.0,2.0 },    {18.369,9.33},} };
     size_t flag=diff(B,B_exact,1e-12);
     if(flag==0) cout<<"D_Multi_Matrix passed complex double test! \n";
     else cout<<"WARNING!!!!!!!!! D_Multi_Matrix failed complex double test! \n"; 
 }


 void pfaffian_test()
 {
     Matrix<complex<double>,2> A={4,4,{ {0.0,0.0},  {1.0,2.0},  {1.5,0.0}, {2.3,0.0}, 
                                        {-1.0,-2.0},{0.0,0.0},  {-3.0,0.0},{1.5,0.0}, 
                                        {-1.5,0.0}, {3.0,0.0},  {0.0,0.0}, {-2.5,-5.3},
                                        {-2.3,0.0}, {-1.5,0.0}, {2.5,5.3}, {0.0,0.0} } };
     check_skew_symmetric(A);
     complex<double> pf=Pfaffian(A);
     complex<double> exact{-1.05,-10.3};
     if( abs(pf-exact)< 1e-13 ) cout<<"Pfaffian passed complex double test! \n";
     else cout<<"WARNING!!!!!!!!! Pfaffian failed complex double test! \n";
 }

 void matrix_2d_common_fun_test()
 {
     LUDecomp_test();
     determinant_test();
     log_determinant_test();
     D_Multi_Matrix_test();
     pfaffian_test();
 }

} //end namespace matrix_hao_lib
