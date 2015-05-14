#ifdef MPI_HAO

#include "matrix_all.h"

using namespace std;

namespace matrix_hao_lib
{

 void MPIBcast_double_one_test()
 {
     Matrix<double,1> A={4,{14.861,12.630129,20.984,23.753129}};
     Matrix<double,1> A_exact(4); A_exact=0.0;
     int rank;MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     if(rank==0) A=0.0;
     MPIBcast(A);
     size_t flag=diff(A,A_exact,1e-12);
     if(flag!=0) cout<<"Warning!!!!Bcast failed the double 1d test! rank: "<<rank<<endl;
 }


 void MPIBcast_double_two_test()
 {
     Matrix<double,2> A={2,2,{14.861,12.630129,20.984,23.753129}};
     Matrix<double,2> A_exact(2,2); A_exact=0.0;
     int rank;MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     if(rank==0) A=0.0;
     MPIBcast(A);
     size_t flag=diff(A,A_exact,1e-12);
     if(flag!=0) cout<<"Warning!!!!Bcast failed the double 2d test! rank: "<<rank<<endl;
 }


 void MPIBcast_complex_double_two_test()
 {
     Matrix<complex<double>,2> A={2,2, 
                                      { {-13.769,40.877}, {-16.551971,38.73806},
                                        {-17.756,56.71},  {-22.838971, 66.77106} }
                                 };
     Matrix<complex<double>,2> A_exact(2,2);A_exact=complex<double>(0,0);
     int rank;MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     if(rank==0) A=complex<double>(0,0);
     MPIBcast(A);
     size_t flag=diff(A,A_exact,1e-12);
     if(flag!=0) cout<<"Warning!!!!Bcast failed the complex double 2d test! rank: "<<rank<<endl;
 }

 void MPIBcast_complex_double_three_test()
 {
     Matrix<complex<double>,3> A={2,3,2,{
                                  complex<double>(1.0,2.3),complex<double>(3.0,4.0),complex<double>(2.123,3.11),
                                  complex<double>(2.0,3.3),complex<double>(4.0,5.0),complex<double>(3.123,4.11),
                                  complex<double>(3.0,2.3),complex<double>(3.0,4.0),complex<double>(2.123,3.11),
                                  complex<double>(2.0,3.3),complex<double>(4.0,5.0),complex<double>(3.123,6.11)
                                      }};
     Matrix<complex<double>,3> A_exact(2,3,2);A_exact=complex<double>(0,0);
     int rank;MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     if(rank==0) A=complex<double>(0,0);
     MPIBcast(A);
     size_t flag=diff(A,A_exact,1e-12);
     if(flag!=0) cout<<"Warning!!!!Bcast failed the complex double 3d test! rank: "<<rank<<endl;
 }


 void matrix_mpi_test()
 {
     int rank;MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     if(rank==0) cout<<"Testing Matrix_hao_mpi version......\n"<<endl;

     MPIBcast_double_one_test();
     MPIBcast_double_two_test();
     MPIBcast_complex_double_two_test();
     MPIBcast_complex_double_three_test();

     MPI_Barrier(MPI_COMM_WORLD);
     if(rank==0) cout<<"\n\nIf these is no warning, we have passed all the test!"<<endl; 
 }

}

#endif
