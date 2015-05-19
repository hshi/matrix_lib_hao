#ifdef USE_MAGMA

#include "matrix_all.h"
#include "random_fill_matrix.h"

using namespace std;
namespace matrix_hao_lib
{
   void gmm_float_timing(int M, int N, int K)
   {
       real_Double_t cpu_time, magma_time;
       Matrix<float,2> a(M, K); fill_random(a);
       Matrix<float,2> b(K, N); fill_random(b);
       Matrix<float,2> c_cpu(M, N), c_magma(M, N);

       cpu_time = magma_wtime();
       gmm_cpu(a,b,c_cpu);
       cpu_time = magma_wtime() - cpu_time;

       magma_time = magma_wtime();
       gmm_magma(a,b,c_magma);
       magma_time = magma_wtime() - magma_time;

       size_t flag=diff(c_cpu,c_magma,1e-5);
       cout<<setw(16)<<M<<setw(16)<<N<<setw(16)<<K<<setw(16)<<cpu_time<<setw(16)<<magma_time<<setw(16)<<flag<<endl;
   }

   void gmm_float_timing_loop()
   {
       cout<<"Timing gmm float:"<<endl;
       cout<<setw(16)<<"M"<<setw(16)<<"N"<<setw(16)<<"K"<<setw(16)<<"cpu_time"<<setw(16)<<"magma_time"<<setw(16)<<"flag"<<endl;
       for (int i = 8; i <= 1087; i *= 2)       gmm_float_timing(i, i, i);
       for (int i = 1088; i <= 5184; i += 1024) gmm_float_timing(i, i, i);
       cout<<"\n\n\n"<<endl;
   }


   void gmm_double_timing(int M, int N, int K)
   {
       real_Double_t cpu_time, magma_time;
       Matrix<double,2> a(M, K); fill_random(a);
       Matrix<double,2> b(K, N); fill_random(b);
       Matrix<double,2> c_cpu(M, N), c_magma(M, N);

       cpu_time = magma_wtime();
       gmm_cpu(a,b,c_cpu);
       cpu_time = magma_wtime() - cpu_time;

       magma_time = magma_wtime();
       gmm_magma(a,b,c_magma);
       magma_time = magma_wtime() - magma_time;

       size_t flag=diff(c_cpu,c_magma,1e-10);
       cout<<setw(16)<<M<<setw(16)<<N<<setw(16)<<K<<setw(16)<<cpu_time<<setw(16)<<magma_time<<setw(16)<<flag<<endl;
   }

   void gmm_double_timing_loop()
   {
       cout<<"Timing gmm double:"<<endl;
       cout<<setw(16)<<"M"<<setw(16)<<"N"<<setw(16)<<"K"<<setw(16)<<"cpu_time"<<setw(16)<<"magma_time"<<setw(16)<<"flag"<<endl;
       for (int i = 8; i <= 1087; i *= 2)       gmm_double_timing(i, i, i);
       for (int i = 1088; i <= 5184; i += 1024) gmm_double_timing(i, i, i);
       cout<<"\n\n\n"<<endl;
   }



   void gmm_complexfloat_timing(int M, int N, int K)
   {
       real_Double_t cpu_time, magma_time;
       Matrix<complex<float>,2> a(M, K); fill_random(a);
       Matrix<complex<float>,2> b(K, N); fill_random(b);
       Matrix<complex<float>,2> c_cpu(M, N), c_magma(M, N);

       cpu_time = magma_wtime();
       gmm_cpu(a,b,c_cpu);
       cpu_time = magma_wtime() - cpu_time;

       magma_time = magma_wtime();
       gmm_magma(a,b,c_magma);
       magma_time = magma_wtime() - magma_time;

       size_t flag=diff(c_cpu,c_magma,1e-5);
       cout<<setw(16)<<M<<setw(16)<<N<<setw(16)<<K<<setw(16)<<cpu_time<<setw(16)<<magma_time<<setw(16)<<flag<<endl;
   }

   void gmm_complexfloat_timing_loop()
   {
       cout<<"Timing gmm complexfloat:"<<endl;
       cout<<setw(16)<<"M"<<setw(16)<<"N"<<setw(16)<<"K"<<setw(16)<<"cpu_time"<<setw(16)<<"magma_time"<<setw(16)<<"flag"<<endl;
       for (int i = 8; i <= 1087; i *= 2)       gmm_complexfloat_timing(i, i, i);
       for (int i = 1088; i <= 5184; i += 1024) gmm_complexfloat_timing(i, i, i);
       cout<<"\n\n\n"<<endl;
   }



   void gmm_complexdouble_timing(int M, int N, int K)
   {
       real_Double_t cpu_time, magma_time;
       Matrix<complex<double>,2> a(M, K); fill_random(a);
       Matrix<complex<double>,2> b(K, N); fill_random(b);
       Matrix<complex<double>,2> c_cpu(M, N), c_magma(M, N);

       cpu_time = magma_wtime();
       gmm_cpu(a,b,c_cpu);
       cpu_time = magma_wtime() - cpu_time;

       magma_time = magma_wtime();
       gmm_magma(a,b,c_magma);
       magma_time = magma_wtime() - magma_time;

       size_t flag=diff(c_cpu,c_magma,1e-10);
       cout<<setw(16)<<M<<setw(16)<<N<<setw(16)<<K<<setw(16)<<cpu_time<<setw(16)<<magma_time<<setw(16)<<flag<<endl;
   }

   void gmm_complexdouble_timing_loop()
   {
       cout<<"Timing gmm complexdouble:"<<endl;
       cout<<setw(16)<<"M"<<setw(16)<<"N"<<setw(16)<<"K"<<setw(16)<<"cpu_time"<<setw(16)<<"magma_time"<<setw(16)<<"flag"<<endl;
       for (int i = 8; i <= 1087; i *= 2)       gmm_complexdouble_timing(i, i, i);
       for (int i = 1088; i <= 5184; i += 1024) gmm_complexdouble_timing(i, i, i);
       cout<<"\n\n\n"<<endl;
   }


   void eigen_double_timing(int N)
   {
       real_Double_t cpu_time, magma_time;
       Matrix<double,2>  a_cpu(N,N), a_magma(N,N); 
       Matrix<double,1>  w_cpu(N),   w_magma(N);

       //Get a real symmetry a_cpu and a_magma
       fill_random(a_cpu);
       for (int i=0; i<N; i++)  {for (int j=i+1; j<N; j++) a_cpu(i,j)=a_cpu(j,i);}
       a_magma=a_cpu;

       cpu_time = magma_wtime();
       eigen_cpu(a_cpu,w_cpu);
       cpu_time = magma_wtime() - cpu_time;

       magma_time = magma_wtime();
       eigen_magma(a_magma,w_magma);
       magma_time = magma_wtime() - magma_time;

       size_t flag=diff(w_cpu,w_magma,1e-10);
       cout<<setw(16)<<N<<setw(16)<<cpu_time<<setw(16)<<magma_time<<setw(16)<<flag<<endl;
   }

   void eigen_double_timing_loop()
   {
       cout<<"Timing eigen double:"<<endl;
       cout<<setw(16)<<"N"<<setw(16)<<"cpu_time"<<setw(16)<<"magma_time"<<setw(16)<<"flag"<<endl;
       for (int i = 210; i <= 1000; i += 200)     eigen_double_timing(i);
       for (int i = 1088; i <= 7232; i += 1024)   eigen_double_timing(i);
       cout<<"\n\n\n"<<endl;
   }


   void eigen_complexdouble_timing(int N)
   {
       real_Double_t cpu_time, magma_time;
       Matrix<complex<double>,2>  a_cpu(N,N), a_magma(N,N);
       Matrix<double,1>           w_cpu(N),   w_magma(N);

       //Get a Hermition matrix a_cpu and a_magma
       fill_random(a_cpu); 
       for (int i=0; i<N; i++)
       {
           a_cpu(i,i)=a_cpu(i,i).real();
           for (int j=i+1; j<N; j++) a_cpu(i,j)=conj(a_cpu(j,i));
       }
       a_magma=a_cpu;
       check_Hermitian(a_cpu); 

       cpu_time = magma_wtime();
       eigen_cpu(a_cpu,w_cpu);
       cpu_time = magma_wtime() - cpu_time;

       magma_time = magma_wtime();
       eigen_magma(a_magma,w_magma);
       magma_time = magma_wtime() - magma_time;

       size_t flag=diff(w_cpu,w_magma,1e-10);
       cout<<setw(16)<<N<<setw(16)<<cpu_time<<setw(16)<<magma_time<<setw(16)<<flag<<endl;
       //if(N==810) cout<<w_cpu-w_magma<<endl;
   }

   void eigen_complexdouble_timing_loop()
   {
       cout<<"Timing eigen complex double:"<<endl;
       cout<<setw(16)<<"N"<<setw(16)<<"cpu_time"<<setw(16)<<"magma_time"<<setw(16)<<"flag"<<endl;
       for (int i = 210; i <= 1000; i += 200)     eigen_complexdouble_timing(i);
       for (int i = 1088; i <= 7232; i += 1024)   eigen_complexdouble_timing(i);
       cout<<"\n\n\n"<<endl;
   }

   void LUconstruct_timing(int N)
   {
       real_Double_t cpu_time, magma_time;
       Matrix<complex<double>,2>  X(N,N); fill_random(X);

       cpu_time = magma_wtime();
       LUDecomp<complex<double>> LU_cpu=LUconstruct_cpu(X);
       cpu_time = magma_wtime() - cpu_time;      

       magma_time = magma_wtime();
       LUDecomp<complex<double>> LU_magma=LUconstruct_magma(X);
       magma_time = magma_wtime() - magma_time;

       size_t flag=diff(LU_cpu.A, LU_magma.A,1e-10);
       cout<<setw(16)<<N<<setw(16)<<cpu_time<<setw(16)<<magma_time<<setw(16)<<flag<<endl; 
   }

   void LUconstruct_timing_loop()
   {
       cout<<"Timing LUconstruct complex double:"<<endl;
       cout<<setw(16)<<"N"<<setw(16)<<"cpu_time"<<setw(16)<<"magma_time"<<setw(16)<<"flag"<<endl;
       for (int i = 16; i <= 1087; i += 128)     LUconstruct_timing(i);
       for (int i = 1088; i <= 7232; i += 1024)  LUconstruct_timing(i);
       cout<<"\n\n\n"<<endl;
   }


   void inverse_timing(int N)
   {
       real_Double_t cpu_time, magma_time;
       Matrix<complex<double>,2> X(N,N); fill_random(X);
       Matrix<complex<double>,2> A_cpu(N,N), A_magma(N,N);
       LUDecomp<complex<double>> LU_cpu=LUconstruct_cpu(X);
       LUDecomp<complex<double>> LU_magma=LUconstruct_magma(X);
       
       cpu_time = magma_wtime();
       A_cpu=inverse_cpu( LU_cpu );
       cpu_time = magma_wtime() - cpu_time;

       magma_time = magma_wtime();
       A_magma=inverse_magma( LU_magma );
       magma_time = magma_wtime() - magma_time;

       size_t flag=diff(A_cpu, A_magma,1e-10);
       cout<<setw(16)<<N<<setw(16)<<cpu_time<<setw(16)<<magma_time<<setw(16)<<flag<<endl;
   }


   void inverse_timing_loop()
   {
       cout<<"Timing inverse complex double:"<<endl;
       cout<<setw(16)<<"N"<<setw(16)<<"cpu_time"<<setw(16)<<"magma_time"<<setw(16)<<"flag"<<endl;
       for (int i = 16; i <= 1087; i += 128)     inverse_timing(i);
       for (int i = 1088; i <= 7232; i += 1024)  inverse_timing(i);
       cout<<"\n\n\n"<<endl;
   }


   void bl_cpu_magma_timing()
   {
      //gmm_float_timing_loop();
      //gmm_double_timing_loop();
      //gmm_complexfloat_timing_loop();
      //gmm_complexdouble_timing_loop();
      //eigen_double_timing_loop();
      //eigen_complexdouble_timing_loop();
      //LUconstruct_timing_loop();
      inverse_timing_loop();
   }

}  //end namespace matrix_hao_lib

#endif
