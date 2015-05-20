#ifndef USE_BL_CPU
#define USE_BL_CPU

//WARNING: it is not safe to use this header
//Another way to swith between cpu and gpu is: 
//#define BL_NAME(x) x##_cpu
//call gmm_cpu by: BL_NAME(gmm)

#define gmm          gmm_cpu
#define eigen        eigen_cpu
#define LUconstruct  LUconstruct_cpu 
#define inverse      inverse_cpu
#define solve_lineq  solve_lineq_cpu
#define QRMatrix     QRMatrix_cpu

#endif
