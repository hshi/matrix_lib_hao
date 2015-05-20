#ifndef USE_BL_MAGMA_H
#define USE_BL_MAGMA_H

//WARNING: it is not safe to use this header
//Another way to swith between cpu and gpu is: 
//#define BL_NAME(x) x##_magma
//call gmm_magma by: BL_NAME(gmm)

#define gmm          gmm_magma
#define eigen        eigen_magma
#define LUconstruct  LUconstruct_magma
#define inverse      inverse_magma
#define solve_lineq  solve_lineq_magma
#define QRMatrix     QRMatrix_magma

#endif
