#ifndef RANDOM_FILL_MATRIX_H
#define RANDOM_FILL_MATRIX_H

namespace matrix_hao_lib
{
    BL_INT lapack_ran_ISEED[4] = { 0, 127, 0, 127 };
    template <int D>
    void fill_random(Matrix<float,D> &A)
    {
        BL_INT itwo = 2; BL_INT size_A = A.L_f();
        FORTRAN_NAME(slarnv)(&itwo, lapack_ran_ISEED, &size_A, (BL_FLOAT*) A.base_array);
    }

    template <int D>
    void fill_random(Matrix<double,D> &A)
    {
        BL_INT itwo = 2; BL_INT size_A = A.L_f();
        FORTRAN_NAME(dlarnv)(&itwo, lapack_ran_ISEED, &size_A, (BL_DOUBLE*) A.base_array);
    }

    template <int D>
    void fill_random(Matrix<std::complex<float>,D> &A)
    {
        BL_INT itwo = 2; BL_INT size_A = A.L_f();
        FORTRAN_NAME(clarnv)(&itwo, lapack_ran_ISEED, &size_A, (BL_COMPLEX8*) A.base_array);
    }

    template <int D>
    void fill_random(Matrix<std::complex<double>,D> &A)
    {
        BL_INT itwo = 2; BL_INT size_A = A.L_f();
        FORTRAN_NAME(zlarnv)(&itwo, lapack_ran_ISEED, &size_A, (BL_COMPLEX16*) A.base_array);
    }

}  //end namespace matrix_hao_lib

#endif
