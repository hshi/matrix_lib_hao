#include "matrix_all.h"

namespace matrix_hao_lib
{
 void matrix_class_test();
 void matrix_2d_common_fun_test();
 void matrix_2d_bl_cpu_test();
#ifdef USE_MAGMA
 void matrix_2d_bl_magma_test();
#endif
#ifdef MPI_HAO
 void matrix_mpi_test();
#endif
}

using namespace std;
using namespace matrix_hao_lib;

int main(int argc, char** argv)
{
    int rank=0;

#ifdef MPI_HAO
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
#ifdef USE_MAGMA
    magma_init();
#endif

    if(rank==0)  cout<<"\n\n\n=======Testing======="<<endl;
    matrix_class_test();
    matrix_2d_common_fun_test(); 
    matrix_2d_bl_cpu_test();
#ifdef USE_MAGMA
    matrix_2d_bl_magma_test();
#endif
#ifdef MPI_HAO
    matrix_mpi_test();
#endif


#ifdef USE_MAGMA
    magma_finalize();
#endif
#ifdef MPI_HAO
    MPI_Finalize();
#endif

    return 0;
}
