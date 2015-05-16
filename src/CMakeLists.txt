set(matrix_srcs matrix_2d_common_fun.cpp matrix_2d_bl_cpu.cpp matrix_mpi.cpp)

add_library(matrixhao STATIC ${matrix_srcs})

install(TARGETS matrixhao  DESTINATION lib)
