set(matrix_test_srcs matrix_2d_bl_cpu.cpp matrix_2d_common_fun_test.cpp matrix_class_test.cpp matrix_mpi_test.cpp test_all.cpp)

add_executable(matrixtest ${matrix_test_srcs})

target_link_libraries(matrixtest matrixhao)

target_link_libraries(matrixtest ${lapackblas_libraries})

install(TARGETS matrixtest DESTINATION test)
