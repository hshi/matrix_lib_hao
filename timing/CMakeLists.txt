set(matrix_timing_srcs bl_cpu_magma_timing.cpp timing_all.cpp)

add_executable(matrixtiming ${matrix_timing_srcs})
target_link_libraries(matrixtiming matrixhao)
target_link_libraries(matrixtiming ${magma_libraries})
target_link_libraries(matrixtiming ${lapackblas_libraries})

install(TARGETS matrixtiming DESTINATION timing)
