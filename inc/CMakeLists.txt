set(header matrix_define.h matrix_base.h matrix_1d.h matrix_2d.h matrix_3d.h matrix_element_wise.h
    matrix_mpi.h matrix_2d_common_fun.h matrix_2d_bl_cpu.h matrix_2d_bl_magma.h pfaffian_utilities.h
    use_bl_cpu.h use_bl_magma.h matrix_all.h acml.h)

install(FILES ${header}  DESTINATION include)
