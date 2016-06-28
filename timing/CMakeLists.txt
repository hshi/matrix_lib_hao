set(time_exe matrixtiming)

file(GLOB time_srcs "*.cpp")

add_executable(${time_exe} ${time_srcs})

target_link_libraries(${time_exe} matrixhao)

target_link_libraries(${time_exe} ${lib_all})

install(TARGETS ${time_exe}  DESTINATION timing)
