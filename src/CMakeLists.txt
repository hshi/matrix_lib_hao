set(src_lib matrixhao)

file(GLOB srcs "*.cpp")

add_library(${src_lib} STATIC ${srcs})

install(TARGETS ${src_lib}  DESTINATION lib)
