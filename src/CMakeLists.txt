# Define the source files
file(GLOB KBLAS_SRC "*.cpp")

include_directories( ${PROJECT_SOURCE_DIR}/include)

# Defines the KBLAS library to build
add_library(kblas STATIC ${KBLAS_SRC})

install( TARGETS kblas
         DESTINATION lib
         )

