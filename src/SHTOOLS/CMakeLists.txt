
include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})
set(SHTOOLS_LIB shtools)

set(SHTOOLS_LIB_SRCS 
    PLegendreA_d1.f90
    PlmBar_d1.f90
    PlmIndex.f90
    )

add_library(${SHTOOLS_LIB} ${SHTOOLS_LIB_SRCS})
install(TARGETS ${SHTOOLS_LIB} ARCHIVE DESTINATION lib)
