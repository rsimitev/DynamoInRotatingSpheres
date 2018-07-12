include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SURRENT_SOURCE_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(nr lubksb.f90 ludcmp.f90 matinv.f90)
install(TARGETS nr ARCHIVE DESTINATION lib)
