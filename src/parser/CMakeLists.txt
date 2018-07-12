include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SURRENT_SOURCE_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(parser parser.F90)
install(TARGETS parser ARCHIVE DESTINATION lib)

if(BUILD_TESTS)
   add_subdirectory(tests)
endif()
