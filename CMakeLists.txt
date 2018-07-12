# CMake project file for DRS

cmake_minimum_required (VERSION 2.8)
project(DRS)
enable_language (Fortran)
set (CMAKE_Fortran_FORMAT FREE)
set (DRS_VERSION_MAJOR 1)
set (DRS_VERSION_MINOR 70)
set (DRS_VERSION_PATCH 0)
set (DRS_VERSION
   "${DRS_VERSION_MAJOR}.${DRS_VERSION_MINOR}.${DRS_VERSION_PATCH}" )

SET(CMAKE_MODULE_PATH ${DRS_SOURCE_DIR}/cmake)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# make sure that the default is a RELEASE
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
   set (CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build, options are: None Debug Release." FORCE)
endif()
# default installation
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   set (CMAKE_INSTALL_PREFIX ${DRS_SOURCE_DIR} CACHE STRING "Choose the installation directory." FORCE)
endif()

# FFLAGS depend on the compiler
# This is a workaround for the intel compiler that puts
# all the arrays in main memory thus exceeding the maximum
# allocation limit very quickly. It may also be specified by hand.
message("Using ${CMAKE_Fortran_COMPILER_ID} compiler")
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
   add_compile_options(-heap-arrays)
endif ()

option(MPI "Compiles the parallel version of the code" OFF)
option(COMP "Compiles the doubly diffusive version of the code" OFF)
option(BUILD_UTILS "Build the utilities associated with the code." OFF)
option(BUILD_TESTS "Build a series of driver programs that test several subroutines of the code." OFF)

string ( TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)
if (BUILD_TYPE STREQUAL "debug")
   add_definitions(-DDEBUG)
endif()

find_package(FFTW 3.3 REQUIRED)
# add a target to generate API documentation with Doxygen.
find_package(Doxygen 1.8)
# Needed to build the user manual.
find_package(LATEX)
# Needed for the plotting scripts.
find_package(Gnuplot 4.4)
# Needed for the IDL plotting scripts.
find_package(IDL)

set (DXMACROS $ENV{HOME}/.dx.d/macros CACHE PATH "Choose the DX macro installation path." )
set (DXPROGS  $ENV{HOME}/.dx.d/programs CACHE PATH "Choose the DX networks installation path." )

if(PDFLATEX_COMPILER)
  option(BUILD_DOCS "Build the user manual from the latex sources." OFF)
endif(PDFLATEX_COMPILER)

if(MPI)
  find_package(MPI 1.6 REQUIRED)
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  add_definitions(-DMPI)
  message("Building with MPI support.")
  LIST(APPEND FEATURES MPI)
endif(MPI)

if(COMP)
  add_definitions(-DCOMP)
  LIST(APPEND FEATURES DD)
  message("Building DRS with support for compositional convection.")
endif(COMP)

if(BUILD_DOCS)
   add_subdirectory(docs)
endif(BUILD_DOCS)
if(BUILD_TESTS)
   enable_testing()
endif(BUILD_TESTS)

# Now that we know all the features we can build a version string
set (DRS_VERSION_STRING "${DRS_VERSION}(${FEATURES}) - ${CMAKE_BUILD_TYPE}")

include_directories(${DRS_BINARY_DIR})

add_subdirectory(src)
add_subdirectory(scripts)
add_subdirectory(DX)

if(DOXYGEN_FOUND)
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
   add_custom_target(APIDOX
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
      )
endif(DOXYGEN_FOUND)

INCLUDE(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/BUILD/;/bin/;test*;/.git*;/.kate*;run-test*")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${DRS_SOURCE_DIR}/README")
set(CPACK_PACKAGE_VENDOR "The University of Glasgow")
set(CPACK_PACKAGE_VERSION_MAJOR ${DRS_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${DRS_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${DRS_VERSION_PATCH})

INCLUDE(CPack)

