include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add the shtools lib to the build
add_subdirectory(SHTOOLS)

add_subdirectory(NR)

# Add the parser library to the build
add_subdirectory(parser)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/parser)

configure_file(drsDefs.F90.in drsDefs.F90 @ONLY)

set(DRS_LIB drs)

set(DRS_LIB_SRCS
    drs_field.F90
    drs_hypDiff.F90
    drs_legendre.F90
    drs_transforms.F90
    drs_real_space.F90
    drs_probes.F90
    drs_temp.F90
    drs_dims.F90
    drs_fftw3.F90
    drs_flow.F90
    drs_io_units.F90
    drs_io_state.F90
    drs_io.F90
    drs_mpi.F90
    drs_params.F90
    drs_Chebyshev.F90
    drs_radial.F90
    drs_time.F90
    drs_debug.F90
    drs_lock.F90
    drs_error_codes.F90
    )

if(COMP)
   set(DRS_LIB_SRCS ${DRS_LIB_SRCS} drs_comp.F90)
endif(COMP)

add_library(${DRS_LIB} ${DRS_LIB_SRCS})
target_link_libraries(${DRS_LIB} ${FFTW_LIBRARIES} shtools parser)
install(TARGETS ${DRS_LIB} ARCHIVE DESTINATION lib)

set(DRS_LINK_LIBRARIES ${DRS_LIB} ${FFTW_LIBRARIES} shtools parser)
if(MPI)
   set(DRS_LINK_LIBRARIES ${DRS_LINK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
endif(MPI)

# build the main program.
set(DRS drs.exe)
set(DRS_SRCS
   drs_heat_equation.F90
   drs_momentum_equation.F90
   drs_induction_equation.F90
   drs_io_conf.F90
   drs.F90)
if(COMP)
   set(DRS_SRCS drs_composition_equation.F90 ${DRS_SRCS})
endif(COMP)

add_executable (${DRS} ${DRS_SRCS})
target_link_libraries(${DRS} ${DRS_LINK_LIBRARIES} nr)
install(TARGETS ${DRS} RUNTIME DESTINATION bin)

add_executable (pureHydroDyn drs_momentum_equation.F90 pureHydroDyn.F90)
target_link_libraries(pureHydroDyn ${DRS_LINK_LIBRARIES} nr)
install(TARGETS pureHydroDyn RUNTIME DESTINATION bin)

if(BUILD_UTILS)
   add_subdirectory(utilities)
else()
   message("- Will not build the utilities.")
endif()

if(BUILD_TESTS)
   add_subdirectory(tests)
endif()
