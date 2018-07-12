
include_directories(${CURRENT_CMAKE_BINARY_DIR})
link_directories(${CURRENT_CMAKE_BINARY_DIR})

# Test the azimuthal Fourier transforms
add_executable(test_drs_fftw test_drs_fftw.F90)
target_link_libraries(test_drs_fftw ${DRS_LIB} ${DRS_LINK_LIBRARIES})
add_test(NAME drs_fftw COMMAND test_drs_fftw)
set_tests_properties(drs_fftw PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")

# Test the drs_radial module
add_executable(test_drs_radial test_drs_radial.F90)
target_link_libraries(test_drs_radial ${DRS_LIB} ${DRS_LINK_LIBRARIES})
# Test the initialization of the collocated points
add_test(NAME collocated_init COMMAND test_drs_radial col)
set_tests_properties( collocated_init PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")
# Test the first radial derivatives
add_test(NAME radial_derivative_r2r COMMAND test_drs_radial dr1)
set_tests_properties( radial_derivative_r2r PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")
# Test the first and second radial derivatives
add_test(NAME radial_dr_ddr_1D_r2r COMMAND test_drs_radial 1D_r2r)
set_tests_properties( radial_dr_ddr_1D_r2r PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")

# Test the Chebyshev module
add_executable(test_drs_Chebyshev test_drs_Chebyshev.F90)
target_link_libraries(test_drs_Chebyshev  ${DRS_LIB} ${DRS_LINK_LIBRARIES})
# Configure the test wrapper
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/test_drs_Chebyshev_wrapper.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/test_drs_Chebyshev_wrapper.sh"
  )
# Test the correct calculation of the polynomials
add_test(NAME chebichev_poly_val COMMAND test_drs_Chebyshev_wrapper.sh che)
set_tests_properties( chebichev_poly_val PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")
# Test the correct calculation of the first derivative of the polynomials
add_test(NAME chebichev_poly_dx COMMAND test_drs_Chebyshev_wrapper.sh  cdx)
set_tests_properties( chebichev_poly_dx PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")
# Test a n -> x ->n roundtrip
add_test(NAME chebichev_roundtrip COMMAND test_drs_Chebyshev_wrapper.sh rdt)
set_tests_properties( chebichev_roundtrip PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")


# Test the cosine transforms
add_executable(test_drs_fftw-r2r test_drs_fftw-r2r.F90)
target_link_libraries(test_drs_fftw-r2r ${DRS_LIB} ${DRS_LINK_LIBRARIES})
add_test(NAME drs_fftw-r2r COMMAND test_drs_fftw-r2r)
set_tests_properties(drs_fftw-r2r PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")

# Test the radial remeshing
add_executable(test_remesh-r2r test_remesh-r2r.F90)
target_link_libraries(test_remesh-r2r ${DRS_LIB} ${DRS_LINK_LIBRARIES})
add_test(NAME remesh-r2r COMMAND test_remesh-r2r)
set_tests_properties(remesh-r2r PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")

# Test parsing of a full and correct config file
#add_executable(test_drs_read_conf_v3 test_drs_read_conf_v3.F90)
#target_link_libraries(test_drs_read_conf_v3 ${DRS_LIB} ${DRS_LINK_LIBRARIES})
#add_test(NAME drs_read_conf_v3 COMMAND
#   ${CMAKE_CURRENT_SOURCE_DIR}/wrap_test_drs_read_conf_v3.sh ${CMAKE_CURRENT_SOURCE_DIR}/data/drs.conf)
#set_tests_properties(drs_read_conf_v3 PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")

set(test_saveDXmeridional_SRCS test_saveDXmeridional.F90)
add_executable(test_saveDXmeridional ${test_saveDXmeridional_SRCS})
target_link_libraries(test_saveDXmeridional ${DRS_LIB} ${DRS_LINK_LIBRARIES})

set(test_logFeature_SRCS test_logFeature.F90)
add_executable(test_logFeature ${test_logFeature_SRCS})
target_link_libraries(test_logFeature ${DRS_LIB} ${DRS_LINK_LIBRARIES})

# test the divergence of a vector field in sperical coordinates and real
# space.
add_executable(test_vectorField2Divergence test_vectorField2Divergence.F90)
target_link_libraries(test_vectorField2Divergence ${DRS_LIB} ${DRS_LINK_LIBRARIES})
add_test(NAME vectorField2Divergence COMMAND test_vectorField2Divergence)
set_tests_properties( vectorField2Divergence PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")
