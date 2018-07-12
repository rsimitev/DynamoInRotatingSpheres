
include_directories(${DRS_BINARY_DIR})
link_directories(${DRS_BINARY_DIR})

set(drs-state-average_SRCS
   drs-state-average.F90
   )
add_executable(drs-state-average ${drs-state-average_SRCS})
target_link_libraries(drs-state-average ${DRS_LIB} ${DRS_LINK_LIBRARIES})
install(TARGETS drs-state-average RUNTIME DESTINATION bin)

set(drs-state-change_SRCS
   drs-state-change.F90
   )
add_executable(drs-state-change ${drs-state-change_SRCS})
target_link_libraries(drs-state-change ${DRS_LIB} ${DRS_LINK_LIBRARIES})
install(TARGETS drs-state-change RUNTIME DESTINATION bin)

set(Benchmark_SRCS
   Benchmark.F90
   )
add_executable(Benchmark ${Benchmark_SRCS})
target_link_libraries(Benchmark ${DRS_LIB} ${DRS_LINK_LIBRARIES})
install(TARGETS Benchmark RUNTIME DESTINATION bin)

set(Benchmark2_SRCS
   Benchmark-v2.F90
   )
add_executable(Benchmark2 ${Benchmark2_SRCS})
target_link_libraries(Benchmark2 ${DRS_LIB} ${DRS_LINK_LIBRARIES})
install(TARGETS Benchmark2 RUNTIME DESTINATION bin)

set(Benchmark-gen_SRCS
   Benchmark-gen.F90
   )
add_executable(Benchmark-gen ${Benchmark-gen_SRCS})
target_link_libraries(Benchmark-gen ${DRS_LIB} ${DRS_LINK_LIBRARIES})
install(TARGETS Benchmark-gen RUNTIME DESTINATION bin)

set(Yokoi_plots_SRCS
   Yokoi-plots.F90
   )
add_executable(Yokoi-plots ${Yokoi_plots_SRCS})
target_link_libraries(Yokoi-plots ${DRS_LIB} ${DRS_LINK_LIBRARIES})
install(TARGETS Yokoi-plots RUNTIME DESTINATION bin)

set(Yokoi_plots_rs_SRCS
   Yokoi-plots-rs.F90
   )
add_executable(Yokoi-plots-rs ${Yokoi_plots_rs_SRCS})
target_link_libraries(Yokoi-plots-rs ${DRS_LIB} ${DRS_LINK_LIBRARIES})
install(TARGETS Yokoi-plots-rs RUNTIME DESTINATION bin)

set(drs2dx_SRCS
   drs2dx.F90
   drs_renderers.F90
   drs_io_DX.F90
   )
add_executable(drs2dx ${drs2dx_SRCS})
target_link_libraries(drs2dx ${DRS_LIB} ${DRS_LINK_LIBRARIES})
install(TARGETS drs2dx RUNTIME DESTINATION bin)

set(drs_spec_SRCS
   drs-spectra.F90
   )
add_executable(drs-spectra ${drs_spec_SRCS})
target_link_libraries(drs-spectra ${DRS_LIB} ${DRS_LINK_LIBRARIES})
install(TARGETS drs-spectra RUNTIME DESTINATION bin)

#
set(drs-get-profile_SRCS drs-get-profile.F90)
add_executable(drs-get-profile ${drs-get-profile_SRCS})
target_link_libraries(drs-get-profile ${DRS_LIB} ${DRS_LINK_LIBRARIES})
install(TARGETS drs-get-profile RUNTIME DESTINATION bin)

set(drs-version_SRCS drs-version.F90)
add_executable(drs-version ${drs-version_SRCS})
install(TARGETS drs-version RUNTIME DESTINATION bin)
