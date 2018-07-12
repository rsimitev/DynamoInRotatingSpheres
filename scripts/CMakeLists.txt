include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

install(PROGRAMS compressdata 
                 uncompressdata 
                 rebuildTimeSeries.sh 
                 resubmitJob.sh
                 drs_cleanRun.sh
                 stopRun.sh
        DESTINATION bin)

# Only install these if the required binary utils were compiled
if(BUILD_UTILS)
   install(PROGRAMS Yokoi-rs-drive.pl DESTINATION bin)
   configure_file(drs_plot.sh.in drs_plot.sh @ONLY)
   install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/drs_plot.sh DESTINATION bin)
endif()

# Only install these if gnuplot is installed on the system
if(GNUPLOT_FOUND)
   install(PROGRAMS latestNormalizedSpectra.sh DESTINATION bin)
   install(PROGRAMS drs_KE_plot.sh DESTINATION bin)
   install(FILES monitor_functions.sh DESTINATION lib/DRS)
endif(GNUPLOT_FOUND)

