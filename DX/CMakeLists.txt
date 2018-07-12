include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

install(FILES ClipCyllinder.net
              ColouredIsolines.net
              MapToSphere.net
              NormMinMax.net
              rtp2xyz.net
              DESTINATION ${DXMACROS})

configure_file(EcutBatch-scal.net.in EcutBatch-scal.net @ONLY)
configure_file(EcutBatch-vect-r.net.in EcutBatch-vect-r.net @ONLY)
configure_file(EcutBatch-vect-p.net.in EcutBatch-vect-p.net @ONLY)
configure_file(McutBatch-scal.net.in McutBatch-scal.net @ONLY)
configure_file(McutBatch-vect-r.net.in McutBatch-vect-r.net @ONLY)
configure_file(McutBatch-vect-p.net.in McutBatch-vect-p.net @ONLY)
configure_file(VcutBatch-vect-outside.net.in VcutBatch-vect-outside.net @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/EcutBatch-scal.net
              ${CMAKE_CURRENT_BINARY_DIR}/EcutBatch-vect-r.net 
              ${CMAKE_CURRENT_BINARY_DIR}/EcutBatch-vect-p.net 
              ${CMAKE_CURRENT_BINARY_DIR}/McutBatch-scal.net 
              ${CMAKE_CURRENT_BINARY_DIR}/McutBatch-vect-r.net 
              ${CMAKE_CURRENT_BINARY_DIR}/McutBatch-vect-p.net 
              ${CMAKE_CURRENT_BINARY_DIR}/VcutBatch-vect-outside.net 
              # EquatorialCut.net 
              # fieldLines.net 
              # flowAndTemperature.net 
              # sequencePlotterEquat.net 
              # sequencePlotterMerid.net 
              # sequencePlotterVolume.net 
              # Yokoi2.net 
              # Yokoi.net         
              DESTINATION ${DXPROGS})
