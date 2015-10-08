IF(NOT MINC_TOOLKIT_BUILD)
  FIND_PACKAGE(LIBMINC REQUIRED)
ENDIF(NOT MINC_TOOLKIT_BUILD)



ADD_EXECUTABLE( volume_patches
                volume_patches.cpp 
                utils.cpp  
                utils.h
              )
TARGET_LINK_LIBRARIES(volume_patches ${EZMINC_LIBRARIES} )


INSTALL(TARGETS 
          volume_patches
        RUNTIME DESTINATION bin )
