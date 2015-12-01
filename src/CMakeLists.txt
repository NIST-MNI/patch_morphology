INCLUDE_DIRECTORIES(
  ${ITK_INCLUDE_DIRS}
  itk
  )

LINK_DIRECTORIES( ${ITK_LIBRARY_DIRS} )
  
ADD_EXECUTABLE(itk_patch_morphology
                itk_patch_morphology.cpp )

TARGET_LINK_LIBRARIES( itk_patch_morphology   
                      ${ITK_LIBRARIES} 
                      )

ADD_EXECUTABLE(itk_patch_morphology_mc
                itk_patch_morphology_mc.cpp )

TARGET_LINK_LIBRARIES( itk_patch_morphology_mc 
                      ${ITK_LIBRARIES} 
                      )

ADD_EXECUTABLE(itk_minc_nonlocal_filter 
                itk_minc_nonlocal_filter.cpp )

TARGET_LINK_LIBRARIES( itk_minc_nonlocal_filter  
                      ${ITK_LIBRARIES} 
                     )

# Experimantal stuff
# ADD_EXECUTABLE(itk_levelset_segment 
#                 itk_levelset_segment.cpp )
# 
# TARGET_LINK_LIBRARIES( itk_levelset_segment  ${ITK_LIBRARIES} )
# 
 
#ADD_EXECUTABLE(itk_patch_normalization
#    itk_patch_normalization.cpp )

#TARGET_LINK_LIBRARIES( itk_patch_normalization  
#                      ${ITK_LIBRARIES} 
#          )



INSTALL(TARGETS 
         itk_patch_morphology 
         itk_patch_morphology_mc 
         itk_minc_nonlocal_filter 
        RUNTIME DESTINATION bin )
