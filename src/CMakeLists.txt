INCLUDE_DIRECTORIES(
  ${ITK_INCLUDE_DIRS}
  itk
  3rd_party/libsvm/
  3rd_party/nnls/
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

ADD_EXECUTABLE( itk_patch_grading
                itk_patch_grading.cpp
                3rd_party/libsvm/svm.cpp
                3rd_party/nnls/nnls.c
                )
                
TARGET_LINK_LIBRARIES( itk_patch_grading  
                ${OpenBLAS_LIBRARIES} 
                ${ITK_LIBRARIES} )
                

ADD_EXECUTABLE(itk_patch_segmentation
                itk_patch_segmentation.cpp
                3rd_party/libsvm/svm.cpp
                3rd_party/nnls/nnls.c
                )

TARGET_LINK_LIBRARIES( itk_patch_segmentation   
                ${OpenBLAS_LIBRARIES} 
                ${ITK_LIBRARIES} )
                
                     
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
         itk_patch_grading 
         itk_patch_segmentation
        RUNTIME DESTINATION bin )
