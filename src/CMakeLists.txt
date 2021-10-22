INCLUDE_DIRECTORIES(
  ${ITK_INCLUDE_DIRS}
  itk
  3rd_party/libsvm/
  )

IF(HAVE_OPENBLAS) # TODO: make this generic for any BLAS
  INCLUDE_DIRECTORIES(
    ${OpenBLAS_INCLUDE_DIRS}
    3rd_party/nnls/
  )
ENDIF(HAVE_OPENBLAS)

LINK_DIRECTORIES( ${ITK_LIBRARY_DIRS} )
  
ADD_EXECUTABLE( itk_patch_morphology
                itk_patch_morphology.cpp )

TARGET_LINK_LIBRARIES( itk_patch_morphology   
                      ${ITK_LIBRARIES} 
                      )

ADD_EXECUTABLE( itk_patch_morphology_mc
                itk_patch_morphology_mc.cpp )

TARGET_LINK_LIBRARIES( itk_patch_morphology_mc 
                      ${ITK_LIBRARIES} 
                      )

ADD_EXECUTABLE( itk_minc_nonlocal_filter 
                itk_minc_nonlocal_filter.cpp )

TARGET_LINK_LIBRARIES( itk_minc_nonlocal_filter  
                      ${ITK_LIBRARIES} 
                     )

SET(  itk_patch_grading_SRC
                itk_patch_grading.cpp
                3rd_party/libsvm/svm.cpp
                )

IF(HAVE_OPENBLAS)
  SET(itk_patch_grading_SRC ${itk_patch_grading_SRC} 3rd_party/nnls/nnls.c) # NOT implemented anyway
ENDIF(HAVE_OPENBLAS)


ADD_EXECUTABLE( itk_patch_grading ${itk_patch_grading_SRC})

TARGET_LINK_LIBRARIES( itk_patch_grading
                ${ITK_LIBRARIES} )
                
IF(HAVE_OPENBLAS)
  TARGET_LINK_LIBRARIES( itk_patch_grading
                  ${OpenBLAS_LIBRARIES} )
ENDIF(HAVE_OPENBLAS)

SET( itk_patch_segmentation_SRC
                itk_patch_segmentation.cpp
                3rd_party/libsvm/svm.cpp
                )
                
IF(HAVE_OPENBLAS)
  SET(itk_patch_segmentation_SRC ${itk_patch_grading_SRC} 3rd_party/nnls/nnls.c)
ENDIF(HAVE_OPENBLAS)

ADD_EXECUTABLE( itk_patch_segmentation ${itk_patch_segmentation_SRC})

TARGET_LINK_LIBRARIES( itk_patch_segmentation
                ${ITK_LIBRARIES} )
                
IF(HAVE_OPENBLAS)
  TARGET_LINK_LIBRARIES( itk_patch_segmentation
                  ${OpenBLAS_LIBRARIES} )
ENDIF(HAVE_OPENBLAS)
                

ADD_EXECUTABLE( itk_split_labels 
                itk_split_labels.cpp
                )

TARGET_LINK_LIBRARIES( itk_split_labels
                ${ITK_LIBRARIES} )
                
ADD_EXECUTABLE( itk_merge_labels
                itk_merge_labels.cpp
                )

TARGET_LINK_LIBRARIES( itk_merge_labels
                ${ITK_LIBRARIES} )
 
ADD_EXECUTABLE( itk_merge_discrete_labels
                itk_merge_discrete_labels.cpp
                )

TARGET_LINK_LIBRARIES( itk_merge_discrete_labels
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
    itk_split_labels
    itk_merge_labels
    itk_merge_discrete_labels
        RUNTIME DESTINATION bin )
