INSTALL(PROGRAMS 
        hcag_segmentation_pipeline.pl
        patch_segmentation_pipeline.pl
        DESTINATION bin )

IF(PATCH_MORPHOLOGY_BUILD_LEGACY)
  INSTALL(PROGRAMS 
          ventricles_segmentation_pipeline.pl
          patch_segmentation_pipeline.pl
          DESTINATION bin )
ENDIF(PATCH_MORPHOLOGY_BUILD_LEGACY)