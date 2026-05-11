# Vladimir S. FONOV - vladimir.fonov@gmail.com
#
PROJECT (PATCH_MORPHOLOGY)
ENABLE_TESTING()

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
if(COMMAND cmake_policy)
 cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)

IF(NOT MINC_TOOLKIT_BUILD)
  # ITK is optional: probe for it but do not require.
  # When ITK isn't available, the ITK-based src/ tree is skipped.
  FIND_PACKAGE(ITK)
  OPTION(USE_BLAS    "Use BLAS" OFF)

  IF(USE_BLAS)
    set(BLA_VENDOR OpenBLAS) # search for OpenBlas for now
    find_package(BLAS)
    if(BLAS_FOUND AND NOT TARGET BLAS::BLAS)
      add_library(BLAS::BLAS INTERFACE IMPORTED)
      set_target_properties(BLAS::BLAS PROPERTIES
        INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
    endif()
  ENDIF()
  OPTION(PATCH_MORPHOLOGY_BUILD_LEGACY    "Build legacy program volume_patchs, requires LIBMINC" OFF)
ELSE()
  SET(PATCH_MORPHOLOGY_BUILD_LEGACY ON)
  # ITK is optional in the superbuild path: probe for it but do not require.
  # When ITK isn't available, the ITK-based src/ tree is skipped.
  FIND_PACKAGE(ITK QUIET)
  # Superbuild path: parent forwards BLAS+LAPACKE state via -D flags.
  # Reconstruct both BLAS::BLAS and LAPACKE::LAPACKE so src/ can link uniformly.
  include(${SUPERBUILD_CMAKE_DIR}/BLASTargetShim.cmake)
  include(${SUPERBUILD_CMAKE_DIR}/LAPACKETargetShim.cmake)
  if(BLAS_LIBRARIES OR BLAS_MKL_MODE STREQUAL "Config")
    if(BLAS_MKL_MODE STREQUAL "Config")
      find_package(MKL CONFIG REQUIRED)
      blas_create_target_shim("MKL::MKL" "" "")
    else()
      blas_create_target_shim("${BLAS_LIBRARIES}" "${BLAS_LINKER_FLAGS}" "${BLAS_INCLUDE_DIRS}")
    endif()
    lapacke_create_target_shim("${LAPACKE_LIBRARIES}" "${LAPACKE_INCLUDE_DIRS}" "${HAVE_LAPACKE}")
    set(USE_BLAS ON)
  endif()
ENDIF()


# nnls.c is the only BLAS-using source path here, and it needs lapacke.h.
# Gate HAVE_OPENBLAS on both BLAS and LAPACKE being available — without
# LAPACKE the NNLS code path can't compile.
IF(TARGET BLAS::BLAS AND USE_BLAS AND HAVE_LAPACKE)
  SET(HAVE_OPENBLAS ON)
  ADD_DEFINITIONS(-DHAVE_OPENBLAS -DUSE_SPAMS)
  message("Building with BLAS+LAPACKE")
ELSE()
 SET(HAVE_OPENBLAS OFF)
 message("Building without BLAS (need both BLAS and LAPACKE for the NNLS path)")
ENDIF()

IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
  INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS})
  LINK_DIRECTORIES(${ITK_LIBRARY_DIRS})
  add_subdirectory(src)
ENDIF(ITK_FOUND)

add_subdirectory(scripts)

IF(PATCH_MORPHOLOGY_BUILD_LEGACY)
  add_subdirectory(legacy)
ENDIF(PATCH_MORPHOLOGY_BUILD_LEGACY)

IF(BUILD_TESTING AND ITK_FOUND)
  ADD_SUBDIRECTORY(testdir)
ENDIF()
