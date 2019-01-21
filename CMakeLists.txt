PROJECT(ITK-CLIs)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11) #2.8.11 for TARGET_INCLUDE_DIRECTORIES

## bash command to check if CMakeLists.txt contains an entry for all *.cxx: 
## for i in *cxx ; do if ! `grep -q ${i%.cxx} CMakeLists.txt`; then echo "$i not found"'!'; fi; done


##from http://www.cmake.org/pipermail/cmake/2008-September/023808.html
#i.e. this is bad: SET(CMAKE_BUILD_TYPE "Release")
# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their definition
# and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()

#${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH} 


IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()


IF(DEFINED CMAKE_INSTALL_PREFIX)
   SET(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "Choose the install directory.")
ELSE()
   SET(CMAKE_INSTALL_PREFIX /opt/itk-simple/ CACHE STRING "Choose the install directory.")
ENDIF()

# IF(DEFINED ITK_DIR)
#    SET(ITK_DIR ${ITK_DIR})
# ELSE()
#    SET(ITK_DIR /opt/itk-4.6.1/lib/cmake/ITK-4.6/)
# ENDIF()

##-Werror for clean programming ;-)
# SET(CMAKE_CXX_FLAGS "") #dynamic linking, comment out for static linking
# IF(DEFINED CMAKE_CXX_FLAGS)
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "NOT compiled with static libc, libgcc and libstdc++")
# ELSE()
#    SET(CMAKE_CXX_FLAGS "-static -static-libgcc -static-libstdc++ -Werror" CACHE STRING "compile with static libc, libgcc and libstdc++")
# ENDIF()




#set search path
SET(CMAKE_PREFIX_PATH "/opt/itk-3.20/lib/InsightToolkit/" ${CMAKE_PREFIX_PATH})
SET(CMAKE_PREFIX_PATH "/opt/itk-4.1.0/lib/cmake/ITK-4.1/" ${CMAKE_PREFIX_PATH})


SET(Required_Packages
  ITK
)


FIND_PACKAGE (ITK REQUIRED)
IF ( ITK_FOUND )
      INCLUDE( ${ITK_USE_FILE} ) # https://stackoverflow.com/questions/33583101/cmake-3-4-0-including-multiple-files#33585398
ENDIF( ITK_FOUND )


INCLUDE_DIRECTORIES(
  # $ENV{HOME}/itk/parabolicMorphology/
  # $ENV{HOME}/itk/ITKbinaryThinningImageFilter3D/Source/
  # $ENV{HOME}/itk/itkWindowImageFilter/
  # $ENV{HOME}/itk/itkLabelSetErodeDilate_IJ_886/
  # filter/internal/
  ${CMAKE_CURRENT_SOURCE_DIR}/tclap/include/
  )

SET(Libraries
  ${ITK_LIBRARIES}
)


SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags") # use RPATH instead of RUNPATH (currently ignored on Debian-9), -WL tells gcc to pass option to the linker

## check if specific modules are enabled
LIST(FIND ITK_MODULES_ENABLED ITKReview ITKReview_FOUND)
LIST(FIND ITK_MODULES_ENABLED LesionSizingToolkit LesionSizingToolkit_FOUND)
LIST(FIND ITK_MODULES_ENABLED MinimalPathExtraction MinimalPathExtraction_FOUND)#in itk-4.8: https://github.com/Kitware/ITK/commit/da17ae4d4955bb92d5259c6ac47bfce6c474be18
LIST(FIND ITK_MODULES_ENABLED ParabolicMorphology ParabolicMorphology_FOUND)#in itk-4.8: https://github.com/InsightSoftwareConsortium/ITK/commit/2b3190575b82e8457c028040c3bfd0f69fa1ef35
LIST(FIND ITK_MODULES_ENABLED LabelErodeDilate LabelErodeDilate_FOUND)#in itk-4.8: https://github.com/Kitware/ITK/commit/5b158490c01b3e3a68816fdbe321751609746444
LIST(FIND ITK_MODULES_ENABLED StreamingSinc StreamingSinc_FOUND)


SET(CurrentExe "file_converter")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "file_converter_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "file-series_reader")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "file-series_reader_SDI")
ADD_EXECUTABLE(${CurrentExe}  "file-series_reader.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_SDI")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "file-series_writer")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "file-series_writer_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "subimage_extract")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "subimage_extract_SDI")
ADD_EXECUTABLE(${CurrentExe}  "subimage_extract.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_SDI")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "subimage_extract_series")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "subimage_extract_series_SDI")
ADD_EXECUTABLE(${CurrentExe}  "subimage_extract_series.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_SDI")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "subimage_mask")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "subimage_mask_SDI")
ADD_EXECUTABLE(${CurrentExe}  "subimage_mask.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_SDI")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

IF(ITKReview_FOUND LESS 0)
  MESSAGE(STATUS "Programs dependent on ITKReview will not be compiled (requires ITK with ITKReview ON)")
ELSE(ITKReview_FOUND LESS 0)
  MESSAGE(STATUS "Programs dependent on ITKReview will be compiled!")

SET(CurrentExe "watershed_morph_f32_UI32")
ADD_EXECUTABLE(${CurrentExe}  "watershed_morph.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT -DUSE_UI32")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "watershed_morph_f32_UI64")
ADD_EXECUTABLE(${CurrentExe}  "watershed_morph.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "watershed_morph_marker_UI32")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "watershed_morph_nX")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "watershed_morph_nX+border-thr")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "watershed_morph_nX_marker_float")
ADD_EXECUTABLE(${CurrentExe}  "watershed_morph_nX_marker.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT" LINK_FLAGS "-lprocps")#ADD_DEFINITIONS(-DUSE_FLOAT) is not target specific!
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "watershed_morph_nX_marker_double")
ADD_EXECUTABLE(${CurrentExe}  "watershed_morph_nX_marker.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES LINK_FLAGS "-lprocps")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "watershed_morph+growing_border")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "watershed_morph_min-max_nX")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

ENDIF(ITKReview_FOUND LESS 0)

SET(CurrentExe "gradient_mag_f32")
ADD_EXECUTABLE(${CurrentExe}  "gradient_mag.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "gradient_mag_f64")
ADD_EXECUTABLE(${CurrentExe}  "gradient_mag.cxx")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "gradient_mag_vec_f32")
ADD_EXECUTABLE(${CurrentExe}  "gradient_mag_vec.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "gradient_mag_vec_f64")
ADD_EXECUTABLE(${CurrentExe}  "gradient_mag_vec.cxx")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "distance_map_signed_maurer_f32")
ADD_EXECUTABLE(${CurrentExe}  "distance_map_signed_maurer.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "distance_map_signed_maurer_f64")
ADD_EXECUTABLE(${CurrentExe}  "distance_map_signed_maurer.cxx")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "distance_map_parabolic_f32")
IF(ParabolicMorphology_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with MinimalPathExtraction ON)")
ELSE(ParabolicMorphology_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  "distance_map_parabolic.cxx")
  SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(ParabolicMorphology_FOUND LESS 0)

SET(CurrentExe "bin")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "bin_SDI")
ADD_EXECUTABLE(${CurrentExe}  "bin.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_SDI")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "bin_series")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "bin_series_SDI")
ADD_EXECUTABLE(${CurrentExe}  "bin_series.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_SDI")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "label_connected_components")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "h-min")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "vec2RGB")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "RGB2vec")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "rgb2RGB")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "rgba2RGBA")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "pad_constant")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "std-mean_ROI_SBS")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "shift-scale_ROI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "thresh-glob")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "thresh-glob_SDI")
ADD_EXECUTABLE(${CurrentExe}  "thresh-glob.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_SDI")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "thresh-erase")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "thresh-erase_SDI")
ADD_EXECUTABLE(${CurrentExe}  "thresh-glob.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_SDI")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "stat")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

## according to the docs itkStatisticsImageFilter cannot stream
## but StreamingStatisticsImageFilter can: https://github.com/blowekamp/itkStreamingSinc/
SET(CurrentExe "stat_SDI")
IF(StreamingSinc_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with Module_StreamingSinc ON)")
ELSE()
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)
MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF()

SET(CurrentExe "label_stack_SBS")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "label_stack")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "LabelShiftImageFilter")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "median")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "median_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "gauss+rescale")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "gauss_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "close_bin_morph")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "close_bin_morph_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "open_bin_morph")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "open_bin_morph_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "erode_bin_morph")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "mask")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "mask-negated")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "resample")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

## itkResampleImageFilter seems to need the whole input no matter what interpolator is used: see GenerateInputRequestedRegion() in Modules/Filtering/ImageGrid/include/itkResampleImageFilter.hxx
##       "So we do the easy thing and request the entire input image."
# SET(CurrentExe "resample_SDI")
# ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
# TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
# INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "resample-iso")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "fill_holes")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "analyse_binary")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "analyse_labels")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "label_stats")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "vessel_enhancing_diffusion")
IF(LesionSizingToolkit_FOUND LESS 0)# -1 if not found!
  #MESSAGE(SEND_ERROR "${CurrentExe} requires ITK to be compiled with LesionSizingToolkit ON")
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with LesionSizingToolkit ON)")
ELSE(LesionSizingToolkit_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(LesionSizingToolkit_FOUND LESS 0)

SET(CurrentExe "change_header")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "min-path_seg_f32")
IF(MinimalPathExtraction_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with MinimalPathExtraction ON)")
ELSE(MinimalPathExtraction_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  "min-path_seg.cxx" "${PROJECT_SOURCE_DIR}/filter/internal/itkConstGradientDescentOptimizer.cxx" "${PROJECT_SOURCE_DIR}/filter/internal/itkInverseGradientDescentOptimizer.cxx")
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  TARGET_INCLUDE_DIRECTORIES(${CurrentExe} PRIVATE "${PROJECT_SOURCE_DIR}/filter/internal/")
  SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(MinimalPathExtraction_FOUND LESS 0)

SET(CurrentExe "skeletonize_SBS")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "skeletonize_3D")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "erode-dilate_dm_f32")
ADD_EXECUTABLE(${CurrentExe}  "erode-dilate_dm.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "open-close_dm_f32")
ADD_EXECUTABLE(${CurrentExe}  "open-close_dm.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "paste_image")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "dilate_parabolic_f32")
IF(ParabolicMorphology_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with MinimalPathExtraction ON)")
ELSE(ParabolicMorphology_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  "dilate_parabolic.cxx")
  SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(ParabolicMorphology_FOUND LESS 0)

SET(CurrentExe "erode_parabolic_f32")
IF(ParabolicMorphology_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with MinimalPathExtraction ON)")
ELSE(ParabolicMorphology_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  "erode_parabolic.cxx")
  SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(ParabolicMorphology_FOUND LESS 0)

SET(CurrentExe "open_parabolic_f32")
IF(ParabolicMorphology_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with MinimalPathExtraction ON)")
ELSE(ParabolicMorphology_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  "open_parabolic.cxx")
  SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(ParabolicMorphology_FOUND LESS 0)

SET(CurrentExe "close_parabolic_f32")
IF(ParabolicMorphology_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with MinimalPathExtraction ON)")
ELSE(ParabolicMorphology_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  "close_parabolic.cxx")
  SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(ParabolicMorphology_FOUND LESS 0)

SET(CurrentExe "fast-marching_f32")
ADD_EXECUTABLE(${CurrentExe}  "fast-marching.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "erode_labels")
IF(LabelErodeDilate_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with LabelErodeDilate ON)")
ELSE(LabelErodeDilate_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(LabelErodeDilate_FOUND LESS 0)

SET(CurrentExe "open_bin_para")
IF(ParabolicMorphology_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with MinimalPathExtraction ON)")
ELSE(ParabolicMorphology_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  "X_bin_para.cxx")
  target_compile_definitions(${CurrentExe} PRIVATE FILTER='O') # https://stackoverflow.com/questions/9017573/define-preprocessor-macro-through-cmake#32255212
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(ParabolicMorphology_FOUND LESS 0)

SET(CurrentExe "close_bin_para")
IF(ParabolicMorphology_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with MinimalPathExtraction ON)")
ELSE(ParabolicMorphology_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  "X_bin_para.cxx")
  target_compile_definitions(${CurrentExe} PRIVATE FILTER='C')
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(ParabolicMorphology_FOUND LESS 0)

SET(CurrentExe "erode_bin_para")
IF(ParabolicMorphology_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with MinimalPathExtraction ON)")
ELSE(ParabolicMorphology_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  "X_bin_para.cxx")
  target_compile_definitions(${CurrentExe} PRIVATE FILTER='E')
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(ParabolicMorphology_FOUND LESS 0)

SET(CurrentExe "dilate_bin_para")
IF(ParabolicMorphology_FOUND LESS 0)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with MinimalPathExtraction ON)")
ELSE(ParabolicMorphology_FOUND LESS 0)
  ADD_EXECUTABLE(${CurrentExe}  "X_bin_para.cxx")
  target_compile_definitions(${CurrentExe} PRIVATE FILTER='D')
  TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
  INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
  MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF(ParabolicMorphology_FOUND LESS 0)

SET(CurrentExe "anti-alias_f32")
ADD_EXECUTABLE(${CurrentExe}  "anti-alias.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "max")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "mean_GPU")
IF(NOT ITKGPUSmoothing_LOADED)# -1 if not found!
  MESSAGE(STATUS "${CurrentExe} will not be compiled (requires ITK with ITK_USE_GPU ON)")
ELSE()
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
MESSAGE(STATUS "${CurrentExe} will be compiled!")
ENDIF()

SET(CurrentExe "count_neighbours")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "vo2ve")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "prune_ends")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "label_uncertainty_float")
ADD_EXECUTABLE(${CurrentExe}  "label_uncertainty.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT" LINK_FLAGS "-lprocps")#ADD_DEFINITIONS(-DUSE_FLOAT) is not target specific!
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "label_uncertainty_double")
ADD_EXECUTABLE(${CurrentExe}  "label_uncertainty.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES LINK_FLAGS "-lprocps")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "toUInt8")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "toUInt16")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "slice")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "add_const")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "add")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "subtract")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "tile")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "dilate_bin")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "anisoDiff-grad_f32")
ADD_EXECUTABLE(${CurrentExe}  "anisoDiff-grad.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "anisoDiff-curv_f32")
ADD_EXECUTABLE(${CurrentExe}  "anisoDiff-curv.cxx")
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DUSE_FLOAT")
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "mean")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "keepNobj")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "open_label-shape")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "hist")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "label_overlap")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "label_measures")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "segmentLung")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "project")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "channel_extract")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
