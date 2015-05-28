PROJECT(ITK-CLIs)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

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
      INCLUDE( ${USE_ITK_FILE} )
ENDIF( ITK_FOUND )
# FIND_PACKAGE ( VTK )
# IF ( VTK_FOUND )
#       INCLUDE( ${USE_VTK_FILE} )
# ENDIF( VTK_FOUND )

INCLUDE_DIRECTORIES(
  # $ENV{HOME}/itk/parabolicMorphology/
  # $ENV{HOME}/itk/ITKbinaryThinningImageFilter3D/Source/
  # $ENV{HOME}/itk/itkWindowImageFilter/
  # $ENV{HOME}/itk/itkLabelSetErodeDilate_IJ_886/
  )

SET(Libraries
  ${ITK_LIBRARIES}
)


SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

## check if specific modules are enabled
LIST(FIND ITK_MODULES_ENABLED ITKReview ITKReview_FOUND)
LIST(FIND ITK_MODULES_ENABLED LesionSizingToolkit LesionSizingToolkit_FOUND)


SET(CurrentExe "file_converter")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "file_converter_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "file-series_reader_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
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

SET(CurrentExe "extract_subimage")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "extract_subimage_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

IF(ITKReview_FOUND LESS 0)
  MESSAGE(STATUS "Programs dependent on ITKReview will not be compiled (requires ITK with ITKReview ON)")
ELSE(ITKReview_FOUND LESS 0)
  MESSAGE(STATUS "Programs dependent on ITKReview will be compiled!")

SET(CurrentExe "watershed_morph")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
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

ENDIF(ITKReview_FOUND LESS 0)

SET(CurrentExe "gradient_mag")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "distance_map_signed_maurer")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "bin")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
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
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "stat")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

## according to the docs itkStatisticsImageFilter cannot stream
## but StreamingStatisticsImageFilter can: https://github.com/blowekamp/itkStreamingSinc/
SET(CurrentExe "stat_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

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

SET(CurrentExe "closing_bin_morph")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "closing_bin_morph_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "opening_bin_morph")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "opening_bin_morph_SDI")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "mask")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "resample")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)

SET(CurrentExe "resample_SDI")
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

SET(CurrentExe "vessel_enhancing_diffusion.1")
ADD_EXECUTABLE(${CurrentExe}  ${CurrentExe})
SET_TARGET_PROPERTIES(${CurrentExe} PROPERTIES COMPILE_FLAGS "-DINTERMEDIATE_OUTPUTS")
TARGET_LINK_LIBRARIES(${CurrentExe}  ${Libraries})
INSTALL(TARGETS  ${CurrentExe}  RUNTIME  DESTINATION  bin)
