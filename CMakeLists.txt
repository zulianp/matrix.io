cmake_minimum_required(VERSION 3.20)

project(matrixio VERSION 0.1.0 LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(MATRIXIO_BUILD_EXTRAS "Build sorting/checks extras library" "${PROJECT_IS_TOP_LEVEL}")
option(MATRIXIO_BUILD_TOOLS "Build driver tools (print_crs, print_array, ...)" "${PROJECT_IS_TOP_LEVEL}")
option(MATRIXIO_ENABLE_PARMETIS "Build ParMETIS decomposition support" OFF)
option(MATRIXIO_INSTALL "Enable install/export rules" "${PROJECT_IS_TOP_LEVEL}")

find_package(MPI REQUIRED COMPONENTS C)

add_library(matrixio STATIC
  matrixio_crs.c
  matrixio_array.c
  matrixio_ndarray.c
  utils.c
  array_dtof.c
  array_ftod.c
)
add_library(matrixio::matrixio ALIAS matrixio)

set_target_properties(matrixio PROPERTIES
  OUTPUT_NAME "matrix.io"
  POSITION_INDEPENDENT_CODE ON
)

target_compile_features(matrixio PUBLIC c_std_99)

target_include_directories(matrixio
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/graphs>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/checks>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sorting>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/graphs>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/checks>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/sorting>
)

target_link_libraries(matrixio PUBLIC MPI::MPI_C)
if(UNIX)
  target_link_libraries(matrixio PUBLIC m)
endif()

if(MATRIXIO_BUILD_EXTRAS)
  find_package(MPI REQUIRED COMPONENTS CXX)

  add_library(matrixio_extras STATIC
    checks/matrixio_checks.cpp
    sorting/sorting.cpp
  )
  add_library(matrixio::extras ALIAS matrixio_extras)

  set_target_properties(matrixio_extras PROPERTIES
    POSITION_INDEPENDENT_CODE ON
  )

  target_compile_features(matrixio_extras PUBLIC cxx_std_14)
  target_link_libraries(matrixio_extras PUBLIC MPI::MPI_CXX)
  target_include_directories(matrixio_extras PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/checks>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sorting>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/checks>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/sorting>
  )
endif()

if(MATRIXIO_ENABLE_PARMETIS)
  find_package(ParMETIS REQUIRED)

  add_library(matrixio_parmetis STATIC graphs/matrixio_parmetis.c)
  add_library(matrixio::parmetis ALIAS matrixio_parmetis)

  set_target_properties(matrixio_parmetis PROPERTIES
    POSITION_INDEPENDENT_CODE ON
  )

  target_link_libraries(matrixio_parmetis PUBLIC matrixio MPI::MPI_C ParMETIS::parmetis)
  target_include_directories(matrixio_parmetis PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/graphs>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/graphs>
  )
endif()

if(MATRIXIO_BUILD_TOOLS)
  add_executable(test drivers/test.c)
  target_link_libraries(test PRIVATE matrixio MPI::MPI_C)

  add_executable(print_crs drivers/print_crs.c)
  target_link_libraries(print_crs PRIVATE matrixio MPI::MPI_C)

  add_executable(print_array drivers/print_array.c)
  target_link_libraries(print_array PRIVATE matrixio MPI::MPI_C)

  if(MATRIXIO_BUILD_EXTRAS)
    add_executable(reorder_crs drivers/reorder_crs.c)
    target_link_libraries(reorder_crs PRIVATE matrixio matrixio_extras MPI::MPI_C)
    set_property(TARGET reorder_crs PROPERTY LINKER_LANGUAGE CXX)

    if(MATRIXIO_ENABLE_PARMETIS)
      add_executable(partition_crs drivers/partition_crs.c)
      target_link_libraries(partition_crs PRIVATE matrixio_extras matrixio_parmetis MPI::MPI_C)
      set_property(TARGET partition_crs PROPERTY LINKER_LANGUAGE CXX)
    endif()
  endif()
endif()

if(MATRIXIO_INSTALL)
  install(
    TARGETS matrixio
    EXPORT matrixioTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  )

  if(MATRIXIO_BUILD_EXTRAS)
    install(
      TARGETS matrixio_extras
      EXPORT matrixioTargets
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
  endif()

  if(MATRIXIO_ENABLE_PARMETIS)
    install(
      TARGETS matrixio_parmetis
      EXPORT matrixioTargets
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
  endif()

  install(
    FILES
      array_dtof.h
      array_ftod.h
      matrixio_array.h
      matrixio_base.h
      matrixio_crs.h
      matrixio_ndarray.h
      utils.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
  install(DIRECTORY graphs/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/graphs" FILES_MATCHING PATTERN "*.h")
  install(DIRECTORY checks/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/checks" FILES_MATCHING PATTERN "*.h")
  install(DIRECTORY sorting/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sorting" FILES_MATCHING PATTERN "*.h")

  install(FILES cmake/FindParMETIS.cmake DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/matrixioConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/matrixioConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  )

  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/matrixioConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
  )

  install(
    EXPORT matrixioTargets
    NAMESPACE matrixio::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  )

  install(
    FILES
      "${CMAKE_CURRENT_BINARY_DIR}/matrixioConfig.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/matrixioConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  )
endif()

