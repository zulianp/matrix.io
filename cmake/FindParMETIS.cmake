# Minimal Find module for ParMETIS.
#
# Provides:
#   ParMETIS_FOUND
#   ParMETIS_INCLUDE_DIR
#   ParMETIS_LIBRARY
#   ParMETIS::parmetis (imported target)

include(FindPackageHandleStandardArgs)

set(_ParMETIS_HINTS "")
if(DEFINED ENV{PARMETIS_DIR})
  list(APPEND _ParMETIS_HINTS "$ENV{PARMETIS_DIR}")
endif()
if(DEFINED PARMETIS_DIR)
  list(APPEND _ParMETIS_HINTS "${PARMETIS_DIR}")
endif()

find_path(ParMETIS_INCLUDE_DIR
  NAMES parmetis.h
  HINTS ${_ParMETIS_HINTS}
  PATH_SUFFIXES include
)

find_library(ParMETIS_LIBRARY
  NAMES parmetis
  HINTS ${_ParMETIS_HINTS}
  PATH_SUFFIXES lib lib64
)

find_package_handle_standard_args(ParMETIS
  REQUIRED_VARS ParMETIS_INCLUDE_DIR ParMETIS_LIBRARY
)

if(ParMETIS_FOUND AND NOT TARGET ParMETIS::parmetis)
  add_library(ParMETIS::parmetis UNKNOWN IMPORTED)
  set_target_properties(ParMETIS::parmetis PROPERTIES
    IMPORTED_LOCATION "${ParMETIS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${ParMETIS_INCLUDE_DIR}"
  )
endif()

