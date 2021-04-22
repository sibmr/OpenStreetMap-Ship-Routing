find_path(osmpbf_INCLUDE_DIRS
  NAMES osmpbf/osmpbf.h
)
mark_as_advanced(osmpbf_INCLUDE_DIRS)

find_library(osmpbf_LIBRARIES
  NAMES osmpbf
)
mark_as_advanced(osmpbf_LIBRARIES)


if(NOT osmpbf_INCLUDE_DIRS)
  message(STATUS "Could NOT find osmpbf/osmpbf.h")
endif()
if(NOT osmpbf_LIBRARIES)
  message(STATUS "Could NOT find osmpbf library")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(osmpbf DEFAULT_MSG osmpbf_INCLUDE_DIRS osmpbf_LIBRARIES)
