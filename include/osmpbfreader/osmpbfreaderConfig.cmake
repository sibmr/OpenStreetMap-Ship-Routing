get_filename_component(osmpbfreader_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

if(NOT TARGET osmpbfreader::osmpbfreader)
  set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${osmpbfreader_CMAKE_DIR}/modules")
  find_package(osmpbf REQUIRED)
  find_package(Protobuf REQUIRED)
  include("${osmpbfreader_CMAKE_DIR}/osmpbfreaderTargets.cmake")
endif()
