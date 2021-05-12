
add_executable(coastlineExtraction costlineExtraction.cpp)
target_link_libraries(coastlineExtraction PRIVATE osmpbfreader)
set_target_properties(coastlineExtraction PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(pointInPolygon pointInPolygon.cpp)
target_link_libraries(pointInPolygon PRIVATE osmpbfreader)
set_target_properties(pointInPolygon PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")