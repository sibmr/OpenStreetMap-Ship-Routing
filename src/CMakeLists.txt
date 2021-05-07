
add_executable(osmrouting main.cpp)
target_link_libraries(osmrouting PRIVATE osmpbfreader)
set_target_properties(osmrouting PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(pointInPolygon pointInPolygon.cpp)
target_link_libraries(pointInPolygon PRIVATE osmpbfreader)
set_target_properties(pointInPolygon PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")