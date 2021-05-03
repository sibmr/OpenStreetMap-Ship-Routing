
# Examples
add_executable(example_counter example_counter.cpp)
target_link_libraries(example_counter PRIVATE osmpbfreader)
set_target_properties(example_counter PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(example_routing example_routing.cpp)
target_link_libraries(example_routing PRIVATE osmpbfreader)
set_target_properties(example_routing PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(osmrouting main.cpp)
target_link_libraries(osmrouting PRIVATE osmpbfreader)
set_target_properties(osmrouting PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(pointInPolygon pointInPolygon.cpp)
target_link_libraries(osmrouting PRIVATE osmpbfreader)
set_target_properties(osmrouting PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")