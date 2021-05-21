
add_executable(server server.cpp)
#target_link_libraries(server PRIVATE osmpbfreader)
target_link_libraries(server httplib::httplib)
set_target_properties(server PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(dijkstra dijkstra.cpp)
target_link_libraries(dijkstra PRIVATE osmpbfreader)
set_target_properties(dijkstra PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(buildGraph buildGraph.cpp)
target_link_libraries(buildGraph PRIVATE osmpbfreader)
set_target_properties(buildGraph PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(coastlineExtraction costlineExtraction.cpp)
target_link_libraries(coastlineExtraction PRIVATE osmpbfreader)
set_target_properties(coastlineExtraction PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(OpenMP REQUIRED)
add_executable(pointInPolygon pointInPolygon.cpp)
target_link_libraries(pointInPolygon PRIVATE osmpbfreader)
target_link_libraries(pointInPolygon PRIVATE OpenMP::OpenMP_CXX)
set_target_properties(pointInPolygon PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")