add_executable(test_single_query test_single_query.cpp)
set_target_properties(test_single_query PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(test test.cpp)
set_target_properties(test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(OpenMP REQUIRED)
add_executable(CH_preprocessing CH_preprocessing.cpp)
target_link_libraries(CH_preprocessing PRIVATE OpenMP::OpenMP_CXX)
set_target_properties(CH_preprocessing PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(server server.cpp)
target_link_libraries(server httplib::httplib)
set_target_properties(server PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

#/add_executable(dijkstra Dijkstra.cpp)
#/target_link_libraries(dijkstra PRIVATE osmpbfreader)
#/set_target_properties(dijkstra PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(benchmarkPathAlgorithm benchmarkPathAlgorithm.cpp)
set_target_properties(benchmarkPathAlgorithm PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(buildGraph buildGraph.cpp)
set_target_properties(buildGraph PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(coastlineExtraction costlineExtraction.cpp)
target_link_libraries(coastlineExtraction PRIVATE osmpbfreader)
set_target_properties(coastlineExtraction PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(OpenMP REQUIRED)
add_executable(pointInPolygon pointInPolygon.cpp)
target_link_libraries(pointInPolygon PRIVATE OpenMP::OpenMP_CXX)
set_target_properties(pointInPolygon PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")