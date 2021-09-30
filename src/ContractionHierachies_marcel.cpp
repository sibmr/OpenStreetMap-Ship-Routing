# include <iostream>
# include <vector>
# include <array>
# include <math.h>

# include "shortestPathUtils.cpp"

/**
 * @brief contraction hierachies preprocessing
 * 
 * @param array takes adjacency array 
 * @param percentage percentage of how many nodes are part of uncontracted core
 */
void contract(AdjacencyArray &array, double percentage){
    std::cout << "contract method" << std::endl;

    std::cout << array.nodes.size() << std::endl;
    std::cout << array.rank.size() << std::endl;
    std::cout << array.offsets.at(array.offsets.size()-1) << std::endl;
    std::cout << array.nodes.at(0) << " ";
    std::cout << array.nodes.at(array.offsets.size()-1) << std::endl;


}