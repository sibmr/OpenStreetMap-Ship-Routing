# include <iostream>
# include <vector>
# include <array>
# include <math.h>

# include "shortestPathUtils.cpp"
# include "Dijkstra.cpp"

/**
 * @brief contraction hierachies preprocessing
 * 
 * @param array takes adjacency array 
 * @param percentage percentage of how many nodes are part of uncontracted core
 */
void contract(AdjacencyArray &array, double percentage){
    std::cout << "--- BEGIN CONTRACH HIERACHIES ---" << std::endl;
    std::cout << "initial number of edges: " << array.edges.size() << std::endl;





    AdjacencyArray workArray (array);

    while (workArray.nodes.size() > (1 - percentage) * array.nodes.size()){   

        // SELECT INDEPENDENT SET OF NODES  C <= V


        // CREATE SET OF OF SHORTCUTS 


        // REMOVE DUPLICATE SHORTCUTS


        // ASSIGN CURRENT RANK TO NODES


        // MARK ADJACENT EDGES FOR REMOVAL 


        // BUILD NEW WORKARRAY G' = (V \ C, E u S \ A)



    }

    // BUILD FINAL GRAPH G'' = (V; E u E')
    
    





    std::cout << "workarray " << workArray.width << std::endl;

    std::cout << "--- END CONTRACH HIERACHIES ---" << std::endl;

}