# include <iostream>
# include <vector>
# include <array>
# include <math.h>
# include <set>

# include "shortestPathUtils.cpp"
# include "Dijkstra.cpp"



void independentNodeSetRandom(
    AdjacencyArray &array, 
    std::vector<bool> &adjacentNodes, 
    std::vector<bool> &independentSetVector, 
    std::set<uint64_t> &independentSet, 
    uint64_t sizeOfSet
    ){
    // use one seedNode on water
    uint64_t seedNode = 0; 

    uint64_t currSetSize = 0;
    while(currSetSize < sizeOfSet){

        do{
            seedNode = std::rand() % array.rank.size();
        }while (isNodeOnLand(array, seedNode) || independentSetVector.at(seedNode) || adjacentNodes.at(seedNode));

        independentSetVector.at(seedNode) = true;
        independentSet.insert(seedNode);
        currSetSize++;

        std::cout << currSetSize  <<  " " << sizeOfSet << std::endl;

        for(uint64_t currEdgeId = array.offsets.at(seedNode); currEdgeId < array.offsets.at(seedNode+1); currEdgeId++){
                uint64_t neighborIdx = array.edges.at(currEdgeId);
                adjacentNodes.at(neighborIdx) = true;
        }

    }



}



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

    // initialize random seed based on size
    std::srand(time(NULL));

    uint64_t numNodes = workArray.rank.size();

    std::set<uint64_t> independentSet;

    SecondDijkstra dijkstra(array);
    
    std::cout << "workarray nodes " << numNodes << std::endl;


    std::cout << workArray.nodes.size() << "\t" << (1 - percentage) * numNodes << std::endl;

    while (workArray.nodes.size() > (1 - percentage) * numNodes){   

        // SELECT INDEPENDENT SET OF NODES  C <= V


        std::vector<bool> adjacentNodes(array.width * array.height, false);
        std::vector<bool> independentSetVector(array.width * array.height, false);
        independentSet.clear();

        independentNodeSetRandom(workArray, adjacentNodes, independentSetVector, independentSet, (array.width * array.height) / 100);



        // VISUALIZE INDEPENDENT SET 
        //for(int i = 0; i < array.height; i++){
        //    for(int j = 0; j < array.width; j++){
        //        std::cout << independentSet.at(array.height * i + j);
        //    }
        //    std::cout << "--" <<  std::endl;

        //}

        //std::cout << "--" <<  std::endl;

        //for(int i = 0; i < array.height; i++){
        //    for(int j = 0; j < array.width; j++){
        //        std::cout << adjacentNodes.at(array.height * i + j);
        //    }
        //    std::cout << "--" <<  std::endl;

        //}


        // CREATE SET OF OF SHORTCUTS 
        uint64_t dijkstraDistance = UINT64_MAX;

        for(uint64_t nodeInIndependentSet : independentSet){
            //std::cout << nodeInIndependentSet << " ";
            for(uint64_t currEdgeId1 = array.offsets.at(nodeInIndependentSet); currEdgeId1 < array.offsets.at(nodeInIndependentSet+1); currEdgeId1++){
                uint64_t sourceNeighbor = array.edges.at(currEdgeId1);
                for(uint64_t currEdgeId2 = array.offsets.at(nodeInIndependentSet); currEdgeId2 < array.offsets.at(nodeInIndependentSet+1); currEdgeId2++){
                    uint64_t targetNeighbor = array.edges.at(currEdgeId2);
                    if(sourceNeighbor != targetNeighbor){

                        //dijkstra.reset();
                        dijkstraDistance = dijkstra.calculateDist(sourceNeighbor, targetNeighbor);



                        if(dijkstraDistance <= array.distances.at(currEdgeId1) + array.distances.at(currEdgeId2)){
                            // add shortcut edge
                            std::cout << sourceNeighbor << "--" << targetNeighbor << "\t" << dijkstraDistance  << "--"  << array.distances.at(currEdgeId1)  + array.distances.at(currEdgeId2) << std::endl;

                        }

                    }

                }
                

            }
        }





        return;
        // REMOVE DUPLICATE SHORTCUTS


        // ASSIGN CURRENT RANK TO NODES


        // MARK ADJACENT EDGES FOR REMOVAL 


        // BUILD NEW WORKARRAY G' = (V \ C, E u S \ A)



    }

    // BUILD FINAL GRAPH G'' = (V; E u E')
    
    





    std::cout << "workarray " << workArray.width << std::endl;

    std::cout << "--- END CONTRACH HIERACHIES ---" << std::endl;

}