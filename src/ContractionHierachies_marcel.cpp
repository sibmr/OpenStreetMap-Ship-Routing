# include <iostream>
# include <vector>
# include <array>
# include <math.h>
# include <set>
# include <time.h>
#include <chrono>

# include "shortestPathUtils.cpp"
# include "Dijkstra.cpp"



/**
 * @brief generate set of independent nodes. Try random nodes and add them if they for an independent set
 *  
 * @param array 
 * @param adjacentNodes 
 * @param independentSetVector 
 * @param independentSet 
 * @param alreadyContractedNodes 
 * @param totalWaterNodes 
 */
void independentNodeSetRandom(
    AdjacencyArray &array, 
    std::set<uint64_t> &independentSet, 
    uint64_t alreadyContractedNodes,
    uint64_t totalWaterNodes
    ){
    std::vector<bool> adjacentNodes(array.rank.size(), false);
    std::vector<bool> independentSetVector(array.rank.size(), false);
    // use one seedNode on water
    uint64_t seedNode = 0; 
    bool correctNode = true;

    uint64_t limitNumberOfTries = totalWaterNodes - alreadyContractedNodes; // if we can not find any more node efficient we stop early

    uint64_t sizeOfSet = (totalWaterNodes - alreadyContractedNodes) / 6; // only one out of 5 nodes can be in independent set

    uint64_t currSetSize = 0;
    while(currSetSize < sizeOfSet && limitNumberOfTries > 0){

        do{
            seedNode = std::rand() % array.rank.size();
            limitNumberOfTries--;
            correctNode = !(isNodeOnLand(array, seedNode) || independentSetVector.at(seedNode) || adjacentNodes.at(seedNode) || array.rank.at(seedNode) > 0);
            if(limitNumberOfTries == 0){
                break;
            }

        }while (!correctNode);

        if(correctNode){
            independentSetVector.at(seedNode) = true;
            independentSet.insert(seedNode);
            currSetSize++;


            // mark adjacent nodes so they can not be in the independent set
            for(uint64_t currEdgeId = array.offsets.at(seedNode); currEdgeId < array.offsets.at(seedNode+1); currEdgeId++){
                    uint64_t neighborIdx = array.edges.at(currEdgeId);
                    adjacentNodes.at(neighborIdx) = true;
            }
        }

    }



}

/**
 * @brief compare two ChEdges sorted by startNode
 * 
 * @param i1 
 * @param i2 
 * @return true 
 * @return false 
 */
bool compareEdge(ChEdge i1, ChEdge i2)
{
    return (i1.startNode < i2.startNode);
}


/**
 * @brief calculate how many percent of nodes have been already contracted
 * 
 * @param numberOfNodesOnWater 
 * @param numberOfNodesContracted 
 * @return double 
 */
double contractedPercent(uint64_t numberOfNodesOnWater, uint64_t numberOfNodesContracted){
    return (((double) (numberOfNodesOnWater - numberOfNodesContracted)) / ((double) numberOfNodesOnWater) * 100);
}




/**
 * @brief contraction hierachies preprocessing
 * 
 * @param array takes adjacency array 
 * @param percentage percentage of how many nodes are part of uncontracted core
 */
void contract(AdjacencyArray &array, double percentage){
    std::cout << "--- BEGIN CONTRACT HIERACHIES ---" << std::endl;
    std::cout << "initial number of edges: " << array.edges.size() << std::endl;

    uint64_t totalIndependentSetTime = 0;
    uint64_t totalDijkstraTime = 0;
    uint64_t totalRebuildTime = 0;
    uint64_t totalTime = 0;

    auto allStart = std::chrono::high_resolution_clock::now();

    AdjacencyArray workArray (array);

    uint64_t numberOfNodesOnWater = 0;
    uint64_t numberOfNodesContracted = 0;


    std::vector<ChEdge> allShortcuts; 

    // maybe non rel
    std::vector<uint64_t> globalEdges = array.edges;

    for(uint64_t nodeId = 0; nodeId < array.rank.size(); nodeId++){
        if(!isNodeOnLand(array, nodeId)){
            numberOfNodesOnWater++;
        }
    }
    std::cout << "Nodes on water: " << numberOfNodesOnWater << std::endl;
    uint64_t sizeOfUncontractedCore =  (uint64_t) ((1 - percentage) * numberOfNodesOnWater);


    // initialize random seed based on time
    std::srand(time(NULL));

    uint64_t numNodes = workArray.rank.size();

    std::set<uint64_t> independentSet;
    std::set<uint64_t> contractedSet;

    std::vector<ChEdge> roundEdges;

    uint64_t currentRank = 1;

    std::cout << "workarray nodes " << numNodes << std::endl;



    // stop if the uncontracted core is smaller than x percent of all nodes
    while ((numberOfNodesOnWater - numberOfNodesContracted) > sizeOfUncontractedCore){   

        std::cout << contractedPercent(numberOfNodesOnWater, numberOfNodesContracted) <<  " percent "  << (numberOfNodesOnWater - numberOfNodesContracted) << " of " << sizeOfUncontractedCore << std::endl;

        AdjacencyArray workArrayNew;

        // SELECT INDEPENDENT SET OF NODES  C <= V


        independentSet.clear();

        auto independentSetStart = std::chrono::high_resolution_clock::now();

        independentNodeSetRandom(workArray, independentSet, numberOfNodesContracted, numberOfNodesOnWater);

        std::cout << "set generated " <<  independentSet.size()   << "\t contraction of " << contractedPercent(numberOfNodesOnWater, numberOfNodesContracted) << " percent" << std::endl;
        auto independentSetStop = std::chrono::high_resolution_clock::now();

        auto independentSetTime = std::chrono::duration_cast<std::chrono::microseconds>(independentSetStop - independentSetStart).count();
        totalIndependentSetTime += independentSetTime;


        SecondDijkstra dijkstra(workArray);
        //dijkstra.removeContractedNodes();

        // CREATE SET OF OF SHORTCUTS 
        uint64_t dijkstraDistance = UINT64_MAX;
        uint64_t arrayDistance = UINT64_MAX;

        
        contractedSet.clear();
        roundEdges.clear();
        uint64_t counter = 0;

        auto dijkstraStart = std::chrono::high_resolution_clock::now();

        for(uint64_t nodeInIndependentSet : independentSet){

            // use 70% of the independent set (ranoom order)
            if((std::rand() % 100) < 70){
                contractedSet.insert(nodeInIndependentSet);

                for(uint64_t currEdgeId1 = workArray.offsets.at(nodeInIndependentSet); currEdgeId1 < workArray.offsets.at(nodeInIndependentSet+1); currEdgeId1++){
                    uint64_t sourceNeighbor = workArray.edges.at(currEdgeId1);

                    // save half of the iterations since all edges are undirected
                    for(uint64_t currEdgeId2 = workArray.offsets.at(nodeInIndependentSet); currEdgeId2 < currEdgeId1; currEdgeId2++){
                        uint64_t targetNeighbor = workArray.edges.at(currEdgeId2);

                        dijkstraDistance = dijkstra.calculateDist(sourceNeighbor, targetNeighbor);
                        arrayDistance = workArray.distances.at(currEdgeId1) + workArray.distances.at(currEdgeId2);

                        if(dijkstraDistance >= arrayDistance){
                            roundEdges.push_back(ChEdge(sourceNeighbor, targetNeighbor, dijkstraDistance, currEdgeId1, currEdgeId2));
                            roundEdges.push_back(ChEdge(targetNeighbor, sourceNeighbor, dijkstraDistance, currEdgeId2, currEdgeId1));
                        }
                    }
                }
            }
        }
        auto dijkstraStop = std::chrono::high_resolution_clock::now();
        auto dijkstraTiming = std::chrono::duration_cast<std::chrono::microseconds>(dijkstraStop - dijkstraStart).count();
        totalDijkstraTime += dijkstraTiming;

        std::cout << "contracted Edges " << roundEdges.size()  << "\t contraction of " << contractedPercent(numberOfNodesOnWater, numberOfNodesContracted) << " percent" << std::endl;

        auto tempGraphBuildStart = std::chrono::high_resolution_clock::now();

        // sort edges in this rank by start node
        std::sort(roundEdges.begin(), roundEdges.end(), compareEdge);

        std::vector<bool> removeShortcut (roundEdges.size(), false);

        if(roundEdges.size() > 0){

            // REMOVE DUPLICATE SHORTCUTS
            uint64_t currentLastStartnode = roundEdges.at(0).startNode;
            std::set<uint64_t> currentLastTargetnodes;
            currentLastTargetnodes.insert(roundEdges.at(0).targetNode);


            for (uint64_t shortcutId = 1; shortcutId < roundEdges.size(); shortcutId++){
                if(currentLastStartnode == roundEdges.at(shortcutId).startNode){
                    if(currentLastTargetnodes.find(roundEdges.at(shortcutId).targetNode) != currentLastTargetnodes.end()){
                        // start and endnode is already in set remove shortcut
                        removeShortcut.at(shortcutId) = true;
                    }else{
                        // new endNode found add to targetSet
                        currentLastTargetnodes.insert(roundEdges.at(shortcutId).targetNode);
                    }
                    
                }else{
                    // new startNode found clear endNode set and set new startNode
                    currentLastStartnode = roundEdges.at(shortcutId).startNode;
                    currentLastTargetnodes.clear();
                    currentLastTargetnodes.insert(roundEdges.at(shortcutId).targetNode);
                }
            }
        }

        std::cout << "shortcuts marked" << std::endl;

        // ASSIGN CURRENT RANK TO NODES

        for(uint64_t nodeInContractedSet : contractedSet){
            workArray.rank.at(nodeInContractedSet) = currentRank;
        }
        currentRank++;
        numberOfNodesContracted += contractedSet.size();

        std::cout << "ranks updated" << std::endl;


        // MARK ADJACENT EDGES FOR REMOVAL 

        std::vector<bool> adjacentEdges (workArray.edges.size(), false);

        for(uint64_t nodeInContractedSet : contractedSet){
            // mark all edges with (v,*)
            for(uint64_t currEdgeId = workArray.offsets.at(nodeInContractedSet); 
                currEdgeId < workArray.offsets.at(nodeInContractedSet+1); 
                currEdgeId++){
                    adjacentEdges.at(currEdgeId) = true;
            }
            
            // mark all edges with (v,*)
            for(uint64_t currEdgeId = 0; currEdgeId < workArray.edges.size(); currEdgeId++){
                if(workArray.edges.at(currEdgeId) == nodeInContractedSet){
                    adjacentEdges.at(currEdgeId) = true;
                }
            }
        }
        std::cout << "adj edges marked" << std::endl;


        // BUILD NEW WORKARRAY G' = (V \ C, E u S \ A)

        // V \ C are all verticies with rank == 0


        workArrayNew.longLow = workArray.longLow;     
        workArrayNew.latLow = workArray.latLow;
        workArrayNew.longHigh = workArray.longHigh;   
        workArrayNew.latHigh = workArray.latHigh;
        workArrayNew.width = workArray.width;         
        workArrayNew.height = workArray.height;

        workArrayNew.nodes = workArray.nodes;    
        workArrayNew.rank = workArray.rank;

        uint64_t roundEdgeIndex = 0;
        uint64_t oldEdgesIndex = 0;
        uint64_t preventedShortcuts = 0;
        ChEdge currEdge;

        workArrayNew.offsets.push_back(0);
        uint64_t currentOffset = 0;
        for(uint64_t nodeId = 0; nodeId < workArray.rank.size(); nodeId++){
            
            // if there is no shortcut found in this iteration we skip this step
            if(roundEdges.size() > 0){
                bool whileAlive = true;
                while(whileAlive && roundEdges.at(roundEdgeIndex).startNode == nodeId){
                    // only add shortcut if there is no identical other one
                    if(!removeShortcut.at(roundEdgeIndex)){
                        currEdge = roundEdges.at(roundEdgeIndex);
                        workArrayNew.edges.push_back(currEdge.targetNode);
                        workArrayNew.distances.push_back(currEdge.distance);

                        allShortcuts.push_back(roundEdges.at(roundEdgeIndex));
                        currentOffset++;
                    }else{
                        preventedShortcuts++;
                    }
                    roundEdgeIndex++;
                    if(roundEdgeIndex == roundEdges.size()){
                        roundEdgeIndex--;
                        whileAlive = false;
                    }

                }
            }
            // add all edges not adjacent to the contracted core
            for(uint64_t currEdge = workArray.offsets.at(nodeId); currEdge < workArray.offsets.at(nodeId +1); currEdge++){
                if(!adjacentEdges.at(currEdge)){
                    workArrayNew.edges.push_back(workArray.edges.at(currEdge));
                    workArrayNew.distances.push_back(workArray.distances.at(currEdge));
                    currentOffset++;
                }

            }

            workArrayNew.offsets.push_back(currentOffset);
        }
        auto tempGraphBuildStop = std::chrono::high_resolution_clock::now();
        auto tempGraphBuildTiming = std::chrono::duration_cast<std::chrono::microseconds>(tempGraphBuildStop - tempGraphBuildStart).count();
        totalRebuildTime += tempGraphBuildTiming;

        totalTime = std::chrono::duration_cast<std::chrono::microseconds>(tempGraphBuildStop - allStart).count();

        std::cout << "intermediate shortcuts for rank " << currentRank-1 << " has been calculated" << std::endl;
        std::cout << "new array has : " << workArrayNew.offsets.at(workArrayNew.offsets.size() - 1)  << " edges old array had " <<  workArray.offsets.at(workArray.offsets.size() - 1) << " edges"<< std::endl;
        std::cout << "prevented shortcuts: " << preventedShortcuts <<   " from " << roundEdges.size() << std::endl;
        std::cout << "independetSet took " << independentSetTime / 1000 << " ms" 
                << " dijkstra took " << dijkstraTiming / 1000000 << " s" 
                << " graph Build took "  << tempGraphBuildTiming / 1000 << " ms" 
                << " overall " << totalTime / 1000000 << " s" << std::endl;

        workArray = workArrayNew;

        std::cout << "--------------------------------------------------------------" << std::endl;






    }

    // BUILD FINAL GRAPH G'' = (V; E u E')

    std::cout << "iterations finished" << std::endl;
    
    
    std::sort(allShortcuts.begin(), allShortcuts.end(), compareEdge);

    AdjacencyArray finalGraph; 

    finalGraph.longLow = workArray.longLow;     
    finalGraph.latLow = workArray.latLow;
    finalGraph.longHigh = workArray.longHigh;   
    finalGraph.latHigh = workArray.latHigh;
    finalGraph.width = workArray.width;         
    finalGraph.height = workArray.height;

    finalGraph.nodes = workArray.nodes;    
    finalGraph.rank = workArray.rank;

    for(uint64_t currNode = 0; currNode < finalGraph.rank.size(); currNode++){
        if(finalGraph.rank.at(currNode) == 0){
            finalGraph.rank.at(currNode) = currentRank;
        }
    }


    std::cout << "end" << std::endl;

    uint64_t allShortcutId = 0;
    uint64_t oldEdgesIndex = 0;
    uint64_t preventedShortcuts = 0;

    finalGraph.offsets.push_back(0);
    uint64_t currentOffset = 0;
    for(uint64_t nodeId = 0; nodeId < workArray.rank.size(); nodeId++){
        

        // add all shortcuts
        bool whileAlive = true;
        if(allShortcuts.size() > 0){
            while(allShortcuts.at(allShortcutId).startNode == nodeId && whileAlive){

                // only add shortcut if there is no identical other one
                finalGraph.edges.push_back(allShortcuts.at(allShortcutId).targetNode);
                finalGraph.distances.push_back(allShortcuts.at(allShortcutId).distance);
                currentOffset++;
                allShortcutId++;
                if(allShortcutId == allShortcuts.size()){
                    allShortcutId--;
                    whileAlive = false;
                }

            }
        }

        // add all initial edges
        for(uint64_t currEdge = array.offsets.at(nodeId); currEdge < array.offsets.at(nodeId +1); currEdge++){
            finalGraph.edges.push_back(array.edges.at(currEdge));
            finalGraph.distances.push_back(array.distances.at(currEdge));
            currentOffset++;

        }

        finalGraph.offsets.push_back(currentOffset);
    }

    std::cout << "Independent Set took " << totalIndependentSetTime / 1000 << "ms in total"  << std::endl;
    std::cout << "Dijkstra took " << totalDijkstraTime / 1000000 << "s in total"  << std::endl;
    std::cout << "Graph rebuild took " << totalRebuildTime / 1000000 << "s in total"  << std::endl;
    std::cout << "In total took " << totalTime / 1000000 << "s in total"  << std::endl;

    std::cout << "finalGraph has " << finalGraph.edges.size() << " edges from initially "  << array.edges.size() << " edges" << std::endl;

    array = finalGraph;

    std::cout << "--- END CONTRACT HIERACHIES ---" << std::endl;

}