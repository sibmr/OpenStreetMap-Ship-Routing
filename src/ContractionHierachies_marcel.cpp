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
    uint64_t alreadyContractedNodes,
    uint64_t totalLandNodes
    ){
    // use one seedNode on water
    uint64_t seedNode = 0; 

    uint64_t limitNumberOfTries = totalLandNodes - alreadyContractedNodes; // if we can not find any more node efficient we stop early

    uint64_t sizeOfSet = (totalLandNodes - alreadyContractedNodes) / 6; // only one out of 5 nodes can be in independent set

    uint64_t currSetSize = 0;
    while(currSetSize < sizeOfSet && limitNumberOfTries > 0){

        do{
            seedNode = std::rand() % array.rank.size();
            limitNumberOfTries--;
            //std::cout << currSetSize  << " of " << sizeOfSet <<  "\t " << limitNumberOfTries << std::endl;
        }while (isNodeOnLand(array, seedNode) || independentSetVector.at(seedNode) || adjacentNodes.at(seedNode) || array.rank.at(seedNode) > 0);

        independentSetVector.at(seedNode) = true;
        independentSet.insert(seedNode);
        currSetSize++;


        for(uint64_t currEdgeId = array.offsets.at(seedNode); currEdgeId < array.offsets.at(seedNode+1); currEdgeId++){
                uint64_t neighborIdx = array.edges.at(currEdgeId);
                adjacentNodes.at(neighborIdx) = true;
        }

    }



}

bool compareEdge(ChEdge i1, ChEdge i2)
{
    return (i1.startNode < i2.startNode);
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


    AdjacencyArray workArray (array);

    uint64_t numberOfNodesOnLand = 0;
    uint64_t numberOfNodesContracted = 0;

    std::vector<ChEdge> allShortcuts; 

    for(uint64_t nodeId = 0; nodeId < array.rank.size(); nodeId++){
        if(isNodeOnLand(array, nodeId)){
            numberOfNodesOnLand++;
        }
    }
    std::cout << "Nodes on land: " << numberOfNodesOnLand << std::endl;


    // initialize random seed based on size
    std::srand(time(NULL));

    uint64_t numNodes = workArray.rank.size();

    std::set<uint64_t> independentSet;
    std::set<uint64_t> contractedSet;

    std::vector<ChEdge> roundEdges;

    uint64_t currentRank = 1;

    SecondDijkstra dijkstra(workArray);
    
    std::cout << "workarray nodes " << numNodes << std::endl;



    // stop if the uncontracted core is smaller than x percent of all nodes

    std::cout << (numberOfNodesOnLand - numberOfNodesContracted) << " " << (1-percentage) * numberOfNodesOnLand << std::endl;

    while ((numberOfNodesOnLand - numberOfNodesContracted) > (1 - percentage) * numberOfNodesOnLand){   

        AdjacencyArray workArrayNew;

        // SELECT INDEPENDENT SET OF NODES  C <= V


        std::vector<bool> adjacentNodes(workArray.rank.size(), false);
        std::vector<bool> independentSetVector(workArray.rank.size(), false);
        independentSet.clear();

        independentNodeSetRandom(workArray, adjacentNodes, independentSetVector, independentSet, numberOfNodesContracted, numberOfNodesOnLand);



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
        uint64_t arrayDistance = UINT64_MAX;


        
        roundEdges.clear();

        for(uint64_t nodeInIndependentSet : independentSet){
            //std::cout << nodeInIndependentSet << " ";

            if((std::rand() % 100) > 70){
                contractedSet.insert(nodeInIndependentSet);

                for(uint64_t currEdgeId1 = workArray.offsets.at(nodeInIndependentSet); currEdgeId1 < workArray.offsets.at(nodeInIndependentSet+1); currEdgeId1++){
                    uint64_t sourceNeighbor = workArray.edges.at(currEdgeId1);
                    for(uint64_t currEdgeId2 = workArray.offsets.at(nodeInIndependentSet); currEdgeId2 < workArray.offsets.at(nodeInIndependentSet+1); currEdgeId2++){
                        uint64_t targetNeighbor = workArray.edges.at(currEdgeId2);
                        if(sourceNeighbor != targetNeighbor){

                            dijkstraDistance = dijkstra.calculateDist(sourceNeighbor, targetNeighbor);
                            arrayDistance = workArray.distances.at(currEdgeId1) + workArray.distances.at(currEdgeId2);

                            if(dijkstraDistance <= arrayDistance){
                              // add shortcut edge
                              roundEdges.push_back(ChEdge(sourceNeighbor, targetNeighbor, 1000, currEdgeId1, currEdgeId2));
                            }
                        }
                    }
                }
            }
        }

        std::cout << "contracted Edges " << roundEdges.size() << std::endl;


        std::sort(roundEdges.begin(), roundEdges.end(), compareEdge);


        std::vector<bool> removeShortcut (roundEdges.size(), false);


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


        // ASSIGN CURRENT RANK TO NODES

        for(uint64_t nodeInContractedSet : contractedSet){
            workArray.rank.at(nodeInContractedSet) = currentRank;
        }
        currentRank++;
        numberOfNodesContracted += contractedSet.size();


        


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


        // BUILD NEW WORKARRAY G' = (V \ C, E u S \ A)


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

        workArrayNew.offsets.push_back(0);
        uint64_t currentOffset = 0;
        for(uint64_t nodeId = 0; nodeId < workArray.rank.size(); nodeId++){
            
            if(roundEdges.size() > 0){
                bool whileAlive = true;
                while(roundEdges.at(roundEdgeIndex).startNode == nodeId && whileAlive){
                    // only add shortcut if there is no identical other one
                    if(!removeShortcut.at(roundEdgeIndex)){
                        workArrayNew.edges.push_back(roundEdges.at(roundEdgeIndex).targetNode);
                        workArrayNew.distances.push_back(roundEdges.at(roundEdgeIndex).distance);
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
            for(uint64_t currEdge = workArray.offsets.at(nodeId); currEdge < workArray.offsets.at(nodeId +1); currEdge++){
                if(!adjacentEdges.at(currEdge)){
                    workArrayNew.edges.push_back(workArray.edges.at(currEdge));
                    workArrayNew.distances.push_back(workArray.distances.at(currEdge));
                    currentOffset++;
                }

            }

            workArrayNew.offsets.push_back(currentOffset);
        }

        std::cout << "workarray new: " << workArrayNew.offsets.at(workArrayNew.offsets.size() - 1)  << " " <<  workArray.offsets.at(workArray.offsets.size() - 1)<< std::endl;
        std::cout << "workarray new: " << workArrayNew.offsets.size()  << " " <<  workArray.offsets.size() << std::endl;
        std::cout << "prevented shortcuts: " << preventedShortcuts <<   " from " << roundEdges.size() << std::endl;

        workArray = workArrayNew;

        std::cout << "--------------------------------------------------------------" << std::endl;

        std::cout << "workarray new: " << workArrayNew.offsets.at(workArrayNew.offsets.size() - 1)  << " " <<  workArray.offsets.at(workArray.offsets.size() - 1)<< std::endl;
        std::cout << "workarray new: " << workArrayNew.offsets.size()  << " " <<  workArray.offsets.size() << std::endl;





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
            finalGraph.rank.at(currNode) == currentRank;
        }
    }

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


    std::cout << "finalGraph " << finalGraph.edges.size() << " from "  << array.edges.size() << std::endl;

    array = finalGraph;

    std::cout << "finalGraph " << array.edges.size() << std::endl;

    std::cout << "--- END CONTRACT HIERACHIES ---" << std::endl;

}