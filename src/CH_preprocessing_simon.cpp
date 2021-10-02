#include <ctime>
#include <numeric>
# include <omp.h>

#include "shortestPathUtils.cpp"
#include "Dijkstra_simon.cpp"


std::vector<std::vector<uint64_t>> createVectorMatrix(uint64_t size_i, uint64_t size_j){
    std::vector<std::vector<uint64_t>> matrix(size_i);
    for(uint64_t i = 0; i<size_i; ++i){
        matrix.at(i).resize(size_j);
    }
    return matrix;
}

int getFilePathInput(int argc, char** argv, std::string &out_path){
    if(argc > 1){
        out_path = std::string(argv[1]);
    }else{
        out_path = "data/planet.graph_2";
        std::cout << "no input file given assume " <<  out_path << std::endl;
    }

    {
        std::ifstream f(out_path);
        if(!f.good()){
            std::cout << "file: " << out_path << " not found\n";
            return 0;
        }
    }
    return 1;
}

void nodeToEdgeIdPathForward(AdjacencyArray &adjArray, std::vector<uint64_t> &in_nodeIdPath, std::vector<uint64_t> &out_edgeIdPath){
    uint64_t prev_nodeId = in_nodeIdPath.at(0);
    for(uint64_t nodeIndex = 1; nodeIndex<in_nodeIdPath.size(); ++nodeIndex){
        uint64_t nodeId = in_nodeIdPath.at(nodeIndex);
        for(uint64_t currEdgeIndex = adjArray.offsets.at(prev_nodeId); currEdgeIndex < adjArray.offsets.at(prev_nodeId+1); ++currEdgeIndex){
            uint64_t nextNode = adjArray.edges.at(currEdgeIndex);
            if(nextNode == nodeId){
                out_edgeIdPath.push_back(adjArray.edgeIds.at(currEdgeIndex));
            }
        }
        prev_nodeId = nodeId;
    }
}

void nodeToEdgeIdPathBackward(AdjacencyArray &adjArray, std::vector<uint64_t> &in_nodeIdPath, std::vector<uint64_t> &out_edgeIdPath){
    uint64_t prev_nodeId = in_nodeIdPath.at(in_nodeIdPath.size()-1);
    for(uint64_t nodeIndex = in_nodeIdPath.size()-2; nodeIndex>=0 && nodeIndex < in_nodeIdPath.size(); --nodeIndex){
        uint64_t nodeId = in_nodeIdPath.at(nodeIndex);
        for(uint64_t currEdgeIndex = adjArray.offsets.at(prev_nodeId); currEdgeIndex < adjArray.offsets.at(prev_nodeId+1); ++currEdgeIndex){
            uint64_t nextNode = adjArray.edges.at(currEdgeIndex);
            if(nextNode == nodeId){
                out_edgeIdPath.push_back(adjArray.edgeIds.at(currEdgeIndex));
            }
        }
        prev_nodeId = nodeId;
    }
}

void extractAllEdges(AdjacencyArray &adjArray, std::vector<Edge> &out_allEdges){
    uint64_t edgeIdCounter = 0;
    uint64_t noEdgesCounter = 0;
    for(uint64_t nodeId = 0; nodeId < adjArray.width*adjArray.height; ++nodeId){
        if(!adjArray.nodes.at(nodeId)){
            for(uint64_t currEdgeIndex = adjArray.offsets.at(nodeId); currEdgeIndex < adjArray.offsets.at(nodeId+1); ++currEdgeIndex){
                uint64_t adjacentId = adjArray.edges.at(currEdgeIndex);
                out_allEdges.push_back(
                    Edge{
                        .edgeId = adjArray.edgeIds.at(currEdgeIndex),
                        .v1 = nodeId, 
                        .v2 = adjacentId, 
                        .edgeDistance = adjArray.distances.at(currEdgeIndex),
                        .shortcutPathEdges = {}
                        }
                );
            }

            if(adjArray.offsets.at(nodeId) == adjArray.offsets.at(nodeId+1)){
                noEdgesCounter++;
            }
        }
    }
    std::cout << "Number of active nodes without edges: " << noEdgesCounter << "\n";
}

bool isNodeAdjacentToSet(AdjacencyArray &adjArray, uint64_t idToCheck, std::vector<uint64_t> &idSet){
    for(uint64_t currentIdx : idSet){
        for(uint64_t currEdgeIndex = adjArray.offsets.at(idToCheck); currEdgeIndex < adjArray.offsets.at(idToCheck+1); currEdgeIndex++){
            uint64_t adjacentId = adjArray.edges.at(currEdgeIndex);
            if(adjacentId == currentIdx){
                return true;
            }
        }
    }
    return false;
}

bool isNodeAdjacentToSet(AdjacencyArray &adjArray, uint64_t idToCheck, std::vector<bool> &isIdInSet){
    for(uint64_t currEdgeIndex = adjArray.offsets.at(idToCheck); currEdgeIndex < adjArray.offsets.at(idToCheck+1); ++currEdgeIndex){
        uint64_t adjacentId = adjArray.edges.at(currEdgeIndex);
        if(isIdInSet.at(adjacentId)){
            return true;
        }
    }
    return false;
}

void randomFillContractionSet(
    AdjacencyArray &adjArray, 
    std::vector<uint64_t> &allContractedIds, 
    std::vector<bool> &isContracted, 
    std::vector<uint64_t> &out_newContractions
    ){

    uint64_t numDraws = 20000;
    uint64_t nodeIdLimit = isContracted.size();
    
    uint64_t nodeIdx = 0;
    for(int i = 0; i<numDraws; ++i){
        
        // random position
        uint64_t draw = std::rand() % nodeIdLimit;
        
        // check if node is in water and not already in contraction set
        if(!isContracted.at(draw) && !adjArray.nodes.at(draw)){
            
            // check if node is not adjacent to nodes in contraction set
            bool isAdjacent = isNodeAdjacentToSet(adjArray, draw, out_newContractions);

            // if node is not adjacent, add to contraction set
            if(!isAdjacent){
                out_newContractions.push_back(draw);
                allContractedIds.push_back(draw);
                isContracted.at(draw) = true;
            }
        
        }
    }
}

void numEdgeFillContractionSet(
    AdjacencyArray &adjArray,
    double fraction, 
    std::vector<uint64_t> &in_out_allContractedIds, 
    std::vector<bool> &in_out_isContracted, 
    std::vector<uint64_t> &out_newContractions
    ){

    uint64_t numDraws = 20000;
    uint64_t nodeIdLimit = in_out_isContracted.size();
    
    // pick first node
    uint64_t initialId = std::rand() % nodeIdLimit;
    std::vector<bool> isInIndependentSet(adjArray.width*adjArray.height, false);
    std::vector<std::pair<uint64_t,uint64_t>> independentSetWithCost;
    do{
        initialId = std::rand() % nodeIdLimit;
    }while(in_out_isContracted.at(initialId) || adjArray.nodes.at(initialId));
    isInIndependentSet.at(initialId) = true;
    independentSetWithCost.push_back(
        std::pair<uint64_t, uint64_t>(initialId, adjArray.offsets.at(initialId+1) - adjArray.offsets.at(initialId))
    );
    
    for(uint64_t nodeId=0; nodeId<adjArray.width*adjArray.height; ++nodeId){
        
        // check if node is in water and not already in contraction set
        if(!in_out_isContracted.at(nodeId) && !adjArray.nodes.at(nodeId)){
            
            // check if node is not adjacent to nodes in contraction set
            bool isAdjacent = isNodeAdjacentToSet(adjArray, nodeId, isInIndependentSet);

            // calculate num edges
            uint64_t numEdges = adjArray.offsets.at(nodeId+1) - adjArray.offsets.at(nodeId);

            // if node is not adjacent, add to independent set
            if(!isAdjacent){
                independentSetWithCost.push_back(std::pair<uint64_t, uint64_t>(nodeId, numEdges));
                isInIndependentSet.at(nodeId) = true;
            }
        }
    }

    std::sort(  independentSetWithCost.begin(), independentSetWithCost.end(),
                [](std::pair<uint64_t,uint64_t> &p1, std::pair<uint64_t,uint64_t> &p2){ return p1.second < p2.second; });
    
    for(uint64_t nodeIndex=0; nodeIndex<((uint64_t)(independentSetWithCost.size()*fraction)); ++nodeIndex){
        uint64_t nodeId = independentSetWithCost.at(nodeIndex).first;
        uint64_t cost = independentSetWithCost.at(nodeIndex).second;
        out_newContractions.push_back(nodeId);
        in_out_allContractedIds.push_back(nodeId);
        in_out_isContracted.at(nodeId) = true;
    }
    std::cout << "num nodes contracted: " << out_newContractions.size() << "\n";
}


void contractNode(  uint64_t contractedNodeId, AdjacencyArray &workArray, uint64_t currentRank, Dijkstra::Dijkstra &dijkstra,
                    std::vector<bool> &out_isEdgeRemoved, std::vector<uint64_t> &out_removedEdgeIndices,
                    std::vector<Edge> &out_allEdges, std::vector<Edge> &out_shortcutEdges){

    // set the rank of contractedNodes to current rank
    workArray.rank.at(contractedNodeId) = currentRank;

    // find all neighbor ids and all edges (to and from contracted node)
    std::vector<uint64_t> adjacentIds;
    std::vector<uint64_t> forwardEdgeIds;
    std::vector<uint64_t> forwardEdgeIndices;
    std::vector<uint64_t> edgeIds;
    std::vector<uint64_t> edgeIndices;
    // process edges away from contracted node
    for(uint64_t currEdgeIndex = workArray.offsets.at(contractedNodeId); currEdgeIndex < workArray.offsets.at(contractedNodeId+1); ++currEdgeIndex){
        uint64_t adjacentId = workArray.edges.at(currEdgeIndex);
        bool isIdNotDuplicated = std::find(adjacentIds.begin(), adjacentIds.end(), adjacentId) == adjacentIds.end();
        if(isIdNotDuplicated){
            adjacentIds.push_back(adjacentId);
            forwardEdgeIds.push_back(workArray.edgeIds.at(currEdgeIndex));
            forwardEdgeIndices.push_back(currEdgeIndex);
        }
        edgeIds.push_back(workArray.edgeIds.at(currEdgeIndex));
        edgeIndices.push_back(currEdgeIndex);
        out_isEdgeRemoved.at(currEdgeIndex) = true;
    }
    
    std::vector<uint64_t> reverseEdgeIds(adjacentIds.size(), UINT64_MAX);
    // process edges to contracted node
    for(uint64_t neighborIndex = 0; neighborIndex<adjacentIds.size(); ++neighborIndex){
        uint64_t neighborId = adjacentIds.at(neighborIndex);
        uint64_t minEdgeId = UINT64_MAX;
        for(uint64_t currEdgeIndex = workArray.offsets.at(neighborId); currEdgeIndex < workArray.offsets.at(neighborId+1); ++currEdgeIndex){
            uint64_t neighborOfNeighborId = workArray.edges.at(currEdgeIndex);
            if(contractedNodeId == neighborOfNeighborId){ // seems like there is not always a backward edge
                edgeIds.push_back(workArray.edgeIds.at(currEdgeIndex));
                edgeIndices.push_back(currEdgeIndex);
                out_isEdgeRemoved.at(currEdgeIndex) = true;
                if(workArray.edgeIds.at(currEdgeIndex) < minEdgeId){
                    minEdgeId = workArray.edgeIds.at(currEdgeIndex);
                    reverseEdgeIds.at(neighborIndex) = workArray.edgeIds.at(currEdgeIndex); // there are still uninitialized entries
                }
            }
        }
    }

    if(edgeIds.size() != 2*adjacentIds.size()){ std::cout << "more/less edges" << edgeIds.size() << " " << adjacentIds.size() << "\n"; }

    //calculate distances
    std::vector<std::vector<uint64_t>> distanceUWMatrix = createVectorMatrix(adjacentIds.size(), adjacentIds.size());
    for(uint64_t uIndex = 0; uIndex<adjacentIds.size(); ++uIndex){
        uint64_t uId = adjacentIds.at(uIndex);
        dijkstra.reset(); // do not disable any nodes
        if(uIndex+1<adjacentIds.size()){
            dijkstra.calculateDist(uId, adjacentIds.at(uIndex+1)); // set correct start point
        } 
        for(uint64_t wIndex = uIndex+1; wIndex<adjacentIds.size(); ++wIndex){
            uint64_t wId = adjacentIds.at(wIndex);
            uint64_t distance = dijkstra.calculateDist(wId); // only query with end point (one-to-many)
            distanceUWMatrix.at(uIndex).at(wIndex) = distance;
            distanceUWMatrix.at(wIndex).at(uIndex) = distance;
        }
    }
    

    // insert shortcuts in both directions
    for(uint64_t i = 0; i<adjacentIds.size(); ++i){
        for(uint64_t j = 0; j<adjacentIds.size(); ++j){
            if(i == j){ continue; }
            
            uint64_t startNode = adjacentIds.at(i); // i is u
            uint64_t endNode = adjacentIds.at(j);   // j is w
            uint64_t startNodeEdgeIndex = forwardEdgeIndices.at(i);
            uint64_t endNodeEdgeIndex = forwardEdgeIndices.at(j);

            uint64_t distanceUW = distanceUWMatrix.at(i).at(j);
            uint64_t distanceUVW = workArray.distances.at(startNodeEdgeIndex) + workArray.distances.at(endNodeEdgeIndex);

            // if deleted edges are potentially part of a shortest path, then add shortcut edge
            if(distanceUW >= distanceUVW){
                // add shortcut to list
                // use distances of edges (contractedNodeId, otherNodeId) - instead of (otherNodeId, contractedNodeId)
                uint64_t edgeDistance = distanceUVW;
                uint64_t shortcutId = out_allEdges.size();
                std::vector<uint64_t> edgePathUVW;
                std::vector<uint64_t> edgePathUV = {reverseEdgeIds.at(i)};
                std::vector<uint64_t> edgePathVW = {forwardEdgeIds.at(j)};
                edgePathUVW.insert(edgePathUVW.end(), edgePathUV.begin(), edgePathUV.end());
                edgePathUVW.insert(edgePathUVW.end(), edgePathVW.begin(), edgePathVW.end());

                Edge newEdge{
                    .edgeId=shortcutId,
                    .v1=startNode, 
                    .v2=endNode, 
                    .edgeDistance=edgeDistance, 
                    .shortcutPathEdges=edgePathUVW
                    };
                out_allEdges.push_back(newEdge);
                out_shortcutEdges.push_back(newEdge);
            }
        }
    }

    //mark adjacent edges of contractedNode for removal
    out_removedEdgeIndices.insert(out_removedEdgeIndices.end(), edgeIndices.begin(), edgeIndices.end());
}


int main(int argc, char** argv){
    std::string path;
    if(getFilePathInput(argc, argv, path)==0){
        return 0;
    }

    std::srand(std::time(nullptr));
    
    AdjacencyArray adjArray(path);

    AdjacencyArray &finalArray = *(new AdjacencyArray(adjArray));

    std::vector<bool> isContracted(adjArray.width*adjArray.height, false);
    std::vector<uint64_t> contractedNodeIds;
    std::vector<Edge> allEdges;
    extractAllEdges(adjArray, allEdges);
    uint64_t currentRank = 1;

    uint64_t numWaterNodes = 0;
    for(bool onLand : adjArray.nodes){
        numWaterNodes += !onLand;
    }

    // copy original adjacency array
    AdjacencyArray workArray = AdjacencyArray(adjArray);

    // while there are still uncontracted nodes
    //for(uint64_t round = 0; round<5; ++round){
    uint64_t round = 1;
    std::cout << workArray.width*workArray.height - std::accumulate(workArray.nodes.begin(), workArray.nodes.end(), 0) << "\n";
    while((workArray.width*workArray.height - std::accumulate(workArray.nodes.begin(), workArray.nodes.end(), 0)) > 7000){
        std::cout << "round " << round++ << "\n";

        Dijkstra::Dijkstra dijkstra(workArray);
        
        std::vector<Edge> shortcutEdges;
        std::vector<uint64_t> removedEdgeIndices;
        std::vector<bool> isEdgeRemoved(workArray.edges.size(), false);

        // select set of independent nodes for contraction
        std::vector<uint64_t> newContractions;
        numEdgeFillContractionSet(workArray, 0.5, contractedNodeIds, isContracted, newContractions);

        uint64_t roundProgress = 0;

        // create shortcuts to replace edges of nodes to contract/remove
        for(uint64_t contractedNodeId : newContractions){
            
            if(roundProgress++%1024==0){
                std::cout << "round " << round-1 << " progress " << ((double)roundProgress)/newContractions.size()  << "\t\r" << std::flush;
            }
            
            contractNode(   contractedNodeId, workArray, currentRank,
                            dijkstra, isEdgeRemoved, removedEdgeIndices, allEdges, shortcutEdges);
        }

        // parallel version
        // uint64_t edgeCount = allEdges.size();
        // #pragma omp parallel num_threads(8)
        // {
        //     std::vector<Edge> shortcutEdges_thread;
        //     std::vector<uint64_t> removedEdgeIndices_thread;
        //     Dijkstra::Dijkstra dijkstra_thread(workArray);
            
        //     #pragma omp for
        //     for(uint64_t contractedNodeId : newContractions){
                
        //         if(roundProgress++%1024==0){
        //             std::cout << "round " << round-1 << " progress " << ((double)roundProgress)/newContractions.size()  << "\t\r" << std::flush;
        //         }
                
        //         contractNodeParallel(   contractedNodeId, workArray, currentRank,
        //                                 dijkstra_thread, removedEdgeIndices_thread, shortcutEdges_thread, edgeCount);
        //     }
            
        //     #pragma omp critical
        //     {
        //         shortcutEdges.insert(shortcutEdges.end(), shortcutEdges_thread.begin(), shortcutEdges_thread.end());
        //         removedEdgeIndices.insert(removedEdgeIndices.end(), removedEdgeIndices_thread.begin(), removedEdgeIndices_thread.end());
        //     }
        // }
        // allEdges.insert(allEdges.end(), shortcutEdges.begin(), shortcutEdges.end());
        // for(uint64_t removedEdgeIndex : removedEdgeIndices){ isEdgeRemoved.at(removedEdgeIndex) = true; };

        std::cout << "rebuild phase\n";

        // build intermediate graph
        AdjacencyArray tempArray = AdjacencyArray(workArray);
        // adjust temp array
        tempArray.edges.clear();
        tempArray.edgeIds.clear();
        tempArray.distances.clear();
        tempArray.offsets.clear();
        tempArray.offsets.push_back(0);
        
        std::vector<Edge> shortcutSortedV1 = copyEdgeVector(shortcutEdges);
        std::sort(  
            shortcutSortedV1.begin(), shortcutSortedV1.end(), 
            [](Edge& a, Edge& b) {return a.v1 < b.v1; });

        uint64_t currentOffset = 0;
        for(uint64_t nodeId = 0; nodeId<tempArray.width*tempArray.height; ++nodeId){
            bool activeNode = !tempArray.nodes.at(nodeId);
            bool wasContracted = isContracted.at(nodeId);
            
            if(activeNode && wasContracted){
                // contracted node is no longer active
                tempArray.nodes.at(nodeId) = true;
            }
            else if(activeNode && !wasContracted){
                // insert old edges (except for removed ones)
                for(uint64_t currEdgeIndex = workArray.offsets.at(nodeId); currEdgeIndex < workArray.offsets.at(nodeId+1); currEdgeIndex++){
                    // check if edge was removed
                    if(!isEdgeRemoved.at(currEdgeIndex)){
                        tempArray.edges.push_back(workArray.edges.at(currEdgeIndex));
                        tempArray.edgeIds.push_back(workArray.edgeIds.at(currEdgeIndex));
                        tempArray.distances.push_back(workArray.distances.at(currEdgeIndex));
                        
                        currentOffset++;
                    }
                }

                // insert shortcut edges (if node is part of any)
                auto resV1 = std::lower_bound(
                    shortcutSortedV1.begin(), shortcutSortedV1.end(), 
                    nodeId, [](Edge& e1, uint64_t targetId){return e1.v1 < targetId;});

                for(auto shortcutEdgeIt = resV1; (*shortcutEdgeIt).v1 <= nodeId && shortcutEdgeIt<shortcutSortedV1.end(); ++shortcutEdgeIt){
                    
                    if(nodeId != shortcutEdgeIt->v1){ std::cout << "nodeID!=shortcutEdge.v1"  << "\n"; }
                    
                    tempArray.edges.push_back(shortcutEdgeIt->v2);
                    tempArray.edgeIds.push_back(shortcutEdgeIt->edgeId);
                    tempArray.distances.push_back(shortcutEdgeIt->edgeDistance);
                    currentOffset++;
                }


            }

            // finalize AdjacencyArray node entry by adding right offset border
            tempArray.offsets.push_back(currentOffset);
        }
        
        workArray = std::move(tempArray);

        // increment rank
        currentRank++;

        // print round stats
        std::cout << "number of edges: " << workArray.edges.size() << "\n";
        std::cout << "total number of edges: " << allEdges.size() << "\n";
        std::cout << "number of offsets: " << workArray.offsets.size() << "\n";
        std::cout << "number of shortcuts: " << shortcutEdges.size() << "\n";
        std::cout << "number of removed edges: " << removedEdgeIndices.size() << "\n";
        std::cout << "remaining nodes: " << workArray.width*workArray.height - std::accumulate(workArray.nodes.begin(), workArray.nodes.end(), 0) << "\n";
        std::cout << "current rank: " << currentRank << "\n";
    
        



        // finally: G'' = (V, E u E')
        
        // prepare array
        
        finalArray.edges.clear();
        finalArray.edgeIds.clear();
        finalArray.distances.clear();
        finalArray.offsets.clear();
        finalArray.offsets.push_back(0);
        finalArray.allEdgeInfo = allEdges;
        finalArray.rank = workArray.rank;
        
        // sort edges
        std::vector<Edge> allEdgesSortedV1 = std::vector<Edge>(allEdges); // shallow copy only
        std::sort(  
            allEdgesSortedV1.begin(), allEdgesSortedV1.end(), 
            [](Edge& a, Edge& b) {return a.v1 < b.v1; });
        
        // add normal and shortcut edges to final adjacency array
        uint64_t currentOffsetIndex = 0;
        for(uint64_t nodeId = 0; nodeId<adjArray.width*adjArray.height; ++nodeId){
            auto resV1 = std::lower_bound(
                allEdgesSortedV1.begin(), allEdgesSortedV1.end(), 
                nodeId, [](Edge& e1, uint64_t targetId){return e1.v1 < targetId;});

            for(auto edgeIt = resV1; (*edgeIt).v1 <= nodeId && edgeIt<allEdgesSortedV1.end(); ++edgeIt){
                finalArray.edges.push_back(edgeIt->v2);
                finalArray.edgeIds.push_back(edgeIt->edgeId);
                finalArray.distances.push_back(edgeIt->edgeDistance);
                currentOffsetIndex++;
            }
            finalArray.offsets.push_back(currentOffsetIndex);
        }
        finalArray.writeToDisk("data/CHAdjArray_" + std::to_string(currentRank-1) + ".graph_2");

        // workArray.allEdgeInfo = allEdges;
        // workArray.writeToDisk("data/CHAdjArray.graph_2");

        

    
    
    }
    
    delete &finalArray;
}


// void contractNodeParallel(  uint64_t contractedNodeId, AdjacencyArray &workArray, uint64_t currentRank, Dijkstra::Dijkstra &dijkstra,
//                             std::vector<uint64_t> &out_removedEdgeIndices,
//                             std::vector<Edge> &out_shortcutEdges,
//                             uint64_t &in_out_edgeCount){

//     // set the rank of contractedNodes to current rank
//     workArray.rank.at(contractedNodeId) = currentRank;

//     // find all neighbor ids and all edges (to and from contracted node)
//     std::vector<uint64_t> adjacentIds; 
//     std::vector<uint64_t> edgeIds;
//     std::vector<uint64_t> edgeIndices;
    
//     // process edges away from contracted node
//     for(uint64_t currEdgeIndex = workArray.offsets.at(contractedNodeId); currEdgeIndex < workArray.offsets.at(contractedNodeId+1); ++currEdgeIndex){
//         uint64_t adjacentId = workArray.edges.at(currEdgeIndex);
//         bool isIdNotDuplicated = std::find(adjacentIds.begin(), adjacentIds.end(), adjacentId) == adjacentIds.end();
//         if(isIdNotDuplicated){
//             adjacentIds.push_back(adjacentId);
//             edgeIds.push_back(workArray.edgeIds.at(currEdgeIndex));
//             edgeIndices.push_back(currEdgeIndex);
//         }
//     }

//     std::vector<uint64_t> reverseEdgeIds(adjacentIds.size(), UINT64_MAX);
//     // process edges to contracted node
//     for(uint64_t neighborIndex = 0; neighborIndex<adjacentIds.size(); ++neighborIndex){
//         uint64_t neighborId = adjacentIds.at(neighborIndex);
//         uint64_t minEdgeId = UINT64_MAX;
//         for(uint64_t currEdgeIndex = workArray.offsets.at(neighborId); currEdgeIndex < workArray.offsets.at(neighborId+1); ++currEdgeIndex){
//             uint64_t neighborOfNeighborId = workArray.edges.at(currEdgeIndex);
//             if(contractedNodeId == neighborOfNeighborId){
//                 edgeIds.push_back(workArray.edgeIds.at(currEdgeIndex));
//                 edgeIndices.push_back(currEdgeIndex);
//                 if(workArray.edgeIds.at(currEdgeIndex) < minEdgeId){
//                     minEdgeId = workArray.edgeIds.at(currEdgeIndex);
//                     reverseEdgeIds.at(neighborIndex) = currEdgeIndex; // there are still uninitialized entries
//                 }
//             }
//         }
//     }

//     if(edgeIds.size() != 2*adjacentIds.size()){ std::cout << "more/less edges" << edgeIds.size() << " " << adjacentIds.size() << "\n"; }

//     //calculate distances
//     std::vector<std::vector<uint64_t>> distanceUWMatrix = createVectorMatrix(adjacentIds.size(), adjacentIds.size());
//     for(uint64_t uIndex = 0; uIndex<adjacentIds.size(); ++uIndex){
//         uint64_t uId = adjacentIds.at(uIndex);
//         dijkstra.reset(); // do not disable any nodes
//         if(uIndex+1<adjacentIds.size()){
//             dijkstra.calculateDist(uId, adjacentIds.at(uIndex+1)); // set correct start point
//         } 
//         for(uint64_t wIndex = uIndex+1; wIndex<adjacentIds.size(); ++wIndex){
//             uint64_t wId = adjacentIds.at(wIndex);
//             uint64_t distance = dijkstra.calculateDist(wId); // only query with end point (one-to-many)
//             distanceUWMatrix.at(uIndex).at(wIndex) = distance;
//             distanceUWMatrix.at(wIndex).at(uIndex) = distance;
//         }
//     }
    

//     // insert shortcuts in both directions
//     for(uint64_t i = 0; i<adjacentIds.size(); ++i){
//         for(uint64_t j = 0; j<adjacentIds.size(); ++j){
//             if(i == j){ continue; }
            
//             uint64_t startNode = adjacentIds.at(i); // i is u
//             uint64_t endNode = adjacentIds.at(j);   // j is w
//             uint64_t startNodeEdgeIndex = edgeIndices.at(i);
//             uint64_t endNodeEdgeIndex = edgeIndices.at(j);

//             uint64_t distanceUW = distanceUWMatrix.at(i).at(j);
//             uint64_t distanceUVW = workArray.distances.at(startNodeEdgeIndex) + workArray.distances.at(endNodeEdgeIndex);

//             // if replaced edges are potentially part of a shortest path, then add shortcut
//             if(distanceUW >= distanceUVW){
//                 // add shortcut to list
//                 // use distances of edges (contractedNodeId, otherNodeId) - instead of (otherNodeId, contractedNodeId)
//                 uint64_t edgeDistance = distanceUVW;
//                 uint64_t shortcutId = 0;
//                 #pragma omp critical
//                 {
//                     shortcutId = in_out_edgeCount++;
//                 }
//                 std::vector<uint64_t> edgePathUVW;
//                 std::vector<uint64_t> edgePathUV = {reverseEdgeIds.at(i)};
//                 std::vector<uint64_t> edgePathVW = {edgeIds.at(       j)};
//                 edgePathUVW.insert(edgePathUVW.end(), edgePathUV.begin(), edgePathUV.end());
//                 edgePathUVW.insert(edgePathUVW.end(), edgePathVW.begin(), edgePathVW.end());

//                 Edge newEdge{
//                     .edgeId=shortcutId,
//                     .v1=startNode, 
//                     .v2=endNode, 
//                     .edgeDistance=edgeDistance, 
//                     .shortcutPathEdges=edgePathUVW
//                     };
//                 out_shortcutEdges.push_back(newEdge);
//             }
//         }
//     }

//     //mark adjacent edges of contractedNode for removal
//     out_removedEdgeIndices.insert(out_removedEdgeIndices.end(), edgeIndices.begin(), edgeIndices.end());
// }