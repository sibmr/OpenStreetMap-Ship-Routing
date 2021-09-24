#include <ctime>
#include <numeric>

#include "shortestPathUtils.cpp"
#include "Dijkstra_simon.cpp"

struct Edge
{
    uint64_t edgeId;
    uint64_t v1;
    uint64_t v2;
    uint64_t edgeDistance;
    std::vector<uint64_t> shortcutPathEdges;
};

std::vector<Edge> copyEdgeVector(std::vector<Edge> &edgeVector){
    std::vector<Edge> copy;
    for (Edge &edge : edgeVector){
        copy.push_back(Edge{
            .edgeId = edge.edgeId,
            .v1 = edge.v1,
            .v2 = edge.v2,
            .edgeDistance = edge.edgeDistance,
            .shortcutPathEdges = std::vector<uint64_t>(edge.shortcutPathEdges)
        });
    }
    return copy;
}

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

void addEdgeIds(AdjacencyArray &adjArray){
    for(uint64_t edgeId = 0; edgeId<adjArray.edges.size(); ++edgeId){
        adjArray.edgeIds.push_back(edgeId);
    }
}

void nodeToEdgeIdPath(AdjacencyArray &adjArray, std::vector<uint64_t> &in_nodeIdPath, std::vector<uint64_t> &out_edgeIdPath){
    uint64_t prev_nodeId = in_nodeIdPath.at(0);
    for(uint64_t nodeIndex = 1; nodeIndex<in_nodeIdPath.size(); ++nodeIndex){
        uint64_t nodeId = in_nodeIdPath.at(nodeIndex);
        for(uint64_t currEdgeIndex = adjArray.offsets.at(prev_nodeId); currEdgeIndex < adjArray.offsets.at(prev_nodeId+1); ++currEdgeIndex){
            uint64_t nextNode = adjArray.edges.at(currEdgeIndex);
            if(nextNode == nodeId){
                out_edgeIdPath.push_back(adjArray.edgeIds.at(currEdgeIndex));
            }
        }
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

void fillContractionSet(
    AdjacencyArray &adjArray, 
    std::vector<uint64_t> &allContractedIds, 
    std::vector<bool> &isContracted, 
    std::vector<uint64_t> &out_newContractions
    ){

    uint64_t numDraws = 200000;
    uint64_t maxNodeId = isContracted.size();
    
    uint64_t nodeIdx = 0;
    for(int i = 0; i<numDraws; ++i){
        
        // random position
        uint64_t draw = std::rand() % maxNodeId;
        
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

int main(int argc, char** argv){
    std::string path;
    if(getFilePathInput(argc, argv, path)==0){
        return 0;
    }

    std::srand(std::time(nullptr));
    
    AdjacencyArray adjArray(path);
    addEdgeIds(adjArray);

    std::vector<bool> isContracted(adjArray.width*adjArray.height, false);
    std::vector<uint64_t> contractedNodeIds;
    std::vector<Edge> allEdges;
    extractAllEdges(adjArray, allEdges);
    uint64_t currentRank = 0;

    uint64_t numWaterNodes = 0;
    for(bool onLand : adjArray.nodes){
        numWaterNodes += !onLand;
    }

    // copy original adjacency array
    AdjacencyArray workArray = AdjacencyArray(adjArray);

    // while there are still uncontracted nodes
    //for(uint64_t round = 0; round<5; ++round){
    uint64_t round = 0;
    std::cout << workArray.width*workArray.height - std::accumulate(workArray.nodes.begin(), workArray.nodes.end(), 0) << "\n";
    while((workArray.width*workArray.height - std::accumulate(workArray.nodes.begin(), workArray.nodes.end(), 0)) > 20000){
        std::cout << "round " << round++ << "\n";

        Dijkstra::Dijkstra dijkstra(workArray);
        
        std::vector<Edge> shortcutEdges;
        std::vector<uint64_t> removedEdgeIndices;
        std::vector<bool> isEdgeRemoved(workArray.edges.size(), false);

        // select set of independent nodes for contraction
        std::vector<uint64_t> newContractions;
        fillContractionSet(workArray, contractedNodeIds, isContracted, newContractions);

        // create shortcuts to replace edges of nodes to contract/remove
        for(uint64_t contractedNodeId : newContractions){
            
            // set the rank of contractedNodes to current rank
            workArray.rank.at(contractedNodeId) = currentRank;

            // find all neighbor ids and all edges (to and from contracted node)
            std::vector<uint64_t> adjacentIds; 
            std::vector<uint64_t> edgeIds;
            std::vector<uint64_t> edgeIndices;
            // process edges away from contracted node
            for(uint64_t currEdgeIndex = workArray.offsets.at(contractedNodeId); currEdgeIndex < workArray.offsets.at(contractedNodeId+1); currEdgeIndex++){
                uint64_t adjacentId = workArray.edges.at(currEdgeIndex);
                adjacentIds.push_back(adjacentId);
                edgeIds.push_back(workArray.edgeIds.at(currEdgeIndex));
                edgeIndices.push_back(currEdgeIndex);
                isEdgeRemoved.at(currEdgeIndex) = true;
            }
            // process edges to contracted node
            for(uint64_t neighborId : adjacentIds){
                for(uint64_t currEdgeIndex = workArray.offsets.at(neighborId); currEdgeIndex < workArray.offsets.at(neighborId+1); currEdgeIndex++){
                    uint64_t neighborOfNeighborId = workArray.edges.at(currEdgeIndex);
                    if(contractedNodeId == neighborOfNeighborId){
                        edgeIds.push_back(workArray.edgeIds.at(currEdgeIndex));
                        edgeIndices.push_back(currEdgeIndex);
                        isEdgeRemoved.at(currEdgeIndex) = true;
                    }
                }
            }

            //calculate distances
            std::vector<std::vector<uint64_t>> distanceUWMatrix = createVectorMatrix(adjacentIds.size(), adjacentIds.size());
            for(uint64_t uIndex = 0; uIndex<adjacentIds.size(); ++uIndex){
                uint64_t uId = adjacentIds.at(uIndex);
                dijkstra.reset();
                dijkstra.disableNode(contractedNodeId); // disable contracted node
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
            std::vector<uint64_t> distanceUV(adjacentIds.size()); // for calculating uvw
            std::vector<uint64_t> shortcutPathPart;
            std::vector<std::vector<uint64_t>> shortcutEdgePathPart(adjacentIds.size());
            if(adjacentIds.size() > 0){
                dijkstra.reset();
                dijkstra.calculateDist(contractedNodeId, adjacentIds.at(0)); // set correct start point
                for(uint64_t uIndex = 0; uIndex<adjacentIds.size(); ++uIndex){
                    uint64_t uId = adjacentIds.at(uIndex);
                    uint64_t distance = dijkstra.calculateDist(uId); // only query with end point (one-to-many)
                    distanceUV.at(uIndex) = distance;
                    shortcutPathPart.clear();
                    dijkstra.getPath(shortcutPathPart); // get shortest path from v to u
                    nodeToEdgeIdPath(workArray, shortcutPathPart, shortcutEdgePathPart.at(uIndex));
                }
            }
            else{
                std::cout << "empty node: " << contractedNodeId << "\n";
            }
            
            

            // insert shortcuts in both directions
            for(uint64_t i = 0; i<adjacentIds.size(); ++i){
                for(uint64_t j = 0; j<adjacentIds.size(); ++j){
                    if(i == j){ continue; }
                    
                    uint64_t startNode = adjacentIds.at(i);
                    uint64_t endNode = adjacentIds.at(j);
                    
                    uint64_t distanceUW = distanceUWMatrix.at(i).at(j);
                    uint64_t distanceUVW = distanceUV.at(i) + distanceUV.at(j); // i is u, j is w

                    // if replaced edges are potentially part of a shortest path, then add shortcut
                    if(distanceUW >= distanceUVW){
                        // add shortcut to list
                        // use distances of edges (contractedNodeId, otherNodeId) - instead of (otherNodeId, contractedNodeId)
                        uint64_t edgeDistance = distanceUVW;
                        uint64_t shortcutId = allEdges.size();
                        std::vector<uint64_t> edgePathUVW;
                        std::vector<uint64_t> edgePathVU = shortcutEdgePathPart.at(i);
                        std::vector<uint64_t> edgePathVW = shortcutEdgePathPart.at(j);
                        for(uint64_t pathIndex = edgePathVU.size()-1; pathIndex >= 0 && pathIndex < edgePathVU.size(); --pathIndex){
                            edgePathUVW.push_back(edgePathVU.at(pathIndex));
                        }
                        for(uint64_t pathIndex = 0; pathIndex < edgePathVW.size(); --pathIndex){
                            edgePathUVW.push_back(edgePathVW.at(pathIndex));
                        }
                        Edge newEdge{
                            .edgeId=shortcutId,
                            .v1=startNode, 
                            .v2=endNode, 
                            .edgeDistance=edgeDistance, 
                            .shortcutPathEdges=edgePathUVW
                            };
                        allEdges.push_back(newEdge);
                        shortcutEdges.push_back(newEdge);
                    }
                }
            }

            //mark adjacent edges of contractedNode for removal
            removedEdgeIndices.insert(removedEdgeIndices.end(), edgeIndices.begin(), edgeIndices.end());
        }
        
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
                    
                    //std::cout << "target: " << nodeId << "  currentId: " << shortcutEdgeIt->v1 << "\n";
                    
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
        std::cout << "number of offsets: " << workArray.offsets.size() << "\n";
        std::cout << "number of offsets: " << workArray.height*workArray.width << "\n";
        std::cout << "number of shortcuts: " << shortcutEdges.size() << "\n";
        std::cout << "number of removed edges: " << removedEdgeIndices.size() << "\n";
        std::cout << "remaining nodes: " << workArray.width*workArray.height - std::accumulate(workArray.nodes.begin(), workArray.nodes.end(), 0) << "\n";
    }
    // finally: G'' = (V, E u E')
}