#include <set>
#include "shortestPathUtils.cpp"
#include "MultiDijkstra.cpp"

void recursiveUnpack(Edge &edge, std::vector<uint64_t> &unpackedShortcut, AdjacencyArray &adjArray){
    std::cout << "num shortcut edges: " << edge.shortcutPathEdges.size() << ", v1: " <<  edge.v1 << ", v2: " <<  edge.v2 << "\n";
    if(edge.shortcutPathEdges.size() > 0){
        for(uint64_t edgeId : edge.shortcutPathEdges){
            recursiveUnpack(adjArray.allEdgeInfo.at(edgeId), unpackedShortcut, adjArray);
        }
    }else{
        unpackedShortcut.push_back(edge.edgeId);
        std::cout << edge.v1 << ", rank " << adjArray.rank.at(edge.v1) << " -> " << edge.v2 << ", rank " << adjArray.rank.at(edge.v2) << "\n"; 
    }
}

/**
 * @brief convert a path of node ids to a path of edges between the nodes
 * 
 * Unused in current version of program (moved from CH_preprocessing)
 * 
 * @param adjArray          AdjacencyArray that the path is occuring in
 * @param in_nodeIdPath     input vector of node ids
 * @param out_edgeIdPath    output vector of edge ids
 */
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


/**
 * @brief convert a path of node ids to a path of edges betwee the nodes in the opposite direction as the original node path
 * 
 * Unused in current version of program (moved from CH_preprocessing)
 * 
 * @param adjArray          AdjacencyArray that the path is occuring in
 * @param in_nodeIdPath     input vector of node ids
 * @param out_edgeIdPath    output vector of edge ids
 */
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

std::vector<bool> checkRedundantEdges(AdjacencyArray &adjArray){
    std::vector<bool> isEdgeRemoved(adjArray.allEdgeInfo.size(), false);
    uint64_t numRemoved = 0;
    
    if(adjArray.allEdgeInfo.size() < 1){
        std::cout << "reduced function since edge info not filled\n";
    }
    
    for(uint64_t nodeId = 0; nodeId < adjArray.width * adjArray.height; ++nodeId){
        std::map<uint64_t, std::vector<uint64_t>> nodesTo;
        for(uint64_t currEdgeIndex = adjArray.offsets.at(nodeId); currEdgeIndex < adjArray.offsets.at(nodeId+1); ++currEdgeIndex){
            uint64_t adjacentNode = adjArray.edges.at(currEdgeIndex);
            uint64_t edgeId = adjArray.edgeIds.at(currEdgeIndex);
            auto entry = nodesTo.find(adjacentNode);
            if(entry == nodesTo.end()){
                nodesTo.insert(std::pair<uint64_t, std::vector<uint64_t>>(adjacentNode, {edgeId}));
            }else{
                entry->second.push_back(edgeId);
            }
        }
        for(auto pair : nodesTo){
            uint64_t lowestDist = UINT64_MAX;
            uint64_t lowestEdgeId = UINT64_MAX;
            for(auto edgeId : pair.second){
                uint64_t edgeDist = UINT64_MAX;
                if(adjArray.allEdgeInfo.size() > 0){
                    edgeDist = adjArray.allEdgeInfo.at(edgeId).edgeDistance;
                }else{
                    edgeDist = 0;
                }
                
                
                if(edgeDist < lowestDist){
                    lowestDist = edgeDist;
                    lowestEdgeId = edgeId;
                }
            }
            for(auto edgeId : pair.second){
                if(edgeId != lowestEdgeId){
                    isEdgeRemoved.at(edgeId) = true;
                    numRemoved++;
                }
            }
        }
    }
    std::cout << numRemoved << " redundant edges removed\n";
    return isEdgeRemoved;
}

void removeRedundantEdges(AdjacencyArray &adjArray){
    std::vector<bool> isEdgeRedundant = checkRedundantEdges(adjArray);
    AdjacencyArray &original = *(new AdjacencyArray(adjArray));
    adjArray.edges.clear();
    adjArray.edgeIds.clear();
    adjArray.offsets.clear();
    adjArray.offsets.push_back(0);
    adjArray.distances.clear();

    uint64_t currentOffset = 0;
    for(uint64_t nodeId = 0; nodeId<adjArray.width*adjArray.height; ++nodeId){
        // insert old edges (except for removed ones)
        for(uint64_t currEdgeIndex = original.offsets.at(nodeId); currEdgeIndex < original.offsets.at(nodeId+1); currEdgeIndex++){
            // check if edge was removed
            if(!isEdgeRedundant.at(original.edgeIds.at(currEdgeIndex))){
                adjArray.edges.push_back(original.edges.at(currEdgeIndex));
                adjArray.edgeIds.push_back(original.edgeIds.at(currEdgeIndex));
                adjArray.distances.push_back(original.distances.at(currEdgeIndex));
                
                currentOffset++;
            }
        }
        adjArray.offsets.push_back(currentOffset);
    }
    delete &original;
}

void checkUndirectedness(AdjacencyArray &adjArray){
    for(uint64_t nodeId = 0; nodeId < adjArray.width*adjArray.height; ++nodeId){
        //if(!adjArray.nodes.at(nodeId)){
            for(uint64_t edgeIndex = adjArray.offsets.at(nodeId); edgeIndex < adjArray.offsets.at(nodeId+1); ++edgeIndex){
                uint64_t adjacentNodeId = adjArray.edges.at(edgeIndex);
                bool foundBackedge = false;
                for(uint64_t backedgeIndex = adjArray.offsets.at(adjacentNodeId); backedgeIndex < adjArray.offsets.at(adjacentNodeId+1); ++backedgeIndex){
                    uint64_t adjacentToNeighborId = adjArray.edges.at(backedgeIndex);
                    if(adjacentToNeighborId == nodeId){
                        foundBackedge = true;
                    }
                }
                if(!foundBackedge){
                    std::cout << "did not find backedge\n";
                }
            }
        //}
    }
}

void checkEdgeIdSorted(AdjacencyArray &adjArray){
    for(uint64_t edgeId = 0; edgeId < adjArray.allEdgeInfo.size(); ++edgeId){
        if(adjArray.allEdgeInfo.at(edgeId).edgeId != edgeId){
            std::cout << "differs: " << edgeId << " " << adjArray.allEdgeInfo.at(edgeId).edgeId << "\n";
        }
    }
}

void rank0ToRankMax(AdjacencyArray &adjArray){
    for(uint64_t nodeId = 0; nodeId < adjArray.width*adjArray.height; ++nodeId){
        if(adjArray.rank.at(nodeId) == 0){
            adjArray.rank.at(nodeId) = UINT64_MAX;
        } 
    }
    for(uint64_t nodeId = 0; nodeId < adjArray.width*adjArray.height; ++nodeId){
        if(adjArray.rank.at(nodeId) == 0){
            std::cout << "nodes with rank 0 still exist\n";
        } 
    }
}

void sortEdgesByRank(AdjacencyArray &adjArray){
    AdjacencyArray &original = *(new AdjacencyArray(adjArray));

    adjArray.edges.clear();
    adjArray.edgeIds.clear();
    adjArray.offsets.clear();
    adjArray.offsets.push_back(0);
    adjArray.distances.clear();

    uint64_t currOffset = 0;
    for(uint64_t nodeId = 0; nodeId < adjArray.width*adjArray.height; ++nodeId){
        for(uint64_t edgeIndex = adjArray.offsets.at(nodeId); edgeIndex < adjArray.offsets.at(nodeId+1); ++edgeIndex){
            //adjArray
            currOffset++;
        }
        adjArray.offsets.push_back(currOffset);
    }

    delete &original;
}

void printContractionPercentage(AdjacencyArray &adjArray){
    uint64_t numNodes = 0;
    uint64_t numContractedNodes = 0;
    for(uint64_t nodeId = 0; nodeId < adjArray.width*adjArray.height; ++nodeId){
        numNodes += !adjArray.nodes.at(nodeId);
        uint64_t rank = adjArray.rank.at(nodeId);
        numContractedNodes += (rank > 0) && (rank < UINT64_MAX);
    }
    std::cout << numNodes << " " << numContractedNodes << "\n";
    std::cout << "contraction percentage: " << ((double)numContractedNodes)/((double)numNodes) << "\n";
}

void checkMultiDijkstra(){
    AdjacencyArray arr;
    arr.width = 2;
    arr.height= 2;
    arr.nodes   =   {true,  true,   true,   true};
    arr.offsets =   {0,     2,      3,      4,      4};
    arr.edges =     {1,2,   3,      3};
    
    arr.distances = std::vector<uint64_t>(arr.edges.size(), 1);

    MultiDijkstra::Dijkstra md(arr);
    md.calculateDist(0, 3);
    std::vector<bool> isIndep(arr.nodes.size(), false);
    isIndep.at(2) = true;
    isIndep.at(1) = true;
    md.checkMultipleShortestPath(isIndep);
}

int main(){
    AdjacencyArray adjArray("data/CHAdjArray_curr_20.graph_2");
    std::cout << "Number of edges: " << adjArray.allEdgeInfo.size() << "\n";
    uint64_t numNodes = 0;
    for(uint64_t i = 0; i<adjArray.nodes.size(); ++i){
        numNodes += !adjArray.nodes.at(i);
    }
    std::cout << "Number of nodes: " << numNodes << "\n";
    std::cout << "Number of edges: " << adjArray.edges.size() << "\n";
    for(uint64_t i = 0; i<adjArray.allEdgeInfo.size(); ++i){
        uint64_t v1 = adjArray.allEdgeInfo.at(i).v1;
        uint64_t v2 = adjArray.allEdgeInfo.at(i).v2;
        if(adjArray.allEdgeInfo.at(i).shortcutPathEdges.size() > 2 || i > adjArray.allEdgeInfo.size()-5){
            std::cout << "EdgeId: " << adjArray.allEdgeInfo.at(i).edgeId << "\n";
            std::cout << "v1: " << v1 << " rank v1: " << adjArray.rank.at(v1) << "\n";
            std::cout << "v2: " << v2 << " rank v2: " << adjArray.rank.at(v2) << "\n";
            std::cout << "edge distance: " << adjArray.allEdgeInfo.at(i).edgeDistance << "\n";
            std::cout << "number of shortcut edges: " << adjArray.allEdgeInfo.at(i).shortcutPathEdges.size() << "\n";
            // test unpacking
            std::vector<uint64_t> unpackedEdgeIds;
            recursiveUnpack(adjArray.allEdgeInfo.at(i), unpackedEdgeIds, adjArray);
        }
    }

    uint64_t contractedNodes = 0;
    for(uint64_t nodeId=0; nodeId<adjArray.width*adjArray.height; ++nodeId){
        if(adjArray.rank.at(nodeId)>0 && adjArray.rank.at(nodeId)<UINT64_MAX){ 
            //std::cout << nodeId << " " << adjArray.rank.at(nodeId) << "\n"; 
            contractedNodes++;
        }
    }
    std::cout << "num contracted nodes: " << contractedNodes<< "\n";

    // checkRedundantEdges(adjArray);

    // checkUndirectedness(adjArray);

    // checkEdgeIdSorted(adjArray);

    // rank0ToRankMax(adjArray);

    // removeRedundantEdges(adjArray);

    // checkRedundantEdges(adjArray);

    // checkUndirectedness(adjArray);

    // checkEdgeIdSorted(adjArray);

    // adjArray.writeToDisk("data/CHAdjArray_currDup_110.graph_2");
    
    //checkMultiDijkstra();
    
    printContractionPercentage(adjArray);
}