#include <set>
# include "shortestPathUtils.cpp"

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
    AdjacencyArray original(adjArray);
    adjArray.edges.clear();
    adjArray.edgeIds.clear();
    adjArray.offsets.clear();
    adjArray.offsets.push_back(0);
    adjArray.distances.clear();

    uint64_t currentOffset = 0;
    for(uint64_t nodeId = 0; nodeId<adjArray.width*adjArray.height; ++nodeId){
        bool activeNode = !original.nodes.at(nodeId);
        
        if(activeNode){
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
        }
        adjArray.offsets.push_back(currentOffset);
    }
}

void checkUndirectedness(AdjacencyArray &adjArray){
    for(uint64_t nodeId = 0; nodeId < adjArray.width*adjArray.height; ++nodeId){
        if(!adjArray.nodes.at(nodeId)){
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
        }
    }
}

void checkEdgeIdSorted(AdjacencyArray &adjArray){
    for(uint64_t edgeId = 0; edgeId < adjArray.allEdgeInfo.size(); ++edgeId){
        if(adjArray.allEdgeInfo.at(edgeId).edgeId != edgeId){
            std::cout << "differs: " << edgeId << " " << adjArray.allEdgeInfo.at(edgeId).edgeId << "\n";
        }
    }
}

int main(){
    AdjacencyArray adjArray("data/CHAdjArray_54.graph_2");
    std::cout << "Number of edges: " << adjArray.allEdgeInfo.size() << "\n";
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
        if(adjArray.rank.at(nodeId)>0){ 
            //std::cout << nodeId << " " << adjArray.rank.at(nodeId) << "\n"; 
            contractedNodes++;
        }
    }
    std::cout << "num contracted nodes: " << contractedNodes<< "\n";

    checkRedundantEdges(adjArray);

    checkUndirectedness(adjArray);

    checkEdgeIdSorted(adjArray);

    removeRedundantEdges(adjArray);

    checkRedundantEdges(adjArray);

    checkUndirectedness(adjArray);

    checkEdgeIdSorted(adjArray);

    adjArray.writeToDisk("data/CHAdjArray_smallNoDuplicate.graph_2");
    
}