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

int main(){
    AdjacencyArray adjArray("data/CHAdjArray.graph_2");
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
    for(uint64_t nodeId; nodeId<adjArray.width*adjArray.height; ++nodeId){
        if(adjArray.rank.at(nodeId)>0){ 
            //std::cout << nodeId << " " << adjArray.rank.at(nodeId) << "\n"; 
            contractedNodes++;
        }
    }
    std::cout << "num contracted nodes: " << contractedNodes<< "\n";


    
}