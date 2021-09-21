#include <ctime>

#include "shortestPathUtils.cpp"

int getPathInput(int argc, char** argv, std::string &out_path){
    if(argc > 1){
        out_path = std::string(argv[1]);
    }else{
        out_path = "data/planet.graph";
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

void fillContractionSet(AdjacencyArray &adjArray, std::vector<uint64_t> &out_contractionSet){
    uint64_t numDraws = 1000;
    uint64_t maxNodeId = adjArray.width*adjArray.height;
    std::vector<bool> inContractionSet(maxNodeId, false);
    
    std::srand(std::time(nullptr));
    
    uint64_t nodeIdx = 0;
    for(int i = 0; i<numDraws; ++i){
        
        // random position
        uint64_t draw = std::rand() % maxNodeId;
        
        // check if node is in water and not already in contraction set
        if(!inContractionSet.at(draw) && adjArray.nodes.at(draw)){
            
            // check if node is not adjacent to nodes in contraction set
            bool isAdjacent = false;
            for(uint64_t currentIdx : out_contractionSet){
                for(uint64_t currEdgeId = adjArray.offsets.at(draw); currEdgeId < adjArray.offsets.at(draw+1); currEdgeId++){
                    uint64_t adjacentId = adjArray.edges.at(currEdgeId);
                    if(adjacentId == currentIdx){
                        isAdjacent = true;
                        break;
                    }
                }
            }
            
            // if node is not adjacent, add to contraction set
            if(!isAdjacent){
                out_contractionSet.push_back(draw);
                inContractionSet.at(draw) = true;
            }
        
        }
    }
}

int main(int argc, char** argv){
    std::string path;
    if(getPathInput(argc, argv, path)==0){
        return 0;
    }

    AdjacencyArray adjArray(path);
    std::vector<uint64_t> nodeRank;
    
    
    // select set of independent nodes for contraction
    std::vector<uint64_t> contractionSet;
    fillContractionSet(adjArray, contractionSet);
    
    // create shortcuts to replace edges of nodes to remove
    // assign current rank to all contracted nodes
    // increment rank
    // mark edges of contracted nodes for removal
    // build intermediate graph and do further contractions
    // finally: G'' = (V, E u E')
}