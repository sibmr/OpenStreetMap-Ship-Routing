#include <ctime>
#include <numeric>

#include "shortestPathUtils.cpp"

struct ShortcutEdge
{
    uint64_t v1;
    uint64_t v2;
    uint64_t edgeDistance;
    std::vector<uint64_t> removedEdgeIds;
};


int getPathInput(int argc, char** argv, std::string &out_path){
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

bool isNodeAdjacentToSet(AdjacencyArray &adjArray, uint64_t idToCheck, std::vector<uint64_t> &idSet){
    for(uint64_t currentIdx : idSet){
        for(uint64_t currEdgeId = adjArray.offsets.at(idToCheck); currEdgeId < adjArray.offsets.at(idToCheck+1); currEdgeId++){
            uint64_t adjacentId = adjArray.edges.at(currEdgeId);
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

    uint64_t numDraws = 1000;
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
    if(getPathInput(argc, argv, path)==0){
        return 0;
    }

    std::srand(std::time(nullptr));
    
    AdjacencyArray adjArray(path);
    std::vector<bool> isContracted(adjArray.width*adjArray.height, false);
    std::vector<uint64_t> contractedNodeIds;
    uint64_t currentRank = 0;

    uint64_t numWaterNodes = 0;
    for(bool onLand : adjArray.nodes){
        numWaterNodes += !onLand;
    }

    // copy original adjacency array
    AdjacencyArray workArray = AdjacencyArray(adjArray);

    // while there are still uncontracted nodes
    for(uint64_t round = 0; round<5; ++round){
        
        std::vector<ShortcutEdge> shortcutEdges;
        std::vector<uint64_t> removedEdges;

        // select set of independent nodes for contraction
        std::vector<uint64_t> newContractions;
        fillContractionSet(workArray, contractedNodeIds, isContracted, newContractions);

        // create shortcuts to replace edges of nodes to contract/remove
        for(uint64_t contractedNodeId : newContractions){
            
            // set the rank of contractedNodes to current rank
            workArray.rank.at(contractedNodeId) = currentRank;

            // find all neighbor ids
            std::vector<uint64_t> adjacentIds; 
            std::vector<uint64_t> edgeIds;
            for(uint64_t currEdgeId = workArray.offsets.at(contractedNodeId); currEdgeId < workArray.offsets.at(contractedNodeId+1); currEdgeId++){
                uint64_t adjacentId = workArray.edges.at(currEdgeId);
                adjacentIds.push_back(adjacentId);
                edgeIds.push_back(currEdgeId);
            }

            //insert shortcuts
            for(uint64_t i = 0; i<adjacentIds.size(); ++i){
                for(uint64_t j = i+1; j<adjacentIds.size(); ++j){
                    uint64_t startNode = adjacentIds.at(i);
                    uint64_t endNode = adjacentIds.at(j);
                    std::vector<uint64_t> replacedEdges{edgeIds.at(i), edgeIds.at(j)};
                    
                    // TODO: check if shortcut is necessary between start and end node
                    
                    // add shortcut to list
                    uint64_t edgeDistance = workArray.distances.at(edgeIds.at(i)) + workArray.distances.at(edgeIds.at(j));
                    shortcutEdges.push_back(ShortcutEdge{startNode, endNode, edgeDistance, replacedEdges});
                }
            }

            //mark adjacent edges of contractedNode for removal
            removedEdges.insert(removedEdges.end(), edgeIds.begin(), edgeIds.end());
        }
        
        // build intermediate graph
        AdjacencyArray tempArray = AdjacencyArray(workArray);
        
        // increment rank
        currentRank++;
    }
    // finally: G'' = (V, E u E')
}