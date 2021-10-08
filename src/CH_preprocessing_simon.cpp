#include <ctime>
#include <numeric>
# include <omp.h>

#include "shortestPathUtils.cpp"
#include "Dijkstra_simon.cpp"
#include "MultiDijkstra_simon.cpp"



/**
 * @brief Create a 2D Vector Matrix object of type T
 * 
 * @tparam T                                type of the matrix elements
 * @param size_i                            size of the first dimenstion
 * @param size_j                            size of the second dimension
 * @return std::vector<std::vector<T>>      Matrix object
 */
template <typename T>
std::vector<std::vector<T>> createVectorMatrix(uint64_t size_i, uint64_t size_j){
    std::vector<std::vector<T>> matrix(size_i);
    for(uint64_t i = 0; i<size_i; ++i){
        matrix.at(i).resize(size_j);
    }
    return matrix;
}


/**
 * @brief Get the File Path Input string from the command line
 * 
 * @param argc      first parameter of main method, containing the number of command line arguments
 * @param argv      second parameter of main method, containing all command line arguments
 * @param out_path  output string reference with file path command line argument
 * @return int 
 */
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


/**
 * @brief From an AdjacencyArray input, generate a vector of all edges in the adjacency array
 * This gives each edge a unique identitiy that can be used to unpack shortcuts
 * 
 * @param adjArray          input AdjacencyArray
 * @param out_allEdges      output vector representing all edges in the AdjacencyArray
 */
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


/**
 * @brief Check if there are redundant edges in the AdjacencyArray and mark such edges for removal
 * 
 * @param adjArray          input AdjacencyArray
 * @param isEdgeRemoved     output binary flags for each edge
 */
void checkRedundantEdges(AdjacencyArray &adjArray, std::vector<bool> isEdgeRemoved){
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
    std::cout << numRemoved << " additional redundant edges removed\n";
}


/**
 * @brief Check if a node (idToCheck) is adjacent to a marked node (isIdInSet) in the AdjacencyArray 
 * 
 * @param adjArray      input AdjacencyArray
 * @param idToCheck     node id that is checked for adjacency to set     
 * @param isIdInSet     vector, where each node is flagged with true if it is in the set
 * @return true         if idToCheck is adjacent to any of the nodes in the set
 * @return false        otherwise
 */
bool isNodeAdjacentToSet(AdjacencyArray &adjArray, uint64_t idToCheck, std::vector<bool> &isIdInSet){
    for(uint64_t currEdgeIndex = adjArray.offsets.at(idToCheck); currEdgeIndex < adjArray.offsets.at(idToCheck+1); ++currEdgeIndex){
        uint64_t adjacentId = adjArray.edges.at(currEdgeIndex);
        if(isIdInSet.at(adjacentId)){
            return true;
        }
    }
    return false;
}


/**
 * @brief create a random set of independent nodes from the adjacency array that does not contain nodes that were already contracted
 * 
 * @param adjArray              input AdjacencyArray
 * @param allContractedIds      vector of node ids that are already contracted - newly marked nodes ids are added
 * @param isContracted          vector of boolean flags for each node, indicating if the node is already contracted - newly marked nodes are marked
 * @param out_newContractions   vector of node ids that are newly marked for contaction by this function
 */
void randomFillContractionSet(
    AdjacencyArray &adjArray, 
    std::vector<uint64_t> &allContractedIds, 
    std::vector<bool> &isContracted, 
    std::vector<uint64_t> &out_newContractions
    ){

    uint64_t numDraws = 200000;
    uint64_t nodeIdLimit = isContracted.size();
    
    uint64_t nodeIdx = 0;
    for(int i = 0; i<numDraws; ++i){
        
        // random position
        uint64_t draw = std::rand() % nodeIdLimit;
        
        // check if node is in water and not already in contraction set
        if(!isContracted.at(draw) && !adjArray.nodes.at(draw)){
            
            // check if node is not adjacent to nodes in contraction set
            bool isAdjacent = isNodeAdjacentToSet(adjArray, draw, isContracted);

            // if node is not adjacent, add to contraction set
            if(!isAdjacent){
                out_newContractions.push_back(draw);
                allContractedIds.push_back(draw);
                isContracted.at(draw) = true;
            }
        
        }
    }
}


/**
 * @brief heuristically create a set of independent nodes containing nodes with few edges from the adjacency array 
 *        the set does not contain nodes that were already contracted
 * 
 * @param adjArray                  input AdjacencyArray
 * @param fraction                  percentage of independent nodes with lowest cost to output
 * @param in_out_allContractedIds   vector of node ids that are already contracted - newly marked nodes ids are added
 * @param in_out_isContracted       vector of boolean flags for each node, indicating if the node is already contracted - newly marked nodes are marked
 * @param out_newContractions       vector of node ids that are newly marked for contaction by this function
 */
void numEdgeFillContractionSet(
    AdjacencyArray &adjArray,
    double fraction, 
    std::vector<uint64_t> &in_out_allContractedIds, 
    std::vector<bool> &in_out_isContracted, 
    std::vector<uint64_t> &out_newContractions
    ){

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


/**
 * @brief contract the specified node and add shortcuts
 * 
 * @param contractedNodeId          id of node to contract
 * @param workArray                 input AdjacencyArray
 * @param currentRank               current rank of the contraction process
 * @param dijkstra                  dijkstra routing instance for shortest path calculation
 * @param multiDijkstra             multiDijkstra routing instance for checking if there exist multiple shortest paths
 * @param isNodeContracted          vector of boolean flags for each node indicating whether the node is contracted (used in multiDijkstra)
 * @param out_isEdgeRemoved         vector to add removed edges of the current node
 * @param out_removedEdgeIndices    vector that edgeIndices removed in this round are added to
 * @param out_allEdges              vector that newly created edges are added to
 * @param out_shortcutEdges         vector that newly created edges are added to
 */
void contractNode(  uint64_t contractedNodeId, AdjacencyArray &workArray, uint64_t currentRank, 
                    Dijkstra::Dijkstra &dijkstra, MultiDijkstra::Dijkstra &multiDijkstra, std::vector<bool> &isNodeContracted,
                    std::vector<bool> &out_isEdgeRemoved, std::vector<uint64_t> &out_removedEdgeIndices,
                    std::vector<Edge> &out_allEdges, std::vector<Edge> &out_shortcutEdges){

    // set the rank of contractedNodes to current rank
    workArray.rank.at(contractedNodeId) = currentRank;


    // find all neighbor ids and all edges (to and from contracted node)
    std::vector<uint64_t> adjacentIds;
    std::vector<uint64_t> forwardEdgeIds;       // used to avoid duplicate edges
    std::vector<uint64_t> forwardEdgeIndices;   // used to avoid duplicate edges
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
    
    // process edges to contracted node
    std::vector<uint64_t> reverseEdgeIds(adjacentIds.size(), UINT64_MAX); // used to avoid duplicate edges
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


    //if(edgeIds.size() != 2*adjacentIds.size()){ std::cout << "more/less edges" << edgeIds.size() << " " << adjacentIds.size() << "\n"; }

    //calculate distances between all adjacent nodes
    std::vector<std::vector<uint64_t>> distanceUWMatrix = createVectorMatrix<uint64_t>(adjacentIds.size(), adjacentIds.size());
    std::vector<std::vector<bool>> multipleShortcutsMatrix = createVectorMatrix<bool>(adjacentIds.size(), adjacentIds.size());
    for(uint64_t uIndex = 0; uIndex<adjacentIds.size(); ++uIndex){
        uint64_t uId = adjacentIds.at(uIndex);
        dijkstra.reset(); // do not disable any nodes
        multiDijkstra.reset();
        if(uIndex+1<adjacentIds.size()){
            dijkstra.calculateDist(uId, adjacentIds.at(uIndex+1)); // set correct start point
            multiDijkstra.calculateDist(uId, adjacentIds.at(uIndex+1));
        } 
        for(uint64_t wIndex = uIndex+1; wIndex<adjacentIds.size(); ++wIndex){
            uint64_t wId = adjacentIds.at(wIndex);
            uint64_t distance = dijkstra.calculateDist(wId); // only query with end point (one-to-many)
            multiDijkstra.calculateDist(wId);
            bool multipleShortcuts = multiDijkstra.checkMultipleShortestPath(isNodeContracted);
            distanceUWMatrix.at(uIndex).at(wIndex) = distance;
            distanceUWMatrix.at(wIndex).at(uIndex) = distance;
            multipleShortcutsMatrix.at(uIndex).at(wIndex) = multipleShortcuts;
            multipleShortcutsMatrix.at(wIndex).at(uIndex) = multipleShortcuts;
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

            // read distances from matrix: order of forward/backward indices and i/j does not matter since graph is undirected 
            // -> always has edges in both directions of same length between nodes 
            uint64_t distanceUW = distanceUWMatrix.at(i).at(j);
            uint64_t distanceUVW = workArray.distances.at(startNodeEdgeIndex) + workArray.distances.at(endNodeEdgeIndex);

            bool multipleShortestPaths = multipleShortcutsMatrix.at(i).at(j); //md.checkMultipleShortestPath(con);
            
            //                              condition for avoiding duplicate shortcuts
            //                       there is a shorter path   OR  there exist shortest paths besides UVW 
            bool noShortcutNeeded = (distanceUW < distanceUVW) || multipleShortestPaths;

            //if((!(distanceUW < distanceUVW)) && multipleShortestPaths){std::cout << "shortcut saved\n";}

            // if deleted edges are potentially part of a shortest path, then add shortcut edge
            if(!noShortcutNeeded){
                // add shortcut to list
                uint64_t edgeDistance = distanceUVW;
                uint64_t shortcutId = out_allEdges.size();
                std::vector<uint64_t> edgePathUVW {reverseEdgeIds.at(i), forwardEdgeIds.at(j)};

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


/**
 * @brief heuristically create a set of independent nodes containing nodes with low edge difference from the adjacency array 
 *        the set does not contain nodes that were already contracted
 * 
 * @param adjArray                  input AdjacencyArray
 * @param fraction                  percentage of independent nodes with lowest cost to output
 * @param in_out_allContractedIds   vector of node ids that are already contracted - newly marked nodes ids are added
 * @param in_out_isContracted       vector of boolean flags for each node, indicating if the node is already contracted - newly marked nodes are marked
 * @param out_newContractions       vector of node ids that are newly marked for contaction by this function
 */
void edgeDifferenceFillContractionSet(
    AdjacencyArray &adjArray,
    double fraction, 
    std::vector<uint64_t> &in_out_allContractedIds, 
    std::vector<bool> &in_out_isContracted, 
    std::vector<uint64_t> &out_newContractions
    ){

    std::vector<uint64_t> independentSet;
    std::vector<bool> isInIndependentSet(adjArray.width*adjArray.height, false);

    uint64_t nodeIdLimit = in_out_isContracted.size();
    uint64_t initialId = 0;
    do{
        initialId = std::rand() % nodeIdLimit;
    }while(in_out_isContracted.at(initialId) || adjArray.nodes.at(initialId));

    independentSet.push_back(initialId);
    isInIndependentSet.at(initialId) = true;

    for(uint64_t nodeId=0; nodeId<adjArray.width*adjArray.height; ++nodeId){
        
        // check if node is in water and not already in contraction set
        if(!in_out_isContracted.at(nodeId) && !adjArray.nodes.at(nodeId) && !isInIndependentSet.at(nodeId)){
            
            // check if node is not adjacent to nodes in contraction set
            bool isAdjacent = isNodeAdjacentToSet(adjArray, nodeId, isInIndependentSet);

            // if node is not adjacent, add to independent set
            if(!isAdjacent){
                independentSet.push_back(nodeId);
                isInIndependentSet.at(nodeId) = true;
            }
        }
    }

    std::cout << "contraction candidates: " << independentSet.size() << "\n";

    std::vector<std::pair<uint64_t, int>> nodeWithCost;

    #pragma omp parallel num_threads(10)
    {
        MultiDijkstra::Dijkstra multiDijkstra(adjArray);
        multiDijkstra.stepLimit = 0;
        std::vector<bool> isNodeInSet(adjArray.width*adjArray.height, true); // no additional shortcuts check (dont know other nodes in final indep set)
        std::vector<bool> isEdgeRemoved(adjArray.edges.size(), false); // this is not read in contract node
        std::vector<uint64_t> removedEdgeIndices;
        std::vector<Edge> allEdges;
        std::vector<Edge> shortcutEdges;
        Dijkstra::Dijkstra dijkstra(adjArray);
        #pragma omp for
        for(uint64_t nodeIndex=0; nodeIndex<independentSet.size(); ++nodeIndex){
            
            uint64_t nodeId = independentSet.at(nodeIndex);

            // check if node is in water and not already in contraction set
            if(!in_out_isContracted.at(nodeId) && !adjArray.nodes.at(nodeId)){

                // calculate criterion
                removedEdgeIndices.clear();
                allEdges.clear();
                shortcutEdges.clear();
                contractNode(nodeId, adjArray, adjArray.rank.at(nodeId), dijkstra, multiDijkstra, isNodeInSet, isEdgeRemoved, removedEdgeIndices, allEdges, shortcutEdges);
                int criterion = shortcutEdges.size() - removedEdgeIndices.size();

                //std::cout << nodeId << "  " << criterion << "\n";
                #pragma omp critical
                {
                    // add to list
                    nodeWithCost.push_back(std::pair<uint64_t, uint64_t>(nodeId, criterion));
                    if(nodeWithCost.size() % 1000 == 0){
                        std::cout << "cost calculation progress " << ((double)nodeWithCost.size())/independentSet.size()  << "\t\r" << std::flush;
                    }
                }
            }
        }
    }
    
    std::cout << "contraction candidates: " << nodeWithCost.size() << "\n";

    std::sort(  nodeWithCost.begin(), nodeWithCost.end(),
                [](std::pair<uint64_t,int> &p1, std::pair<uint64_t,int> &p2){ return p1.second < p2.second; });
    

    // only take the best fraction of the independent set according to the cost
    for(uint64_t nodeIndex=0; nodeIndex<((uint64_t)(nodeWithCost.size()*fraction)); ++nodeIndex){
        uint64_t nodeId = nodeWithCost.at(nodeIndex).first;
        out_newContractions.push_back(nodeId);
        in_out_allContractedIds.push_back(nodeId);
        in_out_isContracted.at(nodeId) = true;
    }
    std::cout << "num nodes contracted: " << out_newContractions.size() << "  best value: " << nodeWithCost.at(0).second <<  "\n";
}

int main(int argc, char** argv){
    std::string path;
    if(getFilePathInput(argc, argv, path)==0){
        return 0;
    }

    std::srand(std::time(nullptr));
    
    // load non-CH AdjacencyArray from disk and set initial rank to UINT_MAX
    AdjacencyArray adjArray(path);
    adjArray.rank = std::vector<uint64_t>(adjArray.width*adjArray.height, UINT64_MAX); // set rank from 0 to UINT64_MAX because uncontracted core has highest rank

    // create final array in dynamic memory
    AdjacencyArray &finalArray = *(new AdjacencyArray(adjArray));

    // create vectors to keep track of which nodes are contracted
    std::vector<bool> isContracted(adjArray.width*adjArray.height, false);
    std::vector<uint64_t> contractedNodeIds;
    
    // create vector containing all edges
    std::vector<Edge> allEdges;
    extractAllEdges(adjArray, allEdges);
    
    // initialize rank to 1
    uint64_t currentRank = 1;
    uint64_t round = 1;

    // copy original adjacency array
    AdjacencyArray workArray = AdjacencyArray(adjArray);

    std::cout << "initial number of uncontracted nodes: " << 
                    workArray.width*workArray.height - std::accumulate(workArray.nodes.begin(), workArray.nodes.end(), 0) << "\n";
    
    // while there is an uncontracted core bigger than the threshold, start new contraction rounds
    while((workArray.width*workArray.height - std::accumulate(workArray.nodes.begin(), workArray.nodes.end(), 0)) > 550){
        std::cout << "round " << round++ << "\n";

        Dijkstra::Dijkstra dijkstra(workArray);
        MultiDijkstra::Dijkstra multiDijkstra(workArray);
        multiDijkstra.stepLimit = 1000;
        
        std::vector<Edge> shortcutEdges;
        std::vector<uint64_t> removedEdgeIndices;
        std::vector<bool> isEdgeRemoved(workArray.edges.size(), false); // flag for edge indices if they are removed

        // select set of independent nodes for contraction
        std::vector<uint64_t> newContractions;
        edgeDifferenceFillContractionSet(workArray, 0.5, contractedNodeIds, isContracted, newContractions);

        uint64_t roundProgress = 0;

        // create shortcuts to replace edges of nodes to contract/remove
        for(uint64_t contractedNodeId : newContractions){
            
            if(roundProgress++%1024==0){
                std::cout << "round " << round-1 << " progress " << ((double)roundProgress)/newContractions.size()  << "\t\r" << std::flush;
            }
            
            contractNode(   contractedNodeId, workArray, currentRank,
                            dijkstra, multiDijkstra, isContracted, isEdgeRemoved, removedEdgeIndices, allEdges, shortcutEdges);
        }

        // check for duplicate shortcut edges and mark them
        std::vector<bool> isShortcutIdRemoved(allEdges.size(), false); // vector for duplicate shortcut edge ids (instead of edge indices in isEdgeRemoved)
        checkRedundantEdges(adjArray, isShortcutIdRemoved);


        std::cout << "rebuild phase\n";
        
        // build intermediate graph
        AdjacencyArray tempArray = AdjacencyArray(workArray);
        
        // adjust temp array
        tempArray.edges.clear();
        tempArray.edgeIds.clear();
        tempArray.distances.clear();
        tempArray.offsets.clear();
        tempArray.offsets.push_back(0);
        
        // sort shortcuts for binary search
        std::vector<Edge> shortcutSortedV1 = copyEdgeVector(shortcutEdges);
        std::sort(  
            shortcutSortedV1.begin(), shortcutSortedV1.end(), 
            [](Edge& a, Edge& b) {return a.v1 < b.v1; });

        // insert edges for all nodes
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
                    if(!isShortcutIdRemoved.at(shortcutEdgeIt->edgeId)){ // ignore marked shortcut edges
                        tempArray.edges.push_back(shortcutEdgeIt->v2);
                        tempArray.edgeIds.push_back(shortcutEdgeIt->edgeId);
                        tempArray.distances.push_back(shortcutEdgeIt->edgeDistance);
                        currentOffset++;
                    }
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
        // final array is constructed and saved in every round
        
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
    }
    
    delete &finalArray;
}