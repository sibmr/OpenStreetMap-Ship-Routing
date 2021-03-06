#include <numeric>

#include "shortestPathUtils.cpp"
#include "PathAlgorithm.cpp"

namespace CH_Astar{

/**
 * @brief HeapElement for A_star implementation
 */
struct HeapElement {
    // heuristic dist is f(x) = g(x) + h(x), dist is g(x)
    uint64_t nodeIdx, prev, heuristic_dist, dist;
    
    /**
     * @brief "reverse" comparison function turning max-heap into min-heap 
     * 
     * @param a         other HeapElement
     * @return true     if own distance is bigger than others distance,
     * @return false    otherwise
     */
    bool operator<(const HeapElement &a){
        return heuristic_dist > a.heuristic_dist;
    }
};

/**
 * @brief Implementation of the A* Algorithm for Contraction Hierarchies
 * 
 */
class A_star: public PathAlgorithm{
    public:
        A_star(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        void reset();
        virtual uint64_t getHeuristic(AdjacencyArray &adjArray, uint64_t firstNodeIdx, uint64_t secondNodeIdx);
        virtual uint64_t getNumNodesPopped();

        uint64_t constLngDist;
        std::vector<uint64_t> constLatDist;
    private:
        uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
        void recursiveUnpackEdge(uint64_t edgeId, std::vector<uint64_t> &unpackedEdges);
        uint64_t getEdgeIdBetween(uint64_t nodeId1, uint64_t nodeId2);
        std::vector<uint64_t> distance;
        std::vector<HeapElement> heap;
        std::vector<uint64_t> visited;
        std::vector<uint64_t> prev;
        AdjacencyArray &adjArray;
        uint64_t startPoint, endPoint, lastCalculatedDistance;
        uint64_t numNodesPopped;
        std::vector<bool> isEdgeMarked;
        std::vector<bool> isNodeMarked;
        std::vector<uint64_t> edgeReverseIndex;
        
};


/**
 * @brief Construct a new A_star object
 * 
 * Initialize datastructures
 * Initialize distances between nodes
 * 
 * @param array 
 */
A_star::A_star(AdjacencyArray &array) : adjArray(array), prev(array.width*array.height, UINT64_MAX){
    reset();

    // calculate distances between nodes:
    // constLngDist is the distance between nodes with the same longitude but different latitude
    // constLatDist is the distance between nodes with the same latitude but different logitude
    // distance between (i,0) and (i,1)
    constLngDist = nodeDistance(adjArray, 0, 1);
    // for each constant latitude "ring" around the globe, the distance is different
    // mirroring is disregarded
    for(uint64_t i = 0; i<adjArray.height; ++i){
        // distance between (0,i) and (1,i)
        constLatDist.push_back(nodeDistance(adjArray, i, adjArray.height+i));
    }

    // fill edge/offset reverse
    edgeReverseIndex = std::vector<uint64_t>(adjArray.edges.size(), UINT64_MAX);
    for(uint64_t nodeId = 0; nodeId < adjArray.width*adjArray.height; ++nodeId){
        for(uint64_t currEdgeId = adjArray.offsets.at(nodeId); currEdgeId < adjArray.offsets.at(nodeId+1); ++currEdgeId){
            uint64_t adjacentId = adjArray.edges.at(currEdgeId);
            for(uint64_t backEdgeIndex = adjArray.offsets.at(adjacentId); backEdgeIndex < adjArray.offsets.at(adjacentId+1); ++backEdgeIndex){
                if(adjArray.edges.at(backEdgeIndex) == nodeId){
                    edgeReverseIndex.at(currEdgeId) = backEdgeIndex;
                }
            }
            if(edgeReverseIndex.at(currEdgeId) == UINT64_MAX){ std::cout << "no back edge \n"; }
            //std::cout << edgeReverseIndex.at(currEdgeId) << " " << currEdgeId << "\n";
        }
    }
}


/**
 * @brief reset datastructures to prepare for next call to calculateDist
 */
void A_star::reset(){
    std::vector<uint64_t> initDist(adjArray.width*adjArray.height, UINT64_MAX);
    distance = std::move(initDist);
    isEdgeMarked = std::vector<bool>(adjArray.edges.size(), false);
    isNodeMarked = std::vector<bool>(adjArray.width*adjArray.height, false);
    // no need to reset prev
    heap.clear();
    // initially no nodes popped
    numNodesPopped = 0;
}


/**
 * @brief for two nodes, calculate the heuristic underapproximation of their path distance
 * 
 * @param adjArray          input AdjacencyArray
 * @param firstNodeIdx      star node
 * @param secondNodeIdx     goal node
 * @return uint64_t         great-circle distance between the two nodes
 */
uint64_t A_star::getHeuristic(AdjacencyArray &adjArray, uint64_t firstNodeIdx, uint64_t secondNodeIdx){
    return nodeDistance(adjArray, firstNodeIdx, secondNodeIdx);
}


/**
 * @brief calculate the distance between two nodes using the astar algorithm for contraction hierarchies
 * 
 * Uses binary Min(Max)-heap for greedy node visitation strategy
 * 
 * For this to work, reset need to be called first
 * 
 * @param startPoint 
 * @param endPoint 
 * @return uint64_t 
 */
uint64_t A_star::calculateDist(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;

    // do depth first search to mark edges on all downward paths towards the goal node
    std::vector<uint64_t> stack;
    stack.push_back(endPoint);
    while(!stack.empty()){
        uint64_t nodeId = stack.back();
        stack.pop_back();
        if(!isNodeMarked.at(nodeId)){
            isNodeMarked.at(nodeId) = true;
            for(uint64_t currEdgeId = adjArray.offsets.at(nodeId); currEdgeId < adjArray.offsets.at(nodeId+1); ++currEdgeId){
                uint64_t adjacentId = adjArray.edges.at(currEdgeId);
                if(     !isEdgeMarked.at(edgeReverseIndex.at(currEdgeId)) 
                    &&  adjArray.rank.at(adjacentId) > adjArray.rank.at(nodeId))
                {
                    isEdgeMarked.at(edgeReverseIndex.at(currEdgeId)) = true;
                    if(adjArray.edges.at(edgeReverseIndex.at(currEdgeId)) != nodeId){
                        std::cout << "is not nodeId\n";
                    }
                    // isEdgeMarked.at(currEdgeId) = true;
                    stack.push_back(adjacentId);
                }
            }
        }
       
    }
    // std::cout << "backward marking finished\n";
    
    // do forward astar

    heap.push_back(HeapElement{startPoint, UINT64_MAX, getHeuristic(adjArray, startPoint, endPoint), 0});

    std::make_heap(heap.begin(), heap.end());

    HeapElement front;

    while(true){

        // no path found
        if(heap.empty()){
            return UINT64_MAX;
        }

        // get heap front
        std::pop_heap(heap.begin(), heap.end());
        front = heap.back();
        heap.pop_back();
        numNodesPopped++;

        // avoid duplicate nodes (nodes that were already visited, indicated by higher distance)
        if(front.dist >= distance.at(front.nodeIdx)){
            continue;
        }

        // update distance and previous node of current node
        distance.at(front.nodeIdx) = front.dist;
        prev.at(front.nodeIdx) = front.prev;
        visited.push_back(front.nodeIdx);

        // iterate over edges of current node
        for(uint64_t currEdgeId = adjArray.offsets.at(front.nodeIdx); currEdgeId < adjArray.offsets.at(front.nodeIdx+1); currEdgeId++){
            
            // get id of adjacent node for current edge (neighboring node)
            uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

            // choose length of edge from precalculated lengths
            uint64_t edgeDist = adjArray.distances.at(currEdgeId);

            // calculate new distances
            uint64_t newNeighborDist = front.dist + edgeDist;
            uint64_t oldNeighborDist = distance.at(neighborIdx);
            uint64_t neighborRank = adjArray.rank.at(neighborIdx);
            uint64_t currentRank = adjArray.rank.at(front.nodeIdx);
            if(neighborRank == 0){ neighborRank = 1000; }
            if(currentRank == 0){ currentRank = 1000; }

            // update node distance if it improves
            if(newNeighborDist<oldNeighborDist && (neighborRank >= currentRank || isEdgeMarked.at(currEdgeId))){
                heap.push_back(
                    HeapElement{
                        neighborIdx, 
                        front.nodeIdx, 
                        newNeighborDist + getHeuristic(adjArray, neighborIdx, endPoint), // f(x) = g(x) + h(x)
                        newNeighborDist}); // g(x)
                std::push_heap(heap.begin(), heap.end());
            }
        }

        if(front.nodeIdx == endPoint){
            lastCalculatedDistance = distance.at(front.nodeIdx);
            return lastCalculatedDistance;
        }
    }

}


/**
 * @brief returns distance of last call to calculateDist
 * 
 * @return uint64_t 
 */
uint64_t A_star::getDist(){
    return lastCalculatedDistance;
}


/**
 * @brief unpack the node-path that is represented by the edgeId 
 * 
 * @param edgeId            id of edge to unpack
 * @param unpackedNodes     unpacked vector of nodes (still contains duplicate nodes)
 */
void A_star::recursiveUnpackEdge(uint64_t edgeId, std::vector<uint64_t> &unpackedNodes){
    Edge &edge = adjArray.allEdgeInfo.at(edgeId);
    if(edge.shortcutPathEdges.size() > 0){
        for(uint64_t edgeId : edge.shortcutPathEdges){
            recursiveUnpackEdge(edgeId, unpackedNodes);
        }
    }else{
        unpackedNodes.push_back(edge.v1);
        unpackedNodes.push_back(edge.v2);
    }
}


/**
 * @brief given two nodes, calculate the id of the edge between them, if there is one
 * 
 * @param nodeId1   first node id
 * @param nodeId2   second node id
 * @return uint64_t the id of the edge between first and second node, if there is none: UINT64_t instead
 */
uint64_t A_star::getEdgeIdBetween(uint64_t nodeId1, uint64_t nodeId2){
    for(uint64_t currEdgeIndex = adjArray.offsets.at(nodeId1); currEdgeIndex < adjArray.offsets.at(nodeId1+1); ++currEdgeIndex){
        uint64_t adjacentNodeId = adjArray.edges.at(currEdgeIndex);
        if(adjacentNodeId == nodeId2){
            return adjArray.edgeIds.at(currEdgeIndex);
        }
    }
    return UINT64_MAX;
}


/**
 * @brief retrieve path calculated by the last call to calculateDist by unpacking shortcuts
 * 
 * @param path  path of nodes from start to goal node, empty vector if there is no path
 */
void A_star::getPath(std::vector<uint64_t> &path){
    if(distance.at(endPoint) < UINT64_MAX){
        
        // build up path
        uint64_t prevNode = endPoint;
        uint64_t currNode = UINT64_MAX;
        while(currNode != startPoint){
            currNode = prev.at(prevNode);
            
            std::cout << "p1: " << startPoint << " " << endPoint << "\n";
            std::cout << currNode << "\n";
            
            uint64_t currEdgeId = getEdgeIdBetween(prevNode, currNode);
            std::vector<uint64_t> unpackedNodeIds;
            recursiveUnpackEdge(currEdgeId, unpackedNodeIds);
            for(uint64_t i = 0; i<unpackedNodeIds.size(); i+=2){
                path.push_back(unpackedNodeIds.at(i));
            }
            path.push_back(unpackedNodeIds.at(unpackedNodeIds.size()-1));
            
            prevNode = currNode;
        }

        // print path
        std::cout << "Path:" << std::endl;
        for (std::vector<uint64_t>::iterator it = path.begin(); it != path.end(); ++it) {
            std::cout << *it << " ";
        }
        std::cout <<  std::endl;
        std::cout << "dist: " << distance.at(endPoint)/1000 << "km" << std::endl;
    }else{
        std::cout << "no path found" << std::endl;
    }
    
}


/**
 * @brief return the number of nodes popped from the heap during shortest path calculation
 * 
 * @return uint64_t number of nodes
 */
uint64_t A_star::getNumNodesPopped(){
    return numNodesPopped;
}



/**
 * @brief A_star with a thighter lower bound on the path distance as heuristic
 * 
 */
class A_star_rectangular : public A_star{
    
    public:
        A_star_rectangular(AdjacencyArray &array) : A_star(array) {}
        uint64_t getHeuristic(AdjacencyArray &adjArray, uint64_t firstNodeIdx, uint64_t secondNodeIdx);
    
};


/**
 * @brief heuristic distance between two nodes similar to the manhattan distance
 * the heuristic calculates min(start-southpole-goal,start-northpole-goal,shortest path around rectange on the surface of the sphere)
 * the calculation considers longitudinal wraparound
 * this only works because of the 4-connected grid with nodes equally distributed over longitude and latitude
 * 
 * @param adjArray          input AdjacencyArray
 * @param firstNodeIdx      node index of the first node
 * @param secondNodeIdx     node index of the second node
 * @return uint64_t         underapproximation of path distance between nodes
 */
uint64_t A_star_rectangular::getHeuristic(AdjacencyArray &adjArray, uint64_t firstNodeIdx, uint64_t secondNodeIdx){
    uint64_t lngFirst, latFirst, lngSecond, latSecond;

    nodeIdToArrayIdx(adjArray, firstNodeIdx, lngFirst, latFirst);
    nodeIdToArrayIdx(adjArray, secondNodeIdx, lngSecond, latSecond);

    // makes sure lngFirst <= lngSecond
    if(lngFirst > lngSecond){
        std::swap(lngFirst, lngSecond);
    }

    // makes sure latFirst <= latSecond
    if(latFirst > latSecond){
        std::swap(latFirst, latSecond);
    }

    return (
        std::min(
            std::min(
            
            // distance over south pole
            (latFirst+latSecond)*constLngDist,  
            
            // distance over north pole
            ((adjArray.height-latFirst)+(adjArray.height-latSecond))*constLngDist  
            ),

            // shortest "manhattan" distance: latitude difference + minimum longitude difference (wraparound) on the smaller ring around the globe
            (latSecond-latFirst)*constLngDist + 
            std::min((lngSecond-lngFirst), (lngFirst+adjArray.width) - lngSecond) * std::min(constLatDist[latFirst], constLatDist[latSecond])
        )
    );
}


}