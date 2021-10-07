#include "shortestPathUtils.cpp"
#include "PathAlgorithm.cpp"

namespace Dijkstra{

/**
 * @brief HeapElement for Dijkstra implementation
 */
struct HeapElement {
    // for normal dijkstra, dist is the current distance to this node
    uint64_t nodeIdx, prev, dist;
    
    /**
     * @brief "reverse" comparison function turning max-heap into min-heap 
     * 
     * @param a         other HeapElement
     * @return true     if own distance is bigger than others distance,
     * @return false    otherwise
     */
    bool operator<(const HeapElement &a){
        return dist > a.dist;
    }
};

/**
 * @brief Corresponds to second more efficient implementation of the dijkstra algorithm in Dijkstra.cpp
 * Modified for one-to-many query
 * Similar implementation to https://github.com/Lesstat/dijkstra-performance-study/
 */
class Dijkstra: public PathAlgorithm{
    public:
        Dijkstra(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        uint64_t calculateDist(uint64_t endPoint);
        void reset();
        uint64_t getNumNodesPopped();
        void disableNode(uint64_t disabledNode);
        void resetDisabledNodes();
    private:
        uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
        uint64_t mainCalculationLoop();
        std::vector<uint64_t> distance;
        std::vector<HeapElement> heap;
        std::vector<uint64_t> visited;
        std::vector<uint64_t> prev;
        std::vector<bool> disabledNodes;
        AdjacencyArray &adjArray;
        uint64_t constLngDist;
        std::vector<uint64_t> constLatDist;
        uint64_t startPoint, endPoint, lastCalculatedDistance;
        uint64_t numNodesPopped;
};


/**
 * @brief Construct a new Dijkstra::Dijkstra object
 * 
 * Initialize datastructures
 * Initialize distances between nodes
 * 
 * @param array 
 */
Dijkstra::Dijkstra(AdjacencyArray &array) : adjArray(array), prev(array.width*array.height, UINT64_MAX){
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
}

/**
 * @brief reset datastructures to prepare for next call to calculateDist
 */
void Dijkstra::reset(){
    std::vector<uint64_t> initDist(adjArray.width*adjArray.height, UINT64_MAX);
    distance = std::move(initDist);
    // no need to reset prev
    heap.clear();
    // initially no nodes popped
    numNodesPopped = 0;
    // no nodes are disabled after reset
    resetDisabledNodes();
}

/**
 * @brief efficient dijkstra shortest-path implementation
 * 
 * Uses binary Min(Max)-heap for greedy node visitation strategy
 * 
 * For this to work, reset need to be called first
 * 
 * @param startPoint 
 * @param endPoint 
 * @return uint64_t 
 */
uint64_t Dijkstra::calculateDist(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;

    heap.push_back(HeapElement{startPoint, UINT64_MAX, 0});

    std::make_heap(heap.begin(), heap.end());
    
    return mainCalculationLoop();
}

/**
 * @brief calculate distance with same start point as the previous call to calculateDist(startPoint, endPoint)
 * 
 * do not call reset before next query
 * 
 * @param endPoint_ 
 * @return uint64_t 
 */
uint64_t Dijkstra::calculateDist(uint64_t endPoint_){
    endPoint = endPoint_;
    // if node is already settled, then return distance directly
    
    if(distance.at(endPoint) != UINT64_MAX){
        return distance.at(endPoint);
    }

    return mainCalculationLoop();
}

/**
 * @brief main dijkstra distance calculation function using min-heap
 * 
 * start node and end node have to be set prior to call
 * 
 * @return uint64_t calculated distance
 */
uint64_t Dijkstra::mainCalculationLoop(){
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

            // ignore disabled nodes
            if(disabledNodes.at(neighborIdx)){ continue; }

            // choose length of edge from precalculated lengths
            uint64_t edgeDist = adjArray.distances.at(currEdgeId);

            // calculate new distances
            uint64_t newNeighborDist = front.dist + edgeDist;
            uint64_t oldNeighborDist = distance.at(neighborIdx);

            // push updated node if distance is improved
            if(newNeighborDist<oldNeighborDist){
                heap.push_back(HeapElement{neighborIdx, front.nodeIdx, newNeighborDist});
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
uint64_t Dijkstra::getDist(){
    return lastCalculatedDistance;
}

/**
 * @brief retrieve path calculated by the last call to calculateDist
 * 
 * @param path 
 */
void Dijkstra::getPath(std::vector<uint64_t> &path){
    if(distance.at(endPoint) < UINT64_MAX){
        // build up path
        uint64_t currNode = endPoint;
        path.push_back(currNode);
        while(currNode != startPoint){
            currNode = prev.at(currNode);
            path.push_back(currNode);
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

uint64_t Dijkstra::getNumNodesPopped(){
    return numNodesPopped;
}

/**
 * @brief disable the specified node for use in queries
 * 
 * not used anymore
 * 
 * @param nodeId disabled node
 */
void Dijkstra::disableNode(uint64_t nodeId){
    disabledNodes.at(nodeId) = true;
}

/**
 * @brief clear disable nodes flag vector
 * 
 */
void Dijkstra::resetDisabledNodes(){
    disabledNodes = std::vector<bool>(adjArray.width*adjArray.height, false);
}


}