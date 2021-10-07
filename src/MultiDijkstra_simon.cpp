#include "shortestPathUtils.cpp"
#include "PathAlgorithm.cpp"

namespace MultiDijkstra{

/**
 * @brief HeapElement for Dijkstra implementation
 */
struct HeapElement {
    // for normal dijkstra, heuristic_dist is the current distance to this node
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
 * @brief modified dijkstra algorithm to check if there exist multiple shortest paths
 * 
 */
class Dijkstra: public PathAlgorithm{
    public:
        Dijkstra(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        bool checkMultipleShortestPath(std::vector<bool> &isNodeInIndepententSet);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        uint64_t calculateDist(uint64_t endPoint);
        void reset();
        uint64_t getNumNodesPopped();
        void disableNode(uint64_t disabledNode);
        void resetDisabledNodes();
        uint64_t stepLimit;
    private:
        uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
        uint64_t mainCalculationLoop();
        std::vector<uint64_t> distance;
        std::vector<HeapElement> heap;
        std::vector<uint64_t> visited;
        std::vector<std::vector<uint64_t>> prev; // previous nodes are now a vector
        std::vector<bool> disabledNodes;
        AdjacencyArray &adjArray;
        uint64_t constLngDist;
        std::vector<uint64_t> constLatDist;
        uint64_t startPoint, endPoint, lastCalculatedDistance;
        uint64_t numNodesPopped;
};


/**
 * @brief Construct a new MultiDijkstra object
 * 
 * Initialize datastructures
 * 
 * @param array 
 */
Dijkstra::Dijkstra(AdjacencyArray &array) : adjArray(array){
    reset();
    stepLimit = UINT64_MAX;
}


/**
 * @brief reset datastructures to prepare for next call to calculateDist
 */
void Dijkstra::reset(){
    prev = std::vector<std::vector<uint64_t>>();
    prev.resize(adjArray.width*adjArray.height);
    std::vector<uint64_t> initDist(adjArray.width*adjArray.height, UINT64_MAX);
    distance = std::move(initDist);
    heap.clear();
    // initially no nodes popped
    numNodesPopped = 0;
    // no nodes are disabled after reset
    resetDisabledNodes();
}


/**
 * @brief dijkstra shortest-path implementation that saves all possible previous candidates for each node
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
    lastCalculatedDistance = UINT64_MAX;
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
        if(heap.empty() || numNodesPopped > stepLimit){
            return lastCalculatedDistance;
        }

        // get heap front
        std::pop_heap(heap.begin(), heap.end());
        front = heap.back();
        heap.pop_back();
        numNodesPopped++;

        // avoid duplicate nodes (nodes that were already visited, indicated by higher distance)
        if(front.heuristic_dist > distance.at(front.nodeIdx)){
            continue;
        }
        // if node has the same distance as previous node, add to set of previous nodes
        else if(front.heuristic_dist == distance.at(front.nodeIdx)){
            prev.at(front.nodeIdx).push_back(front.prev);
            continue;
        }
        // update distance and previous node of current node
        distance.at(front.nodeIdx) = front.heuristic_dist;
        // if there is a new shortest distance for the current node, clear previous, longer path neighbors
        prev.at(front.nodeIdx).clear();
        prev.at(front.nodeIdx).push_back(front.prev);
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
            uint64_t newNeighborDist = front.heuristic_dist + edgeDist;
            uint64_t oldNeighborDist = distance.at(neighborIdx);

            // push updated node if distance is improved
            if(newNeighborDist<oldNeighborDist){
                heap.push_back(HeapElement{neighborIdx, front.nodeIdx, newNeighborDist});
                std::push_heap(heap.begin(), heap.end());
            }
            else if(newNeighborDist == oldNeighborDist){
                prev.at(neighborIdx).push_back(front.nodeIdx);
            }
        }

        if(front.nodeIdx == endPoint){
            lastCalculatedDistance = distance.at(front.nodeIdx);
        }

        // only stop the calculation if all nodes with the shortest path distance lastCalculatedDistance have been been processed
        if(front.dist > lastCalculatedDistance){
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
            currNode = prev.at(currNode).at(0);
            path.push_back(currNode);
        }
    }else{
        std::cout << "no path found" << std::endl;
    }
    std::vector<bool> isNodeInIndependentSet(adjArray.width*adjArray.height, false);
    
}


/**
 * @brief do a backwards depth first search from the end point to the start node to find multiple witness paths
 * 
 * @param isNodeInIndepententSet 
 * @return true 
 * @return false 
 */
bool Dijkstra::checkMultipleShortestPath(std::vector<bool> &isNodeInIndepententSet){
    // check if shortest path exists with no nodes in independent set
    uint64_t limit = 0;
    uint64_t numShortestPaths = 0;
    std::vector<uint64_t> stack;
    stack.push_back(endPoint);

    while(!stack.empty()){
        
        uint64_t currNode = stack.back();
        stack.pop_back();
        limit++;

        if(limit > stepLimit*4){
            break;
        }

        if(currNode == UINT64_MAX){
            continue;
        }

        for(uint64_t nextNode : prev.at(currNode)){
            if(nextNode == UINT64_MAX){
                continue;
            }
            if(nextNode != startPoint && !isNodeInIndepententSet.at(nextNode)){
                stack.push_back(nextNode);
            }
            else if(nextNode == startPoint){
                numShortestPaths++;
                if(numShortestPaths>1){
                    return true;
                }
            }
        }
    }
    // std::cout << "num shortest paths " << numShortestPaths << "\n";
    return false;
};


/**
 * @brief get the number of nodes that were popped from the heap since the last reset
 * 
 * @return uint64_t  number of nodes popped from the heap
 */
uint64_t Dijkstra::getNumNodesPopped(){
    return numNodesPopped;
}


/**
 * @brief disable a node from being used to find a shortest path
 * 
 * @param nodeId the node that is disabled
 */
void Dijkstra::disableNode(uint64_t nodeId){
    disabledNodes.at(nodeId) = true;
}


/**
 * @brief re-enable all nodes
 * 
 */
void Dijkstra::resetDisabledNodes(){
    disabledNodes = std::vector<bool>(adjArray.width*adjArray.height, false);
}


}