#include <set>

#include "shortestPathUtils.cpp"
#include "PathAlgorithm.cpp"

namespace BidirectionalDijkstra{

/**
 * @brief HeapElement for BidirectionalDijkstra implementation
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
 * @brief Our second more efficient implementation of the dijkstra algorithm
 * Similar implementation to https://github.com/Lesstat/dijkstra-performance-study/
 */
class BidirectionalDijkstra: public PathAlgorithm{
    public:
        BidirectionalDijkstra(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        void reset();
    private:
        uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
        std::vector<uint64_t> forwardDistance;
        std::vector<uint64_t> backwardDistance;
        std::vector<HeapElement> forwardHeap;
        std::vector<HeapElement> backwardHeap;
        std::set<uint64_t> visited;
        std::vector<uint64_t> forwardPrev;
        std::vector<uint64_t> backwardPrev;
        AdjacencyArray &adjArray;
        uint64_t startPoint, endPoint, lastCalculatedDistance, forwardMinMeetingNodeId, backwardMinMeetingNodeId;
};


/**
 * @brief Construct a new Second BidirectionalDijkstra:: Second BidirectionalDijkstra object
 * 
 * Initialize datastructures
 * Initialize distances between nodes
 * 
 * @param array 
 */
BidirectionalDijkstra::BidirectionalDijkstra(AdjacencyArray &array) : 
                adjArray(array), 
                forwardPrev(array.width*array.height, UINT64_MAX),
                backwardPrev(array.width*array.height, UINT64_MAX)
{
    reset();
}

/**
 * @brief reset datastructures to prepare for next call to calculateDist
 */
void BidirectionalDijkstra::reset(){
    lastCalculatedDistance = UINT64_MAX;
    std::vector<uint64_t> initDistF(adjArray.width*adjArray.height, UINT64_MAX);
    forwardDistance = std::move(initDistF);
    std::vector<uint64_t> initDistB(adjArray.width*adjArray.height, UINT64_MAX);
    backwardDistance = std::move(initDistB);
    // no need to reset prev
    forwardHeap.clear();
    backwardHeap.clear();
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
uint64_t BidirectionalDijkstra::calculateDist(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;

    forwardHeap.push_back(HeapElement{startPoint, UINT64_MAX, 0});
    backwardHeap.push_back(HeapElement{endPoint, UINT64_MAX, 0});

    std::make_heap(forwardHeap.begin(), forwardHeap.end());
    std::make_heap(backwardHeap.begin(), backwardHeap.end());

    HeapElement forwardFront;
    HeapElement backwardFront;

    uint64_t minDistance = UINT64_MAX;

    while(true){
        
        //// forward step ////

        // get (non-duplicate) heap front
        // avoid duplicate nodes (nodes that were already visited, indicated by higher distance)
        do{
            // no path found
            if(forwardHeap.empty()){
                return UINT64_MAX;
            }
            std::pop_heap(forwardHeap.begin(), forwardHeap.end());
            forwardFront = forwardHeap.back();
            forwardHeap.pop_back();
        }while (forwardFront.heuristic_dist >= forwardDistance.at(forwardFront.nodeIdx));
        do{
            // no path found
            if(backwardHeap.empty()){
                return UINT64_MAX;
            }
            std::pop_heap(backwardHeap.begin(), backwardHeap.end());
            backwardFront = backwardHeap.back();
            backwardHeap.pop_back();
        }while (backwardFront.heuristic_dist >= backwardDistance.at(backwardFront.nodeIdx));

        // update distance and previous node of current forward node
        forwardDistance.at(forwardFront.nodeIdx) = forwardFront.heuristic_dist;
        forwardPrev.at(forwardFront.nodeIdx) = forwardFront.prev;
        // update distance and previous node of current backward node
        backwardDistance.at(backwardFront.nodeIdx) = backwardFront.heuristic_dist;
        backwardPrev.at(backwardFront.nodeIdx) = backwardFront.prev;

        // termination criterion for bidirectional dijkstra
        if(forwardFront.heuristic_dist + backwardFront.heuristic_dist >= minDistance){
            lastCalculatedDistance = minDistance;
            return lastCalculatedDistance;
        }

        // iterate over edges of current forward node
        for(uint64_t currEdgeId = adjArray.offsets.at(forwardFront.nodeIdx); currEdgeId < adjArray.offsets.at(forwardFront.nodeIdx+1); currEdgeId++){
            
            // get id of adjacent node for current edge (neighboring node)
            uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

            // choose length of edge from precalculated lengths
            uint64_t edgeDist = adjArray.distances.at(currEdgeId);

            // calculate new distances
            uint64_t newNeighborDist = forwardFront.heuristic_dist + edgeDist;
            uint64_t oldNeighborDist = forwardDistance.at(neighborIdx);

            // update node distance if it improves (in this case: push new heap node with better distance)
            if(newNeighborDist<oldNeighborDist){
                forwardHeap.push_back(HeapElement{neighborIdx, forwardFront.nodeIdx, newNeighborDist});
                std::push_heap(forwardHeap.begin(), forwardHeap.end());
            }

            // if neighbor has been discovered in "backward" dijkstra then check if min distance improves
            uint64_t neighborGoalDistance = backwardDistance.at(neighborIdx);
            if(neighborGoalDistance < UINT64_MAX){
                uint64_t joinDistance = forwardFront.heuristic_dist + edgeDist + neighborGoalDistance;
                if(joinDistance < minDistance){
                    minDistance = joinDistance;
                    forwardMinMeetingNodeId = forwardFront.nodeIdx;
                    backwardMinMeetingNodeId = neighborIdx;
                }
            }
        }

        // iterate over edges of current backward node
        for(uint64_t currEdgeId = adjArray.offsets.at(backwardFront.nodeIdx); currEdgeId < adjArray.offsets.at(backwardFront.nodeIdx+1); currEdgeId++){
            
            // get id of adjacent node for current edge (neighboring node)
            uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

            // choose length of edge from precalculated lengths
            uint64_t edgeDist = adjArray.distances.at(currEdgeId);

            // calculate new distances
            uint64_t newNeighborDist = backwardFront.heuristic_dist + edgeDist;
            uint64_t oldNeighborDist = backwardDistance.at(neighborIdx);

            // update node distance if it improves (in this case: push new heap node with better distance)
            if(newNeighborDist<oldNeighborDist){
                backwardHeap.push_back(HeapElement{neighborIdx, backwardFront.nodeIdx, newNeighborDist});
                std::push_heap(backwardHeap.begin(), backwardHeap.end());
            }

            // if neighbor has been discovered in "backward" dijkstra then check if min distance improves
            uint64_t neighborStartDistance = forwardDistance.at(neighborIdx);
            if(neighborStartDistance < UINT64_MAX){
                uint64_t joinDistance = backwardFront.heuristic_dist + edgeDist + neighborStartDistance;
                if(joinDistance < minDistance){
                    minDistance = joinDistance;
                    forwardMinMeetingNodeId = neighborIdx;
                    backwardMinMeetingNodeId = backwardFront.nodeIdx;
                }
            }
        }


        // // current node is goal node OR found node visited by other dijkstra
        // if(front.nodeIdx == endPoint || visited.find(front.nodeIdx) != visited.end()){
        //     lastCalculatedDistance = distance.at(front.nodeIdx);
        //     visited.insert(front.nodeIdx);
        //     return lastCalculatedDistance;
        // }else{
        //     visited.insert(front.nodeIdx);
        // }

    }

}

/**
 * @brief returns distance of last call to calculateDist
 * 
 * @return uint64_t 
 */
uint64_t BidirectionalDijkstra::getDist(){
    return lastCalculatedDistance;
}

/**
 * @brief retrieve path calculated by the last call to calculateDist
 * 
 * @param path 
 */
void BidirectionalDijkstra::getPath(std::vector<uint64_t> &path){
    std::vector<uint64_t> forwardPath;
    std::vector<uint64_t> backwardPath;
    if(lastCalculatedDistance < UINT64_MAX){
        // build up path
        uint64_t currNode = forwardMinMeetingNodeId;
        while(currNode != startPoint){
            std::cout << "p1: " << startPoint << " " << endPoint << "\n";
            std::cout << currNode << "\n";
            currNode = forwardPrev.at(currNode);
            forwardPath.push_back(currNode);
        }
        currNode = backwardMinMeetingNodeId;
        while(currNode != endPoint){
            std::cout << "p2: " << startPoint << " " << endPoint << "\n";
            std::cout << currNode << "\n";
            currNode = backwardPrev.at(currNode);
            backwardPath.push_back(currNode);
        }
        for(int i = forwardPath.size()-1; i>=0; --i){
            path.push_back(forwardPath.at(i));
        }
        for(int i = 0; i<backwardPath.size(); ++i){
            path.push_back(backwardPath.at(i));
        }

        // print path
        std::cout << "Path:" << std::endl;
        for (std::vector<uint64_t>::iterator it = path.begin(); it != path.end(); ++it) {
            std::cout << *it << " ";
        }
        std::cout <<  std::endl;
        std::cout << "dist: " << lastCalculatedDistance/1000 << "km" << std::endl;
    }else{
        std::cout << "no path found" << std::endl;
        //path.push_back(startPoint);
        //path.push_back(endPoint);
    }
    
}


}