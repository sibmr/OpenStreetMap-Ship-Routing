#include <numeric>

#include "shortestPathUtils.cpp"
#include "PathAlgorithm.cpp"

namespace A_star{

/**
 * @brief HeapElement for A_star implementation
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
        std::vector<uint64_t> distance;
        std::vector<HeapElement> heap;
        std::vector<uint64_t> visited;
        std::vector<uint64_t> prev;
        AdjacencyArray &adjArray;
        uint64_t startPoint, endPoint, lastCalculatedDistance;
        uint64_t numNodesPopped;
        
};


/**
 * @brief Construct a new Second A_star:: Second A_star object
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
}

/**
 * @brief reset datastructures to prepare for next call to calculateDist
 */
void A_star::reset(){
    std::vector<uint64_t> initDist(adjArray.width*adjArray.height, UINT64_MAX);
    distance = std::move(initDist);
    // no need to reset prev
    heap.clear();
    // initially no nodes popped
    numNodesPopped = 0;
}

uint64_t A_star::getHeuristic(AdjacencyArray &adjArray, uint64_t firstNodeIdx, uint64_t secondNodeIdx){
    return nodeDistance(adjArray, firstNodeIdx, secondNodeIdx);
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
uint64_t A_star::calculateDist(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;

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

            // update node distance if it improves
            if(newNeighborDist<oldNeighborDist){
                // do not update distance array: only update for distances that are final
                // distance.at(neighborIdx) = newNeighborDist;
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
 * @brief retrieve path calculated by the last call to calculateDist
 * 
 * @param path 
 */
void A_star::getPath(std::vector<uint64_t> &path){
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
        //path.push_back(startPoint);
        //path.push_back(endPoint);
    }
    
}

uint64_t A_star::getNumNodesPopped(){
    return numNodesPopped;
}


class A_star_rectangular : public A_star{
    
    public:
        A_star_rectangular(AdjacencyArray &array) : A_star(array) {}
        uint64_t getHeuristic(AdjacencyArray &adjArray, uint64_t firstNodeIdx, uint64_t secondNodeIdx);
    
};


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