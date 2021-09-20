#include "shortestPathUtils.cpp"
#include "PathAlgorithm.cpp"

/**
 * @brief HeapElement for DijkstraImpl implementation
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
 * @brief Our more efficient implementation of the dijkstra algorithm
 * Similar implementation to https://github.com/Lesstat/dijkstra-performance-study/
 */
class DijkstraImpl: public PathAlgorithm{
    public:
        DijkstraImpl(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        uint64_t calculateDistSavedEdges(uint64_t startPoint, uint64_t endPoint);
        void reset();
    private:
        uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
        std::vector<uint64_t> distance;
        std::vector<HeapElement> heap;
        std::vector<uint64_t> visited;
        std::vector<uint64_t> prev;
        AdjacencyArray &adjArray;
        uint64_t constLngDist;
        std::vector<uint64_t> constLatDist;
        uint64_t startPoint, endPoint, lastCalculatedDistance;
};

/**
 * @brief Our more efficient implementation of the dijkstra algorithm
 * Similar implementation to https://github.com/Lesstat/dijkstra-performance-study/
 */
class DijkstraSavedEdges: public PathAlgorithm{
    public:
        DijkstraSavedEdges(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        void reset();
    private:
        uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
        std::vector<uint64_t> distance;
        std::vector<HeapElement> heap;
        std::vector<uint64_t> visited;
        std::vector<uint64_t> prev;
        AdjacencyArray &adjArray;
        uint64_t constLngDist;
        std::vector<uint64_t> constLatDist;
        uint64_t startPoint, endPoint, lastCalculatedDistance;
};

/**
 * @brief Our more efficient implementation of the dijkstra algorithm
 * Similar implementation to https://github.com/Lesstat/dijkstra-performance-study/
 */
class DijkstraBiDirect: public PathAlgorithm{
    public:
        DijkstraBiDirect(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        void reset();
    private:
        uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
        std::vector<uint64_t> distanceSource;
        std::vector<uint64_t> distanceTarget;
        std::vector<HeapElement> heapSource;
        std::vector<HeapElement> heapTarget;
        std::vector<uint64_t> visitedSource;
        std::vector<uint64_t> visitedTarget;
        std::vector<uint64_t> prevSource;
        std::vector<uint64_t> prevTarget;
        AdjacencyArray &adjArray;
        uint64_t constLngDist;
        std::vector<uint64_t> constLatDist;
        uint64_t startPoint, endPoint, lastCalculatedDistance;
};


/**
 * @brief Construct a new Dijkstra:: DijkstraImpl object
 * 
 * Initialize datastructures
 * Initialize distances between nodes
 * 
 * @param array 
 */
DijkstraImpl::DijkstraImpl(AdjacencyArray &array) : adjArray(array), prev(array.width*array.height, UINT64_MAX){
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
void DijkstraImpl::reset(){
    std::vector<uint64_t> initDist(adjArray.width*adjArray.height, UINT64_MAX);
    distance = std::move(initDist);
    // no need to reset prev
    heap.clear();
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
uint64_t DijkstraImpl::calculateDist(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;

    // early stop with same start value
    if(std::find (visited.begin(), visited.end(), endPoint) != visited.end() && distance.at(startPoint) == 0){
        return distance.at(endPoint);
    }else{
        reset();
    }

    heap.push_back(HeapElement{startPoint, UINT64_MAX, 0});

    std::make_heap(heap.begin(), heap.end());

    HeapElement front;
    while(true){
        if(heap.empty()){
            return UINT64_MAX;
        }

        std::pop_heap(heap.begin(), heap.end());
        front = heap.back();
        heap.pop_back();

        // avoid duplicate nodes
        if(front.heuristic_dist >= distance.at(front.nodeIdx)){
            continue;
        }

        distance.at(front.nodeIdx) = front.heuristic_dist;
        prev.at(front.nodeIdx) = front.prev;
        visited.push_back(front.nodeIdx);

        for(uint64_t currEdgeId = adjArray.offsets.at(front.nodeIdx); currEdgeId < adjArray.offsets.at(front.nodeIdx+1); currEdgeId++){
            uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

            // absolute difference of unsigned int
            uint64_t idxDiff = neighborIdx<front.nodeIdx ? front.nodeIdx-neighborIdx : neighborIdx-front.nodeIdx;

            // choose length of edge from precalculated lengths
            // WRONG CALCULATION OF DISTANCE BETWEEN TWO NODES
            uint64_t edgeDist = (idxDiff < 2) ? constLngDist : constLatDist.at(neighborIdx%adjArray.height);
            //uint64_t edgeDist = nodeDistance(adjArray, front.nodeIdx, neighborIdx);

            uint64_t newNeighborDist = front.heuristic_dist + edgeDist;
            uint64_t oldNeighborDist = distance.at(neighborIdx);

            if(newNeighborDist<oldNeighborDist){
                // do not update distance array: only update for distances that are final
                // distance.at(neighborIdx) = newNeighborDist;
                heap.push_back(HeapElement{neighborIdx, front.nodeIdx, newNeighborDist});
                std::push_heap(heap.begin(), heap.end());
            }
        }

        if(front.nodeIdx == endPoint){
            lastCalculatedDistance = distance.at(front.nodeIdx);
            std::cout << "Traversed nodes " << visited.size() << std::endl;
            return lastCalculatedDistance;
        }

    }

}






/**
 * @brief returns distance of last call to calculateDist
 * 
 * @return uint64_t 
 */
uint64_t DijkstraImpl::getDist(){
    return lastCalculatedDistance;
}

/**
 * @brief retrieve path calculated by the last call to calculateDist
 * 
 * @param path 
 */
void DijkstraImpl::getPath(std::vector<uint64_t> &path){
    if(distance.at(endPoint) < UINT64_MAX){
        // build up path
        uint64_t currNode = endPoint;
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


/**
 * @brief Construct a new Dijkstra:: DijkstraImpl object use distances from Memory
 * 
 * Initialize datastructures
 * Initialize distances between nodes
 * 
 * @param array 
 */
DijkstraSavedEdges::DijkstraSavedEdges(AdjacencyArray &array) : adjArray(array), prev(array.width*array.height, UINT64_MAX){
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
void DijkstraSavedEdges::reset(){
    std::vector<uint64_t> initDist(adjArray.width*adjArray.height, UINT64_MAX);
    distance = std::move(initDist);
    // no need to reset prev
    heap.clear();
}



/**
 * @brief efficient dijkstra shortest-path implementation
 * 
 * Uses binary Min(Max)-heap for greedy node visitation strategy
 * 
 * Saved Distances
 * 
 * @param startPoint 
 * @param endPoint 
 * @return uint64_t 
 */
uint64_t DijkstraSavedEdges::calculateDist(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;

    // early stop with same start value
    if(std::find (visited.begin(), visited.end(), endPoint) != visited.end() && distance.at(startPoint) == 0){
        return distance.at(endPoint);
    }else{
        reset();
    }

    heap.push_back(HeapElement{startPoint, UINT64_MAX, 0});

    std::make_heap(heap.begin(), heap.end());

    HeapElement front;
    while(true){
        if(heap.empty()){
            return UINT64_MAX;
        }

        std::pop_heap(heap.begin(), heap.end());
        front = heap.back();
        heap.pop_back();

        // avoid duplicate nodes
        if(front.heuristic_dist >= distance.at(front.nodeIdx)){
            continue;
        }

        distance.at(front.nodeIdx) = front.heuristic_dist;
        prev.at(front.nodeIdx) = front.prev;
        visited.push_back(front.nodeIdx);

        for(uint64_t currEdgeId = adjArray.offsets.at(front.nodeIdx); currEdgeId < adjArray.offsets.at(front.nodeIdx+1); currEdgeId++){
            uint64_t neighborIdx = adjArray.edges.at(currEdgeId);


            // choose length of edge from precalculated lengths
            uint64_t edgeDist = adjArray.distances.at(currEdgeId);

            uint64_t newNeighborDist = front.heuristic_dist + edgeDist;
            uint64_t oldNeighborDist = distance.at(neighborIdx);

            if(newNeighborDist<oldNeighborDist){
                // do not update distance array: only update for distances that are final
                // distance.at(neighborIdx) = newNeighborDist;
                heap.push_back(HeapElement{neighborIdx, front.nodeIdx, newNeighborDist});
                std::push_heap(heap.begin(), heap.end());
            }
        }

        if(front.nodeIdx == endPoint){
            lastCalculatedDistance = distance.at(front.nodeIdx);
            std::cout << "Traversed nodes " << visited.size() << std::endl;
            return lastCalculatedDistance;
        }

    }

}

/**
 * @brief returns distance of last call to calculateDist
 * 
 * @return uint64_t 
 */
uint64_t DijkstraSavedEdges::getDist(){
    return lastCalculatedDistance;
}

/**
 * @brief retrieve path calculated by the last call to calculateDist
 * 
 * @param path 
 */
void DijkstraSavedEdges::getPath(std::vector<uint64_t> &path){
    if(distance.at(endPoint) < UINT64_MAX){
        // build up path
        uint64_t currNode = endPoint;
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


/**
 * @brief Construct a new Dijkstra:: DijkstraImpl object from both sides
 * 
 * Initialize datastructures
 * Initialize distances between nodes
 * 
 * @param array 
 */
DijkstraBiDirect::DijkstraBiDirect(AdjacencyArray &array) : adjArray(array), prevSource(array.width*array.height, UINT64_MAX), prevTarget(array.width*array.height, UINT64_MAX){
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
void DijkstraBiDirect::reset(){
    std::vector<uint64_t> initDist(adjArray.width*adjArray.height, UINT64_MAX);
    std::vector<uint64_t> initDist2(adjArray.width*adjArray.height, UINT64_MAX);
    distanceSource = std::move(initDist);
    distanceTarget = std::move(initDist2);
    // no need to reset prev
    heapSource.clear();
    heapTarget.clear();
}



/**
 * @brief efficient dijkstra shortest-path implementation
 * 
 * Uses binary Min(Max)-heap for greedy node visitation strategy
 * 
 * Saved Distances
 * 
 * @param startPoint 
 * @param endPoint 
 * @return uint64_t 
 */
uint64_t DijkstraBiDirect::calculateDist(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;
    bool sourceTree = true;

    uint64_t sourceDist = 0;
    uint64_t targetDist = 0;
    uint64_t currDist;
    std::vector<uint64_t> * distancePointer;

    // early stop with same start value
    //if(std::find (visited.begin(), visited.end(), endPoint) != visited.end() && distance.at(startPoint) == 0){
    //    return distance.at(endPoint);
    //}else{
    //    reset();
    //}

    heapSource.push_back(HeapElement{startPoint, UINT64_MAX, 0});
    heapTarget.push_back(HeapElement{endPoint, UINT64_MAX, 0});

    std::make_heap(heapSource.begin(), heapSource.end());
    std::make_heap(heapTarget.begin(), heapTarget.end());

    HeapElement front;
    while(true){
        if(heapSource.empty() && heapTarget.empty()){

            std::cout << "Traversed nodes " << visitedSource.size() << " " << visitedTarget.size() << std::endl;

            //return UINT64_MAX;
            if(sourceDist != 0 && targetDist != 0){
                return (sourceDist + targetDist) / 2;
            }else{
                return UINT64_MAX;
            }
        }

        std::pop_heap(heapSource.begin(), heapSource.end());
        std::pop_heap(heapTarget.begin(), heapTarget.end());

        // POTENTIAL SOURCE OF FAILURE
        if(heapTarget.empty() || (!heapSource.empty() && heapSource.back().heuristic_dist <= heapTarget.back().heuristic_dist)){
            sourceTree = true;
            front = heapSource.back();
            heapSource.pop_back();
        }else{
            sourceTree = false;
            front = heapTarget.back();
            heapTarget.pop_back();
        }


        // original dijkstra
        if(sourceTree){
            // avoid duplicate nodes
            if(front.heuristic_dist >= distanceSource.at(front.nodeIdx)){
                continue;
            }
            distanceSource.at(front.nodeIdx) = front.heuristic_dist;
            prevSource.at(front.nodeIdx) = front.prev;
            visitedSource.push_back(front.nodeIdx);
        
            for(uint64_t currEdgeId = adjArray.offsets.at(front.nodeIdx); currEdgeId < adjArray.offsets.at(front.nodeIdx+1); currEdgeId++){
                uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

                // choose length of edge from precalculated lengths
                uint64_t edgeDist = adjArray.distances.at(currEdgeId);

                uint64_t newNeighborDist = front.heuristic_dist + edgeDist;
                uint64_t oldNeighborDist = distanceSource.at(neighborIdx);

                if(newNeighborDist<oldNeighborDist){
                    // do not update distance array: only update for distances that are final
                    // distance.at(neighborIdx) = newNeighborDist;
                    heapSource.push_back(HeapElement{neighborIdx, front.nodeIdx, newNeighborDist});
                    std::push_heap(heapSource.begin(), heapSource.end());
                }
            }
            if(front.nodeIdx == endPoint){
                lastCalculatedDistance = distanceSource.at(front.nodeIdx);
                sourceDist = lastCalculatedDistance;
                heapSource.clear();
                //return lastCalculatedDistance;
            }
        }else{
            // avoid duplicate nodes
            if(front.heuristic_dist >= distanceTarget.at(front.nodeIdx)){
                continue;
            }
            distanceTarget.at(front.nodeIdx) = front.heuristic_dist;
            prevTarget.at(front.nodeIdx) = front.prev;
            visitedTarget.push_back(front.nodeIdx);

            for(uint64_t currEdgeId = adjArray.offsets.at(front.nodeIdx); currEdgeId < adjArray.offsets.at(front.nodeIdx+1); currEdgeId++){
                uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

                // choose length of edge from precalculated lengths
                uint64_t edgeDist = adjArray.distances.at(currEdgeId);

                uint64_t newNeighborDist = front.heuristic_dist + edgeDist;
                uint64_t oldNeighborDist = distanceTarget.at(neighborIdx);

                if(newNeighborDist<oldNeighborDist){
                    // do not update distance array: only update for distances that are final
                    // distance.at(neighborIdx) = newNeighborDist;
                    heapTarget.push_back(HeapElement{neighborIdx, front.nodeIdx, newNeighborDist});
                    std::push_heap(heapTarget.begin(), heapTarget.end());
                }
            }
            if(front.nodeIdx == startPoint){
                lastCalculatedDistance = distanceTarget.at(front.nodeIdx);
                targetDist = lastCalculatedDistance;
                heapTarget.clear();
                //return lastCalculatedDistance;
            }
        }

    }

}

/**
 * @brief returns distance of last call to calculateDist
 * 
 * @return uint64_t 
 */
uint64_t DijkstraBiDirect::getDist(){
    return lastCalculatedDistance;
}

/**
 * @brief retrieve path calculated by the last call to calculateDist
 * 
 * @param path 
 */
void DijkstraBiDirect::getPath(std::vector<uint64_t> &path){
    return;
    //if(distance.at(endPoint) < UINT64_MAX){
    //    // build up path
    //    uint64_t currNode = endPoint;
    //    while(currNode != startPoint){
    //        currNode = prev.at(currNode);
    //        path.push_back(currNode);
    //    }

    //    // print path
    //    std::cout << "Path:" << std::endl;
    //    for (std::vector<uint64_t>::iterator it = path.begin(); it != path.end(); ++it) {
    //        std::cout << *it << " ";
    //    }
    //    std::cout <<  std::endl;
    //    std::cout << "dist: " << distance.at(endPoint)/1000 << "km" << std::endl;
    //}else{
    //    std::cout << "no path found" << std::endl;
    //    //path.push_back(startPoint);
    //    //path.push_back(endPoint);
    //}
    
}


void test() {
    AdjacencyArray adjArray("data/worldGrid_1415_707.save");
    DijkstraBiDirect sd(adjArray);
    PathAlgorithm &pa = sd;

    std::vector<uint64_t> path;
    pa.reset();
    pa.calculateDist(1001, 16001);
    pa.getPath(path);
    std::vector<double> posPath;
    generatePositionPath(posPath, path, adjArray);


    std::string response;
    generateReponse(posPath, response, pa.getDist());
    std::cout << response << std::endl;

    int counter_one = 0;
    int counter_zero = 0;
    for (int i = 0; i < adjArray.nodes.size(); i++){
        if(adjArray.nodes.at(i) == 0){
            counter_zero++;
        }else
        {
            counter_one++;
        }
    }
    std::cout << counter_zero << " " << counter_one << " " << counter_one + counter_zero << std::endl;
    std::cout << adjArray.edges.size() << std::endl;
    testLatLongDistance();
}