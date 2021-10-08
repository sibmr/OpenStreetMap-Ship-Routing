#include "shortestPathUtils.cpp"
#include "PathAlgorithm.cpp"

/**
 * @brief HeapElement for SecondDijkstra implementation
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
 * @brief Our first implementation of the dijkstra algorithm
 * 
 */
class FirstDijkstra: public PathAlgorithm{
    public:
        FirstDijkstra(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        void reset();
    private:
        void generatePath(uint64_t startPoint, uint64_t endPoint, std::vector<uint64_t> &path);
        void fillMaps(uint64_t startPoint, uint64_t endPoint);
        
        AdjacencyArray &array;
        
        std::map<uint64_t,uint64_t> distance;
        std::map<uint64_t,uint64_t> previous;

        uint64_t startPoint, endPoint, lastCalculatedDistance;
};

/**
 * @brief Our second more efficient implementation of the dijkstra algorithm
 * Similar implementation to https://github.com/Lesstat/dijkstra-performance-study/
 */
class SecondDijkstra: public PathAlgorithm{
    public:
        SecondDijkstra(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        void reset();
        uint64_t getNumNodesPopped();
    private:
        uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
        std::vector<uint64_t> distance;
        std::vector<HeapElement> heap;
        std::vector<uint64_t> visited;
        std::vector<uint64_t> prev;
        std::vector<bool> contractedNodes;
        AdjacencyArray &adjArray;
        uint64_t constLngDist;
        std::vector<uint64_t> constLatDist;
        uint64_t startPoint, endPoint, lastCalculatedDistance;
        uint64_t numNodesPopped;
};


FirstDijkstra::FirstDijkstra(AdjacencyArray &array) : array(array){}

/**
 * @brief retrieve a path of nodes from the last call to calculateDist
 * this only works if calculateNodes was called before this
 * 
 * @param path  output: path of nodes
 */
void FirstDijkstra::getPath(std::vector<uint64_t> &path){
    if(distance.at(endPoint) < UINT64_MAX){
        // build up path
        uint64_t currNode = endPoint;
        while(currNode != startPoint){
            currNode = previous.at(currNode);
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
 * @brief reset and prepare the datastructures for the next call to calculateDist
 */
void FirstDijkstra::reset(){
    distance.clear();
    previous.clear();
    // reset map
    for (int i = 0; i < array.nodes.size(); i++){
        distance.insert({i, UINT64_MAX});
        previous.insert({i, UINT64_MAX});
    }
}

/**
 * @brief returns the last distance calculated by calculateDist
 * 
 * @return uint64_t
 */
uint64_t FirstDijkstra::getDist(){
    return lastCalculatedDistance;
}

/**
 * @brief this method computes the dijkstra algorithm to find the shortest path from startPoint_ to endPoint_
 * 
 * @param startPoint_   node id of start node
 * @param endPoint_     node id of goal node
 * @return uint64_t 
 */
uint64_t FirstDijkstra::calculateDist(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;

    distance.at(startPoint) = 0;
    typedef std::pair<uint64_t, uint64_t> pqPair;


        


    std::priority_queue<pqPair, std::vector<pqPair>, std::greater<pqPair>> pq;
    pq.push(std::make_pair(0, startPoint));
    std::pair<uint64_t, uint64_t> currTop; // current top element

    uint64_t currSourceNode;
    uint64_t currTargetNode;
    uint64_t currDistance;

    uint64_t currTargetDistance = UINT64_MAX;
    while (true){
        if(pq.empty()){
            // no node can be expanded anymore
            lastCalculatedDistance = distance.at(endPoint);
            return lastCalculatedDistance;
        }
        currTop = pq.top(); // first element is dist, second in node id
        pq.pop();
        //std::cout << pq.size() << std::endl;
        currSourceNode = currTop.second;
        //std::cout << pq.size() << std::endl;

        for(uint64_t currEdgeId = array.offsets.at(currSourceNode); currEdgeId < array.offsets.at(currSourceNode+1); currEdgeId++){
            currTargetNode = array.edges.at(currEdgeId);
            currDistance = distance.at(currSourceNode) + nodeDistance(array, currSourceNode, currTargetNode);
            
            if(distance.at(currTargetNode) > currDistance){
                if(currTargetNode == endPoint){
                    // save current best distance to endPoint
                    currTargetDistance = currDistance;
                }
                distance.at(currTargetNode) = currDistance;
                previous.at(currTargetNode) = currSourceNode;

                // do not expand node if it has higher distance than the current target distance
                if(currDistance <= currTargetDistance){
                    pq.push(std::make_pair(currDistance, currTargetNode));
                }
            }
        }
    }
}

/**
 * @brief Construct a new Second Dijkstra:: Second Dijkstra object
 * 
 * Initialize datastructures
 * Initialize distances between nodes
 * 
 * @param array 
 */
SecondDijkstra::SecondDijkstra(AdjacencyArray &array) : adjArray(array), prev(array.width*array.height, UINT64_MAX){
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
void SecondDijkstra::reset(){
    std::vector<uint64_t> initDist(adjArray.width*adjArray.height, UINT64_MAX);
    distance = std::move(initDist);
    // no need to reset prev
    heap.clear();
    // initially no nodes are popped
    numNodesPopped = 0;
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
uint64_t SecondDijkstra::calculateDist(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;


    //early stop with same start value
    if(distance.at(startPoint) == 0 && distance.at(endPoint) != UINT64_MAX){
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
        if(front.dist >= distance.at(front.nodeIdx)){
            continue;
        }else{
            numNodesPopped++;
        }

        distance.at(front.nodeIdx) = front.dist;
        prev.at(front.nodeIdx) = front.prev;
        visited.push_back(front.nodeIdx);

        for(uint64_t currEdgeId = adjArray.offsets.at(front.nodeIdx); currEdgeId < adjArray.offsets.at(front.nodeIdx+1); currEdgeId++){
            uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

            // choose length of edge from precalculated lengths
            uint64_t edgeDist = adjArray.distances.at(currEdgeId);

            uint64_t newNeighborDist = front.dist + edgeDist;
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
            return lastCalculatedDistance;
        }

    }

}

/**
 * @brief returns distance of last call to calculateDist
 * 
 * @return uint64_t 
 */
uint64_t SecondDijkstra::getDist(){
    return lastCalculatedDistance;
}

/**
 * @brief retrieve path calculated by the last call to calculateDist
 * 
 * @param path 
 */
void SecondDijkstra::getPath(std::vector<uint64_t> &path){
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

uint64_t SecondDijkstra::getNumNodesPopped(){
    return numNodesPopped;
}



void test() {
    AdjacencyArray adjArray("data/worldGrid_1415_707.save");
    SecondDijkstra sd(adjArray);
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