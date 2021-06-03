
# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <math.h>
# include <map>
# include <queue>
# include <algorithm>

struct AdjacencyArray {
    double longLow, latLow, longHigh, latHigh;
    uint64_t width, height;
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> edges;
    std::vector<bool> nodes;

    /**
     * @brief Construct a new Adjacency Array object
     * 
     * @param path path to input .graph file
     */
    AdjacencyArray(std::string path) : offsets(), edges(), nodes(){

        /*
        * input file
        * longLow       - double
        * latLow        - double
        * longHigh      - double
        * latHigh       - double
        * width         - uint64_t
        * height        - uint64_t
        * offset_size   - uint64_t
        * offsets       - uint64_t
        * edges_size    - uint64_t
        * edges         - uint64_t
        * nodes_size    - uint64_t
        * nodes          - bool
        * */

        std::ifstream adjacency_input_file;

        adjacency_input_file.open(path, std::ios::in);

        // write globe size
        adjacency_input_file.read(reinterpret_cast<char *>(&longLow),     sizeof(longLow));
        adjacency_input_file.read(reinterpret_cast<char *>(&latLow),      sizeof(latLow));
        adjacency_input_file.read(reinterpret_cast<char *>(&longHigh),    sizeof(longHigh));
        adjacency_input_file.read(reinterpret_cast<char *>(&latHigh),     sizeof(latHigh));
        

        // save number of nodes per direction (width, height)
        adjacency_input_file.read(reinterpret_cast<char *>(&width),     sizeof(width));
        adjacency_input_file.read(reinterpret_cast<char *>(&height),     sizeof(height));

        // save offset vectro size
        uint64_t offset_size;
        adjacency_input_file.read(reinterpret_cast<char *>(&offset_size),     sizeof(offset_size));
        offsets.resize(offset_size);


        for(uint64_t i = 0; i < offset_size; i++){
            adjacency_input_file.read(reinterpret_cast<char *>(&offsets.at(i)),     sizeof(offsets.at(i)));
        }


        uint64_t edges_size;
        adjacency_input_file.read(reinterpret_cast<char *>(&edges_size),     sizeof(edges_size));
        edges.resize(edges_size);


        for(uint64_t i = 0; i < edges_size; i++){
            adjacency_input_file.read(reinterpret_cast<char *>(&edges.at(i)), sizeof(edges.at(i)));
        }

        uint64_t nodes_size = nodes.size();
        adjacency_input_file.read(reinterpret_cast<char *>(&nodes_size), sizeof(nodes_size));
        nodes.resize(nodes_size);



        std::vector<bool> data(
            (std::istreambuf_iterator<char>(adjacency_input_file)), 
            std::istreambuf_iterator<char>());
        std::copy(
            data.begin(),
            data.end(),
                nodes.begin());

        adjacency_input_file.close();
    }
};

/**
 * @brief HeapElement for SecondDijkstra implementation
 */
struct HeapElement {
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
 * @brief Distance between two points on the globe
 * Implementation adapted to cpp from: https://www.movable-type.co.uk/scripts/latlong.html
 * 
 * @param p1lng     Input:   Two LongLat Points in Degrees
 * @param p1lat 
 * @param p2lng     Output:  Distance between the Points in Meters
 * @param p2lat 
 * @return double 
 */
double latLongDistance(double p1lng, double p1lat, double p2lng, double p2lat){
    
    // Haversine formula:
    const double R = 6.371e6;
    const double phi1 = p1lat           * M_PI/180.0;
    const double phi2 = p2lat           * M_PI/180.0;
    const double dphi = (p2lat-p1lat)   * M_PI/180.0;
    const double dlam = (p2lng-p1lng)   * M_PI/180.0; 

    const double a = sin(dphi/2) * sin(dphi/2) + cos(phi1) * cos(phi2) * sin(dlam/2) * sin(dlam/2);
    const double c = 2 * atan2(sqrt(a), sqrt(1-a));

    const double d = R * c;

    return d;
}

// return (long, lat) from node id

/**
 * @brief from a node id taken from the adjacency array, get the longitude and latitude of the node 
 * 
 * @param result    output: longitude and latitude of the node
 * @param array     adjacency array
 * @param id        id of node in adjacency array
 */
void nodeIdToLongLat(std::array<double,2> &result, AdjacencyArray &array, uint64_t id){
    const double longStep = (array.longHigh-array.longLow)/(array.width+1);
    const double latStep = (array.latHigh-array.latLow)/(array.height+1);

    result[0] =  array.longLow + std::floor(id/array.height) * longStep + longStep/2;
    result[1] =  array.latLow + (id%array.height) * latStep + latStep/2;
    return;

}

/**
 * @brief given (longitude, latitude), get the id of the closest node to that point
 * 
 * @param array         adjacency array
 * @param n1long        longitude position
 * @param n1lat         latitude position
 * @return uint64_t     node id
 */
uint64_t longLatToNodeId(AdjacencyArray &array, double n1long, double n1lat){
    while(n1long < -180){
        n1long += (array.longHigh - array.longLow);
    }
    uint64_t fastIndex =  uint64_t (std::round(std::fmod((n1lat  - array.latLow) / (array.latHigh - array.latLow), 1.0) * array.height));
    uint64_t slowIndex =  uint64_t (std::round(std::fmod((n1long - array.longLow) / (array.longHigh - array.longLow), 1.0) * array.width));
    
    return fastIndex + slowIndex * array.height;
}

/**
 * @brief calcuate the great-circle distance of two nodes on the globe from their ids
 * 
 * @param array         adjacency array
 * @param node1         first node id
 * @param node2         second node id
 * @return uint64_t     distance rounded to uint64_t in meters
 */
uint64_t nodeDistance(AdjacencyArray &array, uint64_t node1, uint64_t node2){

    std::array<double,2> longLat_1;
    std::array<double,2> longLat_2;
    nodeIdToLongLat(longLat_1, array, node1);
    nodeIdToLongLat(longLat_2, array, node2);

    uint64_t result = uint64_t(round(latLongDistance(longLat_1[0], longLat_1[1], longLat_2[0], longLat_2[1])));

    return result;
}

bool isNodeOnLand(AdjacencyArray &adjArray, uint64_t node){
    return adjArray.nodes.at(node);
}

void testLatLongDistance(){
    std::cout << "Test distance:\n";
    // should result in approx 5929000 and 572500 metres
    std::cout << latLongDistance(40,40,20,-10) << " " << latLongDistance(20, -80, 30, -85) << std::endl;
}

/**
 * @brief convert a path of node ids to a path of positions (long1, lat1, long2, lat2, ...)
 * 
 * @param posPath   path of (long1, lat1, long2, lat2, ...)
 * @param idPath    path of ids
 * @param array     adjacency array
 */
void generatePositionPath(std::vector<double> &posPath, std::vector<uint64_t> &idPath,  AdjacencyArray &array){
    std::array<double,2> tmp_longLat;
    for (std::vector<uint64_t>::iterator it = idPath.begin(); it != idPath.end(); ++it) {
        nodeIdToLongLat(tmp_longLat,array,*it);
        posPath.push_back(tmp_longLat[0]);
        posPath.push_back(tmp_longLat[1]);
    }
}

/**
 * @brief from a position path and distance of the path, generate a json response string
 * 
 * @param posPath   input:  path of positions
 * @param response  output: json response string
 * @param distance  input:  distance of the path
 */
void generateReponse(std::vector<double> &posPath, std::string& response, uint64_t distance){
    // response has to be in [[lat,long],...] format, so (long, lat) is swapped
    response += "{\"path\":[";
    if(distance < UINT64_MAX){
        response += "[" + std::to_string(posPath.at(1)) + "," + std::to_string(posPath.at(0)) + "]";
        for(int i=1; i<posPath.size()/2; ++i){
            response += ",[" + std::to_string(posPath.at(2*i+1)) + "," + std::to_string(posPath.at(2*i)) + "]";
    
        }
    }
    response += "],\n";
    response += "\"dist\": ";
    response += std::to_string(distance); 
    response += ",\n";
    response += "\"route\": ";
    response += std::to_string(distance < UINT64_MAX); 
    response += "}";
}


/**
 * @brief virtual interface for shortest path algorithms
 * 
 */
class PathAlgorithm {
    public:
        virtual void getPath(std::vector<uint64_t> &path) = 0;
        virtual uint64_t getDist() = 0;
        virtual uint64_t calculateDist(uint64_t startPoint_, uint64_t endPoint_) = 0;
        virtual void reset() = 0;
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

    heap.push_back(HeapElement{startPoint, UINT64_MAX, 0});
    heap.push_back(HeapElement{endPoint, UINT64_MAX, UINT64_MAX});

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
        }

        distance.at(front.nodeIdx) = front.dist;
        prev.at(front.nodeIdx) = front.prev;
        visited.push_back(front.nodeIdx);

        for(uint64_t currEdgeId = adjArray.offsets.at(front.nodeIdx); currEdgeId < adjArray.offsets.at(front.nodeIdx+1); currEdgeId++){
            uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

            // absolute difference of unsigned int
            uint64_t idxDiff = neighborIdx<front.nodeIdx ? front.nodeIdx-neighborIdx : neighborIdx-front.nodeIdx;

            // choose length of edge from precalculated lengths
            uint64_t edgeDist = (idxDiff < 2) ? constLngDist : constLatDist.at(neighborIdx%adjArray.height);

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

void test() {
    AdjacencyArray adjArray("data/worldGrid_1415_707.save");
    FirstDijkstra fd(adjArray);
    PathAlgorithm &pa = fd;

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