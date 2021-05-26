
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

    AdjacencyArray(std::string path) : offsets(), edges(), nodes(){
        std::ifstream adjacency_input_file;

        size_t lastindex = path.find_last_of(".");
        std::string adjacency_file_name = path.substr(0, lastindex);
        adjacency_file_name += "_adjacencyarray.save";

        adjacency_input_file.open(adjacency_file_name, std::ios::in);

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

struct HeapElement {
    uint64_t nodeIdx, prev, dist;
    
    // true if own distance is bigger than others distance
    bool operator<(const HeapElement &a){
        return dist > a.dist;
    }
};

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
void loadAdjacencyArray(AdjacencyArray &array, std::string path){
    std::ifstream adjacency_input_file;

    size_t lastindex = path.find_last_of(".");
    std::string adjacency_file_name = path.substr(0, lastindex);
    adjacency_file_name += "_adjacencyarray.save";

    adjacency_input_file.open(adjacency_file_name, std::ios::in);

    // write globe size
    adjacency_input_file.read(reinterpret_cast<char *>(&array.longLow),     sizeof(array.longLow));
    adjacency_input_file.read(reinterpret_cast<char *>(&array.latLow),      sizeof(array.latLow));
    adjacency_input_file.read(reinterpret_cast<char *>(&array.longHigh),    sizeof(array.longHigh));
    adjacency_input_file.read(reinterpret_cast<char *>(&array.latHigh),     sizeof(array.latHigh));
    

    // save number of nodes per direction (width, height)
    adjacency_input_file.read(reinterpret_cast<char *>(&array.width),     sizeof(array.width));
    adjacency_input_file.read(reinterpret_cast<char *>(&array.height),     sizeof(array.height));

    // save offset vectro size
    uint64_t offset_size;
    adjacency_input_file.read(reinterpret_cast<char *>(&offset_size),     sizeof(offset_size));
    array.offsets.resize(offset_size);


    for(uint64_t i = 0; i < offset_size; i++){
        adjacency_input_file.read(reinterpret_cast<char *>(&array.offsets.at(i)),     sizeof(array.offsets.at(i)));
    }


    uint64_t edges_size;
    adjacency_input_file.read(reinterpret_cast<char *>(&edges_size),     sizeof(edges_size));
    array.edges.resize(edges_size);


    for(uint64_t i = 0; i < edges_size; i++){
        adjacency_input_file.read(reinterpret_cast<char *>(&array.edges.at(i)), sizeof(array.edges.at(i)));
    }

    uint64_t nodes_size = array.nodes.size();
    adjacency_input_file.read(reinterpret_cast<char *>(&nodes_size), sizeof(nodes_size));
    array.nodes.resize(nodes_size);



    std::vector<bool> data(
        (std::istreambuf_iterator<char>(adjacency_input_file)), 
        std::istreambuf_iterator<char>());
    std::copy(
        data.begin(),
        data.end(),
            array.nodes.begin());

    adjacency_input_file.close();
}

double latLongDistance(double p1lng, double p1lat, double p2lng, double p2lat){
    /**
     * Input:   Two LongLat Points in Degrees
     * Output:  Distance between the Points in Metres
     **/
    
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
void nodeIdToLongLat(std::array<double,2> &result, AdjacencyArray &array, uint64_t id){
    const double longStep = (array.longHigh-array.longLow)/(array.width+1);
    const double latStep = (array.latHigh-array.latLow)/(array.height+1);

    result[0] =  array.longLow + std::floor(id/array.height) * longStep + longStep/2;
    result[1] =  array.latLow + (id%array.height) * latStep + latStep/2;
    return;

}

uint64_t longLatToNodeId(AdjacencyArray &array, double n1long, double n1lat){
    while(n1long < -180){
        n1long += (array.longHigh - array.longLow);
    }
    uint64_t fastIndex =  uint64_t (std::round(std::fmod((n1lat  - array.latLow) / (array.latHigh - array.latLow), 1.0) * array.height));
    uint64_t slowIndex =  uint64_t (std::round(std::fmod((n1long - array.longLow) / (array.longHigh - array.longLow), 1.0) * array.width));
    //std::cout << fastIndex + slowIndex * array.height << std::endl;
    return fastIndex + slowIndex * array.height;
}

uint64_t nodeDistance(AdjacencyArray &array, uint64_t node1, uint64_t node2){

    std::array<double,2> longLat_1;
    std::array<double,2> longLat_2;
    nodeIdToLongLat(longLat_1, array, node1);
    nodeIdToLongLat(longLat_2, array, node2);

    uint64_t result = uint64_t(latLongDistance(longLat_1[0], longLat_1[1], longLat_2[0], longLat_2[1]));

    //std::cout << "pos: "<< n1_long << " " << n1_lat << "\t" << "pos: "<< n2_long << " " << n2_lat << "\tdist: " << result <<std::endl;
    return result; //result;
}


void testLatLongDistance(){
    std::cout << "Test distance:\n";
    // should result in approx 5929000 and 572500 metres
    std::cout << latLongDistance(40,40,20,-10) << " " << latLongDistance(20, -80, 30, -85) << std::endl;
}

void generateNodePath(std::vector<double> &nodePath, std::vector<uint64_t> &path,  AdjacencyArray &array){
    std::array<double,2> tmp_longLat;
    for (std::vector<uint64_t>::iterator it = path.begin(); it != path.end(); ++it) {
        nodeIdToLongLat(tmp_longLat,array,*it);
        nodePath.push_back(tmp_longLat[0]);
        nodePath.push_back(tmp_longLat[1]);
    }
}

void generateReponse(std::vector<double> &path, std::string& response, uint64_t distance){
    // response has to be in [[lat,long],...] format, so (long, lat) is swapped
    response += "{\"path\":[";
    response += "[" + std::to_string(path.at(1)) + "," + std::to_string(path.at(0)) + "]";
    for(int i=1; i<path.size()/2; ++i){
        response += ",[" + std::to_string(path.at(2*i+1)) + "," + std::to_string(path.at(2*i)) + "]";
    
    }
    response += "],\n";
    response += "\"dist\": ";
    response += std::to_string(distance); 
    response += "}";
}



class PathAlgorithm {
    public:
        virtual void getPath(uint64_t startPoint, uint64_t endPoint, std::vector<uint64_t> &path) = 0;
        virtual uint64_t getDist(uint64_t startPoint, uint64_t endPoint) = 0;
        virtual uint64_t getNewDist(uint64_t startPoint, uint64_t endPoint) = 0;
};

class FirstDijkstra: public PathAlgorithm{
    public:
        FirstDijkstra(AdjacencyArray &array);
        void getPath(uint64_t startPoint, uint64_t endPoint, std::vector<uint64_t> &path);
        uint64_t getDist(uint64_t startPoint, uint64_t endPoint);
        uint64_t getNewDist(uint64_t startPoint, uint64_t endPoint);
    private:
        void generatePath(uint64_t startPoint, uint64_t endPoint, std::vector<uint64_t> &path);
        void fillMaps(uint64_t startPoint, uint64_t endPoint);
        
        AdjacencyArray &array;
        
        std::map<uint64_t,uint64_t> distance;
        std::map<uint64_t,uint64_t> previous;
};

class SecondDijkstra: public PathAlgorithm{
    public:
        SecondDijkstra(AdjacencyArray &array);
        void getPath(uint64_t startPoint, uint64_t endPoint, std::vector<uint64_t> &path);
        uint64_t getDist(uint64_t startPoint, uint64_t endPoint);
        uint64_t getNewDist(uint64_t startPoint, uint64_t endPoint);
        void prepareDatastructures();
    private:
        std::vector<uint64_t> distance;
        std::vector<HeapElement> heap;
        std::vector<uint64_t> visited;
        std::vector<uint64_t> prev;
        AdjacencyArray &adjArray;
        uint64_t constLngDist;
        std::vector<uint64_t> constLatDist;
};


FirstDijkstra::FirstDijkstra(AdjacencyArray &array) : array(array){}

void FirstDijkstra::getPath(uint64_t startPoint, uint64_t endPoint, std::vector<uint64_t> &path){
    distance.clear();
    previous.clear();
    generatePath(startPoint, endPoint, path);
}


uint64_t FirstDijkstra::getNewDist(uint64_t startPoint, uint64_t endPoint){
    return distance.at(endPoint);
}

uint64_t FirstDijkstra::getDist(uint64_t startPoint, uint64_t endPoint){
    return distance.at(endPoint);
}

void FirstDijkstra::generatePath(uint64_t startPoint, uint64_t endPoint, std::vector<uint64_t> &path){
    

    // check if one of the points is on land
    if(array.nodes.at(startPoint) || array.nodes.at(endPoint)){
        std::cout << "one point is on land" << std::endl;
        if(array.nodes.at(startPoint)){
            std::cout << startPoint << " on land" << std::endl;
        }
        if(array.nodes.at(endPoint)){
            std::cout << endPoint << " on land" << std::endl;
        }
        return;
    }    



    // reset map
    for (int i = 0; i < array.nodes.size(); i++){
        distance.insert({i, UINT64_MAX});
        previous.insert({i, UINT64_MAX});
    }
    distance.at(startPoint) = 0;

    std::cout << "start maps fill" << std::endl;
    fillMaps(startPoint, endPoint);


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

void FirstDijkstra::fillMaps(uint64_t startPoint, uint64_t endPoint){

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
            return;
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

                // do not expand node if it has higher distance than 
                if(currDistance <= currTargetDistance){
                    pq.push(std::make_pair(currDistance, currTargetNode));
                }
            }
        }
    }
}

SecondDijkstra::SecondDijkstra(AdjacencyArray &array) : adjArray(array), prev(array.width*array.height, UINT64_MAX){
    prepareDatastructures();
    
    // distance between (i,0) and (i,1)
    constLngDist = nodeDistance(adjArray, 0, 1);
    for(uint64_t i = 0; i<adjArray.height; ++i){
        // distance between (0,i) and (1,i)
        constLatDist.push_back(nodeDistance(adjArray, i, adjArray.height+i));
    }
}

void SecondDijkstra::prepareDatastructures(){
    std::vector<uint64_t> initDist(adjArray.width*adjArray.height, UINT64_MAX);
    distance = std::move(initDist);
    // no need to reset prev
    heap.clear();
}

/**
 * @brief For this to work, prepareDatastructures need to be called first
 * 
 * @param startPoint 
 * @param endPoint 
 * @return uint64_t 
 */
uint64_t SecondDijkstra::getDist(uint64_t startPoint, uint64_t endPoint){

    heap.push_back(HeapElement{startPoint, UINT64_MAX, 0});
    heap.push_back(HeapElement{endPoint, UINT64_MAX, UINT64_MAX});

    std::make_heap(heap.begin(), heap.end());

    HeapElement front;

    //std::cout << "survived1" << "\n";

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

        //std::cout << "survived2" << "\n";

        distance.at(front.nodeIdx) = front.dist;
        // std::cout << prev.size() << "\n";
        // std::cout << distance.size() << "\n";
        // std::cout << "survived2.1" << "\n";
        prev.at(front.nodeIdx) = front.prev;
        //std::cout << "survived2.2" << "\n";
        visited.push_back(front.nodeIdx);

        //std::cout << "survived3" << "\n";

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
            return distance.at(front.nodeIdx);
        }

    }

}

uint64_t SecondDijkstra::getNewDist(uint64_t startPoint, uint64_t endPoint){
    return distance.at(endPoint);
}

void SecondDijkstra::getPath(uint64_t startPoint, uint64_t endPoint, std::vector<uint64_t> &path){
    
    // check if one of the points is on land
    if(adjArray.nodes.at(startPoint) || adjArray.nodes.at(endPoint)){
        std::cout << "one point is on land" << std::endl;
        if(adjArray.nodes.at(startPoint)){
            std::cout << startPoint << " on land" << std::endl;
        }
        if(adjArray.nodes.at(endPoint)){
            std::cout << endPoint << " on land" << std::endl;
        }
        return;
    }    

    
    prepareDatastructures();
    getDist(startPoint, endPoint);


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
//int main() {
//    AdjacencyArray adjArray;
//    // load
//    loadAdjacencyArray(adjArray, "data/worldGrid_1415_707.save");
//
//    std::vector<uint64_t> path;
//    generatePath(path, 1001, 16001, adjArray);
//    std::vector<double> nodePath;
//    generateNodePath(nodePath, path, adjArray);
//
//
//    std::string response;
//    generateReponse(nodePath, response);
//    std::cout << response << std::endl;
//
//    //int counter_one = 0;
//    //int counter_zero = 0;
//    //for (int i = 0; i < adjArray.nodes.size(); i++){
//    //    if(adjArray.nodes.at(i) == 0){
//    //        counter_zero++;
//    //    }else
//    //    {
//    //        counter_one++;
//    //    }
//    //}
//    //std::cout << counter_zero << " " << counter_one << " " << counter_one + counter_zero << std::endl;
//    //std::cout << adjArray.edges.size() << std::endl;
//    //testLatLongDistance();
//}