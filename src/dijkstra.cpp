
# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <math.h>
# include <map>
# include <queue>

struct AdjacencyArray {
    double longLow, latLow, longHigh, latHigh;
    uint64_t width, height;
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> edges;
    std::vector<bool> nodes;
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

void nodeIdToLongLat(std::array<double,2> &result, AdjacencyArray &array, uint64_t id){
    const double longStep = (array.longHigh-array.longLow)/(array.width+1);
    const double latStep = (array.latHigh-array.latLow)/(array.height+1);

    result[0] =  array.longLow + std::floor(id/array.height) * longStep + longStep/2;
    result[1] =  array.latLow + (id%array.height) * latStep + latStep/2;
    return;

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


void fillMaps(
    std::map<uint64_t,uint64_t> &distance, std::map<uint64_t,uint64_t> &previous,
    uint64_t startPoint, uint64_t endPoint, AdjacencyArray &array){



    //auto cmp = [](uint64_t left, uint64_t right, std::map<uint64_t, uint64_t> &dis){return &dis.at(left) < &dis.at(right);};
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


void generatePath(std::vector<uint64_t> &path,  uint64_t startPoint, uint64_t endPoint, AdjacencyArray &array){
    std::map<uint64_t,uint64_t> distance;
    std::map<uint64_t,uint64_t> previous;

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

    fillMaps(distance, previous, 1001, 10001, array);


    if(distance.at(endPoint) < INT64_MAX){
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
    }
}



int main(int argc, char** argv) {
    AdjacencyArray adjArray;
    // load
    loadAdjacencyArray(adjArray, "data/worldGrid_1415_707.save");

    std::vector<uint64_t> result;
    generatePath(result, 1001, 10001, adjArray);

    //int counter_one = 0;
    //int counter_zero = 0;
    //for (int i = 0; i < adjArray.nodes.size(); i++){
    //    if(adjArray.nodes.at(i) == 0){
    //        counter_zero++;
    //    }else
    //    {
    //        counter_one++;
    //    }
    //}
    //std::cout << counter_zero << " " << counter_one << " " << counter_one + counter_zero << std::endl;
    //std::cout << adjArray.edges.size() << std::endl;
    //testLatLongDistance();
}