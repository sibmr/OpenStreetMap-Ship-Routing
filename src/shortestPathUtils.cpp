#ifndef SHORTEST_PATH_UTILS
#define SHORTEST_PATH_UTILS

# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <math.h>
# include <map>
# include <queue>
# include <algorithm>

# include "fileUtils.cpp"

struct Edge
{
    uint64_t edgeId;
    uint64_t v1;
    uint64_t v2;
    uint64_t edgeDistance;
    std::vector<uint64_t> shortcutPathEdges;
};

std::vector<Edge> copyEdgeVector(const std::vector<Edge> &edgeVector){
    std::vector<Edge> copy;
    for (const Edge &edge : edgeVector){
        copy.push_back(Edge{
            .edgeId = edge.edgeId,
            .v1 = edge.v1,
            .v2 = edge.v2,
            .edgeDistance = edge.edgeDistance,
            .shortcutPathEdges = std::vector<uint64_t>(edge.shortcutPathEdges)
        });
    }
    return copy;
}

void writeEdgeVector(std::ofstream &outputFile, std::vector<Edge> &edgeVector){
    uint64_t edgeVectorSize = edgeVector.size();
    writeSimpleValue(outputFile, edgeVectorSize);
    for(Edge &edge : edgeVector){
        writeSimpleValue(outputFile, edge.edgeId);
        writeSimpleValue(outputFile, edge.v1);
        writeSimpleValue(outputFile, edge.v2);
        writeSimpleValue(outputFile, edge.edgeDistance);
        writeVectorOfSimpleValue(outputFile, edge.shortcutPathEdges);
    }
}

void readEdgeVector(std::ifstream &inputFile, std::vector<Edge> &edgeVector){
    uint64_t edgeVectorSize;
    readSimpleValue(inputFile, edgeVectorSize);
    edgeVector.resize(edgeVectorSize);
    for(Edge edge : edgeVector){
        readSimpleValue(inputFile, edge.edgeId);
        readSimpleValue(inputFile, edge.v1);
        readSimpleValue(inputFile, edge.v2);
        readSimpleValue(inputFile, edge.edgeDistance);
        readVectorOfSimpleValue(inputFile, edge.shortcutPathEdges);
    }
}

struct AdjacencyArray {
    double longLow, latLow, longHigh, latHigh;
    uint64_t width, height;
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> edges;
    std::vector<uint64_t> edgeIds;
    std::vector<uint64_t> distances;
    std::vector<uint64_t> rank;
    std::vector<Edge> allEdgeInfo;
    std::vector<bool> nodes;

    /**
     * @brief Construct a new Adjacency Array object for buildGraph
     */
    AdjacencyArray(){};

    /**
     * @brief Deep Copy an existing Adjacency Array object
     * 
     * @param adjArray 
     */
    AdjacencyArray(const AdjacencyArray& adjArray){
        longLow = adjArray.longLow;     latLow = adjArray.latLow;
        longHigh = adjArray.longHigh;   latHigh = adjArray.latHigh;
        width = adjArray.width;         height = adjArray.height;
        offsets = std::vector<uint64_t>(adjArray.offsets);
        edges = std::vector<uint64_t>(adjArray.edges);
        edgeIds = std::vector<uint64_t>(adjArray.edgeIds);
        distances = std::vector<uint64_t>(adjArray.distances);
        rank = std::vector<uint64_t>(adjArray.rank);
        allEdgeInfo = copyEdgeVector(adjArray.allEdgeInfo);
        nodes = std::vector<bool>(adjArray.nodes);
    };



    /**
     * @brief Construct a new Adjacency Array object
     * 
     * @param path path to input .graph file
     */
    AdjacencyArray(std::string path) : offsets(), edges(), distances(), rank(),  nodes(){

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
        readSimpleValue(adjacency_input_file, longLow);
        readSimpleValue(adjacency_input_file, latLow);
        readSimpleValue(adjacency_input_file, longHigh);
        readSimpleValue(adjacency_input_file, latHigh);

        // save number of nodes per direction (width, height)
        readSimpleValue(adjacency_input_file, width);
        readSimpleValue(adjacency_input_file, height);

        // save offset
        readVectorOfSimpleValue(adjacency_input_file, offsets);

        // read edges
        readVectorOfSimpleValue(adjacency_input_file, edges);

        // read distances
        readVectorOfSimpleValue(adjacency_input_file, distances);

        // read rank
        readVectorOfSimpleValue(adjacency_input_file, rank);

        // read edgeIds
        readVectorOfSimpleValue(adjacency_input_file, edgeIds);

        // read allEdgeInfo
        readEdgeVector(adjacency_input_file, allEdgeInfo);

        // read nodes
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

    /**
     * @brief save adjacency array to disk
     * 
     * @param array AdjacencyArray struct that is stored
     * @param path  path to storage location
     * 
     * output file format
     * longLow       - double
     * latLow        - double
     * longHigh      - double
     * latHigh       - double
     * width         - uint64_t
     * height        - uint64_t
     * offset_size   - uint64_t
     * offsets       - uint64_t     (offset_size many)
     * edges_size    - uint64_t
     * edges         - uint64_t     (edges_size many)
     * nodes_size    - uint64_t
     * nodes          - bool        (nodes_size many)
     */
    void writeToDisk(std::string path){
        std::ofstream adjacency_output_file;

        adjacency_output_file.open(path, std::ios::out | std::ios::trunc);
        adjacency_output_file.exceptions(adjacency_output_file.exceptions() | std::ios::failbit | std::ifstream::badbit);

        // write globe size
        writeSimpleValue(adjacency_output_file, longLow);
        writeSimpleValue(adjacency_output_file, latLow);
        writeSimpleValue(adjacency_output_file, longHigh);
        writeSimpleValue(adjacency_output_file, latHigh);

        // save number of nodes per direction (width, height)
        writeSimpleValue(adjacency_output_file, width);
        writeSimpleValue(adjacency_output_file, height);

        // save offset vector
        writeVectorOfSimpleValue(adjacency_output_file, offsets);

        // save edges
        writeVectorOfSimpleValue(adjacency_output_file, edges);

        // save distances
        writeVectorOfSimpleValue(adjacency_output_file, distances);

        // save rank
        writeVectorOfSimpleValue(adjacency_output_file, rank);

        // save edgeIds
        writeVectorOfSimpleValue(adjacency_output_file, edgeIds);

        // save edgeVector
        writeEdgeVector(adjacency_output_file, allEdgeInfo);

        // save nodes
        uint64_t nodes_size = nodes.size();
        writeSimpleValue(adjacency_output_file, nodes_size);
        std::copy(nodes.begin(), nodes.end(), std::ostreambuf_iterator<char>(adjacency_output_file));
        adjacency_output_file.flush();

        adjacency_output_file.close();
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
void nodeIdToLongLat(std::array<double,2> &result, double &longLow, double &latLow, double &longHigh, double &latHigh, uint64_t &width, uint64_t &height, uint64_t id){
    const double longStep = (longHigh-longLow)/(width+1);
    const double latStep = (latHigh-latLow)/(height+1);

    result[0] =  longLow + std::floor(id/height) * longStep + longStep/2;
    result[1] =  latLow + (id%height) * latStep + latStep/2;
    return;

}

void nodeIdToArrayIdx(AdjacencyArray &array, uint64_t id, uint64_t &lng_out, uint64_t &lat_out){
    lng_out = std::floor(id/array.height);
    lat_out = id%array.height;
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
uint64_t nodeDistance(double longLow, double latLow, double longHigh, double latHigh, uint64_t width, uint64_t height, uint64_t node1, uint64_t node2){
    std::array<double,2> longLat_1;
    std::array<double,2> longLat_2;

    nodeIdToLongLat(longLat_1, longLow, latLow, longHigh, latHigh, width, height, node1);
    nodeIdToLongLat(longLat_2, longLow, latLow, longHigh, latHigh, width, height, node2);

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

#endif
