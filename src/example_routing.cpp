/* 
Copyright (c) 2012, Canal TP
This is an example file, do whatever you want with it! (for example if you are in Paris, invite us for a beer)

A slightly more advanced example : we want to use OSM data to extract a graph that represents the road network 

Build the example with (we're to lazy to make it compatible with older standards):
g++ example_routing.cc -O2 -losmpbf -lprotobuf -std=c++0x -o routing 
To run it:
./routing path_to_your_data.osm.pbf
*/

#include <unordered_map>
#include "osmpbfreader.h"
using namespace CanalTP;

// We keep every node and the how many times it is used in order to detect crossings
struct Node {
    public:
        Node(double lon = 0, double lat = 0) : uses(0), lon_m(lon), lat_m(lat){}

        int32_t uses;
        double lon_m;
        double lat_m;
};

struct Routing {
    // Map that stores all the nodes read
    std::unordered_map<uint64_t, Node> nodes;

    // Stores all the nodes of all the ways that are part of the road network
    std::vector< std::vector<uint64_t> > ways;

    // This method is called every time a Node is read
    void node_callback(uint64_t osmid, double lon, double lat, const Tags &/*tags*/){
        this->nodes[osmid] = Node(lon, lat);
    }

    // This method is called every time a Way is read
    void way_callback(uint64_t /*osmid*/, const Tags &tags, const std::vector<uint64_t> &refs){
        // If the way is part of the road network we keep it
        // There are other tags that correspond to the street network, however for simplicity, we don't manage them
        // Homework: read more properties like oneways, bicycle lanes…
        if(tags.find("highway") != tags.end()){
            ways.push_back(refs);
        }
    }

    // Once all the ways and nodes are read, we count how many times a node is used to detect intersections
    void count_nodes_uses() {
        for(std::vector<uint64_t> refs : ways){
            for(uint64_t ref : refs){
                nodes.at(ref).uses++;
            }
            // make sure that the last node is considered as an extremity
            nodes.at(refs.back()).uses++;
        }
    }

    // Returns the source and target node of the edges
    std::vector< std::pair<uint64_t, uint64_t> > edges(){
        std::vector< std::pair<uint64_t, uint64_t> > result;

        for(std::vector<uint64_t> refs : ways){
            if(refs.size() > 0){
                uint64_t source = refs[0];
                for(size_t i = 1; i < refs.size(); ++i){
                    uint64_t current_ref = refs[i];
                    // If a node is used more than once, it is an intersection, hence it's a node of the road network graph
                    if(nodes.at(current_ref).uses > 1){
                        // Homework: measure the length of the edge
                        uint64_t target = current_ref;
                        result.push_back(std::make_pair(source, target));
                        source = target;
                    }
                }
            }
        }
        return result;
    }

    // We don't care about relations
    void relation_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const References & /*refs*/){}
};

int main(int argc, char** argv) {
    if(argc != 2) {
        std::cout << "Usage: " << argv[0] << " file_to_read.osm.pbf" << std::endl;
        return 1;
    }

    // Let's read that file !
    Routing routing;
    read_osm_pbf(argv[1], routing);
    std::cout << "We read " << routing.nodes.size() << " nodes and " << routing.ways.size() << " ways" << std::endl;
    routing.count_nodes_uses();    
    std::cout << "The routing graph has " << routing.edges().size() << " edges" << std::endl;
    return 0;
}

