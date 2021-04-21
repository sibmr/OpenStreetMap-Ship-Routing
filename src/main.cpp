# include <iostream>

# include "lib/osmpbfreader.h"

using namespace CanalTP;




struct Counter {
    // Three integers count how many times each object type occurs
    int nodes;
    int ways;
    int relations;

    Counter() : nodes(0), ways(0), relations(0) {}

    // This method is called every time a Node is read
    void node_callback(uint64_t /*osmid*/, double longitude, double latitude, const Tags &tags){
        if(tags.size() > 0)
        {
            std::cout << longitude << " " << latitude << " " << tags.at(0).at(0) << std::endl;
        }
        else
        {
            std::cout << longitude << " " << latitude << std::endl;
        }
        ++nodes;
    }

    // This method is called every time a Way is read
    // refs is a vector that contains the reference to the nodes that compose the way
    void way_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const std::vector<uint64_t> &/*refs*/){
        ++ways;
    }

    // This method is called every time a Relation is read
    // refs is a vector of pair corresponding of the relation type (Node, Way, Relation) and the reference to the object
    void relation_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const References & /*refs*/){
        ++relations;
    }
};





int main(){
    std::cout << "Hello World\n";

    Counter counter; 
    read_osm_pbf("../data/planet-coastlines.pbf", counter);

    std::cout << "We read " << counter.nodes << " nodes, " << counter.ways << " ways and " << counter.relations << " relations" << std::endl;

    return 0;


}