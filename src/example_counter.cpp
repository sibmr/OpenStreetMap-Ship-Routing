/*
Copyright (c) 2012, Canal TP
This is an example file, do whatever you want with it! (for example if you are in Paris, invite us for a beer)

This shows the simplest way to use the osm.pbf reader. It just counts the number of objects in the file.

To build this file :
g++ -O2 -o counter example_counter.cc -losmpbf -lprotobuf

To run it:
./counter path_to_your_data.osm.pbf
*/

#include "osmpbfreader.h"

using namespace CanalTP;

// We need to define a visitor with three methods that will be called while the file is read
struct Counter {
    // Three integers count how many times each object type occurs
    int nodes;
    int ways;
    int relations;

    Counter() : nodes(0), ways(0), relations(0) {}

    // This method is called every time a Node is read
    void node_callback(uint64_t /*osmid*/, double /*lon*/, double /*lat*/, const Tags &/*tags*/){
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

int main(int argc, char** argv) {
     if(argc != 2) {
         std::cout << "Usage: " << argv[0] << " file_to_read.osm.pbf" << std::endl;
         return 1;
     }


     // Let's read that file !
     Counter counter;
     read_osm_pbf(argv[1], counter);
     std::cout << "We read " << counter.nodes << " nodes, " << counter.ways << " ways and " << counter.relations << " relations" << std::endl;
     return 0;
}
