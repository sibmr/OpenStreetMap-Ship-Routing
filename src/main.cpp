# include <iostream>
# include <string>
# include <vector>

# include "osmpbfreader.h"

using namespace CanalTP;

struct Way {
    uint64_t osmid;
    const Tags &tags;
    const std::vector<uint64_t> &refs;
};

struct CoastlineWaysExtractor {
    
    std::vector<Way> coastline_ways;

    CoastlineWaysExtractor() {}

    // This method is called every time a Way is read
    // refs is a vector that contains the reference to the nodes that compose the way
    void way_callback(uint64_t osmid, const Tags &tags, const std::vector<uint64_t> &refs){
        if(tags.size() > 0)
        {
            if(tags.count("natural") > 0){
                if(tags.at("natural").compare("coastline") == 0){
                    //std::cout << "found stuff\n";
                    coastline_ways.push_back({osmid, tags, refs});
                }
            }
        }
    }
    void relation_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const References & /*refs*/){}
    void node_callback(uint64_t /*osmid*/, double /*longitude*/, double latitude, const Tags &/*tags*/){}
};

struct CoastlineNodesExtractor {
    
    std::vector<Way> coastline_ways;

    CoastlineNodesExtractor() {}

    // This method is called every time a Way is read
    // refs is a vector that contains the reference to the nodes that compose the way
    void way_callback(uint64_t osmid, const Tags &tags, const std::vector<uint64_t> &refs){
        if(tags.size() > 0)
        {
            if(tags.count("natural") > 0){
                if(tags.at("natural").compare("coastline") == 0){
                    //std::cout << "found stuff\n";
                    coastline_ways.push_back({osmid, tags, refs});
                }
            }
        }
    }
    void relation_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const References & /*refs*/){}
    void node_callback(uint64_t /*osmid*/, double /*longitude*/, double latitude, const Tags &/*tags*/){}
};



struct Counter {
    // Three integers count how many times each object type occurs
    int nodes;
    int ways;
    int relations;
    std::vector<Way> coastline_ways;


    Counter() : nodes(0), ways(0), relations(0) {}

    // This method is called every time a Node is read
    void node_callback(uint64_t /*osmid*/, double longitude, double latitude, const Tags &tags){
        if(tags.size() > 0)
        {
            if(tags.count("natural") > 0){
                if(tags.at("natural").compare("coastline") == 0){
                    std::cout << "found stuff\n";
                }
            }
            //std::cout << longitude << " " << latitude << " " << tags.at(0).at(0) << std::endl;
        }
        else
        {
            //std::cout << longitude << " " << latitude << std::endl;
        }
        ++nodes;
    }

    // This method is called every time a Way is read
    // refs is a vector that contains the reference to the nodes that compose the way
    void way_callback(uint64_t osmid, const Tags &tags, const std::vector<uint64_t> &refs){
        if(tags.size() > 0)
        {
            if(tags.count("natural") > 0){
                if(tags.at("natural").compare("coastline") == 0){
                    //std::cout << "found stuff\n";
                    coastline_ways.push_back(Way({osmid, tags, refs}));
                }
            }
            //std::cout << longitude << " " << latitude << " " << tags.at(0).at(0) << std::endl;
        }
        else
        {
            //std::cout << longitude << " " << latitude << std::endl;
        }
        ++ways;
    }

    // This method is called every time a Relation is read
    // refs is a vector of pair corresponding of the relation type (Node, Way, Relation) and the reference to the object
    void relation_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const References & /*refs*/){
        ++relations;
    }
};

int main(int argc, char** argv) {
    std::cout << "Hello World\n";
    if(argc != 2) {
        std::cout << "Usage: " << argv[0] << " file_to_read.osm.pbf" << std::endl;
        return 1;
    }

    //Counter counter; 
    CoastlineWaysExtractor wayExtractor;
    read_osm_pbf(argv[1], wayExtractor);

    //std::cout << "We read " << counter.nodes << " nodes, " << counter.ways << " ways and " << counter.relations << " relations" << std::endl;
    std::cout << "Coastline size: " << wayExtractor.coastline_ways.size() << std::endl;
    const std::vector<uint64_t> &refs = wayExtractor.coastline_ways.at(0).refs;
    std::cout << refs.at(0) << std::endl;
    return 0;

}