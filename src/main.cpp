# include <iostream>
# include <string>
# include <vector>
# include <map>

# include "osmpbfreader.h"

using namespace CanalTP;

struct Way {
    uint64_t osmid;
    const Tags tags;
    const std::vector<uint64_t> refs;
};

struct Node {
    uint64_t osmid;
    double logitude;
    double latitude;
    const Tags tags;
};

struct CoastlineWaysExtractor {
    
    //std::vector<Way> coastline_ways;

    std::map<uint64_t, std::vector<uint64_t>> coastline_list;

    CoastlineWaysExtractor() {}

    // This method is called every time a Way is read
    // refs is a vector that contains the reference to the nodes that compose the way
    void way_callback(uint64_t osmid, const Tags &tags, const std::vector<uint64_t> &refs){

        // IF block not mandatory if file contains only coaslines
        if(tags.size() > 0)
        {
            if(tags.count("natural") > 0){
                if(tags.at("natural").compare("coastline") == 0){
                    for (auto i = 0; i < refs.size(); i++)
                    {
                        // node is already in one of the ways
                        if(coastline_list.count(refs.at(i)) > 0)
                        {
                           coastline_list.at(refs.at(i)).push_back(osmid);
                        }else {
                            coastline_list.insert({refs.at(i), std::vector<uint64_t>{osmid}});
                        }
                    }
                }
            }
        }
    }
    void relation_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const References & /*refs*/){}
    void node_callback(uint64_t /*osmid*/, double /*longitude*/, double latitude, const Tags &/*tags*/){}
};

struct CoastlineNodesExtractor {
    
    std::map<uint64_t, Node*> coastlineNodes;

    CoastlineNodesExtractor(std::map<uint64_t, Node*> &nodeMap) : coastlineNodes(nodeMap) {}

    void node_callback(uint64_t osmid, double longitude, double latitude, const Tags &tags){
        // find node id in the map and link the next node to it
        std::map<uint64_t, Node*>::iterator it = coastlineNodes.find(osmid); 
        if(it != coastlineNodes.end()){
            it->second = new Node({osmid, longitude, latitude, tags});
        }
    }
    void relation_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const References & /*refs*/){}
    void way_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const std::vector<uint64_t> &/*refs*/){}
};



//struct Counter {
//    // Three integers count how many times each object type occurs
//    int nodes;
//    int ways;
//    int relations;
//    std::vector<Way> coastline_ways;
//
//
//    Counter() : nodes(0), ways(0), relations(0) {}
//
//    // This method is called every time a Node is read
//    void node_callback(uint64_t /*osmid*/, double longitude, double latitude, const Tags &tags){
//        //if(tags.size() > 0)
//        //{
//        //    if(tags.count("natural") > 0){
//        //        if(tags.at("natural").compare("coastline") == 0){
//        //            std::cout << "found stuff\n";
//        //        }
//        //    }
//        //    //std::cout << longitude << " " << latitude << " " << tags.at(0).at(0) << std::endl;
//        //}
//        //else
//        //{
//        //    //std::cout << longitude << " " << latitude << std::endl;
//        //}
//        ++nodes;
//    }
//
//    // This method is called every time a Way is read
//    // refs is a vector that contains the reference to the nodes that compose the way
//    void way_callback(uint64_t osmid, const Tags &tags, const std::vector<uint64_t> &refs){
//        if(tags.size() > 0)
//        {
//            if(tags.count("natural") > 0){
//                if(tags.at("natural").compare("coastline") == 0){
//                    std::cout << "found stuff\n";
//                    //coastline_ways.push_back(Way({osmid, tags, refs}));
//                }
//            }
//            //std::cout << longitude << " " << latitude << " " << tags.at(0).at(0) << std::endl;
//        }
//        else
//        {
//            //std::cout << longitude << " " << latitude << std::endl;
//        }
//        ++ways;
//    }
//
//    // This method is called every time a Relation is read
//    // refs is a vector of pair corresponding of the relation type (Node, Way, Relation) and the reference to the object
//    void relation_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const References & /*refs*/){
//        ++relations;
//    }
//};

void coast_polyline_assembly(std::string pbfPath, std::map<uint64_t, Node*> &nodeMap, std::vector<Way> &ways){
    CoastlineWaysExtractor wayExtractor;
    read_osm_pbf(pbfPath, wayExtractor);

    //for(Way way : wayExtractor.coastline_ways){
    //    for(uint64_t ref : way.refs){
    //        //nodeMap.insert(std::pair<uint64_t, Node*>(ref, nullptr));
    //        //std::cout << ref << std::endl;
    //    }
    //}

    //CoastlineNodesExtractor nodeExtractor(nodeMap);
    //read_osm_pbf(pbfPath, nodeExtractor);
}

void test(std::vector<std::string> &testv){
    std::vector<std::string> newV {"changed"};
    testv = newV;
}

int main(int argc, char** argv) {
    std::cout << "Hello World\n";
    if(argc != 2) {
        std::cout << "Usage: " << argv[0] << " file_to_read.osm.pbf" << std::endl;
        return 1;
    }

    //std::vector<std::string> oldV {"ok"};
    //std::cout << oldV.at(0) << std::endl;
    //test(oldV);
    //std::cout << oldV.at(0) << std::endl;

    std::map<uint64_t, Node *> nodes;
    std::vector<Way> ways;
    
    coast_polyline_assembly(argv[1], nodes, ways);

    //Counter counter; 
    //read_osm_pbf(argv[1], counter);
    //CoastlineWaysExtractor wayExtractor;
    //read_osm_pbf(argv[1], wayExtractor);

    //std::cout << "We read " << counter.nodes << " nodes, " << counter.ways << " ways and " << counter.relations << " relations" << std::endl;
    //std::cout << "Coastline size: " << wayExtractor.coastline_ways.size() << std::endl;
    //const std::vector<uint64_t> &refs = wayExtractor.coastline_ways.at(0).refs;
    //std::cout << refs.at(0) << std::endl;
    return 0;

}