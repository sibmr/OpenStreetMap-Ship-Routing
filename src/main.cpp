# include <iostream>
# include <string>
# include <vector>
# include <map>
//# include <ios>

# include "osmpbfreader.h"

using namespace CanalTP;

struct Way {
    uint64_t osmid;
    const Tags tags;
    const std::vector<uint64_t> refs;
};

struct Node {
    uint64_t osmid;
    double longitude;
    double latitude;
    const Tags tags;
};

struct SimpleNode {
    uint64_t osmid;
    double longitude;
    double latitude;
};

struct {
    bool operator()(SimpleNode a, SimpleNode b) const { return a.osmid < b.osmid; }
} compareSimpleNode;

struct CoastlineWaysExtractor {
    
    std::vector<Way> &ways;

    CoastlineWaysExtractor(std::vector<Way> &coastline_ways) : ways(coastline_ways) {}

    // This method is called every time a Way is read
    // refs is a vector that contains the reference to the nodes that compose the way
    void way_callback(uint64_t osmid, const Tags &tags, const std::vector<uint64_t> &refs){
        if(tags.size() > 0)
        {
            if(tags.count("natural") > 0){
                if(tags.at("natural").compare("coastline") == 0){
                    //std::cout << "found stuff\n";
                    ways.push_back({osmid, tags, refs});
                }
            }
        }
    }
    void relation_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const References & /*refs*/){}
    void node_callback(uint64_t /*osmid*/, double /*longitude*/, double latitude, const Tags &/*tags*/){}
};

struct CoastlineNodesExtractor {
    
    static const uint64_t MAX_NODES = 57000000;
    std::vector<SimpleNode> &nodes;

    CoastlineNodesExtractor(std::vector<SimpleNode> &coastlineNodes): nodes(coastlineNodes){}

    void node_callback(uint64_t osmid, double longitude, double latitude, const Tags &tags){
        // binary-search the node with osmid
        std::vector<SimpleNode>::iterator it = std::lower_bound(nodes.begin(), nodes.end(), SimpleNode{osmid}, compareSimpleNode);
        if(it->osmid == osmid){
            it->longitude = longitude;
            it->latitude = latitude;
        }
    }
    void relation_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const References & /*refs*/){}
    void way_callback(uint64_t /*osmid*/, const Tags &/*tags*/, const std::vector<uint64_t> &/*refs*/){}
};

/**
 * Adds around 2.5GB of data to the nodes and ways vector combined 
 **/
void load_coastline_data(std::string pbfPath, std::vector<SimpleNode> &nodes, std::vector<Way> &ways){
    CoastlineWaysExtractor wayExtractor(ways);
    read_osm_pbf(pbfPath, wayExtractor);

    for(Way way : ways){
        for(uint64_t ref : way.refs){
            nodes.push_back(SimpleNode{ref});
        }
    }

    std::sort(nodes.begin(), nodes.end(), compareSimpleNode);

    CoastlineNodesExtractor nodeExtractor(nodes);
    read_osm_pbf(pbfPath, nodeExtractor);
}

void save_coastline_to_geojson(std::string geoJsonPath, std::vector<SimpleNode> &nodes, std::vector<Way> &ways){
    std::ofstream file;

    file.open(geoJsonPath, std::ios::out | std::ios::trunc);
    file.exceptions(file.exceptions() | std::ios::failbit | std::ifstream::badbit);

    file <<     "{ \"type\": \"FeatureCollection\",\n";
    file <<     "  \"features\": [\n";

    uint64_t count = 0;
    bool first_way = true;
    for(Way way : ways){
        if(!first_way){file << ",";}
        file << "   {\"type\": \"Feature\",\n"
                "    \"geometry\": {\n"
                "       \"type\": \"LineString\",\n"
                "       \"coordinates\": [\n";

        bool first_node = true;
        for(SimpleNode node : nodes){
            std::vector<SimpleNode>::iterator it = std::lower_bound(nodes.begin(), nodes.end(), SimpleNode{node.osmid}, compareSimpleNode);
            if(!first_node){file << ",";};
            file << "[" << it->longitude << "," << it->latitude << "]\n";
            first_node=false;
        }

        file << "   ]\n";
        file << "   }}";
        
        count++;
        if(count%1 == 0){file.flush();}
        if(count%1 == 0){break;}
        first_way = false;
    }
    file <<     "]}\n" << std::endl;
    
}

void test(std::vector<std::string> &testv){
    std::vector<std::string> newV {"changed"};
    testv = newV;
}

int main(int argc, char** argv) {
    std::cout << "Hello World\n";
    // if(argc != 2) {
    //     std::cout << "Usage: " << argv[0] << " file_to_read.osm.pbf" << std::endl;
    //     return 1;
    // }

    std::vector<std::string> oldV {"ok"};
    std::cout << oldV.at(0) << std::endl;
    test(oldV);
    std::cout << oldV.at(0) << std::endl;

    std::vector<SimpleNode> nodes;
    std::vector<Way> ways;
    load_coastline_data("/home/simonb/git/OpenStreetMap-Ship-Routing/data/antarctica-latest.osm.pbf", nodes, ways);
    std::cout << "Nodes vector takes up " << (nodes.size()*3.*4.)/1e6 << " MB\n";
    
    save_coastline_to_geojson("/home/simonb/git/OpenStreetMap-Ship-Routing/data/planet-coastlines.json", nodes, ways);
    
    //CoastlineWaysExtractor wayExtractor;
    //read_osm_pbf(argv[1], wayExtractor);

    //std::cout << "We read " << counter.nodes << " nodes, " << counter.ways << " ways and " << counter.relations << " relations" << std::endl;
    //std::cout << "Coastline size: " << wayExtractor.coastline_ways.size() << std::endl;
    //const std::vector<uint64_t> &refs = wayExtractor.coastline_ways.at(0).refs;
    //std::cout << refs.at(0) << std::endl;
    return 0;

}