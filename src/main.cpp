# include <iostream>
# include <string>
# include <vector>
# include <map>
# include <algorithm>
//# include <ios>

# include "osmpbfreader.h"

using namespace CanalTP;

struct Way {
    uint64_t osmid;
    Tags tags;
    std::vector<uint64_t> refs;
    uint64_t firstNode;
    uint64_t lastNode;
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

        // IF block not mandatory if file contains only coaslines
        if(tags.size() > 0)
        {
            if(tags.count("natural") > 0){
                if(tags.at("natural").compare("coastline") == 0){
                    //std::cout << "found stuff\n";
                    ways.push_back({osmid, tags, refs, refs.front(), refs.back()});
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

// add element to map
// if key already exists return other osmid else return same osmid
uint64_t add_element_map(std::map<uint64_t, std::vector<uint64_t>> mymap, uint64_t key, uint64_t value){

    // node is already in one of the ways
    if(mymap.count(key) > 0)
    {
        mymap.at(key).push_back(value);
        return mymap.at(key).front();
    }else {
        mymap.insert({key, std::vector<uint64_t>{value}});
    }
    return value;
}



/**
 * Adds around 2.5GB of data to the nodes and ways vector combined 
 **/
void load_coastline_data(std::string pbfPath, std::vector<SimpleNode> &nodes, std::vector<Way> &ways){
    CoastlineWaysExtractor wayExtractor(ways);
    read_osm_pbf(pbfPath, wayExtractor);

    for (uint64_t i = 0; i < ways.size(); i++)
    {
        for(uint64_t ref : ways.at(i).refs){
            nodes.push_back(SimpleNode{ref});
        }
        //ways.at(lower_wayid).refs.insert(ways.at(lower_wayid).refs.end(), ways.at(i).refs.begin(), ways.at(i).refs.end());
    }
    
    
    
    uint64_t tempLastNode = 0;
    uint64_t tempiLastNode = 0;
    uint64_t tempjLastNode = 0;

    uint64_t iFirstNode = 0;
    uint64_t iLastNode = 0;
    int waysize = 0;
    int nways = ways.size();

    for (uint64_t i = 0; i < ways.size(); i++)
    {
        // if way is not closed
        iFirstNode = ways.at(i).firstNode; 
        iLastNode = ways.at(i).lastNode;
        while (ways.at(i).lastNode != ways.at(i).firstNode)
        {
            //std::cout << ways.at(i).lastNode << " " << ways.at(i).firstNode << std::endl;

            tempLastNode = ways.at(i).lastNode;
            tempiLastNode = ways.at(i).lastNode;

            // FIND NEXT WAY
            for (uint64_t j = 0; j < ways.size(); j++)
            {
                if(tempLastNode == ways.at(j).firstNode)
                {
                    ways.at(i).refs.insert(ways.at(i).refs.end(), ways.at(j).refs.begin(), ways.at(j).refs.end());
                    ways.at(i).lastNode = ways.at(j).lastNode;
                    tempiLastNode = ways.at(i).lastNode;
                    tempjLastNode = ways.at(j).lastNode;
                    tempLastNode = ways.at(j).lastNode;

                    iFirstNode = ways.at(i).firstNode; 
                    iLastNode = ways.at(i).lastNode;
                    waysize = ways.at(i).refs.size();

                    ways.erase(ways.begin() + j);
                    nways = ways.size();
                    std::cout << nways << std::endl;
                    std::cout  << std::endl;
                    break;
                }
            }
            //ways.at(i).refs.insert(ways.at(i).refs.end(), ways.at(i).)

        }
    }


    //for(Way way : ways){
    //    //std::cout << way.refs.front() << std::endl;
    //    //std::cout << way.refs.back() << std::endl;
    //    //std::cout << std::endl;



    //    //front_wayid = add_element_map(node_wayid, way.refs.front(), way.osmid);
    //    //back_wayid = add_element_map(node_wayid, way.refs.back(), way.osmid);
    //    //lower_wayid = std::min(front_wayid, back_wayid);


    //    //origin_osmid.insert({way.osmid, add_element_map(node_wayid, way.refs.front(), way.osmid)})
    //    //if(way.osmid == add_element_map(node_wayid, way.refs.front(), way.osmid))
    //    //{
    //    //    origin_osmid.insert({way.osmid, way.osmid});
    //    //}



    //    for(uint64_t ref : way.refs){
    //        nodes.push_back(SimpleNode{ref});
    //    }
    //}

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

int main(int argc, char** argv) {
    std::cout << "Hello World\n";
    if(argc != 3) {
        std::cout << "Usage: " << argv[0] << " file_to_read.osm.pbf" << " " << "file_to save" << std::endl;
        return 1;
    }

    std::vector<SimpleNode> nodes;
    std::vector<Way> ways;
    load_coastline_data(argv[1], nodes, ways);
    std::cout << "Nodes vector takes up " << (nodes.size()*3.*4.)/1e6 << " MB\n";
    
    save_coastline_to_geojson(argv[2], nodes, ways);
    
    //CoastlineWaysExtractor wayExtractor;
    //read_osm_pbf(argv[1], wayExtractor);

    //std::cout << "We read " << counter.nodes << " nodes, " << counter.ways << " ways and " << counter.relations << " relations" << std::endl;
    //std::cout << "Coastline size: " << wayExtractor.coastline_ways.size() << std::endl;
    //const std::vector<uint64_t> &refs = wayExtractor.coastline_ways.at(0).refs;
    //std::cout << refs.at(0) << std::endl;
    return 0;

}