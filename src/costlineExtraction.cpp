# include <iostream>
# include <string>
# include <vector>
# include <map>
# include <algorithm>
# include <math.h>
# include <chrono>
# include <iomanip>

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
};

/**
 * @brief comparison operator for instances of Node based on their id
 * 
 */
struct {
    bool operator()(Node a, Node b) const { return a.osmid < b.osmid; }
} compareSimpleNode;

/**
 * @brief callback for osmpbfreader library: extracts ways that have the {"natural" : "coastline"} tag
 * 
 */
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

/**
 * @brief callback for osmpbfreader library: after coastline ways are extracted, nodes referenced by those ways are extracted
 * 
 */
struct CoastlineNodesExtractor {
    
    static const uint64_t MAX_NODES = 57000000;
    std::vector<Node> &nodes;

    /**
     * @brief Construct a new Coastline Nodes Extractor object: updates node longitude and latitude based on node-id
     * 
     * @param coastlineNodes this vector needs to be pre-initialized with with empty nodes (containing only node-id) and pre-sorted 
     */
    CoastlineNodesExtractor(std::vector<Node> &coastlineNodes): nodes(coastlineNodes){}

    void node_callback(uint64_t osmid, double longitude, double latitude, const Tags &tags){
        // binary-search the node with osmid
        std::vector<Node>::iterator it = std::lower_bound(nodes.begin(), nodes.end(), Node{osmid}, compareSimpleNode);
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
 * @brief Load coastline data from disk, extract coastlines and prepare coastlines data for saving
 * 
 * @param pbfPath       path where pbf file is located
 * @param nodes         vector of nodes that will be filled with nodes that are part of coastline ways
 * @param ways          vector of ways that will be filled with coastline ways
 * @param export_ways   vector containing indices of ways in the ways vector, specifying the order of export
 * 
 * Adds around 2.5GB of data to the nodes and ways vector combined (on planet_coaslines.pbf)
 */
void load_coastline_data(std::string pbfPath, std::vector<Node> &nodes, std::vector<Way> &ways, std::vector<uint64_t> &export_ways){
    CoastlineWaysExtractor wayExtractor(ways);
    read_osm_pbf(pbfPath, wayExtractor);

    std::map<uint64_t, uint64_t> node_vectid; // store for each first node, which ways starts with it

    for (uint64_t i = 0; i < ways.size(); i++)
    {
        node_vectid.insert({ways.at(i).refs.front(), i});
        for(uint64_t ref : ways.at(i).refs){
            nodes.push_back(Node{ref});
        }

    }
    

    std::cout << "All Nodes pushed" << std::endl;
    
    uint64_t next_way_id = 0;
     
    std::vector<uint64_t> remove_ways = {};
    //std::vector<uint64_t> export_ways = {};
    std::vector<uint64_t>().swap(export_ways);


    for (uint64_t i = 0; i < ways.size(); i++){

        if (std::find(remove_ways.begin(), remove_ways.end(), i) == remove_ways.end()){ // way does not get removed later

            // if way is not closed
            while (ways.at(i).lastNode != ways.at(i).firstNode){
                next_way_id = node_vectid.at(ways.at(i).lastNode);

                ways.at(i).refs.insert(ways.at(i).refs.end(), ways.at(next_way_id).refs.begin(), ways.at(next_way_id).refs.end());
                ways.at(i).lastNode = ways.at(next_way_id).refs.back();
                remove_ways.push_back(next_way_id); // store which ways do not have to be watched further
            }
            export_ways.push_back(i);
        }
    }

    std::cout << "All ways merged. Total ways: " << export_ways.size() <<  std::endl; 


    std::map<uint64_t,uint64_t>().swap(node_vectid); // free up space of temp map


    //std::sort(remove_ways.begin(), remove_ways.end()); // list has to be sorted so the offset (- i) stays consistent

    //for (uint64_t i = 0; i < remove_ways.size(); i++){
    //    ways.erase(ways.begin() + remove_ways.at(i) - i);
    //}

    //std::vector<uint64_t>().swap(remove_ways); // free up space of temp vector

    //std::cout << "Total number of ways: " << ways.size() << std::endl;


    std::sort(nodes.begin(), nodes.end(), compareSimpleNode);

    CoastlineNodesExtractor nodeExtractor(nodes);
    read_osm_pbf(pbfPath, nodeExtractor);
}

/**
 * @brief Save geojson representation of data extracted from pbf file
 * 
 * @param geoJsonPath   path where pbf file is located
 * @param nodes         vector of coastline nodes
 * @param ways          vector of coastline ways
 * @param export_ways   vector specifying ordering of coastline ways
 */
void save_coastline_to_geojson(std::string geoJsonPath, std::vector<Node> &nodes, std::vector<Way> &ways, std::vector<uint64_t> &export_ways){
    std::ofstream file;

    file.open(geoJsonPath, std::ios::out | std::ios::trunc);
    file.exceptions(file.exceptions() | std::ios::failbit | std::ifstream::badbit);

    file <<     "{ \"type\": \"FeatureCollection\",\n";
    file <<     "  \"features\": [\n";

    uint64_t count = 0;
    bool first_way = true;
    for(uint64_t i : export_ways){
        if(!first_way){file << ",";}
        file << "   {\"type\": \"Feature\",\n"
                "    \"geometry\": {\n"
                "       \"type\": \"LineString\",\n"
                "       \"coordinates\": [\n";

        bool first_node = true;
        for(uint64_t osmid : ways.at(i).refs){
            std::vector<Node>::iterator it = std::lower_bound(nodes.begin(), nodes.end(), Node{osmid}, compareSimpleNode);
            if(!first_node){file << ",";};
            file << "[" << it->longitude << "," << it->latitude << "]\n";
            first_node=false;
        }

        file << "   ]\n"
                "   },\n"
                "   \"properties\": {}\n"
                "   }";
        
        count++;
        if(count%1 == 0){file.flush();}
        first_way = false;
    }
    file <<     "]}\n" << std::endl;
    
}

/**
 * @brief save extracted coastline data to disk
 * 
 * @param save_string   path to location where data.save file will be created
 * @param nodes         vector of coastline nodes
 * @param ways          vector of coastline ways
 * @param export_ways   vector specifying index-ordering in which ways are exported
 */
void save_coastline_edges_to_file(std::string save_string, std::vector<Node> &nodes, std::vector<Way> &ways, std::vector<uint64_t> &export_ways){
    std::ofstream textfile;
    size_t lastindex = save_string.find_last_of(".");
    std::string text_file_name = save_string.substr(0, lastindex);
    text_file_name += ".save";
    textfile.open(text_file_name, std::ios::out | std::ios::trunc);
    textfile.exceptions(textfile.exceptions() | std::ios::failbit | std::ifstream::badbit);


    Node SourceNode;
    Node TargetNode;

    uint64_t tmp_osmid_source;
    uint64_t tmp_osmid_target;

    double source_longitude;
    double source_latitude;

    double target_longitude;
    double target_latitude;

    for (uint64_t way_id : export_ways){
        for(uint64_t ref_iterator = 0; ref_iterator < ways.at(way_id).refs.size() - 1; ref_iterator++){
            tmp_osmid_source = ways.at(way_id).refs.at(ref_iterator);
            tmp_osmid_target = ways.at(way_id).refs.at(ref_iterator+1);
            
            std::vector<Node>::iterator it_source = std::lower_bound(nodes.begin(), nodes.end(), Node{tmp_osmid_source}, compareSimpleNode);
            std::vector<Node>::iterator it_target = std::lower_bound(nodes.begin(), nodes.end(), Node{tmp_osmid_target}, compareSimpleNode);

            source_latitude = it_source->latitude;
            source_longitude = it_source->longitude;

            target_latitude = it_target->latitude;
            target_longitude = it_target->longitude;


            double plus_long_bound = 180.0;
            double minus_long_bound = -180.0;

            // cut edge if there is a wraparound on longitude
            if (!(signbit(source_longitude) != signbit(target_longitude) && abs(source_longitude > 100))){                    
                // do not cut the the edges
                textfile.write(reinterpret_cast<const char*>(&source_latitude), sizeof(source_latitude));
                textfile.write(reinterpret_cast<const char*>(&source_longitude), sizeof(source_longitude));
                textfile.write(reinterpret_cast<const char*>(&target_latitude), sizeof(target_latitude));
                textfile.write(reinterpret_cast<const char*>(&target_longitude), sizeof(target_longitude));

            }else{
                std::cout << "nocut:" << source_longitude << " " << target_longitude << std::endl;
                double intermediate_lat = 0.5 * (source_latitude + target_latitude);
                if(signbit(source_longitude)){
                    // first edge
                    textfile.write(reinterpret_cast<const char*>(&source_latitude), sizeof(source_latitude));
                    textfile.write(reinterpret_cast<const char*>(&source_longitude), sizeof(source_longitude));
                    textfile.write(reinterpret_cast<const char*>(&intermediate_lat), sizeof(source_latitude));
                    textfile.write(reinterpret_cast<const char*>(&plus_long_bound), sizeof(source_longitude));

                    // second edge
                    textfile.write(reinterpret_cast<const char*>(&intermediate_lat), sizeof(source_latitude));
                    textfile.write(reinterpret_cast<const char*>(&minus_long_bound), sizeof(source_longitude));
                    textfile.write(reinterpret_cast<const char*>(&target_latitude), sizeof(target_latitude));
                    textfile.write(reinterpret_cast<const char*>(&target_longitude), sizeof(target_longitude));
                }else{
                    // first edge
                    textfile.write(reinterpret_cast<const char*>(&source_latitude), sizeof(source_latitude));
                    textfile.write(reinterpret_cast<const char*>(&source_longitude), sizeof(source_longitude));
                    textfile.write(reinterpret_cast<const char*>(&intermediate_lat), sizeof(source_latitude));
                    textfile.write(reinterpret_cast<const char*>(&minus_long_bound), sizeof(source_longitude));

                    // second edge
                    textfile.write(reinterpret_cast<const char*>(&intermediate_lat), sizeof(source_latitude));
                    textfile.write(reinterpret_cast<const char*>(&plus_long_bound), sizeof(source_longitude));
                    textfile.write(reinterpret_cast<const char*>(&target_latitude), sizeof(target_latitude));
                    textfile.write(reinterpret_cast<const char*>(&target_longitude), sizeof(target_longitude));
                }

            }
        }

        textfile.flush();
    }

}


int main(int argc, char** argv) {
    std::cout << "Hello World\n";
    if(argc != 3) {
        std::cout << "Usage: " << argv[0] << " file_to_read.osm.pbf" << " " << "file_to save" << std::endl;
        return 1;
    }

    // start tracking time
    std::chrono::duration<double> main_program_timing;
    auto startMain = std::chrono::high_resolution_clock::now();


    std::vector<Node> nodes;
    std::vector<Way> ways;
    std::vector<uint64_t> export_ways;
    load_coastline_data(argv[1], nodes, ways, export_ways);
    std::cout << "Nodes vector takes up " << (nodes.size()*3.*4.)/1e6 << " MB\n";
    
    //save_coastline_to_geojson(argv[2], nodes, ways, export_ways);
    save_coastline_edges_to_file(argv[2], nodes, ways, export_ways);
    
    
    // stop tracking time
    auto endMain = std::chrono::high_resolution_clock::now();
    main_program_timing = endMain - startMain;
    std::cout << "Total .save export time: " << main_program_timing.count() << "seconds" << std::endl;

    return 0;

}