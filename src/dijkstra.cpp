
# include <iostream>
# include <fstream>
# include <string>
# include <vector>


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


int main(int argc, char** argv) {
    AdjacencyArray adjArray;
    // load
    loadAdjacencyArray(adjArray, "data/worldGrid_1415_707.save");
    int counter_one = 0;
    int counter_zero = 0;
    for (int i = 0; i < adjArray.nodes.size(); i++){
        if(adjArray.nodes.at(i) == 0){
            counter_zero++;
        }else
        {
            counter_one++;
        }
    }
    std::cout << counter_zero << " " << counter_one << " " << counter_one + counter_zero << std::endl;
}