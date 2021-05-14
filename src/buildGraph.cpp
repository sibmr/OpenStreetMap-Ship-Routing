# include <iostream>
# include <fstream>
# include <string>
# include <vector>

struct GridData {
    double longLow, latLow, longHigh, latHigh;
    uint64_t width;
    uint64_t height;
    std::vector<bool> gridDataList; 
};

struct AdjacencyArray {
    double longLow, latLow, longHigh, latHigh;
    uint64_t width, height;
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> edges;
    std::vector<bool> nodes;
};

void loadGridPoints(GridData &gridData){
    std::ifstream textfile;
    textfile.open("data/worldGrid_1415_707.save", std::ios::in);

    double longLow, latLow, longHigh, latHigh;
    uint64_t width, height;
    bool gridPoint;

    textfile.read(reinterpret_cast<char *>(&longLow),   sizeof(longLow));
    textfile.read(reinterpret_cast<char *>(&latLow),    sizeof(latLow));
    textfile.read(reinterpret_cast<char *>(&longHigh),  sizeof(longHigh));
    textfile.read(reinterpret_cast<char *>(&latHigh),   sizeof(latHigh));
    textfile.read(reinterpret_cast<char *>(&width), sizeof(width));
    textfile.read(reinterpret_cast<char *>(&height), sizeof(height));
    
    gridData.longLow    = longLow;
    gridData.latLow     = latLow;
    gridData.longHigh   = longHigh;
    gridData.latHigh    = latHigh;


    gridData.width = width;
    gridData.height = height;

    uint64_t count = 0;
    uint64_t i = 0;
    uint64_t j = 0;
    while(!textfile.eof()){
        textfile.read(reinterpret_cast<char *>(&gridPoint), sizeof(gridPoint));
        i = count/height;
        j = count%height;
        count += 1;

        gridData.gridDataList.push_back(gridPoint);
        
    }


    std::cout << width << std::endl;
    std::cout << height << std::endl;
    std::cout << gridData.gridDataList.at(0) << std::endl;
    std::cout << gridData.gridDataList.at(1) << std::endl;
    std::cout << gridData.gridDataList.at(2) << std::endl;
    std::cout << latLow << std::endl;
}

void fillAdjacencyArray(GridData &data, AdjacencyArray &array){
    
    array.longLow   = data.longLow;
    array.latLow    = data.latLow;
    array.longHigh  = data.longHigh;
    array.latHigh   = data.latHigh;
    array.width     = data.width;
    array.height    = data.height;
    
    array.nodes = data.gridDataList;

    uint64_t current_offset = 0;
    uint64_t offset_step    = 0;
    array.offsets.push_back(0);
    for(uint64_t i = 0; i<data.width; ++i)
    for(uint64_t j = 0; j<data.height; ++j)
    {
        uint64_t x      = i*data.height + j;
        
        if(data.gridDataList.at(x)){
            uint64_t left   = i*data.height + ((j-1) % data.height);
            uint64_t right  = i*data.height + ((j+1) % data.height);
            uint64_t down   = (i-1)*data.height + j;
            uint64_t up     = (i+1)*data.height + j;

            if(data.gridDataList.at(left)){
                offset_step += 1;
                array.edges.push_back(left);
            }
            if(data.gridDataList.at(right)){
                offset_step += 1;
                array.edges.push_back(right);
            }
            if(i>0              && data.gridDataList.at(down)){
                offset_step += 1;
                array.edges.push_back(down);
            }
            if(i<data.width-1   && data.gridDataList.at(up)){
                offset_step += 1;
                array.edges.push_back(up);
            }
        }

        current_offset += offset_step;
        offset_step = 0;
        array.offsets.push_back(current_offset);
    }
}

void testLoadFill(GridData &dat, AdjacencyArray &adjArray){
    std::cout << "longLow:\n" << adjArray.longLow   << "\n";
    std::cout << "latLow:\n" << adjArray.latLow     << "\n";
    std::cout << "longHigh:\n" << adjArray.longHigh << "\n";
    std::cout << "latHigh:\n" << adjArray.latHigh   << "\n";
    std::cout << "width:\n" << adjArray.width       << "\n";
    std::cout << "height:\n" << adjArray.height     << "\n";
    
    std::cout << "Offsets:\n";
    for(int i = 0; i<10; ++i){
        std::cout << adjArray.offsets.at(i) << "  ";
    }
    std::cout << "\nEdges:\n";
    for(int i = 0; i<10; ++i){
        std::cout << adjArray.edges.at(i) << "  ";
    }
    std::cout << "\nData Nodes:\n";
    for(int i = 0; i<10; ++i){
        std::cout << dat.gridDataList.at(i) << "  ";
    }
    std::cout << "\nArray Nodes:\n";
    for(int i = 0; i<10; ++i){
        std::cout << adjArray.nodes.at(i) << "  ";
    }
    std::cout << "\n";
}

    // double longLow, latLow, longHigh, latHigh;
    // uint64_t width, height;
    // std::vector<uint64_t> offsets;
    // std::vector<uint64_t> edges;
    // std::vector<bool> nodes;
/*
* output file
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
void saveAdjacencyArray(AdjacencyArray &array, std::string path){
    std::ofstream adjacency_output_file;

    size_t lastindex = path.find_last_of(".");
    std::string adjacency_file_name = path.substr(0, lastindex);
    adjacency_file_name += "_adjacencyarray.save";

    adjacency_output_file.open(adjacency_file_name, std::ios::out | std::ios::trunc);
    adjacency_output_file.exceptions(adjacency_output_file.exceptions() | std::ios::failbit | std::ifstream::badbit);

    // write globe size
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.longLow),     sizeof(array.longLow));
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.latLow),      sizeof(array.latLow));
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.longHigh),    sizeof(array.longHigh));
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.latHigh),     sizeof(array.latHigh));
    

    // save number of nodes per direction (width, height)
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.width),     sizeof(array.width));
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.height),     sizeof(array.height));

    // save offset vectro size
    uint64_t offset_size = array.offsets.size();
    adjacency_output_file.write(reinterpret_cast<const char*>(&offset_size),     sizeof(offset_size));

    adjacency_output_file.flush();

    for(uint64_t i = 0; i < offset_size; i++){
        adjacency_output_file.write(reinterpret_cast<const char*>(&array.offsets.at(i)),     sizeof(array.offsets.at(i)));
        if(i % 20 == 0){
            adjacency_output_file.flush();
        }
    }
    adjacency_output_file.flush();


    uint64_t edges_size = array.edges.size();
    adjacency_output_file.write(reinterpret_cast<const char*>(&edges_size),     sizeof(edges_size));

    adjacency_output_file.flush();
    for(uint64_t i = 0; i < edges_size; i++){
        adjacency_output_file.write(reinterpret_cast<const char*>(&array.edges.at(i)), sizeof(array.edges.at(i)));
        if(i % 20 == 0){
            adjacency_output_file.flush();
        }
    }
    adjacency_output_file.flush();

    uint64_t nodes_size = array.nodes.size();
    adjacency_output_file.write(reinterpret_cast<const char*>(&nodes_size), sizeof(nodes_size));

    adjacency_output_file.flush();
    std::copy(array.nodes.begin(), array.nodes.end(), std::ostreambuf_iterator<char>(adjacency_output_file));

    adjacency_output_file.flush();
    adjacency_output_file.close();
}

int main(int argc, char** argv) {
    GridData dat;
    AdjacencyArray adjArray;
    loadGridPoints(dat);
    fillAdjacencyArray(dat, adjArray);
    // save
    std::cout << dat.height << std::endl;
    testLoadFill(dat, adjArray);
    saveAdjacencyArray(adjArray, "data/worldGrid_1415_707.save");
}