# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <array>
# include <math.h>

# include "shortestPathUtils.cpp"
# include "ContractionHierachies_marcel.cpp"


struct GridData {
    double longLow, latLow, longHigh, latHigh;
    uint64_t width;
    uint64_t height;
    std::vector<bool> gridDataList; 
};

//struct AdjacencyArray {
//    double longLow, latLow, longHigh, latHigh;
//    uint64_t width, height;
//    std::vector<uint64_t> offsets;
//    std::vector<uint64_t> edges;
//    std::vector<uint64_t> distances;
//    std::vector<uint16_t> rank;
//    std::vector<bool> nodes;
//};

/**
 * @brief read grid data from disk and store it in struct
 * 
 * @param gridData  out: struct containg the grid data (implicit 2D bool array indicating land/ocean for each point)
 */
void loadGridPoints(GridData &gridData, std::string path){
    std::ifstream textfile;
    textfile.open(path, std::ios::in);

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

    while(!textfile.eof()){
        textfile.read(reinterpret_cast<char *>(&gridPoint), sizeof(gridPoint));
        gridData.gridDataList.push_back(gridPoint);
    }


    std::cout << width << std::endl;
    std::cout << height << std::endl;
    std::cout << gridData.gridDataList.at(0) << std::endl;
    std::cout << gridData.gridDataList.at(1) << std::endl;
    std::cout << gridData.gridDataList.at(2) << std::endl;
    std::cout << latLow << std::endl;
}

/**
 * @brief generates adjacency array from grid data (nodes: grid points, edges: between neighboring nodes with land==false)
 * 
 * @param data  grid data as basis for the adjacency array
 * @param array out: adjacency array corresponding to input grid data
 */
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

    std::cout << "nodes: "  << data.gridDataList.size() << std::endl;


    for(uint64_t i = 0; i<data.width; ++i)
    for(uint64_t j = 0; j<data.height; ++j)
    {
        uint64_t x      = i*data.height + j;

        
        if(!data.gridDataList.at(x)){               // if its water
            uint64_t down   = i*data.height + (j-1);
            uint64_t up     = i*data.height + (j+1);
            uint64_t left   = ((data.width+i-1)%data.width)*data.height + j;
            uint64_t right  = ((i+1) % data.width)*data.height + j;

            if(!data.gridDataList.at(left)){
                offset_step += 1;
                array.edges.push_back(left);
                array.distances.push_back(nodeDistance(array, x, left));
            }
            if(!data.gridDataList.at(right)){
                offset_step += 1;
                array.edges.push_back(right);
                array.distances.push_back(nodeDistance(array, x, right));
            }
            if(j>0              && !data.gridDataList.at(down)){
                offset_step += 1;
                array.edges.push_back(down);
                array.distances.push_back(nodeDistance(array, x, down));
            }
            if(j<data.height-1   && !data.gridDataList.at(up)){
                offset_step += 1;
                array.edges.push_back(up);
                array.distances.push_back(nodeDistance(array, x, up));
            }
        }


        current_offset += offset_step;
        offset_step = 0;
        array.offsets.push_back(current_offset);
        array.rank.push_back(0);
    }
    std::cout << array.distances.size() << std::endl;
    std::cout << array.edges.size() << std::endl;
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

/**
 * @brief save adjacency array to disk
 * 
 * @param array AdjacencyArray struct that is stored
 * @param path  path to storage location
 * 
 * output file format
 * longLow       - double
 * latLow        - double
 * longHigh      - double
 * latHigh       - double
 * width         - uint64_t
 * height        - uint64_t
 * offset_size   - uint64_t
 * offsets       - uint64_t     (offset_size many)
 * edges_size    - uint64_t
 * edges         - uint64_t     (edges_size many)
 * nodes_size    - uint64_t
 * nodes          - bool        (nodes_size many)
 */
void saveAdjacencyArray(AdjacencyArray &array, std::string path){
    std::ofstream adjacency_output_file;

    adjacency_output_file.open(path, std::ios::out | std::ios::trunc);
    adjacency_output_file.exceptions(adjacency_output_file.exceptions() | std::ios::failbit | std::ifstream::badbit);

    // write globe size
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.longLow),     sizeof(array.longLow));
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.latLow),      sizeof(array.latLow));
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.longHigh),    sizeof(array.longHigh));
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.latHigh),     sizeof(array.latHigh));
    

    // save number of nodes per direction (width, height)
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.width),     sizeof(array.width));
    adjacency_output_file.write(reinterpret_cast<const char*>(&array.height),     sizeof(array.height));

    // save offset vector size
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


    // save edges
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

    // save distances
    uint64_t distances_size = array.distances.size();
    adjacency_output_file.write(reinterpret_cast<const char*>(&distances_size),     sizeof(distances_size));

    adjacency_output_file.flush();
    for(uint64_t i = 0; i < distances_size; i++){
        adjacency_output_file.write(reinterpret_cast<const char*>(&array.distances.at(i)), sizeof(array.distances.at(i)));
        if(i % 20 == 0){
            adjacency_output_file.flush();
        }
    }
    adjacency_output_file.flush();

    // save rank
    uint64_t rank_size = array.rank.size();
    adjacency_output_file.write(reinterpret_cast<const char*>(&rank_size),     sizeof(rank_size));

    adjacency_output_file.flush();
    for(uint64_t i = 0; i < rank_size; i++){
        adjacency_output_file.write(reinterpret_cast<const char*>(&array.rank.at(i)), sizeof(array.rank.at(i)));
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
    std::string inputFileName;
    std::string outputFileName;

    int inputFileId = -1;
    int outputFileId = -1;

    bool verboseOutput = false;
    bool ch = false;
    // always go in contract mode
    ch = true;

    // check input parameter
    if(argc < 1 || argc > 4){
        std::cout << "Usage: " << argv[0] << " file_to_read.grid" << " " << "file_to.graph" << std::endl;
        return 1;
    }
    if(argc > 1){
        if(std::string(argv[1]) == "-t"){
            verboseOutput = true;
            if(argc > 2){
                inputFileId = 2;
            }
            if(argc > 3){
                outputFileId = 3;
            }
        }else if(std::string(argv[1]) == "-ch"){
            ch = true;
            if(argc > 2){
                inputFileId = 2;
            }
            if(argc > 3){
                outputFileId = 3;
            }
        }
        else{
            if(argc > 1){
                inputFileId = 1;
            }if(argc > 2){
                outputFileId = 2;
            }
        }
    }

    // generate input filename
    if(inputFileId == -1){
        inputFileName = "data/planet.grid";
        inputFileName = "data/planet_2.grid";
        //inputFileName = "data/antarctica_100K.grid";
        std::cout << "no input file given assume " <<  inputFileName << std::endl;
    }else{
        inputFileName = std::string(argv[inputFileId]);
    }

    // generate output filename
    std::string tmp_oFile;
    if(outputFileId == -1){
        tmp_oFile = inputFileName;
    }else{
        tmp_oFile = argv[outputFileId];
    }
    size_t lastindex = tmp_oFile.find_last_of(".");
    outputFileName = tmp_oFile.substr(0, lastindex);
    // distinguish between graphs
    outputFileName += ".graph_3";

    GridData dat;
    AdjacencyArray adjArray;

    loadGridPoints(dat, inputFileName);
    fillAdjacencyArray(dat, adjArray);

    if(ch){
        contract(adjArray, 0.9);
    }



    // save
    std::cout << dat.height << std::endl;
    if(verboseOutput){
        testLoadFill(dat, adjArray);
    }
    saveAdjacencyArray(adjArray, outputFileName);
}