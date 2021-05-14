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
    textfile.open("data/worldGrid_400_200.save", std::ios::in);

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
    std::cout << gridPoint << std::endl;
}

void fillAdjacencyArray(GridData &data, AdjacencyArray &array){
    for(uint64_t i = 0; i<data.width; ++i)
    for(uint64_t j = 0; j<data.height; ++j)
    {
        uint64_t x = i*data.height + j;
        
    }
}

int main(int argc, char** argv) {
    GridData dat;
    AdjacencyArray adjArray;
    loadGridPoints(dat);
    fillAdjacencyArray(dat, adjArray);
    // save
    std::cout << dat.height << std::endl;
}