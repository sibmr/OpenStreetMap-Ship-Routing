#ifndef FILE_UTILS
#define FILE_UTILS

# include <fstream>
# include <vector>

template<typename T>
void writeSimpleValue(std::ofstream &outputFile, T &value){
    outputFile.write(reinterpret_cast<const char*>(&value),     sizeof(value));
}

template<typename T>
void readSimpleValue(std::ifstream &inputFile, T &value){
    inputFile.read(reinterpret_cast<char *>(&value),     sizeof(value));
}

template<typename T>
void writeVectorOfSimpleValue(std::ofstream &outputFile, std::vector<T> &valueVector){
    // write vector size
    uint64_t vectorSize = valueVector.size();
    writeSimpleValue(outputFile, vectorSize);
    outputFile.flush();
    // write vector values
    for(uint64_t i = 0; i < vectorSize; i++){
        writeSimpleValue(outputFile, valueVector.at(i));
        if(i % 20 == 0){
            outputFile.flush();
        }
    }
    outputFile.flush();
}

template<typename T>
void readVectorOfSimpleValue(std::ifstream &inputFile, std::vector<T> &valueVector){
    // read vector size
    uint64_t vectorSize;
    readSimpleValue(inputFile, vectorSize);
    valueVector.resize(vectorSize);
    // read vector values
    for(uint64_t i = 0; i < vectorSize; i++){
        inputFile.read(reinterpret_cast<char *>(&valueVector.at(i)), sizeof(valueVector.at(i)));
    }
}

#endif