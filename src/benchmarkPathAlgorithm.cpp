#include <iostream>
#include <fstream>
#include <chrono>
#include <numeric>
#include <algorithm>

#include <time.h>

#include "Dijkstra.cpp"
#include "BiDirectDijkstra_marcel.cpp"
#include "ChQuery_marcel.cpp"

double doubleRand(double doubleMin, double doubleMax){
    double d = ((double) rand() / (double) RAND_MAX);
    return doubleMin + d * (doubleMax - doubleMin);
}

uint64_t computeAverage(std::vector<uint64_t> &data){
    if(data.size() > 0){
        return std::accumulate(data.begin(), data.end(), 0) / data.size();
    }else{
        return 1;
    }
}

uint64_t computeStdDev(std::vector<uint64_t> &data, uint64_t dataAvg){
    if(data.size() > 0){
        uint64_t stddev = 0;
        for(uint64_t value : data){
            uint64_t diff = value-dataAvg;
            stddev += diff*diff;
        }
        stddev = sqrt(stddev)/data.size(); 
        return stddev;
    }
    else{
        return 1;
    }
}

void benchmarkDijkstra(PathAlgorithm &pathAlg, AdjacencyArray &adjArray,
    double longStart, double latStart, double longGoal, double latGoal, int numAvg)
{
    std::vector<uint64_t> timings;
    uint64_t queryTiming;

    uint64_t sNode = longLatToNodeId(adjArray, longStart, latStart);
    uint64_t tNode = longLatToNodeId(adjArray, longGoal, latGoal);

    pathAlg.reset();
    for(int i = 0; i<numAvg; ++i){

        auto startQuery = std::chrono::high_resolution_clock::now();

        std::cout << pathAlg.calculateDist(sNode, tNode) << std::endl;

        
        auto endQuery = std::chrono::high_resolution_clock::now();
        queryTiming = std::chrono::duration_cast<std::chrono::microseconds>(endQuery - startQuery).count();
        std::cout << "Run " << i+1 << ": " << queryTiming << " us\n";
        timings.push_back(queryTiming);
    }

    uint64_t stddev = 0; 
    uint64_t avg = std::accumulate(timings.begin(), timings.end(), 0) / numAvg;
    
    stddev = 0;
    for(uint64_t duration : timings){
        uint64_t diff = duration-avg;
        stddev += diff*diff;
    }
    stddev = sqrt(stddev)/numAvg;

    std::cout << "Duration: " << avg << " +/- " << stddev << "us\n";
}



/**
 * @brief Test two PathAlgorithms to be identical and test how long each takes
 * 
 * @param pathAlg           first Path Algorith to test
 * @param adjArray      
 * @param longStart 
 * @param latStart 
 * @param longGoal 
 * @param latGoal 
 * @param numSamplePoints 
 */
void debugDijkstra(PathAlgorithm &pathAlg,  AdjacencyArray &adjArray,
    double longStart, double latStart, double longGoal, double latGoal)
{
    uint64_t sNode = longLatToNodeId(adjArray, longStart, latStart);
    uint64_t tNode = longLatToNodeId(adjArray, longGoal, latGoal);

    pathAlg.reset();

    std::cout << pathAlg.calculateDist(sNode, tNode) << std::endl;
}

/**
 * @brief Test two PathAlgorithms to be identical and test how long each takes
 * 
 * @param pathAlg           first Path Algorith to test
 * @param pathAlg2          second Path Algorith to test
 * @param adjArray      
 * @param longStart 
 * @param latStart 
 * @param longGoal 
 * @param latGoal 
 * @param numSamplePoints 
 */
void testDijkstra(PathAlgorithm &pathAlg, PathAlgorithm &pathAlg2,  AdjacencyArray &adjArray,
    double longStart, double latStart, double longGoal, double latGoal, int numSamplePoints)
{
    srand(time(NULL));

    std::vector<uint64_t> resultPathAlgOne;
    std::vector<uint64_t> resultPathAlgTwo;

    std::vector<uint64_t> timingPathAlgOne;
    std::vector<uint64_t> timingPathAlgTwo;

    std::vector<uint64_t> numNodesPoppedPathAlgOne;
    std::vector<uint64_t> numNodesPoppedPathAlgTwo;    

    std::vector<std::array<double,4>> coordinates;

    uint64_t queryTiming;
    uint64_t temp_result;
     
    uint64_t accumulatedAvgOne = 0;
    uint64_t accumulatedAvgTwo = 0;

    double currLongStart, currLatStart, currLongEnd, currLatEnd;
    auto startQuery = std::chrono::high_resolution_clock::now();
    auto endQuery = std::chrono::high_resolution_clock::now();
    
    int currSample = 0;

    while(currSample < numSamplePoints){

        std::vector<uint64_t> timings;
        
        // get two random samplepoints and check if there exists a path
        currLongStart = doubleRand(longStart, longGoal);
        currLongEnd = doubleRand(longStart, longGoal);
        currLatStart = doubleRand(latStart, latGoal);
        currLatEnd = doubleRand(latStart, latGoal);

        // get start and end node 
        uint64_t sNode = longLatToNodeId(adjArray, currLongStart, currLatStart);
        uint64_t tNode = longLatToNodeId(adjArray, currLongEnd, currLatEnd);

        // check with first algorithm if there is a path between start and endnode
        pathAlg.reset();        
        temp_result = pathAlg.calculateDist(sNode, tNode);
        if(temp_result == UINT64_MAX){
            continue;
        }

        // there exists a path between the two nodes
        currSample++;
        std::cout << "run: " << currSample << " of " << numSamplePoints << std::endl;



        coordinates.push_back(std::array<double,4> {currLongStart, currLatStart, currLongEnd, currLatEnd});

        // start timing of first algorithm
        startQuery = std::chrono::high_resolution_clock::now();

        pathAlg.reset();        
        temp_result = pathAlg.calculateDist(sNode, tNode);

        endQuery = std::chrono::high_resolution_clock::now();
        resultPathAlgOne.push_back(temp_result);
        timingPathAlgOne.push_back(std::chrono::duration_cast<std::chrono::microseconds>(endQuery - startQuery).count());
        numNodesPoppedPathAlgOne.push_back(pathAlg.getNumNodesPopped());


        // start timing of second algorithm
        startQuery = std::chrono::high_resolution_clock::now();

        pathAlg2.reset();        
        temp_result = pathAlg2.calculateDist(sNode, tNode);
        endQuery = std::chrono::high_resolution_clock::now();
        resultPathAlgTwo.push_back(temp_result);
        timingPathAlgTwo.push_back(std::chrono::duration_cast<std::chrono::microseconds>(endQuery - startQuery).count());
        numNodesPoppedPathAlgTwo.push_back(pathAlg2.getNumNodesPopped());


    }
    
    // all queries have been completed output results


    uint64_t avgOne = computeAverage(timingPathAlgOne);
    uint64_t stddevOne = computeStdDev(timingPathAlgOne, avgOne);
    uint64_t avgPoppedOne = computeAverage(numNodesPoppedPathAlgOne);

     
    uint64_t avgTwo = computeAverage(timingPathAlgTwo);
    uint64_t stddevTwo = computeStdDev(timingPathAlgTwo, avgTwo);
    uint64_t avgPoppedTwo = computeAverage(numNodesPoppedPathAlgTwo);

    

    // output the samplepoints and show diff if there is any
    for(int i = 0; i < timingPathAlgOne.size(); i++){
        if(resultPathAlgOne.at(i) != resultPathAlgTwo.at(i)){
            std::cout << "(" << coordinates.at(i)[0] << "," << coordinates.at(i)[1] << ") (" << coordinates.at(i)[2] << "," << coordinates.at(i)[3] << ")\t has wrong results - " << resultPathAlgOne.at(i) << " " << resultPathAlgTwo.at(i) << std::endl;
        }else{
            std::cout << "(" << coordinates.at(i)[0] << "," << coordinates.at(i)[1] << ") (" << coordinates.at(i)[2] << "," << coordinates.at(i)[3] << ")" <<  std::endl;
        }
    }

    // output timing differences
    std::cout << "First alg has from  " << timingPathAlgOne.size() << " Queries an Average of " << avgOne << "us and stddev of: "<< stddevOne <<  "us" << std::endl;
    std::cout << "Second alg has from " << timingPathAlgTwo.size() << " Queries an Average of " << avgTwo << "us and stddev of: "<< stddevTwo << "us" <<  std::endl;
    std::cout << "On average the first algorithm takes  " <<  ((double)((avgOne * 10000)/(avgTwo)) / 10000) << " times longer" << std::endl;
    std::cout << "On average the second algorithm takes " <<  ((double)((avgTwo * 10000)/(avgOne)) / 10000) << " times longer" << std::endl;
    std::cout << "First alg has from  " << numNodesPoppedPathAlgOne.size() << " Queries an Average of " << avgPoppedOne << " nodes\n";
    std::cout << "Second alg has from " << numNodesPoppedPathAlgTwo.size() << " Queries an Average of " << avgPoppedTwo << " nodes\n";
    std::cout << "On average the first algorithm pops  " <<  ((double)((avgPoppedOne * 10000)/(avgPoppedTwo)) / 10000) << " times more nodes" << std::endl;
    std::cout << "On average the second algorithm pops " <<  ((double)((avgPoppedTwo * 10000)/(avgPoppedOne)) / 10000) << " times more nodes" << std::endl;


}

int main(int argc, char** argv) {
    std::cout << "First graph uses dijkstra second graph uses contraction hierachies" << std::endl;
    std::string inputFileName = "data/planet_2.graph_2";
    std::string inputFileName2 = "data/planet_2.graph_3_920p_invT";
    int numberOfTries = 100; 

    bool outputHelp = false;
    if(argc < 1 || argc > 4){
        outputHelp = true;
    }else if(argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")){
        outputHelp = true;
    }


    // check input parameter
    if(outputHelp){
        std::cout << "Usage: " << argv[0]  << std::endl;
        std::cout << "Usage: " << argv[0] << " numberOfTries " << std::endl;
        std::cout << "Usage: " << argv[0] << " numberOfTries " << " file_to.graph" << std::endl;
        std::cout << "Usage: " << argv[0] << " numberOfTries " << " file_to_first.graph" << " " << "file_to_second.graph"<< std::endl;
        return 1;
    }
    if(argc > 1){
        numberOfTries = atoi(argv[1]);
    }
    if(argc > 2){
        inputFileName = std::string(argv[2]);
    }
    if(argc > 3){
        inputFileName2 = std::string(argv[3]);
    }


    AdjacencyArray adjArray(inputFileName);
    AdjacencyArray chArray(inputFileName2);
    SecondDijkstra dijkstra (adjArray);
    //BiDirectDijkstra::BiDirectDijkstra sd2 (adjArray);
    {
        ChQuery::ChQuery chQuery (chArray);

        PathAlgorithm &pa1 = dijkstra;
        PathAlgorithm &pa2 = chQuery;

        std::cout << "First array has " << adjArray.edges.size() << " edges. Second array has " <<  chArray.edges.size() << " edges." << std::endl;
        testDijkstra(pa1, pa2, adjArray, -180, -85, 180, 85, numberOfTries);
    }
}
