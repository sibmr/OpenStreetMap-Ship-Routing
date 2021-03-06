#include <iostream>
#include <fstream>
#include <chrono>
#include <numeric>
#include <algorithm>

#include <time.h>

#include "Dijkstra.cpp"
#include "Bidirectional_Dijkstra.cpp"
#include "A_star.cpp"
#include "CH_query.cpp"
#include "CH_Astar.cpp"



/**
 * @brief benchmark a path algorithm for a single query, multiple times
 * 
 * @param pathAlg       path algorithm object
 * @param adjArray      input AdjacencyArray
 * @param longStart     // 
 * @param latStart      // first point coordinates
 * @param longGoal      //
 * @param latGoal       // second point coordinates
 * @param numAvg        // number of query executions to average
 */
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
 * @brief generate a random double value between (doubleMin, doubleMax)
 * 
 * @param doubleMin     lower bound for random value
 * @param doubleMax     upper bound for random value
 * @return double       random value
 */
double doubleRand(double doubleMin, double doubleMax){
    double d = ((double) rand() / (double) RAND_MAX);
    return doubleMin + d * (doubleMax - doubleMin);
}

/**
 * @brief compute the average of the values in the vector
 * 
 * @param data          vector of values to use in the computation
 * @return uint64_t     mean of the values
 */
uint64_t computeAverage(std::vector<uint64_t> &data){
    return std::accumulate(data.begin(), data.end(), 0) / data.size();
}

/**
 * @brief compute the standart deviation of the values in the vector
 * 
 * @param data          vector of values to use in the computation
 * @return uint64_t     standart deviation of the values
 */
uint64_t computeStdDev(std::vector<uint64_t> &data, uint64_t dataAvg){
    uint64_t stddev = 0;
    for(uint64_t value : data){
        uint64_t diff = value-dataAvg;
        stddev += diff*diff;
    }
    stddev = sqrt(stddev)/data.size(); 
    return stddev;
}

/**
 * @brief Test two PathAlgorithms to be identical and test how long each takes
 * in addition, execution duration and number of heap elements popped mean and standart deviation are printed
 * 
 * @param pathAlg           first Path Algorith to test
 * @param pathAlg2          second Path Algorith to test
 * @param adjArray          input AdjacencyArray
 * @param longStart         lower longitude limit
 * @param latStart          lower latitude limit
 * @param longGoal          upper longitude limit
 * @param latGoal           upper latitude limit
 * @param numSamplePoints   number path calculations to compare
 */
void comparePathAlgorithms(PathAlgorithm &pathAlg, PathAlgorithm &pathAlg2,  AdjacencyArray &adjArray,
    double longStart, double latStart, double longGoal, double latGoal, int numSamplePoints)
{
    srand(2);

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

    uint64_t sNode = UINT64_MAX;
    uint64_t tNode = UINT64_MAX;

    for (int sampleId = 0; sampleId < numSamplePoints; sampleId++){

        std::vector<uint64_t> timings;
        bool found_two_nodes;

        // get two random samplepoints and check if there exists a path
        do{
            currLongStart = doubleRand(longStart, longGoal);
            currLongEnd = doubleRand(longStart, longGoal);
            currLatStart = doubleRand(latStart, latGoal);
            currLatEnd = doubleRand(latStart, latGoal);

            sNode = longLatToNodeId(adjArray, currLongStart, currLatStart);
            tNode = longLatToNodeId(adjArray, currLongEnd, currLatEnd);
            pathAlg.reset();        
            temp_result = pathAlg.calculateDist(sNode, tNode);

            if(sNode >= adjArray.width*adjArray.height || tNode >= adjArray.width*adjArray.height){
                continue;
            }

        }while(temp_result == UINT64_MAX);


        coordinates.push_back(std::array<double,4> {currLongStart, currLatStart, currLongEnd, currLatEnd});

        // benchmark first algorithm
        pathAlg.reset();
        
        startQuery = std::chrono::high_resolution_clock::now();
        
        temp_result = pathAlg.calculateDist(sNode, tNode);

        endQuery = std::chrono::high_resolution_clock::now();
        
        resultPathAlgOne.push_back(temp_result);
        timingPathAlgOne.push_back(std::chrono::duration_cast<std::chrono::microseconds>(endQuery - startQuery).count());
        numNodesPoppedPathAlgOne.push_back(pathAlg.getNumNodesPopped());

        // benchmark second algorithm
        pathAlg2.reset();
        
        startQuery = std::chrono::high_resolution_clock::now();
                
        temp_result = pathAlg2.calculateDist(sNode, tNode);

        endQuery = std::chrono::high_resolution_clock::now();
        
        resultPathAlgTwo.push_back(temp_result);
        timingPathAlgTwo.push_back(std::chrono::duration_cast<std::chrono::microseconds>(endQuery - startQuery).count());
        numNodesPoppedPathAlgTwo.push_back(pathAlg2.getNumNodesPopped());

    }
 
    uint64_t avgOne = computeAverage(timingPathAlgOne);
    uint64_t stddevOne = computeStdDev(timingPathAlgOne, avgOne);
    uint64_t avgPoppedOne = computeAverage(numNodesPoppedPathAlgOne);

     
    uint64_t avgTwo = computeAverage(timingPathAlgTwo);
    uint64_t stddevTwo = computeStdDev(timingPathAlgTwo, avgTwo);
    uint64_t avgPoppedTwo = computeAverage(numNodesPoppedPathAlgTwo);

    for(int i = 0; i < timingPathAlgOne.size(); i++){
        if(resultPathAlgOne.at(i) != resultPathAlgTwo.at(i)){
            std::cout << "(" << coordinates.at(i)[0] << "," << coordinates.at(i)[1] << ") (" << coordinates.at(i)[2] << "," << coordinates.at(i)[3] << ")\t has wrong results - " << resultPathAlgOne.at(i) << " " << resultPathAlgTwo.at(i) << std::endl;
        }else{
            std::cout << "(" << coordinates.at(i)[0] << "," << coordinates.at(i)[1] << ") (" << coordinates.at(i)[2] << "," << coordinates.at(i)[3] << ")" <<  std::endl;
        }
    }

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
    std::string inputFileName;
    std::string filenameCH;
    std::string pathAlg1;
    std::string pathAlg2;
    // check input parameter
    if(argc < 1 || argc > 5){
        std::cout << "Usage: " << argv[0] << " file_to_read.graph ch_file_to_read.graph pathAlg1Name pathAlg2Name" << std::endl;
        return 1;
    }
    if(argc == 4){
        inputFileName = std::string(argv[1]);
        filenameCH =std::string(argv[2]);
        pathAlg1 = "dijkstra";
        pathAlg2 = argv[3];
    }
    else if(argc > 4){
        inputFileName = std::string(argv[1]);
        filenameCH =std::string(argv[2]);
        pathAlg1 = argv[3];
        pathAlg2 = argv[4];
    }
    else{
        inputFileName = "data/planet.graph_2";
        filenameCH = "data/CHAdjArray_54.graph_2";
        pathAlg1 = "dijkstra";
        pathAlg2 = "bidijkstra";
    }

    AdjacencyArray adjArray(inputFileName);
    AdjacencyArray adjArrayCH(filenameCH);
    
    Dijkstra::Dijkstra dijk (adjArray);
    BidirectionalDijkstra::BidirectionalDijkstra bidijk (adjArray);
    A_star::A_star astar (adjArray);
    A_star::A_star_rectangular astar_rect (adjArray);
    CH_query::CH_query chquery (adjArrayCH);
    CH_Astar::A_star chqueryStar (adjArrayCH);
    CH_Astar::A_star_rectangular chqueryRectStar (adjArrayCH);

    std::map<std::string, PathAlgorithm*> algs{
        {"dijkstra", &dijk},
        {"bidijkstra", &bidijk},
        {"astar", &astar},
        {"rectastar", &astar_rect},
        {"chquery", &chquery},
        {"chastar", &chqueryStar},
        {"chrectastar", &chqueryRectStar}
    };
    
    if(algs.find(pathAlg1)==algs.end()){
        std::cout << "path algorithm with name " << pathAlg1 << " does not exist\n";
        return 0;
    }

    if(algs.find(pathAlg2)==algs.end()){
        std::cout << "path algorithm with name " << pathAlg2 << " does not exist\n";
        return 0;
    }
    
    std::cout << inputFileName << "  " << filenameCH << "  " << pathAlg1 << "  " << pathAlg2 << "\n";

    PathAlgorithm &pa = *algs.at(pathAlg1);
    PathAlgorithm &pa_one = *algs.at(pathAlg2);
    comparePathAlgorithms(pa, pa_one, adjArray, -180, -85, 180, 85, 100);
}
