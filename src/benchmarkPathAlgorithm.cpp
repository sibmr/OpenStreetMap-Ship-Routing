#include <iostream>
#include <chrono>
#include <numeric>
#include <algorithm>

#include "Dijkstra.cpp"

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

double doubleRand(double doubleMin, double doubleMax){
    double d = ((double) rand() / (double) RAND_MAX);
    return doubleMin + d * (doubleMax - doubleMin);
}

void testDijkstra(PathAlgorithm &pathAlg, AdjacencyArray &adjArray,
    double longStart, double latStart, double longGoal, double latGoal, int numAvg, int numSamplePoints)
{
    std::vector<uint64_t> timings;
    uint64_t queryTiming;
    uint64_t temp_result;

    uint64_t accumulatedAvgOne = 0;
    uint64_t accumulatedAvgTwo = 0;

    double currLongStart, currLatStart, currLongEnd, currLatEnd;


    for (int sampleId = 0; sampleId < numSamplePoints; sampleId++){
        
        currLongStart = doubleRand(longStart, longGoal);
        currLongEnd = doubleRand(longStart, longGoal);
        currLatStart = doubleRand(latStart, latGoal);
        currLatEnd = doubleRand(latStart, latGoal);

        std::cout << currLongStart << " " << currLatStart << " | " << currLongEnd << " " << currLatEnd << std::endl;

        uint64_t sNode = longLatToNodeId(adjArray, currLongStart, currLatStart);
        uint64_t tNode = longLatToNodeId(adjArray, currLongEnd, currLatEnd);

        uint64_t result1 = 0;

        for(int i = 0; i<numAvg; ++i){

            auto startQuery = std::chrono::high_resolution_clock::now();

            pathAlg.reset();        
            temp_result = pathAlg.calculateDist(sNode, tNode);
            if(temp_result == UINT64_MAX){
                continue;
            }
            result1 += temp_result;

            auto endQuery = std::chrono::high_resolution_clock::now();
            queryTiming = std::chrono::duration_cast<std::chrono::microseconds>(endQuery - startQuery).count();
            //std::cout << "Run " << i+1 << ": " << queryTiming << " us\n";
            timings.push_back(queryTiming);
        }

        uint64_t stddev = 0; 
        uint64_t avg = std::accumulate(timings.begin(), timings.end(), 0) / numAvg;

        accumulatedAvgOne += avg;
    
        stddev = 0;
        for(uint64_t duration : timings){
            uint64_t diff = duration-avg;
            stddev += diff*diff;
        }
        stddev = sqrt(stddev)/numAvg;

        std::cout << "Duration: " << avg << " +/- " << stddev << "us\n";

        std::cout << "--------------------------------------" << std::endl;

        std::vector<uint64_t> timings2;
        uint64_t queryTiming2;

        uint64_t result2 = 0;

        for(int i = 0; i<numAvg; ++i){

            auto startQuery = std::chrono::high_resolution_clock::now();

            pathAlg.reset();
            temp_result = pathAlg.calculateDistSavedEdges(sNode, tNode);
            if(temp_result == UINT64_MAX){
                continue;
            }
            result2 += temp_result;


            auto endQuery = std::chrono::high_resolution_clock::now();
            queryTiming2 = std::chrono::duration_cast<std::chrono::microseconds>(endQuery - startQuery).count();
            //std::cout << "Run " << i+1 << ": " << queryTiming2 << " us\n";
            timings2.push_back(queryTiming2);
        }

        stddev = 0; 
        avg = 0;
        avg = std::accumulate(timings2.begin(), timings2.end(), 0) / numAvg;
    
        accumulatedAvgTwo += avg;
        stddev = 0;
        for(uint64_t duration : timings2){
            uint64_t diff = duration-avg;
            stddev += diff*diff;
        }
        stddev = sqrt(stddev)/numAvg;

        std::cout << "Duration: " << avg << " +/- " << stddev << "us\n";

        std::cout << std::endl;
    
        std::cout << "Result precalculated: " << result1/numAvg  << "\tResult edgesaved: " << result2/numAvg << "\t Results: " << " " << ((result1 == result2) ? "identical" : "NOT identical") << " " << sNode << " " << tNode <<std::endl;
    }

    std::cout << accumulatedAvgOne << "\t" << accumulatedAvgTwo << std::endl;

    std::cout << "To calculate the distance takes " << accumulatedAvgOne/accumulatedAvgTwo << " times the time" << std::endl;

}

int main() {
    AdjacencyArray adjArray("data/planet_2.graph_2");
    {
        DijkstraImpl sd (adjArray);
        PathAlgorithm &pa = sd;
        // across atlantic
        //benchmarkDijkstra(pa, adjArray, -62, 40, -14, 53.5, 3);
        testDijkstra(pa, adjArray, -85, -180, 85, 180, 10, 100);
    }
}
