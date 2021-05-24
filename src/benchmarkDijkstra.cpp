#include <iostream>
#include <chrono>
#include <numeric>
#include <algorithm>

#include "dijkstra.cpp"

void runDijkstra(PathAlgorithm &pathAlg, AdjacencyArray &adjArray,
double longStart, double latStart, double longGoal, double latGoal)
{
    uint64_t sNode = longLatToNodeId(adjArray, longStart, latStart);
    uint64_t tNode = longLatToNodeId(adjArray, longGoal, latGoal);
    std::vector<uint64_t> idPath;
    pathAlg.findPath(sNode, tNode, idPath);
}

void benchmarkDijkstra(PathAlgorithm &pathAlg, AdjacencyArray &adjArray,
    double longStart, double latStart, double longGoal, double latGoal, int numAvg)
{
    std::vector<uint64_t> timings;
    uint64_t queryTiming;
    for(int i = 0; i<numAvg; ++i){
        auto startQuery = std::chrono::high_resolution_clock::now();

        runDijkstra(pathAlg, adjArray, longStart, latStart, longGoal, latGoal);
        
        auto endQuery = std::chrono::high_resolution_clock::now();
        queryTiming = std::chrono::duration_cast<std::chrono::microseconds>(endQuery - startQuery).count();
        std::cout << "Run " << i+1 << ": " << queryTiming << " ms\n";
        timings.push_back(queryTiming);
    }

    uint64_t stddev = 0; 
    uint64_t avg = std::accumulate(timings.begin(), timings.end(), 0) / numAvg;
    
    ;

    stddev = 0;
    for(uint64_t duration : timings){
        uint64_t diff = duration-avg;
        stddev += diff*diff;
    }
    stddev = sqrt(stddev)/numAvg;

    std::cout << "Duration: " << avg << " +/- " << stddev << "ms\n";
}

int main() {
    AdjacencyArray adjArray;
    FirstDijkstra fd (adjArray);
    PathAlgorithm &pa = fd;
    loadAdjacencyArray(adjArray, "data/worldGrid_1415_707.save");
    // across atlantic
    benchmarkDijkstra(pa, adjArray, -62, 40, -14, 53.5, 3);
}
