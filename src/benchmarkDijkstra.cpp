#include <iostream>
#include <chrono>
#include <numeric>
#include <algorithm>

#include "dijkstra.cpp"

void benchmarkDijkstra(PathAlgorithm &pathAlg, AdjacencyArray &adjArray,
    double longStart, double latStart, double longGoal, double latGoal, int numAvg)
{
    std::vector<uint64_t> timings;
    uint64_t queryTiming;

    uint64_t sNode = longLatToNodeId(adjArray, longStart, latStart);
    uint64_t tNode = longLatToNodeId(adjArray, longGoal, latGoal);

    for(int i = 0; i<numAvg; ++i){

        pathAlg.reset();

        auto startQuery = std::chrono::high_resolution_clock::now();

        pathAlg.calculateDist(sNode, tNode);
        
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

int main() {
    AdjacencyArray adjArray("data/planet.graph");
    {
        FirstDijkstra fd (adjArray);
        PathAlgorithm &pa = fd;
        // across atlantic
        benchmarkDijkstra(pa, adjArray, -62, 40, -14, 53.5, 3);
    }
    {
        SecondDijkstra sd (adjArray);
        PathAlgorithm &pa = sd;
        // across atlantic
        benchmarkDijkstra(pa, adjArray, -62, 40, -14, 53.5, 3);
    }
}
