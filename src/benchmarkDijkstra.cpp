#include <iostream>
#include <chrono>
#include <numeric>
#include <algorithm>

#include "dijkstra.cpp"

void runDijkstra(AdjacencyArray &adjArray, double longStart, double latStart, double longGoal, double latGoal){
    uint64_t sNode = longLatToNodeId(adjArray, longStart, latStart);
    uint64_t tNode = longLatToNodeId(adjArray, longGoal, latGoal);
    std::vector<uint64_t> idPath;
    generatePath(idPath, sNode, tNode, adjArray);
}

void benchmarkDijkstra(AdjacencyArray &adjArray, double longStart, double latStart, double longGoal, double latGoal, int numAvg){
    std::vector<double> timings;
    double queryTiming;
    for(int i = 0; i<numAvg; ++i){
        auto startQuery = std::chrono::high_resolution_clock::now();

        runDijkstra(adjArray, longStart, latStart, longGoal, latGoal);
        
        auto endQuery = std::chrono::high_resolution_clock::now();
        queryTiming = std::chrono::duration_cast<std::chrono::microseconds>(endQuery - startQuery).count();
        std::cout << "Run " << i+1 << ": " << queryTiming << " ms\n";
        timings.push_back(queryTiming);
    }

    double stddev = 0; 
    double avg = std::accumulate(timings.begin(), timings.end(), 0) / numAvg;
    
    std::transform(
        timings.begin(), timings.end(), timings.begin(), 
        [avg](double dval) -> double {return std::pow(dval-avg, 2);}
        );
    stddev = sqrt(
        std::accumulate(timings.begin(), timings.end(), 0)/numAvg
        );

    std::cout << "Duration: " << avg << " +/- " << stddev << "ms\n";
}

int main() {
    AdjacencyArray adjArray;
    loadAdjacencyArray(adjArray, "data/worldGrid_1415_707.save");
    // across atlantic
    benchmarkDijkstra(adjArray, -62, 40, -14, 53.5, 3);
}
