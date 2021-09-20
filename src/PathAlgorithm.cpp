#ifndef PATH_ALGORITHM
#define PATH_ALGORITHM

#include <vector>
#include <iostream>

/**
 * @brief virtual interface for shortest path algorithms
 * 
 */
class PathAlgorithm {
    public:
        virtual void getPath(std::vector<uint64_t> &path) = 0;
        virtual uint64_t getDist() = 0;
        virtual uint64_t calculateDist(uint64_t startPoint_, uint64_t endPoint_) = 0;
        virtual void reset() = 0;
};

#endif