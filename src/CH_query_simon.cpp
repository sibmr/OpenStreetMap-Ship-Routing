#include "PathAlgorithm"
#include "shortestPathUtils.cpp"

class CH_query: public PathAlgorithm{
    public:
        CH_query(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        uint64_t calculateDist(uint64_t endPoint);
        void reset();
        uint64_t getNumNodesPopped();
        void disableNode(uint64_t disabledNode);
        void resetDisabledNodes();
    private:
        uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
        uint64_t mainCalculationLoop();
        std::vector<uint64_t> distance;
        std::vector<HeapElement> heap;
        std::vector<uint64_t> visited;
        std::vector<uint64_t> prev;
        std::vector<bool> disabledNodes;
        AdjacencyArray &adjArray;
        uint64_t constLngDist;
        std::vector<uint64_t> constLatDist;
        uint64_t startPoint, endPoint, lastCalculatedDistance;
        uint64_t numNodesPopped;
};