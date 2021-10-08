#include "shortestPathUtils.cpp"
#include "PathAlgorithm.cpp"

namespace ChQuery{

    /**
     * @brief HeapElement for DijkstraImpl implementation
     */
    struct HeapElement {
        // for normal dijkstra, dist is the current distance to this node
        uint64_t nodeIdx, prev, dist;

        /**
         * @brief "reverse" comparison function turning max-heap into min-heap 
         * 
         * @param a         other HeapElement
         * @return true     if own distance is bigger than others distance,
         * @return false    otherwise
         */
        bool operator<(const HeapElement &a){
            return dist > a.dist;
        }
    };


    /**
     * @brief My bidirectional Dijkstra Algorithm
     * Based on our improved Dijkstra similar implementation to https://github.com/Lesstat/dijkstra-performance-study/
     */
    class ChQuery: public PathAlgorithm{
        public:
            ChQuery(AdjacencyArray &array);
            void getPath(std::vector<uint64_t> &path);
            uint64_t getDist();
            uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
            uint64_t calculateDistSavedEdges(uint64_t startPoint, uint64_t endPoint);
            void reset();
            uint64_t getNumNodesPopped();
        private:
            uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
            std::vector<uint64_t> distanceSource;
            std::vector<uint64_t> distanceTarget;
            std::vector<HeapElement> heapSource;
            std::vector<HeapElement> heapTarget;
            std::vector<uint64_t> prevSource;
            std::vector<uint64_t> prevTarget;
            AdjacencyArray &adjArray;
            uint64_t constLngDist;
            std::vector<uint64_t> constLatDist;
            uint64_t startPoint, endPoint, lastCalculatedDistance;
            uint64_t meetingNode;
            uint64_t numNodesPopped;
            uint64_t maxRank;
    };


    /**
     * @brief Construct a new Dijkstra:: DijkstraImpl object from both sides
     * 
     * Initialize datastructures
     * Initialize distances between nodes
     * 
     * @param array 
     */
    ChQuery::ChQuery(AdjacencyArray &array) : adjArray(array), prevSource(array.width*array.height, UINT64_MAX), prevTarget(array.width*array.height, UINT64_MAX){
        maxRank = 0;

        for(uint64_t currNode = 0; currNode < adjArray.rank.size(); currNode++){
            if(adjArray.rank.at(currNode) > maxRank){
                maxRank = adjArray.rank.at(currNode);
            }
        }
        reset();

        // calculate distances between nodes:
        // constLngDist is the distance between nodes with the same longitude but different latitude
        // constLatDist is the distance between nodes with the same latitude but different logitude
        // distance between (i,0) and (i,1)
        constLngDist = nodeDistance(adjArray, 0, 1);
        // for each constant latitude "ring" around the globe, the distance is different
        // mirroring is disregarded
        for(uint64_t i = 0; i<adjArray.height; ++i){
            // distance between (0,i) and (1,i)
            constLatDist.push_back(nodeDistance(adjArray, i, adjArray.height+i));
        }
    }

    /**
     * @brief reset datastructures to prepare for next call to calculateDist
     */
    void ChQuery::reset(){
        std::vector<uint64_t> initDist(adjArray.width*adjArray.height, UINT64_MAX);
        std::vector<uint64_t> initDist2(adjArray.width*adjArray.height, UINT64_MAX);
        distanceSource = std::move(initDist);
        distanceTarget = std::move(initDist2);
        // no need to reset prev
        heapSource.clear();
        heapTarget.clear();

        // initially no nodes are popped
        numNodesPopped = 0;

    }



    /**
     * @brief efficient dijkstra shortest-path implementation
     * 
     * Uses binary Min(Max)-heap for greedy node visitation strategy
     * 
     * Saved Distances
     * 
     * @param startPoint 
     * @param endPoint 
     * @return uint64_t 
     */
    uint64_t ChQuery::calculateDist(uint64_t startPoint, uint64_t endPoint){

        // save start and end Node for path export
        this->startPoint = startPoint;
        this->endPoint = endPoint;
        uint64_t currDist = UINT64_MAX;

        uint64_t rankSourceTree = 0;
        uint64_t rankTargetTree = 0;

        bool sourceTreeReady = false;
        bool targetTreeReady = false;


        // init rank
        rankSourceTree = adjArray.rank.at(startPoint);
        rankTargetTree = adjArray.rank.at(endPoint);

        uint64_t currRank = 0;

        heapSource.push_back(HeapElement{startPoint, UINT64_MAX, 0});
        heapTarget.push_back(HeapElement{endPoint, UINT64_MAX, 0});

        std::make_heap(heapSource.begin(), heapSource.end());
        std::make_heap(heapTarget.begin(), heapTarget.end());

        HeapElement frontSource;
        HeapElement frontTarget;
        while(true){
            if(heapSource.empty() && heapTarget.empty()){
                lastCalculatedDistance = currDist;
                return currDist;
            }

            // find next node in sourceTree not already visited
            if(!sourceTreeReady){
                do{
                    if(heapSource.empty()){
                        sourceTreeReady = true;
                    }else{
                        std::pop_heap(heapSource.begin(), heapSource.end());
                        frontSource = heapSource.back();
                        heapSource.pop_back();
                        // one node will be popped inSourceTree 
                        numNodesPopped++;

                        if(frontSource.dist > currDist){
                            sourceTreeReady = true;
                        }
                    }
                }while(frontSource.dist >= distanceSource.at(frontSource.nodeIdx) && !heapSource.empty() && !sourceTreeReady);
            }
            



            // find next node in targetTree not already visited
            if(!targetTreeReady){
                do{
                    if(heapTarget.empty()){
                        targetTreeReady = true;
                    }else{
                        std::pop_heap(heapTarget.begin(), heapTarget.end());
                        frontTarget = heapTarget.back();
                        heapTarget.pop_back();
                        // one node will be popped in targetTree
                        numNodesPopped++;

                        if(frontTarget.dist > currDist){
                            targetTreeReady = true; 
                        }
                    }
                }while(frontTarget.dist >= distanceTarget.at(frontTarget.nodeIdx) && !heapTarget.empty());


            }


            // ch termination
            if(sourceTreeReady && targetTreeReady){
                lastCalculatedDistance = currDist;
                return lastCalculatedDistance;
            }

            if(!sourceTreeReady){
                rankSourceTree = adjArray.rank.at(frontSource.nodeIdx);

                distanceSource.at(frontSource.nodeIdx) = frontSource.dist;
                prevSource.at(frontSource.nodeIdx) = frontSource.prev;

                // meetingNode found
                if(distanceTarget.at(frontSource.nodeIdx) < UINT64_MAX){
                    uint64_t tmpDist = distanceTarget.at(frontSource.nodeIdx) + frontSource.dist;
                    if(tmpDist < currDist){
                        currDist = tmpDist;
                        meetingNode = frontSource.nodeIdx;
                    }
                }
            }


            if( !targetTreeReady){
                rankTargetTree = adjArray.rank.at(frontTarget.nodeIdx);
                // meetingNode found
                distanceTarget.at(frontTarget.nodeIdx) = frontTarget.dist;
                prevTarget.at(frontTarget.nodeIdx) = frontTarget.prev;
                if(distanceSource.at(frontTarget.nodeIdx) < UINT64_MAX){
                    uint64_t tmpDist = distanceSource.at(frontTarget.nodeIdx) + frontTarget.dist;
                    if(tmpDist < currDist){
                        currDist = tmpDist;
                        meetingNode = frontTarget.nodeIdx;
                    }
                }
            }





            // one step in sourceTree
            if(!sourceTreeReady){
                distanceSource.at(frontSource.nodeIdx) = frontSource.dist;
                prevSource.at(frontSource.nodeIdx) = frontSource.prev;

                for(uint64_t currEdgeId = adjArray.offsets.at(frontSource.nodeIdx); currEdgeId < adjArray.offsets.at(frontSource.nodeIdx+1); currEdgeId++){
                    uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

                    currRank = adjArray.rank.at(neighborIdx);

                    // only push nodes with higer rank
                    if(rankSourceTree < currRank || currRank == maxRank){
    
                        // choose length of edge from precalculated lengths
                        uint64_t edgeDist = adjArray.distances.at(currEdgeId);
    
                        uint64_t newNeighborDist = frontSource.dist + edgeDist;
                        uint64_t oldNeighborDist = distanceSource.at(neighborIdx);

    
                        if(newNeighborDist<oldNeighborDist && currDist == UINT64_MAX){
                            // do not update distance array: only update for distances that are final
                            // distance.at(neighborIdx) = newNeighborDist;
                            heapSource.push_back(HeapElement{neighborIdx, frontSource.nodeIdx, newNeighborDist});
                            std::push_heap(heapSource.begin(), heapSource.end());
                        }
                    }

                }
            }

            // one step in targetTree
            if(!targetTreeReady){
                distanceTarget.at(frontTarget.nodeIdx) = frontTarget.dist;
                prevTarget.at(frontTarget.nodeIdx) = frontTarget.prev;

                for(uint64_t currEdgeId = adjArray.offsets.at(frontTarget.nodeIdx); currEdgeId < adjArray.offsets.at(frontTarget.nodeIdx+1); currEdgeId++){
                    uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

                    currRank = adjArray.rank.at(neighborIdx);

                    if(rankTargetTree < currRank || currRank == maxRank){


                        // choose length of edge from precalculated lengths
                        uint64_t edgeDist = adjArray.distances.at(currEdgeId);

                        uint64_t newNeighborDist = frontTarget.dist + edgeDist;
                        uint64_t oldNeighborDist = distanceTarget.at(neighborIdx);


                        // from 4 times more nodes to 2 times more nodes
                        if(newNeighborDist<oldNeighborDist && currDist == UINT64_MAX){
                            // do not update distance array: only update for distances that are final
                            // distance.at(neighborIdx) = newNeighborDist;
                            heapTarget.push_back(HeapElement{neighborIdx, frontTarget.nodeIdx, newNeighborDist});
                            std::push_heap(heapTarget.begin(), heapTarget.end());
                        }
                    }

                }
            }

        }

    }

    /**
     * @brief returns distance of last call to calculateDist
     * 
     * @return uint64_t 
     */
    uint64_t ChQuery::getDist(){
        return lastCalculatedDistance;
    }

    uint64_t ChQuery::getNumNodesPopped(){
        return numNodesPopped;
    }

    /**
     * @brief retrieve path calculated by the last call to calculateDist
     * 
     * @param path 
     */
    void ChQuery::getPath(std::vector<uint64_t> &path){
        std::vector<uint64_t> tmpVector;
        if(getDist() < UINT64_MAX){
            // build up backward path
            uint64_t currNode = meetingNode;
            tmpVector.push_back(currNode);
            while(currNode != endPoint){
                currNode = prevTarget.at(currNode);
                tmpVector.push_back(currNode);
            }

            std::reverse(tmpVector.begin(), tmpVector.end());
            for (int i = 0; i < tmpVector.size(); i++){
                path.push_back(tmpVector.at(i));
            }

            path.push_back(meetingNode);
            //path.push_back(forwardMeetingNode);

            // build up forward path
            currNode = meetingNode;
            while(currNode != startPoint){
                currNode = prevSource.at(currNode);
                path.push_back(currNode);
            }


            // print path
            std::cout << "Path:" << std::endl;
            for (std::vector<uint64_t>::iterator it = path.begin(); it != path.end(); ++it) {
                std::cout << *it << " ";
            }
            std::cout <<  std::endl;
            std::cout << "dist: " << getDist()/1000 << "km" << std::endl;
        }else{
            std::cout << "no path found" << std::endl;
            //path.push_back(startPoint);
            //path.push_back(endPoint);
        }
        return;

    }


    void test() {
        AdjacencyArray adjArray("data/worldGrid_1415_707.save");
        ChQuery sd(adjArray);
        PathAlgorithm &pa = sd;

        std::vector<uint64_t> path;
        pa.reset();
        pa.calculateDist(1001, 16001);
        pa.getPath(path);
        std::vector<double> posPath;
        generatePositionPath(posPath, path, adjArray);


        std::string response;
        generateReponse(posPath, response, pa.getDist());
        std::cout << response << std::endl;

        int counter_one = 0;
        int counter_zero = 0;
        for (int i = 0; i < adjArray.nodes.size(); i++){
            if(adjArray.nodes.at(i) == 0){
                counter_zero++;
            }else
            {
                counter_one++;
            }
        }
        std::cout << counter_zero << " " << counter_one << " " << counter_one + counter_zero << std::endl;
        std::cout << adjArray.edges.size() << std::endl;
        testLatLongDistance();
    }
}