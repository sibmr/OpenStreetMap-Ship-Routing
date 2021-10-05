#include <set>

#include "PathAlgorithm.cpp"
#include "shortestPathUtils.cpp"


namespace CH_query{

/**
 * @brief HeapElement for Dijkstra implementation
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

class CH_query: public PathAlgorithm{
    public:
        CH_query(AdjacencyArray &array);
        void getPath(std::vector<uint64_t> &path);
        uint64_t getDist();
        uint64_t calculateDist(uint64_t startPoint, uint64_t endPoint);
        uint64_t calculateDist1(uint64_t startPoint_, uint64_t endPoint_);
        void reset();
        uint64_t getNumNodesPopped();
    private:
        uint64_t fillVectors(uint64_t startPoint, uint64_t endPoint);
        void recursiveUnpackEdge(uint64_t edgeId, std::vector<uint64_t> &unpackedEdges);
        uint64_t getEdgeIdBetween(uint64_t nodeId1, uint64_t nodeId2);
        std::vector<uint64_t> forwardDistance;
        std::vector<uint64_t> backwardDistance;
        std::vector<HeapElement> forwardHeap;
        std::vector<HeapElement> backwardHeap;
        std::set<uint64_t> visited;
        std::vector<uint64_t> forwardPrev;
        std::vector<uint64_t> backwardPrev;
        AdjacencyArray &adjArray;
        uint64_t startPoint, endPoint, lastCalculatedDistance, forwardMinMeetingNodeId, backwardMinMeetingNodeId;
        uint64_t numNodesPopped;
};

/**
 * @brief Construct a new CH_query::CH_query object
 * 
 * Initialize datastructures
 * Initialize distances between nodes
 * 
 * @param array 
 */
CH_query::CH_query(AdjacencyArray &array) : 
                adjArray(array), 
                forwardPrev(array.width*array.height, UINT64_MAX),
                backwardPrev(array.width*array.height, UINT64_MAX)
{
    reset();
}

/**
 * @brief reset datastructures to prepare for next call to calculateDist
 */
void CH_query::reset(){
    lastCalculatedDistance = UINT64_MAX;
    std::vector<uint64_t> initDistF(adjArray.width*adjArray.height, UINT64_MAX);
    forwardDistance = std::move(initDistF);
    std::vector<uint64_t> initDistB(adjArray.width*adjArray.height, UINT64_MAX);
    backwardDistance = std::move(initDistB);
    // no need to reset prev
    forwardHeap.clear();
    backwardHeap.clear();
    // initially no nodes popped
    numNodesPopped = 0;
}


/**
 * @brief efficient dijkstra shortest-path implementation
 * 
 * Uses binary Min(Max)-heap for greedy node visitation strategy
 * 
 * For this to work, reset need to be called first
 * 
 * @param startPoint 
 * @param endPoint 
 * @return uint64_t 
 */
uint64_t CH_query::calculateDist(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;

    forwardHeap.push_back(HeapElement{startPoint, UINT64_MAX, 0});
    backwardHeap.push_back(HeapElement{endPoint, UINT64_MAX, 0});

    std::make_heap(forwardHeap.begin(), forwardHeap.end());
    std::make_heap(backwardHeap.begin(), backwardHeap.end());

    HeapElement forwardFront;
    HeapElement backwardFront;

    HeapElement newForwardFront;
    HeapElement newBackwardFront;

    uint64_t minDistance = UINT64_MAX;

    bool forwardStuck = false;
    bool backwardStuck = false;

    

    while(true){

        bool forwardStalled = false;
        bool backwardStalled = false;

        // get (non-duplicate) heap front
        // avoid duplicate nodes (nodes that were already visited, indicated by higher distance)
        do{
            // no path found
            if(forwardHeap.empty()){
                forwardStuck = true;
            }else{
                std::pop_heap(forwardHeap.begin(), forwardHeap.end());
                newForwardFront = forwardHeap.back();
                forwardHeap.pop_back();
                numNodesPopped++;
            }
        }while ((newForwardFront.dist >= forwardDistance.at(newForwardFront.nodeIdx) || newForwardFront.dist >= minDistance) && !forwardStuck);
        do{
            // no path found
            if(backwardHeap.empty()){
                backwardStuck = true;
            }else{
                std::pop_heap(backwardHeap.begin(), backwardHeap.end());
                newBackwardFront = backwardHeap.back();
                backwardHeap.pop_back();
                numNodesPopped++;
            }
        }while ((newBackwardFront.dist >= backwardDistance.at(newBackwardFront.nodeIdx) || newBackwardFront.dist >= minDistance) && !backwardStuck);
        
        if(forwardStuck && backwardStuck){
            lastCalculatedDistance = minDistance;
            return minDistance;
        }
        
        //std::cout << "heap size: " << forwardHeap.size() << " " << backwardHeap.size() << "\n";

        // settle heap fronts
        // if(!forwardStuck){
        //     forwardFront = newForwardFront;
            
        //     uint64_t stallingDistance = UINT64_MAX;
        //     uint64_t stallingPredecessor = UINT64_MAX;

        //     // stall on demand
        //     for(uint64_t currEdgeId = adjArray.offsets.at(forwardFront.nodeIdx); currEdgeId < adjArray.offsets.at(forwardFront.nodeIdx+1); ++currEdgeId){
        //         // get id of adjacent node for current edge (neighboring node)
        //         uint64_t neighborIdx = adjArray.edges.at(currEdgeId);
        //         if(adjArray.rank.at(neighborIdx) < adjArray.rank.at(forwardFront.nodeIdx)){
                    
        //             if(forwardDistance.at(neighborIdx)==UINT64_MAX){ continue; }
                    
        //             stallingDistance = forwardDistance.at(neighborIdx) + adjArray.distances.at(currEdgeId);
        //             stallingPredecessor = neighborIdx;

        //             if(stallingDistance < forwardFront.dist){
        //                 forwardStalled = true;
        //                 break;    
        //             }
        //         }
        //     }

        //     if(!forwardStalled){
        //         // update distance and previous node of current forward node
        //         forwardDistance.at(forwardFront.nodeIdx) = forwardFront.dist;
        //         forwardPrev.at(forwardFront.nodeIdx) = forwardFront.prev;
        //     }
        //     else{
        //         forwardDistance.at(forwardFront.nodeIdx) = stallingDistance;
        //         forwardDistance.at(forwardFront.nodeIdx) = stallingPredecessor;
        //     }
        //     //forwardStalled = false;
        // }
        // if(!backwardStuck){
        //     backwardFront = newBackwardFront;

        //     uint64_t stallingDistance = UINT64_MAX;
        //     uint64_t stallingPredecessor = UINT64_MAX;

        //     // stall on demand
        //     for(uint64_t currEdgeId = adjArray.offsets.at(backwardFront.nodeIdx); currEdgeId < adjArray.offsets.at(backwardFront.nodeIdx+1); ++currEdgeId){
        //         // get id of adjacent node for current edge (neighboring node)
        //         uint64_t neighborIdx = adjArray.edges.at(currEdgeId);
        //         if(adjArray.rank.at(neighborIdx) < adjArray.rank.at(backwardFront.nodeIdx)){
                    
        //             if(backwardDistance.at(neighborIdx)==UINT64_MAX){ continue; }
                    
        //             stallingDistance = backwardDistance.at(neighborIdx) + adjArray.distances.at(currEdgeId);
        //             stallingPredecessor = neighborIdx;
                    
        //             if(stallingDistance < backwardFront.dist){
        //                 backwardStalled = true;
        //                 break;    
        //             }
        //         }
        //     }
            
        //     if(!backwardStalled){
        //         // update distance and previous node of current backward node
        //         backwardDistance.at(backwardFront.nodeIdx) = backwardFront.dist;
        //         backwardPrev.at(backwardFront.nodeIdx) = backwardFront.prev;
        //     }
        //     else{
        //         backwardDistance.at(backwardFront.nodeIdx) = stallingDistance;
        //         backwardPrev.at(backwardFront.nodeIdx) = stallingPredecessor;
        //     }
        //     // backwardStalled = false;
        // }

        if(!forwardStuck){
            forwardFront = newForwardFront;
            
            // update distance and previous node of current forward node
            forwardDistance.at(forwardFront.nodeIdx) = forwardFront.dist;
            forwardPrev.at(forwardFront.nodeIdx) = forwardFront.prev;
        }
        if(!backwardStuck){
            backwardFront = newBackwardFront;
            
            // update distance and previous node of current forward node
            backwardDistance.at(backwardFront.nodeIdx) = backwardFront.dist;
            backwardPrev.at(backwardFront.nodeIdx) = backwardFront.prev;
        }
        forwardStalled = false;
        backwardStalled = false;

        // termination criterion for bidirectional dijkstra
        if(forwardFront.dist + backwardFront.dist >= minDistance){
            lastCalculatedDistance = minDistance;
            return lastCalculatedDistance;
        }

        if(forwardFront.nodeIdx == backwardFront.nodeIdx){
            uint64_t joinDistance = forwardFront.dist + backwardFront.dist;
            if(joinDistance < minDistance){
                minDistance = joinDistance;
                forwardMinMeetingNodeId = forwardFront.nodeIdx;
                backwardMinMeetingNodeId = backwardFront.nodeIdx;
            }
        }

        //// forward step ////
        // iterate over edges of current forward node
        if(!forwardStuck && !forwardStalled){
            for(uint64_t currEdgeId = adjArray.offsets.at(forwardFront.nodeIdx); currEdgeId < adjArray.offsets.at(forwardFront.nodeIdx+1); ++currEdgeId){
                
                // get id of adjacent node for current edge (neighboring node)
                uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

                // choose length of edge from precalculated lengths
                uint64_t edgeDist = adjArray.distances.at(currEdgeId);

                // calculate new distances
                uint64_t newNeighborDist = forwardFront.dist + edgeDist;
                uint64_t oldNeighborDist = forwardDistance.at(neighborIdx);
                uint64_t neighborRank = adjArray.rank.at(neighborIdx);
                uint64_t currentRank = adjArray.rank.at(forwardFront.nodeIdx);
                if(neighborRank == 0){ neighborRank = 1000; }
                if(currentRank == 0){ currentRank = 1000; }

                // update node distance if it improves (in this case: push new heap node with better distance)
                if(newNeighborDist<oldNeighborDist && neighborRank >= currentRank){
                    forwardHeap.push_back(HeapElement{neighborIdx, forwardFront.nodeIdx, newNeighborDist});
                    std::push_heap(forwardHeap.begin(), forwardHeap.end());
                
                }

                // if neighbor has been discovered in "backward" dijkstra then check if min distance improves
                uint64_t neighborGoalDistance = backwardDistance.at(neighborIdx);
                if(neighborGoalDistance < UINT64_MAX){
                    uint64_t joinDistance = newNeighborDist + neighborGoalDistance;
                    if(joinDistance < minDistance){
                        minDistance = joinDistance;
                        forwardMinMeetingNodeId = forwardFront.nodeIdx;
                        backwardMinMeetingNodeId = neighborIdx;
                    }
                }

                
            }
        }

        //// backward step ////
        // iterate over edges of current backward node
        if(!backwardStuck && !backwardStalled){
            for(uint64_t currEdgeId = adjArray.offsets.at(backwardFront.nodeIdx); currEdgeId < adjArray.offsets.at(backwardFront.nodeIdx+1); ++currEdgeId){
                
                // get id of adjacent node for current edge (neighboring node)
                uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

                // choose length of edge from precalculated lengths
                uint64_t edgeDist = adjArray.distances.at(currEdgeId);

                // calculate new distances
                uint64_t newNeighborDist = backwardFront.dist + edgeDist;
                uint64_t oldNeighborDist = backwardDistance.at(neighborIdx);
                uint64_t neighborRank = adjArray.rank.at(neighborIdx);
                uint64_t currentRank = adjArray.rank.at(backwardFront.nodeIdx);
                if(neighborRank == 0){ neighborRank = 1000; }
                if(currentRank == 0){ currentRank = 1000; }

                // update node distance if it improves (in this case: push new heap node with better distance)
                if(newNeighborDist<oldNeighborDist && neighborRank >= currentRank){
                    backwardHeap.push_back(HeapElement{neighborIdx, backwardFront.nodeIdx, newNeighborDist});
                    std::push_heap(backwardHeap.begin(), backwardHeap.end());    
            
                }

                // if neighbor has been discovered in "forward" dijkstra then check if min distance improves
                uint64_t neighborStartDistance = forwardDistance.at(neighborIdx);
                if(neighborStartDistance < UINT64_MAX){
                    uint64_t joinDistance = newNeighborDist + neighborStartDistance;
                    if(joinDistance < minDistance){
                        minDistance = joinDistance;
                        forwardMinMeetingNodeId = neighborIdx;
                        backwardMinMeetingNodeId = backwardFront.nodeIdx;
                    }
                }
                
            }
        }

    }

}

uint64_t CH_query::calculateDist1(uint64_t startPoint_, uint64_t endPoint_){

    startPoint = startPoint_;
    endPoint = endPoint_;

    forwardHeap.push_back(HeapElement{startPoint, UINT64_MAX, 0});
    backwardHeap.push_back(HeapElement{endPoint, UINT64_MAX, 0});

    std::make_heap(forwardHeap.begin(), forwardHeap.end());
    std::make_heap(backwardHeap.begin(), backwardHeap.end());

    HeapElement forwardFront;
    HeapElement backwardFront;

    HeapElement newForwardFront;
    HeapElement newBackwardFront;

    uint64_t forwardRank = 0;
    uint64_t backwardRank = 0;

    uint64_t minDistance = UINT64_MAX;

    while(true){
        //std::cout << "forward heap " << forwardHeap.size() << "\n";
        if(forwardHeap.empty()){
            break;
        }
        std::pop_heap(forwardHeap.begin(), forwardHeap.end());
        newForwardFront = forwardHeap.back();
        forwardHeap.pop_back();
        forwardFront = newForwardFront;

        if(forwardFront.dist >= forwardDistance.at(forwardFront.nodeIdx)){
            continue;
        }else{
            numNodesPopped++;
        }

        // update distance and previous node of current forward node
        forwardDistance.at(forwardFront.nodeIdx) = forwardFront.dist;
        forwardPrev.at(forwardFront.nodeIdx) = forwardFront.prev;

        //// forward step ////
        // iterate over edges of current forward node
        for(uint64_t currEdgeId = adjArray.offsets.at(forwardFront.nodeIdx); currEdgeId < adjArray.offsets.at(forwardFront.nodeIdx+1); ++currEdgeId){
            //std::cout << currEdgeId << " " << adjArray.offsets.at(forwardFront.nodeIdx) << " " << adjArray.offsets.at(forwardFront.nodeIdx+1) << "\n";
            // get id of adjacent node for current edge (neighboring node)
            uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

            // choose length of edge from precalculated lengths
            uint64_t edgeDist = adjArray.distances.at(currEdgeId);

            // calculate new distances
            uint64_t newNeighborDist = forwardFront.dist + edgeDist;
            uint64_t oldNeighborDist = forwardDistance.at(neighborIdx);
            uint64_t neighborRank = adjArray.rank.at(neighborIdx);
            uint64_t currentRank = adjArray.rank.at(forwardFront.nodeIdx);
            if(neighborRank == 0){ neighborRank = 1000; }
            if(currentRank == 0){ currentRank = 1000; }

            // update node distance if it improves (in this case: push new heap node with better distance)
            if(newNeighborDist<oldNeighborDist && neighborRank >= currentRank){
                forwardHeap.push_back(HeapElement{neighborIdx, forwardFront.nodeIdx, newNeighborDist});
                std::push_heap(forwardHeap.begin(), forwardHeap.end());
            }   
        }
    }

    std::vector<uint64_t> backwardSettledNodes;
    while(true){
        //std::cout << "backward heap " << backwardHeap.size() << "\n";
        if(backwardHeap.empty()){
            break;
        }
        std::pop_heap(backwardHeap.begin(), backwardHeap.end());
        newBackwardFront = backwardHeap.back();
        backwardHeap.pop_back();
        backwardFront = newBackwardFront;

        if(backwardFront.dist >= backwardDistance.at(backwardFront.nodeIdx)){
            continue;
        }else{
            numNodesPopped++;
        }


        // update distance and previous node of current forward node
        backwardDistance.at(backwardFront.nodeIdx) = backwardFront.dist;
        backwardPrev.at(backwardFront.nodeIdx) = backwardFront.prev;
        backwardSettledNodes.push_back(backwardFront.nodeIdx);

        //// forward step ////
        // iterate over edges of current forward node
        for(uint64_t currEdgeId = adjArray.offsets.at(backwardFront.nodeIdx); currEdgeId < adjArray.offsets.at(backwardFront.nodeIdx+1); ++currEdgeId){

            // get id of adjacent node for current edge (neighboring node)
            uint64_t neighborIdx = adjArray.edges.at(currEdgeId);

            // choose length of edge from precalculated lengths
            uint64_t edgeDist = adjArray.distances.at(currEdgeId);

            // calculate new distances
            uint64_t newNeighborDist = backwardFront.dist + edgeDist;
            uint64_t oldNeighborDist = backwardDistance.at(neighborIdx);
            uint64_t neighborRank = adjArray.rank.at(neighborIdx);
            uint64_t currentRank = adjArray.rank.at(backwardFront.nodeIdx);
            if(neighborRank == 0){ neighborRank = 1000; }
            if(currentRank == 0){ currentRank = 1000; }

            // update node distance if it improves (in this case: push new heap node with better distance)
            if(newNeighborDist<oldNeighborDist && neighborRank >= currentRank){
                backwardHeap.push_back(HeapElement{neighborIdx, backwardFront.nodeIdx, newNeighborDist});
                std::push_heap(backwardHeap.begin(), backwardHeap.end());
            }


        }


    }

    minDistance = UINT64_MAX;
    uint64_t meetingNode = UINT64_MAX;
    for(uint64_t nodeId : backwardSettledNodes){
        uint64_t bwdist = backwardDistance.at(nodeId);
        uint64_t fwdist = forwardDistance.at(nodeId);
        if(bwdist < UINT64_MAX && fwdist < UINT64_MAX){
            if(minDistance > fwdist + bwdist){
                minDistance = fwdist + bwdist;
                meetingNode = nodeId;
            }
        }
    }

    forwardMinMeetingNodeId = meetingNode;
    backwardMinMeetingNodeId = meetingNode;

    lastCalculatedDistance = minDistance;
    //std::cout << lastCalculatedDistance << "\n";
    return lastCalculatedDistance;
}

/**
 * @brief returns distance of last call to calculateDist
 * 
 * @return uint64_t 
 */
uint64_t CH_query::getDist(){
    return lastCalculatedDistance;
}

void CH_query::recursiveUnpackEdge(uint64_t edgeId, std::vector<uint64_t> &unpackedNodes){
    Edge &edge = adjArray.allEdgeInfo.at(edgeId);
    //std::cout << "num shortcut edges: " << edge.shortcutPathEdges.size() << ", v1: " <<  edge.v1 << ", v2: " <<  edge.v2 << "\n";
    if(edge.shortcutPathEdges.size() > 0){
        for(uint64_t edgeId : edge.shortcutPathEdges){
            recursiveUnpackEdge(edgeId, unpackedNodes);
        }
    }else{
        unpackedNodes.push_back(edge.v1);
        unpackedNodes.push_back(edge.v2);
        //std::cout << edge.v1 << ", rank " << adjArray.rank.at(edge.v1) << " -> " << edge.v2 << ", rank " << adjArray.rank.at(edge.v2) << "\n"; 
    }
}

uint64_t CH_query::getEdgeIdBetween(uint64_t nodeId1, uint64_t nodeId2){
    for(uint64_t currEdgeIndex = adjArray.offsets.at(nodeId1); currEdgeIndex < adjArray.offsets.at(nodeId1+1); ++currEdgeIndex){
        uint64_t adjacentNodeId = adjArray.edges.at(currEdgeIndex);
        if(adjacentNodeId == nodeId2){
            return adjArray.edgeIds.at(currEdgeIndex);
        }
    }
    return UINT64_MAX;
}

/**
 * @brief retrieve path calculated by the last call to calculateDist
 * 
 * @param path 
 */
void CH_query::getPath(std::vector<uint64_t> &path){
    std::vector<uint64_t> forwardPath;
    std::vector<uint64_t> backwardPath;
    if(lastCalculatedDistance < UINT64_MAX){
        // build up path
        uint64_t prevNode = forwardMinMeetingNodeId;
        uint64_t currNode = UINT64_MAX;
        while(currNode != startPoint){
            currNode = forwardPrev.at(prevNode);
            
            std::cout << "p1: " << startPoint << " " << endPoint << "\n";
            std::cout << currNode << "\n";
            
            uint64_t currEdgeId = getEdgeIdBetween(prevNode, currNode);
            std::vector<uint64_t> unpackedNodeIds;
            recursiveUnpackEdge(currEdgeId, unpackedNodeIds);
            for(uint64_t i = 0; i<unpackedNodeIds.size(); i+=2){
                forwardPath.push_back(unpackedNodeIds.at(i));
            }
            forwardPath.push_back(unpackedNodeIds.at(unpackedNodeIds.size()-1));
            
            prevNode = currNode;
        }
        prevNode = backwardMinMeetingNodeId;
        currNode = UINT64_MAX;
        while(currNode != endPoint){
            currNode = backwardPrev.at(prevNode);
            
            std::cout << "p2: " << startPoint << " " << endPoint << "\n";
            std::cout << currNode << "\n";
            
            uint64_t currEdgeId = getEdgeIdBetween(prevNode, currNode);
            std::vector<uint64_t> unpackedNodeIds;
            recursiveUnpackEdge(currEdgeId, unpackedNodeIds);
            for(uint64_t i = 0; i<unpackedNodeIds.size(); i+=2){
                backwardPath.push_back(unpackedNodeIds.at(i));
            }
            backwardPath.push_back(unpackedNodeIds.at(unpackedNodeIds.size()-1));
            
            prevNode = currNode;
        }
        
        for(int i = forwardPath.size()-1; i>=0; --i){
            path.push_back(forwardPath.at(i));
        }

        std::cout << "forward path end: " << path.size();

        if(forwardMinMeetingNodeId != backwardMinMeetingNodeId){
            std::vector<uint64_t> unpackedNodeIds;
            uint64_t edgeId = getEdgeIdBetween(forwardMinMeetingNodeId, backwardMinMeetingNodeId);
            recursiveUnpackEdge(edgeId, unpackedNodeIds);
            for(uint64_t i = 0; i<unpackedNodeIds.size(); i+=2){
                    path.push_back(unpackedNodeIds.at(i));
            }

            path.push_back(unpackedNodeIds.at(unpackedNodeIds.size()-1));

            std::cout << "intermediate path end: " << path.size() << "\n";
        }
       

        for(int i = 0; i<backwardPath.size(); ++i){
            path.push_back(backwardPath.at(i));
        }

        // print path
        std::cout << "Path:" << std::endl;
        for (std::vector<uint64_t>::iterator it = path.begin(); it != path.end(); ++it) {
            std::cout << *it << " ";
        }
        std::cout <<  std::endl;
        std::cout << "dist: " << lastCalculatedDistance/1000 << "km" << std::endl;
    }else{
        std::cout << "no path found" << std::endl;
        //path.push_back(startPoint);
        //path.push_back(endPoint);
    }
    
}

uint64_t CH_query::getNumNodesPopped(){
    return numNodesPopped;
}

}