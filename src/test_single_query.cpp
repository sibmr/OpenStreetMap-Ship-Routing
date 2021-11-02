#include "shortestPathUtils.cpp"
#include "CH_query.cpp"

void addForwardSettled(AdjacencyArray &adjArray, CH_query::CH_query &query, std::vector<uint64_t> &settled){
    for(uint64_t nodeId = 0; nodeId<adjArray.width*adjArray.height; ++nodeId){
        if(query.forwardDistance.at(nodeId) < UINT64_MAX){
            settled.push_back(nodeId);
        }
    }
}

void addBackwardSettled(AdjacencyArray &adjArray, CH_query::CH_query &query, std::vector<uint64_t> &settled){
    for(uint64_t nodeId = 0; nodeId<adjArray.width*adjArray.height; ++nodeId){
        if(query.backwardDistance.at(nodeId) < UINT64_MAX){
            settled.push_back(nodeId);
        }
    }
}

void cutLists(std::vector<uint64_t> &initial, std::vector<uint64_t> &cutaway, std::vector<uint64_t> &result){
    for(uint64_t nodeId : initial){
        if(std::find(cutaway.begin(), cutaway.end(), nodeId) == cutaway.end()){
            result.push_back(nodeId);
        }
    }
}

void printVector(std::vector<uint64_t> &vect, std::string printBefore){
    std::cout << "results size " << vect.size() << "\n";
    for(uint64_t nodeId : vect){
        std::cout << printBefore << nodeId << "\n";
    }
}

int main(){
    AdjacencyArray adjArray("data/CHAdjArray_54.graph_2");
    CH_query::CH_query naive(adjArray);
    CH_query::CH_query bidirect(adjArray);

    uint64_t nodeIdStart = longLatToNodeId(adjArray, 53.0789,140.421);
    uint64_t nodeIdGoal = longLatToNodeId(adjArray, 34.2939,-65.3568);

    uint64_t dist1 = naive.calculateDistNaive(nodeIdStart, nodeIdGoal);
    uint64_t dist2 = bidirect.calculateDist(nodeIdStart, nodeIdGoal);

    std::cout << dist1 << "  " << dist2 << "\n";

    std::vector<uint64_t> settledNaiveForward;
    std::vector<uint64_t> settledBidirForward;
    std::vector<uint64_t> settledNaiveBackward;
    std::vector<uint64_t> settledBidirBackward;

    addForwardSettled(adjArray, naive, settledNaiveForward);
    addBackwardSettled(adjArray, naive, settledNaiveBackward);
    addForwardSettled(adjArray, bidirect, settledBidirForward);
    addBackwardSettled(adjArray, bidirect, settledBidirBackward);

    std::cout << "naive forward settled:  " << settledNaiveForward.size() << " naive backward settled: " << settledNaiveBackward.size() << "\n";
    std::cout << "bidir forward settled: " << settledBidirForward.size() << " bidir backward settled: " << settledBidirBackward.size() << "\n";
    std::cout << "naive meeting nodes: " << naive.forwardMinMeetingNodeId << " " << naive.backwardMinMeetingNodeId << "\n";
    std::cout << "naive meeting node prev: " << naive.forwardPrev.at(naive.forwardMinMeetingNodeId) << " " << naive.backwardPrev.at(naive.backwardMinMeetingNodeId) << "\n";
    std::cout << "bidir meeting nodes: " << bidirect.forwardMinMeetingNodeId << " " << bidirect.backwardMinMeetingNodeId << "\n";

    std::vector<uint64_t> resultsForward;
    cutLists(settledNaiveForward, settledBidirForward, resultsForward);
    std::vector<uint64_t> resultsBackward;
    cutLists(settledNaiveBackward, settledBidirBackward, resultsBackward);
    
    printVector(resultsForward, "fw: ");
    printVector(resultsBackward, "bw: ");
    std::cout << "\n";
}