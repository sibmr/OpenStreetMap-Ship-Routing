#include <httplib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <mutex>
#include <fstream>

#include "Dijkstra.cpp"
#include "Bidirectional_Dijkstra.cpp"
#include "A_star.cpp"
#include "CH_query.cpp"
#include "MultiDijkstra.cpp"
#include "CH_Astar.cpp"

/**
 * @brief load static file (html,js,css) from disk
 * 
 * @param path      in: path to file
 * @param loaded    out: string with file contents 
 */
void loadStatic(std::string& path, std::string& loaded){
    std::ifstream staticFile (path);
    std::string line;
    if(staticFile.is_open())
    {
        while ( getline (staticFile,line) )
            {
                loaded += line + '\n';
            }
    }
}


int main(int argc, char** argv)
{
    std::string inputFileName;
    std::string pathAlg;
    if(argc > 2){
        inputFileName = std::string(argv[1]);
        pathAlg = std::string(argv[2]);
    }else{
        inputFileName = "data/planet.graph";
        pathAlg = "dijkstra";
        std::cout << "no input file given assume " <<  inputFileName << std::endl;
    }

    {
        std::ifstream f(inputFileName);
        if(!f.good()){
            std::cout << "file: " << inputFileName << " not found\n";
            return 0;
        }
    }
    

    using namespace httplib;

    Server svr;

    // global main html file contents, served to every client
    static std::string page;

    static AdjacencyArray adjArray(inputFileName);
    
    Dijkstra::Dijkstra dijk (adjArray);
    BidirectionalDijkstra::BidirectionalDijkstra bidijk (adjArray);
    A_star::A_star astar (adjArray);
    A_star::A_star_rectangular astar_rect (adjArray);
    CH_query::CH_query chquery (adjArray);
    CH_Astar::A_star chqueryStar (adjArray);
    CH_Astar::A_star_rectangular chqueryRectStar (adjArray);

    std::map<std::string, PathAlgorithm*> algs{
        {"dijkstra", &dijk},
        {"bidijkstra", &bidijk},
        {"astar", &astar},
        {"rectastar", &astar_rect},
        {"chquery", &chquery},
        {"chastar", &chqueryStar},
        {"chrectastar", &chqueryRectStar}
    };
    
    
    static PathAlgorithm &pathAlgorithm = *algs.at(pathAlg);
    
    static std::mutex mutex;
    static std::unique_lock<std::mutex> lock (mutex, std::defer_lock);

    // initialize main page from disk
    std::string path = "static/index.html";
    loadStatic(path, page);

    svr.Get("/", [](const Request& req, Response& res) {
        res.set_content(page, "text/html");
    });

    // handle client requests of form {long: double, lat:double}
    // check if node is on water or on land {failure: true/false, water:true/false}
    svr.Post("/testNode", [](const Request& req, Response& res) {
        double lng =  std::stod((*req.params.find("long")).second);
        double lat =   std::stod((*req.params.find("lat")).second);

        uint64_t sNode = longLatToNodeId(adjArray, lng, lat);

        std::string response = "{\"water\": " +  std::to_string(!isNodeOnLand(adjArray, sNode)) + "}";

        // return geojson with of results path
        std::cout << response << std::endl;
        
        // return json with of results path
        res.set_content(response, "application/json");
        
    });


    // handle client requests of form {longStart: double, latStart:double, longGoal:double, latGoal:double}
    // calculate path, return result path {failure: true/false, path:[[long,lat],[long,lat],...]}
    svr.Post("/getRoute", [](const Request& req, Response& res) {
        double longStart =  std::stod((*req.params.find("longStart")).second);
        double latStart =   std::stod((*req.params.find("latStart")).second);
        double longGoal =   std::stod((*req.params.find("longGoal")).second);
        double latGoal =    std::stod((*req.params.find("latGoal")).second);
        std::cout << "Route from (" << longStart << ", " << latStart << ") to (" << longGoal << ", " << latGoal << ")\n";
        
        uint64_t sNode = longLatToNodeId(adjArray, longStart, latStart);
        uint64_t tNode = longLatToNodeId(adjArray, longGoal, latGoal);

        lock.lock();

        std::vector<uint64_t> idPath;
        pathAlgorithm.reset();
        uint64_t distance = pathAlgorithm.calculateDist(sNode, tNode);

        std::vector<double> posPath;

        if(distance < UINT64_MAX){
            pathAlgorithm.getPath(idPath);
            generatePositionPath(posPath, idPath, adjArray);
        }

        lock.unlock();

        std::string response;
        generateReponse(posPath, response, distance);

        // return geojson with of results path
        std::cout << response << std::endl;
        
        // return json with of results path
        res.set_content(response, "application/json");
        
    });

    svr.listen("localhost", 8080);
}