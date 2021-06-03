#include <httplib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <mutex>
#include <fstream>

# include "dijkstra.cpp"

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
    if(argc > 1){
        inputFileName = std::string(argv[1]);
    }else{
        inputFileName = "data/planet.graph";
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
    static SecondDijkstra secondDijkstra(adjArray);
    static PathAlgorithm &pathAlgorithm = secondDijkstra;
    
    static std::mutex mutex;
    static std::unique_lock<std::mutex> lock (mutex, std::defer_lock);

    // initialize main page from disk
    std::string path = "static/index.html";
    loadStatic(path, page);

    //loadAdjacencyArray(adjArray, "data/worldGrid_1415_707.save");

    svr.Get("/", [](const Request& req, Response& res) {
        res.set_content(page, "text/html");
    });

    svr.Post("/testNode", [](const Request& req, Response& res) {
        double longStart =  std::stod((*req.params.find("longStart")).second);
        double latStart =   std::stod((*req.params.find("latStart")).second);

        uint64_t sNode = longLatToNodeId(adjArray, longStart, latStart);

        //lock.lock();

        std::string response;
        if(isNodeOnLand(adjArray, sNode)){
            response = "{\"water\": false}";
        }else{
            response = "{\"water\": true}";
        }

        //lock.unlock();

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
        

        //std::string response;
        //generateReponse(longStart, latStart, longGoal, latGoal, response);
        //std::cout << response << std::endl;

        uint64_t sNode = longLatToNodeId(adjArray, longStart, latStart);
        uint64_t tNode = longLatToNodeId(adjArray, longGoal, latGoal);

        lock.lock();

        std::vector<uint64_t> idPath;
        pathAlgorithm.reset();
        uint64_t distance = pathAlgorithm.calculateDist(sNode, tNode);
        pathAlgorithm.getPath(idPath);

        std::vector<double> posPath;
        generatePositionPath(posPath, idPath, adjArray);

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