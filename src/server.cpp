#include <httplib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

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

void generateReponse(double longStart, double latStart, double longGoal, double latGoal, std::string& response){
    
    // TODO: fill with real routing response
    std::vector<double> path {longStart, latStart, longGoal, latGoal};
    
    // response has to be in [[lat,long],...] format, so (long, lat) is swapped
    response += "{\"path\":[";
    response += "[" + std::to_string(path.at(1)) + "," + std::to_string(path.at(0)) + "]";
    for(int i=1; i<path.size()/2; ++i){
        response += ",[" + std::to_string(path.at(2*i+1)) + "," + std::to_string(path.at(2*i)) + "]";
    
    }
    response += "]}";
}
    

static std::string page;

int main(void)
{

    using namespace httplib;

    Server svr;

    std::string path = "static/index.html";
    loadStatic(path, page);

    svr.Get("/", [](const Request& req, Response& res) {
        res.set_content(page, "text/html");
    });

    svr.Post("/getRoute", [](const Request& req, Response& res) {
        double longStart =  std::stod((*req.params.find("longStart")).second);
        double latStart =   std::stod((*req.params.find("latStart")).second);
        double longGoal =   std::stod((*req.params.find("longGoal")).second);
        double latGoal =    std::stod((*req.params.find("latGoal")).second);
        std::cout << "Route from (" << longStart << ", " << latStart << ") to (" << longGoal << ", " << latGoal << ")\n";
        
        std::string response;
        generateReponse(longStart, latStart, longGoal, latGoal, response);
        std::cout << response << std::endl;
        // return geojson with of results path
        res.set_content(response, "application/json");
    });

    // svr.Get(R"(/numbers/(\d+))", [&](const Request& req, Response& res) {
    //     auto numbers = req.matches[1];
    //     res.set_content(numbers, "text/plain");
    // });

    // svr.Get("/body-header-param", [](const Request& req, Response& res) {
    //     if (req.has_header("Content-Length")) {
    //     auto val = req.get_header_value("Content-Length");
    //     }
    //     if (req.has_param("key")) {
    //     auto val = req.get_param_value("key");
    //     }
    //     res.set_content(req.body, "text/plain");
    // });

    // svr.Get("/stop", [&](const Request& req, Response& res) {
    //     svr.stop();
    // });

    svr.listen("localhost", 8080);
}