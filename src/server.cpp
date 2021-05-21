#include <httplib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

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

/**
 * @brief start calculation of path and create result json file
 * 
 * @param longStart route starting point longitude
 * @param latStart  route starting point latitude
 * @param longGoal  route goal longitude
 * @param latGoal   route goal latitude
 * @param response  out: json response containing calculated route
 */
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
    
// global main html file contents, served to every client
static std::string page;

int main(void)
{

    using namespace httplib;

    Server svr;

    // initialize main page from disk
    std::string path = "static/index.html";
    loadStatic(path, page);

    // serve main page
    svr.Get("/", [](const Request& req, Response& res) {
        res.set_content(page, "text/html");
    });

    // handle client requests of form {longStart: double, latStart:double, longGoal:double, latGoal:double}
    // calculate path, return result path {failure: true/false, path:[[long,lat],[long,lat],...]}
    svr.Post("/getRoute", [](const Request& req, Response& res) {
        double longStart =  std::stod((*req.params.find("longStart")).second);
        double latStart =   std::stod((*req.params.find("latStart")).second);
        double longGoal =   std::stod((*req.params.find("longGoal")).second);
        double latGoal =    std::stod((*req.params.find("latGoal")).second);
        std::cout << "Route from (" << longStart << ", " << latStart << ") to (" << longGoal << ", " << latGoal << ")\n";
        
        std::string response;
        generateReponse(longStart, latStart, longGoal, latGoal, response);
        std::cout << response << std::endl;
        
        // return json with of results path
        res.set_content(response, "application/json");
    });

    svr.listen("localhost", 8080);
}