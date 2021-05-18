#include <httplib.h>
#include <string>
#include <iostream>
#include <fstream>

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