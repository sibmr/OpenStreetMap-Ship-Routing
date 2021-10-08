# Ship Routing using OpenStreetMap

## Dependencies


### cmake:
Install using apt:  

    sudo apt-get install cmake

### g++:
Install using apt:  

    sudo apt-get install g++


### libosmbfreader:
https://github.com/CanalTP/libosmpbfreader

Install using apt:  

    sudo apt-get install libosmpbf-dev

### cpp-httplib:
https://github.com/yhirose/cpp-httplib

is installed as submodule 

## Installation

### Clone repository

    git clone --recurse-submodules https://github.com/sibmr/OpenStreetMap-Ship-Routing.git

### Install Program using CMake

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..  
    make -j4
    make install
    cd .. 

## Usage

Each program will generate an intermediate file with one of the following suffixes:
* .coastline
* .grid
* .graph

It is recommend to create a data folder to store the intermediate files as well as the planet-coastlines.pbf file.
    
    mkdir data
    mv planet-coastlines.pbf data/

Each program can get either one or two files as input parameter.
The first of the parameter will be the input file and the second one the output file.
If only one file is given, the output file will be saved in the same directory with the corresponding suffix.

The programs have to be executed in the following order (**with** the latter described arguments):

    ./build/bin/coastlineExtraction ...
    ./build/bin/pointInPolygon ...
    ./build/bin/buildGraph ...
    ./build/bin/server ...

### Extract Coastlines
The coastlineExtraction program reads a "*.pbf" file and saves it as "*.coastline" file.
From project root:
    ./build/bin/coastlineExtraction {path to pbf file} {path to new coastline file}

Example:

    ./build/bin/coastlineExtraction data/planet-coastlines.pbf data/planet.coastline

### Create File containing PointInPolygon Queries for Grid-Graph creation
The pointInPolygon program reads a "*.coastline" file and saves it as "*.grid" file.
From project root:

    ./build/bin/pointInPolygon {path to coastline file} {path to new grid file}

Example:

    ./build/bin/pointInPolygon data/planet.coastline data/planet.grid

In addition, you can select the number of nodes to be 100.000 Nodes, one million Nodes or ten million Nodes.
With the "-n" argument followed by "100K", "1M", "10M". 

Example:

    ./build/bin/pointInPolygon -n 100K ./build/bin/planet.coastline data/planet.grid 

### Create Adjacency Array from PointInPolygon Queries file
The buildGraph program reads a "*.grid" file and saves it as "*.graph" file.
From project root:

    ./build/bin/buildGraph {path to grid file} {path to new graph file}

Example: 

    ./build/bin/buildGraph data/planet.grid data/planet.graph

### Use Contraction Hierarchies from Grid file
The buildGraph program reads a "*.grid" file and saves it as "*.graph" file.
With a distinguishable filename based on the level of contraction
From project root:

    ./build/bin/buildGraph -ch {percentage} {path to grid file} 

Example: 

    ./build/bin/buildGraph -ch 0.9 data/planet.grid 
    
### Use Benchmarks
From project root:

    ./build/bin/benchmark {number of tries} {file to adjacency array graph} {file to contraction hioerachies graph}

Example: 

    ./build/bin/buildGraph 100 data/planet.graph data/planet_ch_0.900000.graph 



### Start Server
The server program reads a "*.graph" file and runs a server, where the user can interactively select two points on the map.
It uses currently the contraction hierarchies graph (can be chosen over the commented lines)
From project root:

    ./build/bin/server {path to graph file}

Example:

    ./build/bin/server data/planet_ch_0.900000.graph

Visit http://localhost:8080 using the Browser

## Links
OSM:  
http://wiki.openstreetmap.org/wiki/Nodes​  
http://wiki.openstreetmap.org/wiki/Tags  
https://wiki.openstreetmap.org/wiki/Ways  
https://wiki.openstreetmap.org/wiki/Coastline  
Reading Data:  
https://wiki.openstreetmap.org/wiki/PBF_Format#The_code 