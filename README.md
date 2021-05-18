# Ship Routing using OpenStreetMap

## Links
OSM:  
http://wiki.openstreetmap.org/wiki/Nodes​  
http://wiki.openstreetmap.org/wiki/Tags  
https://wiki.openstreetmap.org/wiki/Ways  
https://wiki.openstreetmap.org/wiki/Coastline  
Reading Data:  
https://wiki.openstreetmap.org/wiki/PBF_Format#The_code  
https://osmcode.org/libosmium/

## Dependencies

### libosmiumpbfreader:
https://github.com/CanalTP/libosmpbfreader

Install using apt:  

    sudo apt-get install libosmpbf-dev

### libosmiumpbfreader:
https://github.com/yhirose/cpp-httplib

## Installation

### Install Program using CMake

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..  
    make -j4
    make install
    cd .. 

## Usage

### Extract Coastlines
From project root:

    ./build/bin/coastlineExtraction data/planet-coastlines.pbf data/planet-coastlines

### Create File containing PointInPolygon Queries for Grid-Graph creation
From project root:

    ./build/bin/pointInPolygon {path to extracted coastlines} {result save location} {number of grid points along longitude} {number of grid points along latitude}

Example:

    ./build/bin/pointInPolygon ./build/bin/planet-coastlines.save data/worldGrid_1415_707.save 1415 707

### Create Adjacency Array from PointInPolygon Queries file
From project root:

    ./build/bin/buildGraph data/worldGrid_1415_707.save data/wordGraph_1415_707.save

### Start Server
From project root:

    ./build/bin/server

Visit http://localhost:8080 using the Browser

