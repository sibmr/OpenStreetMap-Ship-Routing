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

## Meeting Notes
Two neighboring Nodes in the water -> We assume there is water in between the nodes

## First Task

Read PBF File and extract coastlines.  
Output them to the console or display them by converting them to geoJSON.

### Installing Osmium
Clone https://github.com/osmcode/libosmium  
Install Dependencies:  

    apt-get install -q -y \
    cmake \
    doxygen \
    g++ \
    git \
    graphviz \
    libboost-dev \
    libbz2-dev \
    libexpat1-dev \
    libgdal-dev \
    libgeos++-dev \
    libproj-dev \
    libsparsehash-dev \
    make \
    ruby \
    ruby-json \
    spatialite-bin \
    zlib1g-dev \
    cmake-curses-gui


Build using:

    mkdir build
    cd build
    cmake ..
    ccmake ..
This will open a curses gui.  
Configure ccmake, then press c, then press g.
Finally, if there were no errors:  
    
    make
Or

    make -j4

### Using Osmium

* Create Reader object of pbf
* Create Handler that reads coastlines
* Save Coastlines to geoJson (or use data directly)