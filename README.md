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

### libosmiumpbfreader:
https://github.com/CanalTP/libosmpbfreader

Dependency:  

    sudo apt-get install libosmpbf-dev


### Approach

* Create Reader object of pbf
* Create Handler that reads coastlines
* Save Coastlines to geoJson (or use data directly)



### install execute program

cmake -DCMAKE_BUILD_TYPE=Release .
make 
make install 

./bin/osmrouting data/planet-coastline.pbf data/planet-coastlin


