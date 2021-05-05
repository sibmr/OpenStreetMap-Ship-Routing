# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <math.h>

struct Node {
    uint64_t id;
    double longitude;
    double latitude;
};

/**
 * to the left there is land, to the right there is ocean
 **/
struct Edge {
    Node &sourceNode;
    Node &destinationNode;  
};


struct Edge2 {
    double sourceLongitude;
    double sourceLatitude;
    double targetLongitude;
    double targetLatitude;
};

struct Vec3D {
    double x;
    double y;
    double z;
};

void printVec(Vec3D vec){
    std::cout << "x=" << vec.x << " y=" << vec.y << " z=" << vec.z << "\n";
}

void printEdge(Edge2 edge){
    std::cout << "source longitude=" << edge.sourceLongitude << " ,latitude=" << edge.sourceLatitude 
    << " ,destination longitude=" << edge.targetLongitude << " ,latitude=" << edge.targetLatitude << "\n";
}

void printEdge(Edge edge){
    std::cout << "source longitude=" << edge.sourceNode.longitude << " ,latitude=" << edge.sourceNode.latitude 
    << " ,destination longitude=" << edge.destinationNode.longitude << " ,latitude=" << edge.destinationNode.latitude << "\n";
}

/**
 * supposed position in 3D on a sphere of radius 1
 **/
Vec3D sphericalToRectangular(double longitude, double latitude){
    double theta = (M_PI*(latitude+90))/180.0;
    double phi = (M_PI*longitude)/180.0;
    double stheta = sin(theta);
    double x = stheta*cos(phi);
    double y = stheta*sin(phi);
    double z = cos(theta);
    return Vec3D{x, y, z};
}

Vec3D rectangularToSpherical(Vec3D coords){
    double r = sqrt(coords.x*coords.x + coords.y*coords.y + coords.z*coords.z);
    double latitude = -atan2(coords.z,sqrt(coords.x*coords.x + coords.y*coords.y))*(180/M_PI);
    double longitude = 0;
    longitude = atan2(coords.y,coords.x)*(180/M_PI);
    return Vec3D{r, longitude, latitude};
}

Vec3D crossProduct(Vec3D v1, Vec3D v2){
    return Vec3D{
        v1.y*v2.z - v1.z*v2.y,
        v1.z*v2.x - v1.x*v2.z,
        v1.x*v2.y - v1.y*v2.x
    };
}

Vec3D add(Vec3D v1, Vec3D v2){
    return Vec3D{
        v1.x+v2.x,
        v1.y+v2.y,
        v1.z+v2.z
    };
}

Vec3D sub(Vec3D v1, Vec3D v2){
    return Vec3D{
        v1.x-v2.x,
        v1.y-v2.y,
        v1.z-v2.z
    };
}

Vec3D div(Vec3D v, double divisor){
    return Vec3D{
        v.x/divisor,
        v.y/divisor,
        v.z/divisor
    };
}

double dotProduct(Vec3D v1, Vec3D v2){
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

Vec3D normalize(Vec3D v){
    double norm = std::sqrt(dotProduct(v,v));
    return div(v, norm); 
}

bool isNodeLeftOfEdge(Node n, Edge e){
    Vec3D c = sphericalToRectangular(n.longitude,                   n.latitude);
    Vec3D a = sphericalToRectangular(e.sourceNode.longitude,        e.sourceNode.latitude);
    Vec3D b = sphericalToRectangular(e.destinationNode.longitude,   e.destinationNode.latitude);
    Vec3D o = Vec3D{0,0,0};
    
    Vec3D ao = sub(o,a);
    Vec3D ab = sub(b,a);
    Vec3D ac = sub(c,a);
    Vec3D aoxab = crossProduct(ao,ab);
    double w = dotProduct(aoxab,ac);
    
    //std::cout << w << std::endl;

    return w > 0;
}

bool isNodeLeftOfEdge(double n_long, double n_lat, Edge2 e){
    Vec3D c = sphericalToRectangular(n_long,                   n_lat);
    Vec3D a = sphericalToRectangular(e.sourceLongitude,        e.sourceLatitude);
    Vec3D b = sphericalToRectangular(e.targetLongitude,   e.targetLatitude);
    Vec3D o = Vec3D{0,0,0};
    
    Vec3D ao = sub(o,a);
    Vec3D ab = sub(b,a);
    Vec3D ac = sub(c,a);
    Vec3D aoxab = crossProduct(ao,ab);
    double w = dotProduct(aoxab,ac);
    
    //std::cout << w << std::endl;
    //printEdge(e);

    return w > 0;
}

bool isArcIntersecting(Edge e1, Edge e2){
    bool e1node1L = isNodeLeftOfEdge(e1.sourceNode, e2);
    bool e1node2L = isNodeLeftOfEdge(e1.destinationNode, e2);
    bool e2node1L = isNodeLeftOfEdge(e2.sourceNode, e1);
    bool e2node2L = isNodeLeftOfEdge(e2.destinationNode, e1);
    
    return (e1node1L != e1node2L) && (e2node1L != e2node2L);
}

bool isArcIntersecting(Edge2 e1, Edge2 e2){
    bool e1node1L = isNodeLeftOfEdge(e1.sourceLongitude, e1.sourceLatitude,  e2);
    bool e1node2L = isNodeLeftOfEdge(e1.targetLongitude, e1.targetLatitude,  e2);
    bool e2node1L = isNodeLeftOfEdge(e2.sourceLongitude, e2.sourceLatitude,  e1);
    bool e2node2L = isNodeLeftOfEdge(e2.targetLongitude, e2.targetLatitude,  e1);
    
    return (e1node1L != e1node2L) && (e2node1L != e2node2L);
}

/**
 * TODO: test this method better
 **/
bool isPointInPolygon(Node n, std::vector<Edge> edges){
    Edge firstEdge = edges.at(0);
    bool nLeftOfEdge = isNodeLeftOfEdge(n, firstEdge);

    // calculate halfway vector beetween edge nodes, then normalize and convert to node with lat-long
    // there seems to be an error with center point calculation or more likely with rect to spherical conversion
    Vec3D centerRect = normalize(
        add(
            sphericalToRectangular(
                firstEdge.sourceNode.longitude, firstEdge.sourceNode.latitude
            ),
            sphericalToRectangular(
                firstEdge.destinationNode.longitude, firstEdge.destinationNode.latitude)
        ));
    Vec3D centerSpherical = rectangularToSpherical(centerRect);
    Node centerNode = Node{0, centerSpherical.y, centerSpherical.z};
    
    // edge representing arc from centerNode of edge to query point
    Edge intersectionEdge {n, centerNode};
    
    // check for all edges (except the first one) if they intersect the intersection edge and count intersections
    uint64_t intersectionCount = 0;
    for(auto edgeIter = edges.begin()+1; edgeIter!=edges.end(); ++edgeIter){
        intersectionCount += isArcIntersecting(*edgeIter, intersectionEdge);
    }

    std::cout << intersectionCount << "\n";
    bool even = (intersectionCount%2) == 0;

    // (left and !even) OR (!left and even) means the node is outside, otherwise it is inside
    // equivalent to left XOR even
    return nLeftOfEdge == even;
}

// new isPointInPolygon less pointers
bool isPointInPolygon(Node n, std::vector<Edge2> edges){
    Edge2 firstEdge = edges.at(0);
    bool nLeftOfEdge = isNodeLeftOfEdge(n.longitude, n.latitude, firstEdge);

    // calculate halfway vector beetween edge nodes, then normalize and convert to node with lat-long
    Vec3D centerRect = normalize(
        add(
            sphericalToRectangular(
                firstEdge.sourceLongitude, firstEdge.sourceLatitude
            ),
            sphericalToRectangular(
                firstEdge.targetLongitude, firstEdge.targetLatitude)
        ));
    Vec3D centerSpherical = rectangularToSpherical(centerRect);
    Node centerNode = Node{0, centerSpherical.y, centerSpherical.z};
    
    // edge representing arc from centerNode of edge to query point
    Edge2 intersectionEdge {n.longitude, n.latitude, centerNode.longitude, centerNode.latitude};
    
    // check for all edges (except the first one) if they intersect the intersection edge and count intersections
    uint64_t intersectionCount = 0;
    for(auto edgeIter = edges.begin()+1; edgeIter!=edges.end(); ++edgeIter){
        intersectionCount += isArcIntersecting(*edgeIter, intersectionEdge);
    }

    std::cout << intersectionCount << "\n";
    bool even = (intersectionCount%2) == 0;

    // (left and !even) OR (!left and even) means the node is outside, otherwise it is inside
    // equivalent to left XOR even
    return nLeftOfEdge == even;
}


void load_ploygon_edges(std::string load_string, std::vector<Edge2> &edges){

    std::ifstream textfile;
    textfile.open(load_string, std::ios::in);
    //textfile.exceptions(textfile.exceptions() | std::ios::failbit | std::ifstream::badbit);

    double sourceLongitude;
    double sourceLatitude;
    double targetLongitude;
    double targetLatitude;

    while(!textfile.eof()){
        textfile.read(reinterpret_cast<char *>(&sourceLatitude), sizeof(sourceLatitude));
        textfile.read(reinterpret_cast<char *>(&sourceLongitude), sizeof(sourceLongitude));
        textfile.read(reinterpret_cast<char *>(&targetLatitude), sizeof(targetLatitude));
        textfile.read(reinterpret_cast<char *>(&targetLongitude), sizeof(targetLongitude));

        edges.push_back(Edge2 {sourceLongitude, sourceLatitude, targetLongitude, targetLatitude});
    }
}

/**
 * This method shows that conversion does not work for the edge case latitude=-90
 * */
void test_conversion(){
    for(double lo=-180; lo<=180; lo+=0.5){
        for(double la=-89.999; la<=90; la+=0.5){
            Vec3D rect = sphericalToRectangular(lo, la);
            Vec3D sph = rectangularToSpherical(rect);
            double lodiff = sph.y-lo;
            double ladiff = sph.z-la;
            if (lodiff > 0.01 || ladiff > 0.01){
                std::cout << "for lo=" << lo << " and la=" << la << " there are higher diffs\n";
                printVec(sph);
            }
        }
    }
}

/**
 * Check if param Edge intersects window defined by bounds on longitude and latitude (should also detect edges whose end points are not in the cell)
 **/
bool isEdgeInWindow(double longLow, double longHigh, double latLow, double latHigh, Edge2 edge){
    return 
        isArcIntersecting(edge, Edge2{longLow,  latLow,     longLow,    latHigh})
    &&                        isArcIntersecting(edge, Edge2{longLow,    latHigh,    longHigh,   latHigh})
    &&                                                isArcIntersecting(edge, Edge2{longHigh,   latHigh, longHigh,  latLow})
    &&                                                                     isArcIntersecting(edge, Edge2{longHigh,  latLow, longLow, latLow});
}


int main(int argc, char** argv) {

    // if(argc != 2){
    //     std::cout << "Usage: " << argv[0] << " file_to_read.save" << std::endl;
    //     return 1;
    // }




    // check rectangular to spherical
    double lo = 0;
    double la = 0;
    Vec3D rect = sphericalToRectangular(lo, la);
    std::cout << "rectangular coordinate of long = " << lo << ", lat = " << la << 
    " is x=" << rect.x << " y=" << rect.y << " z=" << rect.z << "\n";
    Vec3D back = rectangularToSpherical(rect);
    std::cout << "transformed back" 
    " is r=" << back.x << " long=" << back.y << " lat=" << back.z << "\n";

    // check cross product
    Vec3D res = crossProduct(Vec3D{1,2,3},Vec3D{-7,8,9});
    std::cout << "should be -6, -30, 22\n";
    printVec(res);

    // check sub
    Vec3D sres = sub(Vec3D{1,2,3},Vec3D{-7,8,10});
    std::cout << "should be 8, -6, -7\n";
    printVec(sres);

    // check dot product
    double cres = dotProduct(Vec3D{1,2,3},Vec3D{-7,8,9});
    std::cout << "should be 36, is " << cres << std::endl; 

    
    std::vector<Node> nodes {Node{0,-1,0},Node{1,1,0},Node{2,0,-1},Node{2,0,1}};
    std::vector<Edge> edges {Edge{nodes.at(0),nodes.at(1)},Edge{nodes.at(2),nodes.at(3)}};
    
    // read data from binary file
    std::vector<Edge2> edges2;
    load_ploygon_edges("data/antarctica-edges.save", edges2);

    // check on which side of the edge arc node 2 is
    bool bres = isNodeLeftOfEdge(nodes.at(3), edges.at(0));
    std::cout << "is Left? " << bres << std::endl;

    // check if arcs are intersecting
    bool bres2 = isArcIntersecting(edges.at(0),edges.at(1));
    std::cout << "is Intersecting? " << bres2 << std::endl;

    // testing data for polygon
    std::vector<Node> pNodes{
        Node{0,-5.0,-6.5},Node{1, 1.0,-3.0},Node{2,-5.5,-6.0},Node{3,1,-1}
    };
    std::vector<Edge> polygon{
        Edge{pNodes.at(0),pNodes.at(1)},
        Edge{pNodes.at(1),pNodes.at(2)},
        //Edge{pNodes.at(3),pNodes.at(2)},
        //Edge{pNodes.at(0),pNodes.at(3)},
    };

    // check if point is in polygon
    Node toCheck = Node{0,-2.0,-4.2};
    Node midPointApprox = Node{0, -2, -4.75};
    bool inPoly1 = isPointInPolygon(toCheck, polygon);
    std::cout << "in Polygon? " << inPoly1 << std::endl;
    std::cout << "is left of first? " << isNodeLeftOfEdge(toCheck, polygon.at(0)) << std::endl;
    std::cout << "is left of second? " << isNodeLeftOfEdge(toCheck, polygon.at(1)) << std::endl;
    std::cout << "ist arc intersecting? " << isArcIntersecting(Edge{toCheck, midPointApprox},polygon.at(1)) << std::endl;
    printEdge(polygon.at(0));

    test_conversion();

    printEdge(edges2.at(0));
    bool inPoly2 = isPointInPolygon(toCheck, edges2);
    std::cout << "in Polygon? " << inPoly2 << std::endl;
    std::cout << "is left of first? " << isNodeLeftOfEdge(toCheck.longitude, toCheck.latitude, edges2.at(0)) << std::endl;
} 