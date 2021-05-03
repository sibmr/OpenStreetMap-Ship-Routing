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
    double sourceLatitude;
    double sourceLongitude;
    double targetLatitude;
    double targetLongitude;
};

struct Vec3D {
    double x;
    double y;
    double z;
};

void printVec(Vec3D vec){
    std::cout << "x=" << vec.x << " y=" << vec.y << " z=" << vec.z << "\n";
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
    double theta = acos(coords.z/r) - M_PI_2;
    double phi = 0;
    if(coords.x != 0){
        double phi = atan(coords.y/coords.x);
    }
    return Vec3D{r, phi, theta};
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
    
    std::cout << w << std::endl;

    return w < 0;
}
bool isNodeLeftOfEdge(double n_lat, double n_long, Edge2 e){
    Vec3D c = sphericalToRectangular(n_long,                   n_lat);
    Vec3D a = sphericalToRectangular(e.sourceLongitude,        e.sourceLatitude);
    Vec3D b = sphericalToRectangular(e.targetLongitude,   e.targetLatitude);
    Vec3D o = Vec3D{0,0,0};
    
    Vec3D ao = sub(o,a);
    Vec3D ab = sub(b,a);
    Vec3D ac = sub(c,a);
    Vec3D aoxab = crossProduct(ao,ab);
    double w = dotProduct(aoxab,ac);
    
    std::cout << w << std::endl;

    return w < 0;
}

bool isArcIntersecting(Edge e1, Edge e2){
    bool e1node1L = isNodeLeftOfEdge(e1.sourceNode, e2);
    bool e1node2L = isNodeLeftOfEdge(e1.destinationNode, e2);
    bool e2node1L = isNodeLeftOfEdge(e2.sourceNode, e1);
    bool e2node2L = isNodeLeftOfEdge(e2.destinationNode, e1);
    
    return (e1node1L != e1node2L) && (e2node1L != e2node2L);
}

bool isArcIntersecting(Edge2 e1, Edge2 e2){
    bool e1node1L = isNodeLeftOfEdge(e1.sourceLatitude, e1.sourceLongitude, e2);
    bool e1node2L = isNodeLeftOfEdge(e1.targetLatitude, e1.targetLongitude, e2);
    bool e2node1L = isNodeLeftOfEdge(e2.sourceLatitude, e2.sourceLongitude, e1);
    bool e2node2L = isNodeLeftOfEdge(e2.targetLatitude, e2.targetLongitude, e1);
    
    return (e1node1L != e1node2L) && (e2node1L != e2node2L);
}

/**
 * TODO: test this method better
 **/
bool isPointInPolygon(Node n, std::vector<Edge> edges){
    Edge firstEdge = edges.at(0);
    bool nLeftOfEdge = isNodeLeftOfEdge(n, firstEdge);

    // calculate halfway vector beetween edge nodes, then normalize and convert to node with lat-long
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
    return nLeftOfEdge != even;
}

// new isPointInPolygon less pointers
bool isPointInPolygon(Node n, std::vector<Edge2> edges){
    Edge2 firstEdge = edges.at(0);
    bool nLeftOfEdge = isNodeLeftOfEdge(n.latitude, n.longitude, firstEdge);

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
    Edge2 intersectionEdge {n.latitude, n.longitude, centerNode.latitude, centerNode.longitude};
    
    // check for all edges (except the first one) if they intersect the intersection edge and count intersections
    uint64_t intersectionCount = 0;
    for(auto edgeIter = edges.begin()+1; edgeIter!=edges.end(); ++edgeIter){
        intersectionCount += isArcIntersecting(*edgeIter, intersectionEdge);
    }

    std::cout << intersectionCount << "\n";
    bool even = (intersectionCount%2) == 0;

    // (left and !even) OR (!left and even) means the node is outside, otherwise it is inside
    // equivalent to left XOR even
    return nLeftOfEdge != even;
}

void load_ploygon_edges(std::string load_string, std::vector<Edge2> &edges){

    std::ifstream textfile;
    textfile.open(load_string, std::ios::in);
    //textfile.exceptions(textfile.exceptions() | std::ios::failbit | std::ifstream::badbit);

    double sourceLatitude;
    double sourceLongitude;
    double targetLatitude;
    double targetLongitude;

    while(!textfile.eof()){
        textfile.read(reinterpret_cast<char *>(&sourceLatitude), sizeof(sourceLatitude));
        textfile.read(reinterpret_cast<char *>(&sourceLongitude), sizeof(sourceLongitude));
        textfile.read(reinterpret_cast<char *>(&targetLatitude), sizeof(targetLatitude));
        textfile.read(reinterpret_cast<char *>(&targetLongitude), sizeof(targetLongitude));

        edges.push_back(Edge2 {sourceLatitude, sourceLongitude, targetLatitude, targetLongitude});
    }
}


int main(int argc, char** argv) {

    if(argc != 2){
        std::cout << "Usage: " << argv[0] << " file_to_read.save" << std::endl;
        return 1;
    }




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
    load_ploygon_edges(argv[1], edges2);

    // check on which side of the edge arc node 2 is
    bool bres = isNodeLeftOfEdge(nodes.at(2), edges.at(0));
    std::cout << "is Left? " << bres << std::endl;

    // check if arcs are intersecting
    bool bres2 = isArcIntersecting(edges.at(0),edges.at(1));
    std::cout << "is Intersecting? " << bres2 << std::endl;

    // testing data for polygon
    std::vector<Node> pNodes{
        Node{0,-1,-1},Node{1,-1,1},Node{2,1,1},Node{3,1,-1}
    };
    std::vector<Edge> polygon{
        Edge{pNodes.at(1),pNodes.at(0)},
        Edge{pNodes.at(2),pNodes.at(1)},
        Edge{pNodes.at(3),pNodes.at(2)},
        Edge{pNodes.at(0),pNodes.at(3)},
    };

    // check if point is in polygon
    Node toCheck = Node{0,0,0};
    bool inPoly = isPointInPolygon(toCheck, edges2);
    std::cout << "in Polygon? " << inPoly << std::endl;
    std::cout << "is left of first? " << isNodeLeftOfEdge(toCheck, edges.at(0)) << std::endl;
} 