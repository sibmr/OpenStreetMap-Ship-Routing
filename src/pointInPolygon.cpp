/**
 * @file pointInPolygon.cpp
 * @author Marcel Galuschka
 * @author Simon Bihlmaier
 * @brief Point in Polygon check for Coastline Polygons
 * @version 0.1
 * @date 2021-06-02
 * 
 * Implementation is inspired by:
 * 
 * title:           Fast and Robust Point-in-Spherical-Polygon Tests Using Multilevel Spherical Grids
 * authors:         Jing Li, Han Zhang, and Wencheng Wang
 * published in:    Next Generation Computer
 *                  Animation Techniques
 *                  Third International Workshop, AniNex 2017
 *                  Bournemouth, UK, June 22–23, 2017
 *                  Revised Selected Papers
 * pages:           56-66
 * 
 */

# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <array>
# include <math.h>
# include <algorithm>
# include <chrono>
# include <omp.h>

const double globalLongLow = -180;
const double globalLongHigh = 180;
const double globalLatLow = -84;
const double globalLatHigh = 90;

struct Node {
    uint64_t id;
    double longitude;
    double latitude;
};

struct Edge {
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

void printEdge(Edge edge){
    std::cout << "source longitude=" << edge.sourceLongitude << " ,latitude=" << edge.sourceLatitude 
    << " ,destination longitude=" << edge.targetLongitude << " ,latitude=" << edge.targetLatitude << "\n";
}

/**
 * @brief calculate 3d position on a sphere with origin as center
 * 
 * @param longitude 
 * @param latitude 
 * @return Vec3D 
 */
Vec3D sphericalToRectangular(double longitude, double latitude){
    double theta = (M_PI*(latitude+90))/180.0;
    double phi = (M_PI*longitude)/180.0;
    double stheta = sin(theta);
    double x = stheta*cos(phi);
    double y = stheta*sin(phi);
    double z = cos(theta);
    return Vec3D{x, y, z};
}

/**
 * @brief 3d position to longitude, latitude on sphere of radius r
 * 
 * @param coords  3d position
 * @return Vec3D  {r, longitude, latitude}
 */
Vec3D rectangularToSpherical(Vec3D coords){
    double r = sqrt(coords.x*coords.x + coords.y*coords.y + coords.z*coords.z);
    double latitude = -atan2(coords.z,sqrt(coords.x*coords.x + coords.y*coords.y))*(180/M_PI);
    double longitude = 0;
    longitude = atan2(coords.y,coords.x)*(180/M_PI);
    return Vec3D{r, longitude, latitude};
}

/**
 * @brief calculate the cross product of two 3d vectors v1, v2
 * 
 * @param v1
 * @param v2 
 * @return Vec3D orthogonal to v1 and v2
 */
Vec3D crossProduct(const Vec3D &v1, const Vec3D &v2){
    return Vec3D{
        v1.y*v2.z - v1.z*v2.y,
        v1.z*v2.x - v1.x*v2.z,
        v1.x*v2.y - v1.y*v2.x
    };
}

/**
 * @brief add two 3d vectors v1, v2
 * 
 * @param v1 
 * @param v2 
 * @return Vec3D sum of v1 and v2
 */
Vec3D add(const Vec3D &v1, const Vec3D &v2){
    return Vec3D{
        v1.x+v2.x,
        v1.y+v2.y,
        v1.z+v2.z
    };
}

/**
 * @brief substract two 3d vectors v1, v2
 * 
 * @param v1 
 * @param v2 
 * @return Vec3D difference of v1 and v2
 */
Vec3D sub(const Vec3D &v1, const Vec3D &v2){
    return Vec3D{
        v1.x-v2.x,
        v1.y-v2.y,
        v1.z-v2.z
    };
}

/**
 * @brief divide 3d vector v by a scalar divisor
 * 
 * @param v 
 * @param divisor 
 * @return Vec3D v divided by divisor
 */
Vec3D div(const Vec3D &v, double divisor){
    return Vec3D{
        v.x/divisor,
        v.y/divisor,
        v.z/divisor
    };
}

/**
 * @brief calculate dot product between 3d vectors v1 and v2
 * 
 * @param v1 
 * @param v2 
 * @return double 
 */
double dotProduct(Vec3D v1, Vec3D v2){
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

/**
 * @brief normalize 3d vector
 * 
 * @param v 
 * @return Vec3D unit-length vector with same direction as v
 */
Vec3D normalize(Vec3D v){
    double norm = std::sqrt(dotProduct(v,v));
    return div(v, norm); 
}

/**
 * @brief check if node at position (n_long, n_lat) is left of Edge / Line Segment e on a sphere
 * 
 * @param n_long 
 * @param n_lat 
 * @param e 
 * @return true     if the node is left of the edge
 * @return false    otherwise
 */
bool isNodeLeftOfEdge(double n_long, double n_lat, const Edge &e){
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

/**
 * @brief Check if two edges or line segments on a sphere are intersecting
 * 
 * If for each edge e, the source and target node of the other edge e' are on opposite sides of e,
 * then the edges have to cross
 * 
 * In addition we have to check if they are more than half of the maximum distance of two points on a
 * sphere apart by calculating the cos(angle) between them
 * 
 * This prevents detecting intersections with edges on the opposite side of the globe
 * 
 * @param e1 
 * @param e2 
 * @return true     if e1 and e2 are intersecting
 * @return false    otherwise
 */
bool isArcIntersecting(Edge e1, Edge e2){
    bool e1node1L = isNodeLeftOfEdge(e1.sourceLongitude, e1.sourceLatitude,  e2);
    bool e1node2L = isNodeLeftOfEdge(e1.targetLongitude, e1.targetLatitude,  e2);
    bool e2node1L = isNodeLeftOfEdge(e2.sourceLongitude, e2.sourceLatitude,  e1);
    bool e2node2L = isNodeLeftOfEdge(e2.targetLongitude, e2.targetLatitude,  e1);
    
    // need to make sure that arc segments are not on opposite sides of the sphere
    double angle = dotProduct(
        sphericalToRectangular(e1.sourceLongitude, e1.sourceLatitude),
        sphericalToRectangular(e2.sourceLongitude, e2.sourceLatitude)
    );

    return (e1node1L != e1node2L) && (e2node1L != e2node2L) && angle > 0;
}

/**
 * @brief check if node n is inside the polygons defined by a list of edges/line segments
 * 
 * One edge e' is selected and we check if the node is left of that edge.
 * Then the number of intersections of all other edges with 
 * the line segmentfrom the midpoint of e' to node n are counted.
 * Finally we can tell, based on the intersection count and if n is left of e'
 * if n is inside a polygon or not.
 * 
 * this check can be used on all coastline edges for initialization of the first partition center
 * 
 * @param n 
 * @param edges 
 * @return true     if n is inside the polygons defined by edges
 * @return false    otherwise
 */
bool isPointInPolygon(Node n, std::vector<Edge> edges){
    Edge firstEdge = edges.at(0);
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
    Edge intersectionEdge {n.longitude, n.latitude, centerNode.longitude, centerNode.latitude};
    
    // check for all edges (except the first one) if they intersect the intersection edge and count intersections
    uint64_t intersectionCount = 0;
    for(auto edgeIter = edges.begin()+1; edgeIter!=edges.end(); ++edgeIter){
        intersectionCount += isArcIntersecting(*edgeIter, intersectionEdge);
    }

    bool even = (intersectionCount%2) == 0;

    // (left and even) OR (!left and !even) means the node is outside, otherwise it is inside
    // equivalent to left == even
    return nLeftOfEdge == even;
}

/**
 * @brief Load the extracted coastline edges at path load_string into a vector of edges
 * 
 * @param load_string 
 * @param edges 
 */
void load_ploygon_edges(std::string load_string, std::vector<Edge> &edges){

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

        edges.push_back(Edge {sourceLongitude, sourceLatitude, targetLongitude, targetLatitude});
    }
}

/**
 * @brief Check if an edge is fully or partly contained in a rectangle on a sphere
 * 
 * @param longLow   | bottom left corner of rectangle
 * @param latLow    |
 * @param longHigh  | bottom right corner of rectangle
 * @param latHigh   |
 * @param edge      
 * @return true     if edge is contained in the window or intersecting the window
 * @return false    otherwise
 */
bool isEdgeInWindow(double longLow, double latLow, double longHigh, double latHigh, Edge &edge){
    return 
        (edge.sourceLongitude > longLow) && (edge.sourceLongitude < longHigh) && (edge.sourceLatitude > latLow) && (edge.sourceLatitude < latHigh)
    ||  (edge.targetLongitude > longLow) && (edge.targetLongitude < longHigh) && (edge.targetLatitude > latLow) && (edge.targetLatitude < latHigh) //;
    ||      isArcIntersecting(edge, Edge{longLow,  latLow,  longLow,   latHigh} )
    ||      isArcIntersecting(edge, Edge{longLow,  latHigh, longHigh,  latHigh} )
    ||      isArcIntersecting(edge, Edge{longHigh, latHigh, longHigh,  latLow}  )
    ||      isArcIntersecting(edge, Edge{longHigh, latLow,  longLow,   latLow}  );
}
// returns (minI, maxI, minJ, maxJ) 
void getBounds(std::array<int, 4> &array, Edge edge, size_t width, size_t height){
    array[0] = std::max(int ((std::min(edge.sourceLongitude,edge.targetLongitude) - globalLongLow) * width/ (globalLongHigh - globalLongLow)), 0);
    array[1] = std::min(int ((std::max(edge.sourceLongitude,edge.targetLongitude) - globalLongLow) * width/ (globalLongHigh - globalLongLow)), int (width-1));
    array[2] = std::max(int ((std::min(edge.sourceLatitude,edge.targetLatitude) - globalLatLow) * height/ (globalLatHigh - globalLatLow)) , 0);
    array[3] = std::min(int ((std::max(edge.sourceLatitude,edge.targetLatitude) - globalLatLow) * height/ (globalLatHigh - globalLatLow)) , int (height-1));
    return;
};

/**
 * @brief for a grid partitioning, fill each partition cell with the corresponding edges
 * 
 * @tparam width    of the partition grid
 * @tparam height   of the partition grid
 * @param edges         input: list of edges
 * @param partitions    output: partition grid
 */
template<std::size_t width, std::size_t height>
void fillPartitions(std::vector<Edge> &edges, std::vector<Edge*> (&partitions)[width][height]){
    const double longStep = (globalLongHigh-globalLongLow)/width;
    const double latStep = (globalLatHigh-globalLatLow)/height;
    
    uint64_t count = 0;
    Edge *candidate;
    std::array<int,4> boundsArray = {0,0,0,0}; // (minI, maxI, minJ, maxJ) 
    for(uint64_t c=0; c<edges.size(); ++c)
    {
        candidate = &edges.at(c);
        getBounds(boundsArray, *candidate, width, height);

        for(int i = boundsArray[0]; i <= boundsArray[1]; i++){
            for(int j = boundsArray[2]; j <= boundsArray[3]; j++)
            {
                if(isEdgeInWindow(globalLongLow+i*longStep, globalLatLow+j*latStep, globalLongLow+(i+1)*longStep, globalLatLow+(j+1)*latStep, *candidate)){
                    partitions[i][j].push_back(candidate);
                }
            }
        }
        count += 1;
        if(count % 10000 == -1)
            std::cout << (count*100)/edges.size() << "\n";
    }

}

/**
 * @brief for each of the partitions, check if the partition center is in a polygon or outside
 * 
 * this is done in a similar fashion to the article that is referenced at the top of this file
 * 
 * @tparam width            of the partitions grid
 * @tparam height           of the partitions grid
 * @param partitions        input: partitions containing the edges
 * @param partitionCenters  ouput: true if on land, false if in ocean
 */
template<std::size_t width, std::size_t height>
void fillPartitionCenters(std::vector<Edge*> (&partitions)[width][height], bool (&partitionCenters)[width][height]){
    const double longStep = (globalLongHigh-globalLongLow)/width;
    const double latStep = (globalLatHigh-globalLatLow)/height;
    
    // assumption: first partition center is on land (antarctica)
    partitionCenters[0][0] = true;

    int c_i, c_j, p_i, p_j;
    Edge centerEdge;
    // go over all partitions except first one
    for(int x=1; x<width*height; ++x)
    {
        c_i = x/height;
        c_j = x%height;
        p_i = (x-1)/height;
        p_j = (x-1)%height;
        
        // if we jump one latitude step up, then set p_i=c_i to make sure that p_i is beneath c_i
        if(c_j < p_j){
            p_i = c_i -1;
            p_j = c_j;
        }
        //std::cout << "(" << p_i << " " << p_j << ") \t (" << c_i << " " << c_j << ")" << std::endl;


        // edge goes from previous cell center to next cell center
        centerEdge.sourceLongitude  = globalLongLow + p_i*longStep + longStep/2;
        centerEdge.sourceLatitude   = globalLatLow  + p_j*latStep  + latStep/2;
        centerEdge.targetLongitude  = globalLongLow + c_i*longStep + longStep/2;
        centerEdge.targetLatitude   = globalLatLow  + c_j*latStep  + latStep/2;
        
        //std::cout << "----- step " << x << " -----\n";
        //std::cout << "c=(" << c_i << ", " << c_j << " )\n";
        //std::cout << "p=(" << p_i << ", " << p_j << " )\n";
        //printEdge(centerEdge);

        //std::cout << "loop1" << "\n";
        std::vector<Edge*> processedEdges;
        uint64_t intersectionCount = 0;
        for(Edge *edge : partitions[p_i][p_j]){
            processedEdges.push_back(edge);
            intersectionCount += isArcIntersecting(*edge, centerEdge);
        }
        std::sort(processedEdges.begin(), processedEdges.end());
        //std::cout << "loop2" << "\n";
        for(Edge *edge : partitions[c_i][c_j]){
            auto lowerEdgeIt = std::lower_bound(processedEdges.begin(), processedEdges.end(), edge);
            // only increment if edge was not yet processed
            Edge *lowerEdge = (lowerEdgeIt == processedEdges.end()) ? nullptr : *lowerEdgeIt;
            intersectionCount += ((lowerEdge != edge) && isArcIntersecting(*edge, centerEdge));
        }

        // Equivalence operator for intersection count even and previous in polygon
        //  even    p_in    c_in
        //  0       0       1
        //  0       1       0
        //  1       0       0
        //  1       1       1

        //std::cout << "(" << p_i << " " << p_j << ") " <<  partitionCenters[p_i][p_j]  << "\t(" << c_i << " " << c_j << ") " << intersectionCount << std::endl;
        partitionCenters[c_i][c_j] = (intersectionCount%2==0) == partitionCenters[p_i][p_j];
    }
}

/**
 * @brief Given the partitions with predetermined partition centers, calculate if a point is on land or in the ocean
 * 
 * First the partition that the query point (longitude, latitude) is determinde
 * Intersections of edges with the line segment from the query node to the partitionCenter are counted
 * Based on the intersection count and if the partitionCenter is on land we can tell if the node is on land or in the ocean
 * 
 * @tparam width    of partition grid
 * @tparam height   of patition grid
 * @param partitions
 * @param partitionCenters
 * @param longitude of query position
 * @param latitude  of query position
 * @return true     if point (longitude, latitude) is on land
 * @return false    otherwise
 */
template<std::size_t width, std::size_t height>
bool queryPartitions(std::vector<Edge*>(&partitions)[width][height], bool (&partitionCenters)[width][height], double longitude, double latitude){
    const double longStep = (globalLongHigh-globalLongLow)/width;
    const double latStep = (globalLatHigh-globalLatLow)/height;
    
    int i = ((longitude - globalLongLow ) /longStep);
    int j = ((latitude  - globalLatLow ) /latStep);


    Edge queryCenterEdge{longitude, latitude, globalLongLow + i*longStep + longStep/2, globalLatLow + j*latStep + latStep/2};
    uint64_t count = 0;
    for(Edge *edge : partitions[i][j]){
        count += isArcIntersecting(*edge, queryCenterEdge);
    }
    return (count%2==0) == partitionCenters[i][j];
}

/**
 * @brief For a grid of points, determine for each point if it is on land or in the ocean
 * 
 * Call query partitions on each grid point
 * This is parallelized using OpenMP
 * 
 * @tparam partition_width      of partition grid
 * @tparam partition_height     of partition grid
 * @tparam grid_width           of query points grid
 * @tparam grid_height          of query points grid
 * @param partitions            of edges
 * @param partitionCenters      wether the partition center is on land or in the ocean
 * @param gridPoints            wether the grid points are on land or in the ocean
 */
template<std::size_t partition_width, std::size_t partition_height, std::size_t grid_width, std::size_t grid_height>
void determineGridPoints(
    std::vector<Edge*>(&partitions)[partition_width][partition_height], 
    bool (&partitionCenters)[partition_width][partition_height], 
    std::array<std::array<bool, grid_height>, grid_width> &gridPoints
    )
{   
    // add half step offset to borders to get full step width when connecting accross map edge
    // 1 spare step to divide at both map borders
    const double longStep = (globalLongHigh-globalLongLow)/(grid_width+1);
    const double latStep = (globalLatHigh-globalLatLow)/(grid_height+1);

    
    //const uint64_t numLongSteps = (uint64_t)(globalLongHigh-globalLongLow-longStep)/longStep;
    //const uint64_t numLatSteps = (uint64_t) (globalLatHigh-globalLatLow-latStep)/latStep;
    const uint64_t numLongSteps = (uint64_t) grid_width;
    const uint64_t numLatSteps = (uint64_t) grid_height;

    std::cout << "Number of queries: " << numLongSteps*numLatSteps << "\n";
    std::cout << "Step size: " << longStep << " " << latStep << "\n";

    std::chrono::duration<double> query_timing;
    auto startQuery = std::chrono::high_resolution_clock::now();
    uint64_t progress_count = 0;

    #pragma omp parallel num_threads(8)
    {
        #pragma omp for
        for(uint64_t i = 0; i<numLongSteps; ++i)
        for(uint64_t j = 0; j<numLatSteps; ++j)
        {
            gridPoints[i][j] = queryPartitions(partitions, partitionCenters, globalLongLow + longStep/2 + i*longStep, globalLatLow + latStep/2 + j*latStep);
            
            // comment this out when no progress tracking is necessary - critical section slows things down just for the sake of incrementing the progress counter
            #pragma omp critical
            {
            if((progress_count%10000)==0){
                std::cout << "Progess: " << (progress_count)/((double)numLongSteps*numLatSteps) << "\n";
                std::cout << "Coordinate: " << globalLongLow + longStep/2 + i*longStep << " " << globalLatLow + latStep/2 + j*latStep << "\n";
                std::cout << "Result: " << gridPoints[i][j] << "\n";
            }
            progress_count++;
            }
        }
    }

    auto endQuery = std::chrono::high_resolution_clock::now();
    query_timing = endQuery - startQuery;
    std::cout << "Total query time: " << query_timing.count() << "seconds" << std::endl;

}


/**
 * @brief Alternative implementation of queryPartitions, currently unused
 * 
 * @tparam width 
 * @tparam height 
 * @param vectorGrid 
 * @param iNodes 
 * @param jNodes 
 */
template<std::size_t width, std::size_t height>
void generateGrid(std::vector<bool> &vectorGrid, std::vector<Edge*>(&partitions)[width][height], bool (&partitionCenters)[width][height], uint64_t iNodes, uint64_t jNodes){
    const double longStep = ((globalLongHigh - globalLongLow) / iNodes);
    const double latStep = ((globalLatHigh - globalLatLow) / jNodes);

    double curr_long;
    double curr_lat;

    uint64_t maxSize = std::min(vectorGrid.size(), iNodes * jNodes);

    const int size = vectorGrid.size(); 
    for (uint64_t i = 0; i < vectorGrid.size(); i++){
        curr_long = globalLongHigh - ((i % iNodes) * longStep) - longStep/2;
        curr_lat = globalLatHigh - ((i / iNodes) * latStep) - latStep/2;
        vectorGrid.at(i) = queryPartitions(partitions, partitionCenters, curr_long, curr_lat);

        if(i % (size/1000) == 0){
            std::cout << i  << "\t" << size << std::endl;
        }
    }
}

/**
 * @brief method for visualizing partition centers (wether they are on land)
 */
template<std::size_t width, std::size_t height>
bool print_partition_centers(bool (&partitionCenters)[width][height]){
    std::cout << "Printing partition centers" << std::endl;
    for(int i = height-1; i >= 0; i--){
        for (int j = 0; j < width; j++){
            std::cout << partitionCenters[j][i] << " ";
        }
        std::cout << std::endl;
    }
    return true;
}

/**
 * @brief method for printing size of partitions
 */
template<std::size_t width, std::size_t height>
bool print_partitions(std::vector<Edge*>(&partitions)[width][height]){
    std::cout << "Printing partitions" << std::endl;
    for(int i = height-1; i >= 0; i--){
        for (int j = 0; j < width; j++){
            std::cout << partitions[j][i].size() << " ";
        }
        std::cout << std::endl;
    }
    return true;
}

/**
 * @brief method for storing the vector version of the grid
 * currently unused
 */
void saveGrid(std::vector<bool> &vectorGrid, uint64_t iNodes, uint64_t jNodes){
    std::ofstream myfile;

    myfile.open("data/export_graph.txt", std::ios::out | std::ios::trunc);
    myfile.exceptions(myfile.exceptions() | std::ios::failbit | std::ifstream::badbit);

    for(uint64_t i = 0; i < vectorGrid.size(); i++){
        myfile << vectorGrid.at(i) << " ";
        if((i+1)%iNodes == 0){                
            myfile << std::endl;
            myfile.flush();
        }
    }
    myfile <<  "\n" << std::endl;
    myfile.flush();
    myfile.close();

}

/**
 * @brief save grid points to disk: intermediate representation produced by the point in polygon check
 * 
 * creates file at specified path - example: data/worldGrid_width_height.save
 * 
 * @tparam grid_width
 * @tparam grid_height 
 * @param path 
 * @param gridPoints 
 */
template<std::size_t grid_width, std::size_t grid_height>
void saveGridPoints(std::string path,  std::array<std::array<bool, grid_height>, grid_width> &gridPoints){
    std::ofstream textfile;
    textfile.open(path, std::ios::out | std::ios::trunc);
    textfile.exceptions(textfile.exceptions() | std::ios::failbit | std::ifstream::badbit);

    uint64_t width = grid_width;
    uint64_t height = grid_height;
    textfile.write(reinterpret_cast<const char*>(&globalLongLow), sizeof(globalLongLow));
    textfile.write(reinterpret_cast<const char*>(&globalLatLow), sizeof(globalLatLow));
    textfile.write(reinterpret_cast<const char*>(&globalLongHigh), sizeof(globalLongHigh));
    textfile.write(reinterpret_cast<const char*>(&globalLatHigh), sizeof(globalLatHigh));
    textfile.write(reinterpret_cast<const char*>(&width), sizeof(width));
    textfile.write(reinterpret_cast<const char*>(&height), sizeof(height));
    for(uint64_t i = 0; i<grid_width; ++i)
    for(uint64_t j = 0; j<grid_height; ++j)
    {
        textfile.write(reinterpret_cast<const char*>(&gridPoints[i][j]), sizeof(gridPoints[i][j]));
        if(((i*grid_height + j)%2000)==0)
            textfile.flush();
    }
}

/**
 * @brief Create a file containing a geojson representation of edges
 * 
 * @param edges 
 */
void saveEdgesGeoJson(std::vector<Edge> edges){

    std::ofstream file;

    file.open("data/readEdges.json", std::ios::out | std::ios::trunc);
    file.exceptions(file.exceptions() | std::ios::failbit | std::ifstream::badbit);

    file <<     "{ \"type\": \"FeatureCollection\",\n";
    file <<     "  \"features\": [\n";

    uint64_t count = 0;
    bool first_way = true;
    //for(Way way : ways){
    for(Edge edge : edges){
        if(!first_way){file << ",";}
        file << "   {\"type\": \"Feature\",\n"
                "    \"geometry\": {\n"
                "       \"type\": \"LineString\",\n"
                "       \"coordinates\": [\n";

        //bool first_node = true;
        file << "[" << edge.sourceLongitude << "," << edge.sourceLatitude << "],\n";
        file << "[" << edge.targetLongitude << "," << edge.targetLatitude << "]\n";

        file << "   ]\n"
                "   },\n"
                "   \"properties\": {}\n"
                "   }";
        
        if(count%1 == 0){file.flush();}
        first_way = false;
        //if(count%100 == 0){break;}
    }
    file <<     "]}\n" << std::endl;

}

/**
 * @brief Create a file containing a geojson representation of the grid
 */
template<std::size_t grid_width, std::size_t grid_height>
void saveGridPointsGeoJson(int idxLongLow, int idxLatLow, int idxLongHigh, int idxLatHigh, 
    std::array<std::array<bool, grid_height>, grid_width> &gridPoints)
{
    const double longStep = (globalLongHigh-globalLongLow)/(grid_width+1);
    const double latStep = (globalLatHigh-globalLatLow)/(grid_height+1);

    //int idxLongLow = (longLow-globalLongLow)/longStep;
    //int idxLatLow = (longLow-globalLongLow)/longStep;

    std::ofstream file;

    file.open("data/readNodes.json", std::ios::out | std::ios::trunc);
    file.exceptions(file.exceptions() | std::ios::failbit | std::ifstream::badbit);

    file <<     "{ \"type\": \"FeatureCollection\",\n";
    file <<     "  \"features\": [\n";

    uint64_t count = 0;
    bool first_way = true;
    for(int i = idxLongLow; i<idxLongHigh; ++i)
    for(int j = idxLatLow; j<idxLatHigh; ++j)
    {
        if(!first_way){file << ",";}
        file << "   {\"type\": \"Feature\",\n"
                "    \"geometry\": {\n"
                "       \"type\": \"MultiPoint\",\n"
                "       \"coordinates\": [\n";

        if(gridPoints[i][j]){
           file << "[" << globalLongLow + longStep/2 + longStep*i  << "," << globalLatLow + latStep/2 + latStep*j << "]\n";
           first_way = false; 
        }
        

        file << "   ]\n"
                "   },\n"
                "   \"properties\": {}\n"
                "   }";
        
        if(count%1 == 0){file.flush();}
    }
    file <<     "]}\n" << std::endl;

}

/**
 * @brief test sphericalToRectangular and rectangularToSpherical methods
 */
void test_conversion(){
    for(double lo=-180; lo<=180; lo+=0.5){
        for(double la=-89.999; la<=90; la+=0.5){
            Vec3D rect = sphericalToRectangular(lo, la);
            Vec3D sph = rectangularToSpherical(rect);
            double lodiff = abs(sph.y-lo);
            double ladiff = abs(sph.z-la);
            if (lodiff > 0.01 || ladiff > 0.01){
                std::cout << "for lo=" << lo << " and la=" << la << " there are higher diffs\n";
                printVec(sph);
            }
        }
    }
}

/**
 * @brief Test methods on synthetic data
 */
void test_synthetic(){
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

    
    //std::vector<Node> nodes {Node{0,-1,0},Node{1,1,0},Node{2,0,-1},Node{2,0,1}};
    //std::vector<Edge> edges {Edge{nodes.at(0),nodes.at(1)},Edge{nodes.at(2),nodes.at(3)}};
    std::vector<Edge> edges {Edge{-1,0,1,0},Edge{0,-1,0,1}};


    // check on which side of the edge arc node 2 is
    bool bres = isNodeLeftOfEdge(edges.at(1).targetLongitude, edges.at(1).targetLatitude, edges.at(0));
    std::cout << "is Left? " << bres << std::endl;

    // check if arcs are intersecting
    bool bres2 = isArcIntersecting(edges.at(0),edges.at(1));
    std::cout << "is Intersecting? " << bres2 << std::endl;

    // testing data for polygon
    std::vector<Node> pNodes{
        Node{0,-5.0,-6.5},Node{1, 1.0,-3.0},Node{2,-5.5,-6.0},Node{3,1,-1}
    };
    std::vector<Edge> polygon{
        Edge{pNodes.at(0).longitude, pNodes.at(0).latitude,    pNodes.at(1).longitude, pNodes.at(1).latitude},
        Edge{pNodes.at(1).longitude, pNodes.at(1).latitude,    pNodes.at(2).longitude, pNodes.at(2).latitude},
    };

    // check if point is in polygon
    Node toCheck = Node{0,-2.0,-4.2};
    Node midPointApprox = Node{0, -2, -4.75};
    bool inPoly1 = isPointInPolygon(toCheck, polygon);
    std::cout << "in Polygon? " << inPoly1 << std::endl;
    std::cout << "is left of first? " << isNodeLeftOfEdge(toCheck.longitude, toCheck.latitude, polygon.at(0)) << std::endl;
    std::cout << "is left of second? " << isNodeLeftOfEdge(toCheck.longitude, toCheck.latitude, polygon.at(1)) << std::endl;
    std::cout << "is arc intersecting? " << isArcIntersecting(
        Edge{toCheck.longitude, toCheck.latitude, midPointApprox.longitude, midPointApprox.latitude},polygon.at(1)
        ) << std::endl;
    printEdge(polygon.at(0));

    std::vector<Edge> edges2 {
        Edge{10,10,31,-25},
        Edge{25,25,10,10},
        Edge{50,5,25,25},
        Edge{31,-25,50,5},
        //Edge{5,5,6,6},
    };
    std::vector<Edge*> partitions[18][9];
    bool partitionCenters[18][9];
    fillPartitions(edges2, partitions);
    fillPartitionCenters(partitions, partitionCenters);
    print_partitions(partitions);
    print_partition_centers(partitionCenters);

    std::cout << "Test query \n";
    std::cout << "Query ( 31 ,-25.1 ) = " << queryPartitions(partitions, partitionCenters, 31,-25.1) << " =!= 0\n";
    std::cout << "Query ( 50 , 7    ) = " << queryPartitions(partitions, partitionCenters, 50, 7) << " =!= 0\n";
    std::cout << "Query ( 16 , 6    ) = " << queryPartitions(partitions, partitionCenters, 16, 6) << " =!= 1\n";
    std::cout << "Query ( 35 , 8    ) = " << queryPartitions(partitions, partitionCenters, 35, 8) << " =!= 1\n";
    std::cout << "Query ( 36 ,-17   ) = " << queryPartitions(partitions, partitionCenters, 35,-17) << " =!= 1\n";
    std::cout << "Query ( 31 ,-24.9 ) = " << queryPartitions(partitions, partitionCenters, 31,-24.9) << " =!= 1\n";

}

/**
 * @brief Test methods on data non-synthetic data
 * 
 * @param inputFile coastline file
 */
void test_sampleData(std::string inputFile){
    // read data from binary file
    std::vector<Edge> edges2;
    load_ploygon_edges(inputFile, edges2);
    std::cout << "loading done\n";
    Node toCheck = Node{0,-2.0,-4.2};
    std::cout << "is left of first? " << isNodeLeftOfEdge(toCheck.longitude, toCheck.latitude, edges2.at(0)) << std::endl;
    std::vector<Edge*> partitions[201][107];
    bool partitionCenters[201][107];
    std::array<std::array<bool, 1415>, 707> gridPoints;
    fillPartitions(edges2, partitions);
    std::cout << "fill partition centers .." << std::endl;
    fillPartitionCenters(partitions, partitionCenters);
    print_partitions(partitions);
    print_partition_centers(partitionCenters);
    determineGridPoints(partitions, partitionCenters, gridPoints);
    saveGridPointsGeoJson(250, 150, 300, 200, gridPoints);
}

/**
 * @brief reads a .save file at path containing coastlines generated by the coastlineExtraction 
 * and determines for the grid nodes wether they are on land or in the ocean
 * 
 * @tparam grid_width 
 * @tparam grid_height 
 * @param path 
 * @param gridPoints 
 */
template<std::size_t grid_width, std::size_t grid_height>
void prepareGridNodes(std::string path, std::array<std::array<bool, grid_height>, grid_width> &gridPoints){
    std::vector<Edge> edges;
    load_ploygon_edges(path, edges);
    std::cout << "loading done\n";
    std::vector<Edge*> partitions[201][107];
    bool partitionCenters[201][107];
    fillPartitions(edges, partitions);
    fillPartitionCenters(partitions, partitionCenters);
    print_partitions(partitions);
    print_partition_centers(partitionCenters);
    determineGridPoints(partitions, partitionCenters, gridPoints);
    saveGridPointsGeoJson(250, 150, 300, 200, gridPoints);
}

void saveWorldGridPoints_tenmil(std::string inputFile, std::string outputFile){
    const size_t width = 4472, height = 2236;
    std::array<std::array<bool,height>,width>  *gridPoints = new std::array<std::array<bool,height>,width>; 
    std::cout << sizeof(bool)*width*height << "\n";
    prepareGridNodes(inputFile, *gridPoints);
    saveGridPoints(outputFile, *gridPoints);

    delete[] gridPoints;
}
void saveWorldGridPoints_onemil(std::string inputFile, std::string outputFile){
    const size_t width = 1415, height = 707;
    std::array<std::array<bool,height>,width>  *gridPoints = new std::array<std::array<bool,height>,width>; 
    std::cout << sizeof(bool)*width*height << "\n";
    prepareGridNodes(inputFile, *gridPoints);
    saveGridPoints(outputFile, *gridPoints);

    delete[] gridPoints;
}
void saveWorldGridPoints_hundredthousand(std::string inputFile, std::string outputFile){
    const size_t width = 420, height = 240;
    std::array<std::array<bool,height>,width>  *gridPoints = new std::array<std::array<bool,height>,width>; 
    std::cout << sizeof(bool)*width*height << "\n";
    prepareGridNodes(inputFile, *gridPoints);
    saveGridPoints(outputFile, *gridPoints);

    delete[] gridPoints;
}

int main(int argc, char** argv) {
    std::string inputFileName;
    std::string outputFileName;

    if(argc < 1 || argc > 5)
    {
        std::cout << "Usage: " << argv[0] << " file_to_read.coastline" << " " << "file_to.grid" << std::endl;
        return 1;
    }

    int inputFileId = -1;
    int outputFileId = -1;
    int gridSize = 0; // 0 -> 1M, 1 -> 10M, -1 -> 100K
    std::string gridSizeString = "1M";
    // read arguments
    if(argc > 1){
        if(std::string(argv[1])  == "-t"){
            if(argc > 2){
                inputFileId = 2;
            }
        }else if(std::string(argv[1]) == "-n"){
            if(argc > 2){
                if(std::string(argv[2]) == "1M"){
                    gridSize = 0;
                    gridSizeString = "1M";
                }else if(std::string(argv[2]) == "10M"){
                    gridSize = 1;
                    gridSizeString = "10M";
                }else if(std::string(argv[2]) == "100K"){
                    gridSize = -1;
                    gridSizeString = "100K";
                }else{
                    std::cout << "only \"100K\", \"1M\" and \"10M\" nodes are supported" << std::endl;
                    return 1;
                }
            }
            if(argc > 3){
                inputFileId = 3;
            }
            if(argc > 4){
                outputFileId = 4;
            }
        }else{
            if(argc > 1){
                inputFileId = 1;
            }
            if(argc > 2){
                outputFileId = 2;
            }
        }
    }

    // generate input filename
    if(inputFileId == -1){
        inputFileName = "data/planet.coastline";
        std::cout << "no input file given assume " <<  inputFileName << std::endl;
    }else{
        inputFileName = std::string(argv[inputFileId]);
        std::cout << "Use inputfile: " <<  inputFileName << std::endl;
    }

    // generate output filename
    std::string tmp_oFile;
    if(outputFileId == -1){
        tmp_oFile = inputFileName;
    }else{
        tmp_oFile = argv[outputFileId];
    }
    size_t lastindex = tmp_oFile.find_last_of(".");
    outputFileName = tmp_oFile.substr(0, lastindex);
    if(gridSize != 0){
        outputFileName += "_" + gridSizeString;
    }
    outputFileName += ".grid";

    // test mode 
    if(argc > 1){
        if(std::string(argv[1]) == "-t"){    
            std::cout << "Test mode" << std::endl;
            test_conversion();
            test_synthetic();
            std::ifstream ifile(inputFileName);
            if(ifile){
                test_sampleData(inputFileName);
            }else{
                std::cout << "test file not found" << std::endl;
            }
            return 0;
        }
    }
    if(gridSize == 0){
        saveWorldGridPoints_onemil(inputFileName, outputFileName);
    }else if(gridSize == 1){
        saveWorldGridPoints_tenmil(inputFileName, outputFileName);
    }else if(gridSize == -1){
        saveWorldGridPoints_hundredthousand(inputFileName, outputFileName);
    }
    return 0;
} 
