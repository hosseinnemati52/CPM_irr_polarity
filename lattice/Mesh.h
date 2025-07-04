#pragma once
#include <vector>

class RandomGenerator{
public:
    /*----------------Constructor/Destructor---------------*/
    RandomGenerator();

    /*----------------------Methods-------------------------*/
private:
    void setseed();

public:
    /*Generate a random number between in [min, max)*/
    int random(int min, int max);
    double random(double min, double max);
};

struct MeshPoint{
    int meshIndex = -1;                 //Index of this point in the mesh
    double xCoord = 0.0;
    double yCoord = 0.0;
    double voronoiVolume = 0.0;         //The volume of the voronoi cell (area in 2D)
    std::vector<int> neighbors;         //List containing the indices in mesh of the neighbors
    std::vector<double> voronoiWalls;   //The area of the voro cell walls in between neighbors, same indexing as neighbors list
    std::vector<double> neighborDist;   //The distance between neighbors
};

class Mesh{
public:
    int nPoints = -1;           //Number of sites in the mesh
    double lx = 0.0, ly = 0.0;  //Sidelengths of the box
    std::vector<MeshPoint> MeshPoints;


    /*--------------------Constructors/Destructors-----------------*/



    /*---------------------------Methods---------------------------*/

    /*Build a square lattice of size nx * ny points and lattice spacing d*/
    void buildSquareLattice(int nx, int ny, double d);

    /*Reads the data and builds a lattice based on it*/
    void buildRead(double Lx, double Ly, int n, const char* fileName);

    /*Build a random lattice in box of size lx * ly containing n points*/
    void buildRandom(double Lx, double Ly, int n);

    /*Returns the euclidean distance between two sites*/
    double distance(int id1, int id2, double Lx, double Ly);

    /*Generate voronoi cells*/
    void generateVornoi(bool plot);

    /*Make sure if A is neighbour of B, that B is also neighbor of A*/
    unsigned int forceNeighborSymmetry();

    /*Saves mesh to textfile*/
    void save(const char* fileName);

    /*Save mesh and neighbors to textfile*/
    void saveNeighbors(const char* fileName);

    
    void saveVolumes(const char* fileName);

    void saveWalls(const char* fileName);

    void saveDists(const char* fileName);

    void saveVertices(const char* fileName, double Lx, double Ly);

    /*Loads mesh from textfile
    returns true if succesful, false otherwise*/
    bool load(const char* fileName);
};
