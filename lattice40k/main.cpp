#include <stdio.h>
#include "Mesh.h"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;


int main()
{
    Mesh Test;

    double Lx, Ly;
    int NSites;

    // cout<<"Lx: ";
    // cin>>Lx;
    // cout<<"Ly: ";
    // cin>>Ly;
    // cout<<"NumCells: ";
    // cin>>NumCells;
    
    std::ifstream LxLy_file("LxLy.dat");
    if (!LxLy_file)
    {
        std::cerr << "Error opening file: " << "LxLy.dat" << std::endl;
        exit(1);
    }
    std::string line;
    for (int l=0; getline(LxLy_file, line); l++)
    {
        std::stringstream ss(line);
        if (l==0)
        {
            ss >> Lx;
        }
        if (l==1)
        {
            ss >> Ly;
        }
    }
    std::ifstream NSites_file("NSites.dat");
    if (!NSites_file)
    {
        std::cerr << "Error opening file: " << "NSites.dat" << std::endl;
        exit(1);
    }
    getline(NSites_file, line);
    std::stringstream ss(line);
    ss >> NSites;



    Test.buildRead(Lx, Ly, NSites, "test.dat");

    // Test.buildRandom(10.0, 10.0, 100);
    Test.save("lattice.dat");
    Test.generateVornoi(true);
    printf("neighbor mismatches: %u\n", Test.forceNeighborSymmetry());

    Test.saveNeighbors("neighs.dat");

    //By Hossein
    // Test.saveV
    Test.saveVolumes("vols.dat");

    Test.saveWalls("walls.dat");

    Test.saveDists("dists.dat");

    Test.saveVertices("vertices.dat", Lx, Ly);

    system("python3 voroPlotPy_data.py");
    system("python3 voroPlotPy_0.py");
    return 1;
}
