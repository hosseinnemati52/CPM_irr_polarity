#include <iostream>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <string>
#include <string.h>
#include <vector>
#include <sstream>
#include "Mesh.h"
#include "voro++.hh"
#include "mt19937.h"

using namespace std;



RandomGenerator::RandomGenerator()
{
    setseed();
}

void RandomGenerator::setseed()
{
    dsfmt_seed(time(NULL));
}

int RandomGenerator::random(int min, int max)
{
    return (int)((dsfmt_genrand() * ((double)max - min)) + min);
}

double RandomGenerator::random(double min, double max)
{
    return ((dsfmt_genrand() * ((double)max - min)) + min);
}

void Mesh::buildSquareLattice(int nx, int ny, double d)
{
    int ID;
    nPoints = nx * ny;
    lx = nx * d;
    ly = ny * d;

    ID = 0;
    for (int col = 0; col < nx; col++){
        for (int row = 0; row < ny; row++){
            MeshPoints.emplace_back();
            MeshPoints[ID].xCoord = col * d + 0.5 * d;
            MeshPoints[ID].yCoord = row * d + 0.5 * d;
            MeshPoints[ID].meshIndex = ID;
            ID++;
        }
    }
}

void Mesh::buildRead(double Lx, double Ly, int n, const char* fileName)
{

    // cout<<"Hossein Test"<<endl;
    lx = Lx;
    ly = Ly;
    nPoints = n;

    std::ifstream input_file(fileName);
    if (!input_file)
    {
        std::cerr << "Error opening file: " << fileName << std::endl;
        exit(1);
    }

    std::vector<std::vector<double>> data;
    std::string line;
    while (getline(input_file, line))
    {
        std::vector<double> row;
        std::stringstream ss(line);
        double value;
        while (ss >> value)
        {
        row.push_back(value);
        if (ss.peek() == ' ')
        {
            ss.ignore();
        }
        }
        data.push_back(row);
    }

    input_file.close();

    MeshPoints.clear();
    for (int i = 0; i < nPoints; i++){
        MeshPoints.emplace_back();
        MeshPoints[i].meshIndex = i;
        MeshPoints[i].xCoord = data[i][0];
        MeshPoints[i].yCoord = data[i][1];
    }


}

void Mesh::buildIO(double Lx, double Ly, const vector<vector<double>>&xyCom_noVoid)
{

    // cout<<"Hossein Test"<<endl;
    lx = Lx;
    ly = Ly;
    nPoints = xyCom_noVoid.size();
    MeshPoints.clear();
    for (int i = 0; i < nPoints; i++){
        MeshPoints.emplace_back();
        MeshPoints[i].meshIndex = i;
        MeshPoints[i].xCoord = xyCom_noVoid[i][0];
        MeshPoints[i].yCoord = xyCom_noVoid[i][1];
    }


}

void Mesh::buildRandom(double Lx, double Ly, int n)
{
    RandomGenerator Random;
    lx = Lx;
    ly = Ly;
    nPoints = n;

    for (int i = 0; i < nPoints; i++){
        MeshPoints.emplace_back();
        MeshPoints[i].meshIndex = i;
        MeshPoints[i].xCoord = Random.random(0.0, lx);
        MeshPoints[i].yCoord = Random.random(0.0, ly);
    }
}

double Mesh::distance(int id1, int id2, double Lx, double Ly)
{
    double dx, dy;
    dx = MeshPoints[id1].xCoord - MeshPoints[id2].xCoord;
    if (fabs(dx)>Lx/2)
        {dx = fabs(MeshPoints[id1].xCoord - MeshPoints[id2].xCoord)-Lx;}
    dy = MeshPoints[id1].yCoord - MeshPoints[id2].yCoord;
    if (fabs(dy)>Ly/2)
        {dy = fabs(MeshPoints[id1].yCoord - MeshPoints[id2].yCoord)-Ly;}

    return sqrt(dx * dx + dy * dy);
}

void Mesh::generateVornoi(bool plot)
{
    if (MeshPoints.size() == 0) return;

    int nx, ny, nz;     //Number of computational blocks for voro++ routine
    voro::pre_container PreCon(0.0, lx, 0.0, ly, 0.0, 1.0, true, true ,false);

    for (int i = 0; i < MeshPoints.size(); i++){
        PreCon.put(i, MeshPoints[i].xCoord, MeshPoints[i].yCoord, 0.0);
    }
    PreCon.guess_optimal(nx, ny, nz);

    voro::container Con(0.0, lx, 0.0, ly, 0.0, 1.0, nx, ny, nz, true, true, false, 8);
    PreCon.setup(Con);

    voro::c_loop_all Cloop(Con);
    voro::voronoicell_neighbor Cell;
    std::vector<int> NeighborList;
    std::vector<double> FaceAreas;
    int id;

    if (Cloop.start()) do if (Con.compute_cell(Cell, Cloop)){
        Cell.neighbors(NeighborList);
        Cell.face_areas(FaceAreas);
        id = Cloop.pid();
        MeshPoints[id].voronoiVolume = Cell.volume();

        for (int i = 0; i < NeighborList.size(); i++){
            if (NeighborList[i] > -1){ //Is an actual neighbor, not a face due to edge of volume
                MeshPoints[id].neighbors.push_back(NeighborList[i]);
                MeshPoints[id].voronoiWalls.push_back(FaceAreas[i]);
                MeshPoints[id].neighborDist.push_back(distance(id, NeighborList[i], lx, ly));
            }
        }
    } while (Cloop.inc());

    if (plot){
        Con.draw_particles("points_p.gnu");
        Con.draw_cells_gnuplot("points_v.gnu");
    }

}

unsigned int Mesh::forceNeighborSymmetry()
{
	bool neighborMatch;
	unsigned int neighborMisses = 0;
	double faceArea;

	int id, idNeighbor;
	for (id = 0; id < nPoints; id++){
		for (int i = 0; i < MeshPoints[id].neighbors.size(); i++){
			idNeighbor = MeshPoints[id].neighbors[i];
			faceArea = MeshPoints[id].voronoiWalls[i];
			neighborMatch = false;
			for (int j = 0; j < MeshPoints[idNeighbor].neighbors.size(); j++){
				if (MeshPoints[idNeighbor].neighbors[j] == id){
					neighborMatch = true;
					break;
				}

			}

			if (!neighborMatch){ //Add missing particle id to particle idNeighbor
				MeshPoints[idNeighbor].neighbors.push_back(id);
				MeshPoints[idNeighbor].voronoiWalls.push_back(faceArea);
				MeshPoints[idNeighbor].neighborDist.push_back(distance(idNeighbor, id, lx, ly));
				neighborMisses++;
			}
		}
	}

	return neighborMisses;
}

void Mesh::save(const char* fileName)
{
    FILE* f = fopen(fileName, "w");

    fprintf(f, "%lu\t%lf\t%lf\n", MeshPoints.size(), lx, ly);
    for (int i = 0; i < MeshPoints.size(); i++){
        fprintf(f, "%lf\t%lf\n", MeshPoints[i].xCoord, MeshPoints[i].yCoord);
    }
    fclose(f);
}

bool Mesh::load(const char* fileName)
{
    std::ifstream f;
    f.open(fileName);

    if (f.is_open()){
        f >> nPoints;
        f >> lx;
        f >> ly;
        for (int i = 0; i < nPoints; i++){
            MeshPoints.emplace_back();
            MeshPoints[i].meshIndex = i;
            f >> MeshPoints[i].xCoord;
            f >> MeshPoints[i].yCoord;
        }

        f.close();
        return true;
    }

    return false;
}

void Mesh::saveNeighbors(const char* fileName)
{
    FILE* f = fopen(fileName, "w");

    for (int i = 0; i < MeshPoints.size(); i++){
        fprintf(f, "%d\t", MeshPoints[i].meshIndex);
        fprintf(f, "%lf\t%lf\t", MeshPoints[i].xCoord, MeshPoints[i].yCoord);
        for (int j = 0; j < MeshPoints[i].neighbors.size(); j++){
            fprintf(f, "%i\t", MeshPoints[i].neighbors[j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void Mesh::saveVolumes(const char* fileName)
{
    FILE* f = fopen(fileName, "w");

    for (int i = 0; i < MeshPoints.size(); i++){
        fprintf(f, "%d\t%lf\t", MeshPoints[i].meshIndex, MeshPoints[i].voronoiVolume);
        fprintf(f, "\n");
    }
    fclose(f);
}

void Mesh::saveWalls(const char* fileName)
{
    FILE* f = fopen(fileName, "w");

    for (int i = 0; i < MeshPoints.size(); i++){
        fprintf(f, "%d\t", MeshPoints[i].meshIndex);
        fprintf(f, "%lf\t%lf\t", MeshPoints[i].xCoord, MeshPoints[i].yCoord);
        for (int j = 0; j < MeshPoints[i].voronoiWalls.size(); j++){
            fprintf(f, "%lf\t", MeshPoints[i].voronoiWalls[j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void Mesh::saveDists(const char* fileName)
{
    FILE* f = fopen(fileName, "w");

    for (int i = 0; i < MeshPoints.size(); i++){
        fprintf(f, "%d\t", MeshPoints[i].meshIndex);
        fprintf(f, "%lf\t%lf\t", MeshPoints[i].xCoord, MeshPoints[i].yCoord);
        for (int j = 0; j < MeshPoints[i].neighborDist.size(); j++){
            fprintf(f, "%lf\t", MeshPoints[i].neighborDist[j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void Mesh::saveVertices(const char* fileName, double Lx, double Ly)
{
    FILE* f = fopen(fileName, "w");
    int i, j, k;
    int j_counter , k_counter;
    double x_1, y_1, x_2, y_2, x_3, y_3, v_x, v_y;
    double numerator, denominator;
    int vertex_in_set;
    double epsilon = 1.0e-7;
    
    std::vector<double> vertex_x;
    std::vector<double> vertex_y;

    for (int i = 0; i < MeshPoints.size(); i++)
    {
        fprintf(f, "%d:\n", MeshPoints[i].meshIndex);        
        
        vertex_x.clear();
        vertex_y.clear();

        for (int j_counter = 0; j_counter < MeshPoints[i].neighbors.size(); j_counter++){
            
            j = MeshPoints[i].neighbors[j_counter];

            for (int k_counter = 0; k_counter < MeshPoints[j].neighbors.size(); k_counter++)
            {
                k = MeshPoints[j].neighbors[k_counter];
                
                if (k==i)
                {
                    continue;
                }
                else if (std::count(MeshPoints[i].neighbors.begin(), MeshPoints[i].neighbors.end(), k) == 0)
                {
                    continue;
                }
                else
                {
                    x_1 = MeshPoints[i].xCoord;
                    y_1 = MeshPoints[i].yCoord;

                    x_2 = MeshPoints[j].xCoord;
                    if ( x_2-x_1 >  Lx/2)
                        {x_2 = MeshPoints[j].xCoord-Lx;}
                    if ( x_2-x_1 < -Lx/2)
                        {x_2 = MeshPoints[j].xCoord+Lx;}
                    y_2 = MeshPoints[j].yCoord;
                    if ( y_2-y_1 >  Ly/2)
                        {y_2 = MeshPoints[j].yCoord-Ly;}
                    if ( y_2-y_1 < -Ly/2)
                        {y_2 = MeshPoints[j].yCoord+Ly;}

                    x_3 = MeshPoints[k].xCoord;
                    if ( x_3-x_2 >  Lx/2)
                        {x_3 = MeshPoints[k].xCoord-Lx;}
                    if ( x_3-x_2 < -Lx/2)
                        {x_3 = MeshPoints[k].xCoord+Lx;}
                    y_3 = MeshPoints[k].yCoord;
                    if ( y_3-y_2 >  Ly/2)
                        {y_3 = MeshPoints[k].yCoord-Ly;}
                    if ( y_3-y_2 < -Ly/2)
                        {y_3 = MeshPoints[k].yCoord+Ly;}
                    
        

                    if (fabs(y_1-y_2)>epsilon && fabs(y_2-y_3)>epsilon)
                    {
                        numerator = (y_3-y_1)/2. + (x_2*x_2-x_1*x_1)/(2.*(y_1-y_2)) - (x_3*x_3-x_2*x_2)/(2.*(y_2-y_3));
                        denominator = (x_2-x_1)/(y_1-y_2) - (x_3-x_2)/(y_2-y_3);

                        v_x = numerator/denominator;
                        v_y = (y_1+y_2)/2. + ((x_2-x_1)/(y_1-y_2)) * (v_x - (x_1+x_2)/2.);
                    } else
                    {
                        
                        double X_1, Y_1, X_2, Y_2, X_3, Y_3;
                        if (fabs(y_1-y_2)<=epsilon)
                        {
                            X_1 = x_1;
                            Y_1 = y_1;
                            X_2 = x_3;
                            Y_2 = y_3;
                            X_3 = x_2;
                            Y_3 = y_2;
                        }
                        else if(fabs(y_2-y_3)<=epsilon)
                        {
                            X_1 = x_2;
                            Y_1 = y_2;
                            X_2 = x_1;
                            Y_2 = y_1;
                            X_3 = x_3;
                            Y_3 = y_3;
                        }

                        numerator = (Y_3-Y_1)/2. + (X_2*X_2-X_1*X_1)/(2.*(Y_1-Y_2)) - (X_3*X_3-X_2*X_2)/(2.*(Y_2-Y_3));
                        denominator = (X_2-X_1)/(Y_1-Y_2) - (X_3-X_2)/(Y_2-Y_3);

                        v_x = numerator/denominator;
                        v_y = (Y_1+Y_2)/2. + ((X_2-X_1)/(Y_1-Y_2)) * (v_x - (X_1+X_2)/2.);
                    }


                    vertex_in_set = 0;
                    if (vertex_x.size()>0)
                    {
                        for (int v_counter=0; v_counter<vertex_x.size(); v_counter++)
                        {
                            if ((fabs(v_x-vertex_x[v_counter])<epsilon) &&  (fabs(v_y-vertex_y[v_counter])<epsilon))
                            {
                                vertex_in_set = 1;
                                break;
                            }
                        }
                    }
                    if (vertex_in_set == 0)
                    {
                        vertex_x.push_back(v_x);
                        vertex_y.push_back(v_y);
                    }
                }
            }
        }
        for (int v_counter=0; v_counter<vertex_x.size(); v_counter++)
        {
            fprintf(f, "%lf\t%lf\t", vertex_x[v_counter], vertex_y[v_counter]);
            fprintf(f, "\n");
        }
        fprintf(f, "*\n");
    }
    fclose(f);
}