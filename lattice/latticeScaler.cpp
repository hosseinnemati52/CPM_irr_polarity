#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>      // std::setprecision

using namespace std;

int main()
{
    double LxOr, LyOr;
    double LxSc, LySc;
    int NSites;

    // cout<<"Lx: ";
    // cin>>Lx;
    // cout<<"Ly: ";
    // cin>>Ly;
    // cout<<"NumCells: ";
    // cin>>NumCells;
    
    std::ifstream LxLyOr_file("LxLyOriginal.dat");
    if (!LxLyOr_file)
    {
        std::cerr << "Error opening file: " << "LxLyOriginal.dat" << std::endl;
        exit(1);
    }
    std::string line;
    for (int l=0; getline(LxLyOr_file, line); l++)
    {
        std::stringstream ss(line);
        if (l==0)
        {
            ss >> LxOr;
        }
        if (l==1)
        {
            ss >> LyOr;
        }
    }

    std::ifstream LxLySc_file("LxLy.dat");
    if (!LxLySc_file)
    {
        std::cerr << "Error opening file: " << "LxLyScaled.dat" << std::endl;
        exit(1);
    }
    for (int l=0; getline(LxLySc_file, line); l++)
    {
        std::stringstream ss(line);
        if (l==0)
        {
            ss >> LxSc;
        }
        if (l==1)
        {
            ss >> LySc;
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

    // cout<<"Hossein Test"<<endl;

    std::ifstream input_file("testOriginal.dat");
    std::ofstream output_file("test.dat");
    if (!input_file)
    {
        std::cerr << "Error opening file: " << "testOriginal.dat" << std::endl;
        exit(1);
    }

    std::vector<std::vector<double>> data;
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
        output_file << std::setprecision(14) << row[0]*(1.0*LxSc/LxOr);
        output_file << ' ';
        output_file << std::setprecision(14)<< row[1]*(1.0*LySc/LyOr);
        output_file << '\n';
    }

    input_file.close();
    output_file.close();

    // for (int i = 0; i < nPoints; i++){
    //     MeshPoints.emplace_back();
    //     MeshPoints[i].meshIndex = i;
    //     MeshPoints[i].xCoord = data[i][0];
    //     MeshPoints[i].yCoord = data[i][1];
    // }


    return 0;
}