#include <iostream>
#include <iomanip>
#include <math.h>
#include <unistd.h>
// #include <array>
#include <vector>
// #include <cstdlib>
// #include <ctime>
// #include <utility>
// #include <tuple>
// #include <cmath>
// #include <map>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
// #include <cstring>
#include <string>
#include <chrono>
#include <sys/stat.h>
// #include <filesystem>
#include <string.h>
#include <unistd.h>
#include <algorithm>
#include "Mesh.h"
#include "voro++.hh"



///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////#DEFINITIONS///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////-
// #define L 200                                         // lattice size
// #define Kb (1.)                                      // Boltzmann constant
// #define Tem (1.)                                     // Temperature
// #define NumCells 1000                                  // Number of cells
// #define AvgCellArea ((double)((L * L) / (NumCells))) // Average area of cells
// #define Lambda (1.)                                  // Area elasticity strength
// #define SweepLength (L * L)                          // Number of attempts per Monte Carlo sweep
// #define samplesPerWrite 10
#define PI (3.14159265359)
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////





using namespace std;
using namespace std::chrono;
// void dataFolderDeleter()
// {
//     std::filesystem::path subfolder_path = "data"; // Replace with actual path
//     try {
//         std::filesystem::remove_all(subfolder_path);
//         std::cout << "Subfolder successfully deleted.\n";
//     } catch (std::filesystem::filesystem_error& e) {
//         std::cerr << "Error deleting subfolder: " << e.what() << '\n';
//     }
// }
bool directoryExists(const std::string& path);

void dataFolderDeleter();

int totalSnapshotsNumFunction(const int samplesPerWrite);

void simulationDataReader_sq(int* LPtr, double* AlphaPtr, double* KbPtr, double* TemPtr,  int* NumCellsPtr, \
                          double* AvgCellAreaPtr, double* LambdaPtr, long* maxMCStepsPtr, int* samplesPerWritePtr, \
                          int* printingTimeIntervalPtr, string* initConfigPtr);

void simulationDataReader_irr(int* NSitesPtr, double* LxPtr, double* LyPtr, double* AlphaPtr, double* KbPtr, double* TemPtr,  int* NumCellsPtr, \
                          double* AvgCellAreaPtr, double* LambdaPtr, long* maxMCStepsPtr, int* samplesPerWritePtr, \
                          int* printingTimeIntervalPtr, int* numLinkedListPtr, string* initConfigPtr);

std::vector<std::string> parseString(const std::string& inputStr);

void timeMaker(long* time, long length, int samplesPerWrite);

void energyMaker(double* energy, long length, int samplesPerWrite);

void sigmaMatLoader(vector<vector<int>>& mat, const std::string &filename, int rows, int cols);

void doubleMatLoader(vector<vector<double>>& mat, const std::string &filename, int rows, int cols);

void configurationalCalcPerSample(vector<int>& areaSample, vector<int>& periSample,\
vector<double>& isoperiSample, vector<double>& xComSample, vector<double>& yComSample, \
const vector<vector<int>>& sigmaMat, const vector<vector<int>>& sitesX, const vector<vector<int>>& sitesY);

void configurationalCalcPerSample_irr(vector<double>& areaSample, vector<double>& periSample,\
vector<double>& isoperiSample, vector<double>& ARSample, vector<double>& circSample, vector<double>& xComSample, vector<double>& yComSample, \
const std::vector<int>& sigmaMat, const std::vector<double>& sitesX, const std::vector<double>& sitesY,\
const std::vector<int>& neighborNum, const std::vector<vector<int>>& neighborsList, const std::vector<vector<double>>& edges,\
const std::vector<double>& latticeArea);

void sigmaMatExtractor(vector<vector<vector<int>>>& mat, std::ifstream& samplesFile);

void printForTest(const vector<vector<int>>& mat);

void doubleMatSaver(double* matPtr, int Rows, int Cols, const std::string &filename);

void intMatSaver(int* matPtr, int Rows, int Cols, const std::string &filename);

void double2DVecSaver(const vector<vector<double>>& mat, const std::string &filename);

void double1DVecSaver(const vector<double>& mat, const std::string &filename);

void saveInt1DVec(const vector<int>& quantity, const std::string &filename);

void saveDbl1DVec(const vector<double>& quantity, const std::string &filename);

void loadInt1DVec(vector<int>& mat, const std::string &filename);

void loadInt1DVecHor(vector<int>& mat, const std::string &filename);

void loadDbl1DVec(vector<double>& mat, const std::string &filename);

void latticeCreatorPhase1(int* maxNeighPtr, vector<int>& neighborNum, vector<double>& latticeX,  \
                          vector<double>& latticeY,  vector<double>& latticeArea );

void latticeCreatorPhase2(const int NSites, const double Lx, const double Ly, const int maxNeighbors,\
                          vector<vector<int>>& neighborsList, vector<vector<double>>& edges,\
                          vector<vector<double>>& vorDist);

void polygonCom(double* xComPtr, double* yComPtr, const vector<vector<double>> polygonVec);

void latticeCreatorPhase3(const int NSites, const double Lx, const double Ly, const int maxNeighbors, \
                          const vector<int> neighborNum, const vector<vector<int>> neighborsList,\
                          vector<double>& siteComX, vector<double>& siteComY,\
                          vector<vector<double>>& deltaComX, vector<vector<double>>& deltaComY, vector<vector<double>>& comDist);

inline bool existFunc (const std::string& name);

void totComFunc(double& xComTot, double& yComTot, const vector<double>& xComSample, const vector<double>& yComSample);

void saveIntMatCSV(const vector<vector<int>>& sigmaMat, const std::string &filename);

void testPlot(Mesh &vorTesTest, const double Lx, const double Ly, const vector<vector<double>>& xyCom_noVoid);

int main()
{

  int NSites;
  double Lx, Ly;
  double Alpha;
  double Kb;
  double Tem;
  int NumCells;
  double AvgCellArea;
  double Lambda;
  int SweepLength;
  long maxMCSteps;
  int samplesPerWrite;
  int printingTimeInterval;
  int numLinkedList;
  string initConfig; // "r" for rectangular; "c" for crystalline;

  double xShift, yShift;
  
  ////////////// IRR OR SQ //////////////////////////////
  if (directoryExists("lattice")) // for irregular lattice
  // if (0) // for irregular lattice
  {
    simulationDataReader_irr(&NSites, &Lx, &Ly, &Alpha, &Kb, &Tem,  &NumCells, \
                      &AvgCellArea, &Lambda, &maxMCSteps, &samplesPerWrite, &printingTimeInterval, &numLinkedList, &initConfig);
    SweepLength = NSites;
    
    xShift = 0.0;
    yShift = 0.0;
  }
  else // for square lattice
  {
    int L;
    simulationDataReader_sq(&L, &Alpha, &Kb, &Tem,  &NumCells, \
                      &AvgCellArea, &Lambda, &maxMCSteps, &samplesPerWrite, &printingTimeInterval, &initConfig);
    SweepLength = L*L;
    Lx = 1.0 * L;
    Ly = 1.0 * L;

    xShift = 0.5;
    yShift = 0.5;
  }
  ////////////// IRR OR SQ //////////////////////////////
  
  ifstream t_w_file;
  int t_w=0;

  


  /*
  ////////////////////////////////////////////////////
  /////////////////// LATTICE READING ////////////////
  ////////////////////////////////////////////////////

  int maxNeighbors;
  vector<int> neighborNum(NSites);
  vector<double> latticeX(NSites);
  vector<double> latticeY(NSites);
  vector<double> latticeArea(NSites);
  latticeCreatorPhase1(&maxNeighbors, neighborNum, latticeX, latticeY, latticeArea);
  
  vector<vector<int>> neighborsList(NSites, vector<int>(maxNeighbors));
  vector<vector<double>> edges(NSites, vector<double>(maxNeighbors));
  vector<vector<double>> vorDist(NSites, vector<double>(maxNeighbors));
  latticeCreatorPhase2(NSites, Lx, Ly, maxNeighbors,\
                      neighborsList, edges, vorDist);

  vector<double> siteComX(NSites); // center of mass of every single site
  vector<double> siteComY(NSites); // center of mass of every single site
  vector<vector<double>> deltaComX(NSites, vector<double>(maxNeighbors)); // this is between the center of mass of sites, not their voronoi point
  vector<vector<double>> deltaComY(NSites, vector<double>(maxNeighbors)); // this is between the center of mass of sites, not their voronoi point
  vector<vector<double>> comDist(NSites, vector<double>(maxNeighbors)); // this is between the center of mass of sites, not their voronoi point
  
  latticeCreatorPhase3(NSites, Lx, Ly, maxNeighbors, neighborNum,\
                      neighborsList, siteComX, siteComY,\
                      deltaComX, deltaComY, comDist);
  
  // saveSampleE(siteComX, "siteComX.csv");
  // saveSampleE(siteComY, "siteComY.csv");
  // saveDoubleMatCSV(deltaComX, "deltaComX.csv");
  // saveDoubleMatCSV(deltaComY, "deltaComY.csv");
  // saveDoubleMatCSV(comDist, "comDist.csv");

  ////////////////////////////////////////////////////
  /////////////////// LATTICE READING ////////////////
  ////////////////////////////////////////////////////
  */


  string ppDataFolderName = "pp_data";
  string ppDataAddress;
  string ppVorFolderName = "pp_vor";
  string run_ = "run_";
  string folderName;
  string folderAddress;

  int bunchC;
  int snapshotC;
  int snapshotCDum;

  vector<vector<double>> xComBunch(NumCells+1, vector<double>(samplesPerWrite));
  vector<vector<double>> yComBunch(NumCells+1, vector<double>(samplesPerWrite));
  vector<vector<int>> nNeighCellBunch(NumCells+1, vector<int>(samplesPerWrite));
  vector<vector<double>> periCellBunch(NumCells+1, vector<double>(samplesPerWrite));
  vector<vector<double>> areaCellBunch(NumCells+1, vector<double>(samplesPerWrite));
  vector<vector<double>> qCellBunch(NumCells+1, vector<double>(samplesPerWrite));

  
  vector<int> nNeighCell(NumCells+1);
  vector<double> periCell(NumCells+1);
  vector<double> areaCell(NumCells+1);
  vector<double> qCell(NumCells+1);

  int cellC=0;
  vector<vector<double>> xyCom(NumCells+1, vector<double>(2));
  vector<vector<double>> xyCom_noVoid(NumCells, vector<double>(2));

  vector<double> xComInit(NumCells+1);
  vector<double> yComInit(NumCells+1);

  /////////////////// MAKING ppData FOLDER /////////////////
  //// This block is for windows:
  // mkdir(ppDataFolderName.c_str()); //making data folder

  Mesh vorTes;

  /////////////// POLYGONS Q VALS /////////////////
  vector<double> polygons_Q; // p/sqrt(A)
  vector<double> polygons_q; // sqrt(4*pi*A)/p
  
  polygons_Q.push_back(0.0);        // n = 0
  polygons_Q.push_back(0.0);        // n = 1
  polygons_Q.push_back(0.0);        // n = 2
  polygons_Q.push_back(4.559);      // n = 3
  polygons_Q.push_back(4.0);        // n = 4
  polygons_Q.push_back(3.8119);     // n = 5
  polygons_Q.push_back(3.7224);     // n = 6
  polygons_Q.push_back(3.67207);    // n = 7
  polygons_Q.push_back(3.640718);   // n = 8
  polygons_Q.push_back(3.6198);     // n = 9

  for (int i = 0; i < polygons_Q.size(); i++)
  {
    if (i < 3)
    {
      polygons_q.push_back(0.0);
    }else
    {
      polygons_q.push_back(sqrt(4.0*PI)/polygons_Q[i]);
    }
    
  }
  /////////////// POLYGONS Q VALS /////////////////




  ////////////// LOOP ON RUNS /////////////////////
  int runC = 1;
  while (1)
  {
    folderName = run_ + to_string(runC);

    if (! directoryExists(folderName)) ////// Stop if you cannot find the run folder
    {
      break;
    }

    folderAddress = folderName + "/" + ppVorFolderName;

    //// This block is for Linux:
    mkdir(folderAddress.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making backup_resume folder  

    ppDataAddress = folderName + "/" + ppDataFolderName;

    ////////// NUMBER OF SNAPSHOTS /////////////////////
    int totalSnapshotsNum;
    ifstream totalSnapshotsNumFile;
    totalSnapshotsNumFile.open(ppDataAddress + "/" + "totalSnapshotsNum.csv");
    totalSnapshotsNumFile >> totalSnapshotsNum;
    totalSnapshotsNumFile.close(); // random seed saved
    ////////// NUMBER OF SNAPSHOTS /////////////////////

    ////////// T_W /////////////////////
    t_w_file.open(folderName + "/" + "t_w.csv");
    t_w_file >> t_w;
    t_w_file.close();
    ////////// T_W /////////////////////

    ////////// TIME /////////////////////
    vector<int> time(totalSnapshotsNum);
    loadInt1DVecHor(time, ppDataAddress + "/time.csv");

    int t_w_ind = 0;
    // while (time[t_w_ind] != t_w)
    while (time[t_w_ind] != t_w)
    {
      t_w_ind++;
    }

    // double dt_indep = 10 * (time[totalSnapshotsNum-1] - time[totalSnapshotsNum-2]);
    vector<int> eq_snapshots_indices;

    // double dt;
    // snapshotC = t_w_ind;
    // while (time[totalSnapshotsNum - 1] - time[snapshotC] >= dt_indep)
    // {
    //   dt = 0.0;
    //   while (dt < dt_indep)
    //   {
    //     dt = dt + time[snapshotC+1] - time[snapshotC];
    //     snapshotC++;
    //   }
    //   eq_snapshots_indices.push_back(snapshotC);        
    // }
    // saveInt1DVec(eq_snapshots_indices, folderAddress + "/" + "eq_snapshots_indices.csv");
    ifstream eqSamplingTimesFile;
    eqSamplingTimesFile.open( ppDataAddress+ "/" + "eqSamplingTimes.txt");
    if (eqSamplingTimesFile.is_open()){ //checking whether the file is open
        string tp;
        getline(eqSamplingTimesFile, tp);
        while(getline(eqSamplingTimesFile, tp)){ //read data from file object and put it into string.

            eq_snapshots_indices.push_back(std::stoi(tp));  

        }
    }
    eqSamplingTimesFile.close();
    int nSamplesEq = eq_snapshots_indices.size();

    if (nSamplesEq < 10)
    {
      cout<<"##########################################################"<<endl;
      cout<<"Error!: Number of samples for equilibrated state are less than 10!"<<endl;
      cout<<"##########################################################"<<endl;
      exit(0);
    }
    ////////// TIME /////////////////////

    //////////////// these are for equilibration sampling /////////////////////////
    vector<vector<int>> nNeighCellEq(NumCells+1, vector<int>(nSamplesEq));
    vector<vector<double>> qCellEq(NumCells+1, vector<double>(nSamplesEq)); // qVor of each cell in the equilibrated samples
    vector<vector<int>> nNeighCellApproxEq(NumCells+1, vector<int>(nSamplesEq)); // This is an approximation based on q value

    int sampleEqC = 0; // counter of the equilibrated samples
    //////////////// these are for equilibration sampling /////////////////////////


    vector<vector<double>> qCellVorAvgStd(2, vector<double>(totalSnapshotsNum));
    double qMeanVal, qStdVal;


    ////////////// INITIAL CONDITION /////////
    snapshotC = 0;
    loadDbl1DVec(xComInit, ppDataAddress+"/xComInit.csv");
    loadDbl1DVec(yComInit, ppDataAddress+"/yComInit.csv");

    for (int i = 0; i < NumCells+1; i++)
    {
      nNeighCell[i] = 0;
      periCell[i] = 0.0;
      areaCell[i] = 0.0;
      qCell[i] = 0.0;
    }

    for (cellC = 0; cellC <= NumCells; cellC++)
    {
      xyCom[cellC][0] = xComInit[cellC];
      xyCom[cellC][1] = yComInit[cellC];
      
      if (cellC>0)
      {
        xyCom_noVoid[cellC-1][0] = fmod(xyCom[cellC][0] + xShift, Lx);
        if (xyCom_noVoid[cellC-1][0]<-1e-10)
        {
          xyCom_noVoid[cellC-1][0] += Lx;
        }else if (xyCom_noVoid[cellC-1][0] < 0.0)
        {
          xyCom_noVoid[cellC-1][0] = 0.0;
        }
        
        
        xyCom_noVoid[cellC-1][1] = fmod(xyCom[cellC][1] + yShift, Ly);
        if (xyCom_noVoid[cellC-1][1]<-1e-10)
        {
          xyCom_noVoid[cellC-1][1] += Ly;
        }else if (xyCom_noVoid[cellC-1][1] < 0.0)
        {
          xyCom_noVoid[cellC-1][1] = 0.0;
        }
      }
    }
    
    vorTes.buildIO(Lx, Ly, xyCom_noVoid);
    vorTes.generateVornoi(false);

    qMeanVal = 0.0;
    for (cellC = 0; cellC <= NumCells; cellC++)
    {
      if (cellC>0)
      {
        nNeighCell[cellC] = vorTes.MeshPoints[cellC-1].neighbors.size();
        periCell[cellC] = 0.0;
        for (int neighC = 0; neighC < nNeighCell[cellC]; neighC++)
        {
          periCell[cellC] += vorTes.MeshPoints[cellC-1].voronoiWalls[neighC];
        }
        // accumulate(vorTes.MeshPoints[cellC-1].voronoiWalls.begin(), vorTes.MeshPoints[cellC-1].voronoiWalls.end(), 0.0);
        areaCell[cellC] = vorTes.MeshPoints[cellC-1].voronoiVolume;
        qCell[cellC] =sqrt(4.0*PI*areaCell[cellC]) / periCell[cellC];
        
        qMeanVal += (qCell[cellC])/NumCells;

      }
      else
      {
        nNeighCell[cellC] = 0;
        periCell[cellC] = NAN;
        areaCell[cellC] = NAN;
        qCell[cellC] = NAN;
      }
    }
    saveInt1DVec(nNeighCell, folderAddress + "/" + "nNeighCellVorInit.csv");
    saveDbl1DVec(periCell, folderAddress+"/"+"periCellVorInit.csv");
    saveDbl1DVec(areaCell, folderAddress+"/"+"areaCellVorInit.csv");
    saveDbl1DVec(qCell, folderAddress+"/"+"qCellVorInit.csv");

    qStdVal = 0.0;
    for (cellC = 1; cellC <= NumCells; cellC++)
    {
      qStdVal += ((qCell[cellC]-qMeanVal)*(qCell[cellC]-qMeanVal)) / (NumCells - 1);
    }
    qStdVal = sqrt(qStdVal);

    qCellVorAvgStd[0][snapshotC] = qMeanVal;
    qCellVorAvgStd[1][snapshotC] = qStdVal;



    // if (runC==3)
    // {
    //   Mesh vorTesTest;
    //   testPlot(vorTesTest, Lx, Ly, xyCom_noVoid);
    // }

    ////////////// INITIAL CONDITION /////////


    ////////////// for t > 0 samples /////////
    snapshotC++;
    bunchC = 0;
    std::ofstream eq_samples_edges_data;
    std::ofstream eq_samples_distances_data;

    eq_samples_edges_data.open(folderAddress + "/" + "eq_edges_data.txt");
    eq_samples_distances_data.open(folderAddress + "/" + "eq_distances_data.txt");
    while (1)
    {

      std::ifstream file(ppDataAddress+"/xComBunch_"+to_string(bunchC)+".csv");
      if (!file.good())
      {
        file.close();
        break;
      }
      file.close();

      doubleMatLoader(xComBunch, ppDataAddress+"/xComBunch_"+to_string(bunchC)+".csv" , NumCells+1 , samplesPerWrite);
      doubleMatLoader(yComBunch, ppDataAddress+"/yComBunch_"+to_string(bunchC)+".csv" , NumCells+1 , samplesPerWrite);

      for (int i = 0; i < NumCells+1; i++)
      {
        nNeighCell[i] = 0;
        periCell[i] = 0.0;
        areaCell[i] = 0.0;
        qCell[i] = 0.0;

        for (int j = 0; j < samplesPerWrite; j++)
        {
          nNeighCellBunch[i][j] = 0;
          periCellBunch[i][j] = 0.0;
          areaCellBunch[i][j] = 0.0;
          qCellBunch[i][j] = 0.0;
        }
      }
      

      ////// loop on bunch //////
      snapshotCDum = 0;
      while (snapshotCDum<samplesPerWrite)
      {
        
        for (cellC = 0; cellC <= NumCells; cellC++)
        {
          xyCom[cellC][0] = xComBunch[cellC][snapshotCDum];
          xyCom[cellC][1] = yComBunch[cellC][snapshotCDum];
          
          if (cellC>0)
          {
            xyCom_noVoid[cellC-1][0] = fmod(xyCom[cellC][0] + xShift, Lx);
            if (xyCom_noVoid[cellC-1][0]<-1e-10)
            {
              xyCom_noVoid[cellC-1][0] += Lx;
            }else if (xyCom_noVoid[cellC-1][0] < 0.0)
            {
              xyCom_noVoid[cellC-1][0] = 0.0;
            }
            
            xyCom_noVoid[cellC-1][1] = fmod(xyCom[cellC][1] + yShift, Ly);
            if (xyCom_noVoid[cellC-1][1]<-1e-10)
            {
              xyCom_noVoid[cellC-1][1] += Ly;
            }else if (xyCom_noVoid[cellC-1][1] < 0.0)
            {
              xyCom_noVoid[cellC-1][1] = 0.0;
            }
            
          }
        }

        
        vorTes.buildIO(Lx, Ly, xyCom_noVoid); // building the voronoi tessellation
        vorTes.generateVornoi(false);
        
        qMeanVal = 0.0;
        for (cellC = 0; cellC <= NumCells; cellC++)
        {
          if (cellC>0)
          {
            nNeighCell[cellC] = vorTes.MeshPoints[cellC-1].neighbors.size();
            periCell[cellC] = 0.0;
            for (int neighC = 0; neighC < nNeighCell[cellC]; neighC++)
            {
              periCell[cellC] += vorTes.MeshPoints[cellC-1].voronoiWalls[neighC];
            }
            areaCell[cellC] = vorTes.MeshPoints[cellC-1].voronoiVolume;
            qCell[cellC] =sqrt(4.0*PI*areaCell[cellC]) / periCell[cellC];

            qMeanVal += (qCell[cellC])/NumCells;
          }
          else
          {
            nNeighCell[cellC] = NAN;
            periCell[cellC] = NAN;
            areaCell[cellC] = NAN;
            qCell[cellC] = NAN;
          }

          nNeighCellBunch[cellC][snapshotCDum] = nNeighCell[cellC];
          periCellBunch[cellC][snapshotCDum] = periCell[cellC];
          areaCellBunch[cellC][snapshotCDum] = areaCell[cellC];
          qCellBunch[cellC][snapshotCDum] = qCell[cellC];
        }

        
    
        // if (snapshotC == totalSnapshotsNum-1 && runC==2)
        // if (snapshotC == 1 && runC==3)
        // {
        //   Mesh vorTesTest;
        //   testPlot(vorTesTest, Lx, Ly, xyCom_noVoid);
        // }

        // cout<<"runC: "<<runC<<", snapshotC: "<<snapshotC<<endl;

        
        qStdVal = 0.0;
        for (cellC = 1; cellC <= NumCells; cellC++)
        {
          qStdVal += ((qCell[cellC]-qMeanVal)*(qCell[cellC]-qMeanVal)) / (NumCells - 1);
        }
        qStdVal = sqrt(qStdVal);

        qCellVorAvgStd[0][snapshotC] = qMeanVal;
        qCellVorAvgStd[1][snapshotC] = qStdVal;


        ////////////////// EQ SAMPLING /////////////////////
        if (snapshotC == eq_snapshots_indices[sampleEqC])
        {
          std::string veticesFileName = "verticesCellsEq_"+to_string(sampleEqC)+".dat";
          vorTes.saveVertices(veticesFileName.c_str(), Lx, Ly);
          std::string cut_command_str = "mv "+ veticesFileName +" "+folderAddress;
          int result = system(cut_command_str.c_str());


          vector<double> edges_data;
          vector<double> distances_data;
          edges_data.clear();
          distances_data.clear();

          vector<vector<int>> edge_counted(NumCells+1, vector<int>(NumCells+1));
          for (int i = 0; i <= NumCells; i++)
          {
            for (int j = 0; j <= NumCells; j++)
            {
              edge_counted[i][j] = 0;
            }
          }
          int neighInd;

          for (cellC = 0; cellC <= NumCells; cellC++)
          {
            if (cellC>0)
            {
              nNeighCellEq[cellC][sampleEqC] = nNeighCell[cellC];
              qCellEq[cellC][sampleEqC] = qCell[cellC];
              int n_neigh_app = 2;//approximate number of the neighbors (based on the value of qVor of the cell)
              do
              {
                n_neigh_app++; 
                nNeighCellApproxEq[cellC][sampleEqC] = n_neigh_app; // This is an approximation based on q value
              } while (qCell[cellC] > 0.5*(polygons_q[n_neigh_app]+polygons_q[n_neigh_app+1]) );

              for (int neighIndC = 0; neighIndC < vorTes.MeshPoints[cellC-1].neighbors.size(); neighIndC++)
              {

                neighInd = vorTes.MeshPoints[cellC-1].neighbors[neighIndC];
                if (edge_counted[cellC][neighInd+1] == 0)
                {
                  edges_data.push_back(vorTes.MeshPoints[cellC-1].voronoiWalls[neighIndC]);
                  distances_data.push_back(vorTes.MeshPoints[cellC-1].neighborDist[neighIndC]);
                  edge_counted[cellC][neighInd+1] = 1;
                  edge_counted[neighInd+1][cellC] = 1;
                }
              }
            }
            else
            {
              nNeighCellEq[cellC][sampleEqC] = NAN;
              qCellEq[cellC][sampleEqC] = NAN;
              nNeighCellApproxEq[cellC][sampleEqC]=NAN;
            }
          }

          int data_length = edges_data.size();
          for (int elementC = 0; elementC < data_length; elementC++)
          {
              eq_samples_edges_data << edges_data[elementC];
              eq_samples_distances_data << distances_data[elementC];

              if (elementC != data_length-1)
              {
                eq_samples_edges_data << "\t";
                eq_samples_distances_data << "\t";
              }
              else
              {
                eq_samples_edges_data << "\n";
                eq_samples_distances_data << "\n";
              }
          }
          sampleEqC++;
        }
        ////////////////// EQ SAMPLING /////////////////////

        snapshotC++;
        snapshotCDum++;
      }
      ////// loop on bunch //////

      // saveIntMatCSV(nNeighCellBunch, folderAddress + "/" + "nNeighCellVorSamples_" + to_string(bunchC) + ".csv");
      // double2DVecSaver(periCellBunch, folderAddress+"/"+"periCellVorSamples_"+to_string(bunchC)+".csv");
      // double2DVecSaver(areaCellBunch, folderAddress+"/"+"areaCellVorSamples_"+to_string(bunchC)+".csv");
      // double2DVecSaver(qCellBunch, folderAddress+"/"+"qCellVorSamples_"+to_string(bunchC)+".csv");

      bunchC++;
    }
    ////////////// for t > 0 samples /////////
    saveIntMatCSV(nNeighCellEq, folderAddress + "/" + "nNeighCellEq.csv");
    double2DVecSaver(qCellEq, folderAddress+"/"+"qCellEq.csv");
    saveIntMatCSV(nNeighCellApproxEq, folderAddress + "/" + "nNeighCellApproxEq.csv");

    eq_samples_edges_data << "\b";
    eq_samples_distances_data << "\b";
    eq_samples_edges_data.close();
    eq_samples_distances_data.close();
    /////// deriving PDF ///////////
    // ofstream nNeighCellEqPDF;
    // ofstream nNeighCellApproxEqPDF;
    // nNeighCellEqPDF.open(folderAddress + "/" + "nNeighCellEqPDF.csv");
    // nNeighCellApproxEqPDF.open(folderAddress + "/" + "nNeighCellApproxEqPDF.csv");
    int polygons_q_size = polygons_q.size();
    vector<vector<int>> nNeighCellEqPDF(nSamplesEq+1, vector<int>(polygons_q_size));
    vector<vector<int>> nNeighCellApproxEqPDF(nSamplesEq+1, vector<int>(polygons_q_size));
    for (int i = 0; i < polygons_q_size; i++)
    {
      nNeighCellEqPDF[0][i] = i;
      nNeighCellApproxEqPDF[0][i] = i;
    }

    vector<int> PDF;
    vector<int> PDF_approx;
    for (sampleEqC = 0; sampleEqC < nSamplesEq; sampleEqC++)
    {
      PDF.clear();
      PDF_approx.clear();

      for (int i = 0; i < polygons_q_size; i++)
      {
        PDF.push_back(0);
        PDF_approx.push_back(0);
      }
      
      for (int cellC = 1; cellC <= NumCells; cellC++)
      {
        PDF[nNeighCellEq[cellC][sampleEqC]]++;
        PDF_approx[nNeighCellApproxEq[cellC][sampleEqC]]++;
      }

      for (int i = 0; i < polygons_q_size; i++)
      {
        nNeighCellEqPDF[sampleEqC+1][i] = PDF[i];
        nNeighCellApproxEqPDF[sampleEqC+1][i] = PDF_approx[i];
      }
    }

    saveIntMatCSV(nNeighCellEqPDF, folderAddress + "/" + "nNeighCellEqPDF.csv");
    saveIntMatCSV(nNeighCellApproxEqPDF, folderAddress + "/" + "nNeighCellApproxEqPDF.csv");

    // nNeighCellEqPDF.close();
    // nNeighCellApproxEqPDF.close();
    /////// deriving PDF ///////////
    


    cout<<"runC: "<<runC<<endl;

    double2DVecSaver(qCellVorAvgStd, folderAddress+"/"+"qCellVorAvgStd.csv");

    runC++;
  }

  /////////////// POLYGONS Q VALS /////////////////


  
  return 0;
}

//////////////////////////////////////////////////////
/////////////////// FUNCTIONS ////////////////////////
//////////////////////////////////////////////////////

bool directoryExists(const std::string& path) {
    return access(path.c_str(), F_OK | R_OK | X_OK) == 0;
}

void dataFolderDeleter()
{
    // const char* subfolderPath = "data";  // Replace with actual path
    
    const char* command = "rm -r data";  // Replace with your desired command

    int result = system(command);
    
    //if (rmdir(subfolderPath) == 0)
    //{
    //    std::cout << "Subfolder successfully deleted.\n";
    //}
    //else
    //{
    //    perror("Error deleting subfolder");
    //}
}

int totalSnapshotsNumFunction(const int samplesPerWrite)
{
    int totalNum=1; // 1 for the initial snapshot
    int bunchC =0;
    // string namestr = "data/tSamples_"+to_string(bunchC)+".csv";
    ifstream tSamplesFile;
    bool cond;
    do
    {
        string namestr = "data/tSamples_"+to_string(bunchC)+".csv";
        // ifstream tSamplesFile(namestr);
        ifstream tSamplesFile(namestr);
        // cout<<bunchC<<endl;
        // cout<<tSamplesFile.is_open()<<endl;
        cond = tSamplesFile.is_open();
        totalNum+=samplesPerWrite;
        bunchC++;
    }while(cond);
    totalNum-=samplesPerWrite;
    bunchC--;


    return totalNum;
}

void simulationDataReader_sq(int* LPtr, double* AlphaPtr, double* KbPtr, double* TemPtr,  int* NumCellsPtr, \
                          double* AvgCellAreaPtr, double* LambdaPtr, long* maxMCStepsPtr, int* samplesPerWritePtr, \
                          int* printingTimeIntervalPtr, string* initConfigPtr)
{   
    // int L_read;
    // int NumCells_read;
    // int samplesPerWrite_read;

    fstream newfile;
    newfile.open("simulationData_vec.csv",ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()){ //checking whether the file is open
        string tp;
        while(getline(newfile, tp)){ //read data from file object and put it into string.
            // cout << tp << "\n"; //print the data of the string
            // if (strstr((const char)tp, "L = "))
            const char* tpChar = tp.c_str();
            if (strstr(tpChar, "L = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *LPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Alpha = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *AlphaPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Kb = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *KbPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Tem = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *TemPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "NumCells = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *NumCellsPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "AvgCellArea = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *AvgCellAreaPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Lambda = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *LambdaPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "SweepLength = "))
            {
                // Do nothing
                continue;
            }
            if (strstr(tpChar, "maxMCSteps = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *maxMCStepsPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "samplesPerWrite = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *samplesPerWritePtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "printingTimeInterval = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *printingTimeIntervalPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "initConfig = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *initConfigPtr = num_str;
                continue;
            }
        }
    }
    newfile.close(); //close the file object.

    // *L_read_ptr = L_read;
    // *NumCells_read_ptr = NumCells_read;
    // *samplesPerWrite_read_ptr = samplesPerWrite_read;
}

void simulationDataReader_irr(int* NSitesPtr, double* LxPtr, double* LyPtr, double* AlphaPtr, double* KbPtr, double* TemPtr,  int* NumCellsPtr, \
                          double* AvgCellAreaPtr, double* LambdaPtr, long* maxMCStepsPtr, int* samplesPerWritePtr, \
                          int* printingTimeIntervalPtr, int* numLinkedListPtr, string* initConfigPtr)
{   
    // int L_read;
    // int NumCells_read;
    // int samplesPerWrite_read;

    fstream newfile;
    newfile.open("simulationData_vec.csv",ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()){ //checking whether the file is open
        string tp;
        while(getline(newfile, tp)){ //read data from file object and put it into string.
            // cout << tp << "\n"; //print the data of the string
            // if (strstr((const char)tp, "NSites = "))
            const char* tpChar = tp.c_str();
            if (strstr(tpChar, "NSites = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *NSitesPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Lx = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *LxPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Ly = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *LyPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Alpha = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *AlphaPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Kb = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *KbPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Tem = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *TemPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "NumCells = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *NumCellsPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "AvgCellArea = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *AvgCellAreaPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Lambda = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *LambdaPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "SweepLength = "))
            {
                // Do nothing
                continue;
            }
            if (strstr(tpChar, "maxMCSteps = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *maxMCStepsPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "samplesPerWrite = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *samplesPerWritePtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "printingTimeInterval = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *printingTimeIntervalPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "numLinkedList = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *numLinkedListPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "initConfig = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *initConfigPtr = num_str;
                continue;
            }
        }
    }
    newfile.close(); //close the file object.

    // *L_read_ptr = L_read;
    // *NumCells_read_ptr = NumCells_read;
    // *samplesPerWrite_read_ptr = samplesPerWrite_read;
}

std::vector<std::string> parseString(const std::string& inputStr) {
    std::vector<std::string> tokens;
    std::stringstream ss(inputStr);
    std::string token;
    while (std::getline(ss, token, ',')) {
        tokens.push_back(token);
    }
    return tokens;
}

void timeMaker(long* time, long length, int samplesPerWrite)
{
    long* timePtr = time;
    *timePtr = 0;

    std::ofstream outfile("pp_data/time.csv");
    outfile << 0;
    outfile << ",";

    long ind=1;
    int bunchC=0;
    while(1)
    {
        string namestr = "data/tSamples_"+to_string(bunchC)+".csv";
        ifstream tSamplesFile(namestr);
        std::string line;
        std::getline(tSamplesFile, line);
        std::vector<std::string> lineTokens = parseString(line);

        for(int i=0; i<samplesPerWrite; i++)
        {
            timePtr = time+ind;
            *timePtr = (long)(stoi(lineTokens[i]));
            outfile << (long)(stoi(lineTokens[i]));
            if (ind<length-1)
            {outfile << ",";}
            ind++;
        }
        bunchC++;

        if(ind>=length)
        {
            ind--;
            outfile.close();
            break;
        }
    }   
}

void energyMaker(double* energy, long length, int samplesPerWrite)
{
    double* energyPtr = energy;
    
    std::ofstream outfile("pp_data/energy.csv");

    std::ifstream e0File("init/e0Initial.csv");
    e0File >> *energyPtr;
    e0File.close();

    outfile << *energyPtr;
    outfile << ",";

    long ind=1;
    int bunchC=0;
    while(1)
    {
        string namestr = "data/eSamples_"+to_string(bunchC)+".csv";
        ifstream eSamplesFile(namestr);
        std::string line;
        std::getline(eSamplesFile, line);
        std::vector<std::string> lineTokens = parseString(line);

        for(int i=0; i<samplesPerWrite; i++)
        {
            energyPtr = energy+ind;
            *energyPtr = (double)(stod(lineTokens[i]));
            outfile << (double)(stod(lineTokens[i]));
            if (ind<length-1)
            {outfile << ",";}
            ind++;
        }
        bunchC++;

        if(ind>=length)
        {
            ind--;
            outfile.close();
            break;
        }
    }   
}

void sigmaMatLoader(vector<vector<int>>& mat, const std::string &filename, int rows, int cols)
{
  std::ifstream input_file(filename);
  if (!input_file)
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }

  std::vector<std::vector<int>> data;
  std::string line;
  while (getline(input_file, line))
  {
    std::vector<int> row;
    std::stringstream ss(line);
    int value;
    while (ss >> value)
    {
      row.push_back(value);
      if (ss.peek() == ',')
      {
        ss.ignore();
      }
    }
    data.push_back(row);
  }

  input_file.close();

  // rows = data.size();
  // cols = data[0].size();

  // int** array = new int*[rows];
  for (int i = 0; i < rows; ++i)
  {
    // array[i] = new int[cols];
    for (int j = 0; j < cols; ++j)
    {
      mat[i][j] = data[i][j];
    }
  }

  // return array;
}

void doubleMatLoader(vector<vector<double>>& mat, const std::string &filename, int rows, int cols)
{
  std::ifstream input_file(filename);
  if (!input_file)
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }

  std::vector<std::vector<double>> data;
  std::string line;
  while (getline(input_file, line))
  {
    std::vector<double> row;
    size_t found = line.find("nan");
    if (found != std::string::npos) 
    {
      for (int j = 0; j < cols; j++)
      {
        row.push_back(NAN);
      }  
    }
    else{
      std::istringstream ss(line);
      double value;
      while (ss >> value)
      {
        row.push_back(value);
        if (ss.peek() == ',')
        {
          ss.ignore();
        }
      }
    }
    data.push_back(row);
  }

  input_file.close();

  // rows = data.size();
  // cols = data[0].size();

  // int** array = new int*[rows];
  for (int i = 0; i < rows; ++i)
  {
    // array[i] = new int[cols];
    for (int j = 0; j < cols; ++j)
    {
      mat[i][j] = data[i][j];
    }
  }

  // return array;
}

void configurationalCalcPerSample(vector<int>& areaSample, vector<int>& periSample,\
vector<double>& isoperiSample, vector<double>& xComSample, vector<double>& yComSample, \
const vector<vector<int>>& sigmaMat, const vector<vector<int>>& sitesX, const vector<vector<int>>& sitesY)
{
    int L = sigmaMat.size();
    int NumCells = areaSample.size() -1;

    int row, col, cellIndex;
    int row_r, col_r, row_d, col_d; // for right and down neighbors
    int cellIndex_r, cellIndex_d;

    for (cellIndex=0; cellIndex<=NumCells; cellIndex++) //must start from 0, because we must make everything zero for index 0.
    {
        areaSample[cellIndex] = 0;
        periSample[cellIndex] = 0;
        isoperiSample[cellIndex] = 0;
        xComSample[cellIndex] = 0;
        yComSample[cellIndex] = 0;
    }
    
    for (row=0; row<L; row++)
    {
        for (col=0; col<L; col++)
        {
            cellIndex = sigmaMat[row][col];
            areaSample[cellIndex]++;

            //right edge
            row_r = row;
            col_r = (col+1)%L;
            cellIndex_r = sigmaMat[row_r][col_r];
            if(cellIndex != cellIndex_r)
            {
                periSample[cellIndex]++;
                periSample[cellIndex_r]++;
            }

            //down edge
            row_d = (row+1)%L;
            col_d = col;
            cellIndex_d = sigmaMat[row_d][col_d];
            if(cellIndex != cellIndex_d)
            {
                periSample[cellIndex]++;
                periSample[cellIndex_d]++;
            }

            xComSample[cellIndex] += (double)(sitesX[row][col]);
            yComSample[cellIndex] += (double)(sitesY[row][col]);
        }
    }

    for (cellIndex=1; cellIndex<=NumCells; cellIndex++) //must start from 1, because for 0, everything is already 0
    {
        isoperiSample[cellIndex] = 2.*sqrt(PI)*sqrt((double)(areaSample[cellIndex])) / ((double)(periSample[cellIndex]));
        xComSample[cellIndex] = xComSample[cellIndex] / ((double)(areaSample[cellIndex]));
        yComSample[cellIndex] = yComSample[cellIndex] / ((double)(areaSample[cellIndex]));
    }
}

void configurationalCalcPerSample_irr(vector<double>& areaSample, vector<double>& periSample,\
vector<double>& isoperiSample, vector<double>& ARSample, vector<double>& circSample, vector<double>& xComSample, vector<double>& yComSample, \
const std::vector<int>& sigmaMat, const std::vector<double>& sitesX, const std::vector<double>& sitesY,\
const std::vector<int>& neighborNum, const std::vector<vector<int>>& neighborsList, const std::vector<vector<double>>& edges,\
const std::vector<double>& latticeArea)
{
    int NSites = sigmaMat.size();
    int NumCells = areaSample.size() -1;

    
    int siteC, cellIndex, neighInd;
    // int row_r, col_r, row_d, col_d; // for right and down neighbors
    // int cellIndex_r, cellIndex_d;

    for (cellIndex=0; cellIndex<=NumCells; cellIndex++) //must start from 0, because we must make everything zero for index 0.
    {
        areaSample[cellIndex] = 0;
        periSample[cellIndex] = 0;
        isoperiSample[cellIndex] = 0;
        ARSample[cellIndex] = 0;
        circSample[cellIndex] = 0;
        xComSample[cellIndex] = 0;
        yComSample[cellIndex] = 0;
    }
    
    double filledArea =0.;
    
    for (siteC = 0; siteC < NSites; siteC++)
    {
        cellIndex = sigmaMat[siteC];
        areaSample[cellIndex]+=latticeArea[siteC];
        filledArea += (latticeArea[siteC])*(cellIndex>0);

        for (int neighIndC = 0; neighIndC < neighborNum[siteC]; neighIndC++)
        {
            neighInd = neighborsList[siteC][neighIndC];
            if (cellIndex != sigmaMat[neighInd])
            {
                periSample[cellIndex] += 0.5 * edges[siteC][neighIndC];
                periSample[sigmaMat[neighInd]] +=  0.5 * edges[siteC][neighIndC];
            }
        }

        xComSample[cellIndex] += (double)(sitesX[siteC]*latticeArea[siteC]);
        yComSample[cellIndex] += (double)(sitesY[siteC]*latticeArea[siteC]);
    }
    
    // totXCom = 0.;
    // totYCom = 0.;
    for (cellIndex=0; cellIndex<=NumCells; cellIndex++) //must start from 1, because for 0, everything is already 0
    {
        isoperiSample[cellIndex] = 2.*sqrt(PI)*sqrt((double)(areaSample[cellIndex])) / ((double)(periSample[cellIndex]));
        xComSample[cellIndex] = xComSample[cellIndex] / ((double)(areaSample[cellIndex]));
        yComSample[cellIndex] = yComSample[cellIndex] / ((double)(areaSample[cellIndex]));

        // if (cellIndex>0)
        // {
        //   totXCom += (xComSample[cellIndex] * areaSample[cellIndex] / filledArea);
        //   totYCom += (yComSample[cellIndex] * areaSample[cellIndex] / filledArea);
        // }
    }

    // Aspect Ratio calculation
    vector<double> IxxSample(NumCells + 1);
    vector<double> IyySample(NumCells + 1);
    vector<double> IxySample(NumCells + 1);
    double I_max, I_min, I_1, I_2, Ixx, Ixy, Iyy;
    for (cellIndex=0; cellIndex<=NumCells; cellIndex++)
    {
        IxxSample[cellIndex] = 0.0;
        IyySample[cellIndex] = 0.0;
        IxySample[cellIndex] = 0.0;
    }
    for (siteC = 0; siteC < NSites; siteC++)
    {
        cellIndex = sigmaMat[siteC];
        // IxxSample[cellIndex] += pow((sitesX[siteC]-xComSample[cellIndex]),2) * latticeArea[siteC] / areaSample[cellIndex];
        // IyySample[cellIndex] += pow((sitesY[siteC]-yComSample[cellIndex]),2) * latticeArea[siteC] / areaSample[cellIndex];
        // IxySample[cellIndex] += (sitesX[siteC]-xComSample[cellIndex])*(sitesY[siteC]-yComSample[cellIndex]) * latticeArea[siteC] / areaSample[cellIndex];

        IxxSample[cellIndex] += pow((sitesX[siteC]-xComSample[cellIndex]),2) * latticeArea[siteC];
        IyySample[cellIndex] += pow((sitesY[siteC]-yComSample[cellIndex]),2) * latticeArea[siteC];
        IxySample[cellIndex] += (sitesX[siteC]-xComSample[cellIndex])*(sitesY[siteC]-yComSample[cellIndex]) * latticeArea[siteC];
    }

    for (cellIndex=0; cellIndex<=NumCells; cellIndex++)
    {
        Ixx = IxxSample[cellIndex];
        Iyy = IyySample[cellIndex];
        Ixy = IxySample[cellIndex];
        I_min = 0.5*((Ixx+Iyy)-sqrt((Ixx-Iyy)*(Ixx-Iyy)+4.0*Ixy*Ixy));
        I_max = 0.5*((Ixx+Iyy)+sqrt((Ixx-Iyy)*(Ixx-Iyy)+4.0*Ixy*Ixy));
        ARSample[cellIndex] = sqrt(I_max/I_min);
        circSample[cellIndex] = (areaSample[cellIndex]*areaSample[cellIndex])/(2*PI*(Ixx+Iyy));
    }

    // for (cellIndex=1; cellIndex<=NumCells; cellIndex++)
    // {
    //     xComSample[cellIndex] -= (totXCom-97.117);
    //     yComSample[cellIndex] -= (totYCom-96.8557);
    // }
    
}

void sigmaMatExtractor(vector<vector<vector<int>>>& mat, std::ifstream& samplesFile)
{
    
    int L = mat.size();
    int samplesPerWrite = mat[0][0].size();

    // cout<< "entered into function"<<endl<<endl;
    int sampleC = 0;
    int row, col;
    string line;
    string teststr="---";

    for (int i =0; i<(samplesPerWrite*L+samplesPerWrite); i++) {
        getline(samplesFile, line);
        if(strstr(line.c_str(), "---"))
        {   
            sampleC++;
            continue;
        }
        else
        {
            // std::cout << tp << std::endl;
            std::vector<std::string> lineTokens = parseString(line);
            row = i%(1+L);
            for (col=0; col<L; col++)
            {
                mat[row][col][sampleC] = (int)(stoi(lineTokens[col]));
            }
        }
    }
}

void printForTest(const vector<vector<int>>& mat)
{
    int L = mat.size();

    for(int row =0; row<L; row++)
    {
        for(int col =0; col<L; col++)
        {
            cout<<mat[row][col];
            if(col==L-1)
            {cout<<"\n";}
            else
            {cout<<", ";}
        }
    }
}

void doubleMatSaver(double* matPtr, int Rows, int Cols, const std::string &filename)
{

  std::ofstream outfile(filename);
  double* ptr;

  for (int i = 0; i < Rows; i++)
  {
    for (int j = 0; j < Cols; j++)
    {
        ptr = matPtr+(i*Cols+j);

        outfile << *ptr;

        if (j != Cols - 1)
        {
        outfile << ",";
        }
    }
    outfile << std::endl;
  }
  outfile.close();
}

void intMatSaver(int* matPtr, int Rows, int Cols, const std::string &filename)
{

  std::ofstream outfile(filename);
  int* ptr;

  for (int i = 0; i < Rows; i++)
  {
    for (int j = 0; j < Cols; j++)
    {
        ptr = matPtr+(i*Cols+j);

        outfile << *ptr;

        if (j != Cols - 1)
        {
        outfile << ",";
        }
    }
    outfile << std::endl;
  }
  outfile.close();
}

void double2DVecSaver(const vector<vector<double>>& mat, const std::string &filename)
{
    int Rows = mat.size();
    int Cols = mat[0].size();

  std::ofstream outfile(filename);

  for (int i = 0; i < Rows; i++)
  {
    for (int j = 0; j < Cols; j++)
    {
        outfile << mat[i][j];

        if (j != Cols - 1)
        {
        outfile << ",";
        }
    }
    outfile << std::endl;
  }
  outfile.close();
}

void double1DVecSaver(const vector<double>& mat, const std::string &filename)
{
    int Rows = mat.size();

  std::ofstream outfile(filename);

  for (int i = 0; i < Rows; i++)
  {
    outfile << mat[i];
    if (i != Rows - 1)
    {
    outfile << ",";
    }
  }
  outfile.close();
}

void saveInt1DVec(const vector<int>& quantity, const std::string &filename)
{
  
  int length = quantity.size();

  std::ofstream outfile(filename);

  for (int i = 0; i < length; i++)
  {
    outfile << quantity[i];
    if (i < length - 1)
    {
      outfile << "\n";
    }
  }
  outfile.close();
}

void saveDbl1DVec(const vector<double>& quantity, const std::string &filename)
{
  
  int length = quantity.size();

  std::ofstream outfile(filename);

  for (int i = 0; i < length; i++)
  {
    outfile << std::setprecision(10)<< quantity[i];
    if (i < length - 1)
    {
      outfile << "\n";
    }
  }
  outfile.close();
}

void loadInt1DVec(vector<int>& mat, const std::string &filename)
{
  int length = mat.size();

  std::ifstream input_file(filename);
  if (!input_file)
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }

  std::vector<int> data;
  std::string line;
  while (getline(input_file, line))
  {
    // std::vector<int> row;
    std::stringstream ss(line);
    int value;
    ss >> value;
    data.push_back(value);
  }
  input_file.close();

  for (int i = 0; i < length; i++)
  {
    mat[i] = data[i];
  }

  // return array;
}

void loadInt1DVecHor(vector<int>& mat, const std::string &filename)
{
  int length = mat.size();

  std::ifstream input_file(filename);
  if (!input_file)
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }

  std::vector<int> data;
  std::string line;
  while(getline(input_file, line))
  {
    // std::vector<int> row;
    std::stringstream ss(line);
    // cout<<line<<endl;
    // cout<<ss<<endl;

    int value;
    // ss >> value;
    // data.push_back(value);
    // ss >> value;
    // ss >> value;
    // ss >> value;
    // ss >> value;
    while (ss >> value)
    {
      data.push_back(value);
      if (ss.peek() == ',')
      {
        ss.ignore();
      }
    }
  }

  input_file.close();

  for (int i = 0; i < length; i++)
  {
    mat[i] = data[i];
  }

  // return array;
}

void loadDbl1DVec(vector<double>& mat, const std::string &filename)
{
  int length = mat.size();

  std::ifstream input_file(filename);
  if (!input_file)
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }

  std::vector<double> data;
  std::string line;
  while (getline(input_file, line))
  {
    // std::vector<int> row;
    std::stringstream ss(line);
    double value;
    ss >> value;
    data.push_back(value);
  }
  input_file.close();

  for (int i = 0; i < length; i++)
  {
    mat[i] = data[i];
  }

  // return array;
}

void latticeCreatorPhase1(int* maxNeighPtr, vector<int>& neighborNum, vector<double>& latticeX,  \
                          vector<double>& latticeY,  vector<double>& latticeArea )
{
    int maxNeigh = 0;
    int neighC = 0;
    int siteC = 0;
    ifstream newfile;
    newfile.open("lattice/neighs.dat",ios::in);
    string tp;
    while(getline(newfile, tp))
    {
      if(siteC != stoi(tp))
      {
        cout<<"Indices inconsistency in reading lattice! (neighs)"<<endl;
        cout<<"Program ended!"<<endl;
        exit(1);
      }

      for (int i = 0; i < 3; i++)
      {
        std::size_t pos = tp.find('\t');  // find position of '='
        tp = tp.substr(pos+1);  // extract substring starting from position after '='

        if (i==0)
        {
          latticeX[siteC] = stod(tp);
        }
        if (i==1)
        {
          latticeY[siteC] = stod(tp);
        }
      }
      neighC=0;

      // const char* tpChar = tp.c_str();
      while(strcmp(tp.c_str(), ""))
      {
        std::size_t pos = tp.find('\t');
        tp = tp.substr(pos+1);  // extract substring starting from position after '='
        neighC++;
        // const char* tpChar = tp.c_str();
      }

      neighborNum[siteC] = neighC;

      if (neighC>maxNeigh)
      {
        maxNeigh = neighC;
      }
      siteC++;
    }
    newfile.close();

    *maxNeighPtr = maxNeigh;

    siteC = 0;
    newfile.open("lattice/vols.dat",ios::in);
    while(getline(newfile, tp))
    {
      if(siteC != stoi(tp))
      {
        cout<<"Indices inconsistency in reading lattice! (vols)"<<endl;
        cout<<"Program ended!"<<endl;
        exit(1);
      }

      std::size_t pos = tp.find('\t');
      tp = tp.substr(pos+1);
      latticeArea[siteC] = stod(tp);
      siteC++;
    }
}

void latticeCreatorPhase2(const int NSites, const double Lx, const double Ly, const int maxNeighbors,\
                          vector<vector<int>>& neighborsList, vector<vector<double>>& edges,\
                          vector<vector<double>>& vorDist)
{
  
  

  for (int i = 0; i < NSites; i++)
  {
    for (int j = 0; j < maxNeighbors; j++)
    {
      neighborsList[i][j] = -1;
      edges[i][j] = -1;
      vorDist[i][j] = -1;
    }
    
  }

  
  ///////////////////NEIGH, DX, DY//////////////////////
  ifstream neighFile;
  neighFile.open("lattice/neighs.dat",ios::in);
  string tp;
  int neighC;
  int siteC = 0;
  while(getline(neighFile, tp))
  {
    if(siteC != stoi(tp))
    {
      cout<<"Indices inconsistency in reading lattice! (neighs, phase 2)"<<endl;
      cout<<"Program ended!"<<endl;
      exit(1);
    }

    for (int i = 0; i < 3; i++)
    {
      std::size_t pos = tp.find('\t');  // find position of '='
      tp = tp.substr(pos+1);  // extract substring starting from position after '='
    }
    
    // const char* tpChar = tp.c_str();
    neighC=0;
    while(strcmp(tp.c_str(), ""))
    {
      neighborsList[siteC][neighC] = stoi(tp);

      std::size_t pos = tp.find('\t');
      tp = tp.substr(pos+1);  // extract substring starting from position after '='
      neighC++;
      // const char* tpChar = tp.c_str();
    }
    siteC++;
  }
  neighFile.close();
  ///////////////////NEIGH, DX, DY//////////////////////


  ///////////////////EDGES//////////////////////
  ifstream wallsFile;
  wallsFile.open("lattice/walls.dat",ios::in);
  siteC = 0;
  while(getline(wallsFile, tp))
  {
    if(siteC != stoi(tp))
    {
      cout<<"Indices inconsistency in reading lattice! (walls, phase 2)"<<endl;
      cout<<"Program ended!"<<endl;
      exit(1);
    }

    for (int i = 0; i < 3; i++)
    {
      std::size_t pos = tp.find('\t');  // find position of '='
      tp = tp.substr(pos+1);  // extract substring starting from position after '='
    }

    // const char* tpChar = tp.c_str();
    neighC=0;
    while(strcmp(tp.c_str(), ""))
    {
      edges[siteC][neighC] = stod(tp);
      std::size_t pos = tp.find('\t');
      tp = tp.substr(pos+1);  // extract substring starting from position after '='
      neighC++;
      // const char* tpChar = tp.c_str();
    }
    siteC++;
  }
  wallsFile.close();
  ///////////////////EDGES//////////////////////

  ///////////////////DISTS//////////////////////
  ifstream distFile;
  distFile.open("lattice/dists.dat",ios::in);
  siteC = 0;
  while(getline(distFile, tp))
  {
    if(siteC != stoi(tp))
    {
      cout<<"Indices inconsistency in reading lattice! (dists, phase 2)"<<endl;
      cout<<"Program ended!"<<endl;
      exit(1);
    }

    for (int i = 0; i < 3; i++)
    {
      std::size_t pos = tp.find('\t');  // find position of '='
      tp = tp.substr(pos+1);  // extract substring starting from position after '='
    }

    // const char* tpChar = tp.c_str();
    neighC=0;
    while(strcmp(tp.c_str(), ""))
    {
      vorDist[siteC][neighC] = stod(tp);
      std::size_t pos = tp.find('\t');
      tp = tp.substr(pos+1);  // extract substring starting from position after '='
      neighC++;
      // const char* tpChar = tp.c_str();
    }
    siteC++;
  }
  distFile.close();
  ///////////////////DISTS//////////////////////

}

void polygonCom(double* xComPtr, double* yComPtr, const vector<vector<double>> polygonVec)
{
  int nNeigh = polygonVec.size();

  vector<vector<double>> orderedPolygonVec;
  vector<int> includedVertices(nNeigh);
  for (int i = 0; i < nNeigh; i++)
  {
    includedVertices[i] = 0;
  }
  int condition;

  double x1, y1, x2, y2, x3, y3;

  int headInd = 0;
  double headX = polygonVec[headInd][0];
  double headY = polygonVec[headInd][1];
  orderedPolygonVec.push_back({headX, headY});
  includedVertices[headInd] = 1;
  int passedVertices = 1;
  while(passedVertices < nNeigh)
  {

    for (int i = 0; i < nNeigh; i++) // i shows the candidate of new head
    {
      if (includedVertices[i]>0)
      {
        continue;
      }
      else
      {
        
        int verticalSwitch =0;
        if (polygonVec[headInd][0] == polygonVec[i][0]) // if x_head == x_head_candidate
        {
          verticalSwitch =1;
        }

        condition=0;
        vector<double> distToLine;

        x1 = headX;
        y1 = headY;
        x2 = polygonVec[i][0];
        y2 = polygonVec[i][1];

        // y = mx+b
        double m;
        double b;
        if (verticalSwitch==0)
        {
          m = (y2-y1)/(x2-x1);
          b = (y1*x2-y2*x1)/(x2-x1);
        }
        
        double dist;
        for (int j = 0; j < nNeigh; j++)
        {
          if (j == headInd || j == i)
          {
            continue;
          }
          else
          {
            x3 = polygonVec[j][0];
            y3 = polygonVec[j][1];
            if (verticalSwitch)
            {
              dist = x3 - x1;
            }
            else
            {
              dist = y3 - (m*x3 + b);
            }
            distToLine.push_back(dist);
          }
        }
        auto maxDistPtr = std::max_element(distToLine.begin(), distToLine.end());
        auto minDistPtr = std::min_element(distToLine.begin(), distToLine.end());
        double maxDist = *maxDistPtr;
        double minDist = *minDistPtr;


        if(maxDist*minDist > 0)
        {
          condition = 1;
        }
        else if (maxDist*minDist < 0)
        {
          condition = 0;
        }
        else
        {
          cout<<"***************************"<<endl;
          cout<<"Problem: maxDist*minDist == 0"<<endl;
          cout<<"***************************"<<endl;
          exit(5);
        }

        if (condition)
        {
          headInd = i;
          headX = polygonVec[headInd][0];
          headY = polygonVec[headInd][1];
          orderedPolygonVec.push_back({headX, headY});
          includedVertices[headInd] = 1;
          passedVertices++;
        }
      }
    }
  }

  orderedPolygonVec.push_back({polygonVec[0][0], polygonVec[0][1]});
  
  double A = 0.;
  double term;
  for (int i = 0; i < nNeigh; i++)
  {
    term = 0.5*(orderedPolygonVec[i][0]*orderedPolygonVec[i+1][1]-orderedPolygonVec[i+1][0]*orderedPolygonVec[i][1]);
    A = A + term;
  }
  
  double xCom =0.;
  double yCom =0.;
  double termX = 0;
  double termY = 0;
  for (int i = 0; i < nNeigh; i++)
  {
    termX = (1/(6*A))*(orderedPolygonVec[i][0]+orderedPolygonVec[i+1][0])*(orderedPolygonVec[i][0]*orderedPolygonVec[i+1][1]-orderedPolygonVec[i+1][0]*orderedPolygonVec[i][1]);
    termY = (1/(6*A))*(orderedPolygonVec[i][1]+orderedPolygonVec[i+1][1])*(orderedPolygonVec[i][0]*orderedPolygonVec[i+1][1]-orderedPolygonVec[i+1][0]*orderedPolygonVec[i][1]);
    xCom = xCom + termX;
    yCom = yCom + termY;
  }
  
  *xComPtr = xCom;
  *yComPtr = yCom;
  
}

void latticeCreatorPhase3(const int NSites, const double Lx, const double Ly, const int maxNeighbors, \
                          const vector<int> neighborNum, const vector<vector<int>> neighborsList,\
                          vector<double>& siteComX, vector<double>& siteComY,\
                          vector<vector<double>>& deltaComX, vector<vector<double>>& deltaComY, vector<vector<double>>& comDist)
{
  
  

  for (int i = 0; i < NSites; i++)
  {
    siteComX[i] = 0;
    siteComY[i] = 0;

    for (int j = 0; j < maxNeighbors; j++)
    {
      comDist[i][j] = -1;
      deltaComX[i][j] = -1;
      deltaComY[i][j] = -1;
    }
  }

  double x1, y1, x2, y2, dx, dy, x_v, y_v, xCom, yCom;

  ///////////////////// siteComX, siteComY ////////////////////
  ifstream verticesFile;
  verticesFile.open("lattice/vertices.dat",ios::in);
  string tp;
  // int neighC;
  int siteC = 0;
  vector<vector<double>> polygonVec;

  while(getline(verticesFile, tp))
  {

    if(tp.find(':') != std::string::npos) // if ":" exists, go to the next line in the file
    {
      continue;
    }
    else if(tp.find('*') != std::string::npos) // if "*" exists, increase siteC by one, and go to the next line in the file
    {
      polygonCom(&xCom, &yCom, polygonVec);
      if (xCom > Lx)
      {
        xCom = xCom - Lx;
      }
      else if (xCom < 0)
      {
        xCom = xCom + Lx;
      }

      if (yCom > Ly)
      {
        yCom = yCom - Ly;
      }
      else if (yCom < 0)
      {
        yCom = yCom + Ly;
      }
      
      siteComX[siteC] = xCom;
      siteComY[siteC] = yCom;

      polygonVec.clear();
      siteC++;
      continue;
    }

    else
    {
      x_v = stod(tp);
      std::size_t pos = tp.find('\t');
      tp = tp.substr(pos+1);  // extract substring starting from position after '='
      y_v = stod(tp);
      polygonVec.push_back({x_v, y_v});
    }
  }
  verticesFile.close();
  ///////////////////// siteComX, siteComY ////////////////////

  ///////////////////// deltaComX, deltaComY, comDist ////////////////////
  for ( siteC = 0; siteC < NSites; siteC++)
  {
    x1 = siteComX[siteC];
    y1 = siteComY[siteC];

    for (int neighC = 0; neighC < neighborNum[siteC]; neighC++)
    {
      x2 = siteComX[neighborsList[siteC][neighC]];
      y2 = siteComY[neighborsList[siteC][neighC]];

      dx = x2-x1;
      dy = y2-y1;
      if(dx > Lx/2.)
      {
        dx = (x2-Lx)-x1;
      }
      if(dx < -Lx/2.)
      {
        dx = (x2+Lx)-x1;
      }
      if(dy > Ly/2.)
      {
        dy = (y2-Ly)-y1;
      }
      if(dy < -Ly/2.)
      {
        dy = (y2+Ly)-y1;
      }

      deltaComX[siteC][neighC] = dx;
      deltaComY[siteC][neighC] = dy;
      comDist[siteC][neighC] = sqrt(pow(dx,2)+pow(dy,2));
    }
  }
//   // saving siteComXY
//   vector<vector<double>> siteComXY(siteComX.size(), vector<double>(2));
//   for (int i = 0; i < NSites; i++)
//   {
//     siteComXY[i][0] = siteComX[i];
//     siteComXY[i][1] = siteComY[i];
//   }
//   saveDoubleMatCSV(siteComXY, "lattice/siteComXY.csv");
  ///////////////////// deltaComX, deltaComY, comDist ////////////////////
}

inline bool existFunc (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

void totComFunc(double& xComTot, double& yComTot, const vector<double>& xComSample, const vector<double>& yComSample)
{
    int NumCells = xComSample.size()-1;
    int cellC;

    xComTot =0.0;
    yComTot =0.0;
    
    for (cellC = 1; cellC <= NumCells; cellC++)
    {
        xComTot += xComSample[cellC]/NumCells;
        yComTot += yComSample[cellC]/NumCells;
    }
}

void saveIntMatCSV(const vector<vector<int>>& sigmaMat, const std::string &filename)
{
  
  int Nx = sigmaMat.size();
  int Ny = sigmaMat[0].size();

  std::ofstream outfile(filename);

  for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {
      outfile << sigmaMat[i][j];

      if (j != Ny - 1)
      {
        outfile << ",";
      }
    }
    outfile << std::endl;
  }
  outfile.close();
}

void testPlot(Mesh &vorTesTest, const double Lx, const double Ly, const vector<vector<double>>& xyCom_noVoid)
{
  vorTesTest.buildIO(Lx, Ly, xyCom_noVoid);
        
  FILE* f = fopen("test.dat", "w");
  for (int i = 0; i < vorTesTest.MeshPoints.size(); i++){
      fprintf(f, "%lf\t%lf\n", vorTesTest.MeshPoints[i].xCoord, vorTesTest.MeshPoints[i].yCoord);
  }
  fclose(f);

  int NumCells = xyCom_noVoid.size();

  vorTesTest.buildRead(Lx, Ly, NumCells, "test.dat");

  // Test.buildRandom(10.0, 10.0, 100);
  vorTesTest.save("lattice.dat");
  vorTesTest.generateVornoi(true);
  printf("neighbor mismatches: %u\n", vorTesTest.forceNeighborSymmetry());

  vorTesTest.saveNeighbors("neighs.dat");

  //By Hossein
  // Test.saveV
  vorTesTest.saveVolumes("vols.dat");

  vorTesTest.saveWalls("walls.dat");

  vorTesTest.saveDists("dists.dat");

  vorTesTest.saveVertices("vertices.dat", Lx, Ly);
  system("python3 voroPlotPy_data.py");
  system("python3 voroPlotPy_0.py");
}
/////////////////////////////////////////////////////////////////
/////////////////////////FUNCTIONS///////////////////////////////
/////////////////////////////////////////////////////////////////