#include <iostream>
#include <iomanip>
#include <math.h>
// #include <array>
#include <vector>
// #include <cstdlib>
// #include <ctime>
// #include <utility>
// #include <tuple>
#include <cmath>
// #include <map>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
// #include <cstring>
#include <string>
#include <string.h>
#include <random>
#include <chrono>
#include <sys/stat.h>
//#include <filesystem>
#include <unistd.h>
#include <algorithm>

using namespace std;
using namespace std::chrono;

///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////#DEFINITIONS///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////-
// #define L 200                                         // lattice size
// //#define Alpha (1.0)                                 // interaction strength
// #define Kb (1.)                                      // Boltzmann constant
// #define Tem (1.)                                     // Temperature
// #define NumCells 1000                                  // Number of cells
// #define AvgCellArea ((double)((L * L) / (NumCells))) // Average area of cells
// #define Lambda (1.)                                  // Area elasticity strength
// #define SweepLength (L * L)                          // Number of attempts per Monte Carlo sweep
// // #define MaxMCSteps 100000                              // Maximum number of MC steps to simulate
// // #define samplingInterval 100                          // Sweeps between smplings
// #define samplesPerWrite 10
// // #define writingInterval (samplesPerWrite * samplingInterval) // Sweeps between writing operations
// #define printingTimeInterval 1000                              // Printing time interval
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////GLOBAL MT RANDOM GENERATOR///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
/*
random_device rd;
int mt_rand_seed = rd();
mt19937 mt_rand(mt_rand_seed);
long MT_MAX = mt_rand.max();
long MT_MIN = mt_rand.min();
*/
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////STRUCTS DEFINITION////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
struct Energy
{
  double inter; // border interaction energy
  double area;  // area elasticity energy
  double act;   // activity energy
  double total; // total energy
};

struct LinkedListElement
{
  // int sitesPair[2];
  // std::pair<int, int> sitesPair;
  int site;
  struct LinkedListElement *next;
};

struct LinkedListElementPair
{
  // int sitesPair[2];
  std::pair<int, int> sitesPair;
  // int site;
  struct LinkedListElementPair *next;
};
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////FUNCTION PROTOTYPES///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
bool directoryExists(const std::string& path) {
    return access(path.c_str(), F_OK | R_OK | X_OK) == 0;
}
//
void sitesXsitesYInit(const vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                    const vector<double>& siteComX, const vector<double>& siteComY, \
                    const int NumCells, const double AvgCellArea, const double Lx, const double Ly);
//
void energyEpsilonsFinder(double &E_test_eps, double &dE_eps,\
                          const vector<vector<double>>& J_int,\
                          const vector<vector<double>>& J_ext,\
                          const double Lambda, const double p);
//
struct Energy wholeSystemEnergyCalc(const vector<int>& sigmaMat,\ 
                                    const vector<vector<int>>& neighborsList, const vector<int>& neighborNum,\ 
                                    const vector<double>& cellArea,\
                                    const double Alpha, const double Lambda, const double AvgCellArea, const int NumCells);
//
struct Energy wholeSystemEnergyCalc_actPlus(const vector<int>& sigmaMat,\ 
                                    const vector<vector<int>>& neighborsList, const vector<int>& neighborNum,\ 
                                    const vector<double>& cellArea,\
                                    const double Alpha, const double Lambda, const double AvgCellArea, const int NumCells, \
                                    const double p, const vector<double>& theta, \
                                    const vector<double>& comX, const vector<double>& comY);
//
struct Energy wholeSystemEnergyCalc_actPlus_W(const vector<int>& sigmaMat,\ 
                                    const vector<vector<int>>& neighborsList, const vector<int>& neighborNum,\ 
                                    const vector<vector<double>>& edges, const double avgEdge,\ 
                                    const vector<double>& cellArea,\
                                    const double Alpha, const double Lambda, const double AvgCellArea, const int NumCells, \
                                    const double p, const vector<double>& theta, \
                                    const vector<double>& comX, const vector<double>& comY);
//
struct Energy energyCalcCompart(const vector<int>& sigmaMat, const vector<int>& compartMat,\
                                const vector<vector<double>>& J_int, const vector<vector<double>>& J_ext,\
                                const vector<vector<int>>& neighborsList, const vector<int>& neighborNum,\ 
                                const vector<vector<double>>& edges, const double avgEdge,\ 
                                const vector<vector<double>>& compartArea,\
                                const double Lambda, const vector<double>& avgArea, const int NumCells, \
                                const double p, const vector<double>& theta, \
                                const vector<double>& comX, const vector<double>& comY);
//
void LinkedListCreation(vector<LinkedListElementPair*>& LLheadVec, vector<LinkedListElementPair*>& LLtailVec, vector<int>& LLNumVec,\
                        const vector<int>& sigmaMat, const vector<int>& cellSiteNum, \
                        const vector<int>& neighborNum, const vector<vector<int>>& neighborsList,\
                        const vector<double>& siteComX, const double Lx);
//
void readCompartData(
    const std::string& filename,
    std::vector<double>& avgAreaFrac,
    std::vector<std::vector<double>>& J_int,
    std::vector<std::vector<double>>& J_ext
);
//
void initial_LS_saver(const string FolderName, const std::mt19937 &mt_rand);
//
void LinkedListSiteCreation(vector<LinkedListElement*>& LLheadVec, vector<LinkedListElement*>& LLtailVec, vector<int>& LLNumVec,\
                            const vector<int>& sigmaMat, vector<int>& isBorder,\
                            const vector<int>& neighborNum, const vector<vector<int>>& neighborsList,\
                            const vector<double>& siteComX, const double Lx);
//
void initializer(vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                 const vector<double>& siteComX, const vector<double>& siteComY, \
                 int NumCells, double AvgCellArea, double Lx, double Ly);
//
void vacInitializer(vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                    const vector<double>& siteComX, const vector<double>& siteComY, \
                    int NumCells, double AvgCellArea, double Lx, double Ly);
//
void recConfluentInitializer(vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                    const vector<double>& siteComX, const vector<double>& siteComY, \
                    int NumCells, double AvgCellArea, double Lx, double Ly);
//
void hexInitializer(vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                    const vector<double>& siteComX, const vector<double>& siteComY, \
                    int NumCells, double AvgCellArea, double Lx, double Ly);
//
void singleInitializer(std::mt19937 &mt_rand, vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                    const vector<double>& siteComX, const vector<double>& siteComY, \
                    int NumCells, double AvgCellArea, double Lx, double Ly);
//
void crystalInitializer(vector<vector<int>>& sigmaMat, vector<vector<int>>& sitesX, vector<vector<int>>& sitesY, int NumCells, double AvgCellArea);
//
void compartInitializer(std::mt19937 &mt_rand, const vector<int>& sigmaMat, vector<int>& compartMat, const vector<double>& avgAreaFrac);
//
void areaCalc(const vector<int>& sigmaMat, vector<double>& cellArea, vector<int>& cellSiteNum, const vector<double>& latticeArea);
//
void compartAreaCalc(const vector<int>& sigmaMat, const vector<int>& compartMat, vector<vector<double>>& compartArea, const vector<double>& latticeArea);
//
void saveSigmaMatCSV(const vector<vector<int>>& sigmaMat, const std::string &filename);
//
void saveDoubleMatCSV(const vector<vector<double>>& sigmaMat, const std::string &filename);
//
void saveIntMatCSV(const vector<vector<int>>& sigmaMat, const std::string &filename);
//
void saveSampleAP(const vector<vector<int>>& AP, const std::string &filename);
//
void saveSampleT(const vector<int>& quantityTemporal, const std::string &filename);
//
void saveSampleE(const vector<double>& quantityTemporal, const std::string &filename);
//
void saveSnapshots(const vector<vector<vector<int>>>& snapshots, const std::string &filename);
//
// void comCalc(const int sigmaMat[L][L], const int sitesX[L][L], const int sitesY[L][L], const int cellArea[NumCells + 1],
//              double comXDummy[NumCells + 1], double comYDummy[NumCells + 1]);
//
void comCalcIrr(const vector<int>& sigmaMat, const vector<double>& sitesX, const vector<double>& sitesY, const vector<double>& cellArea, \
                const vector<double>& latticeArea, \
                vector<double>& cellComX, vector<double>& cellComY);
//
void sigmaMatLoader(vector<vector<int>>& mat, const std::string &filename, int rows, int cols);
//
void samplingPatternLoader(vector<vector<long>>& mat, const std::string &filename, int rows, int cols);
//
inline bool existFunc (const std::string& name);
//
void simulationDataReader(int* NSitesPtr, double* LxPtr, double* LyPtr, double* AlphaPtr, double* KbPtr, double* TemPtr,  int* NumCellsPtr, \
                          double* AvgCellAreaPtr, double* LambdaPtr, long* maxMCStepsPtr, int* samplesPerWritePtr, \
                          int* printingTimeIntervalPtr, int* numLinkedListPtr, string* initConfigPtr);
//
void activityDataReader(double* DrPtr, double* pPtr);
//
int samplingPatternLenghFunc();
//
void latticeCreatorPhase1(int* maxNeighPtr, vector<int>& neighborNum, vector<double>& latticeX,  \
                          vector<double>& latticeY,  vector<double>& latticeArea );
//
void latticeCreatorPhase2(const int NSites, const double Lx, const double Ly, const int maxNeighbors,\
                          vector<vector<int>>& neighborsList, vector<vector<double>>& edges, double* avgEdgePtr, \
                          vector<vector<double>>& vorDist);
//
void polygonCom(double* xComPtr, double* yComPtr, const vector<vector<double>> polygonVec);
//
void latticeCreatorPhase3(const int NSites, const double Lx, const double Ly, const int maxNeighbors, \
                          const vector<int> neighborNum, const vector<vector<int>> neighborsList,\
                          vector<double>& siteComX, vector<double>& siteComY,\
                          vector<vector<double>>& deltaComX, vector<vector<double>>& deltaComY, vector<vector<double>>& comDist);
//
void saveInt1DVec(const vector<int>& quantity, const std::string &filename);
//
void saveDbl1DVec(const vector<double>& quantity, const std::string &filename);
//
void loadInt1DVec(vector<int>& mat, const std::string &filename);
//
void loadDbl1DVec(vector<double>& mat, const std::string &filename);
//
void LL_based();
//
void no_LL();
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
  #define OS "Windows"
#elif __APPLE__
  #define OS "MacOS"
#elif __linux__
  #define OS "Linux"
#else
  #define OS "Unknown"
#endif

int main()
{
  int LL_switch;
  
  std::ifstream LL_switch_file("LL_switch.csv");
  LL_switch_file >> LL_switch;
  LL_switch_file.close();

  if (LL_switch)
  {
    LL_based();
  }
  else
  {
    no_LL();
  }
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// FUNCTIONS ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void no_LL()
{

  auto start_time = high_resolution_clock::now();

  int NSites;
  double Lx, Ly;
  double Alpha = 0.0 / 0.0; // I want it to be NaN, because J matrices include the information
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
  
  simulationDataReader(&NSites, &Lx, &Ly, &Alpha, &Kb, &Tem,  &NumCells, \
                       &AvgCellArea, &Lambda, &maxMCSteps, &samplesPerWrite, &printingTimeInterval, &numLinkedList, &initConfig);
                     
  SweepLength = NSites;

  // Reading compartments data //
  std::vector<double> avgAreaFrac;
  std::vector<std::vector<double>> J_int, J_ext;
  readCompartData("compart_params.txt", avgAreaFrac, J_int, J_ext);
  int NComparts = avgAreaFrac.size();
  vector<double> avgArea(NComparts);
  for (int i = 0; i < NComparts; i++)
  {
    // average area of each compartment
    // 0= cytoplasm ; 1= Fz; 2= Vang;
    avgArea[i] = avgAreaFrac[i] * AvgCellArea;
  }
  // Reading compartments data //

  int writePerZip=0;
  std::ifstream writePerZipFile("writePerZip.csv");
  writePerZipFile >> writePerZip;
  writePerZipFile.close();


  /////////// ACTIVITY DATA //////////////////
  double Dr;
  double p;
  activityDataReader(&Dr, &p);
  /////////// ACTIVITY DATA //////////////////

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
  double avgEdge;
  latticeCreatorPhase2(NSites, Lx, Ly, maxNeighbors,\
                       neighborsList, edges, &avgEdge, vorDist);

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

  
  
  ////////////////////////////////////////////////////
  /////////////////// DEFENITIONS ////////////////////
  vector<int> sigmaMat(NSites); // The main field of spins
  vector<double> sitesX(NSites); // The value of the X(t) of each site
  vector<double> sitesY(NSites); // The value of the Y(t) of each site

  vector<int> compartMat(NSites); // The compartments of cells

  vector<vector<int> > sigmaSamples(NSites, vector<int>(samplesPerWrite));
  vector<vector<int> > compartSamples(NSites, vector<int>(samplesPerWrite));
  // vector<vector<double> > xSamples(NSites, vector<double>(samplesPerWrite));
  // vector<vector<double> > ySamples(NSites, vector<double>(samplesPerWrite));
  vector<vector<int> > xSamples(NSites, vector<int>(samplesPerWrite));
  vector<vector<int> > ySamples(NSites, vector<int>(samplesPerWrite));

  vector<double> cellArea(NumCells + 1); // The array showing the area of each cell (NOTE: values are from
                                         // the index 1 to NumCells, inclusive)
  vector<vector<double>> compartArea(NumCells + 1, vector<double>(NComparts));// The array showing the area of each compartment
  for (int compart_c = 0; compart_c < NComparts; compart_c++) // index 0 is for the void space (which we dont have here)
  {
    compartArea[0][compart_c] = 0;
  }
  cellArea[0] = 0;
  
  vector<vector<double>> areaSamples(NumCells + 1, vector<double>(samplesPerWrite)); // samples for
  for (int i=0; i<samplesPerWrite; i++)
  {
    areaSamples[0][i]=0;
  }
  vector<int> cellSiteNum(NumCells + 1); // The array showing the number of occupied sites by each cell (NOTE: values are from
  cellSiteNum[0] = 0;                    // the index 1 to NumCells, inclusive)

  vector<double> cellAreaTest(NumCells + 1); // The array showing the area of each cell (NOTE: values are from
                                            // the index 1 to NumCells, inclusive)
  vector<vector<double>> compartAreaTest(NumCells + 1, vector<double>(NComparts)); // The array showing the area of each compartment (for test)
  for (int compart_c = 0; compart_c < NComparts; compart_c++) // index 0 is for the void space (which we dont have here)
  {
    compartAreaTest[0][compart_c] = 0;
  }
  cellAreaTest[0] = 0; 

  vector<int> cellSiteNumTest(NumCells + 1); // The array showing the number of occupied sites by each cell (NOTE: values are from
  cellSiteNumTest[0] = 0;                    // the index 1 to NumCells, inclusive)


  struct Energy e0;
  double E, dE;
  vector<double> eSamples(samplesPerWrite);

  int t;
  vector<int> tSamples(samplesPerWrite);


  // activity vectors of the cells:
  vector<double> theta(NumCells + 1); // activity angle
  vector<double> n_x(NumCells + 1); // = cos(theta)
  vector<double> n_y(NumCells + 1); // = sin(theta)
  // activity vectors of the cells.

  /////////////////// DEFENITIONS ////////////////////
  ////////////////////////////////////////////////////

  /////////////////// FOLDERS NAMES ////////////////////
  std::string dataFolderName = "data";
  std::string initFolderName = "init";
  std::string mainResumeFolderName = "main_resume";
  std::string backupResumeFolderName = "backup_resume";
  std::string loadFolderName;
  /////////////////// FOLDERS NAMES ////////////////////


  ////////////////////// SAMPLING PATTERN DEFINITION ///////////
  int SPMLLength;
  SPMLLength = samplingPatternLenghFunc(); // SPML: sampling Pattern Milestone List
  unsigned long samplingPatternMilestoneList[SPMLLength];
  unsigned long samplingIntervalList[SPMLLength];
  unsigned long MAXMILESTONE;
  std::ifstream samplingPatternFile("samplingPattern_vec.csv");
  if (samplingPatternFile.good())
  {
    vector<vector<long>> samplingPatternTotal(2, vector<long>(SPMLLength));
    samplingPatternLoader(samplingPatternTotal, "samplingPattern_vec.csv", 2, SPMLLength);

    for(int colC=0; colC<SPMLLength; colC++)
    {
      samplingPatternMilestoneList[colC] = samplingPatternTotal[0][colC];
      samplingIntervalList[colC] = samplingPatternTotal[1][colC];
    }
    MAXMILESTONE = samplingPatternMilestoneList[SPMLLength-1];
  }else
  {
    cout<<"Problem in openning samplingPattern_vec.csv"<<endl;
    cout<<"Program ended!"<<endl;
    exit(0);
  }
  //////////////////// SAMPLING PATTERN DEFINITION ///////////
  


  ////////////////////// DETERMINING LOADING SWITCH ///////////
  unsigned long mt_rand_seed;
  int tOld;
  int tLSValue;
  int writeCounter;
  int loadSwitch;
  std::ifstream tLS(mainResumeFolderName+"/"+"tLS.csv"); // All the files with LS at the end, are being saved in the direct folder
  if (!(tLS.good()))
  {
    //std::cout << "tLS File does not exist." << std::endl;
    loadSwitch = 0;
  }
  else
  {
    //std::cout << "tLS File exists." << std::endl;
    tLS >> tLSValue;
    tLS.close();

    if (tLSValue >= maxMCSteps)
    {
      cout<<"The whole simulation has been already done!"<<endl;
      cout<<"Program ended!"<<endl;
      exit(0);
    }
    

    if (tLSValue == 0)
    {
      loadSwitch = 0;
    }
    else
    {
      std::ifstream tLSCheck(mainResumeFolderName+"/"+"tLSCheck.csv");
      int tOldCheck;
      tLSCheck >> tOldCheck;
      tLSCheck.close();

      if(tOldCheck == tLSValue)
      {
        loadSwitch = 1;
      }else{

        // go to backup folder
        std::ifstream tLS(backupResumeFolderName+"/"+"tLS.csv");
        tLS >> tLSValue;
        tLS.close();

        std::ifstream tLSCheck(backupResumeFolderName+"/"+"tLSCheck.csv");
        tLSCheck >> tOldCheck;
        tLSCheck.close();

        if(!(tOldCheck == tLSValue))
        {
          cout<<"Something is wrong with Resume and Backup folders!"<<endl<<endl;
          exit(0);
        }else{

          if(tLSValue==0)
          {
            loadSwitch = 0;
          }else{
            loadSwitch = 2;
          }

        }
      }


      if (loadSwitch>0)
      {
        char currentPath[500];
        getcwd(currentPath, 500);
        string currentPathStr = currentPath;
        std::string directoryToCheck =currentPathStr+"/"+dataFolderName;

        // std::filesystem::path currentPath = std::filesystem::current_path();
        // std::filesystem::path dataDirectory = currentPath / dataFolderName;

        // if (std::filesystem::is_directory(dataDirectory)) {
        if (directoryExists(directoryToCheck)) {
          cout << "data Folder exists!"<<endl;
          cout << "Data will be saved in this folder!"<<endl;
          cout << "Resuming..."<<endl;

        } else if(existFunc("data.zip"))
        {
          cout << "data Folder does NOT exist!"<<endl;
          cout << "data.zip exists."<<endl;
          cout << "Extracting..."<<endl;
          system("python3 dataUnzipper.py data");
          cout << "Resuming..."<<endl;
        }else
        {
          cout << "Neither data.zip exists nor data folder"<<endl;
          cout << "Reset all the settings and start the simulation from beginning!"<<endl;
          cout << "Program ended!"<<endl;
          exit(0);
        }

        // cout<<"\n"<<"\n"<<directoryToCheck<<endl<<endl;
        // exit(0);
      }

    }
  }
  ////////////////////// DETERMINING LOADING SWITCH ///////////

  ////////////////////// DEFINING RANDOM GENERATOR ///////////////
  std::mt19937 mt_rand;
  unsigned long MT_MAX = mt_rand.max();
  unsigned long MT_MIN = mt_rand.min();
  unsigned long randState;
  ////////////////////// DEFINING RANDOM GENERATOR ///////////////
  
  ////////////////////// LOADING OR INITIATING ////////////////
  if (loadSwitch == 0) // initialization must be done
  {
    /////////////////// MAKING SUB DIRECTORIES /////////////////
    // This block is for windows:
    // mkdir(dataFolderName.c_str()); //making data folder
    // mkdir(initFolderName.c_str()); //making init folder
    // mkdir(mainResumeFolderName.c_str()); //making main_resume folder
    // mkdir(backupResumeFolderName.c_str()); //making backup_resume folder

    
    // This block is for Linux:
    mkdir(dataFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making data folder
    mkdir(initFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making init folder
    mkdir(mainResumeFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making main_resume folder
    mkdir(backupResumeFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making backup_resume folder
    /////////////////// MAKING SUB DIRECTORIES /////////////////
    
    
    /////////////////// RANDOM GENERATOR SEEDING /////////////////
    random_device rd; // random seed creation
    mt_rand_seed = rd();
    // mt_rand_seed = 3758017832;

    //Seeding
    mt_rand.seed(mt_rand_seed);

    // saving initial random seed
    ofstream randSeedInit;
    randSeedInit.open(initFolderName + "/" + "randSeedInit.csv");
    randSeedInit << mt_rand_seed;
    randSeedInit.close(); // random seed saved
    
    // mt19937 mt_rand_init;
    // mt_rand_init.seed(mt_rand_seed);

    // saving initial random generator
    std::ofstream randStateInit(initFolderName + "/" + "randStateInit.csv");
    randStateInit << mt_rand;
    randStateInit.close();
    /////////////////// RANDOM GENERATOR SEEDING /////////////////


    /////////////////// SAVING INITIAL VALUES /////////////////

    // Initializing the sigmaMat
    if (initConfig.compare("c")==0)
    {
      // crystalInitializer(sigmaMat, sitesX, sitesY, NumCells, AvgCellArea);
      cout<<"#####################"<<endl;
      cout<<"No crystalized init for irr yet!"<<endl;
      cout<<"#####################"<<endl;
      exit(0);
    }
    else if (initConfig.compare("h")==0)
    {
      cout<<"********************"<<endl;
      cout<<"Hex initialization"<<endl;
      cout<<"********************"<<endl;
      hexInitializer(sigmaMat, sitesX, sitesY, \
                  siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
      compartInitializer(mt_rand, sigmaMat, compartMat, avgAreaFrac); // initialization of compartments
    }
    
    else if (initConfig.compare("r")==0)
    {
      // if (NumCells == 1)
      // {
      //   vacInitializer(sigmaMat, sitesX, sitesY, \
      //             siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
      // }
      // else
      // {
      //   initializer(sigmaMat, sitesX, sitesY, \
      //             siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
      // }

      recConfluentInitializer(sigmaMat, sitesX, sitesY, \
                  siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
      compartInitializer(mt_rand, sigmaMat, compartMat, avgAreaFrac); // initialization of compartments
      
    }
    else if (initConfig.compare("s")==0)
    {
      if (NumCells>1)
      {
        cout<<"#####################"<<endl;
        cout<<"If NumCells>1, singleInit is not accepted!"<<endl;
        cout<<"#####################"<<endl;
        exit(0);
      }

      singleInitializer(mt_rand, sigmaMat, sitesX, sitesY, \
                  siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
    }
    else if (initConfig.compare("p")==0) // provided sigmaMat and compartMat
    {
      loadInt1DVec(sigmaMat,   initFolderName+"/"+"sigma_init.csv");
      loadInt1DVec(compartMat, initFolderName+"/"+"compart_init.csv");
      sitesXsitesYInit(sigmaMat, sitesX, sitesY, \
                       siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
    }
    else
    {
      cout << "Please correct initConfig in simulationData_vec.csv"<<endl;
      cout << "Program ended!"<<endl;
      exit(0);
    }

    areaCalc(sigmaMat, cellArea, cellSiteNum, latticeArea);
    compartAreaCalc(sigmaMat, compartMat, compartArea, latticeArea);

    // for (int i = 0; i < NumCells; i++)
    // {
    //   if(abs(cellArea[i]-compartArea[i][0]-compartArea[i][1]-compartArea[i][2])>0.00001)
    //   {
    //     cout<<" insonsistency";
    //   }
    // }
    

    ///////////////// uncommented this part for activity /////////////////
    // Activity initialization
    // mt_rand.seed(1111);
    theta[0] = 0;
    for (int cellIndex =1; cellIndex<=NumCells; cellIndex++)
    {
      // randState = mt_rand();
      theta[cellIndex] = 2.0*M_PI*( ((long double)(mt_rand())-MT_MIN)/((long double)MT_MAX-MT_MIN) ) - M_PI;
      // theta[cellIndex] = 0.;
    }
    // saveThetaCSV(theta, initFolderName+"/"+"theta_init.csv");
    saveDbl1DVec(theta, initFolderName+"/"+"theta_init.csv");
    ///////////////// uncommented this part for activity /////////////////

    // saving the initial configuration
    saveInt1DVec(sigmaMat, initFolderName + "/" + "sigma_init.csv");
    saveInt1DVec(compartMat, initFolderName + "/" + "compart_init.csv");
    saveDbl1DVec(sitesX, initFolderName + "/" + "sitesX_init.csv");
    saveDbl1DVec(sitesY, initFolderName + "/" + "sitesY_init.csv");

    ///////////////// uncommented this part for activity /////////////////
    // initial com
    vector<double> comX_init(NumCells + 1);
    vector<double> comY_init(NumCells + 1);
    // comCalc(sigmaMat, sitesX, sitesY, cellArea, comX_init, comY_init);
    comCalcIrr(sigmaMat, sitesX, sitesY, cellArea, latticeArea, comX_init, comY_init);
    ///////////////// uncommented this part for activity /////////////////

    // e0 = wholeSystemEnergyCalc(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea, NumCells);
    // e0 = wholeSystemEnergyCalc_actPlus(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea, NumCells, p, theta, comX_init, comY_init);
    // e0 = wholeSystemEnergyCalc_actPlus_W(sigmaMat, neighborsList, neighborNum, edges, avgEdge, cellArea, Alpha, Lambda, AvgCellArea, NumCells, p, theta, comX_init, comY_init);
    e0 = energyCalcCompart(sigmaMat, compartMat, J_int, J_ext,
                           neighborsList, neighborNum, edges, avgEdge, 
                           compartArea, Lambda, avgArea, NumCells, 
                           p, theta, comX_init, comY_init);
    ofstream e0Initial;
    e0Initial.open(initFolderName + "/" + "e0Initial.csv");
    e0Initial << e0.total;
    e0Initial.close();
    E = e0.total;
    /////////////////// SAVING INITIAL VALUES /////////////////

    
    // initial LS saving
    ofstream tLSCheck;
    ofstream randStateLS;
    ofstream roundCounterLS;
    ofstream tLS;
    initial_LS_saver(mainResumeFolderName, mt_rand); //SAVING IN MAIN RESUME//
    initial_LS_saver(backupResumeFolderName, mt_rand); //SAVING IN BACKUP RESUME//
    tLSValue = 0;
    // initial LS saving

    // loadFolderName = mainResumeFolderName;

    writeCounter = 0;

  }
  else // loadSwitch != 0
  {
    if(loadSwitch == 1)
    {
      loadFolderName = mainResumeFolderName;
    }
    else if(loadSwitch == 2)
    {
      loadFolderName = backupResumeFolderName;
    }

    // loading oldRoundCounter to be zero
    int roundCounterOld;
    std::ifstream roundCounterLS(loadFolderName+"/"+"roundCounterLS.csv");
    roundCounterLS >> roundCounterOld;
    roundCounterLS.close();
    writeCounter = roundCounterOld+1;

    // LOADING RANDOM GENERATOR STATE LAST SAVED
    std::ifstream randStateLS(loadFolderName+"/"+"randStateLS.csv");
    randStateLS >> mt_rand;
    randStateLS.close();

    loadInt1DVec(sigmaMat, loadFolderName+"/"+"sigmaMatLS.csv");
    loadDbl1DVec(sitesX, loadFolderName+"/"+"sitesX_LS.csv");
    loadDbl1DVec(sitesY, loadFolderName+"/"+"sitesY_LS.csv");

    // reproducing the final last saved areas amd COM coordinates
    areaCalc(sigmaMat, cellArea, cellSiteNum, latticeArea);
    compartAreaCalc(sigmaMat, compartMat, compartArea, latticeArea);


    // calculate comXY based on loaded configuration; I need them for activity energy term
    vector<double> comX_load(NumCells + 1);
    vector<double> comY_load(NumCells + 1);
    // comCalc(sigmaMat, sitesX, sitesY, cellArea, comX_load, comY_load);
    comCalcIrr(sigmaMat, sitesX, sitesY, cellArea, latticeArea, comX_load, comY_load);

    // loading thetaLS
    // thetaLoader(theta, loadFolderName+"/"+"thetaLS.csv");
    loadDbl1DVec(theta, loadFolderName+"/"+"thetaLS.csv");

    // e0 = wholeSystemEnergyCalc(sigmaMat, cellArea, Alpha, Lambda, AvgCellArea); // Calculationg final last saved energy
    // e0 = wholeSystemEnergyCalc(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea, NumCells);
    // e0 = wholeSystemEnergyCalc_actPlus(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea, NumCells, p, theta, comX_load, comY_load);
    // e0 = wholeSystemEnergyCalc_actPlus_W(sigmaMat, neighborsList, neighborNum, edges, avgEdge, cellArea, Alpha, Lambda, AvgCellArea, NumCells, p, theta, comX_load, comY_load);
    e0 = energyCalcCompart(sigmaMat, compartMat, J_int, J_ext,
                           neighborsList, neighborNum, edges, avgEdge, 
                           compartArea, Lambda, avgArea, NumCells, 
                           p, theta, comX_load, comY_load);
    E = e0.total;
    

    tLS.open(loadFolderName+"/"+"tLS.csv");
    tLS >> tLSValue;
    tLS.close();
  }
  ////////////////////// LOADING OR INITIATING ////////////////

  int sampleCounter = 0;
  int cellIndex;
  tOld = tLSValue;
  t = tOld + 1;

  long attemptCounter;
  long indAttFrom; // randomly selected positions to be copied FROM
  long indAttInto; // randomly selected positions to be copied INTO
  // long pairIndex, pairIndexTemp, elementInd, elementC;
  long siteIndex, siteIndexTemp, elementInd, elementC;
  int sigmaFrom, sigmaInto; // value in the site FROM and INTO
  int flag;
  int dFNInto;               // delta FOREIGN NEIGHBOR INTO
  vector<int> variation_vec(maxNeighbors); // stores the variation of spin equality
  double dEInter, dEArea , dEAct; // INTERACTION and ELASTIC parts of energy
  int siteC, NNeigh, neighIndC, neighInd, neighIndCFromInto, NNeighNeigh, neighNeighIndC;
  


  int samplingPatternCounter = 0;
  for (int i=1; i<SPMLLength; i++)
  {
    int condition = ((t>=samplingPatternMilestoneList[i-1]) && (t<samplingPatternMilestoneList[i]));
    if(condition)
    {
      samplingPatternCounter=i;
      break;
    }
  }
  long samplingInterval = samplingIntervalList[samplingPatternCounter];
  long tLastSamplingPatternMilestone = 0;
  unsigned long tNextSamplingPatternMilestone = samplingPatternMilestoneList[samplingPatternCounter];

  int timeLoopFlag;
  timeLoopFlag=0;

  struct Energy E_test;
  struct Energy E_test_2;
  // double E_test_eps = 0.01*Alpha;
  double E_test_eps;
  double dE_eps;
  energyEpsilonsFinder(E_test_eps, dE_eps, J_int, J_ext, Lambda, p);

  

  ///////// RELATED TO CTIVITY
  vector<double> comX(NumCells + 1); // center of mass of the cells (X)
  vector<double> comY(NumCells + 1); // center of mass of the cells (Y)
  // comCalc(sigmaMat, sitesX, sitesY, cellArea, comX, comY);
  comCalcIrr(sigmaMat, sitesX, sitesY, cellArea, latticeArea, comX, comY);
  double cellAreaFromUpdated, cellAreaIntoUpdated;

  vector<double> comX_test(NumCells + 1);
  vector<double> comY_test(NumCells + 1);
  double com_test_eps = 0.0001;
  int com_test_flag = 0;

  int cellC;
  double thetaUpdateTerm;
  double eta;
  double SIG = sqrt(2.0*Dr); // standard dev of theta fluctuation
  double u1, u2; // for Box-Muller transformation

  double comX_old_Into, comY_old_Into, area_old_Into;
  double comX_new_Into, comY_new_Into, area_new_Into;
  double comX_old_From, comY_old_From, area_old_From;
  double comX_new_From, comY_new_From, area_new_From;
  double area_att; // the area of the site whose index is going to be changed in the current copying attempt
  double newX, newY;

  for (int cellC=1; cellC<=NumCells; cellC++)
  {
    n_x[cellC] = cos(theta[cellC]);
    n_y[cellC] = sin(theta[cellC]);
  }
  // double theta_eps = 1e-6;
  ///////// RELATED TO CTIVITY


  int compartOld, compartNew;
  std::vector<std::vector<std::vector<double>>> J_tot{ J_int, J_ext }; // combine all J matrices together
  int diff_cells_key_att;
  int diff_cells_key_old, diff_cells_key_new;


  while (tLSValue < maxMCSteps)
  {
    timeLoopFlag=1;
    for (attemptCounter = 0; attemptCounter < NSites; attemptCounter++)
    {
      do
      {        
        indAttInto = (mt_rand()) % NSites;
        
        neighIndC = (mt_rand()) % neighborNum[indAttInto];
        indAttFrom = neighborsList[indAttInto][neighIndC];

      } while (cellSiteNum[sigmaMat[indAttInto]]<=1);
      
      sigmaFrom = sigmaMat[indAttFrom];
      sigmaInto = sigmaMat[indAttInto];
      neighIndCFromInto = neighIndC; //BE CAREFUL: this is the index of 'indAttFrom' among the list of the neighbors of 'indAttInto'.

      diff_cells_key_att = bool(sigmaFrom - sigmaInto); // 1 if different, 0 if the same

      compartOld = compartMat[indAttInto];
      compartNew = (mt_rand()) % NComparts; // candidate for the compart of Into, in copying

      if (sigmaInto==sigmaFrom && compartOld==compartNew) 
      {// the same cell, the same compartment. Nothing to update
        break;
      } // Otherwise, updates must be done
      

      // for dEInter
      // dFNInto = 0;
      NNeigh = neighborNum[indAttInto];
      dEInter = 0.0;
      for (neighIndC = 0; neighIndC < NNeigh; neighIndC++)
      {
        neighInd = neighborsList[indAttInto][neighIndC];
        diff_cells_key_old = bool(sigmaMat[neighInd] - sigmaInto); // 1 if different, 0 if the same
        diff_cells_key_new = bool(sigmaMat[neighInd] - sigmaFrom); // 1 if different, 0 if the same
        dEInter += (J_tot[diff_cells_key_new][compartNew][compartMat[neighInd]] \
                  - J_tot[diff_cells_key_old][compartOld][compartMat[neighInd]] )\
                 * (edges[indAttInto][neighIndC] / avgEdge);
      }
      
      // dEInter = Alpha * dFNInto;
      // for dEInter
      

      area_att = latticeArea[indAttInto];
      
      // for dEAct
      
      comX_old_Into = comX[sigmaInto];
      comY_old_Into = comY[sigmaInto];
      area_old_Into = cellArea[sigmaInto];

      comX_old_From = comX[sigmaFrom];
      comY_old_From = comY[sigmaFrom];
      area_old_From = cellArea[sigmaFrom];
      
      if (diff_cells_key_att) { // different cells
        comX_new_Into = (comX_old_Into * area_old_Into - area_att * sitesX[indAttInto]) / (area_old_Into - area_att);
        comY_new_Into = (comY_old_Into * area_old_Into - area_att * sitesY[indAttInto]) / (area_old_Into - area_att);

        newX = sitesX[indAttFrom] - deltaComX[indAttInto][neighIndCFromInto]; //BE CAREFUL: 'neighIndCFromInto' is the index of 'indAttFrom' among the list of the neighbors of 'indAttInto'.
        newY = sitesY[indAttFrom] - deltaComY[indAttInto][neighIndCFromInto]; //BE CAREFUL: 'neighIndCFromInto' is the index of 'indAttFrom' among the list of the neighbors of 'indAttInto'.

        comX_new_From = (comX_old_From * area_old_From + area_att * newX) / (area_old_From + area_att);
        comY_new_From = (comY_old_From * area_old_From + area_att * newY) / (area_old_From + area_att);

        dEAct = - p * ( \
        n_x[sigmaFrom] * (comX_new_From - comX_old_From) + \
        n_y[sigmaFrom] * (comY_new_From - comY_old_From) + \
        n_x[sigmaInto] * (comX_new_Into - comX_old_Into) + \
        n_y[sigmaInto] * (comY_new_Into - comY_old_Into) );

      } else { // the same cell
        comX_new_Into = comX[sigmaInto];
        comY_new_Into = comY[sigmaInto];
        // although sigmaInto == sigmaFrom here
        
        newX = sitesX[indAttInto]; //BE CAREFUL: 'neighIndCFromInto' is the index of 'indAttFrom' among the list of the neighbors of 'indAttInto'.
        newY = sitesY[indAttInto]; //BE CAREFUL: 'neighIndCFromInto' is the index of 'indAttFrom' among the list of the neighbors of 'indAttInto'.

        comX_new_From = comX[sigmaFrom];
        comY_new_From = comY[sigmaFrom];

        dEAct = 0;
        }
      // for dEAct
      
      // for dEArea
      if (sigmaFrom*sigmaInto > 0)
        {
          // dEArea = 2.0 * Lambda * area_att * (area_att + area_old_From - area_old_Into); // without compartment

          // // with compartment; detailed:
          // dEArea = - Lambda * pow(            compartArea[sigmaInto][compartOld] - avgArea[compartOld], 2) \
          //          + Lambda * pow(-area_att + compartArea[sigmaInto][compartOld] - avgArea[compartOld], 2) \
          //          - Lambda * pow(            compartArea[sigmaFrom][compartNew] - avgArea[compartNew], 2) \
          //          + Lambda * pow(+area_att + compartArea[sigmaFrom][compartNew] - avgArea[compartNew], 2) ;

          // with compartment; simplified: (calculations in the notebook)
          dEArea = 2 * Lambda * area_att * ( area_att \
          -compartArea[sigmaInto][compartOld] + avgArea[compartOld] \
          +compartArea[sigmaFrom][compartNew] - avgArea[compartNew] );
        }
      else if (sigmaInto==0)
        {
          // dEArea =  Lambda * area_att * (area_att + 2 * ( area_old_From - AvgCellArea));
          // correct this case later
        }
      else if (sigmaFrom==0)
        {
          // dEArea =  Lambda * area_att * (area_att - 2 * ( area_old_Into - AvgCellArea));
          // correct this case later
        }
      // for dEArea



      // dE = dE_inter + dE_area + dE_act;
      dE = dEInter + dEArea + dEAct;

      // E_test = energyCalcCompart(sigmaMat, compartMat, J_int, J_ext,
      //                       neighborsList, neighborNum, edges, avgEdge, 
      //                       compartArea, Lambda, avgArea, NumCells, 
      //                       p, theta, comX, comY);

      if (  (dE < dE_eps) || ((((long double)(mt_rand())-MT_MIN)/((long double)MT_MAX-MT_MIN)) < 
                  exp(-dE/(Tem*Kb)) )   )
      {

        E += dE;
        

        // Area update
        cellArea[sigmaFrom] += area_att;
        cellArea[sigmaInto] -= area_att;
        cellSiteNum[sigmaFrom]++;
        cellSiteNum[sigmaInto]--;
        compartArea[sigmaFrom][compartNew] += area_att;
        compartArea[sigmaInto][compartOld] -= area_att;
        
        //sigmaMat update
        sigmaMat[indAttInto] = sigmaFrom;
        compartMat[indAttInto] = compartNew;

        // E_test = wholeSystemEnergyCalc(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea);

        // comXY update INTO
        comX[sigmaInto] = comX_new_Into;
        comY[sigmaInto] = comY_new_Into;

        // siteXY update
        sitesX[indAttInto] = newX;
        sitesY[indAttInto] = newY;

        // comXY update FROM
        comX[sigmaFrom] = comX_new_From;
        comY[sigmaFrom] = comY_new_From;

        // E_test_2 = energyCalcCompart(sigmaMat, compartMat, J_int, J_ext,
        //                     neighborsList, neighborNum, edges, avgEdge, 
        //                     compartArea, Lambda, avgArea, NumCells, 
        //                     p, theta, comX, comY);

        // if (abs(E_test_2.inter-E_test.inter - dEInter) > (1e-8))
        // {
        //   int ttt =0;
        // }
        // if (abs(E_test_2.area-E_test.area - dEArea) > (1e-8))
        // {
        //   int ttt =0;
        // }
        
        
      }
    }

    //////////////////// UPDATING THETA //////////////////
    for (int cellC=1; cellC<=NumCells; cellC++)
    {

      randState = MT_MIN;
      while(randState == MT_MIN)
      {
        randState = mt_rand();
      }
      // randState = mt_rand();
      u1 = ((long double)(randState) - MT_MIN) / ((long double)MT_MAX - MT_MIN);

      randState = MT_MIN;
      while(randState == MT_MIN)
      {
        randState = mt_rand();
      }
      // randState = mt_rand();
      u2 = ((long double)(randState) - MT_MIN) / ((long double)MT_MAX - MT_MIN);

      eta = SIG * sqrt(-2.0*log(u1)) * cos(2.0*M_PI*u2);
      thetaUpdateTerm = eta * 1.00; // dt = 1.00

      E -= ((-p) * (n_x[cellC] * comX[cellC] + n_y[cellC] * comY[cellC]) );

      theta[cellC] = theta[cellC] + thetaUpdateTerm;

      while (theta[cellC] > M_PI)
      {
        theta[cellC] = theta[cellC] - 2.0*M_PI;
      }
      while (theta[cellC] < -M_PI)
      { 
        theta[cellC] = theta[cellC] + 2.0*M_PI;
      }

      n_x[cellC] = cos(theta[cellC]);
      n_y[cellC] = sin(theta[cellC]);

      E += ((-p) * (n_x[cellC] * comX[cellC] + n_y[cellC] * comY[cellC]) );
    }
    //////////////////// UPDATING THETA //////////////////


  /////////////////////// SAMPLING OPERATION ////////////////////////
    if ((t-tLastSamplingPatternMilestone) % samplingInterval == 0)
    {
      tSamples[sampleCounter] = t;
      eSamples[sampleCounter] = E;
      for (cellIndex = 0; cellIndex <= NumCells; cellIndex++)
      {
        areaSamples[cellIndex][sampleCounter] = cellArea[cellIndex];
      }
      for (siteC = 0; siteC < NSites; siteC++)
      {
        sigmaSamples[siteC][sampleCounter] = sigmaMat[siteC];
        compartSamples[siteC][sampleCounter] = compartMat[siteC];
        xSamples[siteC][sampleCounter] = floor(sitesX[siteC]/Lx);
        ySamples[siteC][sampleCounter] = floor(sitesY[siteC]/Ly);
      }
      sampleCounter++;

      /////////////////////// WRITING OPERATION ////////////////////////
      if (sampleCounter == samplesPerWrite)
      {
        
        // cout << "\033[2J\033[1;1H";
        // cout << "Saving! Do not terminate!"<<endl;

        ///////////////////////// DATA SAVING ///////////////////////
        // writing the center of mass of the cells
        saveIntMatCSV(sigmaSamples, dataFolderName + "/" + "sigmaSamples_" + to_string(writeCounter) + ".csv");
        saveIntMatCSV(compartSamples, dataFolderName + "/" + "compartSamples_" + to_string(writeCounter) + ".csv");
        // saveDoubleMatCSV(xSamples, dataFolderName + "/" + "xSamples_" + to_string(writeCounter) + ".csv");
        // saveDoubleMatCSV(ySamples, dataFolderName + "/" + "ySamples_" + to_string(writeCounter) + ".csv");
        saveIntMatCSV(xSamples, dataFolderName + "/" + "xSamples_" + to_string(writeCounter) + ".csv");
        saveIntMatCSV(ySamples, dataFolderName + "/" + "ySamples_" + to_string(writeCounter) + ".csv");

        // writing the area and perimiter of the cells
        // saveDoubleMatCSV(areaSamples, dataFolderName + "/" + "areaSamples_" + to_string(writeCounter) + ".csv");
        // saveSampleAP(sampleP, dataFolderName+"/"+"sampleP_"+to_string(writeCounter)+".csv");

        // writing the sampling times and energies
        saveSampleT(tSamples, dataFolderName + "/" + "tSamples_" + to_string(writeCounter) + ".csv");
        saveSampleE(eSamples, dataFolderName + "/" + "eSamples_" + to_string(writeCounter) + ".csv");
        ///////////////////////// DATA SAVING ///////////////////////


        ///////////// NUMBERED DATA ZIPPING ///////////////////
        if ((writeCounter+1)%writePerZip ==0)
        {
          // string argv1 = 'data_'+to_string((int)(writeCounter/writePerZip))+'.zip';
          std::string argv1 = "data_" + std::to_string(static_cast<int>(writeCounter / writePerZip)+1) + ".zip";
          std::string argv2 = dataFolderName;
          std::string command = "python3 dataZipperNum.py " + argv1 + " " + argv2;
          system(command.c_str());
        }
        ///////////// NUMBERED DATA ZIPPING ///////////////////

        
        ///////////////////////// LS SAVING ///////////////////////
        for(int resumeFolderCounter=1; resumeFolderCounter<=2; resumeFolderCounter++)
        {
          if (resumeFolderCounter==1)
          {
            loadFolderName = mainResumeFolderName;
          }
          else if (resumeFolderCounter==2)
          {
            loadFolderName = backupResumeFolderName;
          }

          // writing the final time
          ofstream tLSCheck;
          tLSCheck.open(loadFolderName+"/"+"tLSCheck.csv");
          tLSCheck << t;
          tLSCheck.close();

          // writing the final configuration of the sites
          saveInt1DVec(sigmaMat, loadFolderName+"/"+"sigmaMatLS.csv");
          saveInt1DVec(compartMat, loadFolderName+"/"+"compartMatLS.csv");
          saveDbl1DVec(sitesX, loadFolderName+"/"+"sitesX_LS.csv");
          saveDbl1DVec(sitesY, loadFolderName+"/"+"sitesY_LS.csv");

          saveDbl1DVec(theta, loadFolderName+"/"+"thetaLS.csv");

          // writing the final state of random generator
          std::ofstream randStateLS(loadFolderName+"/"+"randStateLS.csv");
          randStateLS << mt_rand;
          randStateLS.close();

          // writing round counter
          ofstream roundCounterLS;
          roundCounterLS.open(loadFolderName+"/"+"roundCounterLS.csv");
          roundCounterLS << writeCounter;
          roundCounterLS.close();

          // writing the final time
          ofstream tLS;
          tLS.open(loadFolderName+"/"+"tLS.csv");
          tLSValue = t;
          tLS << tLSValue;
          tLS.close();

          // cout << "\033[2J\033[1;1H";
        }
        ///////////////////////// LS SAVING ///////////////////////

        ///////////// TESTING (ENERGY & COM) ///////////////////
        // E_test = wholeSystemEnergyCalc(sigmaMat, cellArea, Alpha, Lambda, AvgCellArea); // Calculationg final last saved energy
        // E_test = wholeSystemEnergyCalc_actPlus(sigmaMat, cellArea, Alpha, Lambda, AvgCellArea, p, theta, comX, comY);
        // E_test = wholeSystemEnergyCalc(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea, NumCells);
        // E_test = wholeSystemEnergyCalc_actPlus(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea, NumCells, p, theta, comX, comY);
        // E_test = wholeSystemEnergyCalc_actPlus_W(sigmaMat, neighborsList, neighborNum, edges, avgEdge, cellArea, Alpha, Lambda, AvgCellArea, NumCells, p, theta, comX, comY);
        E_test = energyCalcCompart(sigmaMat, compartMat, J_int, J_ext,
                            neighborsList, neighborNum, edges, avgEdge, 
                            compartArea, Lambda, avgArea, NumCells, 
                            p, theta, comX, comY);

        if (fabs(E_test.total-E)>E_test_eps)
        {
          ofstream EnergyError;
          EnergyError.open("EnergyError.csv");
          EnergyError << "E: ";
          EnergyError << E;
          EnergyError << endl;
          EnergyError << "E_test: ";
          EnergyError << E_test.total;
          EnergyError.close();

          cout<<"Energy calculation inconsistency!"<<endl;
          cout<<"Look at the error file."<<endl;
          cout<<"Program ended!"<<endl;
          exit(0);
        } else
        {
          // cout<<"Energy successfully checked!"<<endl;
        }

        
        areaCalc(sigmaMat, cellAreaTest, cellSiteNumTest, latticeArea);
        compartAreaCalc(sigmaMat, compartMat, compartAreaTest, latticeArea);
        int areaFlag = 0;
        int areaErrorCondition = 0;

        for (int cellCTest = 0; cellCTest <= NumCells; cellCTest++)
        {
          areaErrorCondition += bool(fabs(cellAreaTest[cellCTest]-cellArea[cellCTest])>(1e-6));
          areaErrorCondition += bool((fabs(cellSiteNumTest[cellCTest]-cellSiteNum[cellCTest])>(1e-6)));
          for (int compartC = 0; compartC < NComparts; compartC++)
          {
            areaErrorCondition += bool(fabs(compartAreaTest[cellCTest][compartC]-compartArea[cellCTest][compartC])>(1e-6));
          }

          if (areaErrorCondition)
          {
            areaFlag = 1;

            ofstream AreaError;
            AreaError.open("AreaError.csv");
            // AreaError << "Area["<<cellCTest<<"]= "<<cellArea[cellCTest];
            // AreaError << endl;
            // AreaError << "AreaTest["<<cellCTest<<"]= "<<cellAreaTest[cellCTest];
            // AreaError.close();
            AreaError << "Area inconsistency in cell "<<cellCTest;
            AreaError << endl;
            AreaError.close();

            cout<<"Area calculation inconsistency!"<<endl;
            cout<<"Look at the error file."<<endl;
            cout<<"Program ended!"<<endl;
            exit(0);
          }
          
        }
        if (areaFlag==0)
        {
          // cout<<"Area successfully checked!"<<endl;
        }
        
        

        
        // comXY
        // comCalc(sigmaMat, sitesX, sitesY, cellArea, comX_test, comY_test);
        comCalcIrr(sigmaMat, sitesX, sitesY, cellArea, latticeArea, comX_test, comY_test);
        for (int dum_counter = 0; dum_counter <= NumCells; dum_counter++)
        {
          if( (fabs(comX[dum_counter]-comX_test[dum_counter])>com_test_eps) || (fabs(comY[dum_counter]-comY_test[dum_counter])>com_test_eps))
          {
            ofstream comError;
            comError.open("comError.csv");
            comError << "cell Index: ";
            comError << dum_counter;
            comError << endl;
            comError << "comX (update): ";
            comError << comX[dum_counter];
            comError << endl;
            comError << "comX_test (function): ";
            comError << comX_test[dum_counter];
            comError << endl;
            comError << "comY (update): ";
            comError << comY[dum_counter];
            comError << endl;
            comError << "comY_test (function): ";
            comError << comY_test[dum_counter];
            comError << endl;
            comError.close();

            cout<<"comXY calculation inconsistency!"<<endl;
            cout<<"Look at the error file."<<endl;
            cout<<"Program ended!"<<endl;
            com_test_flag = 1;
            exit(0);
          }
        }
        if(com_test_flag==0)
        {
          // cout<<"comXY successfully checked!"<<endl;
        }
        /////////// TESTING (ENERGY & COM) ///////////////////

        writeCounter++;
        sampleCounter = 0;
      }
      /////////////////////// WRITING OPERATION ////////////////////////


      /////////////////////// CHANGING SAMPLING PATTERN ////////////////////////
      if(t == tNextSamplingPatternMilestone)
      {
        samplingPatternCounter++;
        samplingInterval = samplingIntervalList[samplingPatternCounter];
        tLastSamplingPatternMilestone = samplingPatternMilestoneList[samplingPatternCounter-1];
        tNextSamplingPatternMilestone = samplingPatternMilestoneList[samplingPatternCounter];
      }
      /////////////////////// CHANGING SAMPLING PATTERN ////////////////////////
    }
    /////////////////////// SAMPLING OPERATION ////////////////////////

    if (!(t % printingTimeInterval))
    {
      // cout << "\033[2J\033[1;1H";
      cout << "t : " << t << "/" << maxMCSteps << endl;
    }

    t++;
  }
  t--;

  
  if(timeLoopFlag)
  {

    auto stop_time = high_resolution_clock::now();
    // Calculate the duration of the code execution in minutes
    auto duration = duration_cast<seconds>(stop_time - start_time);
    // Output the duration in minutes
    cout << "Time taken by function: " << duration.count() << " seconds" << endl;
    std::ofstream Time_outfile("time.dat");
    if (Time_outfile.is_open())
    {
      Time_outfile << duration.count();
      Time_outfile << std::endl;
      Time_outfile.close();
    }
    else
    {
      std::cerr << "Could not open file for writing" << std::endl;
    }
  }

  system("python3 dataZipper.py");

  // system("./pp_irr_v1");

  // system("python3 ppDataZipper.py");

  cout << "Program done!" << endl<< endl;

  // return 0;
}
//
void LL_based()
{

  auto start_time = high_resolution_clock::now();

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
  
  simulationDataReader(&NSites, &Lx, &Ly, &Alpha, &Kb, &Tem,  &NumCells, \
                       &AvgCellArea, &Lambda, &maxMCSteps, &samplesPerWrite, &printingTimeInterval, &numLinkedList, &initConfig);
  SweepLength = NSites;

  int writePerZip;
  std::ifstream writePerZipFile("writePerZip.csv");
  writePerZipFile >> writePerZip;
  writePerZipFile.close();

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
  double avgEdge;
  latticeCreatorPhase2(NSites, Lx, Ly, maxNeighbors,\
                       neighborsList, edges, &avgEdge, vorDist);

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

  
  
  ////////////////////////////////////////////////////
  /////////////////// DEFENITIONS ////////////////////
  vector<int> sigmaMat(NSites); // The main field of spins
  vector<double> sitesX(NSites); // The value of the X(t) of each site
  vector<double> sitesY(NSites); // The value of the Y(t) of each site

  vector<vector<int> > sigmaSamples(NSites, vector<int>(samplesPerWrite));
  // vector<vector<double> > xSamples(NSites, vector<double>(samplesPerWrite));
  // vector<vector<double> > ySamples(NSites, vector<double>(samplesPerWrite));
  vector<vector<int> > xSamples(NSites, vector<int>(samplesPerWrite));
  vector<vector<int> > ySamples(NSites, vector<int>(samplesPerWrite));

  vector<double> cellArea(NumCells + 1); // The array showing the area of each cell (NOTE: values are from
  cellArea[0] = 0;                       // the index 1 to NumCells, inclusive)
  vector<vector<double>> areaSamples(NumCells + 1, vector<double>(samplesPerWrite));
  for (int i=0; i<samplesPerWrite; i++)
  {
    areaSamples[0][i]=0;
  }
  vector<int> cellSiteNum(NumCells + 1); // The array showing the number of occupied sites by each cell (NOTE: values are from
  cellSiteNum[0] = 0;                    // the index 1 to NumCells, inclusive)

  vector<double> cellAreaTest(NumCells + 1); // The array showing the area of each cell (NOTE: values are from
  cellAreaTest[0] = 0;                       // the index 1 to NumCells, inclusive)

  vector<int> cellSiteNumTest(NumCells + 1); // The array showing the number of occupied sites by each cell (NOTE: values are from
  cellSiteNumTest[0] = 0;                    // the index 1 to NumCells, inclusive)

  struct Energy e0;
  double E, dE;
  vector<double> eSamples(samplesPerWrite);

  int t;
  vector<int> tSamples(samplesPerWrite);
  /////////////////// DEFENITIONS ////////////////////
  ////////////////////////////////////////////////////

  /////////////////// FOLDERS NAMES ////////////////////
  std::string dataFolderName = "data";
  std::string initFolderName = "init";
  std::string mainResumeFolderName = "main_resume";
  std::string backupResumeFolderName = "backup_resume";
  std::string loadFolderName;
  /////////////////// FOLDERS NAMES ////////////////////


  ////////////////////// SAMPLING PATTERN DEFINITION ///////////
  int SPMLLength;
  SPMLLength = samplingPatternLenghFunc(); // SPML: sampling Pattern Milestone List
  unsigned long samplingPatternMilestoneList[SPMLLength];
  unsigned long samplingIntervalList[SPMLLength];
  unsigned long MAXMILESTONE;
  std::ifstream samplingPatternFile("samplingPattern_vec.csv");
  if (samplingPatternFile.good())
  {
    vector<vector<long>> samplingPatternTotal(2, vector<long>(SPMLLength));
    samplingPatternLoader(samplingPatternTotal, "samplingPattern_vec.csv", 2, SPMLLength);

    for(int colC=0; colC<SPMLLength; colC++)
    {
      samplingPatternMilestoneList[colC] = samplingPatternTotal[0][colC];
      samplingIntervalList[colC] = samplingPatternTotal[1][colC];
    }
    MAXMILESTONE = samplingPatternMilestoneList[SPMLLength-1];
  }else
  {
    cout<<"Problem in openning samplingPattern_vec.csv"<<endl;
    cout<<"Program ended!"<<endl;
    exit(0);
  }
  //////////////////// SAMPLING PATTERN DEFINITION ///////////
  


  ////////////////////// DETERMINING LOADING SWITCH ///////////
  unsigned long mt_rand_seed;
  int tOld;
  int tLSValue;
  int writeCounter;
  int loadSwitch;
  std::ifstream tLS(mainResumeFolderName+"/"+"tLS.csv"); // All the files with LS at the end, are being saved in the direct folder
  if (!(tLS.good()))
  {
    //std::cout << "tLS File does not exist." << std::endl;
    loadSwitch = 0;
  }
  else
  {
    //std::cout << "tLS File exists." << std::endl;
    tLS >> tLSValue;
    tLS.close();

    if (tLSValue >= maxMCSteps)
    {
      cout<<"The whole simulation has been already done!"<<endl;
      cout<<"Program ended!"<<endl;
      exit(0);
    }
    

    if (tLSValue == 0)
    {
      loadSwitch = 0;
    }
    else
    {
      std::ifstream tLSCheck(mainResumeFolderName+"/"+"tLSCheck.csv");
      int tOldCheck;
      tLSCheck >> tOldCheck;
      tLSCheck.close();

      if(tOldCheck == tLSValue)
      {
        loadSwitch = 1;
      }else{

        // go to backup folder
        std::ifstream tLS(backupResumeFolderName+"/"+"tLS.csv");
        tLS >> tLSValue;
        tLS.close();

        std::ifstream tLSCheck(backupResumeFolderName+"/"+"tLSCheck.csv");
        tLSCheck >> tOldCheck;
        tLSCheck.close();

        if(!(tOldCheck == tLSValue))
        {
          cout<<"Something is wrong with Resume and Backup folders!"<<endl<<endl;
          exit(0);
        }else{

          if(tLSValue==0)
          {
            loadSwitch = 0;
          }else{
            loadSwitch = 2;
          }

        }
      }


      if (loadSwitch>0)
      {
        char currentPath[500];
        getcwd(currentPath, 500);
        string currentPathStr = currentPath;
        std::string directoryToCheck =currentPathStr+"/"+dataFolderName;

        // std::filesystem::path currentPath = std::filesystem::current_path();
        // std::filesystem::path dataDirectory = currentPath / dataFolderName;

        // if (std::filesystem::is_directory(dataDirectory)) {
        if (directoryExists(directoryToCheck)) {
          cout << "data Folder exists!"<<endl;
          cout << "Data will be saved in this folder!"<<endl;
          cout << "Resuming..."<<endl;

        } else if(existFunc("data.zip"))
        {
          cout << "data Folder does NOT exist!"<<endl;
          cout << "data.zip exists."<<endl;
          cout << "Extracting..."<<endl;
          system("python3 dataUnzipper.py");
          cout << "Resuming..."<<endl;
        }else
        {
          cout << "Neither data.zip exists nor data folder"<<endl;
          cout << "Reset all the settings and start the simulation from beginning!"<<endl;
          cout << "Program ended!"<<endl;
          exit(0);
        }
      }

    }
  }
  ////////////////////// DETERMINING LOADING SWITCH ///////////

  ////////////////////// DEFINING RANDOM GENERATOR ///////////////
  std::mt19937 mt_rand;
  unsigned long MT_MAX = mt_rand.max();
  unsigned long MT_MIN = mt_rand.min();
  unsigned long randState;
  ////////////////////// DEFINING RANDOM GENERATOR ///////////////
  
  ////////////////////// LOADING OR INITIATING ////////////////
  if (loadSwitch == 0)
  {
    /////////////////// MAKING SUB DIRECTORIES /////////////////

    // This block is for windows:
    // mkdir(dataFolderName.c_str()); //making data folder
    // mkdir(initFolderName.c_str()); //making init folder
    // mkdir(mainResumeFolderName.c_str()); //making main_resume folder
    // mkdir(backupResumeFolderName.c_str()); //making backup_resume folder

    // This block is for Linux:
    mkdir(dataFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making data folder
    mkdir(initFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making init folder
    mkdir(mainResumeFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making main_resume folder
    mkdir(backupResumeFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making backup_resume folder
    /////////////////// MAKING SUB DIRECTORIES /////////////////

    
    /////////////////// RANDOM GENERATOR SEEDING /////////////////
    random_device rd; // random seed creation
    mt_rand_seed = rd();
    // mt_rand_seed = 18;

    //Seeding
    mt_rand.seed(mt_rand_seed);

    // saving initial random seed
    ofstream randSeedInit;
    randSeedInit.open(initFolderName + "/" + "randSeedInit.csv");
    randSeedInit << mt_rand_seed;
    randSeedInit.close(); // random seed saved
    
    // mt19937 mt_rand_init;
    // mt_rand_init.seed(mt_rand_seed);

    // saving initial random generator
    std::ofstream randStateInit(initFolderName + "/" + "randStateInit.csv");
    randStateInit << mt_rand;
    randStateInit.close();
    /////////////////// RANDOM GENERATOR SEEDING /////////////////


    /////////////////// SAVING INITIAL VALUES /////////////////

    // Initializing the sigmaMat
    if (initConfig.compare("c")==0)
    {
      // crystalInitializer(sigmaMat, sitesX, sitesY, NumCells, AvgCellArea);
      cout<<"#####################"<<endl;
      cout<<"No crystalized init for irr yet!"<<endl;
      cout<<"#####################"<<endl;
      exit(0);
    }
    if (initConfig.compare("h")==0)
    {
      // crystalInitializer(sigmaMat, sitesX, sitesY, NumCells, AvgCellArea);
      cout<<"*********************"<<endl;
      cout<<"Hex initialization"<<endl;
      cout<<"*********************"<<endl;
      hexInitializer(sigmaMat, sitesX, sitesY, \
                  siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
    }
    else if (initConfig.compare("r")==0)
    {
      // if (NumCells == 1)
      // {
      //   vacInitializer(sigmaMat, sitesX, sitesY, \
      //             siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
      // }
      // else
      // {
      //   initializer(sigmaMat, sitesX, sitesY, \
      //             siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
      // }

      recConfluentInitializer(sigmaMat, sitesX, sitesY, \
                  siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
      
    }
    else if (initConfig.compare("s")==0)
    {
      if (NumCells>1)
      {
        cout<<"#####################"<<endl;
        cout<<"If NumCells>1, singleInit is not accepted!"<<endl;
        cout<<"#####################"<<endl;
        exit(0);
      }

      singleInitializer(mt_rand, sigmaMat, sitesX, sitesY, \
                  siteComX, siteComY, NumCells, AvgCellArea, Lx, Ly);
    }
    else
    {
      cout << "Please correct initConfig in simulationData_vec.csv"<<endl;
      cout << "Program ended!"<<endl;
      exit(0);
    }

    areaCalc(sigmaMat, cellArea, cellSiteNum, latticeArea);

    // // Activity initialization
    // // mt_rand.seed(1111);
    // theta[0] = 0;
    // for (int cellIndex =1; cellIndex<=NumCells; cellIndex++)
    // {
    //   // randState = mt_rand();
    //   theta[cellIndex] = 2.0*M_PI*( ((long double)(mt_rand())-MT_MIN)/((long double)MT_MAX-MT_MIN) );
    //   // theta[cellIndex] = 0.;
    // }
    // saveThetaCSV(theta, initFolderName+"/"+"theta_init.csv");

    // saving the initial configuration
    saveInt1DVec(sigmaMat, initFolderName + "/" + "sigma_init.csv");
    saveDbl1DVec(sitesX, initFolderName + "/" + "sitesX_init.csv");
    saveDbl1DVec(sitesY, initFolderName + "/" + "sitesY_init.csv");

    // // initial com
    // vector<double> comX_init(NumCells + 1);
    // vector<double> comY_init(NumCells + 1);
    // comCalc(sigmaMat, sitesX, sitesY, cellArea, comX_init, comY_init);

    e0 = wholeSystemEnergyCalc(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea, NumCells);
    ofstream e0Initial;
    e0Initial.open(initFolderName + "/" + "e0Initial.csv");
    e0Initial << e0.total;
    e0Initial.close();
    E = e0.total;
    /////////////////// SAVING INITIAL VALUES /////////////////

    
    //////////////////////////SAVING IN MAIN RESUME/////////////////////////////
    // saving the initial time
    ofstream tLSCheck;
    tLSCheck.open(mainResumeFolderName+"/"+"tLSCheck.csv");
    tLSCheck << 0;
    tLSCheck.close();

    // saving LS random generator
    //std::ofstream randStateLS(mainResumeFolderName+"/"+"randStateLS.csv");
    ofstream randStateLS;
    randStateLS.open(mainResumeFolderName+"/"+"randStateLS.csv");
    randStateLS << mt_rand;
    randStateLS.close();


    // saving oldRoundCounter to be zero
    ofstream roundCounterLS;
    roundCounterLS.open(mainResumeFolderName+"/"+"roundCounterLS.csv");
    roundCounterLS << 0;
    roundCounterLS.close();

    // saving the initial time
    ofstream tLS;
    tLS.open(mainResumeFolderName+"/"+"tLS.csv");
    tLS << 0;
    tLSValue=0;
    tLS.close();
    //////////////////////////SAVING IN MAIN RESUME/////////////////////////////
    

    //////////////////////////SAVING IN BACKUP RESUME///////////////////////////
    // saving the initial time
    //ofstream tLSCheck;
    tLSCheck.open(backupResumeFolderName+"/"+"tLSCheck.csv");
    tLSCheck << 0;
    tLSCheck.close();

    // saving LS random generator
    //ofstream randStateLS;
    randStateLS.open(backupResumeFolderName+"/"+"randStateLS.csv");
    randStateLS << mt_rand;
    randStateLS.close();

    // saving oldRoundCounter to be zero
    //ofstream roundCounterLS;
    roundCounterLS.open(backupResumeFolderName+"/"+"roundCounterLS.csv");
    roundCounterLS << 0;
    roundCounterLS.close();

    // saving the initial time
    //ofstream tLS;
    tLS.open(backupResumeFolderName+"/"+"tLS.csv");
    tLS << 0;
    tLS.close();
    //////////////////////////SAVING IN BACKUP RESUME///////////////////////////

    // loadFolderName = mainResumeFolderName;

    writeCounter = 0;

  }
  else
  {
    if(loadSwitch == 1)
    {
      loadFolderName = mainResumeFolderName;
    }
    else if(loadSwitch == 2)
    {
      loadFolderName = backupResumeFolderName;
    }

    // loading oldRoundCounter to be zero
    int roundCounterOld;
    std::ifstream roundCounterLS(loadFolderName+"/"+"roundCounterLS.csv");
    roundCounterLS >> roundCounterOld;
    roundCounterLS.close();
    writeCounter = roundCounterOld+1;

    // LOADING RANDOM GENERATOR STATE LAST SAVED
    std::ifstream randStateLS(loadFolderName+"/"+"randStateLS.csv");
    randStateLS >> mt_rand;
    randStateLS.close();

    loadInt1DVec(sigmaMat, loadFolderName+"/"+"sigmaMatLS.csv");
    loadDbl1DVec(sitesX, loadFolderName+"/"+"sitesX_LS.csv");
    loadDbl1DVec(sitesY, loadFolderName+"/"+"sitesY_LS.csv");

    // reproducing the final last saved areas amd COM coordinates
    areaCalc(sigmaMat, cellArea, cellSiteNum, latticeArea);

    // // calculate comXY based on loaded configuration
    // vector<double> comX_load(NumCells + 1);
    // vector<double> comY_load(NumCells + 1);
    // comCalc(sigmaMat, sitesX, sitesY, cellArea, comX_load, comY_load);

    // // loading thetaLS
    // thetaLoader(theta, loadFolderName+"/"+"thetaLS.csv");

    // e0 = wholeSystemEnergyCalc(sigmaMat, cellArea, Alpha, Lambda, AvgCellArea); // Calculationg final last saved energy
    e0 = wholeSystemEnergyCalc(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea, NumCells);
    E = e0.total;
    

    tLS.open(loadFolderName+"/"+"tLS.csv");
    tLS >> tLSValue;
    tLS.close();
  }
  ////////////////////// LOADING OR INITIATING ////////////////

  int sampleCounter = 0;
  int cellIndex;
  tOld = tLSValue;
  t = tOld + 1;

  long attemptCounter;
  long indAttFrom; // randomly selected positions to be copied FROM
  long indAttInto; // randomly selected positions to be copied INTO
  // long pairIndex, pairIndexTemp, elementInd, elementC;
  long siteIndex, siteIndexTemp, elementInd, elementC;
  int sigmaFrom, sigmaInto; // value in the site FROM and INTO
  int flag;
  int dFNInto;               // delta FOREIGN NEIGHBOR INTO
  vector<int> variation_vec(maxNeighbors); // stores the variation of spin equality
  double dEInter, dEArea , dEAct; // INTERACTION and ELASTIC parts of energy
  int siteC, NNeigh, neighIndC, neighInd, neighIndCFromInto, NNeighNeigh, neighNeighIndC;
  


  int samplingPatternCounter = 0;
  for (int i=1; i<SPMLLength; i++)
  {
    int condition = ((t>=samplingPatternMilestoneList[i-1]) && (t<samplingPatternMilestoneList[i]));
    if(condition)
    {
      samplingPatternCounter=i;
      break;
    }
  }
  long samplingInterval = samplingIntervalList[samplingPatternCounter];
  long tLastSamplingPatternMilestone = 0;
  unsigned long tNextSamplingPatternMilestone = samplingPatternMilestoneList[samplingPatternCounter];

  int timeLoopFlag;
  timeLoopFlag=0;

  struct Energy E_test;
  // double E_test_eps = 0.01*Alpha;
  double E_test_eps = std::min(0.01*Alpha, 0.01*Lambda);
  double dE_eps = 1.0 * (std::min(fabs(Alpha), fabs(Lambda))) * (1e-10);

  ///////// RELATED TO CTIVITY
  // vector<double> comX(NumCells + 1);
  // vector<double> comY(NumCells + 1);
  // comCalc(sigmaMat, sitesX, sitesY, cellArea, comX, comY);
  // double cellAreaFromUpdated, cellAreaIntoUpdated;

  // vector<double> comX_test(NumCells + 1);
  // vector<double> comY_test(NumCells + 1);
  // double com_test_eps = 0.0001;
  // int com_test_flag = 0;

  // int cellC;
  // double thetaUpdateTerm;
  // double eta;
  // double SIG = sqrt(2.0*Dr); // standard dev of theta fluctuation
  // double u1, u2; // for Box-Muller transformation

  // double comX_old_Into, comY_old_Into, area_old_Into;
  // double comX_new_Into, comY_new_Into, area_new_Into;
  // double comX_old_From, comY_old_From, area_old_From;
  // double comX_new_From, comY_new_From, area_new_From;
  // int newX, newY;

  // for (int cellC=1; cellC<=NumCells; cellC++)
  // {
  //   n_x[cellC] = cos(theta[cellC]);
  //   n_y[cellC] = sin(theta[cellC]);
  // }
  // // double theta_eps = 1e-6;
  ///////// RELATED TO CTIVITY



  struct LinkedListElement *LLhead, *LLtemp, *LLtempPrev, *LLtail;
  if (numLinkedList==0)
  {
    cout<<"####################################################\n";
    cout<<"Error! Number of linked lists must be higher than 0!\n";
    cout<<"####################################################\n";
    exit(0);
  }
  vector<LinkedListElement*> LLheadVec(numLinkedList);
  vector<LinkedListElement*> LLtailVec(numLinkedList);
  vector<int> LLNumVec(numLinkedList);
  vector<int> isBorder(NSites);
  vector<int> LLUpdateSwitch(maxNeighbors); // stores the switch of LL update action
  // LinkedListCreation(LLheadVec, LLtailVec, LLNumVec, sigmaMat, cellSiteNum, neighborNum, neighborsList, siteComX, Lx);
  LinkedListSiteCreation(LLheadVec, LLtailVec, LLNumVec, sigmaMat, isBorder, neighborNum, neighborsList, siteComX, Lx);
  // int totalBorderPairs=0;
  int totalBorderSites=0;
  int into_border_switch;
  int LLInd, LLflag;
  for (int i = 0; i < numLinkedList; i++)
  {
    // totalBorderPairs += LLNumVec[i];
    totalBorderSites += LLNumVec[i];
  }
  

  while (tLSValue < maxMCSteps)
  {
    timeLoopFlag=1;
    for (attemptCounter = 0; attemptCounter < NSites; attemptCounter++)
    {
      do
      {
        siteIndexTemp = (mt_rand()) % totalBorderSites;
        // siteIndexTemp = siteIndex;
        LLInd = 0;

        while (siteIndexTemp >= LLNumVec[LLInd])
        {
          siteIndexTemp -= LLNumVec[LLInd];
          LLInd++;
        }
        // elementInd = siteIndexTemp;

        LLtemp = LLheadVec[LLInd];
        // LLtemp = LLhead;
        for (elementC = 0; elementC < siteIndexTemp; elementC++)
        {
          LLtemp = LLtemp->next;
        }
        // indAttFrom = LLtemp->site;
        
        // neighIndC = (mt_rand()) % neighborNum[indAttFrom];
        // indAttInto = neighborsList[indAttFrom][neighIndC];
        
        indAttInto = LLtemp->site;
        
        neighIndC = (mt_rand()) % neighborNum[indAttInto];
        indAttFrom = neighborsList[indAttInto][neighIndC];

      } while ((sigmaMat[indAttInto]==sigmaMat[indAttFrom]) || (cellSiteNum[sigmaMat[indAttInto]]<=1));
      
      sigmaFrom = sigmaMat[indAttFrom];
      sigmaInto = sigmaMat[indAttInto];
      neighIndCFromInto = neighIndC; //BE CAREFUL: this is the index of 'indAttFrom' among the list of the neighbors of 'indAttInto'.

      
      dFNInto = 0;
      NNeigh = neighborNum[indAttInto];
      for (neighIndC = 0; neighIndC < NNeigh; neighIndC++)
      {
        neighInd = neighborsList[indAttInto][neighIndC];
        // dFNInto = dFNInto + (int)(sigmaFrom != sigmaMat[neighInd]) - (int)(sigmaInto != sigmaMat[neighInd]);
        // variation_vec[neighIndC] = (int)(sigmaFrom != sigmaMat[neighInd]) - (int)(sigmaInto != sigmaMat[neighInd]);
        variation_vec[neighIndC] = (sigmaFrom != sigmaMat[neighInd]) - (sigmaInto != sigmaMat[neighInd]);
        dFNInto = dFNInto + variation_vec[neighIndC];
      }
      
      // comX_old_Into = comX[sigmaInto];
      // comY_old_Into = comY[sigmaInto];
      // area_old_Into = cellArea[sigmaInto];

      // comX_old_From = comX[sigmaFrom];
      // comY_old_From = comY[sigmaFrom];
      // area_old_From = cellArea[sigmaFrom];

      // comX_new_Into = (1. / (area_old_Into-1.)) * (comX_old_Into * area_old_Into - sitesX[rowAttInto][colAttInto]);
      // comY_new_Into = (1. / (area_old_Into-1.)) * (comY_old_Into * area_old_Into - sitesY[rowAttInto][colAttInto]);
      
      // newX = sitesX[rowAttFrom][colAttFrom] + intoSelectorMat[randIntInto][0]; // 0 is for x axis
      // newY = sitesY[rowAttFrom][colAttFrom] + intoSelectorMat[randIntInto][1]; // 1 is for y axis

      // comX_new_From = (1. / (area_old_From+1.)) * (comX_old_From * area_old_From + newX);
      // comY_new_From = (1. / (area_old_From+1.)) * (comY_old_From * area_old_From + newY);

      dEInter = Alpha * dFNInto;
      if (sigmaFrom*sigmaInto > 0)
        {dEArea = 2.0 * Lambda * latticeArea[indAttInto] * (latticeArea[indAttInto] + cellArea[sigmaFrom] - cellArea[sigmaInto]);}
      else if (sigmaInto==0)
        {dEArea =  Lambda * latticeArea[indAttInto] * (latticeArea[indAttInto] + 2 * ( cellArea[sigmaFrom] - AvgCellArea));}
      else if (sigmaFrom==0)
        {dEArea =  Lambda * latticeArea[indAttInto] * (latticeArea[indAttInto] - 2 * ( cellArea[sigmaInto] - AvgCellArea));}
      // dEAct = ;

      // dE = dE_inter + dE_area + dE_act;
      dE = dEInter + dEArea;

      if (  (dE < dE_eps) || ((((long double)(mt_rand())-MT_MIN)/((long double)MT_MAX-MT_MIN)) < 
                  exp(-dE/(Tem*Kb)) )   )
      {

        E += dE;
        

        // Area update
        cellArea[sigmaFrom] += latticeArea[indAttInto];
        cellArea[sigmaInto] -= latticeArea[indAttInto];
        cellSiteNum[sigmaFrom]++;
        cellSiteNum[sigmaInto]--;
        
        //sigmaMat update
        sigmaMat[indAttInto] = sigmaFrom;

        // E_test = wholeSystemEnergyCalc(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea);

        // comXY update INTO
        // comX[sigmaInto] = comX_new_Into;
        // comY[sigmaInto] = comY_new_Into;

        // siteXY update
        // neighIndC = 0;
        // while (indAttInto != neighborsList[indAttFrom][neighIndC])
        // {
        //   neighIndC++;
        // }
        // sitesX[indAttInto] = sitesX[indAttFrom] + deltaComX[indAttFrom][neighIndC];
        // sitesY[indAttInto] = sitesY[indAttFrom] + deltaComY[indAttFrom][neighIndC];

        // sitesX[indAttInto] = sitesX[indAttFrom] + deltaComX[indAttFrom][neighIndCFromInto];
        // sitesY[indAttInto] = sitesY[indAttFrom] + deltaComY[indAttFrom][neighIndCFromInto];

        sitesX[indAttInto] = sitesX[indAttFrom] - deltaComX[indAttInto][neighIndCFromInto];
        sitesY[indAttInto] = sitesY[indAttFrom] - deltaComY[indAttInto][neighIndCFromInto];

        // comXY update FROM
        // comX[sigmaFrom] = comX_new_From;
        // comY[sigmaFrom] = comY_new_From;


        // SWITCHES OF UPDATING LINKED LISTS (SITES)
        into_border_switch = 0;
        for (neighIndC = 0; neighIndC < maxNeighbors; neighIndC++)
        {
          LLUpdateSwitch[neighIndC] = 0;
          // 0 menas unknown (must be checked)
          // 1 means: it was not border before. It is border now. (must be added to LL)
          //-1 means: it was border before. It is not now. (must be removed from LL)
          // 2 means: Do nothing. Its status is exactly like before. (It was border, and it is still border)
          // it is not possible that a site was not border, and is not a border yet, because one of its neighbors has changed.
        }
        
        for (neighIndC = 0; neighIndC < NNeigh; neighIndC++)
        {
          neighInd = neighborsList[indAttInto][neighIndC];

          if (variation_vec[neighIndC]==0) // (1-1) They were inequal. They were both inside LL. They are still inequal. They must remain in LL.
          { // So, do not touch anything. Just turn `into_border_flag' on.
            into_border_switch = 2;
            LLUpdateSwitch[neighIndC] = 2;
          }
          else if (variation_vec[neighIndC]==1) //(1-0) They were equal. Must be checked if they were in LL. They are now inequal. They must be in LL now.
          {
            into_border_switch = 2; // indAttInto was definitely border. It is still border. So, do nothing about it.

            if(isBorder[neighInd]==0)
            {
              // add [neighInd] to proper LL;
              LLUpdateSwitch[neighIndC] = 1;
              // isBorder[neighInd]=1;
            }
          }
          else if(variation_vec[neighIndC]==-1) //(0-1) They were inequal. They were both in LL. They are now equal.
          {
            NNeighNeigh = neighborNum[neighInd];
            for (neighNeighIndC = 0; neighNeighIndC < NNeighNeigh; neighNeighIndC++)
            {
              if (sigmaMat[neighInd] != sigmaMat[neighborsList[neighInd][neighNeighIndC]])
              {
                LLUpdateSwitch[neighIndC] = 2;
                break;
              }
            }
            if (LLUpdateSwitch[neighIndC] == 0)
            {
              LLUpdateSwitch[neighIndC] = -1;
            }
          }
        }

        if (into_border_switch==0)
        {
          into_border_switch = -1;
        }
        // SWITCHES OF UPDATING LINKED LISTS (SITES)


        // UPDATE LINKED LISTS (SITES)
        if (into_border_switch==-1) //delete the indAttInto from the LL
        {
          LLInd = (int)(floor(1.0*siteComX[indAttInto]/(1.0*Lx/numLinkedList)));
          LLhead = LLheadVec[LLInd];
          LLtail = LLtailVec[LLInd];
          if(LLhead->site != indAttInto)
          {
            LLtempPrev = LLhead;
            LLtemp     = LLhead->next;
            while (LLtemp->site != indAttInto)
            {
              LLtempPrev = LLtempPrev -> next;
              LLtemp     = LLtemp -> next;
            }
            LLtempPrev -> next = LLtemp -> next;
            if (LLtemp == LLtail)
              {
                LLtail = LLtempPrev;
                LLtailVec[LLInd] = LLtail;
              }
            delete LLtemp;
          }
          else //LLhead->site == indAttInto
          {
            LLtemp = LLhead->next;  
            delete LLhead;
            LLhead = LLtemp;
            LLheadVec[LLInd] = LLhead;
            
            if (LLNumVec[LLInd]==1) //the list is going to become empty
            {
              LLtailVec[LLInd] = NULL;
            }
          }
          isBorder[indAttInto]=0;
          LLNumVec[LLInd]--;
        }
        
        //neigbors of Into
        for (neighIndC = 0; neighIndC < NNeigh; neighIndC++)
        {
          switch (LLUpdateSwitch[neighIndC])
          {
            case (-1): //delete the neighInd from the LL
              neighInd = neighborsList[indAttInto][neighIndC];
              LLInd = (int)(floor(1.0*siteComX[neighInd]/(1.0*Lx/numLinkedList)));
              LLhead = LLheadVec[LLInd];
              LLtail = LLtailVec[LLInd];
              if(LLhead->site != neighInd)
              {
                LLtempPrev = LLhead;
                LLtemp     = LLhead->next;
                while (LLtemp->site != neighInd)
                {
                  LLtempPrev = LLtempPrev -> next;
                  LLtemp     = LLtemp -> next;
                }
                LLtempPrev -> next = LLtemp -> next;
                if (LLtemp == LLtail)
                  {
                    LLtail = LLtempPrev;
                    LLtailVec[LLInd] = LLtail;
                  }
                delete LLtemp;
              }else //LLhead->site == neighInd
              {
                LLtemp = LLhead->next;  
                delete LLhead;
                LLhead = LLtemp;
                LLheadVec[LLInd] = LLhead;
                if (LLNumVec[LLInd]==1) //the list is going to become empty
                {
                  LLtailVec[LLInd] = NULL;
                }
              }
              isBorder[neighInd]=0;
              LLNumVec[LLInd]--;
              break;
              

            case 1: //add the neighInd to the LL
              neighInd = neighborsList[indAttInto][neighIndC];
              LLInd = (int)(floor(1.0*siteComX[neighInd]/(1.0*Lx/numLinkedList)));
              LLhead = LLheadVec[LLInd];
              LLtail = LLtailVec[LLInd];
              if (LLtailVec[LLInd]) // the LL is NOT initially empty.
              {
                // LLhead = LLheadVec[LLInd];
                LLtail = LLtailVec[LLInd];
                LLtail -> next = new LinkedListElement;
                (LLtail -> next)->site  = neighInd;
                (LLtail -> next)->next = NULL;
                LLtail = LLtail->next;
                LLtailVec[LLInd] = LLtail;
              }
              else //the list is initially empty (head and tail, both are NULL)
              {
                LLhead = new LinkedListElement;
                LLheadVec[LLInd] = LLhead;
                LLtail = LLhead;
                LLtail->site  = neighInd;
                LLtail->next = NULL;
                LLtailVec[LLInd] = LLtail;
              }
              isBorder[neighInd]=1;
              LLNumVec[LLInd]++;
              break;

            case 2: //do nothing
              break;
          }
        }

        totalBorderSites=0;
        for (int i = 0; i < numLinkedList; i++)
        {
          totalBorderSites += LLNumVec[i];
        }
        // UPDATE LINKED LISTS (SITES)

      }
    }

    // //////////////////// UPDATING THETA //////////////////
    // for (int cellC=1; cellC<=NumCells; cellC++)
    // {

    //   randState = MT_MIN;
    //   while(randState == MT_MIN)
    //   {
    //     randState = mt_rand();
    //   }
    //   // randState = mt_rand();
    //   u1 = ((long double)(randState) - MT_MIN) / ((long double)MT_MAX - MT_MIN);

    //   randState = MT_MIN;
    //   while(randState == MT_MIN)
    //   {
    //     randState = mt_rand();
    //   }
    //   // randState = mt_rand();
    //   u2 = ((long double)(randState) - MT_MIN) / ((long double)MT_MAX - MT_MIN);

    //   eta = SIG * sqrt(-2.0*log(u1)) * cos(2.0*M_PI*u2);
    //   thetaUpdateTerm = eta * 1.00; // dt = 1.00

    //   E -= ((-p) * (n_x[cellC] * comX[cellC] + n_y[cellC] * comY[cellC]) );

    //   theta[cellC] = theta[cellC] + thetaUpdateTerm;

    //   while (theta[cellC] > 2.0*M_PI)
    //   {
    //     theta[cellC] = theta[cellC] - 2.0*M_PI;
    //   }
    //   while (theta[cellC] < 0.0)
    //   { 
    //     theta[cellC] = theta[cellC] + 2.0*M_PI;
    //   }

    //   n_x[cellC] = cos(theta[cellC]);
    //   n_y[cellC] = sin(theta[cellC]);

    //   E += ((-p) * (n_x[cellC] * comX[cellC] + n_y[cellC] * comY[cellC]) );
    // }
    // //////////////////// UPDATING THETA //////////////////


  /////////////////////// SAMPLING OPERATION ////////////////////////
    if ((t-tLastSamplingPatternMilestone) % samplingInterval == 0)
    {
      tSamples[sampleCounter] = t;
      eSamples[sampleCounter] = E;
      for (cellIndex = 0; cellIndex <= NumCells; cellIndex++)
      {
        areaSamples[cellIndex][sampleCounter] = cellArea[cellIndex];
      }
      for (siteC = 0; siteC < NSites; siteC++)
      {
        sigmaSamples[siteC][sampleCounter] = sigmaMat[siteC];
        xSamples[siteC][sampleCounter] = floor(sitesX[siteC]/Lx);
        ySamples[siteC][sampleCounter] = floor(sitesY[siteC]/Ly);
      }
      sampleCounter++;

      /////////////////////// WRITING OPERATION ////////////////////////
      if (sampleCounter == samplesPerWrite)
      {
        
        // cout << "\033[2J\033[1;1H";
        // cout << "Saving! Do not terminate!"<<endl;


        // writing the center of mass of the cells
        saveIntMatCSV(sigmaSamples, dataFolderName + "/" + "sigmaSamples_" + to_string(writeCounter) + ".csv");
        // saveDoubleMatCSV(xSamples, dataFolderName + "/" + "xSamples_" + to_string(writeCounter) + ".csv");
        // saveDoubleMatCSV(ySamples, dataFolderName + "/" + "ySamples_" + to_string(writeCounter) + ".csv");
        saveIntMatCSV(xSamples, dataFolderName + "/" + "xSamples_" + to_string(writeCounter) + ".csv");
        saveIntMatCSV(ySamples, dataFolderName + "/" + "ySamples_" + to_string(writeCounter) + ".csv");

        // writing the area and perimiter of the cells
        // saveDoubleMatCSV(areaSamples, dataFolderName + "/" + "areaSamples_" + to_string(writeCounter) + ".csv");
        // saveSampleAP(sampleP, dataFolderName+"/"+"sampleP_"+to_string(writeCounter)+".csv");

        // writing the sampling times and energies
        saveSampleT(tSamples, dataFolderName + "/" + "tSamples_" + to_string(writeCounter) + ".csv");
        saveSampleE(eSamples, dataFolderName + "/" + "eSamples_" + to_string(writeCounter) + ".csv");


        ///////////// NUMBERED DATA ZIPPING ///////////////////
        if ((writeCounter+1)%writePerZip ==0)
        {
          // string argv1 = 'data_'+to_string((int)(writeCounter/writePerZip))+'.zip';
          std::string argv1 = "data_" + std::to_string(static_cast<int>(writeCounter / writePerZip)+1) + ".zip";
          std::string argv2 = dataFolderName;
          std::string command = "python3 dataZipperNum.py " + argv1 + " " + argv2;
          system(command.c_str());
        }
        ///////////// NUMBERED DATA ZIPPING ///////////////////


        ///////////////////////// LS SAVING ///////////////////////
        for(int resumeFolderCounter=1; resumeFolderCounter<=2; resumeFolderCounter++)
        {
          if (resumeFolderCounter==1)
          {
            loadFolderName = mainResumeFolderName;
          }
          else if (resumeFolderCounter==2)
          {
            loadFolderName = backupResumeFolderName;
          }

          // writing the final time
          ofstream tLSCheck;
          tLSCheck.open(loadFolderName+"/"+"tLSCheck.csv");
          tLSCheck << t;
          tLSCheck.close();

          // writing the final configuration of the sites
          saveInt1DVec(sigmaMat, loadFolderName+"/"+"sigmaMatLS.csv");
          saveDbl1DVec(sitesX, loadFolderName+"/"+"sitesX_LS.csv");
          saveDbl1DVec(sitesY, loadFolderName+"/"+"sitesY_LS.csv");

          // saveThetaCSV(theta, loadFolderName+"/"+"thetaLS.csv");

          // writing the final state of random generator
          std::ofstream randStateLS(loadFolderName+"/"+"randStateLS.csv");
          randStateLS << mt_rand;
          randStateLS.close();

          // writing round counter
          ofstream roundCounterLS;
          roundCounterLS.open(loadFolderName+"/"+"roundCounterLS.csv");
          roundCounterLS << writeCounter;
          roundCounterLS.close();

          // writing the final time
          ofstream tLS;
          tLS.open(loadFolderName+"/"+"tLS.csv");
          tLSValue = t;
          tLS << tLSValue;
          tLS.close();

          // cout << "\033[2J\033[1;1H";
        }
        ///////////////////////// LS SAVING ///////////////////////


        ///////////// TESTING (ENERGY & COM) ///////////////////
        // E_test = wholeSystemEnergyCalc(sigmaMat, cellArea, Alpha, Lambda, AvgCellArea); // Calculationg final last saved energy
        // E_test = wholeSystemEnergyCalc_actPlus(sigmaMat, cellArea, Alpha, Lambda, AvgCellArea, p, theta, comX, comY);
        E_test = wholeSystemEnergyCalc(sigmaMat, neighborsList, neighborNum, cellArea, Alpha, Lambda, AvgCellArea, NumCells);
        

        if (fabs(E_test.total-E)>E_test_eps)
        {
          ofstream EnergyError;
          EnergyError.open("EnergyError.csv");
          EnergyError << "E: ";
          EnergyError << E;
          EnergyError << endl;
          EnergyError << "E_test: ";
          EnergyError << E_test.total;
          EnergyError.close();

          cout<<"Energy calculation inconsistency!"<<endl;
          cout<<"Look at the error file."<<endl;
          cout<<"Program ended!"<<endl;
          exit(0);
        } else
        {
          // cout<<"Energy successfully checked!"<<endl;
        }
        

        areaCalc(sigmaMat, cellAreaTest, cellSiteNumTest, latticeArea);
        int areaFlag = 0;
        for (int cellCTest = 0; cellCTest <= NumCells; cellCTest++)
        {
          if ( (fabs(cellAreaTest[cellCTest]-cellArea[cellCTest])>(1e-6)) || (fabs(cellSiteNumTest[cellCTest]-cellSiteNum[cellCTest])>(1e-6)))
          {
            areaFlag = 1;

            ofstream AreaError;
            AreaError.open("AreaError.csv");
            AreaError << "Area["<<cellCTest<<"]= "<<cellArea[cellCTest];
            AreaError << endl;
            AreaError << "AreaTest["<<cellCTest<<"]= "<<cellAreaTest[cellCTest];
            AreaError.close();

            cout<<"Area calculation inconsistency!"<<endl;
            cout<<"Look at the error file."<<endl;
            cout<<"Program ended!"<<endl;
            exit(0);
          }
        }
        if (areaFlag==0)
        {
          // cout<<"Area successfully checked!"<<endl;
        }



        // // comXY
        // comCalc(sigmaMat, sitesX, sitesY, cellArea, comX_test, comY_test);
        // for (int dum_counter = 0; dum_counter <= NumCells; dum_counter++)
        // {
        //   if( (fabs(comX[dum_counter]-comX_test[dum_counter])>com_test_eps) || (fabs(comY[dum_counter]-comY_test[dum_counter])>com_test_eps))
        //   {
        //     ofstream comError;
        //     comError.open("comError.csv");
        //     comError << "cell Index: ";
        //     comError << dum_counter;
        //     comError << endl;
        //     comError << "comX (update): ";
        //     comError << comX[dum_counter];
        //     comError << endl;
        //     comError << "comX_test (function): ";
        //     comError << comX_test[dum_counter];
        //     comError << endl;
        //     comError << "comY (update): ";
        //     comError << comY[dum_counter];
        //     comError << endl;
        //     comError << "comY_test (function): ";
        //     comError << comY_test[dum_counter];
        //     comError << endl;
        //     comError.close();

        //     cout<<"comXY calculation inconsistency!"<<endl;
        //     cout<<"Look at the error file."<<endl;
        //     cout<<"Program ended!"<<endl;
        //     com_test_flag = 1;
        //     exit(0);
        //   }
        // }
        // if(com_test_flag==0)
        // {
        //   cout<<"comXY successfully checked!"<<endl;
        // }
        ///////////// TESTING (ENERGY & COM) ///////////////////

        writeCounter++;
        sampleCounter = 0;
      }
      /////////////////////// WRITING OPERATION ////////////////////////


      /////////////////////// CHANGING SAMPLING PATTERN ////////////////////////
      if(t == tNextSamplingPatternMilestone)
      {
        samplingPatternCounter++;
        samplingInterval = samplingIntervalList[samplingPatternCounter];
        tLastSamplingPatternMilestone = samplingPatternMilestoneList[samplingPatternCounter-1];
        tNextSamplingPatternMilestone = samplingPatternMilestoneList[samplingPatternCounter];
      }
      /////////////////////// CHANGING SAMPLING PATTERN ////////////////////////
    }
    /////////////////////// SAMPLING OPERATION ////////////////////////

    if (!(t % printingTimeInterval))
    {
      // cout << "\033[2J\033[1;1H";
      cout << "t : " << t << "/" << maxMCSteps << endl;
    }

    t++;
  }
  t--;

  
  if(timeLoopFlag)
  {

    auto stop_time = high_resolution_clock::now();
    // Calculate the duration of the code execution in minutes
    auto duration = duration_cast<seconds>(stop_time - start_time);
    // Output the duration in minutes
    cout << "Time taken by function: " << duration.count() << " seconds" << endl;
    std::ofstream Time_outfile("time.dat");
    if (Time_outfile.is_open())
    {
      Time_outfile << duration.count();
      Time_outfile << std::endl;
      Time_outfile.close();
    }
    else
    {
      std::cerr << "Could not open file for writing" << std::endl;
    }
  }

  system("python3 dataZipper.py");

  // system("./pp_irr_v1");

  // system("python3 ppDataZipper.py");

  cout << "Program done!" << endl<< endl;

  // return 0;
}
//
void energyEpsilonsFinder(double &E_test_eps, double &dE_eps,\
                          const vector<vector<double>>& J_int,\
                          const vector<vector<double>>& J_ext,\
                          const double Lambda, const double p)
{
  double min_J;
  int NCompart = J_int.size();

  min_J = J_ext[0][0]; // I chose this as initial value, because this is always positive.

  for (int i = 0; i < NCompart; i++) // sweep on J_ext
  {
    for (int j = 0; j < NCompart; j++)
    {
      if (J_ext[i][j] < min_J && J_ext[i][j] > (1e-8))
      {
        min_J = J_ext[i][j];
      }
    }
  }

  for (int i = 0; i < NCompart; i++) // sweep on J_int
  {
    for (int j = 0; j < NCompart; j++)
    {
      if (J_int[i][j] < min_J && J_int[i][j] > (1e-10))
      {
        min_J = J_int[i][j];
      }
    }
  }
  

  if (fabs(p) > (1e-10)) // p is practically zero, so, we do not consider it
  {
    E_test_eps = std::min(std::min(0.01*min_J, 0.01*Lambda), 0.01*p );
    dE_eps = 1.0 * std::min( (std::min(fabs(min_J), fabs(Lambda))), fabs(p)) * (1e-8); 
  }else // p is non-zero, so, we have to consider it
  {
    E_test_eps = std::min(0.01*min_J, 0.01*Lambda);
    dE_eps = 1.0 * (std::min(fabs(min_J), fabs(Lambda))) * (1e-8);
  }
}
//
struct Energy wholeSystemEnergyCalc(const vector<int>& sigmaMat,\ 
                                    const vector<vector<int>>& neighborsList, const vector<int>& neighborNum,\ 
                                    const vector<double>& cellArea,\
                                    const double Alpha, const double Lambda, const double AvgCellArea, const int NumCells)
{
  int NSites = sigmaMat.size();
  // int NumCells = cellArea.size() - 1;

  struct Energy wholeSystemEnergy;
  
  int siteC, NNeigh, neighIndC, neighInd;

  wholeSystemEnergy.inter = 0.;
  for (siteC = 0; siteC < NSites; siteC++)
  {
    NNeigh = neighborNum[siteC];
    for (neighIndC = 0; neighIndC < NNeigh; neighIndC++)
    {
      neighInd = neighborsList[siteC][neighIndC];
      if (sigmaMat[siteC] != sigmaMat[neighInd])
      {
        wholeSystemEnergy.inter += 0.5* Alpha;
      }
    }
  }
  

  wholeSystemEnergy.area = 0.;
  for (int cell_c = 1; cell_c <= NumCells; cell_c++)
  {
    wholeSystemEnergy.area += (1.0 * Lambda * pow((cellArea[cell_c] - AvgCellArea), 2.));
  }

  wholeSystemEnergy.total = wholeSystemEnergy.inter + wholeSystemEnergy.area;

  return wholeSystemEnergy;
}
//
struct Energy energyCalcCompart(const vector<int>& sigmaMat, const vector<int>& compartMat,\
                                const vector<vector<double>>& J_int, const vector<vector<double>>& J_ext,\
                                const vector<vector<int>>& neighborsList, const vector<int>& neighborNum,\ 
                                const vector<vector<double>>& edges, const double avgEdge,\ 
                                const vector<vector<double>>& compartArea,\
                                const double Lambda, const vector<double>& avgArea, const int NumCells, \
                                const double p, const vector<double>& theta, \
                                const vector<double>& comX, const vector<double>& comY)
{

  int NSites = sigmaMat.size();
  int NCompart = avgArea.size();
  // int NumCells = cellArea.size() - 1;
  struct Energy wholeSystemEnergy;
  
  int siteC, NNeigh, neighIndC, neighInd;
  int cellInd_1, compartInd_1;
  int cellInd_2, compartInd_2;

  //inter
  wholeSystemEnergy.inter = 0.;
  for (siteC = 0; siteC < NSites; siteC++)
  {
    cellInd_1 = sigmaMat[siteC];
    compartInd_1 = compartMat[siteC];

    NNeigh = neighborNum[siteC];
    for (neighIndC = 0; neighIndC < NNeigh; neighIndC++)
    {
      neighInd = neighborsList[siteC][neighIndC];
      
      cellInd_2 = sigmaMat[neighInd];
      compartInd_2 = compartMat[neighInd];

      if (cellInd_1==cellInd_2) // of the same cell (J_int)
      {
        wholeSystemEnergy.inter += 0.5 * J_int[compartInd_1][compartInd_2] * (edges[siteC][neighIndC] / avgEdge);
      }
      else // of different cells (J_ext)
      {
        wholeSystemEnergy.inter += 0.5 * J_ext[compartInd_1][compartInd_2] * (edges[siteC][neighIndC] / avgEdge);
      }
    }
  }
  //inter
  
  //area
  wholeSystemEnergy.area = 0.;
  for (int cell_c = 1; cell_c <= NumCells; cell_c++)
  {
    // double total_A=0;
    for (int compartC = 0; compartC < NCompart; compartC++)
    {
      wholeSystemEnergy.area += (1.0 * Lambda * pow((compartArea[cell_c][compartC] - avgArea[compartC]), 2.)); 
      // total_A += compartArea[cell_c][compartC];
    }
    // wholeSystemEnergy.area += (1.0 * Lambda * pow((total_A - 40.0), 2.)); 
  }
  //area

  //act
  wholeSystemEnergy.act = 0.;
  for (int cell_c = 1; cell_c <= NumCells; cell_c++)
  {
    wholeSystemEnergy.act +=  ((-p) * (cos(theta[cell_c])) * comX[cell_c] );
    wholeSystemEnergy.act +=  ((-p) * (sin(theta[cell_c])) * comY[cell_c] );
  }
  //act

  wholeSystemEnergy.total = wholeSystemEnergy.inter + wholeSystemEnergy.area + wholeSystemEnergy.act;

  return wholeSystemEnergy;

}
//
struct Energy wholeSystemEnergyCalc_actPlus_W(const vector<int>& sigmaMat,\ 
                                    const vector<vector<int>>& neighborsList, const vector<int>& neighborNum,\ 
                                    const vector<vector<double>>& edges, const double avgEdge,\ 
                                    const vector<double>& cellArea,\
                                    const double Alpha, const double Lambda, const double AvgCellArea, const int NumCells, \
                                    const double p, const vector<double>& theta, \
                                    const vector<double>& comX, const vector<double>& comY)
{
  int NSites = sigmaMat.size();
  // int NumCells = cellArea.size() - 1;

  struct Energy wholeSystemEnergy;
  
  int siteC, NNeigh, neighIndC, neighInd;

  //inter
  wholeSystemEnergy.inter = 0.;
  for (siteC = 0; siteC < NSites; siteC++)
  {
    NNeigh = neighborNum[siteC];
    for (neighIndC = 0; neighIndC < NNeigh; neighIndC++)
    {
      neighInd = neighborsList[siteC][neighIndC];
      if (sigmaMat[siteC] != sigmaMat[neighInd])
      {
        wholeSystemEnergy.inter += 0.5 * Alpha * (edges[siteC][neighIndC] / avgEdge);
      }
    }
  }
  //inter
  
  //area
  wholeSystemEnergy.area = 0.;
  for (int cell_c = 1; cell_c <= NumCells; cell_c++)
  {
    wholeSystemEnergy.area += (1.0 * Lambda * pow((cellArea[cell_c] - AvgCellArea), 2.));
  }
  //area

  //act
  wholeSystemEnergy.act = 0.;
  for (int cell_c = 1; cell_c <= NumCells; cell_c++)
  {
    wholeSystemEnergy.act +=  ((-p) * (cos(theta[cell_c])) * comX[cell_c] );
    wholeSystemEnergy.act +=  ((-p) * (sin(theta[cell_c])) * comY[cell_c] );
  }
  //act

  wholeSystemEnergy.total = wholeSystemEnergy.inter + wholeSystemEnergy.area + wholeSystemEnergy.act;

  return wholeSystemEnergy;
}
//
struct Energy wholeSystemEnergyCalc_actPlus(const vector<int>& sigmaMat,\ 
                                    const vector<vector<int>>& neighborsList, const vector<int>& neighborNum,\ 
                                    const vector<double>& cellArea,\
                                    const double Alpha, const double Lambda, const double AvgCellArea, const int NumCells, \
                                    const double p, const vector<double>& theta, \
                                    const vector<double>& comX, const vector<double>& comY)
{
  int NSites = sigmaMat.size();
  // int NumCells = cellArea.size() - 1;

  struct Energy wholeSystemEnergy;
  
  int siteC, NNeigh, neighIndC, neighInd;

  //inter
  wholeSystemEnergy.inter = 0.;
  for (siteC = 0; siteC < NSites; siteC++)
  {
    NNeigh = neighborNum[siteC];
    for (neighIndC = 0; neighIndC < NNeigh; neighIndC++)
    {
      neighInd = neighborsList[siteC][neighIndC];
      if (sigmaMat[siteC] != sigmaMat[neighInd])
      {
        wholeSystemEnergy.inter += 0.5* Alpha;
      }
    }
  }
  //inter
  
  //area
  wholeSystemEnergy.area = 0.;
  for (int cell_c = 1; cell_c <= NumCells; cell_c++)
  {
    wholeSystemEnergy.area += (1.0 * Lambda * pow((cellArea[cell_c] - AvgCellArea), 2.));
  }
  //area

  //act
  wholeSystemEnergy.act = 0.;
  for (int cell_c = 1; cell_c <= NumCells; cell_c++)
  {
    wholeSystemEnergy.act +=  ((-p) * (cos(theta[cell_c])) * comX[cell_c] );
    wholeSystemEnergy.act +=  ((-p) * (sin(theta[cell_c])) * comY[cell_c] );
  }
  //act

  wholeSystemEnergy.total = wholeSystemEnergy.inter + wholeSystemEnergy.area + wholeSystemEnergy.act;

  return wholeSystemEnergy;
}
//
void LinkedListCreation(vector<LinkedListElementPair*>& LLheadVec, vector<LinkedListElementPair*>& LLtailVec, vector<int>& LLNumVec,\
                        const vector<int>& sigmaMat, const vector<int>& cellSiteNum, \
                        const vector<int>& neighborNum, const vector<vector<int>>& neighborsList,\
                        const vector<double>& siteComX, const double Lx)
{

  // vector<std::pair<int, int>> test_vec;
  // std::pair<int, int> test_element;
  
  int NSites = sigmaMat.size();
  int numLinkedList = LLheadVec.size();

  vector<LinkedListElementPair*> LLTempVec(numLinkedList);
  LinkedListElementPair *LLTemp;
  for (int i = 0; i < numLinkedList; i++)
  {
    LLNumVec[i] = 0;
    LLheadVec[i] = NULL;
    LLTempVec[i] = NULL;
  }
  

  int siteC, NNeigh, neighIndC, neighInd;
  int LLIndex;
  for (siteC = 0; siteC < NSites; siteC++)
  {
    LLIndex = (int)(floor(1.0*siteComX[siteC]/(1.0*Lx/numLinkedList)));

    NNeigh = neighborNum[siteC];
    for (neighIndC = 0; neighIndC < NNeigh; neighIndC++)
    {
      neighInd = neighborsList[siteC][neighIndC];

      // if ((sigmaMat[siteC]!=sigmaMat[neighInd]) && (cellSiteNum[neighInd]>1))
      if ( sigmaMat[siteC]!=sigmaMat[neighInd] )
      {
        
        if (LLNumVec[LLIndex]>0)
        {
          LLTemp = LLTempVec[LLIndex];
          LLTemp -> next = new LinkedListElementPair;
          LLTemp = LLTemp -> next;
          LLTemp -> sitesPair.first  = siteC;
          LLTemp -> sitesPair.second = neighInd;
          LLTemp -> next = NULL;
          LLTempVec[LLIndex] = LLTemp;
          LLtailVec[LLIndex] = LLTemp;
        }
        else if (LLNumVec[LLIndex]==0)
        {
          LLTemp = new LinkedListElementPair;
          LLTemp -> sitesPair.first  = siteC;
          LLTemp -> sitesPair.second = neighInd;
          LLTemp -> next = NULL;
          LLTempVec[LLIndex] = LLTemp;
          LLtailVec[LLIndex] = LLTemp;
          LLheadVec[LLIndex] = LLTemp;
        }
        
        // test_element.first  = siteC;
        // test_element.second = neighInd;
        // test_vec.push_back(test_element);
        LLNumVec[LLIndex]++;
      }
    }
  }

  // LinkedListElement* Ptr = LLheadVec[3];
  // while (Ptr)
  // {
  //   int first  = Ptr->sitesPair.first;
  //   int second = Ptr->sitesPair.second;
  //   LinkedListElement* next = Ptr->next;
  //   cout<<"----------------------------"<<endl;
  //   cout<<"address: "<<Ptr<<endl;
  //   cout<<"first:   "<<first<<endl;
  //   cout<<"second:  "<<second<<endl;
  //   cout<<"next:    "<<next<<endl;
  //   cout<<"----------------------------"<<endl;
  //   Ptr = Ptr->next;
  // }

  // for (int i = 0; i < numLinkedList; i++)
  // {
  //   if (LLheadVec[i] == NULL)
  //   {
  //     LLheadVec[i] = new LinkedListElement;
  //     LLTempVec[i] = LLheadVec[i];
  //   }
  // }

  
  // int r;
  // r=5;
}
//
void readCompartData(
    const std::string& filename,
    std::vector<double>& avgAreaFrac,
    std::vector<std::vector<double>>& J_int,
    std::vector<std::vector<double>>& J_ext
) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Cannot open " + filename);
    }

    std::string line;

    // 1) Read avgAreaFrac until "#"
    avgAreaFrac.clear();
    std::getline(in, line); // skip the name of matrix
    while (std::getline(in, line)) {
        if (line == "#") break;
        if (line.empty()) continue;
        avgAreaFrac.push_back(std::stod(line));
    }
    int N = int(avgAreaFrac.size());
    if (N == 0) {
        throw std::runtime_error("No data in first block");
    }

    // prepare matrices
    J_int.assign(N, std::vector<double>(N, 0.0));
    J_ext.assign(N, std::vector<double>(N, 0.0));

    // 2) Read J_int (lower‐triangle ragged form)
    // skip blank or "#"
    while (std::getline(in, line)) {
        if (line.empty() || line == "#") continue;
        else break;
    }
    std::getline(in, line); // skip the name of matrix
    // line now has row 0
    {
        std::istringstream ss(line);
        double v0;
        ss >> v0;
        J_int[0][0] = v0;
    }
    for (int i = 1; i < N; ++i) {
        std::getline(in, line);
        std::istringstream ss(line);
        for (int j = 0; j <= i; ++j) {
            double v;
            ss >> v;
            J_int[i][j] = v;
            J_int[j][i] = v;
        }
    }

    // 3) Read J_ext
    // skip blank or "#"
    while (std::getline(in, line)) {
        if (line.empty() || line == "#") continue;
        else break;
        
    }
    std::getline(in, line); // skip the name of matrix
    // line now has row 0
    {
        std::istringstream ss(line);
        double v0;
        ss >> v0;
        J_ext[0][0] = v0;
    }
    for (int i = 1; i < N; ++i) {
        std::getline(in, line);
        std::istringstream ss(line);
        for (int j = 0; j <= i; ++j) {
            double v;
            ss >> v;
            J_ext[i][j] = v;
            J_ext[j][i] = v;
        }
    }
}
//
void initial_LS_saver(const string FolderName, const std::mt19937 &mt_rand)
{
  // saving the initial time
    ofstream tLSCheck;
    tLSCheck.open(FolderName+"/"+"tLSCheck.csv");
    tLSCheck << 0;
    tLSCheck.close();

    // saving LS random generator
    //std::ofstream randStateLS(mainResumeFolderName+"/"+"randStateLS.csv");
    ofstream randStateLS;
    randStateLS.open(FolderName+"/"+"randStateLS.csv");
    randStateLS << mt_rand;
    randStateLS.close();


    // saving oldRoundCounter to be zero
    ofstream roundCounterLS;
    roundCounterLS.open(FolderName+"/"+"roundCounterLS.csv");
    roundCounterLS << 0;
    roundCounterLS.close();

    // saving the initial time
    ofstream tLS;
    tLS.open(FolderName+"/"+"tLS.csv");
    tLS << 0;
    // tLSValue=0;
    tLS.close();
}
//
void LinkedListSiteCreation(vector<LinkedListElement*>& LLheadVec, vector<LinkedListElement*>& LLtailVec, vector<int>& LLNumVec,\
                            const vector<int>& sigmaMat, vector<int>& isBorder,\
                            const vector<int>& neighborNum, const vector<vector<int>>& neighborsList,\
                            const vector<double>& siteComX, const double Lx)
{

  // vector<std::pair<int, int>> test_vec;
  // std::pair<int, int> test_element;
  
  int NSites = sigmaMat.size();
  int numLinkedList = LLheadVec.size();

  vector<LinkedListElement*> LLTempVec(numLinkedList);
  LinkedListElement *LLTemp;
  for (int i = 0; i < numLinkedList; i++)
  {
    LLNumVec[i] = 0;
    LLheadVec[i] = NULL;
    LLTempVec[i] = NULL;
  }
  
  int siteC, NNeigh, neighIndC, neighInd;
  for (siteC = 0; siteC < NSites; siteC++)
  {
    isBorder[siteC] = 0;
    NNeigh = neighborNum[siteC];
    for (neighIndC = 0; neighIndC < NNeigh; neighIndC++)
    {
      neighInd = neighborsList[siteC][neighIndC];
      if ( sigmaMat[siteC]!=sigmaMat[neighInd] )
      {
        isBorder[siteC] = 1;
      }
    }
  }


  int LLIndex;
  for (siteC = 0; siteC < NSites; siteC++)
  {
    if ( isBorder[siteC] )
    {
      LLIndex = (int)(floor(1.0*siteComX[siteC]/(1.0*Lx/numLinkedList)));
      if (LLNumVec[LLIndex]>0)
      {
        LLTemp = LLTempVec[LLIndex];
        LLTemp -> next = new LinkedListElement;
        LLTemp = LLTemp -> next;
        LLTemp -> site  = siteC;
        LLTemp -> next = NULL;
        LLTempVec[LLIndex] = LLTemp;
        LLtailVec[LLIndex] = LLTemp;
      }
      else if (LLNumVec[LLIndex]==0)
      {
        LLTemp = new LinkedListElement;
        LLTemp -> site  = siteC;
        LLTemp -> next = NULL;
        LLTempVec[LLIndex] = LLTemp;
        LLtailVec[LLIndex] = LLTemp;
        LLheadVec[LLIndex] = LLTemp;
      }
      LLNumVec[LLIndex]++;
    }
  }
}
//
void initializer(vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                 const vector<double>& siteComX, const vector<double>& siteComY, \
                 int NumCells, double AvgCellArea, double Lx, double Ly)
{
  // int L = sigmaMat.size();
  int NSites = sigmaMat.size();

  int numCellsXdirection = sqrt(NumCells); 
  int numCellsYdirection = sqrt(NumCells); 

  int row, col;

  double cellXLength = Lx/numCellsXdirection;
  double cellYLength = Ly/numCellsYdirection;
  
  int xSite, ySite;

  int index = 1;
  for (int siteC=0; siteC<NSites; siteC++)
  {
    xSite = siteComX[siteC];
    ySite = siteComY[siteC];
    
    sitesX[siteC] = siteComX[siteC];
    sitesY[siteC] = siteComY[siteC];

    row = (int)(floor(xSite/cellXLength));
    col = (int)(floor(ySite/cellYLength));

    index = (row * numCellsYdirection + col) + 1;

    sigmaMat[siteC] = index;
  }
}
//
void hexInitializer(vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                    const vector<double>& siteComX, const vector<double>& siteComY, \
                    int NumCells, double AvgCellArea, double Lx, double Ly)
{
  // int L = sigmaMat.size();
  
  if ( fabs(AvgCellArea*NumCells-(Lx*Ly)) > (1e-6) )
  {
    cout<<"#################################"<<endl;
    cout<<"A_0 mismatch in hexagonal initialization!"<<endl;
    cout<<"Program ended!"<<endl;
    cout<<"#################################"<<endl;
    exit(0);
  }

  
  int n_row = 34;
  int n_col = 29;

  if ( n_row%2 == 1 )
  {
    cout<<"#################################"<<endl;
    cout<<"n_row must be even in hexagonal initialization!"<<endl;
    cout<<"Program ended!"<<endl;
    cout<<"#################################"<<endl;
    exit(0);
  }

  if ( fabs(n_row*n_col-NumCells) > (1e-6) )
  {
    cout<<"#################################"<<endl;
    cout<<"(n_row, n_col) mismatch in hexagonal initialization!"<<endl;
    cout<<"Program ended!"<<endl;
    cout<<"#################################"<<endl;
    exit(0);
  }

  // cout<<"**************************"<<endl;
  // cout<<""<<endl;
  // cout<<"**************************"<<endl;

  int NSites = sigmaMat.size();
  int site_c;
  int row_c, col_c;
  double a, h;
  // a = sqrt(2*AvgCellArea/(3*sqrt(3.0)));
  // h = a * sqrt(3.0)/2;
  a = Lx/(3*n_row/2);
  h = Ly/(2*n_col);

  double y_min_tile, y_max_tile;
  double x_center, y_center;
  vector<double> line_ru(2); // [m , d] s.t. "y = m x + d"
  vector<double> line_lu(2);
  vector<double> line_rd(2);
  vector<double> line_ld(2);


  for (site_c = 0; site_c < NSites; site_c++)
  {
    sigmaMat[site_c] = -1;
  }

  int cell_ind = 1;
  bool condition;
  double x_row;
  for (row_c = 0; row_c < n_row; row_c++)
  {

    if (row_c==0) // upper row
    {
      for (col_c = 0; col_c < n_col; col_c++)
      {
        y_min_tile = col_c*(2.0*h);
        y_max_tile = (col_c+1)*(2.0*h);
        y_center = 0.5* (y_min_tile + y_max_tile);

        x_center = 0.0;
        line_rd[0] = -2.0*h/a; // [m , d] s.t. "y = m x + d"
        line_rd[1] = y_center - (line_rd[0])*(x_center+a);
        line_ld[0] = +2.0*h/a;
        line_ld[1] = y_center - (line_ld[0])*(x_center+a);

        x_center = Lx;
        line_ru[0] = +2.0*h/a;
        line_ru[1] = y_center - (line_ru[0])*(x_center-a);
        line_lu[0] = -2.0*h/a;
        line_lu[1] = y_center - (line_lu[0])*(x_center-a);

        for (site_c = 0; site_c < NSites; site_c++)
        {
          if (sigmaMat[site_c]>=0)
          {
            continue;
          }
          condition =((siteComX[site_c]<=a) \
                   && (siteComY[site_c]<=y_max_tile) && (siteComY[site_c]>=y_min_tile) \
                   && (siteComY[site_c]<=line_rd[0]*siteComX[site_c]+line_rd[1]) && (siteComY[site_c]>=line_ld[0]*siteComX[site_c]+line_ld[1])) \
                   || \
                     ((siteComX[site_c]>=Lx-a) \
                   && (siteComY[site_c]<=y_max_tile) && (siteComY[site_c]>=y_min_tile) \
                   && (siteComY[site_c]<=line_ru[0]*siteComX[site_c]+line_ru[1]) && (siteComY[site_c]>=line_lu[0]*siteComX[site_c]+line_lu[1]));
          if (condition)
          {
            sigmaMat[site_c]=cell_ind;
          }
        }
        cell_ind++;        
      }
    }


    else // row_c>0
    {

      x_center = (3*a/2 + (3*a) * (row_c-1)/2) * (row_c%2==1) \ 
               + (        (3*a) * (row_c  )/2) * (row_c%2==0) ;

      if (row_c%2==1)
      {
        // y_min_tile = 0.0;
        // y_max_tile = h;
        y_center = 0.0;
        line_ru[0] = +2.0*h/a;
        line_ru[1] = y_center - (line_ru[0])*(x_center-a);
        line_rd[0] = -2.0*h/a;
        line_rd[1] = y_center - (line_rd[0])*(x_center+a);

        // y_min_tile = Ly-h;
        // y_max_tile = Ly;
        y_center = Ly;
        line_lu[0] = -2.0*h/a;
        line_lu[1] = y_center - (line_lu[0])*(x_center-a);
        line_ld[0] = +2.0*h/a;
        line_ld[1] = y_center - (line_ld[0])*(x_center+a);

        for (site_c = 0; site_c < NSites; site_c++)
        {
          if (sigmaMat[site_c]>=0)
          {
            continue;
          }
          condition =((siteComY[site_c]<=h) \
                   && (siteComY[site_c]<=line_ru[0]*siteComX[site_c]+line_ru[1]) && (siteComY[site_c]<=line_rd[0]*siteComX[site_c]+line_rd[1])) \
                   || \
                     ((siteComY[site_c]>=Ly-h) \
                   && (siteComY[site_c]>=line_lu[0]*siteComX[site_c]+line_lu[1]) && (siteComY[site_c]>=line_ld[0]*siteComX[site_c]+line_ld[1]));
          if (condition)
          {
            sigmaMat[site_c]=cell_ind;
          }
        }
        cell_ind++;
      }
      
      // ordinary tiles
      y_min_tile = (h)*(row_c%2==1)+(0.0)*(row_c%2==0);
      y_max_tile = y_min_tile + 2*h;
      y_center = 0.5* (y_min_tile + y_max_tile);

      while ( y_max_tile-Ly < (1e-6))
      {
        
        line_ru[0] = +2.0*h/a;
        line_ru[1] = y_center - (line_ru[0])*(x_center-a);
        line_rd[0] = -2.0*h/a;
        line_rd[1] = y_center - (line_rd[0])*(x_center+a);

        line_lu[0] = -2.0*h/a;
        line_lu[1] = y_center - (line_lu[0])*(x_center-a);
        line_ld[0] = +2.0*h/a;
        line_ld[1] = y_center - (line_ld[0])*(x_center+a);

        for (site_c = 0; site_c < NSites; site_c++)
        {
          if (sigmaMat[site_c]>=0)
          {
            continue;
          }
          condition = (siteComY[site_c]<=y_max_tile) \
                   && (siteComY[site_c]<=line_ru[0]*siteComX[site_c]+line_ru[1]) \
                   && (siteComY[site_c]<=line_rd[0]*siteComX[site_c]+line_rd[1]) \
                   && (siteComY[site_c]>=y_min_tile) \
                   && (siteComY[site_c]>=line_lu[0]*siteComX[site_c]+line_lu[1]) \
                   && (siteComY[site_c]>=line_ld[0]*siteComX[site_c]+line_ld[1]);

          if (condition)
          {
            sigmaMat[site_c]=cell_ind;
          }
        }
        cell_ind++;

        y_min_tile += (2*h);
        y_max_tile += (2*h);
        y_center   += (2*h);
      }
    }
  }
  
  for (site_c = 0; site_c < NSites; site_c++)
  {
    sitesX[site_c] = siteComX[site_c];
    sitesY[site_c] = siteComY[site_c];

    if (sigmaMat[site_c]<0)
    {
      cout<<"#################################"<<endl;
      cout<<"There is a -1 in hex initialization"<<endl;
      cout<<"Program ended!"<<endl;
      cout<<"#################################"<<endl;
      exit(0);
    }
    

    if ( (sigmaMat[site_c] <= n_col) && siteComX[site_c] > 0.5*Lx )
    {
      sitesX[site_c] = siteComX[site_c] - Lx;
    }

    if ( (sigmaMat[site_c] % (2*n_col) == (n_col+1)) && siteComY[site_c] > 0.5*Ly)
    {
      sitesY[site_c] = siteComY[site_c] - Ly;
    }
  }
}
//
void vacInitializer(vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                    const vector<double>& siteComX, const vector<double>& siteComY, \
                    int NumCells, double AvgCellArea, double Lx, double Ly)
{
  // int L = sigmaMat.size();
  int NSites = sigmaMat.size();

 
  
  int NumCellMax = (int)floor((Lx*Ly+(1e-6))/AvgCellArea);

  int numCellsXdirection = sqrt(NumCellMax); 
  int numCellsYdirection = sqrt(NumCellMax);
  if (numCellsXdirection*numCellsYdirection != NumCellMax)
  {
    numCellsXdirection = (int)floor(sqrt(1.0*NumCellMax));
    while (NumCellMax % numCellsXdirection)
    {
      numCellsXdirection++;
    }
    numCellsYdirection = NumCellMax / numCellsXdirection;
  }

  int row, col;

  double cellXLength = Lx/numCellsXdirection;
  double cellYLength = Ly/numCellsYdirection;
  
  int xSite, ySite;

  int index;
  for (int siteC=0; siteC<NSites; siteC++)
  {
    xSite = siteComX[siteC];
    ySite = siteComY[siteC];
    
    sitesX[siteC] = siteComX[siteC];
    sitesY[siteC] = siteComY[siteC];

    row = (int)(floor(xSite/cellXLength));
    col = (int)(floor(ySite/cellYLength));

    index = (row * numCellsYdirection + col) + 1;
    if (index <= NumCells)
    {
      sigmaMat[siteC] = index;
    } else
    {
      sigmaMat[siteC] = 0;
    }
  }

  // for (int siteC=0; siteC<NSites; siteC++)
  // {
  //   xSite = siteComX[siteC];
  //   ySite = siteComY[siteC];

  //   if (xSite>0 && xSite<sqrt(AvgCellArea) && ySite>0 && ySite<sqrt(AvgCellArea))
  //   {
  //     sitesX[siteC] = siteComX[siteC];
  //     sitesY[siteC] = siteComY[siteC];
  //     sigmaMat[siteC] = 1;
  //   }
  //   else
  //   {
  //     sitesX[siteC] = siteComX[siteC];
  //     sitesY[siteC] = siteComY[siteC];
  //     sigmaMat[siteC] = 0;
  //   }
  // }
}
//
void recConfluentInitializer(vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                    const vector<double>& siteComX, const vector<double>& siteComY, \
                    int NumCells, double AvgCellArea, double Lx, double Ly)
{
  // int L = sigmaMat.size();
  int NSites = sigmaMat.size();

 

  int numCellsXdirection = sqrt(NumCells); 
  int numCellsYdirection = sqrt(NumCells);
  if (numCellsXdirection*numCellsYdirection != NumCells)
  {
    numCellsXdirection = (int)floor(sqrt(1.0*NumCells)) + 1;
    while (NumCells % numCellsXdirection)
    {
      numCellsXdirection++;
    }
    numCellsYdirection = NumCells / numCellsXdirection;
  }

  int row, col;

  double cellXLength = Lx/numCellsXdirection;
  double cellYLength = Ly/numCellsYdirection;
  
  int xSite, ySite;

  int index;
  for (int siteC=0; siteC<NSites; siteC++)
  {
    xSite = siteComX[siteC];
    ySite = siteComY[siteC];
    
    sitesX[siteC] = siteComX[siteC];
    sitesY[siteC] = siteComY[siteC];

    row = (int)(floor(xSite/cellXLength));
    col = (int)(floor(ySite/cellYLength));

    index = (row * numCellsYdirection + col) + 1;
    if (index <= NumCells && index>0)
    {
      sigmaMat[siteC] = index;
    } else
    {
      cout<< "Error in recConfluentInitializer! " <<endl;
      exit(9);
    }
  }

  // for (int siteC=0; siteC<NSites; siteC++)
  // {
  //   xSite = siteComX[siteC];
  //   ySite = siteComY[siteC];

  //   if (xSite>0 && xSite<sqrt(AvgCellArea) && ySite>0 && ySite<sqrt(AvgCellArea))
  //   {
  //     sitesX[siteC] = siteComX[siteC];
  //     sitesY[siteC] = siteComY[siteC];
  //     sigmaMat[siteC] = 1;
  //   }
  //   else
  //   {
  //     sitesX[siteC] = siteComX[siteC];
  //     sitesY[siteC] = siteComY[siteC];
  //     sigmaMat[siteC] = 0;
  //   }
  // }
}
//
void sitesXsitesYInit(const vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                    const vector<double>& siteComX, const vector<double>& siteComY, \
                    const int NumCells, const double AvgCellArea, const double Lx, const double Ly)
{
  // int L = sigmaMat.size();
  int NSites = sigmaMat.size();

  int xSite, ySite;


  vector<double> comXEstimate(NumCells + 1);
  vector<double> comYEstimate(NumCells + 1);
  vector<int> cellSiteNum(NumCells + 1);

  for (int cell_c = 0; cell_c <= NumCells; cell_c++)
  {
    comXEstimate[cell_c] = 0.0; // estimation of cumulative COM of cells
    comYEstimate[cell_c] = 0.0;
    cellSiteNum[cell_c]  = 0;
  }
  
  int cell_ind;
  for (int siteC=0; siteC<NSites; siteC++)
  {
    cell_ind = sigmaMat[siteC];

    if (cellSiteNum[cell_ind] == 0)
    {
      sitesX[siteC] = siteComX[siteC];
      sitesY[siteC] = siteComY[siteC];

      comXEstimate[cell_ind] = siteComX[siteC];
      comYEstimate[cell_ind] = siteComY[siteC];

      cellSiteNum[cell_ind]++;
    }
    else // cellSiteNum[cell_ind] > 0
    {
      if (siteComX[siteC] - comXEstimate[cell_ind] > Lx/2)
      {
        sitesX[siteC] = siteComX[siteC] - Lx;
      } else if (comXEstimate[cell_ind] - siteComX[siteC] > Lx/2)
      {
        sitesX[siteC] = siteComX[siteC] + Lx;
      } else
      {
        sitesX[siteC] = siteComX[siteC];
      }

      if (siteComY[siteC] - comYEstimate[cell_ind] > Ly/2)
      {
        sitesY[siteC] = siteComY[siteC] - Ly;
      } else if (comYEstimate[cell_ind] - siteComY[siteC] > Ly/2)
      {
        sitesY[siteC] = siteComY[siteC] + Ly;
      } else
      {
        sitesY[siteC] = siteComY[siteC];
      }

      comXEstimate[cell_ind] = (comXEstimate[cell_ind]*cellSiteNum[cell_ind]+sitesX[siteC])/(1+cellSiteNum[cell_ind]);
      comYEstimate[cell_ind] = (comYEstimate[cell_ind]*cellSiteNum[cell_ind]+sitesY[siteC])/(1+cellSiteNum[cell_ind]);

      cellSiteNum[cell_ind]++;
    }
  }
}
//
void compartInitializer(std::mt19937 &mt_rand, const vector<int>& sigmaMat, vector<int>& compartMat, const vector<double>& avgAreaFrac)
{
  int N_compart = avgAreaFrac.size();
  int NSites = sigmaMat.size();

  int minRange = 0;
  int maxRange = 10000;
  
  double randUniform;
  double accum;
  for (int siteC = 0; siteC < NSites; siteC++)
  {
    randUniform = 1.0*(minRange+(mt_rand())%(maxRange-minRange+1))/(maxRange-minRange);

    accum = 0; // accumulation of probabilities
    for (int compartC = 0; compartC < N_compart; compartC++)
    {
      if (randUniform <= accum + avgAreaFrac[compartC])
      {
        compartMat[siteC] = compartC;
        break;
      } else
      {
        accum = accum + avgAreaFrac[compartC];
      }
    }    
  }
}
//
void singleInitializer(std::mt19937 &mt_rand, vector<int>& sigmaMat, vector<double>& sitesX, vector<double>& sitesY, \
                    const vector<double>& siteComX, const vector<double>& siteComY, \
                    int NumCells, double AvgCellArea, double Lx, double Ly)
{
  // int L = sigmaMat.size();
  int NSites = sigmaMat.size();
  std::srand(static_cast<unsigned int>(std::time(nullptr))); // Seed with current time
  int minRange = 0;
  int maxRange = 1000;
  
  double x_0, y_0;
  double randUniform;
  randUniform = 1.0*(minRange+(mt_rand())%(maxRange-minRange+1))/(maxRange-minRange);
  // x_0 = 0.5*Lx + (2*randUniform-1.0)*sqrt(AvgCellArea)/2;
  x_0 = randUniform*(Lx-sqrt(AvgCellArea));
  randUniform = 1.0*(minRange+(mt_rand())%(maxRange-minRange+1))/(maxRange-minRange);
  // y_0 = 0.5*Ly + (2*randUniform-1.0)*sqrt(AvgCellArea)/2;
  y_0 = randUniform*(Ly-sqrt(AvgCellArea));

  double x_upper, y_upper, x_lower, y_lower;
  // x_upper = x_0 + 0.5*sqrt(AvgCellArea);
  // x_lower = x_0 - 0.5*sqrt(AvgCellArea); 
  // y_upper = y_0 + 0.5*sqrt(AvgCellArea);
  // y_lower = y_0 - 0.5*sqrt(AvgCellArea); 

  x_upper = x_0 + sqrt(AvgCellArea);
  x_lower = x_0 ; 
  y_upper = y_0 + sqrt(AvgCellArea);
  y_lower = y_0 ; 

  int xSite, ySite;

  int index;
  for (int siteC=0; siteC<NSites; siteC++)
  {
    xSite = siteComX[siteC];
    ySite = siteComY[siteC];
    
    sitesX[siteC] = siteComX[siteC];
    sitesY[siteC] = siteComY[siteC];

    bool condition;
    condition = (xSite<x_upper) && (xSite>x_lower) && (ySite<y_upper) && (ySite>y_lower);

    if (condition)
    {
      sigmaMat[siteC] = 1;
    } else
    {
      sigmaMat[siteC] = 0;
    }
  }

  // for (int siteC=0; siteC<NSites; siteC++)
  // {
  //   xSite = siteComX[siteC];
  //   ySite = siteComY[siteC];

  //   if (xSite>0 && xSite<sqrt(AvgCellArea) && ySite>0 && ySite<sqrt(AvgCellArea))
  //   {
  //     sitesX[siteC] = siteComX[siteC];
  //     sitesY[siteC] = siteComY[siteC];
  //     sigmaMat[siteC] = 1;
  //   }
  //   else
  //   {
  //     sitesX[siteC] = siteComX[siteC];
  //     sitesY[siteC] = siteComY[siteC];
  //     sigmaMat[siteC] = 0;
  //   }
  // }
}
//
void crystalInitializer(vector<vector<int>>& sigmaMat, vector<vector<int>>& sitesX, vector<vector<int>>& sitesY, int NumCells, double AvgCellArea)
{
  int L = sigmaMat.size();

  double eps = 1.0e-6;
  int tileXLength;
  int tileYLength;
  if (L==200 && NumCells==1000)
  {
    tileXLength = 10;
    tileYLength = 8;
  }else if (L==40 && NumCells==40)
  {
    tileXLength = 10;
    tileYLength = 8;
  }else
  {
    cout<<"crystal initializer must be corrected!"<<endl;
    cout<<"Program ended!"<<endl;
    exit(0);
  }

  int tileXCounter;
  int tileYCounter;

  int tilesNumX = (L/tileXLength);
  int tilesNumY = (L/tileYLength);

  int xCounter; //from 0 to L
  int yCounter; //from 0 to L
  int xCounterInTile; //from 0 to tileXLength
  int yCounterInTile; //from 0 to tileYLength
  
  int ind;

  int correctionSiteX, correctionSiteY;

  int correctionSiteX_NW;
  int correctionSiteX_NE;
  int correctionSiteX_OO;
  int correctionSiteX_SW;
  int correctionSiteX_SE;

  int correctionSiteY_NW;
  int correctionSiteY_NE;
  int correctionSiteY_OO;
  int correctionSiteY_SW;
  int correctionSiteY_SE;

  int ind_NW;
  int ind_NE;
  int ind_OO;
  int ind_SW;
  int ind_SE;

  ind_NW = 0;
  ind_NE = ind_NW+1;
  ind_OO = tilesNumY;
  ind_SW = 2*tilesNumY;
  ind_SE = ind_SW+1;

  int tileCounter =0;

  for (tileXCounter = 0; tileXCounter < tilesNumX; tileXCounter++)
  {
    for (tileYCounter = 0; tileYCounter < tilesNumY; tileYCounter++)
    {

      if (tileCounter==0)
      {
        // Values are already correct
          correctionSiteX_NW = 0;
          correctionSiteX_NE = 0;
          correctionSiteX_OO = 0;
          correctionSiteX_SW = 0;
          correctionSiteX_SE = 0;

          correctionSiteY_NW = 0;
          correctionSiteY_NE = 0;
          correctionSiteY_OO = 0;
          correctionSiteY_SW = 0;
          correctionSiteY_SE = 0;

      } else if (tileCounter<tilesNumY-1)
      {
        ind_NW = ind_NE;
        ind_SW = ind_SE;
        ind_OO ++;
        ind_NE ++;
        ind_SE ++;

        correctionSiteX_NW = 0;
        correctionSiteX_NE = 0;
        correctionSiteX_OO = 0;
        correctionSiteX_SW = 0;
        correctionSiteX_SE = 0;

        correctionSiteY_NW = 0;
        correctionSiteY_NE = 0;
        correctionSiteY_OO = 0;
        correctionSiteY_SW = 0;
        correctionSiteY_SE = 0;

      } else if (tileCounter==tilesNumY-1)
      {
        ind_NW = ind_NE;
        ind_SW = ind_SE;
        ind_OO ++;
        ind_NE = (ind_NE+1)-tilesNumY;
        ind_SE = (ind_SE+1)-tilesNumY;

        correctionSiteX_NW = 0;
        correctionSiteX_NE = 0;
        correctionSiteX_OO = 0;
        correctionSiteX_SW = 0;
        correctionSiteX_SE = 0;

        correctionSiteY_NW = 0;
        correctionSiteY_NE = -L;
        correctionSiteY_OO = 0;
        correctionSiteY_SW = 0;
        correctionSiteY_SE = -L;

      } else if (tileYCounter==0 && tileXCounter<tilesNumX-1)
      {
        ind_NW = ind_NW+tilesNumY+1;
        ind_NE = ind_NW+1;
        ind_OO = ind_OO+tilesNumY+1;
        ind_SW = ind_OO+tilesNumY;
        ind_SE = ind_SW+1;

        correctionSiteX_NW = 0;
        correctionSiteX_NE = 0;
        correctionSiteX_OO = 0;
        correctionSiteX_SW = 0;
        correctionSiteX_SE = 0;

        correctionSiteY_NW = 0;
        correctionSiteY_NE = 0;
        correctionSiteY_OO = 0;
        correctionSiteY_SW = 0;
        correctionSiteY_SE = 0;

      } else if (tileYCounter==0 && tileXCounter==tilesNumX-1)
      {
        ind_NW = ind_NW+tilesNumY+1;
        ind_NE = ind_NW+1;
        ind_OO = ind_OO+tilesNumY+1;
        ind_SW = 0;
        ind_SE = 1;

        correctionSiteX_NW = 0;
        correctionSiteX_NE = 0;
        correctionSiteX_OO = 0;
        correctionSiteX_SW = -L;
        correctionSiteX_SE = -L;

        correctionSiteY_NW = 0;
        correctionSiteY_NE = 0;
        correctionSiteY_OO = 0;
        correctionSiteY_SW = 0;
        correctionSiteY_SE = 0;

      } else if (tileYCounter==tilesNumY-1 && tileXCounter<tilesNumX-1)
      {
        ind_NW = ind_NE;
        ind_SW = ind_SE;
        ind_OO ++;
        ind_NE =(ind_NE+1)-tilesNumY;
        ind_SE =(ind_SE+1)-tilesNumY;

        correctionSiteX_NW = 0;
        correctionSiteX_NE = 0;
        correctionSiteX_OO = 0;
        correctionSiteX_SW = 0;
        correctionSiteX_SE = 0;

        correctionSiteY_NW = 0;
        correctionSiteY_NE = -L;
        correctionSiteY_OO = 0;
        correctionSiteY_SW = 0;
        correctionSiteY_SE = -L;

      } else if (tileYCounter==tilesNumY-1 && tileXCounter==tilesNumX-1)
      {
        ind_NW = ind_NE;
        ind_SW = ind_SE;
        ind_OO ++;
        ind_NE =(ind_NE+1)-tilesNumY;
        ind_SE =(ind_SE+1)-tilesNumY;

        correctionSiteX_NW = 0;
        correctionSiteX_NE = 0;
        correctionSiteX_OO = 0;
        correctionSiteX_SW = -L;
        correctionSiteX_SE = -L;

        correctionSiteY_NW = 0;
        correctionSiteY_NE = -L;
        correctionSiteY_OO = 0;
        correctionSiteY_SW = 0;
        correctionSiteY_SE = -L;

      } else if (tileXCounter==tilesNumX-1)
      {
        ind_NW = ind_NE;
        ind_SW = ind_SE;
        ind_OO ++;
        ind_NE ++;
        ind_SE = 1+tileYCounter;

        correctionSiteX_NW = 0;
        correctionSiteX_NE = 0;
        correctionSiteX_OO = 0;
        correctionSiteX_SW = -L;
        correctionSiteX_SE = -L;

        correctionSiteY_NW = 0;
        correctionSiteY_NE = 0;
        correctionSiteY_OO = 0;
        correctionSiteY_SW = 0;
        correctionSiteY_SE = 0;

      } else
      {
        ind_NW = ind_NE;
        ind_SW = ind_SE;
        ind_OO ++;
        ind_NE ++;
        ind_SE ++;

        correctionSiteX_NW = 0;
        correctionSiteX_NE = 0;
        correctionSiteX_OO = 0;
        correctionSiteX_SW = 0;
        correctionSiteX_SE = 0;

        correctionSiteY_NW = 0;
        correctionSiteY_NE = 0;
        correctionSiteY_OO = 0;
        correctionSiteY_SW = 0;
        correctionSiteY_SE = 0;
      }
      
      

      for (xCounterInTile = 0; xCounterInTile < tileXLength; xCounterInTile++)
      {
        for (yCounterInTile = 0; yCounterInTile < tileYLength; yCounterInTile++)
        {
          
          if      (yCounterInTile < ((tileYLength/2.0)-0.5) - xCounterInTile)
            {
              ind = ind_NW;
              correctionSiteX = correctionSiteX_NW;
              correctionSiteY = correctionSiteY_NW;
            }
          else if (yCounterInTile > ((tileYLength/2.0)-0.5) + (tileXLength-1.0) - xCounterInTile)
            {
              ind = ind_SE;
              correctionSiteX = correctionSiteX_SE;
              correctionSiteY = correctionSiteY_SE;
            }
          else if (yCounterInTile < ((tileYLength/2.0)-0.5) - (tileXLength-1.0) + xCounterInTile)
            {
              ind = ind_SW;
              correctionSiteX = correctionSiteX_SW;
              correctionSiteY = correctionSiteY_SW;
            }
          else if (yCounterInTile > ((tileYLength/2.0)-0.5) + xCounterInTile)
            {
              ind = ind_NE;
              correctionSiteX = correctionSiteX_NE;
              correctionSiteY = correctionSiteY_NE;
            }
          else
            {
              ind = ind_OO;
              correctionSiteX = correctionSiteX_OO;
              correctionSiteY = correctionSiteY_OO;
            }
          

          xCounter = tileXCounter*tileXLength+xCounterInTile;
          yCounter = tileYCounter*tileYLength+yCounterInTile;

          sigmaMat[xCounter][yCounter] = ind+1; // +1 because I was wrong and started the indices from 0
          sitesX[xCounter][yCounter] = xCounter+correctionSiteX;
          sitesY[xCounter][yCounter] = yCounter+correctionSiteY;

        }
        
      }
      tileCounter++;
    }
  }
}
//
void areaCalc(const vector<int>& sigmaMat, vector<double>& cellArea, vector<int>& cellSiteNum, const vector<double>& latticeArea)
{
  
  // int NumCells = cellArea.size() - 1;
  int NSites = sigmaMat.size();

  // cellArea[0] = 0;
  // cellSiteNum[0] = 0;
  int cellIndex;
  for (cellIndex = 0; cellIndex < cellArea.size(); cellIndex++)
  {
    cellArea[cellIndex] = 0;
    cellSiteNum[cellIndex] = 0;
  }


  for (int siteC = 0; siteC < NSites; siteC++)
  {
    cellIndex = sigmaMat[siteC];
    cellArea[cellIndex] += latticeArea[siteC];
    cellSiteNum[cellIndex] += 1;
  }
}
//
void compartAreaCalc(const vector<int>& sigmaMat, const vector<int>& compartMat, vector<vector<double>>& compartArea, const vector<double>& latticeArea)
{
  
  
  int NSites = sigmaMat.size();
  int NumCells = compartArea.size() - 1;
  int NComparts = compartArea[0].size();

  // cellArea[0] = 0;
  // cellSiteNum[0] = 0;
  int cellIndex;
  int compartIndex;
  for (cellIndex = 0; cellIndex <= NumCells; cellIndex++)
  {
    for (compartIndex = 0; compartIndex < NComparts; compartIndex++)
    {
      compartArea[cellIndex][compartIndex] = 0.0;
    }
  }


  for (int siteC = 0; siteC < NSites; siteC++)
  {
    cellIndex = sigmaMat[siteC];
    compartIndex = compartMat[siteC];

    compartArea[cellIndex][compartIndex] += latticeArea[siteC];
    // cellSiteNum[cellIndex] += 1;
  }
}
//
void saveSigmaMatCSV(const vector<vector<int>>& sigmaMat, const std::string &filename)
{
  
  int L = sigmaMat.size();

  std::ofstream outfile(filename);

  for (int i = 0; i < L; i++)
  {
    for (int j = 0; j < L; j++)
    {
      outfile << sigmaMat[i][j];

      if (j != L - 1)
      {
        outfile << ",";
      }
    }
    outfile << std::endl;
  }
  outfile.close();
}

void saveDoubleMatCSV(const vector<vector<double>>& sigmaMat, const std::string &filename)
{
  
  int Nx = sigmaMat.size();
  int Ny = sigmaMat[0].size();

  std::ofstream outfile(filename);

  for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {
      outfile << std::setprecision(10) << sigmaMat[i][j];

      if (j != Ny - 1)
      {
        outfile << ",";
      }
    }
    outfile << std::endl;
  }
  outfile.close();
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

void saveSampleAP(const vector<vector<int>>& AP, const std::string &filename)
{

  int NumCells = AP.size() - 1;
  int samplesPerWrite = AP[0].size();

  std::ofstream outfile(filename);

  for (int i = 0; i < NumCells + 1; i++)
  {
    for (int j = 0; j < samplesPerWrite; j++)
    {
      outfile << AP[i][j];

      if (j != samplesPerWrite - 1)
      {
        outfile << ",";
      }
    }
    outfile << std::endl;
  }
  outfile.close();
}

void saveSampleT(const vector<int>& quantityTemporal, const std::string &filename)
{

  int samplesPerWrite = quantityTemporal.size();

  std::ofstream outfile(filename);

  for (int i = 0; i < samplesPerWrite; i++)
  {
    outfile << quantityTemporal[i];
    if (i < samplesPerWrite - 1)
    {
      outfile << ",";
    }
  }
  outfile.close();
}

void saveSampleE(const vector<double>& quantityTemporal, const std::string &filename)
{
  
  int samplesPerWrite = quantityTemporal.size();

  std::ofstream outfile(filename);

  for (int i = 0; i < samplesPerWrite; i++)
  {
    outfile << quantityTemporal[i];
    if (i < samplesPerWrite - 1)
    {
      outfile << ",";
    }
  }
  outfile.close();
}

void saveSnapshots(const vector<vector<vector<int>>>& snapshots, const std::string &filename)
{
  int L = snapshots.size();
  int samplesPerWrite = snapshots[0][0].size();

  std::ofstream outfile(filename);
  int sampleCounter;
  for (sampleCounter = 0; sampleCounter < samplesPerWrite; sampleCounter++)
  {
    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        outfile << snapshots[i][j][sampleCounter];

        if (j != L - 1)
        {
          outfile << ",";
        }
      }
      outfile << std::endl;
    }

    outfile << "------------------------------------------------------";
    outfile << std::endl;
  }
  outfile.close();
}

// void comCalc(const int sigmaMat[L][L], const int sitesX[L][L], const int sitesY[L][L], const int cellArea[NumCells + 1],
//              double comXDummy[NumCells + 1], double comYDummy[NumCells + 1])
// {

//   int row, col, cellIndex;
//   // double xSite, ySite;

//   for (cellIndex = 0; cellIndex < NumCells + 1; cellIndex++)
//   {
//     comXDummy[cellIndex] = 0.;
//     comYDummy[cellIndex] = 0.;
//   }
//   ////////// The origin of the coordinates is considered to be at the //////////
//   ////////// center of the site [row=0, col=0] //////////
//   ////////// so, the 'xStie' and 'ySite' are simply equal to 'row' and 'col' //////////
//   for (row = 0; row < L; row++)
//   {
//     for (col = 0; col < L; col++)
//     {
//       // xSite = (double)(sitesX[row][col]); //It is very important to cast them to double
//       // ySite = (double)(sitesY[row][col]); //It is very important to cast them to double
//       cellIndex = sigmaMat[row][col];
//       // comXDummy[cellIndex]+=xSite/cellArea[cellIndex];
//       // comYDummy[cellIndex]+=ySite/cellArea[cellIndex];
//       comXDummy[cellIndex] += ((double)(sitesX[row][col])) / cellArea[cellIndex];
//       comYDummy[cellIndex] += ((double)(sitesY[row][col])) / cellArea[cellIndex];
//     }
//   }
// }
//
void comCalcIrr(const vector<int>& sigmaMat, const vector<double>& sitesX, const vector<double>& sitesY, const vector<double>& cellArea, \
                const vector<double>& latticeArea, \
                vector<double>& cellComX, vector<double>& cellComY)
{
  int NumCells = cellComX.size() - 1;
  int NSites = sigmaMat.size();

  for (int cellC = 0; cellC <= NumCells; cellC++)
  {
    cellComX[cellC] = 0.0;
    cellComY[cellC] = 0.0;
  }

  int cellIdx;

  for (int siteC = 0; siteC < NSites; siteC++)
  {
    cellIdx = sigmaMat[siteC];
    cellComX[cellIdx] += (latticeArea[siteC]*sitesX[siteC]);
    cellComY[cellIdx] += (latticeArea[siteC]*sitesY[siteC]);
  }

  for (int cellC = 1; cellC <= NumCells; cellC++)
  {
    cellComX[cellC] /= cellArea[cellC];
    cellComY[cellC] /= cellArea[cellC];
  }
  
}
//
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

void samplingPatternLoader(vector<vector<long>>& mat, const std::string &filename, int rows, int cols)
{
  std::ifstream input_file(filename);
  if (!input_file)
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }

  std::vector<std::vector<unsigned long>> data;
  std::string line;
  while (getline(input_file, line))
  {
    std::vector<unsigned long> row;
    std::stringstream ss(line);
    unsigned long value;
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

inline bool existFunc (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

void simulationDataReader(int* NSitesPtr, double* LxPtr, double* LyPtr, double* AlphaPtr, double* KbPtr, double* TemPtr,  int* NumCellsPtr, \
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
//
void activityDataReader(double* DrPtr, double* pPtr)
{  
  fstream newfile;
  newfile.open("activityData_vec.csv",ios::in); //open a file to perform read operation using file object
  if (newfile.is_open()){ //checking whether the file is open
    string tp;
    while(getline(newfile, tp)){ //read data from file object and put it into string.
      // cout << tp << "\n"; //print the data of the string
      // if (strstr((const char)tp, "L = "))
      const char* tpChar = tp.c_str();
      if (strstr(tpChar, "Dr = "))
      {
          std::size_t pos = tp.find('=');  // find position of '='
          std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
          *DrPtr = std::stod(num_str);  // convert substring to int
          continue;
      }
      if (strstr(tpChar, "p = "))
      {
          std::size_t pos = tp.find('=');  // find position of '='
          std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
          *pPtr = std::stod(num_str);  // convert substring to int
          continue;
      }
    }
  }
  else
  {
    cout << "Error: Activity file not found."<<endl;
    exit(0);
  }
  newfile.close();
}
//
int samplingPatternLenghFunc()
{
  int length=0;
  // char c;
  fstream newfile;
  newfile.open("samplingPattern_vec.csv",ios::in); //open a file to perform read operation using file object
  if (newfile.is_open()){ //checking whether the file is open
      string tp;
      getline(newfile, tp);

      // Count the number of commas in the string
      for (char c : tp) {
        if (c == ',') {
          length++;
        }
      }
      

  }else
  {
    cout<<"Problem in openning samplingPattern_vec.csv"<<endl;
    cout<<"Program ended!"<<endl;
    exit(0);
  }
  length++;
  newfile.close();
  return length;
}

void latticeCreatorPhase1(int* maxNeighPtr, vector<int>& neighborNum, vector<double>& latticeX,  \
                          vector<double>& latticeY,  vector<double>& latticeArea )
{
    int maxNeigh = 0;
    int neighC = 0;
    int siteC = 0;
    ifstream newfile;
    newfile.open("../lattice/neighs.dat",ios::in);
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
    newfile.open("../lattice/vols.dat",ios::in);
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
                          vector<vector<int>>& neighborsList, vector<vector<double>>& edges, double* avgEdgePtr,\
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
  neighFile.open("../lattice/neighs.dat",ios::in);
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
  ifstream avgEdgeFile;
  avgEdgeFile.open("../lattice/edges_mean.csv",ios::in);
  getline(avgEdgeFile, tp);
  *avgEdgePtr = stod(tp);
  avgEdgeFile.close();

  ifstream wallsFile;
  wallsFile.open("../lattice/walls.dat",ios::in);
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
  distFile.open("../lattice/dists.dat",ios::in);
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
  verticesFile.open("../lattice/vertices.dat",ios::in);
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
  // saving siteComXY
  vector<vector<double>> siteComXY(siteComX.size(), vector<double>(2));
  for (int i = 0; i < NSites; i++)
  {
    siteComXY[i][0] = siteComX[i];
    siteComXY[i][1] = siteComY[i];
  }
  saveDoubleMatCSV(siteComXY, "siteComXY.csv");
  ///////////////////// deltaComX, deltaComY, comDist ////////////////////
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
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// FUNCTIONS ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////


