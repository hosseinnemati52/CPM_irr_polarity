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
#include <random>
#include <chrono>
#include <sys/stat.h>
// #include <filesystem>
#include <string.h>
#include <algorithm>


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

//////////////////////////////////////////////////////////////////////////
//////////////////////////// PROTOTYPES //////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void readCompartData(
    const std::string& filename,
    std::vector<double>& avgAreaFrac,
    std::vector<std::vector<double>>& J_int,
    std::vector<std::vector<double>>& J_ext);

void dataFolderDeleter();

int totalSnapshotsNumFunction(const int samplesPerWrite);

void simulationDataReader(int* NSitesPtr, double* LxPtr, double* LyPtr, double* AlphaPtr, double* KbPtr, double* TemPtr,  int* NumCellsPtr, \
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
const std::vector<double>& latticeArea,\
const std::vector<double>& siteIxx, const std::vector<double>& siteIyy, const std::vector<double>& siteIxy);

void sigmaMatExtractor(vector<vector<vector<int>>>& mat, std::ifstream& samplesFile);

void printForTest(const vector<vector<int>>& mat);

void doubleMatSaver(double* matPtr, int Rows, int Cols, const std::string &filename);

void intMatSaver(int* matPtr, int Rows, int Cols, const std::string &filename);

void double2DVecSaver(const vector<vector<double>>& mat, const std::string &filename);

void double1DVecSaver(const vector<double>& mat, const std::string &filename);

void saveInt1DVec(const vector<int>& quantity, const std::string &filename);

void saveDbl1DVec(const vector<double>& quantity, const std::string &filename);

void loadInt1DVec(vector<int>& mat, const std::string &filename);

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

void S_h_calculator(const vector<double>& xComSample, const vector<double>& yComSample, const double Lx, const double Ly, const vector<double>& h_vec, const vector<double>&  h_ang_vec, \
                    vector<double>& S_h_real_term, vector<double>& S_h_imag_term);

//////////////////////////////////////////////////////////////////////////
//////////////////////////// PROTOTYPES //////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main()
{

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
    
    vector<double> siteIxx(NSites); // Ixx of every lattice site around its COM
    vector<double> siteIyy(NSites); // Iyy of every lattice site around its COM
    vector<double> siteIxy(NSites); // Ixy of every lattice site around its COM
    vector<vector<double>> I_info(NSites, vector<double>(5));
    std::string INamestr = "../lattice/I_info_array.csv";
    doubleMatLoader(I_info, INamestr, NSites, 5);
    for (int siteC = 0; siteC < NSites; siteC++)
    {
      siteIxx[siteC] = I_info[siteC][2];
      siteIyy[siteC] = I_info[siteC][3];
      siteIxy[siteC] = I_info[siteC][4];
      // if (fabs(I_info[siteC][0]-siteComX[siteC])>0.00001 || fabs(I_info[siteC][1]-siteComY[siteC])>0.00001)
      // {
      //   int r;
      //   r=5;
      // }
      // cout<< fabs(I_info[siteC][0]-siteComX[siteC]);
      // cout<< fabs(I_info[siteC][1]-siteComY[siteC]);
    }
    

    // saveSampleE(siteComX, "siteComX.csv");
    // saveSampleE(siteComY, "siteComY.csv");
    // saveDoubleMatCSV(deltaComX, "deltaComX.csv");
    // saveDoubleMatCSV(deltaComY, "deltaComY.csv");
    // saveDoubleMatCSV(comDist, "comDist.csv");

    ////////////////////////////////////////////////////
    /////////////////// LATTICE READING ////////////////
    ////////////////////////////////////////////////////



    string ppDataFolderName = "pp_data";
    /////////////////// MAKING ppData FOLDER /////////////////
    //// This block is for windows:
    // mkdir(ppDataFolderName.c_str()); //making data folder

    //// This block is for Linux:
    mkdir(ppDataFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making backup_resume folder
    /////////////////// MAKING ppData FOLDER /////////////////

    ///// Unzipping the data.zip file /////
    system("python3 dataUnzipper.py all");
    ///// Unzipping the data.zip file /////

    ////////// NUMBER OF SNAPSHOTS /////////////////////
    int totalSnapshotsNum = totalSnapshotsNumFunction(samplesPerWrite);
    ofstream totalSnapshotsNumFile;
    totalSnapshotsNumFile.open(ppDataFolderName + "/" + "totalSnapshotsNum.csv");
    totalSnapshotsNumFile << totalSnapshotsNum;
    totalSnapshotsNumFile.close(); // random seed saved
    ////////// NUMBER OF SNAPSHOTS /////////////////////

    ////////// TIME MAKER /////////////////////
    long time[totalSnapshotsNum];
    // vector<long> time(totalSnapshotsNum);
    timeMaker(time, totalSnapshotsNum, samplesPerWrite);
    ////////// TIME MAKER /////////////////////

    int t_w;
    ifstream t_w_File;
    t_w_File.open("t_w.csv");
    t_w_File >> t_w;
    t_w_File.close();
    
    ////////// EQ SAMPLING TIMES /////////////////////
    int dt_indep = 10 * (time[totalSnapshotsNum-1] - time[totalSnapshotsNum-2]);
    // vector<int> eqSamplingTimesInverse;
    vector<int> eqSamplingTimes;
    eqSamplingTimes.push_back(time[totalSnapshotsNum-1]); // The final config is in the eq samples

    int tUpper = time[totalSnapshotsNum-1];
    int indDummy = totalSnapshotsNum-1;
    int tDummy = time[indDummy];
    while (tDummy>=t_w)
    {
      if (tUpper-tDummy >= dt_indep)
      {
        eqSamplingTimes.insert(eqSamplingTimes.begin(), tDummy);
        tUpper = tDummy;
      }

      indDummy--;
      tDummy = time[indDummy];
    }
    int eqSamplingTimesC = 0; // counter
    ////////// EQ SAMPLING TIMES /////////////////////


    ////////// ENERGY MAKER /////////////////////
    double energy[totalSnapshotsNum];
    energyMaker(energy, totalSnapshotsNum, samplesPerWrite);
    ////////// ENERGY MAKER /////////////////////

    

    /////////////////area, peri, isoperi////////////////

    int snapshotC; //from 0 to (totalSnapshotsNum-1)
    int sampleC; //from 0 to (samplesPerWrite-1)
    int cellInd;
    int bunchC;

    // int sigmaMat[L][L]; // The main field of spins
    // int sitesX[L][L];   // The value of the X(t) of each site
    // int sitesY[L][L];   // The value of the Y(t) of each site
    vector<int> sigmaMat(NSites); // The main field of spins
    vector<double> sitesX(NSites); // The value of the X(t) of each site
    vector<double> sitesY(NSites); // The value of the Y(t) of each site
    vector<int> compartMat(NSites); // The compartments of cells

    // int areaSample[NumCells+1];
    // int periSample[NumCells+1];
    // double isoperiSample[NumCells+1];
    // double xComSample[NumCells+1];
    // double yComSample[NumCells+1];
    // double R2Sample[NumCells+1] = {0} ;
    // double x_w[NumCells+1];
    // double y_w[NumCells+1];
    vector<double> areaSample(NumCells+1);
    vector<double> periSample(NumCells+1);
    vector<double> isoperiSample(NumCells+1);
    vector<double> ARSample(NumCells+1);
    vector<double> circSample(NumCells+1);
    vector<double> xComSample(NumCells+1);
    vector<double> yComSample(NumCells+1);
    vector<double> R2Sample(NumCells+1);
    for(int i=0; i<NumCells+1; i++)
    {
        R2Sample[i]=0;
    }
    vector<double> x_w(NumCells+1);
    vector<double> y_w(NumCells+1);

    double xComTot, yComTot, RComTot, xComTotInit, yComTotInit, xComTotW, yComTotW;

    // double areaAvgStd[2][totalSnapshotsNum];
    // double periAvgStd[2][totalSnapshotsNum];
    // double isoperiAvgStd[2][totalSnapshotsNum];
    // double xComInit[NumCells+1];
    // double yComInit[NumCells+1];
    // double xComBunch[NumCells+1][samplesPerWrite];
    // double yComBunch[NumCells+1][samplesPerWrite];
    // double msdAvgStd[2][totalSnapshotsNum];
    vector<vector<double>> areaAvgStd(2, vector<double>(totalSnapshotsNum));
    vector<double> kurt(totalSnapshotsNum);
    double r2Avg, r4Avg;
    vector<vector<double>> periAvgStd(2, vector<double>(totalSnapshotsNum));
    vector<vector<double>> isoperiAvgStd(2, vector<double>(totalSnapshotsNum));
    vector<vector<double>> ARAvgStd(2, vector<double>(totalSnapshotsNum));
    vector<vector<double>> circAvgStd(2, vector<double>(totalSnapshotsNum));
    
    vector<double> xComInit(NumCells+1);
    vector<double> yComInit(NumCells+1);

    vector<double> xComTotVec(totalSnapshotsNum);
    vector<double> yComTotVec(totalSnapshotsNum);
    vector<double> RComTotVec(totalSnapshotsNum);
    
    vector<vector<double>> xComBunch(NumCells+1, vector<double>(samplesPerWrite));
    vector<vector<double>> yComBunch(NumCells+1, vector<double>(samplesPerWrite));
    vector<vector<double>> msdAvgStd(2, vector<double>(totalSnapshotsNum));

    

    loadInt1DVec(sigmaMat, "init/sigma_init.csv");
    loadInt1DVec(compartMat, "init/compart_init.csv");
    loadDbl1DVec(sitesX,   "init/sitesX_init.csv");
    loadDbl1DVec(sitesY,   "init/sitesY_init.csv");

    // printForTest(sigmaMat);
    // printForTest(sitesX);
    // printForTest(sitesY);

    configurationalCalcPerSample_irr(areaSample, periSample, isoperiSample, ARSample, circSample, xComSample, yComSample, sigmaMat, sitesX, sitesY,\
                                     neighborNum, neighborsList, edges, latticeArea, siteIxx, siteIyy, siteIxy);
    double meanValArea = 0;
    double meanValPeri = 0;
    double meanValIsoperi = 0;
    double meanValAR = 0;
    double meanValCirc = 0;
    for (cellInd = 0; cellInd <= NumCells; cellInd++)
    {   
        if (cellInd > 0)
        {
            meanValArea += ((double)areaSample[cellInd])/NumCells;
            meanValPeri += ((double)periSample[cellInd])/NumCells;
            meanValIsoperi += ((double)isoperiSample[cellInd])/NumCells;
            meanValAR      += ((double)ARSample[cellInd])/NumCells;
            meanValCirc    += ((double)circSample[cellInd])/NumCells;
        }

        xComInit[cellInd] = xComSample[cellInd];
        yComInit[cellInd] = yComSample[cellInd];
    }
    totComFunc(xComTot, yComTot, xComInit, yComInit); // TOTAL COM (INIT)
    xComTotInit = xComTot;
    yComTotInit = yComTot;

    snapshotC = 0;
    areaAvgStd[0][snapshotC] = meanValArea;
    kurt[snapshotC] = -0.5; //By definition
    periAvgStd[0][snapshotC] = meanValPeri;
    isoperiAvgStd[0][snapshotC] = meanValIsoperi;
    ARAvgStd[0][snapshotC] = meanValAR;
    circAvgStd[0][snapshotC] = meanValCirc;
    msdAvgStd[0][snapshotC] = 0;
    doubleMatSaver(&xComInit[0], NumCells+1, 1, ppDataFolderName+"/"+"xComInit.csv");
    doubleMatSaver(&yComInit[0], NumCells+1, 1, ppDataFolderName+"/"+"yComInit.csv");
    xComTotVec[0] = xComTot;
    yComTotVec[0] = yComTot;
    RComTotVec[0] = 0.0;

    double stdValArea = 0;
    double stdValPeri = 0;
    double stdValIsoperi = 0;
    double stdValAR = 0;
    double stdValCirc = 0;
    for (cellInd = 0; cellInd <= NumCells; cellInd++)
    {   
        if (cellInd > 0)
        {
            stdValArea += pow((areaSample[cellInd]-areaAvgStd[0][snapshotC]), 2)/(NumCells);
            stdValPeri += pow((periSample[cellInd]-periAvgStd[0][snapshotC]), 2)/(NumCells);
            stdValIsoperi += pow((isoperiSample[cellInd]-isoperiAvgStd[0][snapshotC]), 2)/(NumCells);
            stdValAR   += pow((ARSample[cellInd]-ARAvgStd[0][snapshotC]), 2)/(NumCells);
            stdValCirc += pow((circSample[cellInd]-circAvgStd[0][snapshotC]), 2)/(NumCells);
        }
    }
    stdValArea = sqrt(stdValArea);
    stdValPeri = sqrt(stdValPeri);
    stdValIsoperi = sqrt(stdValIsoperi);
    stdValAR = sqrt(stdValAR);
    stdValCirc = sqrt(stdValCirc);

    areaAvgStd[1][snapshotC] = stdValArea;
    periAvgStd[1][snapshotC] = stdValPeri;
    isoperiAvgStd[1][snapshotC] = stdValIsoperi;
    ARAvgStd[1][snapshotC] = stdValAR;
    circAvgStd[1][snapshotC] = stdValCirc;
    msdAvgStd[1][snapshotC] = 0;
    



    // // int sigmaMatBunch[L][L][samplesPerWrite]; // The main field of spins
    // // int sitesXBunch[L][L][samplesPerWrite];   // The value of the X(t) of each site
    // // int sitesYBunch[L][L][samplesPerWrite];   // The value of the Y(t) of each site
    // vector<vector<vector<int> > > sigmaMatBunch(L, vector<vector<int> >(L, vector<int>(samplesPerWrite)));
    // vector<vector<vector<int> > >     sitesXBunch(L, vector<vector<int> >(L, vector<int>(samplesPerWrite)));
    // vector<vector<vector<int> > >     sitesYBunch(L, vector<vector<int> >(L, vector<int>(samplesPerWrite)));
    vector<vector<int> >  sigmaMatBunch(NSites, vector<int> (samplesPerWrite));
    // vector<vector<double> >  sitesXBunch(NSites, vector<double> (samplesPerWrite));
    // vector<vector<double> >  sitesYBunch(NSites, vector<double> (samplesPerWrite));
    vector<vector<int> >  sitesXBunch(NSites, vector<int> (samplesPerWrite));
    vector<vector<int> >  sitesYBunch(NSites, vector<int> (samplesPerWrite));

    string namestr, sigmaNamestr, xNamestr, yNamestr;
    bool cond;
    bunchC = 0;
    snapshotC = 1;
    vector<vector<double>> Area_dist_eq(NumCells+1, vector<double>(eqSamplingTimes.size()));
    vector<vector<double>> Peri_dist_eq(NumCells+1, vector<double>(eqSamplingTimes.size()));
    vector<vector<double>> isoperi_dist_eq(NumCells+1, vector<double>(eqSamplingTimes.size()));
    vector<vector<double>> AR_dist_eq(NumCells+1, vector<double>(eqSamplingTimes.size()));
    vector<vector<double>> circ_dist_eq(NumCells+1, vector<double>(eqSamplingTimes.size()));

    /////////// STUFF FOR S(h) AND SISF ////////////////
    double max_h = 2.0 * (2.0 * PI / sqrt(AvgCellArea));
    double dh    = 0.1 * (2.0 * PI / sqrt(AvgCellArea));
    vector<double> h_vec;
    double h_val = dh;
    while (h_val <= max_h)
    {
      h_vec.push_back(h_val);
      h_val += dh;
    }
    
    double dh_ang = 0.13 * (2.0 * PI);
    vector<double> h_ang_vec;
    double h_ang_val = - PI;
    while (h_ang_val <= PI)
    {
      h_ang_vec.push_back(h_ang_val);
      h_ang_val += dh_ang;
    }

    vector<vector<double>> S_h_real(eqSamplingTimes.size(), vector<double>(h_vec.size()));
    vector<vector<double>> S_h_imag(eqSamplingTimes.size(), vector<double>(h_vec.size()));
    vector<double> S_h_real_term(h_vec.size());
    vector<double> S_h_imag_term(h_vec.size());

    for (int j = 0; j < h_vec.size(); j++)
    {
      S_h_real_term[j] = 0.0;
      S_h_imag_term[j] = 0.0;
      for (int i = 0; i < eqSamplingTimes.size(); i++)
      {
        S_h_real[i][j] = 0.0;
        S_h_imag[i][j] = 0.0;
      }
    }
    /////////// STUFF FOR S(h) AND SISF ////////////////

    while (1)
    {
        if (snapshotC >= totalSnapshotsNum)
        {
            break;
        }

        // namestr = "data/sigmaSamples_"+to_string(bunchC)+".csv";
        // ifstream sigmaSamplesFile(namestr);

        // namestr = "data/xSamples_"+to_string(bunchC)+".csv";
        // ifstream xSamplesFile(namestr);

        // namestr = "data/ySamples_"+to_string(bunchC)+".csv";
        // ifstream ySamplesFile(namestr);

        // sigmaMatExtractor(sigmaMatBunch, sigmaSamplesFile);
        // sigmaMatExtractor(sitesXBunch, xSamplesFile);
        // sigmaMatExtractor(sitesYBunch, ySamplesFile);
        
        sigmaNamestr = "data/sigmaSamples_"+to_string(bunchC)+".csv";
        xNamestr = "data/xSamples_"+to_string(bunchC)+".csv";
        yNamestr = "data/ySamples_"+to_string(bunchC)+".csv";

        sigmaMatLoader(sigmaMatBunch, sigmaNamestr, NSites, samplesPerWrite);
        // doubleMatLoader(sitesXBunch, xNamestr, NSites, samplesPerWrite);
        // doubleMatLoader(sitesYBunch, yNamestr, NSites, samplesPerWrite);
        sigmaMatLoader(sitesXBunch, xNamestr, NSites, samplesPerWrite);
        sigmaMatLoader(sitesYBunch, yNamestr, NSites, samplesPerWrite);

        for(sampleC=0; sampleC<samplesPerWrite; sampleC++)
        {
            ////////// Making sigmaMat, sitesX, sitesY for this snapshot
            // for (int row=0; row<L; row++)
            // {
            //     for (int col=0; col<L; col++)
            //     {
            //         sigmaMat[row][col] = sigmaMatBunch[row][col][sampleC];
            //         sitesX[row][col]   = sitesXBunch[row][col][sampleC];
            //         sitesY[row][col]   = sitesYBunch[row][col][sampleC];
            //     }
            // }
            for (int siteC = 0; siteC < NSites; siteC++)
            {
                sigmaMat[siteC] = sigmaMatBunch[siteC][sampleC];
                sitesX[siteC]   = 1.0* Lx* sitesXBunch[siteC][sampleC] + siteComX[siteC] ;
                sitesY[siteC]   = 1.0* Ly* sitesYBunch[siteC][sampleC] + siteComY[siteC] ;
            }
            
            ////////// Making sigmaMat, sitesX, sitesY for this snapshot


            ////////// Making sample 1d vectors for this snapshot
            // configurationalCalcPerSample(areaSample, periSample, isoperiSample, xComSample, yComSample, sigmaMat, sitesX, sitesY);
            configurationalCalcPerSample_irr(areaSample, periSample, isoperiSample, ARSample, circSample, xComSample, yComSample, sigmaMat, sitesX, sitesY,\
                                     neighborNum, neighborsList, edges, latticeArea, siteIxx, siteIyy, siteIxy);
            ////////// Making sample 1d vectors for this snapshot
            totComFunc(xComTot, yComTot, xComSample, yComSample);

            //////////////// equilibrium distributions ////////////////////
            if (time[snapshotC] == eqSamplingTimes[eqSamplingTimesC])
            {
              ////// Writing the eq indices ////////
              if (eqSamplingTimesC == 0)
              {
                std::fstream fs;
                fs.open (ppDataFolderName+"/"+"eqSamplingTimes.txt", std::fstream::out);
                fs << "snapshotC"<<"\t\t"<<"time"<<"\t\t"<<"sampleC"<<"\t\t"<<"bunchC";
                fs.close();
              }
              std::fstream fs;
              fs.open (ppDataFolderName+"/"+"eqSamplingTimes.txt", std::fstream::app);
              fs <<"\n"<< snapshotC<<"\t\t"<<time[snapshotC]<<"\t\t"<<sampleC<<"\t\t"<<bunchC;
              fs.close();
              ////// Writing the eq indices ////////

              ////////// Writing pp of eq samples ////////////
              for (int cellC = 0; cellC <= NumCells; cellC++)
              {
                if (cellC==0)
                {
                  Area_dist_eq[cellC][eqSamplingTimesC] = 0.;
                  Peri_dist_eq[cellC][eqSamplingTimesC] = 0.;
                  isoperi_dist_eq[cellC][eqSamplingTimesC] = 0.;
                  AR_dist_eq[cellC][eqSamplingTimesC] = 0.;
                  circ_dist_eq[cellC][eqSamplingTimesC] = 0.;
                }
                else
                {
                  Area_dist_eq[cellC][eqSamplingTimesC] = areaSample[cellC];
                  Peri_dist_eq[cellC][eqSamplingTimesC] = periSample[cellC];
                  isoperi_dist_eq[cellC][eqSamplingTimesC] = isoperiSample[cellC];
                  AR_dist_eq[cellC][eqSamplingTimesC] = ARSample[cellC];
                  circ_dist_eq[cellC][eqSamplingTimesC] = circSample[cellC];
                }
              }
              ////////// Writing pp of eq samples ////////////


              ////////// Calc S(h), structure factor of this snapshot////////////

              S_h_calculator(xComSample, yComSample, Lx, Ly, h_vec, h_ang_vec, S_h_real_term, S_h_imag_term);
              S_h_real[eqSamplingTimesC] = S_h_real_term;
              S_h_imag[eqSamplingTimesC] = S_h_imag_term;
              ////////// Calc S(h), structure factor of this snapshot////////////



              eqSamplingTimesC++;
            }
            

            ////////// Mean //////////

            xComTotVec[snapshotC] = xComTot;
            yComTotVec[snapshotC] = yComTot;
            RComTotVec[snapshotC] = sqrt(pow(xComTot-xComTotVec[0],2)+pow(yComTot-yComTotVec[0],2));

            meanValArea = 0;
            meanValPeri = 0;
            meanValIsoperi = 0;
            meanValAR = 0;
            meanValCirc = 0;

            r2Avg = 0;
            r4Avg = 0;
            for (cellInd = 0; cellInd <= NumCells; cellInd++)
            {   
                if (cellInd > 0)
                {
                    meanValArea += ((double)areaSample[cellInd])/NumCells;
                    meanValPeri += ((double)periSample[cellInd])/NumCells;
                    meanValIsoperi += ((double)isoperiSample[cellInd])/NumCells;
                    meanValAR      += ((double)ARSample[cellInd])/NumCells;
                    meanValCirc    += ((double)circSample[cellInd])/NumCells;
                    r2Avg += (pow(xComSample[cellInd]-xComInit[cellInd],2)+pow(yComSample[cellInd]-yComInit[cellInd],2))/NumCells;
                    r4Avg += pow( (pow(xComSample[cellInd]-xComInit[cellInd],2)+pow(yComSample[cellInd]-yComInit[cellInd],2)) ,2)/NumCells;
                }
                xComBunch[cellInd][sampleC] = xComSample[cellInd];
                yComBunch[cellInd][sampleC] = yComSample[cellInd];
            }
            areaAvgStd[0][snapshotC] = meanValArea;
            periAvgStd[0][snapshotC] = meanValPeri;
            isoperiAvgStd[0][snapshotC] = meanValIsoperi;
            ARAvgStd[0][snapshotC] = meanValAR;
            circAvgStd[0][snapshotC] = meanValCirc;
            kurt[snapshotC] = ((0.5*r4Avg)/(pow(r2Avg,2)) - 1.);

            if (time[snapshotC]<t_w)
            {
                msdAvgStd[0][snapshotC] = 0;
            }else if (time[snapshotC]==t_w)
            {
                msdAvgStd[0][snapshotC] = 0;
                x_w[0]=0;
                y_w[0]=0;
                for (cellInd = 1; cellInd <= NumCells; cellInd++)
                {
                    x_w[cellInd]=xComSample[cellInd];
                    y_w[cellInd]=yComSample[cellInd];
                }
                doubleMatSaver(&x_w[0], NumCells+1, 1, ppDataFolderName+"/"+"x_w.csv");
                doubleMatSaver(&y_w[0], NumCells+1, 1, ppDataFolderName+"/"+"y_w.csv");

                xComTotW = xComTot;
                yComTotW = yComTot;

                std::ofstream xyComTotW;
                xyComTotW.open("xyComTotW.csv");
                xyComTotW << xComTotW;
                xyComTotW << "\n";
                xyComTotW << yComTotW;
                xyComTotW.close();

            }else if (time[snapshotC]>t_w)
            {
                msdAvgStd[0][snapshotC] = 0;
                for (cellInd = 1; cellInd <= NumCells; cellInd++)
                {
                    // R2Sample[cellInd] = (pow(xComSample[cellInd]-x_w[cellInd], 2) + pow(yComSample[cellInd]-y_w[cellInd], 2));
                    R2Sample[cellInd] = (pow((xComSample[cellInd]-xComTot)-(x_w[cellInd]-xComTotW), 2) + pow((yComSample[cellInd]-yComTot)-(y_w[cellInd]-yComTotW), 2));
                    msdAvgStd[0][snapshotC] += (R2Sample[cellInd] / NumCells);
                }
            }
            ////////// Mean //////////


            ////////// STD //////////
            stdValArea = 0;
            stdValPeri = 0;
            stdValIsoperi = 0;
            stdValAR = 0;
            stdValCirc = 0;

            for (cellInd = 0; cellInd <= NumCells; cellInd++)
            {   
                if (cellInd > 0)
                {
                    stdValArea += pow(areaSample[cellInd]-areaAvgStd[0][snapshotC], 2)/(NumCells);
                    stdValPeri += pow(periSample[cellInd]-periAvgStd[0][snapshotC], 2)/(NumCells);
                    stdValIsoperi += pow(isoperiSample[cellInd]-isoperiAvgStd[0][snapshotC], 2)/(NumCells);
                    stdValAR      += pow(ARSample[cellInd]-ARAvgStd[0][snapshotC], 2)/(NumCells);
                    stdValCirc    += pow(circSample[cellInd]-circAvgStd[0][snapshotC], 2)/(NumCells);
                }
            }
            stdValArea = sqrt(stdValArea);
            stdValPeri = sqrt(stdValPeri);
            stdValIsoperi = sqrt(stdValIsoperi);
            stdValAR = sqrt(stdValAR);
            stdValCirc = sqrt(stdValCirc);

            areaAvgStd[1][snapshotC] = stdValArea;
            periAvgStd[1][snapshotC] = stdValPeri;
            isoperiAvgStd[1][snapshotC] = stdValIsoperi;
            ARAvgStd[1][snapshotC] = stdValAR;
            circAvgStd[1][snapshotC] = stdValCirc;

            if (time[snapshotC]<t_w)
            {
                msdAvgStd[1][snapshotC] = 0;
            }else if (time[snapshotC]==t_w)
            {
                msdAvgStd[1][snapshotC] = 0;
            }else if (time[snapshotC]>t_w)
            {
                msdAvgStd[1][snapshotC] = 0;
                for (cellInd = 1; cellInd <= NumCells; cellInd++)
                {
                    msdAvgStd[1][snapshotC] += pow(R2Sample[cellInd]-msdAvgStd[0][snapshotC], 2) / (NumCells-1);
                }
                msdAvgStd[1][snapshotC] = sqrt(msdAvgStd[1][snapshotC]);
            }
            ////////// STD //////////
            snapshotC++;
        }
        // doubleMatSaver(&xComBunch[0][0], NumCells+1, samplesPerWrite, ppDataFolderName+"/"+"xComBunch_"+to_string(bunchC)+".csv");
        // doubleMatSaver(&yComBunch[0][0], NumCells+1, samplesPerWrite, ppDataFolderName+"/"+"yComBunch_"+to_string(bunchC)+".csv");
        double2DVecSaver(xComBunch, ppDataFolderName+"/"+"xComBunch_"+to_string(bunchC)+".csv");
        double2DVecSaver(yComBunch, ppDataFolderName+"/"+"yComBunch_"+to_string(bunchC)+".csv");
        bunchC++;
    }
    bunchC--;
    snapshotC--;

    // doubleMatSaver(&areaAvgStd[0][0], 2, totalSnapshotsNum, ppDataFolderName+"/"+"areaAvgStd.csv");
    // doubleMatSaver(&periAvgStd[0][0], 2, totalSnapshotsNum, ppDataFolderName+"/"+"periAvgStd.csv");
    // doubleMatSaver(&isoperiAvgStd[0][0], 2, totalSnapshotsNum, ppDataFolderName+"/"+"isoperiAvgStd.csv");
    // doubleMatSaver(&msdAvgStd[0][0], 2, totalSnapshotsNum, ppDataFolderName+"/"+"msdAvgStd.csv");
    double2DVecSaver(areaAvgStd, ppDataFolderName+"/"+"areaAvgStd.csv");
    double2DVecSaver(periAvgStd, ppDataFolderName+"/"+"periAvgStd.csv");
    double2DVecSaver(isoperiAvgStd, ppDataFolderName+"/"+"isoperiAvgStd.csv");
    double2DVecSaver(ARAvgStd, ppDataFolderName+"/"+"ARAvgStd.csv");
    double2DVecSaver(circAvgStd, ppDataFolderName+"/"+"circAvgStd.csv");
    double2DVecSaver(msdAvgStd, ppDataFolderName+"/"+"msdAvgStd.csv");
    double1DVecSaver(kurt, ppDataFolderName+"/"+"kurt.csv");
    double1DVecSaver(xComTotVec, ppDataFolderName+"/"+"xComTotVec.csv");
    double1DVecSaver(yComTotVec, ppDataFolderName+"/"+"yComTotVec.csv");
    double1DVecSaver(RComTotVec, ppDataFolderName+"/"+"RComTotVec.csv");

    /// equilibrium sampling distributions
    double2DVecSaver(Area_dist_eq, ppDataFolderName+"/"+"Area_dist_eq.csv");
    double2DVecSaver(Peri_dist_eq, ppDataFolderName+"/"+"Peri_dist_eq.csv");
    double2DVecSaver(isoperi_dist_eq, ppDataFolderName+"/"+"isoperi_dist_eq.csv");
    double2DVecSaver(AR_dist_eq, ppDataFolderName+"/"+"AR_dist_eq.csv");
    double2DVecSaver(circ_dist_eq, ppDataFolderName+"/"+"circ_dist_eq.csv");
    



    // // configurationalCalcPerBunch(&areaAvgStd[0][0], &periAvgStd[0][0], &isoperiAvgStd[0][0], totalSnapshotsNum);



    // /////////////////area, peri, isoperi////////////////

    // /////////////////xCom, yCom, MSD////////////////
    // /////////////////xCom, yCom, MSD////////////////






    /////////////// Python part of pp ////////////////
    // system("python3 pp_irr_v16.py");
    // system("python3 pp_py_gen_v0.py");
    /////////////// Python part of pp ////////////////
    
    ///// Re-deleting the data folder /////
    if (existFunc("data.zip"))
      {dataFolderDeleter();}
    ///// Re-deleting the data folder /////
    
    
    return 0;
}

//////////////////////////////////////////////////////////////////////////
////////////////////////////// FUNCTIONS /////////////////////////////////
//////////////////////////////////////////////////////////////////////////
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
    std::stringstream ss(line);
    double value;
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
const std::vector<double>& latticeArea,\
const std::vector<double>& siteIxx, const std::vector<double>& siteIyy, const std::vector<double>& siteIxy)
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

        IxxSample[cellIndex] += pow((sitesX[siteC]-xComSample[cellIndex]),2) * latticeArea[siteC] + siteIxx[siteC];
        IyySample[cellIndex] += pow((sitesY[siteC]-yComSample[cellIndex]),2) * latticeArea[siteC] + siteIyy[siteC];
        IxySample[cellIndex] += (sitesX[siteC]-xComSample[cellIndex])*(sitesY[siteC]-yComSample[cellIndex]) * latticeArea[siteC]  + siteIxy[siteC];

        // if (cellIndex==2 && (fabs(pow((sitesX[siteC]-xComSample[cellIndex]),2) * latticeArea[siteC] + siteIxx[siteC])>100 || fabs(pow((sitesY[siteC]-yComSample[cellIndex]),2) * latticeArea[siteC] + siteIyy[siteC])>100 ) )
        // {
        //   int r;
        //   r=5;
        //   cout<<siteC<<endl;
        //   cout<<sitesX[siteC]<<endl;
        //   cout<<sitesY[siteC]<<endl;
        //   cout<<xComSample[cellIndex]<<endl;
        //   cout<<yComSample[cellIndex]<<endl;
        //   cout<<latticeArea[siteC]<<endl;
        //   cout<<siteIxx[siteC]<<endl;
        //   cout<<siteIyy[siteC]<<endl;
        //   cout<<siteIxy[siteC]<<endl;
        // }
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
//
void S_h_calculator(const vector<double>& xComSample, const vector<double>& yComSample, const double Lx, const double Ly, const vector<double>& h_abs_vec, const vector<double>&  h_ang_vec, \
                    vector<double>& S_h_real_term, vector<double>& S_h_imag_term)
{
    int NumCells = xComSample.size()-1;
    int cellC_1, cellC_2, cellC;

    vector<double> xComSampleFolded(NumCells+1);
    vector<double> yComSampleFolded(NumCells+1);

    for (cellC = 1; cellC <= NumCells; cellC++)
    {
      /* code */
    }
    

    int N_h = h_abs_vec.size();
    int N_ang = h_ang_vec.size();

    double h_abs, h_ang;

    double S_h_real_term_ang, S_h_imag_term_ang;

    double h_vec_x, h_vec_y;
    
    vector<vector<double>> deltaXMat(NumCells+1, vector<double>(NumCells+1));
    vector<vector<double>> deltaYMat(NumCells+1, vector<double>(NumCells+1));

    for ( cellC_1 = 1; cellC_1 <= NumCells; cellC_1++)
    {
      for (cellC_2 = cellC_1; cellC_2 <= NumCells; cellC_2++)
      {
        
      }
    }
    

    for (int hC = 0; hC < N_h; hC++)
    {

      h_abs = h_abs_vec[hC];

      S_h_real_term[hC] = 0.0;
      S_h_imag_term[hC] = 0.0;

      for (int h_ang_C = 0; h_ang_C < N_ang; h_ang_C++)
      {
        
        h_ang = h_ang_vec[h_ang_C];
        
        h_vec_x = h_abs * cos(h_ang);
        h_vec_y = h_abs * sin(h_ang);

        S_h_real_term_ang = 0.0;
        S_h_imag_term_ang = 0.0;





        
        S_h_real_term[hC] += S_h_real_term_ang/N_ang;
        S_h_imag_term[hC] += S_h_imag_term_ang/N_ang;
      }

      
    }
    

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

//////////////////////////////////////////////////////////////////////////
////////////////////////////// FUNCTIONS /////////////////////////////////
//////////////////////////////////////////////////////////////////////////
