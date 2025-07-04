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
#include <map>
#include <set>
#include <utility> // for std::pair


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

void newMSD();

void newMSD_non_corr();

//////////////////////////////////////////////////////////////////////////
//////////////////////////// PROTOTYPES //////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main()
{
    newMSD_non_corr();
    return 0;
}

//////////////////////////////////////////////////////////////////////////
////////////////////////////// FUNCTIONS /////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void newMSD()
{
    struct interval {
    int basic_dt;
    int begin_time;
    int begin_ind;
    int end_time;
    int end_ind;
    int num_steps;
    };

    //////////////////// time read /////////////////
    int t_w;
    ifstream t_w_File;
    t_w_File.open("t_w.csv");
    t_w_File >> t_w;
    t_w_File.close();

    std::vector<int> time;
    std::ifstream file("pp_data/time.csv");
    if (!file.is_open()) {
        std::cerr << "Error opening file time.csv" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;

        while (std::getline(ss, item, ',')) {
            try {
                time.push_back(std::stoi(item));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number in file: " << item << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Number out of range in file: " << item << std::endl;
            }
        }
    }
    file.close();

    int t_w_ind=0;
    while(time[t_w_ind] != t_w)
    {
        t_w_ind++;
    }
    if (t_w_ind == time.size())
    {
        std::cerr << "tw is not in time.csv" << std::endl;
        return;
    }
    //////////////////// time read /////////////////

    //////////////////// smpling pattern read //////////////////
    std::vector<int> sampling_steps;
    std::vector<int> sampling_milestones;
    file.open ("samplingPattern_vec.csv", std::ifstream::in);
    if (!file.is_open()) {
        std::cerr << "Error opening file samplingPattern_vec.csv" << std::endl;
        return;
    }

    int lineNumber = 0;
    while (std::getline(file, line)) {
        std::vector<int> values;
        values.clear();
        std::stringstream ss(line);
        std::string item;

        while (std::getline(ss, item, ',')) {
            try {
                values.push_back(std::stoi(item));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number in file: " << item << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Number out of range in file: " << item << std::endl;
            }
        }

        for (int i = 0; i < values.size(); i++)
        {
            if (lineNumber == 0) {
                sampling_milestones.push_back(values[i]);
            } else if (lineNumber == 1) {
                sampling_steps.push_back(values[i]);
            }
        }

        lineNumber++;
    }
    file.close();

    // correction of milestones, by  maxMCSteps
    int milestone_c = sampling_milestones.size() - 1;
    int elimination_number = -1; // it sould not be zero
    while(1)
    {
        if (sampling_milestones[milestone_c] > time[time.size()-1])
        {
            sampling_milestones[milestone_c] = time[time.size()-1];
            elimination_number++;
        }
        else
        {break;}
        milestone_c--;
    }
    for (int elim_c = 0; elim_c < elimination_number; elim_c++)
    {
        sampling_milestones.pop_back();
        sampling_steps.pop_back();
    }

    milestone_c = 0;
    elimination_number = 0; // it sould be zero
    while(1)
    {
        if (sampling_milestones[milestone_c] <= t_w)
        {
            sampling_milestones[milestone_c] = t_w;
            elimination_number++;
        }
        else
        {break;}
        milestone_c++;
    }
    for (int elim_c = 0; elim_c < elimination_number; elim_c++)
    {
        sampling_milestones.erase(sampling_milestones.begin());
        sampling_steps.erase(sampling_steps.begin());
    }
    // correction of milestones, by maxMCSteps
    

    std::map<int, std::vector<int>>  sampling_intervals;
    std::map<int, std::vector<int>>  sampling_intervals_ind;
    std::vector<struct interval>  intvl_vec;
    struct interval intvl;
    // intvl.clear();
    
    for (int i = 0; i < sampling_steps.size(); i++)
    {
        // if (t_w >= sampling_milestones[i])
        // {continue;}
        
        sampling_intervals[sampling_steps[i]] = {0, 0};

        if (i == 0)
        {
            // sampling_intervals[sampling_steps[i]][0] = 0 ;
            sampling_intervals[sampling_steps[i]][0] = t_w ;
        } else
        {
            sampling_intervals[sampling_steps[i]][0] =  sampling_milestones[i-1] ;
        }

        if (i == sampling_steps.size() - 1)
        {
            sampling_intervals[sampling_steps[i]][1] =  time[time.size()-1] ;
        }
        else
        {
            sampling_intervals[sampling_steps[i]][1] =  sampling_milestones[i] ;
        }

        // if (i == sampling_steps.size()-1)
        // {
        //     int b = 6;
        // }

        int ind = 0;
        while ( time[ind] !=  sampling_intervals[sampling_steps[i]][1] )
        {
            if ( time[ind] == sampling_intervals[sampling_steps[i]][0] )
            {
                sampling_intervals_ind[sampling_steps[i]].push_back(ind);
            }
            ind++;
        }
        sampling_intervals_ind[sampling_steps[i]].push_back(ind);

        // intvl.clear();
        intvl.basic_dt    = sampling_steps[i];
        intvl.begin_time  = sampling_intervals[sampling_steps[i]][0];
        intvl.begin_ind   =  sampling_intervals_ind[sampling_steps[i]][0];
        intvl.end_time    =  sampling_intervals[sampling_steps[i]][1];
        intvl.end_ind     =  sampling_intervals_ind[sampling_steps[i]][1];
        intvl.num_steps = intvl.end_ind - intvl.begin_ind + 1;

        intvl_vec.push_back(intvl);
        // continue;
    }

    //////////////////// smpling pattern read //////////////////

    //////////////////xCom, yCom read////////////////////////
    std::vector<std::vector<double>>  xCom;
    std::vector<std::vector<double>>  yCom;

    /// init
    std::ifstream xfile("pp_data/xComInit.csv");
    std::ifstream yfile("pp_data/yComInit.csv");

    if (!xfile.is_open()) { // Check if the file is open
        std::cerr << "Error: Could not open the file xComInit" << std::endl;
        return;
    }
    if (!yfile.is_open()) { // Check if the file is open
        std::cerr << "Error: Could not open the file yComInit" << std::endl;
        return;
    }

    
    while (std::getline(xfile, line)) { // Read each line from the file
        try {
            double value = std::stod(line); // Convert the line to a double
            xCom.push_back({value}); // Add the value to the vector
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid number format in file." << std::endl;
            return;
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: Number out of range in file." << std::endl;
            return;
        }
    }
    while (std::getline(yfile, line)) { // Read each line from the file
        try {
            double value = std::stod(line); // Convert the line to a double
            yCom.push_back({value}); // Add the value to the vector
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid number format in file." << std::endl;
            return;
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: Number out of range in file." << std::endl;
            return;
        }
    }
    yfile.close(); // Close the file
    xfile.close(); // Close the file
    /// init

    /// bunches
    int bunchC = 0;
    int NumCells = xCom.size()-1;
    string Xnamestr, Ynamestr;
    int cellInd;

    while (1)
    {
        Xnamestr = "pp_data/xComBunch_"+to_string(bunchC)+".csv";
        Ynamestr = "pp_data/yComBunch_"+to_string(bunchC)+".csv";
        xfile.open (Xnamestr , std::ifstream::in);
        yfile.open (Ynamestr , std::ifstream::in);

        if ( (!xfile.is_open()) || (!yfile.is_open()) )
        {
            break; // if the bunches are over
        }

        cellInd = 0;
        while (std::getline(xfile, line))
        {
            std::vector<double> xCom_row;
            std::stringstream xCom_ss(line);
            std::string xCom_value;
            
            while (std::getline(xCom_ss, xCom_value, ',')) {
                try {
                    double xCom_number = std::stod(xCom_value);
                    xCom[cellInd].push_back(xCom_number);
                    // xCom_row.push_back(xCom_number);
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid number format: " << xCom_value << std::endl;
                } catch (const std::out_of_range& e) {
                    std::cerr << "Number out of range: " << xCom_value << std::endl;
                }
            }
            cellInd++;
        }
        cellInd = 0;
        while (std::getline(yfile, line))
        {
            std::vector<double> yCom_row;
            std::stringstream yCom_ss(line);
            std::string yCom_value;
            
            while (std::getline(yCom_ss, yCom_value, ',')) {
                try {
                    double yCom_number = std::stod(yCom_value);
                    yCom[cellInd].push_back(yCom_number);
                    // yCom_row.push_back(yCom_number);
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid number format: " << yCom_value << std::endl;
                } catch (const std::out_of_range& e) {
                    std::cerr << "Number out of range: " << yCom_value << std::endl;
                }
            }
            cellInd++;
            // yCom.push_back(yCom_row);
        }
        
        xfile.close();
        yfile.close();

        bunchC++;
    }
    /// bunches

    //////////////////xCom, yCom read////////////////////////


    //////////////////// MSD calc //////////////////
    std::set<int> MSD_data_points_keys; // These will be the keys for a map that stores all the datapoints
    // t_minus_t_w.insert(10);
    std::map<int, std::vector<double>>  MSD_data_points;

    int n;
    int dt_base;
    int dt;
    int intvl_total_DT;
    int dt_beg_ind, dt_end_ind;
    double avg_MSD;
    double x_term, y_term, sum;
    double xComTot_beg, yComTot_beg;
    double xComTot_end, yComTot_end;
    
    // within each sampling interval
    for (int intvl_c = 0; intvl_c < intvl_vec.size(); intvl_c++)
    {
        dt_base = intvl_vec[intvl_c].basic_dt;
        intvl_total_DT = intvl_vec[intvl_c].end_time - intvl_vec[intvl_c].begin_time;
        n = 1; // dt_base is multiplied by this
        dt = n * dt_base;
        

        while (dt <= intvl_total_DT)
        {
            // if (dt > intvl_total_DT/2.0+(1e-6))
            // {
            //     dt = intvl_total_DT;
            //     dt_beg_ind = intvl_vec[intvl_c].begin_ind;
            //     dt_end_ind = intvl_vec[intvl_c].end_ind;
            // }
            // else
            // {
            //     dt_beg_ind = intvl_vec[intvl_c].begin_ind;
            //     dt_end_ind = dt_beg_ind + n;
            // }

            dt_beg_ind = intvl_vec[intvl_c].begin_ind;
            dt_end_ind = dt_beg_ind + n;
            
            MSD_data_points_keys.insert(dt);

            while (dt_end_ind <= intvl_vec[intvl_c].end_ind)
            {

                xComTot_beg = 0.0;
                yComTot_beg = 0.0;

                xComTot_end = 0.0;
                yComTot_end = 0.0;
                for (cellInd = 1; cellInd <= NumCells; cellInd++)
                {
                    xComTot_beg += xCom[cellInd][dt_beg_ind]/NumCells;
                    yComTot_beg += yCom[cellInd][dt_beg_ind]/NumCells;

                    xComTot_end += xCom[cellInd][dt_end_ind]/NumCells;
                    yComTot_end += yCom[cellInd][dt_end_ind]/NumCells;
                }

                sum = 0.0;
                for (cellInd = 1; cellInd <= NumCells; cellInd++)
                {
                    x_term = (xCom[cellInd][dt_end_ind] - xComTot_end) - (xCom[cellInd][dt_beg_ind] - xComTot_beg);
                    y_term = (yCom[cellInd][dt_end_ind] - yComTot_end) - (yCom[cellInd][dt_beg_ind] - yComTot_beg);
                    sum = sum + pow( x_term , 2) + pow( y_term , 2);
                    // sum = sum + pow(( ( xCom[cellInd][dt_end_ind] ) - ( xCom[cellInd][dt_beg_ind]) ) , 2) + pow(( ( yCom[cellInd][dt_beg_ind] ) - ( yCom[cellInd][dt_end_ind] ) ), 2);
                }
                avg_MSD = sum / NumCells;

                MSD_data_points[dt].push_back(avg_MSD);

                dt_beg_ind += n;
                dt_end_ind += n;
            }
            
            n++;
            // n *= 2;
            dt = n * dt_base;
        }
    }

    // across multiple sampling intervals
    for (int intvl_c_1 = 0; intvl_c_1 < intvl_vec.size(); intvl_c_1++)
    {
        for (int intvl_c_2 = intvl_c_1 + 1 ; intvl_c_2 < intvl_vec.size(); intvl_c_2++)
        {
            dt = intvl_vec[intvl_c_2].end_time - intvl_vec[intvl_c_1].begin_time;
            MSD_data_points_keys.insert(dt);

            dt_beg_ind = intvl_vec[intvl_c_1].begin_ind;
            dt_end_ind = intvl_vec[intvl_c_2].end_ind;

            xComTot_beg = 0.0;
            yComTot_beg = 0.0;

            xComTot_end = 0.0;
            yComTot_end = 0.0;
            for (cellInd = 1; cellInd <= NumCells; cellInd++)
            {
                xComTot_beg += xCom[cellInd][dt_beg_ind]/NumCells;
                yComTot_beg += yCom[cellInd][dt_beg_ind]/NumCells;

                xComTot_end += xCom[cellInd][dt_end_ind]/NumCells;
                yComTot_end += yCom[cellInd][dt_end_ind]/NumCells;
            }

            sum = 0.0;
            for (cellInd = 1; cellInd <= NumCells; cellInd++)
            {
                x_term = (xCom[cellInd][dt_end_ind] - xComTot_end) - (xCom[cellInd][dt_beg_ind] - xComTot_beg);
                y_term = (yCom[cellInd][dt_end_ind] - yComTot_end) - (yCom[cellInd][dt_beg_ind] - yComTot_beg);
                sum = sum + pow( x_term , 2) + pow( y_term , 2);
                // sum = sum + pow(( ( xCom[cellInd][dt_end_ind] ) - ( xCom[cellInd][dt_beg_ind]) ) , 2) + pow(( ( yCom[cellInd][dt_beg_ind] ) - ( yCom[cellInd][dt_end_ind] ) ), 2);
            }
            avg_MSD = sum / NumCells;

            MSD_data_points[dt].push_back(avg_MSD);
        }
    }
    //////////////////// MSD calc //////////////////

    //////////////////// MSD writing in file //////////////////
    std::vector<int> MSD_data_points_keys_vec(MSD_data_points_keys.begin(), MSD_data_points_keys.end());
    std::sort(MSD_data_points_keys_vec.begin(), MSD_data_points_keys_vec.end());
    
    std::ofstream MSD_file_write;
    MSD_file_write.open("pp_data/newMSD.csv");
    for (int dt_c = 0; dt_c < MSD_data_points_keys_vec.size(); dt_c++)
    {
        dt = MSD_data_points_keys_vec[dt_c];

        MSD_file_write << dt;
        MSD_file_write << ": ";

        for (int datapoint_c = 0; datapoint_c < MSD_data_points[dt].size(); datapoint_c++)
        {
            MSD_file_write << MSD_data_points[dt][datapoint_c] << std::setprecision(8);
            if (datapoint_c < MSD_data_points[dt].size() - 1)
            {
                MSD_file_write << ", ";
            }
        }
        
        if (dt_c < MSD_data_points_keys_vec.size() - 1)
        {
            MSD_file_write << "\n";
        }
    }
    MSD_file_write.close();
    //////////////////// MSD writing in file //////////////////
}

void newMSD_non_corr()
{
    struct interval {
    int basic_dt;
    int begin_time;
    int begin_ind;
    int end_time;
    int end_ind;
    int num_steps;
    };

    //////////////////// time read /////////////////
    int t_w;
    ifstream t_w_File;
    t_w_File.open("t_w.csv");
    t_w_File >> t_w;
    t_w_File.close();

    std::vector<int> time;
    std::ifstream file("pp_data/time.csv");
    if (!file.is_open()) {
        std::cerr << "Error opening file time.csv" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;

        while (std::getline(ss, item, ',')) {
            try {
                time.push_back(std::stoi(item));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number in file: " << item << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Number out of range in file: " << item << std::endl;
            }
        }
    }
    file.close();

    int t_w_ind=0;
    while(time[t_w_ind] != t_w)
    {
        t_w_ind++;
    }
    if (t_w_ind == time.size())
    {
        std::cerr << "tw is not in time.csv" << std::endl;
        return;
    }
    //////////////////// time read /////////////////

    //////////////////// smpling pattern read //////////////////
    std::vector<int> sampling_steps;
    std::vector<int> sampling_milestones;
    file.open ("samplingPattern_vec.csv", std::ifstream::in);
    if (!file.is_open()) {
        std::cerr << "Error opening file samplingPattern_vec.csv" << std::endl;
        return;
    }

    int lineNumber = 0;
    while (std::getline(file, line)) {
        std::vector<int> values;
        values.clear();
        std::stringstream ss(line);
        std::string item;

        while (std::getline(ss, item, ',')) {
            try {
                values.push_back(std::stoi(item));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number in file: " << item << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Number out of range in file: " << item << std::endl;
            }
        }

        for (int i = 0; i < values.size(); i++)
        {
            if (lineNumber == 0) {
                sampling_milestones.push_back(values[i]);
            } else if (lineNumber == 1) {
                sampling_steps.push_back(values[i]);
            }
        }

        lineNumber++;
    }
    file.close();

    // correction of milestones, by  maxMCSteps
    int milestone_c = sampling_milestones.size() - 1;
    int elimination_number = -1; // it sould not be zero
    while(1)
    {
        if (sampling_milestones[milestone_c] > time[time.size()-1])
        {
            sampling_milestones[milestone_c] = time[time.size()-1];
            elimination_number++;
        }
        else
        {break;}
        milestone_c--;
    }
    for (int elim_c = 0; elim_c < elimination_number; elim_c++)
    {
        sampling_milestones.pop_back();
        sampling_steps.pop_back();
    }

    milestone_c = 0;
    elimination_number = 0; // it sould be zero
    while(1)
    {
        if (sampling_milestones[milestone_c] <= t_w)
        {
            sampling_milestones[milestone_c] = t_w;
            elimination_number++;
        }
        else
        {break;}
        milestone_c++;
    }
    for (int elim_c = 0; elim_c < elimination_number; elim_c++)
    {
        sampling_milestones.erase(sampling_milestones.begin());
        sampling_steps.erase(sampling_steps.begin());
    }
    // correction of milestones, by maxMCSteps
    

    std::map<int, std::vector<int>>  sampling_intervals;
    std::map<int, std::vector<int>>  sampling_intervals_ind;
    std::vector<struct interval>  intvl_vec;
    struct interval intvl;
    // intvl.clear();
    
    for (int i = 0; i < sampling_steps.size(); i++)
    {
        // if (t_w >= sampling_milestones[i])
        // {continue;}
        
        sampling_intervals[sampling_steps[i]] = {0, 0};

        if (i == 0)
        {
            // sampling_intervals[sampling_steps[i]][0] = 0 ;
            sampling_intervals[sampling_steps[i]][0] = t_w ;
        } else
        {
            sampling_intervals[sampling_steps[i]][0] =  sampling_milestones[i-1] ;
        }

        if (i == sampling_steps.size() - 1)
        {
            sampling_intervals[sampling_steps[i]][1] =  time[time.size()-1] ;
        }
        else
        {
            sampling_intervals[sampling_steps[i]][1] =  sampling_milestones[i] ;
        }

        // if (i == sampling_steps.size()-1)
        // {
        //     int b = 6;
        // }

        int ind = 0;
        while ( time[ind] !=  sampling_intervals[sampling_steps[i]][1] )
        {
            if ( time[ind] == sampling_intervals[sampling_steps[i]][0] )
            {
                sampling_intervals_ind[sampling_steps[i]].push_back(ind);
            }
            ind++;
        }
        sampling_intervals_ind[sampling_steps[i]].push_back(ind);

        // intvl.clear();
        intvl.basic_dt    = sampling_steps[i];
        intvl.begin_time  = sampling_intervals[sampling_steps[i]][0];
        intvl.begin_ind   =  sampling_intervals_ind[sampling_steps[i]][0];
        intvl.end_time    =  sampling_intervals[sampling_steps[i]][1];
        intvl.end_ind     =  sampling_intervals_ind[sampling_steps[i]][1];
        intvl.num_steps = intvl.end_ind - intvl.begin_ind + 1;

        intvl_vec.push_back(intvl);
        // continue;
    }

    //////////////////// smpling pattern read //////////////////

    //////////////////xCom, yCom read////////////////////////
    std::vector<std::vector<double>>  xCom;
    std::vector<std::vector<double>>  yCom;

    /// init
    std::ifstream xfile("pp_data/xComInit.csv");
    std::ifstream yfile("pp_data/yComInit.csv");

    if (!xfile.is_open()) { // Check if the file is open
        std::cerr << "Error: Could not open the file xComInit" << std::endl;
        return;
    }
    if (!yfile.is_open()) { // Check if the file is open
        std::cerr << "Error: Could not open the file yComInit" << std::endl;
        return;
    }

    
    while (std::getline(xfile, line)) { // Read each line from the file
        try {
            double value = std::stod(line); // Convert the line to a double
            xCom.push_back({value}); // Add the value to the vector
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid number format in file." << std::endl;
            return;
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: Number out of range in file." << std::endl;
            return;
        }
    }
    while (std::getline(yfile, line)) { // Read each line from the file
        try {
            double value = std::stod(line); // Convert the line to a double
            yCom.push_back({value}); // Add the value to the vector
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid number format in file." << std::endl;
            return;
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: Number out of range in file." << std::endl;
            return;
        }
    }
    yfile.close(); // Close the file
    xfile.close(); // Close the file
    /// init

    /// bunches
    int bunchC = 0;
    int NumCells = xCom.size()-1;
    string Xnamestr, Ynamestr;
    int cellInd;

    while (1)
    {
        Xnamestr = "pp_data/xComBunch_"+to_string(bunchC)+".csv";
        Ynamestr = "pp_data/yComBunch_"+to_string(bunchC)+".csv";
        xfile.open (Xnamestr , std::ifstream::in);
        yfile.open (Ynamestr , std::ifstream::in);

        if ( (!xfile.is_open()) || (!yfile.is_open()) )
        {
            break; // if the bunches are over
        }

        cellInd = 0;
        while (std::getline(xfile, line))
        {
            std::vector<double> xCom_row;
            std::stringstream xCom_ss(line);
            std::string xCom_value;
            
            while (std::getline(xCom_ss, xCom_value, ',')) {
                try {
                    double xCom_number = std::stod(xCom_value);
                    xCom[cellInd].push_back(xCom_number);
                    // xCom_row.push_back(xCom_number);
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid number format: " << xCom_value << std::endl;
                } catch (const std::out_of_range& e) {
                    std::cerr << "Number out of range: " << xCom_value << std::endl;
                }
            }
            cellInd++;
        }
        cellInd = 0;
        while (std::getline(yfile, line))
        {
            std::vector<double> yCom_row;
            std::stringstream yCom_ss(line);
            std::string yCom_value;
            
            while (std::getline(yCom_ss, yCom_value, ',')) {
                try {
                    double yCom_number = std::stod(yCom_value);
                    yCom[cellInd].push_back(yCom_number);
                    // yCom_row.push_back(yCom_number);
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid number format: " << yCom_value << std::endl;
                } catch (const std::out_of_range& e) {
                    std::cerr << "Number out of range: " << yCom_value << std::endl;
                }
            }
            cellInd++;
            // yCom.push_back(yCom_row);
        }
        
        xfile.close();
        yfile.close();

        bunchC++;
    }
    /// bunches

    //////////////////xCom, yCom read////////////////////////


    //////////////////// MSD calc //////////////////
    std::set<int> MSD_data_points_keys; // These will be the keys for a map that stores all the datapoints
    // t_minus_t_w.insert(10);
    std::map<int, std::vector<double>>  MSD_data_points;

    int n;
    int dt_base;
    int dt;
    int intvl_total_DT;
    int dt_beg_ind, dt_end_ind;
    double avg_MSD;
    double x_term, y_term, sum;
    double xComTot_beg, yComTot_beg;
    double xComTot_end, yComTot_end;
    
    // simple MSD calc
    dt_beg_ind = intvl_vec[0].begin_ind;
    xComTot_beg = 0.0;
    yComTot_beg = 0.0;
    for (cellInd = 1; cellInd <= NumCells; cellInd++)
    {
        xComTot_beg += xCom[cellInd][dt_beg_ind]/NumCells;
        yComTot_beg += yCom[cellInd][dt_beg_ind]/NumCells;
    }

    for (int sample_c = dt_beg_ind+1; sample_c < intvl_vec[intvl_vec.size()-1].end_ind; sample_c++)
    {
        dt_end_ind = sample_c;
        dt = time[dt_end_ind] - time[dt_beg_ind];
        MSD_data_points_keys.insert(dt);

        xComTot_end = 0.0;
        yComTot_end = 0.0;
        for (cellInd = 1; cellInd <= NumCells; cellInd++)
        {
            xComTot_end += xCom[cellInd][dt_end_ind]/NumCells;
            yComTot_end += yCom[cellInd][dt_end_ind]/NumCells;
        }

        sum = 0.0;
        for (cellInd = 1; cellInd <= NumCells; cellInd++)
        {
            x_term = (xCom[cellInd][dt_end_ind] - xComTot_end) - (xCom[cellInd][dt_beg_ind] - xComTot_beg);
            y_term = (yCom[cellInd][dt_end_ind] - yComTot_end) - (yCom[cellInd][dt_beg_ind] - yComTot_beg);
            sum = sum + pow( x_term , 2) + pow( y_term , 2);
            // sum = sum + pow(( ( xCom[cellInd][dt_end_ind] ) - ( xCom[cellInd][dt_beg_ind]) ) , 2) + pow(( ( yCom[cellInd][dt_beg_ind] ) - ( yCom[cellInd][dt_end_ind] ) ), 2);
        }
        avg_MSD = sum / NumCells;

        MSD_data_points[dt].push_back(avg_MSD);
    }
    // simple MSD calc

    //////////////////// MSD writing in file //////////////////
    std::vector<int> MSD_data_points_keys_vec(MSD_data_points_keys.begin(), MSD_data_points_keys.end());
    std::sort(MSD_data_points_keys_vec.begin(), MSD_data_points_keys_vec.end());
    
    std::ofstream MSD_file_write;
    MSD_file_write.open("pp_data/newMSD.csv");
    for (int dt_c = 0; dt_c < MSD_data_points_keys_vec.size(); dt_c++)
    {
        dt = MSD_data_points_keys_vec[dt_c];

        MSD_file_write << dt;
        MSD_file_write << ": ";

        for (int datapoint_c = 0; datapoint_c < MSD_data_points[dt].size(); datapoint_c++)
        {
            MSD_file_write << MSD_data_points[dt][datapoint_c] << std::setprecision(8);
            if (datapoint_c < MSD_data_points[dt].size() - 1)
            {
                MSD_file_write << ", ";
            }
        }
        
        if (dt_c < MSD_data_points_keys_vec.size() - 1)
        {
            MSD_file_write << "\n";
        }
    }
    MSD_file_write.close();
    //////////////////// MSD writing in file //////////////////


}

//////////////////////////////////////////////////////////////////////////
////////////////////////////// FUNCTIONS /////////////////////////////////
//////////////////////////////////////////////////////////////////////////