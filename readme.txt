- Folders for different runs for each alpha value.
- The name of the folders must be like "run_1", "run_2" and so on.
- For "hex" and "irr", there must be a "lattice" folder which includes the data of the related lattice.
- There must be a "simulationData_Vec.csv" file here with the right amount of alpha. You can simply copy it from the folder "run_1".
- Once all the individual runs in this folder (and their individual post-process) are done, the MSD, beta, and D_eff can be calculated by running "MSD_calc.sh", and are saved in a folder called "msd_data". Corrent the number of folders in this file.
- After MSD calculation, one needs to do a post-process for all of the runs for this alpha.
To do so, first, run the file "compile.sh". You need to correct it for the right address of voro++ package in your system.
Once the file "pp_vor_v4.cpp" is compiled, run the file "pp_each_alpha.sh".


