- This folder contains an individual run for the parameters specified in the file "simulationData_vec.csv".
- The file "simulationData_vec.csv" contains some lattice properties (Lx, Ly, NSites) , Alpha,
total simulation time (maxMCSteps) and some other data that are easy to understand.
- The file "t_w.csv" contains waiting time.
- To compile all the codes in this folder, run the file "compile.sh"
- To run the simulation, run the file "run.bat".
- When the simulation is done (when the folder data is automatically deleted), run 
the file "pp_individual_run.bat" for post process. The post-processed data will be saved in
folder "pp_data", and the plots will be saved in the folder "pp_plot". You can turn the plot switches
off (on zero) in the file "plot_switches.txt".
There are other steps of post-process in higher folders. Read their own readme files.
