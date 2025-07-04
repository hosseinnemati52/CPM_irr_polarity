#!/bin/bash
declare -i N=2


for i in $(seq 1 $N); do cd run_${i}; ./MSD_2x; cd ..; done
python3 MSD_perA_v4.py "_2x"

for i in $(seq 1 $N); do cd run_${i}; ./MSD_a10b; cd ..; done
python3 MSD_perA_v4.py "_a10b"

for i in $(seq 1 $N); do cd run_${i}; ./MSD_simple; cd ..; done
python3 MSD_perA_v4.py "_simple"
