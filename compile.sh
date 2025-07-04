#!/bin/bash
g++ -fdiagnostics-color=always -g pp_vor_v4.cpp Mesh.cpp -o pp_vor_v4 -I/usr/local/include/voro++ -L/usr/local/lib -lvoro++
