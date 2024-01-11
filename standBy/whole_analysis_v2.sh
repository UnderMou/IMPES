#!/bin/bash

# rm plotdata_*
# rm data_*

g++ -o impes impes.cpp
./impes

# python3 post_proc_plots.py

# python3 post_proc_errorL2.py
