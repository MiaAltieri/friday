#!/bin/bash

# rm -rf build
mkdir build
cd build

cmake .. -DPYTHON_LIBRARY=/Applications/Anaconda/anaconda/bin/python3 && make
cd ..
python3 main.py