#!/bin/bash
sudo apt install build-essential libssl-dev
wget https://github.com/Kitware/CMake/releases/download/v3.16.5/cmake-3.16.5.tar.gz
tar -zxvf cmake-3.16.5.tar.gz
cd cmake-3.16.5
./bootstrap
make 
sudo make install