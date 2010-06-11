#!/bin/bash
#export VES3D_DIR=`pwd`

make all
cd src;
make all
cd ..
cd test
make all
cd ..

