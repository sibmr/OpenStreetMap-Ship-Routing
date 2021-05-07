#!/bin/sh

cd build
cmake -DCMAKE_BUILD_TYPE=Release .
make -j4
make install
cd ..