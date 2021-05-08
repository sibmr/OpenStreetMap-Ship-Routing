#!/bin/sh

cmake -DCMAKE_BUILD_TYPE=Release .
make -j4
make install
