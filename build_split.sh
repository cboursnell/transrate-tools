#! /bin/sh

g++ -std=c++11 -o bam-split -L bamtools/lib -lbamtools -I bamtools/include -lz transrate-bam-split.cpp
