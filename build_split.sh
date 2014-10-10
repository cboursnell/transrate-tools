#! /bin/sh

g++-4.9 -std=c++11 -g -o bam-split -L bamtools/lib -lbamtools -I bamtools/include -lz transrate-bam-split.cpp
